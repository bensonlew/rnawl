# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# 通过代理爬取多个网页

import re
import os
import sys
import xml.etree.ElementTree as ET
import lxml.html
import chardet
import random
import requests
from pymongo import MongoClient
from multiprocessing import Pool
import time
from lxml import etree
import json

requests.adapters.DEFAULT_RETRIES = 5
proxys = list()
proxys_fail = set()


def user_agent():
    """
    return an User-Agent at random
    :return:
    """
    ua_list = [
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/30.0.1599.101',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/38.0.2125.122',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.71',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95',
        'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.71',
        'Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1; QQDownload 732; .NET4.0C; .NET4.0E)',
        'Mozilla/5.0 (Windows NT 5.1; U; en; rv:1.8.1) Gecko/20061208 Firefox/2.0.0 Opera 9.50',
        'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:34.0) Gecko/20100101 Firefox/34.0',
    ]
    return random.choice(ua_list)

def header():
    """
    basic header
    :return:
    """
    return {'User-Agent': user_agent(),
            'Accept': '*/*',
            'Connection': 'close',
            'Accept-Language': 'zh-CN,zh;q=0.8'}

def get_one_html(u_in):
    '''
    获取单一路径
    '''
    url = u_in
    # proxy = u_in[1]
    retry_count = 5
    # proxies={
    #     "http": "socks4://{}".format(proxy),
    #     "https": "socks4://{}".format(proxy)
    # }


    result = ""
    while retry_count > 0:
        try:
            headers = header()
            s = requests.session()
            s.keep_alive = False # 关闭多余连接
            data = s.get(url, headers=headers, timeout=100)
            # data = s.get(url, headers=headers, timeout=100, proxies=proxies)
            # data = requests.get(url, headers=headers, proxies=proxies)
            result = data.text
            break

        except Exception, e:
            time.sleep(2)
            print "获取信息错误{} {}".format(u_in, e)
            retry_count -= 1

    return result
    '''
    data = requests.get(url, timeout=1000).text

    with open(out, 'w') as fo:
        fo.write(data)
    '''

def wiki_xml2dict(xml_text):
    '''
    convert wiki html to json
    '''
    des_dict = dict()
    html = lxml.html.parse(xml_text)
    root = html.getroot()
    body = root.find('body')
    content =  body.xpath("//div[@id='mw-content-text']")
    content1 = content[0]
    cs = content1.getchildren()
    c1 = cs[0]
    ds = c1.getchildren()
    att = "default"
    for d in ds:
        if d.tag == "h2":
            att = d.getchildren()[0].text
        elif d.tag == "p":
            es = d.getchildren()
            if d.text:
                des = d.text
            else:
                des = ""
            for e in es:
                if e.text:
                    des += str(e.text)
                if e.tail:
                    des += str(e.tail)
            des_dict[att] = des
    return des_dict

def wiki_all(gene_list, name="gene_list"):
    if os.path.exists("tmp_wiki"):
        pass
    else:
        os.mkdir("tmp_wiki")
    all_dict = dict()
    for gene in gene_list:
        if "/" in gene:
            continue
        gene = gene.upper()
        if os.path.exists("tmp_wiki/" + gene + ".xml"):
            pass
        else:
            xml_text = get_one_html("https://en.wikipedia.org/wiki/{}".format(gene.upper()))
            with open("tmp_wiki/" + gene + ".xml", 'w') as fo:
                fo.write(xml_text.encode('utf8'))
        try:
            a_dict = wiki_xml2dict("tmp_wiki/" + gene + ".xml")

            all_dict[gene] = a_dict
        except:
            pass

    with open(name + ".json", 'w') as f:
        f.write(json.dumps(all_dict, indent=4, ensure_ascii=False).encode('utf8'))



def para_get_html(multi_url_file, fa_file):
    '''
    批量下载
    '''

    global proxys
    global proxys_fail

    client = MongoClient("mongodb://{}:{}@{}/sanger_denovo_assembly?authMechanism=SCRAM-SHA-1".format("rna", "y6a3t5n1w0y7", "10.8.0.23"), connect=True)
    db = client['sanger_denovo_assembly']
    proxys = [p['proxy'] for p in db["raw_proxy"].find({"source": "freeProxy30", "type": "socks4"})]
    proxy = proxys[random.randint(0, len(proxys)-1)]

    url2file = list()
    fa_num = 0
    with open(multi_url_file, 'rb') as f, open(fa_file + ".fna", 'w') as fw1, open(fa_file + ".nofaa", 'w') as fw2:
        for line in f:
            url = line.strip()
            text = get_one_html([url, proxy])
            if text == "":
                pass
                # print "文件为空， 代理{}可能失效"
            else:
                fa_num += 1
                tree = etree.HTML(text)
                a = tree.xpath('//pre')[0]
                if a.getchildren()[0].text.startswith("No"):
                    fw2.write(url + "\n")
                else:
                    fw1.write(a.getchildren()[1].tail)

    if fa_num == 0:
        db["raw_proxy"].update({"proxy": proxy}, {"$set": {"fail_count": 1}},  upsert=True)

if __name__ == "__main__":
    # url_file = sys.argv[1]
    # fa_file = sys.argv[2]
    # para_get_html(url_file, fa_file)
    # html_file = sys.argv[1]
    # a = wiki_xml2dict(html_file)
    # print json.dumps(a, indent=4, ensure_ascii=False).encode('utf8')
    gene_list = list()
    with open(sys.argv[1], 'r') as f:
        [gene_list.append(g.strip()) for g in f]
    wiki_all(gene_list, name=sys.argv[1])
