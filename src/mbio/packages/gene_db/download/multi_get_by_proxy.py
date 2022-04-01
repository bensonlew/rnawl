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

def get_html(u_in):
    '''
    获取单一路径
    '''
    url = u_in[0]
    out = u_in[1]
    global proxys
    global proxys_fail
    retry_count = 5
    while retry_count > 0:
        proxy = proxys[random.randint(0, len(proxys)-1)]
        proxies={
            "http": "socks4://{}".format("88.199.82.66:54194"),
            "https": "socks4://{}".format("88.199.82.66:54194"),
        }
        '''
        proxies={
            "http": "socks4://{}".format(proxy),
            "https": "socks4://{}".format(proxy)
        }
        '''
        try:
            with open(out, 'w') as fo:
                headers = header()
                s = requests.session()
                s.keep_alive = False # 关闭多余连接
                data = s.get(url, headers=headers, timeout=100, proxies=proxies)
                # data = requests.get(url, headers=headers, proxies=proxies)
                fo.write(data.text)
                break

        except Exception, e:
            time.sleep(0.1)
            print "获取信息错误{} {}".format(proxy, e)
            proxys_fail.add(proxy)
            retry_count -= 1

    '''
    data = requests.get(url, timeout=1000).text

    with open(out, 'w') as fo:
        fo.write(data)
    '''

def para_get_html(multi_url_file, thread=10):
    '''
    批量下载
    '''

    global proxys
    global proxys_fail

    client = MongoClient("mongodb://{}:{}@{}/sanger_denovo_assembly?authMechanism=SCRAM-SHA-1".format("rna", "y6a3t5n1w0y7", "10.8.0.23"), connect=True)
    db = client['sanger_denovo_assembly']
    proxys = [p['proxy'] for p in db["raw_proxy"].find({"source": "freeProxy30", "type": "socks4"})]
    print proxys

    url2file = list()
    with open(multi_url_file, 'rb') as f:
        for line in f:
            url = line.strip()
            print url
            out = url.split("/")[-1] + ".ko2gene"
            url2file.append([url, out])
            get_html([url, out])
    # p = Pool(thread)
    # p.map(get_html, url2file)
    # p.close()
    # p.join()
    print "failed ip {}".format(proxys_fail)

if __name__ == "__main__":
    url_files = sys.argv[1]
    para_get_html(url_files)
