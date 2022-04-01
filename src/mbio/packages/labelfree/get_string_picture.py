#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/13 14:44
@file    : get_string_picture.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


"""这个脚本是基于 https://string-db.org/cgi/help.pl?sessionId=7dsPPkRbfFV3&subpage=api上面的内容写的，string网站提供了一系列的api可以
方便大家爬虫，之前都不知道。我看了一下，大体可以做到爬图片，爬对应关系，爬id转化等，这些大体可以满足我们的要求，但是其中线下分析中提供的domain信息
我没有找到，以后找到可以再添加，这个脚本只能本地运行，基本上可以实现不去string官网，直接运行脚本就得到string图片等信息，如果有啥要求想加的话，可以
联系yitong.feng@majorbio.com"""


import os
import shutil
import psutil
import signal
import re
import sys
import argparse
import time
import random
import glob
from Bio.Blast import NCBIXML
import requests
import urllib
import urllib2
import ssl
import pandas as pd
import functools
from collections import Counter
ssl._create_default_https_context = ssl._create_unverified_context   #处理了ssl证书出错的问题

def time_out(func):
    def kill(pid):
        p = psutil.Process(os.getpid())
        if p.name == 'kill':
            return
        cc = p.children()
        if not cc:
            os.kill(pid, signal.SIGKILL)
        else:
            for child in cc:
                kill(child.pid)
    @functools.wraps(func)
    def wrapper(*args, **kw):
        func_name = func.__name__
        pid = os.fork()
        rel = None
        if pid == 0:
            time.sleep(600)
            print('函数%s运行了10分钟后没有响应，即将被杀掉'%func_name)
            parent = psutil.Process(os.getppid())
            # os.kill(parent.pid, signal.SIGKILL)
            # os._exit()
            kill(parent.pid)
        else:
            start = time.time()
            start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
            print('Run ' + func_name + ' at ' + start_time)
            rel = func(*args, **kw)
            end = time.time()
            end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
            print('End ' + func_name + ' at ' + end_time)
            print("{}函数执行时间约为{}s".format(func.__name__, end - start))
            # child = psutil.Process(os.getpid()).children()[0]
            for child in psutil.Process(os.getpid()).children():
                kill(child.pid)
            #     print(child.pid,child.status())
            #     child.kill
            #     print(child.status())
            #     child.terminal()
            #     print(child.status())
            #     os.kill(child.pid, signal.SIGKILL)
            # print(child.status())
        return rel
    return wrapper

class string_api(object):
    def __init__(self, specie, list_path, vsstring, identity, useblast, max_num, out):
        self.specie = specie
        self.list_files = glob.glob(os.path.join(list_path,'*.list'))
        if not self.list_files:
            exit('路径下需要存在用来做string的accession号列表文件')
        self.vstringfile = ''
        if not self.specie and not os.path.exists(vsstring):
            exit('%s在string数据库里面没有收录，请传入string的blast文件'%str(specie))
        self.vstringfile = os.path.abspath(vsstring)
        self.identity = float(identity)
        self.useblast = useblast.lower()
        self.max_num = max_num
        if not os.path.exists(out):
            exit('输出目录不存在啊')
        self.out = out

    def deal_blast(self):
        acc2string_identity = dict()
        if self.vstringfile.endswith('.xls'):
            with open(self.vstringfile, 'r') as blast_r:
                for line in blast_r:
                    if line.strip():
                        items = line.strip().split('\t')
                        if float(items[3]) >= self.identity:
                            if items[0] not in acc2string_identity:
                                acc2string_identity[items[0]] = list()
                            acc2string_identity[items[0]].append((items[1],items[3]))
        else:
            with open(self.vstringfile, 'r') as blast_r:
                records = NCBIXML.parse(blast_r)
                for rec in records:
                    query = re.split(' ', rec.query, maxsplit=1)[0]
                    for align in rec.alignments:
                        for hsp in align.hsps:
                            ident = float(hsp.identities)
                            hit_len = float(hsp.align_length)
                            identity = ident / hit_len * 100
                            if identity < self.identity:
                                continue
                            hit = align.hit_def
                            if u'(' in hit or u'|' in hit:
                                hit = align.hit_id
                            if hit.endswith('-P'):
                                hit = hit.split('-P')[0]
                            if query not in acc2string_identity:
                                acc2string_identity[query] = list()
                            acc2string_identity[query].append((hit, identity))
        return acc2string_identity

    def get_stringid_from_file(self, list_file, acc2string_identity):
        with open(list_file) as lf:
            acc_list = [line.strip().split('\t')[0] for line in lf if line.strip()]
        if self.useblast != 'yes':
            acc_dict = {acc:acc for acc in acc_list}
            return acc_dict
        # string_dict = dict()
        def filter_dict():
            string_dict = dict()
            for acc in acc_list:
                if acc in acc2string_identity:
                    string,identity = acc2string_identity[acc][0]
                    if identity < self.identity:
                        continue
                    for s,i in acc2string_identity[acc]:
                        if i < self.identity:
                            continue
                        if i > identity:
                            identity = i
                            string = s
                    string_dict[string] = acc
            return string_dict
        string_dict = filter_dict()
        while len(string_dict) > self.max_num and self.identity < 99.9:
            # print(self.identity)
            # print(len(string_dict))
            self.identity += 0.01
            string_dict = filter_dict()
        return string_dict

    # @time_out
    def get_protein_info_from_net(self, string_dict, out):
        string2acc = dict()
        id_lists = list()
        l = len(string_dict.keys())
        for i in range(0, l, 50):
            if i + 50 <= l:
                id_lists.append(string_dict.keys()[i:i+50])
            else:
                id_lists.append(string_dict.keys()[i:l])
        time.sleep(1+random.random())
        string_api_url = "http://string-db.org/api"
        output_format = "tsv-no-header"
        method = "get_string_ids"
        identifiers = list()
        for n,ids in enumerate(id_lists):
            time.sleep(1 + random.random())
            params = {
                "identifiers": "\r".join(ids),
            # your protein list
                # "species" : 9606, # species NCBI identifier
                "limit": 1,  # only one (best) identifier per input protein
                "echo_query": 1,  # see your input identifiers in the output
                # "caller_identity": "www.awesome_app.org"  # your app name
            }
            print('https://string-db.org/api/tsv/get_string_ids?identifiers=' + '%0d'.join(params['identifiers'].split('\r')))
            if self.specie:
                params.update({'species':self.specie})
            ## contruct method URL
            request_url = string_api_url + "/" + output_format + "/" + method
            ## Call STRING
            response = requests.post(request_url, params=params)
            print("https://string-db.org/api/tsv/get_string_ids request url and params: "+str(request_url)+"\n"+str(params))
            print("https://string-db.org/api/tsv/get_string_ids response: "+str(response.text))

            if response.status_code == 400:
                continue
            else:
                pass

            if n == 0:
                ow = open(out, 'w')
                ow.write('#node\taccession\tstring_id\tannotation\n')
            else:
                ow = open(out, 'a')
                if os.path.getsize(out) == 0:   # dealing with empty result for the first set of ids
                    ow.write('#node\taccession\tstring_id\tannotation\n')
            for line in response.text.strip().split("\n"):
                if not line:
                    continue
                l = line.split("\t")
                input_identifier, string_identifier, annotation, node = l[0], l[2], l[-1], l[-2]
                string2acc[string_identifier] = string_dict[input_identifier]
                ow.write(node + '\t' + string_dict[input_identifier] + '\t' + string_identifier + '\t' + annotation + '\n')
                if string_identifier not in identifiers:
                    identifiers.append(string_identifier)
            ow.close()
        print(string2acc)
        return identifiers, string2acc

    # @time_out
    def get_identifiers_bitscore_and_filter(self, identifiers, string2acc, out):
        id_lists = list()
        l = len(identifiers)
        for i in range(0, l, 50):
            if i + 50 <= l:
                id_lists.append(identifiers[i:i+50])
            else:
                id_lists.append(identifiers[i:l])
        time.sleep(1+random.random())
        string_api_url = "http://string-db.org/api"
        output_format = "tsv-no-header"
        method = "homology"
        for n,ids in enumerate(id_lists):
            time.sleep(1 + random.random())
            params = {
                "identifiers": "\r".join(ids),
            # your protein list
                # "species" : 9606, # species NCBI identifier
                "limit": 1,  # only one (best) identifier per input protein
                "echo_query": 1,  # see your input identifiers in the output
                # "caller_identity": "www.awesome_app.org"  # your app name
            }
            print('https://string-db.org/api/tsv/homology?identifiers=' + '%0d'.join(params['identifiers'].split('\r')))
            if self.specie:
                params.update({'species':self.specie})
            ## contruct method URL
            request_url = string_api_url + "/" + output_format + "/" + method
            ## Call STRING
            response = requests.post(request_url, params=params)
            print("string-db.org request url and params: "+str(request_url)+"\n"+str(params))
            print("string-db.org response: "+str(response.text))

            if n == 0:
                ow = open(out, 'w')
                ow.write('ncbiTaxonId_A\tstringId_A\tncbiTaxonId_B\tstringId_B\tbitscore\n')
            else:
                ow = open(out, 'a')
            for line in response.text.strip().split("\n"):
                if not line:
                    continue
                ow.write(line + '\n')
            ow.close()
        score_df = pd.read_csv(out, sep='\t')
        for ind in score_df.index:
            sa = score_df.loc[ind, 'stringId_A']
            sb = score_df.loc[ind, 'stringId_B']
            pa, pb = '', ''
            for string, acc in string2acc.items():
                if '.' + sa in string:
                    pa = acc
                if '.' + sb in string:
                    pb = acc
            score_df.loc[ind, 'accession_A'] = pa
            score_df.loc[ind, 'accession_B'] = pb
        score_df = score_df[['ncbiTaxonId_A', 'stringId_A', 'accession_A', 'ncbiTaxonId_B', 'stringId_B', 'accession_B', 'bitscore']]
        score_df.to_csv(out, sep='\t', index=False, header=True)
        while len(list(set(score_df['stringId_A'].tolist()+score_df['stringId_B'].tolist()))) > 100:
            min_ = score_df['bitscore'].min()
            score_df = score_df[score_df['bitscore']>(min_+1)]
        sa = [str(a) + '.' + b for a,b in zip(score_df['ncbiTaxonId_A'].tolist(), score_df['stringId_A'].tolist())]
        sb = [str(a) + '.' + b for a,b in zip(score_df['ncbiTaxonId_B'].tolist(), score_df['stringId_B'].tolist())]
        identifiers = list(set(sa + sb))
        # 从202012起string网发生了变化，不需要再额外加taxid了，很惆怅
        for i in range(len(identifiers)):
            if identifiers[i].count('.') > 1:
                identifiers[i] = identifiers[i].split('.', 1)[1]
        return identifiers

    # @time_out
    def get_picture_from_net(self, identifiers, out):
        # print(identifiers)
        time.sleep(1+random.random())
        string_api_url = "https://string-db.org/api"
        output_format = "svg"
        method = "network"

        ## Construct the request
        my_app = "www.awesome_app.org"
        request_url = string_api_url + "/" + output_format + "/" + method + "?"
        request_url += "identifiers=%s"
        if self.specie:
            request_url += "&" + "species=" + str(self.specie)
        request_url += "&" + "network_flavor=" + 'confidence'
        request_url += "&" + "required_score=30"
        # request_url += "&" + "add_white_nodes=5"
        # request_url += "&" + "caller_identity=" + my_app

        ## For each gene call STRING
        print(request_url % "%0d".join(identifiers))
        file_name, _ = urllib.urlretrieve(request_url % "%0d".join(identifiers), "%s.network.svg" % out)
        # time.sleep(1 + random.random())
        # file_name, _ = urllib.urlretrieve(request_url % "%0d".join(identifiers).replace('/svg/', '/highres_image/'), "%s.png" % out)
        # cmd = 'convert -flatten -quality 100 -density 130 -background white %s %s' % (file_name, "%s.png" % out)
        # os.system(cmd)

    # @time_out
    def get_interaction_from_net(self, identifiers, string2acc, out):
        time.sleep(1 + random.random())
        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "network"
        ## Construct the request
        request_url = string_api_url + "/" + output_format + "/" + method + "?"
        request_url += "identifiers=%s" % "%0d".join(identifiers)
        if self.specie:
            request_url += "&" + "species=" + str(self.specie)
        # request_url += "&" + "caller_identity=" + my_app

        try:
            response = urllib2.urlopen(request_url)
        except urllib2.HTTPError as err:
            error_message = err.read()
            print(error_message)
            sys.exit()

        ## Read and parse the results

        with open(out, 'w') as ow:
            ow.write('Accession (protein A)\tAccession (protein B)\tSTRING identifier (protein A)\tSTRING identifier (protein B)\tcommon protein name (protein A)\tcommon protein name (protein B)\tNCBI taxon identifier\tcombined score\tgene neighborhood score\tgene fusion score\tphylogenetic profile score\tcoexpression score\texperimental score\tdatabase score\ttextmining score\n')
            # ow.write(response.read())
            for line in response.readlines():
                if not line.strip():
                    continue

                tmp = line.strip().split('\t')
                sa = tmp[0]
                sb = tmp[1]
                pa, pb = '', ''
                for string, acc in string2acc.items():
                    if '.' + sa in string or string == sa:
                        pa = acc
                    if '.' + sb in string or string == sb:
                        pb = acc
                ow.write(pa + '\t' + pb + '\t' + line)
                # try:
                #     ow.write(string2acc[tmp[0]]+ '\t' + string2acc[tmp[1]] + '\t' +line)
                # except:
                #     print(line)

    def run(self):
        acc2string_identity = dict()
        if self.useblast == 'yes':
            acc2string_identity = self.deal_blast()
        for file in self.list_files:
            if self.list_files.index(file) != 0:
                time.sleep(60+random.randint(1,18))
            prefix = os.path.basename(file).split('_protein.list')[0]
            out = os.path.join(self.out, prefix)
            if os.path.exists(out):
                shutil.rmtree(out)
            os.makedirs(out)
            out = os.path.join(out, prefix)
            print("1 out is:"+str(out))
            print("2 file is:"+str(file))
            print("3 acc2string_identity is:"+str(acc2string_identity))
            string2acc = self.get_stringid_from_file(file, acc2string_identity)
            print("4 string2acc is:"+str(string2acc))
            identifiers,string2acc = self.get_protein_info_from_net(string2acc, out+'.annotation.xls')
            print("5 identifiers is:"+str(identifiers))
            print("6 string2acc is:"+str(string2acc))
            identifiers = self.get_identifiers_bitscore_and_filter(identifiers, string2acc, out+'.bitscore.xls')
            print("7 identifiers is:"+str(identifiers))
            self.get_picture_from_net(identifiers, out)
            self.get_interaction_from_net(identifiers, string2acc, out+'.interaction.xls')


if __name__ == '__main__':
    # def __init__(self, specie, list_path, vsstring, identity, useblast):
    parser = argparse.ArgumentParser(description="use the string net api to get the picture and no need to go internet")
    parser.add_argument("-specie", type=int, default=0, help="your species taxno number")
    parser.add_argument("-list_path", type=str, required=True, help="the dir path contains your list files which would be used for string picture")
    parser.add_argument("-vsstring", type=str, default='', help="the string blast result,if not use,it can be space")
    parser.add_argument("-identity", type=float, default=98,help='the filter para to filter the string blast result')
    parser.add_argument("-max_num", type=int, default=300,help='the filter para to filter the string blast result')
    parser.add_argument("-useblast", type=str, default='yes',help='you should judge whether to use the blast result,default yes')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'string_pictures'),help='the output dir of string pictures')

    args = parser.parse_args()
    STRING = string_api(args.specie, args.list_path, args.vsstring, args.identity, args.useblast, args.max_num, args.out)
    STRING.run()
