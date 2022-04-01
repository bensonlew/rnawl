# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modified: 20190402


import sys
import pandas as pd
from biocluster.config import Config
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time

'''
统计go对应的gene 列表
查mongo 获取go id的ontology和描述信息
给ontology排在前面的go id做标记，标记为1
'''



def get_go_genes(infile,num,out):
    if infile == out:
        raise('输出文件和输入文件名一样，请修改输出文件名')
    client = Config().get_mongo_client(mtype="metagenomic", ref=True)
    mongodb = client[Config().get_mongo_dbname("metagenomic", ref=True)]
    global process_rerun
    process_rerun = 0

    try:
        ontology_container = {}
        fw = open(out, 'w')
        fw.write('GO (Lev1)\tGO Term\tGO ID\tSeq Number\tLs_Draw\tSeq List\n')
        data = pd.read_table(infile, sep='\t')
        data.columns = ['gene_id', 'go_id', 'gene_desc']
        group = dict(list(data.groupby(by=['go_id'])))
        for k in group.keys():
            search = mongodb['GO'].find_one({'go_id': k})
            ontology = search['ontology']
            name = search['name']
            gene_set = set(list(group[k]['gene_id']))
            if ontology not in ontology_container.keys():
                ontology_container[ontology] = [[name, k, len(gene_set), ';'.join(gene_set)]]
            else:
                ontology_container[ontology].append([name, k, len(gene_set), ';'.join(gene_set)])

        for k in sorted(ontology_container.keys()):
            s = sorted(ontology_container[k], key=lambda a: a[2], reverse=True)
            n = 0
            for i in s:
                n += 1
                if n < num:
                    fw.write('\t'.join([k, i[0], i[1], str(i[2]), '1', i[3]]) + '\n')
                else:
                    fw.write('\t'.join([k, i[0], i[1], str(i[2]), '0', i[3]]) + '\n')
    except (ServerSelectionTimeoutError, NetworkTimeout):  # 捕获因为mongo服务器问题导致的异常后重运行此方法
        if process_rerun < 5:
            process_rerun += 1
            time.sleep(5)
            get_go_genes(infile,num,out)
        else:
            raise Exception("重运行5次仍未成功连接mongo")


get_go_genes(sys.argv[1],15,sys.argv[2])


