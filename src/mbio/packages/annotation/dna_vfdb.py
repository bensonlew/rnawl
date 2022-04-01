# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'


import sys
import os
from pymongo import MongoClient
import argparse
from biocluster.config import Config
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time


class meta_vfdb_anno(object):
    """
    dna vfdb数据详细注释信息
    "Gene\tType\tIdentity\tEvalue\tResistance\tRequirement\tClass\tClass_Description\n"
    """

    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.client = MongoClient('mongodb://192.168.10.187:27017')
        # self.client = MongoClient(Config().MONGO_BIO_URI)
        # self.mongodb = self.client.sanger_biodb
        self.process_rerun = 0
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.vfdb = self.mongodb['vfdb_v20200703']

    def vfdb_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('Gene ID\tVFDB ID\tVFs\tSpecies\tDescription\tLevel1\tLevel2\tIdentity(%)\tEvalue\tScore\tCoverage(%)\n')
            l2 = {}
            l2_l1={}
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    query = line[5]
                    # print line[0],"\n",query
                    Subject_ID = line[10]
                    align_len = line[2]
                    VFG_ID = Subject_ID.split("(")[0]
                    # print ardb_ID
                    evalue = line[1]
                    act_iden = line[3]
                    coverge = round(abs(float(line[8])-float(line[7]))/float(line[6]),3)
                    coverge = coverge*100
                    try:
                        detail = self.vfdb.find_one({"vfg_id": VFG_ID})
                    except (ServerSelectionTimeoutError, NetworkTimeout):  # 捕获因为mongo服务器问题导致的异常后重运行此方法
                        if self.process_rerun < 5:
                            self.process_rerun += 1
                            # self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                            time.sleep(5)
                            self.vfdb_by_mongo(align_table, anno_table)
                        else:
                            raise Exception("重运行5次仍未成功连接mongo")
                    if detail:
                        Gi = detail["gi_number"]
                        gene = detail["gene"]
                        # if isinstance(resistance,list) == True:
                        #   resistance = ','.join(resistance)
                        gene_des = detail["gene_des"]
                        species = detail["origin"]
                        VF = detail["vfs"]
                        Function = detail["vf_function"]
                        level1_list = detail["level1"].strip().split(";") ## fix by qingchen.zhang
                        level2_list = detail["level2"].strip().split(";") ## fix by qingchen.zhang
                        for level2 in level2_list:## vfdb改成了一对多关系，要改正此对应关系
                            index_num = level2_list.index(level2)
                            level1 = level1_list[index_num]
                            if level2 in l2:
                                l2[level2] += 1
                            else:
                                l2[level2] = 1
                                l2_l1[level2]=level1
                        outfile.write(
                            '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, Subject_ID, VF, species, Function,
                                                                              ";".join(level1_list),
                                                                                  ";".join(level2_list), act_iden,
                                                                                  line[1], line[0], coverge))
                    else:
                        print
                        "wrong ID"
        level_file = os.path.dirname(anno_table) + '/vfdb_level.xls'
        with open(level_file, 'wb') as outfile1:
            outfile1.write('Level1\tLevel2\tGene No.\n')
            for key in l2:
                if key != '-':
                    outfile1.write('{}\t{}\t{}\n'.format(l2_l1[key],key,l2[key]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    # parser.add_argument('-iden',metavar='[identity]',help='input ardb identity cutoff')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_vfdb_anno = meta_vfdb_anno()
    meta_vfdb_anno.vfdb_by_mongo(align_table, anno_table)
