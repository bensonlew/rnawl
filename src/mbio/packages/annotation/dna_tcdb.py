# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'


import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time


class meta_tcdb_anno(object):
    """
    dna tcdb数据详细注释信息
    ardb.parse.anno.xls
    """

    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.mongodb = self.client.sanger_biodb
        # self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.process_rerun = 0
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.tcdb = self.mongodb.TCDB_v20200917

    def tcdb_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write(
                'Gene ID\tTCDB ID\tTCDB Description\tTCDB Family\tTCDB Subclass\tTCDB Class\tIdentity(%)\tEvalue\tScore\tAlign_len\n')
            #guanqing.zou 20181016  宏基因组结果增加比对长度 Align_len ，注：细菌也调用此package，故前面的列的顺序不能轻易改变
            # query_list = []
            # done_list = []
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    # query = line[0]
                    query = line[5]
                    hit = line[10]
                    id = hit.split("|")[-1]
                    evalue = line[1]
                    score = line[0]
                    align_len = line[11]
                    act_iden = line[3]
                    try:
                        detail = self.tcdb.find_one({"_id": id})
                    except (ServerSelectionTimeoutError, NetworkTimeout):  # 捕获因为mongo服务器问题导致的异常后重运行此方法
                        if self.process_rerun < 5:
                            self.process_rerun += 1
                            # self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                            time.sleep(5)
                            self.tcdb_by_mongo(align_table, anno_table)
                        else:
                            raise Exception("重运行5次仍未成功连接mongo")
                    if detail:
                        subclass = detail["subclass"]
                        class_name = detail["class"]
                        des = detail["des"]
                        family = detail["family"]
                        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, id, des,
                                                                                        family, subclass,
                                                                                        class_name, act_iden,
                                                                                        evalue, score, align_len))
                    else:
                        print "wrong ID"


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_tcdb_anno = meta_tcdb_anno()
    meta_tcdb_anno.tcdb_by_mongo(align_table, anno_table)
