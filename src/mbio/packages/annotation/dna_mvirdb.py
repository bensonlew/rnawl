# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'


import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time



class meta_mvirdb_anno(object):
    """
    dna mvirdb数据详细注释信息
    ardb.parse.anno.xls
    """

    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.mongodb = self.client.sanger_biodb
        # self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.process_rerun = 0
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.mvirdb = self.mongodb.mvirdb

    def mvirdb_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write(
             'Gene ID\tVirulence Factor ID\tDatabase Source\tUser Designated Name\tGene Name\tShort Description\tVirulence Factor Type\tstatus\tIdentity(%)\tEvalue\tScore\tAlign_len\n')

            # query_list = []
            # done_list = []
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    # query = line[0]
                    query = line[5]
                    hit = line[10]
                    id = hit.split("|")[0]
                    evalue = line[1]
                    score = line[0]
                    align_len = line[11]
                    act_iden = line[3]
                    try:
                        detail = self.mvirdb.find_one({"virulence_id": id})
                    except (ServerSelectionTimeoutError, NetworkTimeout):  # 捕获因为mongo服务器问题导致的异常后重运行此方法
                        if self.process_rerun < 5:
                            self.process_rerun += 1
                            # self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                            time.sleep(5)
                            self.mvirdb_by_mongo(align_table, anno_table)
                        else:
                            raise Exception("重运行5次仍未成功连接mongo")
                    if detail:
                        data_source = detail["source"]
                        designated_name=detail['designated_name']
                        gene_name=detail['gene_name']
                        type=detail['type']
                        description=detail['description']
                        status=detail['status']


                        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, id, data_source,designated_name, gene_name,
                                                         description, type, status,act_iden, evalue, score, align_len))
                    else:
                        print "wrong ID"



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_mvirdb_anno = meta_mvirdb_anno()
    meta_mvirdb_anno.mvirdb_by_mongo(align_table, anno_table)
