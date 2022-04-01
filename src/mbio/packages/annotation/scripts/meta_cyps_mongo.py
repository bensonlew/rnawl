# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang"

import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config

class meta_cyps_anno(object):
    """
    宏基因组p450数据详细注释信息
    Gene\tIdentity\tEvalue\tsid\thfid\tsfid\tgi\tsp\thomo_family\tsuper_family
    """
    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        #self.mongodb = self.client.sanger_biodb
        #self.mongodb = Config().mongo_client.sanger_biodb
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.cyps = self.mongodb.cyps

    def cyps_by_mongo(self, align_table, anno_table):
        """
        从mongo参考库cyps注释
        """
        with open(align_table, "r") as infile, open(anno_table, "wb") as outfile:
            outfile.write("#Query\tSid\tIdentity(%)\tAlign_len\thfam_id\tsfam_id\tnr_hit\tspecies\thomo_family\tsuper_family\n")
            for line in infile:
                line = line.strip().split('\t')
                if not "Score" in line[0]:
                    query = ('_').join(line[5].split('_')[:-1])
                    sid = line[10].split('_')[1]
                    act_len = line[2]
                    act_iden = line[3]
                    detail = self.cyps.find_one({"sid": sid})
                    if detail:
                        hfid = detail["hfid"]
                        sfid = detail["sfid"]
                        gi = "gi|"+ str(detail["gi"])
                        sp = detail.get("sp")
                        homo = detail["homo_family"]
                        super = detail["super_family"]
                        if sp:
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, sid, act_iden, act_len, hfid, sfid, gi, sp, homo, super))
                        else:
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, sid, act_iden, act_len, hfid, sfid, gi, "-", homo, super))

                    else:
                        print "wrong ID"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_cyps_anno = meta_cyps_anno()
    meta_cyps_anno.cyps_by_mongo(align_table, anno_table)


