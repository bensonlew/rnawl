# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config

class meta_cog_anno(object):
    """
    meta宏基因cog数据详细注释信息
    """
    def __init__(self):
        #self.client = MongoClient('mongodb://10.100.200.129:27017')
        #self.mongodb = self.client.sanger_biodb
        #self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.eggNOG_ID = self.mongodb.eggNOG4_seqID
        self.eggNOG = self.mongodb.eggNOG4
        self.seq_nog = dict()

    def cog_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tNOG\tNOG_description\tFunction\tFun_description\tCategory\tIdentity(%)\tAlign_len\n')
            NOGseq_list = {}
            NOG_list = {} 
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    #query = line[0]
                    #query = ('_').join(line[5].split("_")[:-1])
                    if line[5].endswith("_1"):
                        query = line[5].rsplit("_1",1)[0]
                    else:
                        query = line[5]
                    NOG_seqID = line[10]
                    align_len = line[2]
                    iden = line[3]
                    if not NOGseq_list.has_key(NOG_seqID):
                        detail = self.eggNOG_ID.find_one({"nog_seq":NOG_seqID})
                        if detail:
                            NOG = detail["nog"]
                            NOGseq_list[NOG_seqID] = NOG
                            detail2 = self.eggNOG.find_one({"nog":NOG})
                            if not NOG_list.has_key(NOG):
                                if detail2:
                                    NOG_des = detail2["nog_des"]
                                    if NOG_des == "NA":
                                        NOG_des = "-"
                                    Function = detail2["function"]
                                    if len(Function) > 1:
                                        Function = ";".join([i for i in Function])
                                        #Function1 = Function.split("")
                                        #Function = ";".join(Function1)
                                    Function_des = detail2["function_des"]
                                    Category = detail2["category"]
                                    all_des = '{}\t{}\t{}\t{}'.format(NOG_des, Function, Function_des, Category)
                                    NOG_list[NOG] = all_des
                                    if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                                        pass
                                    else:
                                        outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, all_des, iden, align_len))
                                        if self.seq_nog.has_key(query):
                                            self.seq_nog[query].append(NOG)
                                        else:
                                            self.seq_nog.update({query: [NOG]})
                                else:
                                    print "wrong NOG ID",NOG
                            else:
                                if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                                    pass
                                else:
                                    all_des = NOG_list[NOG]
                                    outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, all_des, iden, align_len))
                                    if self.seq_nog.has_key(query):
                                        self.seq_nog[query].append(NOG)
                                    else:
                                        self.seq_nog.update({query: [NOG]})
                        else:
                            print "wrong seq ID",NOG_seqID
                    else:
                        NOG = NOGseq_list[NOG_seqID]
                        all_des = NOG_list[NOG]
                        if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                            pass
                        else:
                            outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, all_des, iden, align_len))
                            if self.seq_nog.has_key(query):
                                self.seq_nog[query].append(NOG)
                            else:
                                self.seq_nog.update({query: [NOG]})

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_cog_anno = meta_cog_anno()
    meta_cog_anno.cog_by_mongo(align_table, anno_table)
