# -*- coding: utf-8 -*-
# __author__ = 'yuanshaohua'

import argparse
import os
import time
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
from biocluster.config import Config

parser = argparse.ArgumentParser(description='Get cog annotation from mongo by reading input BLAST tabular file')
parser.add_argument('-i', metavar='[xml_table]', required=True, help='input BLAST tabular file')
parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='Output cog annotation tabular file')
args = parser.parse_args()


class MetagenomicCogAnnotation(object):
    '''
    Obtain detailed cog annotation information of metagenomic
    '''

    def __init__(self):
        self.process_rerun = 0
        try:
            self.connect_mongo()
            print '成功链接mongo库'
        except (ServerSelectionTimeoutError, NetworkTimeout):  # 捕获因为mongo服务器问题导致的异常后重运行此方法
            if self.process_rerun < 5:
                self.process_rerun += 1
                print "检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun)
                time.sleep(5)
                self.connect_mongo()
            else:
                raise Exception("链接失败")
        self.seq_nog = dict()

    def connect_mongo(self):
        self.client = Config().get_mongo_client(mtype='metagenomic', ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname('metagenomic', ref=True)]
        self.eggNOG_ID = self.mongodb.eggNOG4_seqID
        self.eggNOG = self.mongodb.eggNOG4

    def main(self, align_table, anno_table):
        print 'INFO: start building dictionary'
        nog_seq2nog_dict = {document['nog_seq']: document['nog'] for document in self.eggNOG_ID.find({})}
        nog2detail2_dict = {document['nog']: document for document in self.eggNOG.find({})}
        print 'INFO: start processing {}'.format(align_table)
        with open(align_table) as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tNOG\tNOG_description\tFunction\tFun_description\tCategory\tIdentity(%)\tAlign_len\n')
            NOGseq_list = dict()
            NOG_list = dict()
            for line in infile:
                line = line.strip().split('\t')
                if not 'Score' in line[0]:
                    # if line[5].endswith('_1'):
                    #     query = line[5].rsplit('_1',1)[0]
                    # else:
                    query = line[5]
                    NOG_seqID = line[10]
                    align_len = line[2]
                    iden = line[3]
                    if not NOGseq_list.has_key(NOG_seqID):
                        # detail = self.eggNOG_ID.find_one({'nog_seq':NOG_seqID})
                        nog = nog_seq2nog_dict.get(NOG_seqID)
                        # if detail:
                        #     NOG = detail['nog']
                        if nog:
                            NOG = nog
                            NOGseq_list[NOG_seqID] = NOG
                            # detail2 = self.eggNOG.find_one({'nog':NOG})
                            detail2 = nog2detail2_dict.get(NOG)
                            if not detail2:
                                continue
                            if not NOG_list.has_key(NOG):
                                if detail2:
                                    NOG_des = detail2['nog_des']
                                    if NOG_des == 'NA':
                                        NOG_des = '-'
                                    Function = detail2['function']
                                    if len(Function) > 1:
                                        Function = ';'.join([i for i in Function])
                                    Function_des = detail2['function_des']
                                    Category = detail2['category']
                                    all = '{}\t{}\t{}\t{}'.format(NOG_des, Function, Function_des, Category)
                                    NOG_list[NOG] = all
                                    if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                                        pass
                                    else:
                                        outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, all, iden, align_len))
                                        if self.seq_nog.has_key(query):
                                            self.seq_nog[query].append(NOG)
                                        else:
                                            self.seq_nog.update({query: [NOG]})
                                else:
                                    print 'wrong NOG ID -> {}'.format(NOG)
                            else:
                                if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                                    pass
                                else:
                                    outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, all, iden, align_len))
                                    if self.seq_nog.has_key(query):
                                        self.seq_nog[query].append(NOG)
                                    else:
                                        self.seq_nog.update({query: [NOG]})
                        else:
                            print 'wrong seq ID -> {}'.format(NOG_seqID)
                    else:
                        NOG = NOGseq_list[NOG_seqID]
                        all = NOG_list[NOG]
                        if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                            pass
                        else:
                            outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, all, iden, align_len))
                            if self.seq_nog.has_key(query):
                                self.seq_nog[query].append(NOG)
                            else:
                                self.seq_nog.update({query: [NOG]})
        if os.path.getsize(anno_table) > 0:
            print 'INFO: succeed in exporting {}'.format(anno_table)


if __name__ == '__main__':
    if args.i and args.o:
        instance = MetagenomicCogAnnotation()
        instance.main(args.i, args.o)
    elif args.h:
        parser.print_help()
