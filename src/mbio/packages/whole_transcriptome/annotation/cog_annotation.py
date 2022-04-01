# -*- coding: utf-8 -*-
# __author__ = 'yuanshaohua'

import argparse
import logging
import os

from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Get cog annotation from mongo with reading input BLAST tabular file')
parser.add_argument('-i', metavar='[xml_table]', required=True, help='input BLAST tabular file')
parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='Output cog annotation tabular file')
args = parser.parse_args()


class MetagenomicCogAnnotation(object):
    '''
    Obtain detailed cog annotation information of metagenomic
    '''

    def __init__(self):
        self.client = Config().get_mongo_client(mtype='metagenomic', ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname('metagenomic', ref=True)]
        self.eggNOG_ID = self.mongodb.eggNOG4_seqID
        self.eggNOG = self.mongodb.eggNOG4
        self.seq_nog = dict()

    def main(self, align_table, anno_table):
        logger.info('start building dictionary')
        nog_seq2nog_dict = {document['nog_seq']: document['nog'] for document in self.eggNOG_ID.find({})}
        logger.info('succeed in building nog dict with {} items'.format(len(nog_seq2nog_dict)))
        nog2details_dict = {document['nog']: document for document in self.eggNOG.find({})}
        logger.info('succeed in building detail dict with {} items'.format(len(nog2details_dict)))
        logger.info('start processing {}'.format(align_table))
        with open(align_table) as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tNOG\tNOG_description\tFunction\tFun_description\tCategory\tIdentity(%)\tAlign_len\n')
            NOGseq_list = dict()
            NOG_list = dict()
            for line in infile:
                line = line.strip().split('\t')
                if not 'Score' in line[0]:
                    query = line[5]
                    NOG_seqID = line[10]
                    align_len = line[2]
                    iden = line[3]
                    if not NOGseq_list.has_key(NOG_seqID):
                        nog = nog_seq2nog_dict.get(NOG_seqID)
                        if nog:
                            NOG = nog
                            NOGseq_list[NOG_seqID] = NOG
                            detail = nog2details_dict.get(NOG)
                            if not detail:
                                continue
                            if not NOG_list.has_key(NOG):
                                if detail:
                                    NOG_des = detail['nog_des']
                                    if NOG_des == 'NA':
                                        NOG_des = '-'
                                    Function = detail['function']
                                    if len(Function) > 1:
                                        Function = ';'.join([i for i in Function])
                                    Function_des = detail['function_des']
                                    Category = detail['category']
                                    All = '{}\t{}\t{}\t{}'.format(NOG_des, Function, Function_des, Category)
                                    NOG_list[NOG] = All
                                    if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                                        pass
                                    else:
                                        outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, All, iden, align_len))
                                        if self.seq_nog.has_key(query):
                                            self.seq_nog[query].append(NOG)
                                        else:
                                            self.seq_nog.update({query: [NOG]})
                                else:
                                    logger.debug('wrong NOG ID -> {}'.format(NOG))
                            else:
                                if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                                    pass
                                else:
                                    outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, All, iden, align_len))
                                    if self.seq_nog.has_key(query):
                                        self.seq_nog[query].append(NOG)
                                    else:
                                        self.seq_nog.update({query: [NOG]})
                        else:
                            logger.debug('wrong seq ID -> {}'.format(NOG_seqID))
                    else:
                        NOG = NOGseq_list[NOG_seqID]
                        All = NOG_list[NOG]
                        if self.seq_nog.has_key(query) and NOG in self.seq_nog[query]:
                            pass
                        else:
                            outfile.write('{}\t{}\t{}\t{}\t{}\n'.format(query, NOG, All, iden, align_len))
                            if self.seq_nog.has_key(query):
                                self.seq_nog[query].append(NOG)
                            else:
                                self.seq_nog.update({query: [NOG]})
        if os.path.getsize(anno_table) > 0:
            logger.info('succeed in exporting {}'.format(anno_table))


if __name__ == '__main__':
    if args.i and args.o:
        instance = MetagenomicCogAnnotation()
        instance.main(args.i, args.o)
    elif args.h:
        parser.print_help()
