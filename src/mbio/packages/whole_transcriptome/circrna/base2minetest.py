# -*- coding: utf-8 -*-
import logging
import re
import sys

import pandas as pd
from Bio import SeqIO

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def basetomine(database, minecirc, detail, detail_new):
    circbase_seq_list = list()
    circmine_seq_list = list()
    circmine_id_seq = dict()
    circbase_id_seq = dict()
    logging.debug('start calling {}'.format(sys._getframe().f_code.co_name))
    for seq_record_circbase in SeqIO.parse(database, "fasta"):  #读取数据库中的fasta文件
        circbase_seq = seq_record_circbase.seq     #读取fasta文件中的序列
        circbase_id = seq_record_circbase.id       #读取fasta文件中的序列id
        circbase_seq_upper = str(circbase_seq.upper())   #将序列的字母全部大写
        if circbase_seq_upper not in circbase_id_seq:      #创建字典,{序列:id}
            circbase_id_seq[circbase_seq_upper] = [circbase_id]
        else:
            circbase_id_seq[circbase_seq_upper].append(circbase_id)
        circbase_seq_list.append(circbase_seq_upper)    #将所有序列组成一个列表
    logging.debug('succeed in building base dict with {} items'.format(len(circbase_id_seq)))

    for seq_record_circmine in SeqIO.parse(minecirc, "fasta"):
        circmine_seq = seq_record_circmine.seq
        circmine_id = seq_record_circmine.id
        # (chr, start, end) = re.split('[:-]', circmine_id)
        # circmine_name = chr + '_' + start + '_' + end
        circmine_name = circmine_id
        circmine_seq_upper = str(circmine_seq.upper())
        if circmine_seq_upper not in circmine_id_seq:
            circmine_id_seq[circmine_seq_upper] = [circmine_name]
        else:
            circmine_id_seq[circmine_seq_upper].append(circmine_name)
        circmine_seq_list.append(circmine_seq_upper)
    logging.debug('succeed in building mine dict with {} items'.format(len(circmine_id_seq)))

    circrna_id = list()
    database_id = list()
    circrna_databae_id = dict()
    inter_seq_set = set(circbase_id_seq.keys()) & set(circmine_id_seq.keys())
    for seq in inter_seq_set:
        circrna_id.extend(circmine_id_seq[seq])
        database_id.extend(circbase_id_seq[seq])
        circrna_databae_id[circmine_id_seq[seq][0]] = circbase_id_seq[seq][0]
    circrna_id = list(set(circrna_id))

    logging.debug('succeed in getting intersection with {} ids'.format(len(circrna_id)))
    circbase = list()
    database = list()
    df = pd.read_table(detail)
    for l in range(len(df)):
        if df['circrna_id'][l] in circrna_id:
            circbase.append('yes')
            database.append(circrna_databae_id[df['circrna_id'][l]])
        else:
            circbase.append('no')
            database.append('')
    df['circbase'] = circbase
    df['database_id'] = database
    df_none = df.fillna("")
    df_none.to_csv(detail_new, index=False, sep='\t', header=True)
    logging.debug('succeed in exporting {}'.format(detail_new))


def main(args):
    basetomine(args.database, args.minecirc, args.detail, args.detail_new)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='circBase')
    parser.add_argument('-i', action='store', required=True,
                        help='detail', dest='detail')
    parser.add_argument('-d', action='store', required=True,
                        help='database', dest='database')
    parser.add_argument('-m', action='store', required=True,
                        help='minecirc', dest='minecirc')
    parser.add_argument('-o', action='store', required=True,
                        help='output ', dest='detail_new')
    args = parser.parse_args()

    main(args)
