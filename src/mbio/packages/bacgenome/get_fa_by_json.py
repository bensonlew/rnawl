# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os,re
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
from biocluster.file import download
from bson.objectid import ObjectId
from pymongo import MongoClient
from biocluster.api.database.base import Base
from biocluster.config import Config

my_output_fa_record = []
my_output_log_record = pd.DataFrame()

class FromMongo(Base):
    def __init__(self, chr):
        super(FromMongo, self).__init__()
        self._project_type = "bac_assem"
        self.chr = chr

    def parse_mongo(self, id):
        id = ObjectId(id)
        coll = self.db["gap_newseq"]
        result = coll.find_one({"_id": id})
        if "/" in result["seq"]:
            self.parse_file(result)
        else:
            self.parse_atgc(result)
        self.log = {"loc": self.chr, "seq": result["seq_id"], "from": result["seq_id"], "is_circle": "False"}

    def parse_file(self, result):
        remote_path = result["seq"]
        file_name = os.path.basename(remote_path)
        fa_path = download(remote_path, file_name)
        seq_tmp = SeqIO.parse(fa_path, "fasta")
        for record in seq_tmp:
            record.id = result["seq_id"]
            record.name = result["seq_id"]
            record.description=result["seq_id"]
            self.seq = record
            break

    def parse_atgc(self, result):
        seq_str = result["seq"]
        print seq_str
        match_id = re.match(r">.*\n", seq_str)
        if match_id:
            print "matched"
            seq_str = seq_str[len(match_id.group()):]
        self.seq = SeqRecord(Seq(seq_str), id=result["seq_id"], description=result["seq_id"])

def parse(fasta, genome_json, log, seq_prefix, output_prefix):
    global my_output_fa_record,my_output_log_record
    seq_records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    log_data = pd.read_table(log, header=None, index_col=0)
    for chr in genome_json:
        id_list = genome_json[chr]
        for id in id_list:
            if len(str(id)) == 24:
                from_mongo = FromMongo(chr)
                from_mongo.parse_mongo(id)
                my_output_fa_record.append(from_mongo.seq)
                my_output_log_record = my_output_log_record.append(from_mongo.log, ignore_index=True)
            else:
                seq_id = "%s%s" % (seq_prefix, id)
                my_output_fa_record.append(seq_records[seq_id])
                if log_data.loc[seq_id, 2] in ['circle', 'Circle']:
                    is_circle = "True"
                elif log_data.loc[seq_id, 2] in ['Linear', 'linear']:
                    is_circle = "False"
                my_output_log_record = my_output_log_record.append({"loc": chr, "seq": seq_id, "from": seq_id, "is_circle": is_circle}, ignore_index=True)
    SeqIO.write(my_output_fa_record, output_prefix + ".fa", "fasta")
    my_output_log_record.reindex(columns=["loc", "seq", "from", "is_circle"]).to_csv(output_prefix + ".log", sep="\t", header=False, index=False)

def parse2(fasta, genome_json, seq_prefix, output_prefix):
    global my_output_fa_record,my_output_log_record
    seq_records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    for chr in genome_json:
        id_list = genome_json[chr]
        for id in id_list:
            if len(str(id)) == 24:
                from_mongo = FromMongo(chr)
                from_mongo.parse_mongo(id)
                my_output_fa_record.append(from_mongo.seq)
                my_output_log_record = my_output_log_record.append(from_mongo.log, ignore_index=True)
            else:
                seq_id = "%s%s" % (seq_prefix, id)
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": chr, "seq": seq_id, "from": seq_id, "is_circle": "False"}, ignore_index=True)
    SeqIO.write(my_output_fa_record, output_prefix + ".fa", "fasta")
    my_output_log_record.reindex(columns=["loc", "seq", "from", "is_circle"]).to_csv(output_prefix + ".log", sep="\t", header=False, index=False)

def _main():
    parser = argparse.ArgumentParser(description='filter fa by json file')
    parser.add_argument('-fa', '--fa', help="fasta file")
    parser.add_argument('-json', '--json', help="json string")
    parser.add_argument('-log', '--log', help="log file")
    parser.add_argument('-seq_prefix', '--seq_prefix', help="seq name prefix in fasta file")
    parser.add_argument('-o_prefix', '--o_prefix', help="output prefix, output a new fasta file and a log file")
    parser.add_argument('-mongo', help="mongo verison", type=int, default=None)
    args = parser.parse_args()
    if args.mongo:
        Config().DBVersion = args.mongo
    if args.log:
        parse(args.fa, eval(args.json),args.log, args.seq_prefix, args.o_prefix)
    else:
        parse2(args.fa, eval(args.json), args.seq_prefix, args.o_prefix)

if __name__ == "__main__":
    _main()
    '''
        python get_fa_by_json.py -fa ~/app/bioinfo/bin/download/test.test.fa -json "{\"chr1\": [1,2,3,5,9,55], \"plasmidB\": [4,6,8]}" -seq_prefix Scaffold -o_prefix test
    '''