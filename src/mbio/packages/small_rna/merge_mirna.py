# -*- coding: utf-8 -*-
# __author__ : 'shicaiping'
# __date__: 20200715

from Bio import SeqIO
import os
from biocluster.config import Config


class Merge(object):
    """
   This script is used to merge and remove redundancy of mirna mature and hairpin seqs.
    """

    def __init__(self, priority_list, mature2hairpin, database='mirbase', prefix='all'):
        if not os.path.exists(priority_list):
            raise Exception("{} not exist!".format(priority_list))
        if not os.path.exists(mature2hairpin):
            raise Exception("{} not exist!".format(mature2hairpin))
        if database.lower() not in ['mirbase', 'pmiren']:
            raise Exception("{} not supported!".format(database))
        self.mature2hairpin = mature2hairpin
        self.priority_list = priority_list
        self.prefix = prefix
        self.database = database
        self.file_dir = Config().SOFTWARE_DIR + '/database'
        self.output_mature = prefix + "_mature.fa"
        self.output_hairpin = prefix + "_hairpin.fa"

    def run(self):
        mature2hairpin = self.get_mature2hairpin()
        seqs = set()
        merge_mature = dict()
        merge_hairpin = dict()
        with open(self.priority_list, "r") as f:
            for line in f:
                organism = line.strip().split()[0]
                if self.database.lower() == "mirbase":
                    mature_fa = self.file_dir + "/mirbase/species/" + organism + "/" + organism + ".mature.fa"
                    hairpin_fa = self.file_dir + "/mirbase/species/" + organism + "/" + organism + ".hairpin.fa"
                else:
                    mature_fa = self.file_dir + "/PmiREN/species/" + organism + "/" + organism + "_mature.fa"
                    hairpin_fa = self.file_dir + "/PmiREN/species/" + organism + "/" + organism + "_hairpin.fa"
                if not os.path.exists(mature_fa) or not os.path.exists(hairpin_fa):
                    continue
                else:
                    mature_fa_dict = self.parse_fa(mature_fa)
                    hairpin_fa_dict = self.parse_fa(hairpin_fa)
                    for mature in mature_fa_dict:
                        hairpin = mature2hairpin[mature]
                        seq = mature_fa_dict[mature] + hairpin_fa_dict[hairpin]
                        if seq not in seqs:
                            merge_mature[mature] = mature_fa_dict[mature]
                            merge_hairpin[hairpin] = hairpin_fa_dict[hairpin]
                            seqs.add(seq)
                        else:
                            continue
        with open(self.output_mature, "w") as w:
            for record in sorted(merge_mature):
                w.write(">" + record + "\n")
                w.write(str(merge_mature[record]) + "\n")
        with open(self.output_hairpin, "w") as w:
            for record in sorted(merge_hairpin):
                w.write(">" + record + "\n")
                w.write(str(merge_hairpin[record]) + "\n")

    def get_mature2hairpin(self):
        mature2hairpin = dict()
        with open(self.mature2hairpin, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if item[0] not in mature2hairpin:
                    mature2hairpin[item[0]] = item[1]
        return mature2hairpin

    def parse_fa(self, fa):
        fa_dict = dict()
        for seq_record in SeqIO.parse(fa, "fasta"):
            mirna_id = seq_record.id
            mirna_seq = seq_record.seq
            fa_dict[mirna_id] = mirna_seq
        return fa_dict


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='This script is used to merge and remove redundancy of mirna mature and hairpin seqs.')
    parser.add_argument('-priority_list', type=str, help='file of organism priority list')
    parser.add_argument('-mature2hairpin', type=str, help='file of organism priority list')
    parser.add_argument('-database', type=str, default="mirbase", help='can be mirbase or pmiren')
    parser.add_argument('-prefix', type=str,
                        help='prefix of output file, default will be all_mature.fa and all_hairpin.fa')
    args = parser.parse_args()

    merge = Merge(args.priority_list, args.mature2hairpin, args.database, args.prefix)
    merge.run()
