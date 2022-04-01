# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
from optparse import OptionParser
import types
from Bio import SeqIO
from bson.objectid import ObjectId
from biocluster.api.database.base import Base
from biocluster.config import Config

class Process(Base):
    def __init__(self, genome_id):
        super(Process, self).__init__()
        self._project_type = "bac_assem"
        self.genome_id = genome_id

    def check_id(self, main_id):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
        return main_id

    def get_genome_info(self):
        self.genome_id = self.check_id(self.genome_id)
        genome_info = self.db["genome"].find_one({"_id": self.genome_id})
        return genome_info["gap_detail_ids"]

    def get_seq_info(self):
        genome_info = self.get_genome_info()
        print ">>>genome_info"
        print genome_info
        seq_info = {}
        chr_count_info = {}
        for comp_id in genome_info:
            comp_id = self.check_id(comp_id)
            seq_detail = self.db["gap_fill_detail"].find_one({"_id": comp_id})
            scaf_name = seq_detail["scaf_name"]
            location = seq_detail["location"]
            if location not in chr_count_info.keys():
                chr_count_info[location] = 1
            else:
                chr_count_info[location] += 1
                location = location + ".%s" % chr_count_info[location]
            seq_info[scaf_name] = location
        return seq_info,chr_count_info

    def get_seq_file(self, fa, out):
        seq_list = []
        seq_info,chr_count_info = self.get_seq_info()
        print ">>>seq_info"
        print seq_info
        seq_set = set(seq_info)
        for seq_record in SeqIO.parse(fa, "fasta"):
            if seq_record.id in seq_set:
                location = seq_info[seq_record.id]
                yuanshi = location.split(".")[0]
                if chr_count_info[yuanshi] == 1:
                    seq_record.id = seq_info[seq_record.id]
                elif "." in location:
                    seq_record.id = location
                else:
                    seq_record.id = location + ".1"
                seq_list.append(seq_record)
        SeqIO.write(seq_list, out, "fasta")

def main():
    parser = OptionParser()
    parser.add_option('--f', dest ='fa',metavar='[fasta file]')
    parser.add_option('--id', dest ='id',metavar='[genome id]')
    parser.add_option("--o", dest="output", metavar="[out fa file]")
    parser.add_option('--m', dest="mongo", help="mongo verison", type=int, default=None)
    (options, args) = parser.parse_args()
    if options.mongo:
        Config().DBVersion = options.mongo
    if not options.fa or not options.id or not options.output:
        print "python complete_get_seq.py --f fasta --id genome_id --o out.fa"
        return
    Process(options.id).get_seq_file(options.fa, options.output)

if __name__=='__main__':
    main()