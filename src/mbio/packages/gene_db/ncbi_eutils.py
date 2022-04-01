# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from BCBio import GFF
import urllib
import sys
import argparse
from eutils import Client


class NcbiEutils(object):
    def __init__(self, genome_acc):

        self.ec = Client(api_key=os.environ.get("NCBI_API_KEY", None))

    def pipe(self, gff, out_pre):
        genome_dict = self.get_genome_acc(self.genome_acc)
        # print genome_dict
        self.gff2table(gff, out_pre, genome_dict)

    def get_genome_acc(self, genome_acc):


    def gff2table(self, in_file, out_pre, genome_dict):



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-acc', type=str, required=True,
                        help='genome acc')
    parser.add_argument('-gff', type=str, required=True,
                        help='refseq gff')
    parser.add_argument('-out_pre', type=str, required=False, default="table",
                        help="output pre")

    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    mg = NcbiGff(args.acc)
    mg.pipe(args.gff, args.out_pre)
