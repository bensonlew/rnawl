# -*- coding: utf-8 -*-
# author: "shicaiping"
# __date__: 20200306

from Bio import SeqIO
import logging
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import os, re, math, time
from biocluster.config import Config
import subprocess


class SeqKit(object):
    def __init__(self, args):
        self.args = args
        if 'output_dir' in args and args.output_dir:
            if not os.path.exists(args.output_dir):
                os.makedirs(args.output_dir)
            dir = args.output_dir
        else:
            dir = os.getcwd()
        if 'output_name' in args and args.output_name:
            name = args.output_name
        else:
            if args.subparsers == "length":
                name = os.path.basename(args.input_file) + ".length.txt"
            elif args.subparsers == "translate":
                name = os.path.basename(args.input_file) + ".cds2pep.faa"
            elif args.subparsers == "filter":
                name = os.path.basename(args.input_file) + ".filtered.fa"
            else:
                name = ""
        self.output = os.path.join(dir, name)
        self.func = {
            'length': self.length,
            'translate': self.translate,
            'extract': self.extract,
            'filter': self.filter,
        }[args.subparsers]

    def run(self):
        logging.info('begin of the function ({})'.format(self.func.__name__))
        self.func(self.args)
        logging.info('final of the function ({})'.format(self.func.__name__))

    def length(self, args):
        with open(self.output, "w") as w:
            for seq_record in SeqIO.parse(args.input_file, args.file_format):
                seq_id = seq_record.id
                seq_length = len(seq_record.seq)
                w.write("{}\t{}\n".format(seq_id, seq_length))

    def translate(self, args):
        output_handle = open(self.output, "w")
        for rec in SeqIO.parse(args.input_file, "fasta"):
            rec.seq=rec.seq.translate(table=args.transl_table)
            SeqIO.write(rec, output_handle, "fasta")
        output_handle.close()

    def extract(self, args):
        self.samtools_path = Config().SOFTWARE_DIR + '/bioinfo/rna/miniconda2/bin/'
        header = os.path.join(os.getcwd(), "header.txt")
        cmd = "{}samtools view -H {} > {}".format(self.samtools_path, args.input_file, header)
        logging.info(cmd)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            raise Exception("运行{}出错！".format(cmd))
        with open(header, "r") as f:
            for line in f:
                if line.startswith("@HD"):
                    # if 'SO:unsorted' in line or "SO:queryname" in line:
                    if 'SO:coordinate' in line:
                        input_file_sort = args.input_file
                    else:
                        input_file_sort = os.path.join(os.getcwd(), os.path.basename(args.input_file) + ".sorted")
                        cmd = "{}samtools sort -n {} -o {}".format(self.samtools_path, args.input_file, input_file_sort)
                        logging.info(cmd)
                        try:
                            subprocess.check_call(cmd, shell=True)
                        except subprocess.CalledProcessError:
                            raise Exception("运行{}出错！".format(cmd))
        if 'input_file_sort' not in locals().keys():
            input_file_sort = args.input_file
        if args.type == "mapped":
            if args.seq_type == "PE":
                cmd = "{}samtools fastq -N -f 2 -1 {} -2 {} -t {}".format(self.samtools_path, args.read1, args.read2, input_file_sort)
            else:
                cmd = "{}samtools fastq -f 2 -s {} -t {}".format(self.samtools_path, args.single, input_file_sort)
            logging.info(cmd)
            try:
                subprocess.check_call(cmd, shell=True)
            except subprocess.CalledProcessError:
                raise Exception("运行{}出错！".format(cmd))
        elif args.type == "unmapped":
            if args.seq_type == "PE":
                cmd = "{}samtools fastq -N -f 12 -1 {} -2 {} -t {}".format(self.samtools_path, args.read1, args.read2, input_file_sort)
            else:
                cmd = "{}samtools fastq -f 12 -s {} -t {}".format(self.samtools_path, args.single, input_file_sort)
            logging.info(cmd)
            try:
                subprocess.check_call(cmd, shell=True)
            except subprocess.CalledProcessError:
                raise Exception("运行{}出错！".format(cmd))
        elif args.type == "duplicate":
            if args.seq_type == "PE":
                cmd = "{}samtools fastq -f 1024 -1 {} -2 {} -t {}".format(self.samtools_path, args.read1, args.read2, input_file_sort)
            else:
                cmd = "{}samtools fastq -f 1024 -s {} -t {}".format(self.samtools_path, args.single, input_file_sort)
            logging.info(cmd)
            try:
                subprocess.check_call(cmd, shell=True)
            except subprocess.CalledProcessError:
                raise Exception("运行{}出错！".format(cmd))


    def filter(self, args):
        pass


if __name__ == "__main__":
    import argparse
    import os
    import sys
    parser = argparse.ArgumentParser(description='Tool kit for sequence file processing.')
    subparsers = parser.add_subparsers(dest='subparsers', help='sub-command help')

    parser_l = subparsers.add_parser('length', help='This script is used to calculate the sequence length.')
    parser_l.add_argument('-i', '--input_file', type=str, required=True, help='input file path')
    parser_l.add_argument('-f', '--file_format', type=str, required=True, choices=['fasta', 'fastq'], help='file type')
    parser_l.add_argument('-o', '--output_dir', type=str, help='output directory, default is current dir.')
    parser_l.add_argument('-n', '--output_name', type=str, help='output file name, default is add suffix ".length.txt" to the input file name')

    parser_t = subparsers.add_parser('translate', help='This script is used to translate cds to pep.')
    parser_t.add_argument('-i', '--input_file', type=str, required=True, help='input file path')
    parser_t.add_argument('-t', '--transl_table', type=int, default=1, help='Genetic codes, 1 for standard code, 2 for vertebrate mitochondrial code, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi to get more details.')
    parser_t.add_argument('-o', '--output_dir', type=str, help='output directory, default is current dir.')
    parser_t.add_argument('-n', '--output_name', type=str, help='output file name, default is add suffix ".pep.faa" to the input file name')

    parser_e = subparsers.add_parser('extract', help='This script is used to extracts seqs from input bam/sam file.')
    parser_e.add_argument('-i', '--input_file', type=str, required=True, help='input file path')
    #parser_l.add_argument('-f', '--file_format', type=str, required=True, choices=['bam', 'sam'], help='input file format')
    parser_e.add_argument('-t', '--type', type=str, required=True, choices=['mapped', 'unmapped', 'duplicate'], help='mapped: read mapped in proper pair; unmapped: reads unmapped; duplicate: read is PCR or optical duplicate')
    parser_e.add_argument('-q', '--seq_type', type=str, required=True, choices=['SE', 'PE'], help='SE: single-end; PE: paired-end')
    parser_e.add_argument('-o', '--output_dir', type=str, help='output directory, default is current dir.')
    parser_e.add_argument('-s', '--single', type=str, help='write singleton reads to FILE', default="output.fq")
    parser_e.add_argument('-r1', '--read1', type=str, help='write paired reads flagged READ1 to FILE', default="output_r1.fq")
    parser_e.add_argument('-r2', '--read2', type=str, help='write paired reads flagged READ2 to FILE', default="output_r2.fq")

    parser_f = subparsers.add_parser('filter', help='This script is used to filter sequence.')
    parser_f.add_argument('-i', '--input_file', type=str, required=True, help='input file path')
    parser_f.add_argument('-o', '--output_dir', type=str, help='output directory, default is current dir.')
    parser_f.add_argument('-n', '--output_name', type=str, help='output file name, default is add suffix ".output" to the input file name')

    logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

    if len(sys.argv) > 1 and sys.argv[1] not in ['-h', '--help']:
        if sys.argv[1] in ['length', 'translate', 'extract', 'filter']:
            args = parser.parse_args()
        else:
            logging.error('unrecognized command ({})'.format(sys.argv[1]))
            sys.exit(-2)
    else:
        parser.print_help()
        sys.exit(-1)

    inst = SeqKit(args)
    inst.run()