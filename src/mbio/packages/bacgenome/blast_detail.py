# -*- coding: utf-8 -*-
import pandas as pd
from biocluster.config import Config
import os
import argparse

def blast_detail(sample, blast_table, fasta_file, outdir):
    location_dict = add_location(fasta_file)
    outfile = os.path.join(outDir, sample + ".xls")
    with open(blast_table, 'r') as infile, open(outfile, 'wb') as out:
        out.write(
            'Sample Name\tgene ID\tLocation\tquery start\tquery end\tsubject start\tsubject end\talign length\tIdentity (%)\tCoverage(%)\tEvalue\tScore\n')
        head = infile.next()
        if "Score" in head:
            for line in infile:
                line = line.strip().split("\t")
                #gene_id = line[5]
                gene_id = line[10]
                location = location_dict[gene_id]
                qstat = line[7]
                qend = line[8]
                sstat = line[12]
                send = line[13]
                align_len = abs(float(line[8]) - float(line[7]))
                iden = line[3]
                evalue = line[1]
                score = line[0]
                qlen = line[6]
                coverge = round(float(align_len) / int(qlen), 3)
                coverge = coverge * 100
                out.write(
                    '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sample, gene_id, location, qstat, qend,
                                                                              sstat, send, align_len, iden,
                                                                              coverge, evalue, score))
        else:
            print "file error"

def add_location(fasta):
    location = {}
    with open(fasta, 'rb') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                line = line.strip().split(' ')
                gene_name = line[0][1:]
                location_info = line[3].split('_')
                location[gene_name] = location_info[0]
    return location

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[blast_table]', required=True, help='Input blast table')
    parser.add_argument('-fa', metavar='[faa_file]', required=True, help='Input faa_file')
    parser.add_argument('-s', metavar='[sample]', required=True, help='input sample name')
    parser.add_argument('-o', metavar='[outDir]', required=True, help='output Dir')
    args = parser.parse_args()
    blast_table = args.i
    sample = args.s
    faafile = args.fa
    outDir = args.o
    blast_detail(sample, blast_table, faafile, outDir)
