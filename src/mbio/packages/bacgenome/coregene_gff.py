# -*- coding: utf-8 -*-
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', metavar='[balst table]', required=True, help='blast table')
parser.add_argument("-g", metavar="[gff file]", required=False, help="gff file")
parser.add_argument('-o', metavar='[output]', required=True, help='out file')
args = parser.parse_args()

blastfile = args.b
gfffile = args.g
out = args.o

def run():
    genes ={}
    with open(blastfile,"r") as f:
        lines =f.readlines()
        for line in lines:
            lin = line.strip().split("\t")
            covage = (abs(int(lin[7]) - int(lin[6]))+1)/float(lin[3])*100
            genes[lin[0]] = "\t".join([lin[1], lin[2], str(covage), lin[10], lin[11]])
    with open(gfffile, "r") as k,open(out, "w") as g:
        g.write("Gene id\tGene Name\tLocation\tStrand\tStart\tEnd\tIdentity\tCoverage\tEvalue\tScore\n")
        lines = k.readlines()
        for line in lines[1:]:
            lin = line.strip().split("\t")
            if lin[0] in genes.keys():
                location = lin[1].split("_ORF")[0]
                dd =genes[lin[0]].split("\t")
                g.write("\t".join([lin[0], dd[0], location, lin[4], lin[2], lin[3], dd[1], dd[2], dd[3], dd[4]])+"\n")

if __name__ == '__main__':
    run()