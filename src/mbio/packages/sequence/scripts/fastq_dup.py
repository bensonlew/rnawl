# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
import argparse
import sys

usage = "USAGE\nSE: fastq_dup.py -s se.fastq - o out.fastq\nPE: fastq_dup.py -l fq_l.fastq -r fq_r.fastq - o out.fastq"
parser = argparse.ArgumentParser(description="calculation of fastq Duplicate sequences")
parser.add_argument("-s", "--fastq_s", type=str, help="se fastq")
parser.add_argument("-l", "--fastq_l", type=str, help="left fastq(PE)")
parser.add_argument("-r", "--fastq_r", type=str, help="right fastq(PE)")
# parser.add_argument("-t", "--type", type=str, choices=["PE", "SE"], help="fastq type, PE OR SE")
parser.add_argument("-o", "--output", type=str, help="output")
args = parser.parse_args()

# if args.fastq_s: 

if not args.fastq_l and not args.fastq_s:
    print usage
    # print("Error:must give the left fastq when the type is PE!")
    sys.exit(1)

if args.fastq_r and not args.fastq_l:
    print usage
    print("Error:must give the left fastq when the type is PE!")
    sys.exit(1)

if args.fastq_l and not args.fastq_r:
    print usage
    print("Error:must give the left fastq when the type is PE!")
    sys.exit(1)

if not args.output:
    print("Error:must give the out file name")
    sys.exit(1)

# if not args.type:
#     print("Error:must give the fastq type")
#     sys.exit(1)

fq_type = ""
if args.fastq_s:
    fq_type = "SE"
elif args.fastq_l and args.fastq_r:
    fq_type = "PE"


def fastq_dup(fastq):
    with open(fastq, "r") as f:
        line_count = {}
        base_line = 1
        dup = 0
        total_base = 0
        lines = set()
        for n, line in enumerate(f):
           if n % 4 == 1:
            total_base += 1
            lines.add(line.strip())
        un_dup = len(lines)
        dup = total_base - un_dup
    return dup/total_base
        
def fastq_pe_dup(fq_l, fq_r):
    with open(fq_l, "r") as l, open(fq_r, "r") as r:
        total_base = 0
        lines_lr = set()
        line_l = set()
        line_r = set()
        fql = {}
        for n, line in enumerate(l):
           if n % 4 == 1:
               total_base += 1
               line = line.strip()
               line_l.add(line)
               fql[n] = line
        for n, line in enumerate(r):
           if n % 4 == 1:
               line = line.strip()
               line_r.add(line)
               line_lr = fql[n] + line
               lines_lr.add(line_lr)
        dup_l = total_base - len(line_l)
        dup_r = total_base - len(line_r)
        un_dup = len(lines_lr)
        dup = total_base - un_dup
    return dup_l/total_base, dup_r/total_base, dup/total_base

def pe_dup(fq_l, fq_r):
    with open(fq_l, "r") as l,open(fq_r, "r") as r:
        lines_set = set()
        set_lineL = set()
        set_lineR = set()
        total_base = 0
        for n, (leftline, rightline) in enumerate(zip(l, r)):
            if n % 4 == 1:
                total_base += 1
                line_l = leftline.strip()
                set_lineL.add(line_l)
                line_r = rightline.strip()
                set_lineR.add(line_r)
                line_lr = line_l + line_r
                lines_set.add(line_lr)
        dup_lr = total_base - len(lines_set)
        dup_l = total_base - len(set_lineL)
        dup_r = total_base - len(set_lineR)
        # print(dup_l/total_base, dup_r/total_base, dup_lr/total_base)
        return  dup_l/total_base, dup_r/total_base, dup_lr/total_base

    # print dup

def write_out(dup_out):
    with open(dup_out, "w") as w:
        if fq_type == "SE":
            w.write("dup%\n")
            dup_rate = fastq_dup(args.fastq_s)
            sample = args.fastq_s.split("/")[-1]
            w.write("{}\n".format(dup_rate))
        else:
            w.write("read1Dup%\tread2Dup%\tPairedDup%\n")
            dupRate = fastq_pe_dup(args.fastq_l, args.fastq_r)
            # dupRate = pe_dup(args.fastq_l, args.fastq_r)
            w.write("{}\t{}\t{}\n".format(dupRate[0], dupRate[1], dupRate[2]))


if __name__ == "__main__":
    write_out(args.output)

# pe_dup("part_r1.fq", "part_r2.fq")
