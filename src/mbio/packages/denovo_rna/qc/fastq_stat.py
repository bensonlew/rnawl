# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
import math


def fastq_qual_stat(qual_stat):
    with open(qual_stat, 'rb') as f, open('{}.base'.format(qual_stat), 'wb') as b, open('{}.err'.format(qual_stat), 'wb') as e, \
            open('{}.qual'.format(qual_stat), 'wb') as q:
        b.write('position\tA\tT\tC\tG\tN\n')
        e.write('position\terror\n')
        q.write('position\tLW\tQ1\tmed\tQ3\tRW\n')
        f.readline()
        for line in f:
            ln = line.strip().split('\t')
            A_base = int(ln[12])/int(ln[17]) * 100
            T_base = int(ln[15])/int(ln[17]) * 100
            C_base = int(ln[13])/int(ln[17]) * 100
            G_base = int(ln[14])/int(ln[17]) * 100
            N_base = int(ln[16])/int(ln[17]) * 100
            err = 10 ** (float(ln[5])/(-10)) * 100
            b_line = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(ln[0], A_base, T_base, C_base, G_base, N_base)
            e_line = '{}\t{}\n'.format(ln[0], err)
            q_line = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(ln[0], ln[10], ln[6], ln[7], ln[8], ln[11])
            b.write(b_line)
            q.write(q_line)
            e.write(e_line)


def get_fastq_info(fastq):
    with open(fastq, "r") as f:
        total_reads = 0
        # total_bases = 0
        total_Reads_with_Ns = 0
        bases = ""
        qualities = []
        more_than_20 = 0
        more_than_30 = 0
        base_line = 1
        quality_line = 3
        for n, line in enumerate(f):
            if n == base_line:
                total_reads += 1
                bases += line.strip()
                if "N" in line:
                    total_Reads_with_Ns += 1
                base_line += 4
            elif n == quality_line:
                for i in line.strip():
                    quality = ord(i) - 33
                    qualities.append(quality)
                    if quality > 20:
                        more_than_20 += 1
                    if quality > 30:
                        more_than_30 += 1
                quality_line += 4
        total_bases = len(bases)
        reads_with_N = total_Reads_with_Ns/total_reads
        A_Rate = bases.count("A")/total_bases
        T_Rate = bases.count("T")/total_bases
        C_Rate = bases.count("C")/total_bases
        G_Rate = bases.count("G")/total_bases
        N_Rate = bases.count("N")/total_bases
        Q20_Rate = more_than_20/total_bases
        Q30_Rate = more_than_30/total_bases
        error_Rate = math.pow(10.0, (sum(qualities)/total_bases) * -0.1) * 100.0
        return {"total_reads": total_reads, "total_bases": total_bases, "Total_Reads_with_Ns": total_Reads_with_Ns,
                "N_Reads%": reads_with_N, "A%": A_Rate, "T%": T_Rate, "C%": C_Rate, "G%": G_Rate, "N%": N_Rate,
                "Error%": error_Rate, "Q20%": Q20_Rate, "Q30%": Q30_Rate, "CG%": C_Rate+G_Rate}
