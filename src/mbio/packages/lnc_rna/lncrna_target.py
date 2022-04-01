# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
import subprocess
import shlex
import argparse
import re
from multiprocessing.dummy import Pool
from Bio import SeqIO

def run_cmd(cmd):
    os.system(cmd)

class TargetPrediction():
    def __init__(self, query_fa, target_fa,seq_num=100, method='rnaplex',
                 output=None, rnaplex=None, riblast=None,
                 intarna=None, lnctar=None, parafly=None,  pool=10, other_param = ""):
        self.method = method
        self.seq_num = seq_num
        self.query = query_fa
        self.target = target_fa
        self.split_fasta = []
        self.rnaplex = rnaplex
        self.riblast = riblast
        self.intarna = intarna
        self.lnctar = lnctar
        self.pool = pool
        self.output = output
        self.other_param = other_param
        self.max_e = -20
        self.max_e_riblast = -10
        self.parafly = parafly

    def splitfasta(self):
        '''
        按子线程数拆分文件
        '''
        split_fasta = self.query
        if self.method == "riblast":
            split_fasta = self.target

        file_split = dict()

        for i in range(self.pool):
            split_file = os.path.basename(split_fasta) + '_' + str(i)
            self.split_fasta.append(split_file)
            file_split[i] = open(split_file, 'w')

        num = 0
        for seq in SeqIO.parse(split_fasta, "fasta"):
            file_num = num % self.pool
            num += 1
            file_split[file_num].write('>{}\n{}\n'.format(seq.name, seq.seq))
        for i in range(self.pool):
            file_split[i].close()

    def rnaplex_list(self):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} '.format(self.rnaplex)
            cmd += '-f 2 '
            cmd += '-q {} '.format(split)
            cmd += '-t {} '.format(self.target)
            cmd += ' {} '.format(self.other_param)
            cmd += '>{} '.format(split + self.output)
            cmd_list.append(cmd)
        return cmd_list

    def riblast_list(self, db_file):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} ris '.format(self.riblast)
            cmd += '-i {} '.format(split)
            cmd += '-d {} '.format(db_file)
            cmd += ' {} '.format(self.other_param)
            cmd += '-o {} '.format(split + self.output)
            cmd_list.append(cmd)
        return cmd_list

    def riblast_db(self):
        db = os.path.basename(self.query) + "_db"
        cmd = '{} db '.format(self.riblast)
        cmd += '-i {} '.format(self.query)
        cmd += '-o {} '.format(db)
        os.system(cmd)
        return db

    def intarna2_list(self):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} --outMode C --threads 4 '.format(self.intarna)
            cmd += '-q {} '.format(split)
            cmd += '-t {} '.format(self.target)
            cmd += ' {} '.format(self.other_param)
            cmd += '> {} '.format(split + self.output)
            cmd_list.append(cmd)
        return cmd_list

    def intarna_list(self):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} -o '.format(self.intarna)
            cmd += '-t {} '.format(split)
            cmd += '-m {} '.format(self.target)
            cmd += ' {} '.format(self.other_param)
            cmd += '> {} '.format(split + self.output)
            cmd_list.append(cmd)
        return cmd_list

    def lnctar_list(self):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} -p 1 -d -0.1 -s F '.format(self.lnctar)
            cmd += '-l {} '.format(split)
            cmd += '-m {} '.format(self.target)
            cmd += ' {} '.format(self.other_param)
            cmd += '-o {} '.format(split + self.output)
            cmd_list.append(cmd)
        return cmd_list

    def para_run(self, cmd_list):
        with open("para_cmd.sh", 'w') as para_f:
            for cmd in cmd_list:
                para_f.write(cmd + "\n")
        para_cmd = "{} -c {} -CPU {}".format(self.parafly, "para_cmd.sh", self.pool)
        print para_cmd
        os.system(para_cmd)

    def run_predict(self):
        if self.method.lower() == 'rnaplex':
            cmd_list = self.rnaplex_list()
        elif self.method.lower() == 'riblast':
            db_file = self.riblast_db()
            cmd_list = self.riblast_list(db_file)
        elif self.method.lower() == 'intarna':
            cmd_list = self.intarna2_list()
        elif self.method.lower() == 'lnctar':
            cmd_list = self.lnctar_list()
        else:
            raise Exception(self.method + ' is not supported')
        '''
        if len(self.split_fasta) > self.pool:
            pool = Pool(self.pool)
        else:
            pool = Pool(len(self.split_fasta))
        print(cmd_list)
        pool.map(run_cmd, cmd_list)
        pool.close()
        pool.join()
        '''
        self.para_run(cmd_list)
        self.merge_generate()

    def merge_generate(self):
        if self.method.lower() == 'rnaplex':
            result_list = []
            last_list = []
            with open('rnaplex_merge_out', 'w') as f:
                f.write("lncRNA_id\ttarget_mRNA_id\tlncRNA_start\tlncRNA_end\tmRNA_start\tmRNA_end\tEnergy\n")
                #Seq\tTarget\tAlignment\tEnergy\tQ_start\tQ_end\tT_start\tT_end\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                for single in result_list:
                    split_f = open(single,'r')
                    while split_f:
                        target = split_f.readline().strip().lstrip(">")
                        if target == "":
                            break
                        lnc_rna =  split_f.readline().strip().lstrip(">")
                        align = split_f.readline().strip()
                        eles = align.split()
                        try:
                            if float(eles[4].rstrip(")").lstrip("(")) < self.max_e:
                                f.write("\t".join([
                                    lnc_rna, target, eles[3].split(",")[0], eles[3].split(",")[1], eles[1].split(",")[0], eles[1].split(",")[1], eles[4].rstrip(")").lstrip("(")
                                ]) + "\n")
                        except:
                            pass
                    split_f.close()

        elif self.method.lower() == 'riblast':
            result_list = []
            last_list = []
            with open('riblast_merge_out', 'w') as f:
                f.write("lncRNA_id\ttarget_mRNA_id\tlncRNA_start\tlncRNA_end\tmRNA_start\tmRNA_end\tEnergy\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                for single in result_list:
                    if not os.path.exists(single):
                        continue
                    with open(single, 'r') as f_in:
                        pair = ""
                        low_energy = ""
                        last_energy = 0
                        f_in.readline()
                        f_in.readline()
                        f_in.readline()
                        for line in f_in:
                            cols = line.strip().split(",")
                            pair_new = cols[1] + cols[3]
                            if pair_new != pair and low_energy != "":
                                cols = low_energy.strip().split(",")
                                # print "**" + cols[6]
                                pos = re.split(r"[:-]", cols[6].rstrip(")").lstrip("("))
                                if float(cols[5]) < self.max_e_riblast:
                                    f.write("\t".join([
                                        cols[3], cols[1]
                                    ] + [pos[2], pos[3], pos[0], pos[1]] + [cols[5]]) + "\n")
                                last_energy = 0
                            else:
                                if float(cols[5]) < last_energy:
                                    last_energy = cols[5]
                                    low_energy = line
                            pair = pair_new
        elif self.method.lower() == 'intarna':
            result_list = []
            with open('intarna_merge_out', 'w') as f:
                f.write("lncRNA_id\ttarget_mRNA_id\tlncRNA_start\tlncRNA_end\tmRNA_start\tmRNA_end\tEnergy\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                for single in result_list:
                    if not os.path.exists(single):
                        continue
                    with open(single, 'r') as f_in:
                        f_in.readline()
                        for line in f_in:
                            cols = line.strip().split(";")
                            if float(cols[-1]) < self.max_e:
                                f.write("\t".join([
                                    cols[3],
                                    cols[0],
                                    cols[4],
                                    cols[5],
                                    cols[1],
                                    cols[2],
                                    cols[-1]
                                ]) + "\n")
        elif self.method.lower() == 'intarna1':
            result_list = []
            last_list = []
            with open('intarna_merge_out', 'w') as f:
                f.write("lncRNA_id\ttarget_mRNA_id\tlncRNA_start\tlncRNA_end\tmRNA_start\tmRNA_end\tEnergy\n")
                # Seq\tTarget\tEnergy\tQ_start\tQ_end\tT_start\tT_end\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                target_dict = dict()
                for single in result_list:
                    lnc_rna = ""
                    target = ""
                    start = True
                    split_f = open(single,'r')
                    for line in split_f:
                        if line.startswith("======"):
                            start = True
                        if line.startswith(">"):
                            if start:
                                lnc_rna = line.strip().lstrip(">")
                                start = False
                            else:
                                target = line.strip().lstrip(">")
                        if line.startswith("positions(target)"):
                            poss = line.strip().split()
                            target_dict[lnc_rna + "||" + target] = {"target_start" : poss[2], "target_end" : poss[4]}
                        if line.startswith("positions(ncRNA)"):
                            poss = line.strip().split()
                            target_dict[lnc_rna + "||" + target]["lnc_start"] = poss[2]
                            target_dict[lnc_rna + "||" + target]["lnc_end"] = poss[4]

                        if line.startswith("energy:"):
                            poss = line.strip().split()
                            target_dict[lnc_rna + "||" + target]["energy"] = poss[1]

                for key,value in target_dict.items():
                    f.write("\t".join([
                        key.split("||")[0],
                        key.split("||")[1],
                        value["lnc_start"],
                        value["lnc_end"],
                        value["target_start"],
                        value["target_end"],
                        value["energy"],
                    ]) + "\n")

        elif self.method.lower() == 'lnctar':
            result_list = []
            with open('lnctar_merge_out', 'w') as f:
                f.write("lncRNA_id\ttarget_mRNA_id\tlncRNA_start\tlncRNA_end\tmRNA_start\tmRNA_end\tEnergy\n")
                #Target\tQuery\tdG\tndG\tQ_start\tQ_end\tT_start\tT_end\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                for single in result_list:
                    split_f = open(single,'r')
                    header = split_f.readline()
                    for line in split_f:
                        cols = line.strip().split("\t")
                        if float(cols[4]) < 0:
                            f.write("\t".join([
                                cols[0], cols[2],  cols[6], cols[7], cols[8], cols[9], cols[4]
                            ]) + "\n")
                    split_f.close()
        else:
            pass



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', type=str, metavar="query", required=True,
                        help="miRNA fasta file")
    parser.add_argument('-parafly', type=str, metavar="parafly", required=True,
                        help="parafly dir")
    parser.add_argument('-t', type=str, metavar="target", required=True,
                        help="target fasta file")
    parser.add_argument('-p', type=int, metavar="process", default=10,
                        help="max num of processes")
    parser.add_argument('-m', type=str, metavar="method", default="miranda",
                        help="miranda or PITA or RNAhybrid or psRobot or targetfinder")
    parser.add_argument('-n', type=str, metavar="split_num",default=100,
                        help="split query num")
    parser.add_argument('-o', type=str, metavar="output", default="target_predict",
                        help="prefix to the output")
    parser.add_argument('-rnaplex', type=str, metavar="rnaplex_path",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/lnc_rna/RNAplex")
    parser.add_argument('-riblast', type=str, metavar="riblast_path",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/lnc_rna/RIblast-1.1.3/RIblast")
    parser.add_argument('-intarna', type=str, metavar="intarna_path",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/lnc_rna/IntaRNA/bin/IntaRNA")
    parser.add_argument('-lnctar', type=str, metavar="lnctar",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/lnc_rna/LncTar/LncTar.pl")

    args = parser.parse_args()
    toolbox=TargetPrediction(
        query_fa=args.q,
        target_fa=args.t,
        seq_num=args.n,
        method=args.m,
        pool=args.p,
        output=args.o,
        rnaplex=args.rnaplex,
        intarna=args.intarna,
        riblast=args.riblast,
        lnctar=args.lnctar,
        parafly=args.parafly
    )
    toolbox.splitfasta()
    toolbox.run_predict()
