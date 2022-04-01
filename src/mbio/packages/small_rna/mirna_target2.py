# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:2018.11.30
import os
import subprocess
import shlex
import argparse
from multiprocessing.dummy import Pool
from Bio import SeqIO
import gzip

def run_cmd(cmd):
    subprocess.call(cmd)

class TagertPrediction():
    def __init__(self, query_fa, target_fa,seq_num=100, method='psRobot',
                 output=None, psRobot=None, targetfinder=None, spiece='animal',
                 miranda=None, pita=None, rnahybrid=None, targetscan=None, target_species=None, target_context=None, pool=10):
        self.method = method
        self.seq_num = seq_num
        self.query = query_fa
        self.target = target_fa
        self.split_fasta = []
        self.psRobot = psRobot
        self.targetfinder = targetfinder
        self.rnahybrid = rnahybrid
        self.miranda = miranda
        self.pita = pita
        self.pool = pool
        self.output = output
        self.spiece = spiece
        self.targetscan = targetscan
        self.target_species = target_species
        self.target_context = target_context

    def set_options(self, opt_dict):
        self.opts = opt_dict

    def splitfasta(self):
        with open(self.query) as f:
            seq_id = 0
            split_num = 0
            self.split_fasta.append(os.path.basename(f.name) + '_' + str(split_num))
            out_name = os.path.basename(f.name) + '_' + str(split_num)
            file_split = {split_num: open(out_name, 'w')}
            for line in f:
                if line.startswith('>'):
                    split_num = seq_id//self.seq_num
                    seq_id = seq_id + 1
                    if split_num not in file_split.keys():
                        self.split_fasta.append(os.path.basename(f.name) + '_' + str(split_num))
                        file_split[split_num-1].close()
                        out_name = os.path.basename(f.name) + '_' + str(split_num)
                        file_split[split_num] = open(out_name,'w')
                file_split[split_num].write(line)
            else:
                file_split[split_num].close()

    def miranda_list(self,sc="165",en="-25",strict="yes"):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} '.format(self.miranda)
            cmd += '{} '.format(split)
            cmd += '{} '.format(self.target)
            cmd += '-sc {} '.format(sc)
            cmd += '-en {} '.format(en)
            if strict == "yes":
                cmd += '-strict '                 #animal demand 5' seed region pairing strictly
            cmd += '-out {} '.format(split + self.output)
            cmd_list.append(shlex.split(cmd))
        return cmd_list

    def pita_list(self,gu="6;0,7;0,8;1",m="6;0,7;0,8;1"):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} '.format(self.pita)
            cmd += '-utr {} '.format(self.target)
            cmd += '-mir {} '.format(split)
            cmd += '-prefix {} '.format(split + self.output)
            cmd += '-gu {} '.format(gu)
            cmd += '-m {} '.format(m)
            cmd_list.append(shlex.split(cmd))
        return cmd_list

    def RNAhybrid_list(self, b=1, e=-20, pvalue=0.01,distribution="2,0.2"):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} '.format(self.rnahybrid)
            cmd += '-d {} '.format(distribution)
            cmd += '-t {} '.format(self.target)
            cmd += '-q {} '.format(split)
            cmd += '-b {} '.format(b)
            cmd += '-e {} '.format(e)
            cmd += '-m {} '.format(10000)
            cmd += '-p {} '.format(pvalue)
            cmd += '-f {} '.format("2,7")
            cmd += '> {}'.format(split + self.output)
            with open(split+'_hybrid.bash', 'w') as hy_h:
                hy_h.write(cmd)
            cmd_sh = 'bash ' + split+'_hybrid.bash'
            cmd_list.append(shlex.split(cmd_sh))
        return cmd_list

    # 旧版target_score 默认值为4
    def psRobot_list(self,target_score=2.5,):
        cmd_list = []
        if os.path.exists("target.fa"):
            os.remove("target.fa")
        os.link(self.target, "target.fa")
        for split in self.split_fasta:
            cmd = '{} '.format(self.psRobot)
            cmd += '-s {} '.format(split)
            cmd += '-t {} '.format("target.fa")
            cmd += '-o {} '.format(split+self.output)
            cmd += '-ts {}'.format(target_score)
            cmd_list.append(shlex.split(cmd))
        return cmd_list

    # 旧版score默认值为8
    def targetfinder_list(self,fmt='table',score=4):
        cmd_list = []
        for split in self.split_fasta:
            cmd = '{} '.format(self.targetfinder)
            cmd += '-f {} '.format(split)
            cmd += '-d {} '.format(self.target)
            cmd += '-o {} '.format(split + self.output)
            cmd += '-p {} '.format(fmt)
            cmd += '-c {}'.format(score)
            cmd_list.append(shlex.split(cmd))
        return cmd_list

    def run_predict(self):
        if self.method.lower() == 'miranda' and self.spiece.lower() == 'animal':
            cmd_list = self.miranda_list(sc=self.opts["miranda_score"],
                                         en=self.opts["miranda_energy"],
                                         strict=self.opts["miranda_strict"])
        elif self.method.lower() == 'pita' and self.spiece.lower() == 'animal':
            cmd_list = self.pita_list()
        elif self.method.lower() == 'rnahybrid':
            cmd_list = self.RNAhybrid_list(e=self.opts["rnahybird_energy"],
                                           pvalue=self.opts["rnahybird_pvalue"],
                                           b=self.opts["rnahybird_num"])
        elif self.method.lower() == 'psrobot' and self.spiece.lower() == 'plant':
            cmd_list = self.psRobot_list(target_score=self.opts["ps_robot_score"])
        elif self.method.lower() == 'targetfinder' and self.spiece.lower() == 'plant':
            cmd_list = self.targetfinder_list(score=self.opts["targetfinder_score"])
        elif self.method.lower() == 'targetscan':
            self.merge_generate()
            return
        else:
            raise Exception(self.method + ' is not supported')
        if len(self.split_fasta) > self.pool:
            pool = Pool(self.pool)
        else:
            pool = Pool(len(self.split_fasta))
        print(cmd_list)
        pool.map(run_cmd, cmd_list)
        pool.close()
        pool.join()
        self.merge_generate()

    def merge_generate(self):
        if self.method.lower() == 'miranda' and self.spiece.lower() == 'animal':
            result_list = []
            last_list = []
            with open('miranda_merge_out', 'w') as f, gzip.open("miranda_detail.txt.gz", 'w') as f_d:
                f.write("Query\tTarget\tScore\tEnergy\tQ_start\tQ_end\tT_start\tT_end\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                pair = ""
                for single in result_list:
                    detail = False
                    split_f = open(single,'r')
                    for line in split_f:
                        if not line.startswith('>>') and line.startswith('>'):
                            line_list = line.split()
                            if line_list[0:2] == last_list:
                                pass
                            else:
                                line_str = "\t".join(line_list[:-3])
                                last_list = line_list[0:2]
                                f.write(line_str.strip('>') + "\n")
                        elif line.startswith("Performing Scan:"):
                            pair = line
                        elif line.startswith("   Forward"):
                            f_d.write(pair)
                            detail = True
                        elif line.startswith("   Energy"):
                            detail = False

                        if detail:
                            f_d.write(line)

                    split_f.close()

        elif self.method.lower() == 'targetscan' and self.spiece.lower() == 'animal':
            query_dict = dict()
            target_list = list()
            for seq in SeqIO.parse(self.target, "fasta"):
                target_list.append(seq.id)
            for seq in SeqIO.parse(self.query, "fasta"):
                if len(seq.id.split("-miR-")) == 2:
                    mi_id = seq.id.split("-miR-")[1]
                    print mi_id
                else:
                    mi_id = seq.id
                    print mi_id
                query_dict[mi_id] = seq.id

            hsa2id = dict()

            if self.target_species != "9606":
                with open(os.path.dirname(self.targetscan) + '/' + self.target_species  + '_2hsa.txt', 'r') as id_convert:
                    for line in id_convert:
                        cols = line.split("\t")
                        if cols[2] == "":
                            pass
                        else:
                            if cols[2] in hsa2id:
                                hsa2id[cols[2]].append(cols[1])
                            else:
                                hsa2id[cols[2]] = [cols[1]]

            tar2context = dict()
            with open(self.target_context, 'r') as f:
                f.readline()
                for line in f:
                    cols = line.split("\t")
                    tar = cols[2] + "|" + cols[4] + "|" + cols[6] + "|" + cols[7]
                    tar2context[tar] = cols[8]


            target_nonconserved = os.path.join(os.path.dirname(self.target_context), "Nonconserved_Site_Context_Scores.txt")
            with open(target_nonconserved, 'r') as f:
                f.readline()
                for line in f:
                    cols = line.split("\t")
                    tar = cols[2] + "|" + cols[4] + "|" + cols[6] + "|" + cols[7]
                    tar2context[tar] = cols[8]

            target_nonconserved2 = os.path.join(os.path.dirname(self.target_context), "Nonconserved_Family_Info.txt")
            with open(self.targetscan, 'r') as f1, open(target_nonconserved2, 'r') as f1_nc,  open('targetscan_merge_out', 'w') as f2:
                f1.readline()
                f1_nc.readline()
                f2.write("#small_rna\ttarget\tUTR start\tUTR end\tMSA start\tMSA end\tSeed match\tPCT\tcontext score\n")
                for line in f1:
                    cols = line.strip().split("\t")
                    if cols[4] == self.target_species:
                        for mi in cols[0][4:].split("/"):
                            if mi in query_dict:
                                tar = "|".join([cols[3], query_dict[mi], cols[5], cols[6]])
                                if tar in tar2context:
                                    cscore = tar2context[tar]
                                else:
                                    cscore = ""

                                if self.target_species == "9606":
                                    if cols[3].split(".")[0] in target_list:
                                        f2.write("\t".join([query_dict[mi], cols[3].split(".")[0], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10], cscore]) + "\n")
                                else:
                                    gene = cols[1].split(".")[0]
                                    if gene in hsa2id:
                                        for transcript in hsa2id[gene]:
                                            if transcript in target_list:
                                                f2.write("\t".join([query_dict[mi], transcript, cols[5], cols[6], cols[7], cols[8], cols[9], cols[10], cscore]) + "\n")
                for line in f1_nc:
                    cols = line.strip().split("\t")
                    if cols[4] == self.target_species:
                        for mi in cols[0][4:].split("/"):
                            if mi in query_dict:
                                tar = "|".join([cols[3], query_dict[mi], cols[5], cols[6]])
                                if tar in tar2context:
                                    cscore = tar2context[tar]
                                else:
                                    cscore = ""

                                if self.target_species == "9606":
                                    if cols[3].split(".")[0] in target_list:
                                        f2.write("\t".join([query_dict[mi], cols[3].split(".")[0], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10], cscore]) + "\n")
                                else:
                                    gene = cols[1].split(".")[0]
                                    if gene in hsa2id:
                                        for transcript in hsa2id[gene]:
                                            if transcript in target_list:
                                                f2.write("\t".join([query_dict[mi], transcript, cols[5], cols[6], cols[7], cols[8], cols[9], cols[10], cscore]) + "\n")

        elif self.method.lower() == 'pita' and self.spiece.lower() == 'animal':
            result_list = []
            tmp_dict = dict()
            with open('pita_merge_out', 'w') as f:
                f.write("Target\tQuery\tT_start\tT_end\tSeed:Mismatch:Wobble\tScore\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output + "_pita_results_targets.tab")
                for single in result_list:
                    split_f = open(single,'r')
                    _ = split_f.readline()
                    for line in split_f:
                        line_list = line.split()
                        col1_col2='_'.join(line_list[:2])
                        if col1_col2 not in tmp_dict:
                            tmp_dict[col1_col2] = 1
                            line_str = "\t".join(line_list[0:5] + line_list[-1:])
                            f.write(line_str.strip('>') + "\n")
                    split_f.close()

        elif self.method.lower() == 'rnahybrid':
            result_list = []
            with open('rnahybrid_merge_out', 'w') as f, gzip.open("rnahybrid_detail.txt.gz", 'w') as f_d:
                f.write("Target\tQuery\tStart\tEnd\tEnergy\tPvalue\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                for single in result_list:
                    split_f = open(single,'r')
                    target = mirna = mfe = pvalue = '-'
                    start = 0
                    length = 0
                    for line in split_f:
                        # line_list = line.split(":")
                        # line_str = "\t".join(line_list[0:1] + line_list[2:3] + line_list[4:6])
                        if line.startswith("target too long"):
                            pass
                        else:
                            f_d.write(line)
                        if u'target:' in line:
                            target = line.split(':')[1].strip()
                        if u'miRNA :' in line:
                            mirna = line.split(':')[1].strip()
                        if u'mfe:' in line:
                            mfe = line.split(':')[1].strip().rstrip(' kcal/mol')
                        if u'length:' in line:
                            length = int(line.split(':')[1].strip())
                        if u'p-value:' in line:
                            pvalue = line.split(':')[1].strip()
                        if u'position ' in line:
                            start = int(line.split()[1].strip())
                            f.write(target + '\t' + mirna + '\t' + str(start) + '\t' + str(start+length-1) + '\t' + mfe + '\t' + pvalue + "\n")
                            target = mirna = mfe = pvalue = '-'
                    split_f.close()

        elif self.method.lower() == 'psrobot' and self.spiece.lower() == 'plant':
            result_list = []
            with open('psrobot_merge_out', 'w') as f, gzip.open("psrobot_detail.txt.gz", 'w') as f_d:
                f.write("Query\tScore\tTarget\tQ_start\tQ_end\tT_start\tT_end\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                for single in result_list:
                    split_f = open(single,'r')
                    for line in split_f:
                        f_d.write(line)
                        if line.startswith('>'):
                            line_list = line.split()[0:1] + line.split()[2:4]
                        elif line.startswith('Query'):
                            line_list += line.split()[1:2]  # + line.split()[3:4]
                            line_list += line.split()[3:4]
                        elif line.startswith('Sbjct'):
                            line_list += line.split()[1:2]
                            line_list += line.split()[3:4]
                            line_str = "\t".join(line_list)
                            f.write(line_str.strip('>') + "\n")
                    split_f.close()

        elif self.method.lower() == 'targetfinder' and self.spiece.lower() == 'plant':
            result_list = []
            with open('targetfinder_merge_out', 'w') as f, gzip.open("targetfinder_detail.txt.gz", 'w') as f_d:
                f.write("Query\tTarget\tQ_start\tQ_end\tScore\n")
                for split in self.split_fasta:
                    result_list.append(split + self.output)
                for single in result_list:
                    split_f = open(single,'r')
                    for line in split_f:
                        if line.startswith('No'):
                            pass
                        else:
                            cols = line.strip("\n").split("\t")
                            cols[1] = cols[1].split()[0]
                            line_list = cols[0:4] + cols[5:6]
                            line_str = "\t".join(line_list)
                            f.write(line_str + "\n")
                            f_d.write(cols[0] + " vs " + cols[1] + "\n\n")
                            f_d.write("target  5' " + cols[6] + " 3'\n")
                            f_d.write("           " + cols[7].replace(":", "|") + "   \n")
                            f_d.write("query   3' " + cols[8] + " 5'\n\n")
                    split_f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', type=str, metavar="query", required=True,
                        help="miRNA fasta file")
    parser.add_argument('-t', type=str, metavar="target", required=True,
                        help="target fasta file")
    parser.add_argument('-p', type=str, metavar="process", default=10,
                        help="max num of processes")
    parser.add_argument('-m', type=str, metavar="method", default="miranda",
                        help="miranda or PITA or RNAhybrid or psRobot or targetfinder")
    parser.add_argument('-n', type=str, metavar="split_num",default=100,
                        help="split query num")
    parser.add_argument('-o', type=str, metavar="output", default="target_predict",
                        help="prefix to the output")
    parser.add_argument('-miranda', type=str, metavar="miranda_path",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/miRNA/miRanda-3.3a/bin/miranda")
    parser.add_argument('-miranda_s', type=str, metavar="miranda_score",
                        default="140"),
    parser.add_argument('-miranda_e', type=str, metavar="miranda_energy",
                        default="-20"),
    parser.add_argument('-miranda_strict', type=str, metavar="miranda_strict",
                        default="on"),
    parser.add_argument('-pita', type=str, metavar="pita_path",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/miRNA/pita_v6/pita_prediction.pl")
    parser.add_argument('-rnahybrid', type=str, metavar="hybrid_path",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/miRNA/RNAhybrid-2.1.2/src/RNAhybrid")
    parser.add_argument('-rnahybrid_b', type=str, metavar="rnahybrid_num",
                        default="1"),
    parser.add_argument('-rnahybrid_e', type=str, metavar="rnahybrid_energy",
                        default="-20"),
    parser.add_argument('-rnahybrid_p', type=str, metavar="rnahybrid_pvalue",
                        default="0.01"),
    parser.add_argument('-rnahybrid_s', type=str, metavar="rnahybrid_species",
                        default="none"),
    parser.add_argument('-psrobot', type=str, metavar="psrobot",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/miRNA/psRobot_v1.2/psRobot_tar")
    parser.add_argument('-psrobot_ts', type=str, metavar="psrobot_score",
                        default="2.5")
    parser.add_argument('-targetfinder', type=str, metavar="targetfinder",
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/miRNA/targetfinder/targetfinder_threads.pl")
    parser.add_argument('-targetfinder_c', type=str, metavar="targetfinder_score",
                        default="4")
    parser.add_argument('-targetscan', type=str, metavar="targetscan",
                        default="/mnt/ilustre/users/sanger-dev/app/database/Predicted_Targets_Info.default_predictions.txt")
    parser.add_argument('-target_context', type=str, metavar="target_context",
                        default="/mnt/ilustre/users/sanger-dev/app/database/Predicted_Targets_Info.default_predictions.txt")
    parser.add_argument('-target_species', type=str, metavar="target_species",
                        default="Homo sapiens")
    parser.add_argument('-spiece', type=str,metavar="spiece_type",required=True,
                        help="type of spiece, animal or plant")
    args = parser.parse_args()
    toolbox=TagertPrediction(
        query_fa=args.q,
        target_fa=args.t,
        seq_num=args.n,
        method=args.m,
        pool=args.p,
        output=args.o,
        miranda=args.miranda,
        pita=args.pita,
        rnahybrid=args.rnahybrid,
        psRobot=args.psrobot,
        targetfinder=args.targetfinder,
        spiece=args.spiece,
        targetscan=args.targetscan,
        target_species=args.target_species,
        target_context=args.target_context
    )
    toolbox.set_options({
        "miranda_score": args.miranda_s,
        "miranda_energy": args.miranda_e,
        "miranda_strict": args.miranda_strict,
        "rnahybird_num": args.rnahybrid_b,
        "rnahybird_energy": args.rnahybrid_e,
        "rnahybird_pvalue": args.rnahybrid_p,
        "ps_robot_score": args.psrobot_ts,
        "targetfinder_score": args.targetfinder_c,
    })
    toolbox.splitfasta()
    toolbox.run_predict()
