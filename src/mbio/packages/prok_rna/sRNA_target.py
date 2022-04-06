## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
# 20180724


import os
import re
import sys
from collections import defaultdict
from biocluster.config import Config
from mako.template import Template


class srna_target(object):
    def __init__(self, predict_fa, genome_bed, genome_fa):
        self.predict_fa = predict_fa
        self.genome_bed = genome_bed
        self.genome_fa = genome_fa
        self.bedtool_path = Config().SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
        self.parafly = Config().SOFTWARE_DIR + '/program/parafly-r2013-01-21/src/ParaFly'
        self.rnahybrid = Config().SOFTWARE_DIR + '/bioinfo/rna_pro/RNAhybrid-2.1.2/src/RNAhybrid'
        self.rnaplex = Config().SOFTWARE_DIR + '/bioinfo/rna_pro/RNAplex'
        self.fastalength = Config().SOFTWARE_DIR + '/bioinfo/annotation/scripts/fastalength'

    def get_target_fa(self):
        seq_cmd = r'''awk '{if($7-35<0){$7=35};if($8-20<0){$8=20}; if($6=="+"){print $1"\t"$7-35"\t"$7+21"\t"$4"\t0\t+"}else{print $1"\t"$8-20"\t"$8+36"\t"$4"\t0\t-"}}' ${genome_bed} >genome.tar.bed
${fastalength} ${genome_fa} | awk '{print $2"\t0\t"$1}' > genome.bed
${bedtools} intersect -a genome.tar.bed -b genome.bed > genome.tar.range.bed
${bedtools} getfasta -name -fi ${genome_fa} -bed genome.tar.range.bed -fo genome.tar.fa
'''
        f = Template(seq_cmd)
        seq_cmd = f.render(genome_bed=self.genome_bed,
                           genome_fa=self.genome_fa,
                           fastalength=self.fastalength,
                           bedtools=self.bedtool_path,
                             )
        with open('get_target.bash', 'w') as bash:
            bash.write(seq_cmd)
        os.system('bash get_target.bash')

    def split_fasta(self, fa):
        sp_fas = list()
        n = 1
        na2fa = defaultdict(int)
        with open(fa, 'r') as fa_r:
            fa_info = fa_r.read().lstrip('>').split('\n>')
            limit_num = len(fa_info)/10
            for i in range(len(fa_info)):
                with open('splited' + str(n) + '.fa', 'a') as sp_w:
                    if n not in na2fa:
                        na2fa[n] = 0
                    if na2fa[n] == int(limit_num) - 1:
                        sp_w.write('>' + fa_info[i] + '\n')
                        na2fa[n] += 1
                        sp_fas.append('splited' + str(n) + '.fa')
                        n += 1
                    else:
                        sp_w.write('>' + fa_info[i] + '\n')
                        na2fa[n] += 1
            if 'splited' + str(n) + '.fa' not in sp_fas and os.path.exists('splited' + str(n) + '.fa'):
                sp_fas.append('splited' + str(n) + '.fa')
            return sp_fas

    def run_target(self, fa_list):
        target_cmd = ''
        for fa in fa_list:
            target_cmd += self.rnahybrid + ' -t genome.tar.fa -q ' + fa + ' -c -p 0.05 -m 50000 -n 5000 -d 0.05,1 > ' + fa + 'RNAhybrid.txt\n'
            target_cmd += self.rnaplex + ' -t genome.tar.fa -q ' + fa + ' > ' + fa + 'RNAplex.txt\n'
        with open('run_target.bash', 'w') as tar_sh:
            tar_sh.write(target_cmd)
        cmd = self.parafly + ' -c run_target.bash -CPU 35'
        os.system(cmd)

    def merge_result(self, fa_list):
        with open('RNAhybrid_merge', 'w') as hy_merge:
            hy_merge.write('sRNA ID\tTarget Gene ID\tLen miR\tLen Tar\tEnergy\tP-value\n')
            for fa in fa_list:
                with open(fa + 'RNAhybrid.txt', 'r') as hy_r:
                    for line in hy_r.readlines():
                        line = line.strip().split(':')
                        if len(line) >= 6 and float(line[4]) < -40 and float(line[5]) < 0.01:
                            hy_merge.write(line[2] + '\t'+ line[0] + '\t' + line[3] + '\t' + line[1] + '\t' + line[4] + '\t' + line[5] + '\n')
        with open('RNAplex_merge', 'w') as plex_merge:
            plex_merge.write('sRNA ID\tTarget Gene ID\tAlignment\tEnergy\tQuery Start\tQuery End\tTarget Start\tTarget End\n')
            for fa in fa_list:
                with open(fa + 'RNAplex.txt', 'r') as plex_r:
                    plex_info = plex_r.readlines()
                    i = 0
                    while i+2 < len(plex_info):
                        target = plex_info[i].strip().lstrip('>')
                        srna = plex_info[i+1].strip().lstrip('>')
                        re_group = re.match('(\S+)\s+(\d+),(\d+)\s+:\s+(\d+),(\d+)\s+\((\S+)\)$', plex_info[i+2].strip())
                        if re_group:
                            align = re_group.group(1)
                            energy = re_group.group(6)
                            qstart = re_group.group(4)
                            qend = re_group.group(5)
                            tstart = re_group.group(2)
                            tend = re_group.group(3)
                            if float(energy) < -40:
                                plex_merge.write(srna + '\t' + target + '\t' + align + '\t' + energy + '\t' + qstart + '\t' + qend + '\t' + tstart + '\t' + tend + '\n')
                        i += 3

    def combine_two(self):
        hybrid_s2t = defaultdict(set)
        plex_s2t = defaultdict(set)
        hy_plex_s2t = defaultdict(set)
        with open('RNAplex_merge', 'r') as plex_r:
            _ = plex_r.readline()
            for line in plex_r.readlines():
                line = line.strip().split('\t')
                plex_s2t[line[0]].add(line[0] + '\t' + line[1])
        with open('RNAhybrid_merge', 'r') as hy_r:
            _ = hy_r.readline()
            for line in hy_r.readlines():
                line = line.strip().split('\t')
                hybrid_s2t[line[0]].add(line[0] + '\t' + line[1])
        all_srna = set(list(hybrid_s2t.keys()) + list(hybrid_s2t.keys()))
        for srna in all_srna:
            if srna in hybrid_s2t and srna in plex_s2t:
                if hybrid_s2t[srna] & plex_s2t[srna]:
                    hy_plex_s2t[srna] = hybrid_s2t[srna] & plex_s2t[srna]
        with open('combine_RNAplex_RNAhybrid', 'w') as com_w:
            com_w.write('sRNA ID\tTarget Gene ID\tDatabase\n')
            for srna in sorted(all_srna):
                if srna in hy_plex_s2t:
                    for i in hy_plex_s2t[srna]:
                        com_w.write(i + '\t' + 'RNAphybrid, RNAplex' + '\n')
                    for i in hybrid_s2t[srna]:
                        if i not in hy_plex_s2t[srna]:
                            com_w.write(i + '\t' + 'RNAphybrid' + '\n')
                    for i in plex_s2t[srna]:
                        if i not in hy_plex_s2t[srna]:
                            com_w.write(i + '\t' + 'RNAplex' + '\n')
                else:
                    if srna in hybrid_s2t:
                        for i in hybrid_s2t[srna]:
                            com_w.write(i + '\t' + 'RNAphybrid' + '\n')
                    if srna in plex_s2t:
                        for i in plex_s2t[srna]:
                            com_w.write(i + '\t' + 'RNAplex' + '\n')
    def run(self):
        self.get_target_fa()
        splited_list = self.split_fasta(self.predict_fa)
        self.run_target(splited_list)
        self.merge_result(splited_list)
        self.combine_two()


if __name__ == '__main__':
    if len(sys.argv) != 4:
        exit('USAGE %s <predict_fa> <genome_bed> <genome_fa>' %sys.argv[0])

    predict_fa = sys.argv[1]
    genome_bed = sys.argv[2]
    genome_fa = sys.argv[3]
    TARGET = srna_target(predict_fa, genome_bed, genome_fa)
    TARGET.run()
