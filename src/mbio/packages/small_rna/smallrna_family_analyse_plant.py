## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
#last_modify:20180706

import argparse
from collections import defaultdict


class family_analyse(object):
    def __init__(self, miRfamily, miR_pre_mature, mir_list_tab, mature_fa, denovo_mature_fa):
        self.miRfamily = miRfamily
        self.miR_pre_mature = miR_pre_mature
        self.mir_list_tab = mir_list_tab
        self.mature_fa = mature_fa
        self.denovo_mature_fa = denovo_mature_fa

    def tr_perl(self, target_str, old, after):
        trans = ''
        for i in target_str:
            if i in old:
                i = after[old.index(i)]
            trans += i
        return trans

    def fasta_dict(self, fa_file):
        id2seq = {}
        with open(fa_file) as fa_r:
            fas =  fa_r.read().lstrip('>').split('\n>')
            for fa in fas:
                fa = fa.split('\n')
                id = fa[0].strip()
                seq = ''.join(fa[1:]).replace('\n', '')
                id2seq[id] = seq
        return id2seq

    def fam2pre_dict(self, mifam):
        fam2pre = defaultdict(list)
        with open(mifam, 'r') as mifam_r:
            mifams = mifam_r.readlines()
            for line in mifams:
                line = line.strip().split('\t')
                fam2pre[line[2]].append(line[0])
        return fam2pre

    def spe_pre_mi_dicts(self, mature_pre):
        pre2mi = defaultdict(list)
        spe2pre = defaultdict(list)
        with open(mature_pre, 'r') as mat_r:
            mats = mat_r.readlines()
            for line in mats:
                line = line.strip().split('\t')
                spe2pre['%s\t%s'%(line[0], line[1])].append(line[2])
                for mi in line[3:]:
                    pre2mi[line[2]].append(mi)
        return pre2mi, spe2pre

    def creat_known_txt(self, mir_list_tab, fam2pre, pre2mi, spe2pre):
        def sub_in(target, test_list):
            result = False
            for i in test_list:
                if target in i:
                    result = True
                    break
            return result
        fam_txt = ''
        spe_txt = ''
        families = list()
        with open(mir_list_tab, 'r') as mir_r:
            _ = mir_r.readline()
            mirs = mir_r.readlines()
            for line in mirs:
                tmp_pre = list()
                tmp_fam = list()
                mir = line.strip().split('\t')[0]
                for key in pre2mi:
                    if mir in pre2mi[key] or sub_in(mir, pre2mi[key]):
                        tmp_pre.append(key)
                        for fam in fam2pre:
                            if key in fam2pre[fam]:
                                tmp_fam.append(fam)
                if tmp_pre:
                    fam_txt = fam_txt + mir + '\t' + ';'.join(set(tmp_pre)) + '\t'
                else:
                    fam_txt = fam_txt + mir + '-' + '\t'
                if tmp_fam:
                    families += tmp_fam
                    fam_txt = fam_txt + ';'.join(set(tmp_fam)) + '\n'
                else:
                    fam_txt = fam_txt + '-' + '\n'
        species = sorted(spe2pre.keys())
        families = sorted(set(families))
        spe_txt = 'Fullname\tAbbreviation\t%s\n'%('\t'.join(families))
        for spe in species:
            spe_txt += '%s\t'%spe
            for fam in families:
                if set(fam2pre[fam]) & set(spe2pre[spe]):
                    spe_txt += '+\t'
                else:
                    spe_txt += '-\t'
            spe_txt += '\n'
            # print(fam_txt, spe_txt)
        return fam_txt, spe_txt

    def creat_novel_txt(self, mature_fa, novo_fa, fam_txt):
        # print(fam_txt)
        def matct(fam2seed, seq):
            score2fam = defaultdict(list)
            for fam in fam2seed:
                score = 0
                for n, i in enumerate(seq):
                    if i == fam2seed[fam][n]:
                        score += 1
                        score2fam[score].append(fam)
            if max(score2fam) >= 6:
                return score2fam[max(score2fam)][0]
            else:
                return 0
        novel_txt = ''
        mature_dict = self.fasta_dict(mature_fa)
        novo_dict = self.fasta_dict(novo_fa)
        fam2seed = dict()
        for line in fam_txt.split('\n'):
            line = line.strip().split('\t')
            if line[-1] != '-' and line[-1] not in fam2seed and line[0]:
                fam2seed[line[-1]] = mature_dict[line[0]][1:7]
            # print(fam2seed)
        for nov in novo_dict:
            seq = self.tr_perl(novo_dict[nov], 'aucg', 'ATCG')
            seq_z = seq[1:7]
            seq_f = self.tr_perl(seq[::-1], 'ATCGU', 'TAGCT')[1:7]
            matct_z = matct(fam2seed, seq_z)
            matct_f = matct(fam2seed, seq_f)
            is_z = 0
            # print(nov, matct_z, matct_f)
            if matct_z:
                novel_txt += '%s\t%s\n' % (nov, matct_z)
                is_z += 1
            if matct_f and is_z:
                novel_txt += '%s\t%s\n' % (nov, matct_f)
            if not matct_f and not matct_z:
                novel_txt += '%s\t%s\n' % (nov, '-')
            # print(novel_txt)
        return novel_txt

    def run(self):
        fam2pre = self.fam2pre_dict(self.miRfamily)
        pre2mi, spe2pre = self.spe_pre_mi_dicts(self.miR_pre_mature)
        fam_txt, spe_txt = self.creat_known_txt(self.mir_list_tab, fam2pre, pre2mi, spe2pre)
        novel_txt = self.creat_novel_txt(self.mature_fa, self.denovo_mature_fa, fam_txt)
        return fam_txt, spe_txt, novel_txt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
    small_rna定量之后的家族分析
    """)
    parser.add_argument(
        "-fam",
        dest="family",
        required=True,
        type=str,
        help="前体和家族的对应关系文件")
    parser.add_argument(
        "-pre",
        dest="pre",
        required=True,
        type=str,
        help="前体和成熟体的对应关系文件")
    parser.add_argument(
        "-mir",
        dest="mir",
        type=str,
        required=True,
        help="含有全部smallRNA ID的文件")
    parser.add_argument(
        "-matfa",
        dest="mature_fa",
        type=str,
        required=True,
        help="含有全部已知mature_smallrna的序列文件")
    parser.add_argument(
        "-novofa",
        dest="novo_fa",
        type=str,
        required=True,
        help="含有全部新的mature_smallrna的序列文件")
    args = parser.parse_args()

    FA = family_analyse(args.family, args.pre, args.mir, args.mature_fa, args.novo_fa)
    fam_txt, spe_txt, novel_txt = FA.run()
    with open('known_miR_family.xls', 'w') as known_w:
        known_w.write('miRNA\tmiR_pre\tfamily\n')
        known_w.write(fam_txt)
    with open('family.species.xls', 'w') as fam_w:
        fam_w.write(spe_txt)
    with open('novel_miR_family.xls', 'w') as novel_w:
        novel_w.write('miRNA\tFamily\n')
        novel_w.write(novel_txt)




