## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
# 20180725


import os
import sys
from biocluster.config import Config
import shutil


predict_fa = sys.argv[1]
rnafold = Config().SOFTWARE_DIR + '/bioinfo/rna_pro/RNAfold'
parafly = Config().SOFTWARE_DIR + '/program/parafly-r2013-01-21/src/ParaFly'
relplot = Config().SOFTWARE_DIR + '/bioinfo/rna_pro/ViennaRNA-2.3.1/src/Utils/relplot.pl'
perl_path = Config().SOFTWARE_DIR + '/program/perl-5.24.0/bin/perl'
os.environ['PATH'] += ':' + perl_path

fold_cmd = rnafold + ' ' + '-p < ' + predict_fa +' > RNAfold.str.txt'
with open('fold.bash', 'w') as fold_sh:
    fold_sh.write(fold_cmd)

os.system('/bin/bash fold.bash')
if os.path.exists('RNAfold_pdf'):
    shutil.rmtree('RNAfold_pdf')
os.mkdir('RNAfold_pdf')

with open(predict_fa, 'r') as fa_r, open('parafly_list', 'w') as para_w:
    for line in fa_r.readlines():
        if line.startswith('>'):
            id = line.strip().lstrip('>')
            para_w.write(perl_path + ' ' + relplot + ' ' + id + '_ss.ps ' + id + '_dp.ps >RNAfold_pdf/' + id + '.ps'+ ' && ' + 'ps2pdf RNAfold_pdf/' + id + '.ps RNAfold_pdf/' + id + '.pdf\n')

os.system(parafly +' -c parafly_list -CPU 20')

with open('RNAfold.str.txt', 'r') as fold_r, open('RNAfold_pdf/RNAfold_stat.txt', 'w') as fold_stat:
    fold_stat.write('rRNA ID\tMFE Structure in Ensemble\tEnsemble diversity\n')
    for block in fold_r.read().split('\n>'):
        block = block.split('\n')
        id = block[0].lstrip('>')
        for l in block:
            if u'ensemble' in l:
                mfe = diver = '-'
                try:
                    mfe = l.split(';')[0].strip().split(' ')[-1]
                    diver = l.split(';')[1].strip().split(' ')[-1]
                except:
                    pass
                fold_stat.write(id + '\t' + mfe + '\t' + diver + '\n')



