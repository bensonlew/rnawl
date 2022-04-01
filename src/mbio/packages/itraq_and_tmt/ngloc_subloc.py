#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/8/1 17:40
@file    : ngloc_subloc.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import os
import sys
import shutil
from mako.template import Template
from biocluster.config import Config
import pandas as pd
import glob


this_file_dir = os.path.dirname(os.path.realpath(__file__))
ngloc = Config().SOFTWARE_DIR + '/bioinfo/itraq_and_tmt/ngLoc_v1.0/'
try:
    _, expfasta, gram = sys.argv
except:
    gram = 'neg'
    expfasta = os.path.join(ngloc, 'DataSets', 'TestSet', 'test.fa')
expfasta = os.path.realpath(expfasta)

if os.path.exists('ngloc_rundir'):
    shutil.rmtree('ngloc_rundir')
os.makedirs('ngloc_rundir')
os.chdir('ngloc_rundir')

shutil.copy(os.path.join(ngloc, 'ngLOC'), '.')
shutil.copy(expfasta, '.')
sed_cmd = r'''sed -i 's/\r//g' %s''' % os.path.basename(expfasta)
os.system(sed_cmd)
print(sed_cmd)
if gram.lower() == 'neg':
    trainfasta = os.path.join(ngloc, 'DataSets', 'TrainSet', 'bac_gram_neg_unq.fa')
else:
    trainfasta = os.path.join(ngloc, 'DataSets', 'TrainSet', 'bac_gram_pos_unq.fa')
shutil.copy(trainfasta, '.')

if os.path.exists('Log'):
    shutil.rmtree('Log')
os.makedirs('Log')

f = Template(filename= this_file_dir + '/config.ini')
text = f.render(expfasta=os.path.basename(expfasta), trainfasta=os.path.basename(trainfasta))
with open('config.ini', 'w') as iw:
    iw.write(text)

os.system('./ngLOC')

def convert_ngloc2subloc(ngloc_res, out='multiloc.xls'):
    convert_dict = {
        'CYT': 'cytoplasmic',
        'PLA': 'plasma membrane',
        'EXC': 'excellur',
        'NUC': 'nuclear',
    }
    ng_df = pd.read_csv(ngloc_res, sep='\t', index_col=0)
    ng_df.index = ng_df.index.map(lambda x:x.strip('>'))
    for ind in ng_df.index:
        try:
            pre1 = ng_df.loc[ind, 'Pred1']
        except:
            pre1 = '_'
        try:
            pro1 = float(ng_df.loc[ind, 'Prob1'])/100
        except:
            pro1 = 0
        ng_df.loc[ind, 'col1'] = pre1 + ': ' + str(round(pro1, 2))
        try:
            pre2 = ng_df.loc[ind, 'Pred2']
        except:
            pre2 = '_'
        try:
            pro2 = float(ng_df.loc[ind, 'Prob2'])/100
        except:
            pro2 = 0
        ng_df.loc[ind, 'col2'] = pre2 + ': ' + str(round(pro2, 2))
        try:
            pre3 = ng_df.loc[ind, 'Pred3']
        except:
            pre3 = '_'
        try:
            pro3 = float(ng_df.loc[ind, 'Prob3'])/100
        except:
            pro3 = 0
        ng_df.loc[ind, 'col3'] = pre3 + ': ' + str(round(pro3, 2))
    subloc_df = ng_df[['col1', 'col2', 'col3']]
    subloc_df.to_csv(out, sep='\t', index=True, header=False)

ng_file = glob.glob('./Log/*Predictions.txt')[0]
convert_ngloc2subloc(ng_file, '../multiloc.xls')
os.chdir('..')