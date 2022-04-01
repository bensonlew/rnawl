# -*- coding: utf-8 -*-

import glob
import logging
import os

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

work_dir = os.getcwd()

feed_dict = {'T.together': 'T.together.output', 'T.split': 'T.split.output',
             'G.together': 'G.together.output', 'G.split': 'G.split.output'}

data = list()
for key, _dir in feed_dict.items():
    document = {'method': key}
    output_dir = os.path.join(work_dir, _dir)
    for table in glob.glob(os.path.join(output_dir, '*.deseq2.txt')):
        df = pd.read_table(table)
        df = df[pd.notna(df.padj)]
        df = df[(abs(df.log2FoldChange) > 2) & (df.padj < 0.05)]
        sig_num = df.shape[0]
        vs_pair = os.path.basename(table)[:-11]
        document[vs_pair] = sig_num
        logger.info('method: {}; pair: {}'.format(key, vs_pair))
    data.append(document)

df = pd.DataFrame(data).set_index('method').T
df.index.name = 'Compare'
df.to_csv(os.path.join('stat.txt'), sep='\t')
