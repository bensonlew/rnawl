# -*- coding:utf-8 -*-

import glob
import gzip
import logging
import os
import sys

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: {} <last.result DIR> <mapper.name STR>'.format(os.path.basename(sys.argv[0])))
    dirs_path = sys.argv[1]
    mapper_name = sys.argv[2]
    dir_paths = glob.glob(os.path.join(dirs_path, '*'))
    if mapper_name.lower() == 'bismark':
        for path in dir_paths:
            if os.path.basename(path) == 'output':
                continue
            sp_name = os.path.basename(path)
            fp = os.path.join(path, 'result/{}_pe.deduplicated.bismark.cov'.format(sp_name))
            if not os.path.isfile(fp) and os.path.isfile('{}.gz'.format(fp)):
                logger.debug('Unzip {}.gz'.format(fp))
                with gzip.GzipFile('{}.gz'.format(fp)) as gz_file:
                    open(fp, 'wb+').write(gz_file.read())
            logger.debug('Process {}'.format(fp))
            df = pd.read_table(fp, header=None)
            df[6] = df[4] + df[5]
            df.rename({0: 'chr', 1: 'pos', 4: 'X', 6: 'N'}, axis=1, inplace=True)
            df = df.reindex(['chr', 'pos', 'N', 'X'], axis=1)
            output_dir = os.path.join(path, 'output')
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            df.to_csv(os.path.join(output_dir, '{}.dss.input.txt'.format(sp_name)), sep='\t', index=False)
    elif mapper_name.lower() == 'bsmap':
        for path in dir_paths:
            if os.path.basename(path) == 'output':
                continue
            sp_name = os.path.basename(path)
            fp = os.path.join(path, '{}.txt'.format(sp_name))
            logger.debug('Parse {}'.format(sp_name))
            df = pd.read_table(fp)
            df.rename({'C_count': 'X', 'CT_count': 'N'}, axis=1, inplace=True)
            df = df.reindex(['chr', 'pos', 'N', 'X'], axis=1)
            output_dir = os.path.join(path, 'output')
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            df.to_csv(os.path.join(output_dir, '{}.dss.input.txt'.format(sp_name)), sep='\t', index=False)
