# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import pandas as pd
import numpy as np
from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class MatrixFile(File):
    def __init__(self):
        super(MatrixFile, self).__init__()

    def check(self):
        super(MatrixFile, self).check()

    def export_group_matrix(self, out_file, group_file):
        df = pd.read_table(self.path)
        samples = list()
        groups = list()
        for n, line in enumerate(open(group_file)):
            if n != 0:
                sample, group = line.strip().split('\t')
                samples.append(sample)
                groups.append(group)
        df = df.reindex(['seq_id'] + samples, axis=1).set_index('seq_id')
        columns = pd.MultiIndex.from_arrays([groups, samples], names=['group', 'sample'])
        df = pd.DataFrame(np.array(df), index=df.index, columns=columns)
        df = df.groupby(level='group', axis=1).mean()
        df.to_csv(out_file, sep='\t')

    def rename_by_map(self, out_file, map_file):
        df = pd.read_table(self.path)
        columns = dict(line.strip().split('\t') for line in open(map_file) if line[0] != '#')
        df = df.rename(columns=columns)
        df.to_csv(out_file, sep='\t', index=False)
