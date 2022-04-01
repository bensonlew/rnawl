# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os
import pandas as pd

parser = OptionParser(description='Process results of lncRNA sequence conservation analysis')
parser.add_option('-l', '--liftover', dest='liftover', help='input LIFTOVER result BED file')
parser.add_option('-b', '--bedtools', dest='bedtools', help='input BEDTOOLs intersect out file')
parser.add_option('-t', '--type', dest='type', help='input file containing name and its corresponding type')
parser.add_option('-o', '--output', dest='output', help='output tabular file')
(opts, args) = parser.parse_args()

def main(liftover_out, bedtools_out, type_in, file_out):
    if os.path.getsize(liftover_out):
        print 'INFO: start reading {}'.format(liftover_out)
        dfl = pd.read_table(liftover_out, header=None)
        dfl = dfl.reindex([3, 0, 5, 1, 2, 'location'], axis=1)
        dfl = dfl.rename(columns={3: 'lncrna_id', 0: 'chromosome', 5: 'strand', 1: 'start', 2: 'end'})
        for i in dfl.index:
            dfl.loc[i, 'location'] = '{}[{}]{}-{}'.format(*dfl.loc[i, ['chromosome', 'strand', 'start', 'end']])
        dfl = dfl.drop(['chromosome', 'strand', 'start', 'end'], axis=1)
        if os.path.getsize(bedtools_out):
            print 'INFO: start reading {}'.format(bedtools_out)
            dfb = pd.read_table(bedtools_out, header=None)
            dfb = dfb.reindex([3, 15, 24], axis=1)
            dfb = dfb.rename(columns={3: 'lncrna_id', 15: 'transcript_id', 24: 'overlap_length'})
            df = dfl.merge(dfb, on='lncrna_id', how='right')
        else:
            dfl['transcript_id'] = '-'
            dfl['overlap_length'] = 0
            df = dfl
        df = df.merge(pd.read_table(type_in, names=['lncrna_id', 'lncrna_type'], usecols=[0, 3]), on='lncrna_id', how='left')
        df.to_csv(file_out, sep='\t', index=False)
    else:
        print 'INFO: start generating blank result because {} is empty'.format(liftover_out)
        open(file_out, 'w').write('lncrna_id\tlocation\ttranscript_id\toverlap_length\tlncrna_type\n')
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 4, ['liftover', 'bedtools', 'type', 'output'])):
        main(opts.liftover, opts.bedtools, opts.type, opts.output)
    else:
        parser.print_help()
