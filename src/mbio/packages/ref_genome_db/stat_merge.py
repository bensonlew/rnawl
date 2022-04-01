# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd

parser = OptionParser(description='Merge stat results from stat_gtf and seqkit_fx2tab')
parser.add_option('-t', '--tabular', dest='tabular' ,help='Input tabular file')
parser.add_option('-f', '--fx2tab', dest='fx2tab', help='Input fx2tab file')
parser.add_option('-o', '--output', dest='output', help='Output genome_stat file')
(opts, args) = parser.parse_args()

def main(fx2tab, tabular, output):
    df1 = pd.DataFrame([line2items(line) for line in open(fx2tab)], columns=['Chr', 'Size(Mb)', 'GC%'])
    df1.dropna(axis=1, inplace=True)
    df1.iloc[:, 1] = df1.iloc[:, 1].map(lambda x: round(float(x)/10.0**6, 2))
    df2 = pd.read_table(tabular)
    if 'Chr' in df2 and df2['Chr'].dtypes == 'int64':
        df2['Chr'] = df2['Chr'].apply(str)
    df = pd.merge(df1, df2, how='left')
    df = df.reindex(columns=['Chr', 'Size(Mb)', 'GC%', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.iloc[:, 1] = df.iloc[:, 1].map(lambda x: '{:.2f}'.format(x))
    df.iloc[:, 3:7] = df.iloc[:, 3:7].applymap(lambda x: '_' if pd.isnull(x) else int(x))
    df.to_csv(output, sep='\t', index=None)

def line2items(line):
    items = line.strip().split('\t')
    return items[0], items[-2], items[-1]

if __name__ == '__main__':
    if opts.tabular and opts.fx2tab and opts.output:
        main(opts.fx2tab, opts.tabular, opts.output)
    else:
        parser.print_help()
