# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Merge CMSCAN TSV files')
parser.add_option('-i', '--input', dest='input', help='input CMSCAN TSV list file')
parser.add_option('-k', '--known', dest='known', help='input known lncRNA id list file')
parser.add_option('-n', '--novel', dest='novel', help='input novel lncRNA id list file')
parser.add_option('-t', '--t2g', dest='t2g', help='input TXPT2GENE file')
parser.add_option('-o', '--output', dest='output', help='output tabular file')
(opts, args) = parser.parse_args()

def main(file_in, known_id_file, novel_id_file, t2g_file, file_out):
    print 'INFO: start reading {}'.format(file_in)
    df = pd.concat([pd.read_table(line.strip()) for line in open(file_in)], ignore_index=True)
    df['e_value'] = df['e_value'].apply(float)
    df = df.sort_values(by=['family_id', 'e_value', 'lncrna_id'])
    print 'INFO: start reading {} and {}'.format(known_id_file, novel_id_file)
    df = df.merge(pd.concat([
        pd.DataFrame(
            data=[(i, 'known') for i in [line.strip() for line in open(known_id_file)]],
            columns=['lncrna_id', 'lncrna_type']
        ),
        pd.DataFrame(
            data=[(i, 'novel') for i in [line.strip() for line in open(novel_id_file)]],
            columns=['lncrna_id', 'lncrna_type']
        )
    ], ignore_index=True), on='lncrna_id', how='left')
    print 'INFO: start reading {}'.format(t2g_file)
    df = df.merge(pd.DataFrame(
        data=[line.strip().split('\t') for line in open(t2g_file)], columns=['lncrna_id', 'lncrna_gene_id']
    ), on='lncrna_id', how='left')
    df = df.reindex(columns=[
        'lncrna_id', 'lncrna_gene_id', 'family_name', 'family_id', 'lncrna_start', 'lncrna_end',
        'e_value', 'score', 'lncrna_type'
    ])
    df.to_csv(file_out, sep='\t', index=False)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 5, ['input', 'known', 'novel', 't2g', 'output'])):
        main(opts.input, opts.known, opts.novel, opts.t2g, opts.output)
    else:
        parser.print_help()
