# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Process CMSCAN tblout')
parser.add_option('-i', '--input', dest='input', help='input CMSCAN tblout')
parser.add_option('-o', '--output', dest='output', help='output TSV file')
(opts, args) = parser.parse_args()

def main(file_in, file_out):
    print 'INFO: start reading {}'.format(file_in)
    data = list()
    for line in open(file_in):
        if line[0] != '#':
            items = line.strip().split()
            if '!' in items:
                data.append(items[:items.index('!')])
    columns = [
        'idx', 'family_name', 'family_id', 'lncrna_id', 'accession', 'clan_name', 'mdl', 'mdl_from', 'mdl_to',
        'lncrna_start', 'lncrna_end', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'e_value'
    ]
    df = pd.DataFrame(data=data, columns=columns)
    df = df.reindex(['lncrna_id', 'family_name', 'family_id', 'lncrna_start', 'lncrna_end', 'e_value', 'score'], axis=1)
    df['e_value'] = df['e_value'].apply(float)
    df = df.reindex(columns=['lncrna_id', 'family_name', 'family_id', 'lncrna_start', 'lncrna_end', 'e_value', 'score'])
    df.to_csv(file_out, sep='\t', index=False)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 2, ['input', 'output'])):
        main(opts.input, opts.output)
    else:
        parser.print_help()
