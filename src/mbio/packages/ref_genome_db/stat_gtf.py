# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import re
import pandas as pd

parser = OptionParser(description='Evaluate gtf file')
parser.add_option('-i', '--input', dest='input' ,help='Input gtf file')
parser.add_option('-o', '--output', dest='output', help='Output tabular file')
parser.add_option('-s', '--source', dest='source', help='Source of gtf file [ensemble]')
(opts, args) = parser.parse_args()

def process_ensemble_gtf(gtf, tabular):
    dct = dict()
    with open(gtf) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                feature = items[2]
                attributes= items[8]
                if feature == 'gene':
                    if seqname in dct:
                        dct[seqname]['Gene'] += 1
                    else:
                        dct[seqname] = {'Chr': seqname, 'Gene': 1, 'ProteinCoding': 0, 'Pseudogene': 0, 'OtherRNA': 0}
                    m = re.match(r'.*biotype\s"(\S+)";', attributes)
                    if m:
                        biotype = m.group(1)
                        if biotype == 'protein_coding':
                            dct[seqname]['ProteinCoding'] += 1
                        elif biotype == 'pseudogene':
                            dct[seqname]['Pseudogene'] += 1
                        elif 'RNA' in biotype:
                            dct[seqname]['OtherRNA'] += 1
    df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.sort_values(by='Chr', inplace=True)
    df.to_csv(tabular, sep='\t', index=None)

if __name__ == '__main__':
    if opts.input and opts.output:
        if opts.source.lower() == 'ensemble':
            process_ensemble_gtf(opts.input, opts.output)
        else:
            parser.print_help()
    else:
        parser.print_help()