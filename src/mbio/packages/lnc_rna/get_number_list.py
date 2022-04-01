# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os
import pandas as pd
import re

parser = OptionParser(description='Export number_list.txt from input directory')
parser.add_option('-i', '--input', dest='input', help='Input directory containing specified GTF files')
parser.add_option('-o', '--output', dest='output', help='Output number_list.txt')
(opts, args) = parser.parse_args()

def main(dir_in, file_out):
    df = pd.DataFrame(list(), columns=['#file_names', 'trans', 'genes'])
    for file_in in os.listdir(dir_in):
        if file_in.endswith('_out.gtf') or file_in.endswith('merged.gtf'):
            transcripts_set = set()
            genes_set = set()
            print 'INFO: start processing {}'.format(os.path.join(dir_in, file_in))
            for line in open(os.path.join(dir_in, file_in)):
                if line[0] != '#' and 'gene_id' in line and 'transcript_id' in line:
                    transcripts_set.add(re.search(r'transcript_id\s*"(\S+)";', line).group(1))
                    genes_set.add(re.search(r'gene_id\s*"(\S+)";', line).group(1))
            df = df.append(pd.Series(dict([
                ('#file_names', file_in), ('trans', len(transcripts_set)), ('genes', len(genes_set))
            ])), ignore_index=True)
    df.to_csv(file_out, sep='\t', index=None)
    print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(opts.input, opts.output)
    else:
        parser.print_help()