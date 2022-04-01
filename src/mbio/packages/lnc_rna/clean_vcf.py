# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Convert bam path or variant name to sample name in header of VCF file')
parser.add_option('-i', '--input', dest='input', help='Input raw VCF file')
parser.add_option('-o', '--output', dest='output', help='Output clean VCF file')
(opts, args) = parser.parse_args()

def main(file_in, file_out):
    print 'INFO: start reading {}'.format(file_in)
    with open(file_out, 'w') as f:
        for line in open(file_in):
            if line.startswith('#CHROM'):
                items = line.strip().split('\t')
                header = items[:9]
                for item in items[9:]:
                    if item.endswith('.bam'):
                        header.append(os.path.basename(item).rstrip('.bam'))
                    elif '.variant' in item:
                        header.append('.'.join(item.split('.')[:-1]))
                    else:
                        raise Exception('ERROR: find unvalid item in header -> {}'.format(item))
                f.write('{}\n'.format('\t'.join(header)))
            elif line[0] != '#':
                items = line.strip().split('\t')
                if 'AD' in items[8] and 'DP' in items[8]:
                    f.write(line)
            else:
                f.write(line)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(opts.input, opts.output)
    else:
        parser.print_help()