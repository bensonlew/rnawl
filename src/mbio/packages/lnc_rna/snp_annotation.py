# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from optparse import OptionParser
import re
import os

parser = OptionParser(description='Integrate annovar results of SNP annotation')
parser.add_option('-i', '--input', dest='input', help='Input VCF file')
parser.add_option('-e', '--evf', dest='evf', help='The exonic_variant_function file')
parser.add_option('-v', '--vf', dest='vf', help='Input variant_function file')
parser.add_option('-o', '--output', dest='output', help='Output tabular file')
(opts, args) = parser.parse_args()

def main(vcf_in, exonic_variant_function, variant_function, xls_out):
    vcf_anno_dict = {}
    header_items = [
        'CHROM', 'START', 'END', 'REF', 'ALT', 'QUAL', 'ANNO', 'GENE(in or nearby)', 'Depth', 'MUT_type', 'MUT_info'
    ]
    print 'INFO: start reading {}'.format(vcf_in)
    with open(vcf_in) as vcf:
        for n, line in enumerate(vcf):
            if line[:2] == '##':
                continue
            elif re.match(r'#CHROM', line):
                samples = [i for i in line.strip().split('\t')[9:]]
            else:
                items = line.strip().split('\t')
                ln_mark = '{}_{}'.format(items[0], items[1])
                m = re.search(r'DP=([\d]*)', items[7])
                if m:
                    info_dp = m.group(1)
                else:
                    info_dp = '0'
                freqs = []
                qual = items[5]
                try:
                    index_ad = items[8].split(':').index('AD')
                    index_dp = items[8].split(':').index('DP')
                except:
                    raise Exception('ERROR: check line {}'.format(n + 1))
                for s in items[9:]:
                    if s == './.':
                        freqs.append('./.')
                    else:
                        fmat = s.split(':')
                        # problem needed to be processed
                        freqs.append('{}/{}'.format(fmat[index_ad].split(',')[-1], fmat[index_dp]))
                vcf_anno_dict[ln_mark] = [info_dp, freqs, qual, []]

    print 'INFO: start reading {}'.format(exonic_variant_function)
    with open(exonic_variant_function, 'r') as evf:
        for line in evf:
            if line[0] != '#':
                line = line.split('\t')
                ln_mark = '{}_{}'.format(line[3], line[4])
                mut_type = line[1]
                mut_info = line[2]
                if ln_mark in vcf_anno_dict:
                    vcf_anno_dict[ln_mark][-1] = [mut_type, mut_info]

    print 'INFO: start reading {}'.format(variant_function)
    with open(variant_function) as vf, open(xls_out, 'w') as xls:
        xls.write('\t'.join(header_items) + '\t')
        xls.write('\t'.join(samples) + '\n')
        for n, line in enumerate(vf):
            if line[0] != '#':
                items = line.split('\t')
                ln_mark = '{}_{}'.format(items[2], items[3])
                if ln_mark in vcf_anno_dict:
                    vcf_anno = vcf_anno_dict[ln_mark]
                    qual = vcf_anno[-2]
                    info_dp = vcf_anno[0]
                    if vcf_anno[-1]:
                        mut_type = vcf_anno[-1][0]
                        mut_info = vcf_anno[-1][1]
                    else:
                        mut_type = '.'
                        mut_info = '.'
                    xls.write('{}\t{}\t{}\n'.format(
                        '\t'.join([items[2], items[3], items[4], items[5], items[6], qual, items[0], items[1], info_dp]),
                        '\t'.join([mut_type, mut_info]),
                        '\t'.join(vcf_anno[1])
                    ))

    if os.path.getsize(xls_out):
        print 'INFO: succeed in exporting {}'.format(xls_out)

if __name__ == '__main__':
    if opts.input and opts.evf and opts.vf and opts.output:
        main(opts.input, opts.evf, opts.vf, opts.output)
    else:
        parser.print_help()