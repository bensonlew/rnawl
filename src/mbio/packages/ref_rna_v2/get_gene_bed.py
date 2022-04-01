# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import re
import os

parser = OptionParser(description='Generate gene BED file')
parser.add_option('--gtf', dest='gtf', type=str, help='input GTF file')
parser.add_option('--bed', dest='bed', type=str, help='input transcript BED file')
parser.add_option('--output', dest='output', type=str, help='output directory')
(opts, args) = parser.parse_args()

def main(opts):
    print 'INFO: start reading {}'.format(opts.gtf)
    g2t_dict, trans_pos = get_g2t_dict(opts.gtf)
    print 'INFO: start reading {}'.format(opts.bed)
    t2b_dict = get_t2b_dict(opts.bed)
    print 'INFO: start building relation dict between gene and bed'
    g2b_dict = get_g2b_dict(g2t_dict, t2b_dict, trans_pos)
    export_bed(g2b_dict, opts.output)
    if os.path.getsize(opts.output):
        print 'INFO: succeed in exporting {}'.format(opts.output)

def get_g2t_dict(gtf):
    dct = dict()
    trans_pos =dict()
    for line in open(gtf):
        if line[0] != '#':
            items = line.strip().split('\t')
            if len(items) >= 8 and 'gene_id' in items[8] and 'transcript_id' in items[8]:
                mg = re.search(r'gene_id "(\S+)";', items[8])
                mt = re.search(r'transcript_id "(\S+)";', items[8])
                if mg and mt:
                    g = mg.group(1)
                    t = mt.group(1)
                    if g in dct:
                        dct[g].append(t)
                    else:
                        dct[g] = [t]

                    if items[2] == "transcript":
                        trans_pos[t] = [items[0], items[3], items[4], items[5], items[6], items[7]]
    else:
        return dct, trans_pos

def get_t2b_dict(bed):
    dct = dict()
    for line in open(bed):
        items = line.strip().split('\t')
        if len(items) >= 6:
            dct[items[3]] = items[:6]
    else:
        return dct

def get_g2b_dict(g2t, t2b, trans_pos):
    dct = dict()
    for g, ts in g2t.items():
        for t in ts:
            if g in dct:
                if int(t2b[t][1]) < dct[g]['chromStart']:
                    dct[g]['chromStart'] = int(t2b[t][1])
                if int(t2b[t][2]) > dct[g]['chromEnd']:
                    dct[g]['chromEnd'] = int(t2b[t][2])
            else:
                if t in t2b:
                    dct[g] = {'chrom': t2b[t][0], 'chromStart': int(t2b[t][1]), 'chromEnd': int(t2b[t][2]),
                              'name': g, 'score': t2b[t][4], 'strand': t2b[t][5]}
                elif t in trans_pos:
                    dct[g] = {'chrom': trans_pos[t][0], 'chromStart': int(trans_pos[t][1]), 'chromEnd': int(trans_pos[t][2]),
                              'name': g, 'score': trans_pos[t][3], 'strand': trans_pos[t][5]}
                else:
                    raise Exception("找不到基因的坐标{}".format(g))

    else:
        return dct

def export_bed(g2b, bed):
    fields = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    open(opts.output, 'w').writelines('{}\n'.format('\t'.join(str(v[k]) for k in fields)) for v in g2b.values())

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['gtf', 'bed', 'output'])):
        main(opts)
    else:
        parser.print_help()
