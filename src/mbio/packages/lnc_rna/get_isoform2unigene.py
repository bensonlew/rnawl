# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os
import re

parser = OptionParser(description='Generate isoform2unigene file from input GTF file')
parser.add_option('-i', '--gtf', dest='gtf', help='Input GTF file')
parser.add_option('-r', '--g2t2p', dest='g2t2p', help='Input g2t2p file')
parser.add_option('-o', '--i2u', dest='i2u', help='Output isoform2unigene file')
(opts, args) = parser.parse_args()

def main(gtf_in, g2t2p_in, i2u_out):
    if os.path.isfile(gtf_in):
        print 'INFO: start reading {}'.format(gtf_in)
        if os.path.isfile(g2t2p_in):
            print 'INFO: start reading {}'.format(g2t2p_in)
            gtf2isoform2unigene(gtf_in, g2t2p_in, i2u_out)
        else:
            print 'INFO: start processing without g2t2p file'
            gtf2isoform2unigene(gtf_in, None, i2u_out)
        if os.path.getsize(i2u_out) > 0:
            print 'INFO: succeed in exporting {}'.format(i2u_out)
    else:
        raise Exception('ERROR: fail to read {}'.format(gtf_in))

def gtf2isoform2unigene(gtf_in, g2t2p_in, i2u_out):
    print 'INFO: start building transcript to gene dict'
    t2g = get_transcript_id_gene_id_dict(gtf_in)
    print 'INFO: start building transcript to exon length dict'
    t2el = get_transcript_id_exon_len_dict(gtf_in)
    print 'INFO: start building transcript to cds length dict'
    t2cl = get_transcript_id_cds_len_dict(gtf_in)
    print 'INFO: start mapping keys of three built dict'
    for t in t2g:
        if t2el.has_key(t):
            pass
        elif t in t2cl:
            t2el[t] = t2cl[t]
        else:
            t2el[t] = 0
    print 'INFO: start sorting transcript to gene to length dict'
    t2g2l = [(t, t2g[t], t2el[t]) for t in t2g]
    t2g2l_sort = sorted(t2g2l, key=lambda x: x[2], reverse=True)
    if g2t2p_in:
        print 'INFO: start building transcript to protein dict'
        t2p = dict([tuple(line.strip().split('\t')[1:]) for line in open(g2t2p_in) if len(line.strip().split('\t')) >= 3])
    print 'INFO: start iterating for writing line'
    with open(i2u_out, 'wb') as f:
        gene_list = list()
        for t, g, l in t2g2l_sort:
            if g2t2p_in:
                if t2p.has_key(t):
                    p = t2p[t]
                else:
                    p = str()
            else:
                p = t
            if g.startswith('MSTRG') or g.startswith('XLOC') or g.startswith('XLOC'):
                if g in gene_list:
                    f.write('{}\t{}\t{}\t{}\t{}\n'.format(t, g, 'no', l, p))
                else:
                    f.write('{}\t{}\t{}\t{}\t{}\n'.format(t, g, 'yes', l, p))
                    gene_list.append(g)
            elif not (t.startswith('MSTRG') or t.startswith('XLOC') or t.startswith('TCONS')):
                if g in gene_list:
                    f.write('{}\t{}\t{}\t{}\t{}\n'.format(t, g, 'no', l, p))
                else:
                    f.write('{}\t{}\t{}\t{}\t{}\n'.format(t, g, 'yes', l, p))
                    gene_list.append(g)
            else:
                f.write('{}\t{}\t{}\t{}\t{}\n'.format(t, g, 'no', l, p))

def get_transcript_id_gene_id_dict(gtf):
    dct = dict()
    for line in open(gtf):
        if line[0] != '#':
            mg = re.match(r'.*gene_id\s"(\S+)";.*', line)
            mt = re.match(r'.*transcript_id\s"(\S+)";.*', line)
            if mg and mt:
                dct[mt.group(1)] = mg.group(1)
    return dct

def get_transcript_id_exon_len_dict(gtf):
    dct = dict()
    for line in open(gtf):
        if line[0] != '#':
            items = line.strip().split('\t')
            if len(items) > 8 and items[2] == 'exon':
                mt = re.match(r'.*transcript_id\s"(\S+)";.*', items[8])
                if mt:
                    transcript_id = mt.group(1)
                    exon_len = abs(int(items[4]) - int(items[3])) + 1
                    if dct.has_key(transcript_id):
                        dct[transcript_id] += exon_len
                    else:
                        dct[transcript_id] = exon_len
    return dct

def get_transcript_id_cds_len_dict(gtf):
    dct = dict()
    for line in open(gtf):
        if line[0] != '#':
            items = line.strip().split('\t')
            if len(items) > 8 and items[2] == 'cds':
                mt = re.match(r'.*transcript_id\s"(\S+)";.*', items[8])
                if mt:
                    transcript_id = mt.group(1)
                    cds_len = abs(int(items[4]) - int(items[3])) + 1
                    if dct.has_key(transcript_id):
                        dct[transcript_id] += cds_len
                    else:
                        dct[transcript_id] = cds_len
    return dct

if __name__ == '__main__':
    if opts.gtf and opts.i2u:
        if opts.g2t2p:
            main(opts.gtf, opts.g2t2p, opts.i2u)
        else:
            main(opts.gtf, str(), opts.i2u)
    else:
        parser.print_help()