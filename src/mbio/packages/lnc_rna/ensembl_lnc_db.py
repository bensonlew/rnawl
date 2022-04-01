# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import re
import os
import subprocess

parser = OptionParser(description='Export lncRNA database from ensembl reference')
parser.add_option('-g', '--gtf', dest='gtf', help='input ENSEMBL GTF file')
parser.add_option('-t', '--type', dest='type', help='transcript biotypes of lncRNA, separated by comma')
parser.add_option('-f', '--fasta', dest='fasta', help='input ENSEMBL FASTA file')
parser.add_option('-b', '--gffread', dest='gffread', help='program path of gffread')
parser.add_option('-o', '--output', dest='output', help='output directory')
(opts, args) = parser.parse_args()

def main(gtf_in, biotypes, gffread, fasta_in, dir_out):
    print 'INFO: start searching lncRNA in {}'.format(gtf_in)
    lncrna_ids, gene_ids = search_lncrna(gtf_in, biotypes.split(','))
    if lncrna_ids:
        print 'INFO: succeed in finding {} lncRNAs'.format(len(lncrna_ids))
        gtf_out = os.path.join(dir_out, 'lncrna.gtf')
        ret = select_lncrna(gtf_in, lncrna_ids, gene_ids, gtf_out)
        if ret and os.path.getsize(gtf_out):
            print 'INFO: succeed in exporting {}'.format(gtf_out)
            fasta_out = os.path.join(dir_out, 'lncrna.fa')
            ret = export_fasta(gffread, gtf_in, fasta_in, fasta_out)
            if ret == True and os.path.getsize(fasta_out):
                print 'INFO: succeed in exporting {}'.format(fasta_out)
                ids_map = os.path.join(dir_out, 'ids_map.tsv')
                export_ids_map(gtf_out, ids_map)
                if os.path.getsize(ids_map):
                    print 'INFO: succeed in exporting {}'.format(ids_map)
            else:
                print 'WARN: {}'.format(ret)

def search_lncrna(gtf_in, lncrna_biotypes):
    lncrna_ids = set()
    gene_ids = set()
    for line in open(gtf_in):
        if line[0] != '#' and line.strip():
            items = line.strip().split('\t')
            if len(items) == 9:
                m = re.match(r'.*gene_id "(\S+)";.*transcript_id "(\S+)";.*transcript_biotype "(\S+)";', items[8])
                if m:
                    gene_id, transcript_id, transcript_biotype = m.groups()
                    if transcript_biotype in lncrna_biotypes:
                        lncrna_ids.add(transcript_id)
                        gene_ids.add(gene_id)
    else:
        return lncrna_ids, gene_ids

def select_lncrna(gtf_in, lncrna_ids, gene_ids, gtf_out):
    lines = list()
    for line in open(gtf_in):
        if line[0] != '#' and line.strip():
            items = line.strip().split('\t')
            if len(items) == 9:
                if 'transcript_id' in items[8] and 'gene_id' in items[8]:
                    m = re.match(r'.*transcript_id "(\S+)";', items[8])
                    if m and m.group(1) in lncrna_ids:
                        lines.append(line)
                elif 'transcript_id' not in items[8] and 'gene_id' in items[8]:
                    m = re.match(r'.*gene_id "(\S+)";', items[8])
                    if m and m.group(1) in gene_ids:
                        lines.append(line)
    else:
        open(gtf_out, 'w').writelines(lines)
        return True

def export_fasta(gffread, gtf_in, fasta_in, fasta_out):
    cmd = '{} {}'.format(gffread, gtf_in)
    cmd += ' -g {}'.format(fasta_in)
    cmd += ' -w {}'.format(fasta_out)
    spc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = spc.wait()
    if ret:
        return cmd
    else:
        return True

def export_ids_map(gtf_out, ids_map):
    header = 'transcript_id\tensembl_transcript_id\tgene_id\tgene_name\n'
    lines = set()
    for line in open(gtf_out):
        attributes = line.split('\t')[8]
        m = re.match(r'.*gene_id "(\S+)";.*transcript_id "(\S+)";.*gene_name "(\S+)";', attributes)
        if m:
            gene_id, transcript_id, gene_name = m.groups()
            lines.add('{}\t{}\t{}\t{}\n'.format(transcript_id, transcript_id, gene_id, gene_name))
    else:
        open(ids_map, 'w').writelines([header] + list(lines))

if __name__ == '__main__':
    if opts.gtf and opts.type and opts.gffread and opts.fasta and opts.output:
        main(opts.gtf, opts.type, opts.gffread, opts.fasta, opts.output)
    else:
        parser.print_help()
