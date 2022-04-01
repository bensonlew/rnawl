# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng,liubinxu'

import os
import pickle
import re
import subprocess
from optparse import OptionParser

import numpy as np
from Bio import SeqIO

parser = OptionParser(description='Generate relation map files from specified input')
parser.add_option('--gffread', dest='gffread', help='program path of gffread')
parser.add_option('--gtf', dest='gtf', type=str, help='input GTF file')
parser.add_option('--fasta', dest='fasta', type=str, help='input reference genome FASTA file')
parser.add_option('--g2t2p', dest='g2t2p', type=str, help='input G2T2P file')
parser.add_option('--output', dest='output', type=str, help='output directory containing files')
(opts, args) = parser.parse_args()


def main(gffread, gtf, fasta, g2t2p, dir_out):
    print 'INFO: start building T2G dict'
    t2g_dict = build_t2g(gtf)
    print 'INFO: start building T2L dict'
    t2l_dict = build_t2l(gffread, gtf, fasta, dir_out)
    print 'INFO: start building T2PL dict'
    t2pl_dict = build_t2pl(gffread, gtf, fasta, dir_out)
    print 'INFO: start building T2R dict'
    t2r_dict = build_t2r(g2t2p, t2g_dict, t2l_dict)
    print 'INFO: start building T2P dict'
    t2p_dict = build_t2p(g2t2p)
    pickle.dump(t2p_dict, open(os.path.join(dir_out, 't2p.pkl'), 'w'))
    print 'INFO: start merging dict into dataframe'
    print 't2g: {}, t2l: {}, t2r: {}, t2p: {}'.format(len(t2g_dict), len(t2l_dict), len(t2r_dict), len(t2p_dict))
    all_map_out = os.path.join(dir_out, 't2g2r2l2p.tsv')
    longest_out = os.path.join(dir_out, 'longest.t2g.tsv')
    export_result(all_map_out, longest_out, t2g_dict, t2r_dict, t2l_dict, t2p_dict, t2pl_dict)
    if os.path.getsize(all_map_out) > 0:
        print 'INFO: succeed in exporting {}'.format(all_map_out)
    if os.path.getsize(longest_out) > 0:
        print 'INFO: succeed in exporting {}'.format(longest_out)


def build_t2g(gtf):
    t2g_dict = dict()
    for line in open(gtf):
        if line[0] != '#' and len(line.strip().split('\t')) == 9:
            mt = re.match(r'.*transcript_id\s"(\S+)";.*', line)
            mg = re.match(r'.*gene_id\s"(\S+)";.*', line)
            if mt and mg:
                t2g_dict[mt.group(1)] = mg.group(1)
    else:
        return t2g_dict


def build_t2l(gffread, gtf, fasta, dir_out):
    txpt_fa = os.path.join(dir_out, 'transcripts.fasta')
    if run_gffread(gffread, gtf, fasta, txpt_fa):
        return dict((seq_record.id, len(seq_record)) for seq_record in SeqIO.parse(txpt_fa, 'fasta'))


def build_t2pl(gffread, gtf, fasta, dir_out):
    pep_fa = os.path.join(dir_out, 'pep.fasta')
    if run_gffread_pep(gffread, gtf, fasta, pep_fa):
        return dict((seq_record.id, len(seq_record)) for seq_record in SeqIO.parse(pep_fa, 'fasta'))


def build_t2r(g2t2p_file, t2g_dict, t2l_dict):
    t2r_dict = dict()
    ref_gene_set = set(t2g_dict.values())
    rep_gene_set = set()
    txpt_arr = np.intersect1d(t2g_dict.keys(), t2l_dict.keys())
    for i in sorted([(t, t2g_dict[t], t2l_dict[t]) for t in txpt_arr], key=lambda x: x[2], reverse=True):
        t, g = i[:2]
        if g in rep_gene_set:
            t2r_dict[t] = 'no'
        elif g in ref_gene_set:
            t2r_dict[t] = 'yes'
            rep_gene_set.add(g)
    else:
        return t2r_dict


def build_t2p(g2t2p_file):
    t2p_dict = dict()
    for line in open(g2t2p_file):
        items = line.strip().split('\t')
        t = items[1]
        if len(items) == 3:
            t2p_dict[t] = items[2]
        else:
            t2p_dict[t] = str()
    else:
        return t2p_dict


def export_result(all_map_out, longest_out, t2g_dict, t2r_dict, t2l_dict, t2p_dict, t2pl_dict):
    g2t_dict = dict()
    # 只判断包含蛋白的序列
    for trans, gene in t2g_dict.items():
        if trans in t2pl_dict:
            if gene in g2t_dict:
                if t2pl_dict[trans] > g2t_dict[gene]["len"]:
                    g2t_dict[gene] = {"r_trans": trans, "len": t2pl_dict[trans]}
                else:
                    pass
            else:
                g2t_dict[gene] = {"r_trans": trans, "len": t2pl_dict[trans]}

    pep_genes = g2t_dict.keys()
    for trans, gene in t2g_dict.items():
        if gene not in pep_genes:
            if gene in g2t_dict:
                if t2l_dict[trans] > g2t_dict[gene]["len"]:
                    g2t_dict[gene] = {"r_trans": trans, "len": t2l_dict[trans]}
                else:
                    pass
            else:
                g2t_dict[gene] = {"r_trans": trans, "len": t2l_dict[trans]}

    with open(all_map_out, 'w') as f1_o, open(longest_out, 'w') as f_o:
        for trans, gene in t2g_dict.items():
            if trans in t2p_dict:
                pep = t2p_dict[trans]
            else:
                pep = ""

            if g2t_dict[gene]["r_trans"] == trans:
                f1_o.write("\t".join([trans, gene, "yes", str(t2l_dict[trans]), pep]) + "\n")
                f_o.write("\t".join([trans, gene]) + "\n")
            else:
                f1_o.write("\t".join([trans, gene, "no", str(t2l_dict[trans]), pep]) + "\n")


'''
def export_result(all_map_out, longest_out, *t2x_dicts):
    df = pd.concat([pd.Series(t2x_dict) for t2x_dict in t2x_dicts], axis=1)
    df = df[df[1].isna().isin([False])]
    df[2] = df[2].apply(np.int)
    df = df.sort_values(by=2, ascending=False)
    df.to_csv(all_map_out, sep='\t', header=None)
    df[df[1] == 'yes'].reindex([0], axis=1).to_csv(longest_out, sep='\t', header=None)
'''


def run_gffread(gffread, gtf_in, fasta_in, fasta_out):
    cmd = '{} {}'.format(gffread, gtf_in)
    cmd += ' -g {}'.format(fasta_in)
    cmd += ' -w {}'.format(fasta_out)
    spc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = spc.wait()
    if ret:
        print 'WARN: check the following command:\n{}'.format(cmd)
        raise Exception('ERROR: fail to export {}'.format(fasta_out))
    else:
        print 'INFO: succeed in exporting {}'.format(fasta_out)
        return True


def run_gffread_pep(gffread, gtf_in, fasta_in, fasta_out):
    cmd = '{} {}'.format(gffread, gtf_in)
    cmd += ' -g {}'.format(fasta_in)
    cmd += ' -y {}'.format(fasta_out)
    spc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = spc.wait()
    if ret:
        print 'WARN: check the following command:\n{}'.format(cmd)
        raise Exception('ERROR: fail to export {}'.format(fasta_out))
    else:
        print 'INFO: succeed in exporting {}'.format(fasta_out)
        return True


if __name__ == '__main__':
    if all(map(hasattr, [opts] * 5, ['gffread', 'gtf', 'fasta', 'g2t2p', 'output'])):
        main(opts.gffread, opts.gtf, opts.fasta, opts.g2t2p, opts.output)
    else:
        parser.print_help()
