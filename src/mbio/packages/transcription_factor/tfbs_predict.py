# coding=utf-8
import os
import argparse
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
# import json
# import re
__author__ = "gdq"


def parse_args():
    parser = argparse.ArgumentParser(description="""
    A simple Wrapper for MEMEM/fimo. 
    ArgumentSpecification:
    'gtf, genome, gtfParser(-450_Tss_+50), bedtools, genes(optional), backward, forward' are used together to get promoter fasta.
    'seqdb' is an alternative for input promoter fasta directly. promoter fasta will be used by FIMO.
    'thresh, qv_thresh, motifs' are another three args for FIMO.
    'corr*' are used for correlation calculation between genes.
    """)
    parser.add_argument('-motifs', required=True, help="motif file")
    parser.add_argument('-gtf', default=None, help="genome annotation file in gtf format")
    parser.add_argument('-genome', default=None, help="genome fasta file")
    parser.add_argument('-gtfParser', default="/mnt/ilustre/users/sanger-dev/sg-users/deqing/mbio/packages/transcription_factor/parseGTF.pl",
                        help="A script to get promoter bed file based on gtf file, result file is 'promoter.bed'")
    parser.add_argument('-bedtools', default="/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bedtools-2.25.0/bin/bedtools",
                        help="The bedtools to get promoter fasta based on promoter.bed")
    parser.add_argument('-backward', default=450, type=int, help="backward distance from Tss, used for promoter extract ")
    parser.add_argument('-forward', default=50, type=int, help="forward distance from Tss, used for promoter extract ")
    parser.add_argument('-genes', default=None, help="A file with each row contains one gene of genome, these genes' promoter fasta will be extracted.")
    parser.add_argument('-fimo', default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/meme_4.12.0/bin/fimo",
                        help="The MEME/fimo to scan motif in promoter fasta")
    parser.add_argument('-seqdb', default=None, help="A file containing a collection of sequences in FASTA format. "
                                                     "If used, promoter fasta from genome will not be used")
    parser.add_argument('-thresh', type=float, default=0.0001, help=""" 
     The output threshold for displaying search results.
     Only search results with a p-value less than the threshold will be output. 
     The threshold can be set to use q-values rather than p-values via the -qv_thresh option.""")
    parser.add_argument('-qv_thresh', type=int, default=0, help="if not 0, use q-value threshold")
    parser.add_argument('-exp_matrix', default=None, help="expression matrix")
    parser.add_argument('-geneid2tfid', default=None, help="Header required. File with two columns, tfid(used in motif file) -> gene_id(sed in exp_matrix)")
    parser.add_argument('-corr_cutoff', default=1, type=float, help="absolute correlation cutoff")
    parser.add_argument('-corr_pvalue', default=0, type=float, help="correlation pvalue/padjust cutoff")
    parser.add_argument('-corr_use_padjust', default=0, type=int, help="if use padjust correlation cutoff")
    args = parser.parse_args()
    genome_info = args.gtf is None or args.genome is None
    seqdb_info = args.seqdb is None
    if genome_info and seqdb_info:
        raise Exception("seqdb and Genome info cannot be both empty!")
    if (not genome_info) and (not seqdb_info):
        raise Exception("seqdb and Genome info are both provide, but which to use for motif scan?")
    return args


def generate_promoter_fasta(args):
    cmd = "perl {} {} tss -gid -backward {} -forward {} ".format(args.gtfParser, args.gtf, args.backward, args.forward)
    print("get promoter.bed: ", cmd)
    subprocess.check_call(cmd, shell=True)
    if args.genes is not None:
        with open(args.genes) as f:
            genes = [x.strip().split()[0] for x in f]
        with open('new_promoter.bed','w') as fw, open("promoter.bed") as fr:
            all_genes = list()
            for line in fr:
                tmp_gene = line.strip().split()[3]
                all_genes.append(tmp_gene)
                if tmp_gene in genes:
                    fw.write(line)
        diff_genes = set(genes) - set(all_genes)
        if len(diff_genes) == len(set(genes)):
            raise Exception("None of specified genes are found in genome, please check input!")
        if len(diff_genes) > 1:
            print(diff_genes, 'they are not found in genome!')
        os.rename("promoter.bed", "all_promoter.bed")
        os.rename("new_promoter.bed", "promoter.bed")
    cmd2 = "{} getfasta -name -s -fi {} -bed promoter.bed -fo promoter.fa".format(args.bedtools, args.genome)
    print("get promoter fasta: ", cmd2)
    if os.path.exists("promoter.fa"):
        pass
    else:
        subprocess.check_call(cmd2, shell=True)


def multtest_correct(p_values, method=3):
    """
    1. Bonferroni. ---> bonferroni
    2. Bonferroni Step-down(Holm) ---> Holm
    3. Benjamini and Hochberg False Discovery Rate ---> BH
    :param pvalue_list:
    :param method:
    :return: np.array
    """
    pvalue_list = list(p_values)
    n = len(pvalue_list)
    if method == 1:
        fdr = [eachP * n for eachP in pvalue_list]
    elif method == 2:
        sorted_pvalues = sorted(pvalue_list)
        fdr = [eachP * (n - sorted_pvalues.index(eachP)-1) for eachP in pvalue_list]
    elif method == 3:
        sorted_pvalues = sorted(pvalue_list)
        fdr = [eachP * n / (sorted_pvalues.index(eachP) + 1) for eachP in pvalue_list]
    else:
        raise Exception("the method is not supported")
    fdr = np.array(fdr)
    fdr[fdr > 1] = 1.
    return fdr


def run_fimo(args):
    cmd = "{} ".format(args.fimo)
    cmd += "--bgfile --motif-- "  # Use '--bgfile --motif--' to read the background from the motif file.
    cmd += "--max-strand "
    cmd += "--thresh {} ".format(args.thresh)
    if args.qv_thresh != 0:
        cmd += "--qv-thresh ".format()
    cmd += "{} ".format(args.motifs)
    if args.seqdb is None:
        generate_promoter_fasta(args)
        cmd += 'promoter.fa '
    else:
        cmd += args.seqdb
    print("run fimo: ", cmd)
    if os.path.exists("fimo_out/fimo.txt"):
        pass
    else:
        subprocess.check_call(cmd, shell=True)

    # modify result
    fimo_pd = pd.read_table("fimo_out/fimo.txt", header=0)
    fimo_pd.columns = [u'tf_id', u'motif_id', u'target_id', u'start', u'stop',
                       u'strand', u'score', u'p-value', u'q-value', u'matched_sequence']
    if args.geneid2tfid:
        with open(args.geneid2tfid) as f:
            _ = f.readline()
            tf_id2gene_id = dict()
            for line in f:
                if not line.strip():
                    continue
                if len(line.strip().split()) < 2:
                    continue
                k, v = line.strip().split()[0:2]
                tf_id2gene_id.setdefault(v, set())
                tf_id2gene_id[v].add(k)
                if v[-2] == '.' or v[-3] == '.':
                    v_else = v[:v.rfind('.')]
                    if v_else in tf_id2gene_id:
                        tf_id2gene_id[v_else].add(k)
                    else:
                        tf_id2gene_id[v_else] = set(k)
    else:
        # assume tf_id is also gene_id
        tf_id2gene_id = dict()
    new_fimo_result = list()
    if tf_id2gene_id:
        for _, each_row in fimo_pd.iterrows():
            each_tf = each_row['tf_id']
            if each_tf in tf_id2gene_id:
                for each_gene in tf_id2gene_id[each_tf]:
                    tmp_row = each_row.copy()
                    tmp_row['tf_geneid'] = each_gene
                    new_fimo_result.append(tmp_row)
            else:
                print('{} not found in {}. Thus itself will be used as gene_id'.format(each_tf, args.geneid2tfid))
                each_row['tf_geneid'] = each_tf
                new_fimo_result.append(each_row)
        fimo_pd = pd.DataFrame(new_fimo_result, columns=list(fimo_pd.columns) + ['tf_geneid'])
    else:
        fimo_pd['tf_geneid'] = fimo_pd['tf_id']

    # calculate corr
    if args.exp_matrix is not None:
        source = fimo_pd['tf_geneid']
        target = fimo_pd['target_id']
        exp_matrix_pd = pd.read_table(args.exp_matrix, header=0, index_col=0)
        corr_value_list = list()
        corr_pvalue_list = list()
        for s, t in zip(source, target):
            try:
                # s 和 t可能不在表达矩阵中，因此无法提取表达信息
                corr_value, corr_pvalue = pearsonr(exp_matrix_pd.loc[s, :], exp_matrix_pd.loc[t, :])
            except Exception as e:
                print("Warning", e)
                corr_value, corr_pvalue = 0, 1
            corr_pvalue_list.append(corr_pvalue)
            corr_value_list.append(corr_value)
        corr_fdr_list = multtest_correct(corr_pvalue_list)
        fimo_pd['corr'] = corr_value_list
        fimo_pd['corr_pvalue'] = corr_pvalue_list
        fimo_pd['corr_padjust'] = corr_fdr_list

        # filter
        target_rows = list()
        num = -1
        # print fimo_pd
        for corr, pval, fdr in zip(corr_value_list, corr_pvalue_list, corr_fdr_list):
            num += 1
            if abs(corr) >= args.corr_cutoff:
                if not args.corr_use_padjust:
                    if pval <= args.corr_pvalue:
                        target_rows.append(num)
                else:
                    if fdr <= args.corr_pvalue:
                        target_rows.append(num)

        # print target_rows
        fimo_pd = fimo_pd.iloc[target_rows, :]

    fimo_pd=fimo_pd.drop_duplicates(subset=None, keep='first', inplace=False)
    fimo_pd.to_csv("tfbs_predict.xls", header=True, index=False, sep='\t')
    # stat
    geneid2tfid = dict(zip(fimo_pd['tf_geneid'], fimo_pd['tf_id']))
    target_stat = pd.DataFrame(fimo_pd['tf_geneid'].value_counts())
    target_stat.index.name = 'tf_geneid'
    target_stat.columns = ['target_number']
    tmp_list = list()
    for each in target_stat.index:
        if each in geneid2tfid:
            tmp_list.append(geneid2tfid[each])
        else:
            tmp_list.append(each)
    target_stat['tf_id'] = tmp_list
    target_stat.to_csv("tf_target_number.xls", header=True, index=True, sep='\t')


if __name__ == '__main__':
    run_fimo(parse_args())

