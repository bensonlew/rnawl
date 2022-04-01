# -*- coding:utf-8 -*-
from __future__ import print_function
from collections import defaultdict
from scipy.stats import hypergeom
import numpy as np
import argparse
import glob
import re
from textwrap import wrap
import matplotlib as mpl
import matplotlib.pyplot as plt
#plt.style.use('ggplot')
__author__ = "gdq"

def parse_gene_annot(map_file, header=True):
    """
    parse gene_classification file such as:
    ----------------------
    GOLocus Tag Accession
    PA0001  GO:0006260
    PA0001  GO:0003677
    ...
    PA0009  GO:0003688
    **Note, also support line likes: "PA0001  GO:0006260;GO:0003677"
    ----------------------
    """
    gene_class = defaultdict(set)
    class_gene = defaultdict(set)
    d_r = defaultdict(set)
    f = open(map_file)
    if header:
        head = f.readline()
    for line in f:
        if line.strip():
            if len(line.strip().split())<=1:
                continue
            a, b = line.strip().split()
            if ";" in b:
                cls = b.split(";")
                gene_class[a].update(cls)
                for each_class in cls:
                    class_gene[each_class].add(a)
            else:
                gene_class[a].add(b)
                class_gene[b].add(a)
    f.close()
    return class_gene, gene_class, len(gene_class.keys())


def read_diff_genes(deg_file):
    """ get diff gene list"""
    with open(deg_file) as f:
        deg_list = [x.strip().split()[0] for x in f if x.strip()]
    return set(deg_list)


def prepare_hypergeom_data(class_gene_dict, gene_class_dict, deg_set, total_gene_number):
    pop_number = total_gene_number
    study_number = len(deg_set)
    # get all DE gene associated classification, named considered_classes 
    considered_classes = set()
    for gene in deg_set:
        if gene not in gene_class_dict:
            pass
        else:
            considered_classes.update(gene_class_dict[gene])
    # print(considered_classes)

    for each_class in considered_classes:
        associated_genes = class_gene_dict[each_class]
        pop_hitnumber = len(associated_genes)
        associated_diff_genes = associated_genes.intersection(deg_set)
        associated_diff_info = '|'.join(list(associated_diff_genes))
        study_hitnumber = len(associated_diff_genes)
        yield study_hitnumber, pop_number, pop_hitnumber, study_number, each_class, associated_diff_info


def multtest_correct(p_values, methods=3):
    """
    1. Bonferroni
    2. Bonferroni Step-down(Holm)
    3. Benjamini and Hochberg False Discovery Rate
    4. FDR Benjamini-Yekutieli
    :param pvalue_list:
    :param methods:
    :return: np.array
    """
    def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):
        '''pvalue correction for false discovery rate
        This covers Benjamini/Hochberg for independent or positively correlated and
        Benjamini/Yekutieli for general or negatively correlated tests. Both are
        available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
        Parameters
        ----------
        pvals : array_like
            set of p-values of the individual tests.
        alpha : float
            error rate
        method : {'indep', 'negcorr')
        Returns
        -------
        rejected : array, bool
            True if a hypothesis is rejected, False if not
        pvalue-corrected : array
            pvalues adjusted for multiple hypothesis testing to limit FDR
        Notes
        -----
        If there is prior information on the fraction of true hypothesis, then alpha
        should be set to alpha * m/m_0 where m is the number of tests,
        given by the p-values, and m_0 is an estimate of the true hypothesis.
        (see Benjamini, Krieger and Yekuteli)
        The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
        of false hypotheses will be available (soon).
        Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
        fdr_by.
        '''
        def _ecdf(x):
            '''
            no frills empirical cdf used in fdrcorrection
            '''
            nobs = len(x)
            return np.arange(1, nobs+1)/float(nobs)

        pvals = np.asarray(pvals)
        if not is_sorted:
            pvals_sortind = np.argsort(pvals)
            pvals_sorted = np.take(pvals, pvals_sortind)
        else:
            pvals_sorted = pvals  # alias

        if method in ['i', 'indep', 'p', 'poscorr']:
            ecdffactor = _ecdf(pvals_sorted)
        elif method in ['n', 'negcorr']:
            cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
            ecdffactor = _ecdf(pvals_sorted) / cm
        else:
            raise ValueError('only indep and negcorr implemented')
        reject = pvals_sorted <= ecdffactor*alpha
        if reject.any():
            rejectmax = max(np.nonzero(reject)[0])
            reject[:rejectmax] = True

        pvals_corrected_raw = pvals_sorted / ecdffactor
        pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
        del pvals_corrected_raw
        pvals_corrected[pvals_corrected>1] = 1
        if not is_sorted:
            pvals_corrected_ = np.empty_like(pvals_corrected)
            pvals_corrected_[pvals_sortind] = pvals_corrected
            del pvals_corrected
            reject_ = np.empty_like(reject)
            reject_[pvals_sortind] = reject
            return reject_, pvals_corrected_
        else:
            return reject, pvals_corrected

    pvalue_list = list(p_values)
    n = len(pvalue_list)
    if methods == 1:
        fdr = [eachP*n for eachP in pvalue_list]
    elif methods == 2:
        sorted_pvalues = sorted(pvalue_list)
        fdr = [eachP*(n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
    elif methods == 3:
        sorted_pvalues = sorted(pvalue_list)
        fdr = [eachP*n/(sorted_pvalues.index(eachP)+1) for eachP in pvalue_list]
    elif methods == 4:
        _, fdr = fdrcorrection(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
    fdr = np.array(fdr)
    fdr[fdr > 1] = 1.
    return fdr


def hypergeom_test(data, sort_fdr=True, correct_method=3):
    p_value_list = []
    ratio_in_study_list = []
    ratio_in_pop_list = []
    classes = []
    hit_genes = []
    hit_links = []
    study_number_list = []
    for study_hitnumber, pop_number, pop_hitnumber, study_number, each_class, associated_diff_info in data:
        study_number_list.append(study_hitnumber)
        # p_value = 1 - hypergeom.cdf(study_hitnumber-1,pop_number, pop_hitnumber, study_number)
        #20191121 modify by fwy 改为sf（文档介绍精读或许更高）
        p_value = hypergeom.sf(study_hitnumber - 1, pop_number, pop_hitnumber, study_number)
        ratio_in_study = str(study_hitnumber)+'/'+str(study_number)
        ratio_in_pop = str(pop_hitnumber)+'/'+str(pop_number)
        p_value_list.append(p_value)
        ratio_in_study_list.append(ratio_in_study)
        ratio_in_pop_list.append(ratio_in_pop)
        classes.append(each_class)
        hit_genes.append(associated_diff_info)


    q_value_list = multtest_correct(p_value_list, methods=correct_method)
    number = len(q_value_list)
    #print(number)

    result = zip(study_number_list, classes,ratio_in_study_list,ratio_in_pop_list, p_value_list, q_value_list, hit_genes)
    if sort_fdr:
        sorted_result = sorted(result, key=lambda x: (x[5],x[4]))
    else:
        sorted_result = sorted(result, key=lambda x: (x[4],x[5]))

    return sorted_result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-deg', required=True, help='file with two columns: gene\tup/down')
    parser.add_argument('-g2p', required=True, help='file with two columns: gene\tpath:konumber')

    parser.add_argument('--FDR', default=False, action='store_true', help='if used, FDR will be used for plotting')
    parser.add_argument('-correct', metavar='correct_method', default=3, type=int, help='multi-test method')
    args = parser.parse_args()

    g2p_file = args.g2p
    deg_file = args.deg

    deg_files = [deg_file]
    # calculating result
    p_gene, gene_p, gene_number = parse_gene_annot(g2p_file)

    for deg_file in deg_files:
        deg_set = read_diff_genes(deg_file)
        deg_set = deg_set.intersection(set(gene_p.keys()))
        data = prepare_hypergeom_data(p_gene, gene_p, deg_set, gene_number)
        result = hypergeom_test(data, sort_fdr=args.FDR, correct_method=args.correct)
        f = open(deg_file + '.do_enrichment.xls', 'w')
        # print(deg_file + '.kegg_enrichment.xls')
        # zip(study_number_list, classes,ratio_in_study_list,ratio_in_pop_list, p_value_list, q_value_list, hit_genes)
        f.write('#Study_num\tDO_ID\tRatio_in_study\tRatio_in_pop\tPvalue\tPadjust\tGenes\n')
        for each in result:
            f.write('\t'.join([str(x) for x in each])+'\n')
        f.close()
