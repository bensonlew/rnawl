# -*- coding: utf-8 -*-
import os
import math
import subprocess
import numpy as np
import pandas as pd
import argparse
from pprint import pprint
from multiprocessing import Pool
import glob
import matplotlib
matplotlib.use('Agg')
from collections import defaultdict
from collections import OrderedDict
from matplotlib import pyplot as plt
import dask
# import pathos.pools as pp
import datetime
import re
import psutil
__author__ = 'gdq'
matplotlib.style.use('ggplot')

def run_script(codes):
    subprocess.call("Rscript {}".format(codes), shell=True)


class PvalueCorrect(object):
    def multtest_correct(self, p_values, method=3):
        """
        1. Bonferroni. ---> bonferroni
        2. Bonferroni Step-down(Holm) ---> Holm
        3. Benjamini and Hochberg False Discovery Rate ---> BH
        4. FDR Benjamini-Yekutieli --->BY
        :param pvalue_list:
        :param method:
        :return: np.array
        """
        pvalue_list = list(p_values)
        n = len(pvalue_list)
        if method == 1:
            fdr = [eachP*n for eachP in pvalue_list]
        elif method == 2:
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP*(n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
        elif method == 3:
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP*n/(sorted_pvalues.index(eachP)+1) for eachP in pvalue_list]
        elif method == 4:
            _, fdr = self.fdr_correction(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
        fdr = np.array(fdr)
        fdr[fdr > 1] = 1.
        return fdr

    @staticmethod
    def fdr_correction(pvals, alpha=0.05, method='indep', is_sorted=False):
        """pvalue correction for false discovery rate
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
        """
        def _ecdf(x):
            """
            no frills empirical cdf used in fdrcorrection
            """
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
        pvals_corrected[pvals_corrected > 1] = 1
        if not is_sorted:
            pvals_corrected_ = np.empty_like(pvals_corrected)
            pvals_corrected_[pvals_sortind] = pvals_corrected
            del pvals_corrected
            reject_ = np.empty_like(reject)
            reject_[pvals_sortind] = reject
            return reject_, pvals_corrected_
        else:
            return reject, pvals_corrected


class DiffExpToolbox(PvalueCorrect):
    """
    A toolbox contains several differential analysis tools for RNAseq data.
    Currently, only 3 tools supported: edgeR, DEseq2, DEGseq
    Note: only parts of the functions of each tool are implemented.
    """
    def __init__(self, count_matrix, group_info, cmp_info, batch_info=None, exp_matrix=None, exp_type='fpkm',
                 sig_type='pvalue', stat_cutoff=0.05, fc_cutoff=2, padjust_way=3, pool_size=5,filter_method=None,tpm_filter_threshold=0,
                 normalized='tmm', is_duplicate=True, prob=0.9):
        """
        initial inputs
        :param count_matrix: path of raw count table, '\t' as separator. No duplicated row header !
        :param exp_matrix: path of normalized expression value table,'\t' as separator.
                           if None, the second column of count_matrix must be gene length which
                           will be used to calculate fpkm using edgeR.
        :param group_info: path of group info, file with at least two columns. if no replicate
               exist, just use sample name as group name. Header line starts with '#'.
               --------------------
               #sample group_name  group_name
               s1   group1
               s2   group1  group3
               s3   group2
               s4   group2  group3
               s5   s5
               s6   s6
               --------------------
        :param cmp_info: path of cmp info, file with only two columns. Header line starts with '#'.
               -----------------
               #ctrl    test
               group1   group2
               group2   group3
               s5       s6
               -----------------
        :return: Results will be generated in current directory.
                 tmp/ contain raw results of diff analysis
                 *_vs*.{edgeR, deseq, ...}.diffexp.xls
        """
        self.pool_size = pool_size
        self.count = count_matrix
        self.exp = exp_matrix
        self.batch = batch_info
        self.exp_type = exp_type
        if sig_type not in ['pvalue', 'padjust']:
            raise Exception('sig_type is not pvalue/padjust')
        self.sig_type = sig_type
        self.fc_cutoff = fc_cutoff
        self.stat_cutoff = stat_cutoff
        self.padjust_way = padjust_way
        self.tpm_filter_threshold = tpm_filter_threshold
        self.count_filtered = None
        self.filtered_seqs = []
        self.normalized = normalized
        self.sample = is_duplicate
        self.prob = prob
        self.filter_method = filter_method
        # group_info -> dict, group_name as key, list of sample names as values. {group:[s1,s2,]}
        used_sample_list = list()
        sample_list = list()
        with open(group_info) as f:
            group_dict = dict()
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                tmp_list = line.strip().split()
                sample_list.append(tmp_list[0])
                for g in tmp_list[1:]:
                    group_dict.setdefault(g, list())
                    group_dict[g].append(tmp_list[0])
            for g in group_dict.keys():
                group_dict[g] = sorted(list(set(group_dict[g])))
        self.group_dict = group_dict
        self.samples = sorted(list(set(sample_list)))
        self.sample_dict = dict()
        for key,value in self.group_dict.items():
            for sample in value:
                self.sample_dict[sample] = key
        # batch_info -> dict, sample_name as key, batch as values. {s1:batch1}
        if self.batch is not None:
            with open(batch_info) as f:
                batch_dict = dict()
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    tmp_list = line.strip().split()
                    batch_dict[tmp_list[0]]=tmp_list[1]
                #     for b in tmp_list[1:]:
                #         batch_dict.setdefault(b,list())
                #         batch_dict[b].append(tmp_list[0])
                # for b in batch_dict.keys():
                #     batch_dict[b] = sorted(list(set(batch_dict[b])))
            self.batch_dict = batch_dict
        # comparison info -> list. [(ctrl, test), ...]
        with open(cmp_info) as f:
            cmp_list = list()
            error_names = list()
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                tmp_ctrl, tmp_test = line.strip().split()
                if tmp_ctrl not in self.group_dict and tmp_ctrl not in sample_list:
                    error_names.append(tmp_ctrl)
                if tmp_test not in self.group_dict and tmp_test not in sample_list:
                    error_names.append(tmp_test)
                cmp_list.append((tmp_ctrl, tmp_test))
                used_sample_list.extend([sample for sample in self.group_dict[tmp_ctrl]])
                used_sample_list.extend([sample for sample in self.group_dict[tmp_test]])
        cmp_list = sorted(list(set(cmp_list)))
        used_sample_list = sorted(list(set(used_sample_list)))
        self.use_sample_list = used_sample_list
        self.cmp_list = cmp_list
        # print sample info
        pprint('group_dict is: ')
        pprint(self.group_dict)
        pprint('comparison list is (ctrl, test): ')
        pprint(self.cmp_list)
        pprint('sample_dict is: ')
        pprint(self.sample_dict)
        if error_names:
            raise Exception('Each group name of {} is not in {}!'.format(error_names, group_info))

        # check the consistency between group info and count table
        with open(self.count) as f:
            count_samples = f.readline().strip('\n').split('\t')[1:]
        diff = set(self.samples).difference(set(count_samples))
        if diff:
            raise Exception('samples: {} are not contained in count table file'.format(diff))

        # transform count_table and exp_table to python dict
        if self.exp is None:
            with open(self.count) as f:
                count_samples = f.readline().strip('\n').split('\t')[2:]
            if not set(count_samples).difference(set(self.samples)):
                raise Exception('gene length column missing')
            self.exp_calculator_with_count(self.count, exp_type=exp_type)
            self.exp = self.count+'.{}.xls'.format(exp_type)
            self.count = str(self.count) + '.count.xls'
        df = pd.read_table(self.count, index_col=0)
        # self.count_dicts = df.to_dict('index')
        self.total_count_df = df
        df = pd.read_table(self.exp, index_col=0)
        # self.exp_dicts = df.to_dict('index')
        self.total_exp_df = df
        # if sorted(self.count_dicts.keys()) != sorted(self.exp_dicts.keys()):
        if sorted(list(self.total_count_df.index)) != sorted(list(self.total_exp_df.index)):
            raise Exception("The first id column of count table and exp table are different !")


    @staticmethod
    def exp_calculator_with_count(count_table_file, exp_type='both'):
        """
        calculate fpkm and tpm based on count table with second column containing gene length.
        :param count_table_file: example:
        -----------
        gene_id gene_length sample1 sample2
        gene1   1001    29  50
        gene2   1300    30  14
        -----------
        :param exp_type: expression type, fpkm, tpm, or 'both'. default:'both'.
        :return: rpkm_dict, tpm_dict
        """
        if exp_type not in ['fpkm', 'tpm', 'both']:
            raise Exception('exp_type should be fpkm or tpm or both')
        count_table = pd.read_table(count_table_file, index_col=0)
        columns = count_table.columns
        gene_len = count_table[columns[0]]
        if gene_len.min() < 11 or gene_len.max() > 200000:
            print('The minimum gene length and maximum gene length is abnormal!')
        rpkm_dict = dict()
        tpm_dict = dict()
        for sample in columns[1:]:
            # Divide the read counts by the length of each gene in kilobases.
            # This gives you reads per kilobase (RPK)
            rpk = count_table[sample]/gene_len
            # get rpkm/fpkm
            if exp_type == 'fpkm' or exp_type == 'both':
                total_counts = sum(count_table[sample])
                rpkm = rpk/total_counts*1000000*1000
                rpkm_dict[sample] = rpkm
            # get tpm
            if exp_type == 'tpm' or exp_type == 'both':
                norm_gene_len_total_counts = sum(rpk)
                tpm = rpk/norm_gene_len_total_counts*1000000
                tpm_dict[sample] = tpm
        # save results
        if exp_type == 'fpkm' or exp_type == 'both':
            df_rpkm = pd.DataFrame(rpkm_dict)
            df_rpkm.to_csv(count_table_file+'.fpkm.xls', sep='\t')
        if exp_type == 'tpm' or exp_type == 'both':
            df_tpm = pd.DataFrame(tpm_dict)
            df_tpm.to_csv(count_table_file+'.tpm.xls', sep='\t')
        df_count = count_table.iloc[:, 1:]
        df_count.to_csv(count_table_file+'.count.xls', sep='\t')

    def run(self, script_list):
        a=  psutil.swap_memory().total
        b= psutil.swap_memory().free/1024
        c = psutil.swap_memory().used/1024
        print("all{}".format(a))
        print("free{}".format(b))
        print("used{}".format(c))
        pool = Pool(self.pool_size)
        pool.map(run_script, script_list)
        pool.close()
        pool.join()

    def summary_stat(self, detail_dfs,out_stat):
        print(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        summary_ids_set=set()
        sig_id_dict = OrderedDict()
        cmp_deg_sum = dict()
        all_stat_dicts = dict()
        for cmp, detail_df in detail_dfs:
              ctrl, case = cmp.split('_vs_')
              sig_id_set = set(detail_df[detail_df['significant'] == 'yes'].index)
              sig_id_dict[cmp] = sig_id_set
              summary_ids_set.update(sig_id_set)
              cmp_deg_sum[cmp] = len(sig_id_set)
              all_stat_dicts[cmp] = detail_df
        summary_data = list()
        for seq_id in summary_ids_set:
            summary_dict = {'seq_id': seq_id, 'sum': 0}
            for pair, sig_id_set in sig_id_dict.items():
                 ctrl, case = pair.split('_vs_')
                 if seq_id in sig_id_set:
                     summary_dict[pair] = 'yes|{}'.format(all_stat_dicts[pair].loc[seq_id, 'regulate'])
                     summary_dict['sum'] += 1
                 else:
                     if seq_id in all_stat_dicts[pair].index:
                         summary_dict[pair] = 'no|{}'.format(all_stat_dicts[pair].loc[seq_id, 'regulate'])
                     else:
                         summary_dict[pair] = 'no|no test'
            if summary_dict['sum']:
                 summary_data.append(summary_dict)

        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df = summary_df.reindex(['seq_id'] + sig_id_dict.keys() + ['sum'], axis=1)
            stat_num = [self.total_count_df.shape[0]] + [cmp_deg_sum[i] for i in sig_id_dict.keys()] + [
                len(summary_ids_set)]
            df_stat = pd.DataFrame(columns=list(summary_df.columns))
            df_stat.loc[0] = stat_num
            summary_final = pd.concat([df_stat, summary_df], axis=0)
            summary_final.to_csv(out_stat, sep='\t', index=False)
        #add by fwy #20200820 如果出现了没有差异基因的情况
        else:
            summary_dict ={'seq_id': 0, 'sum': 0}
            for cmp, detail_df in detail_dfs:
                summary_dict[cmp] = 0
            summary_data=[summary_dict]
            summary_df = pd.DataFrame(summary_data)
            summary_df.set_index("seq_id",inplace=True)
            summary_df.to_csv(out_stat, sep='\t')
            print("buyaobuyaobuyao")

    def get_filter_counts_df(self,exp_total,used_samples):
        exp_cmp_total = pd.DataFrame(exp_total, columns=used_samples)
        if self.filter_method is None:
            exp_cmp_use = exp_cmp_total
        elif self.filter_method == "sum":
            exp_cmp_use = exp_cmp_total.loc[exp_cmp_total.sum(axis=1) > self.tpm_filter_threshold]
        elif self.filter_method == "average":
            exp_cmp_use = exp_cmp_total.loc[exp_cmp_total.sum(axis=1) / len(used_samples) > self.tpm_filter_threshold]
        elif self.filter_method == "max":
            exp_cmp_use = exp_cmp_total.loc[exp_cmp_total.max(axis=1) > self.tpm_filter_threshold]
        elif self.filter_method == "min":
            exp_cmp_use = exp_cmp_total.loc[exp_cmp_total.min(axis=1) > self.tpm_filter_threshold]
        counts_total = pd.read_table(self.count, index_col="seq_id")
        counts_total = pd.DataFrame(counts_total, columns=self.use_sample_list)
        counts_cmp_use = counts_total.loc[exp_cmp_use.index]
        return counts_cmp_use

    def DEGseq(self, sep='\t', method='MARS', threshold_kind=0, output=None):
        """
        Differential Analysis with DEGseq. Currently, Only MARS method are Supported.
        :param sep: separator of count_table
        :param method: "LRT", "CTR", "FET", "MARS", "MATR", "FC"
        :param stat_value: pvalue or qvalue cutoff
        :param fold_change: fold change cutoff
        :param threshold_kind: possible kinds are:
            • ‘1’: pValue threshold,
            • ‘2’: zScore threshold,
            • ‘3’: qValue threshold (Benjamini et al. 1995),
            • ‘4’: qValue threshold (Storey et al. 2003),
            • ‘5’: qValue threshold (Storey et al. 2003) and
              Fold-Change threshold on MA-plot are both required (can
              be used only when ‘method="MARS"’).
        :param output: output directory. If None, current directory used.
        :return: Results will be in output directory
        """

        if int(threshold_kind) not in range(6):
            raise NameError("threshold_kind must be one of [0,1,2,3,4,5]")
        #with open(count_table) as f:
        #    header = f.readline().strip('\n').split(sep)

        # using R package
        no_sig_final_use = defaultdict(list)

        script_list = list()
        if output is None:
            output = os.getcwd()
        for ctrl, test in self.cmp_list:
            used_samples = self.group_dict[ctrl] + self.group_dict[test]
            exp_total = pd.read_table(self.exp, index_col="seq_id")
            counts_cmp_use = self.get_filter_counts_df(exp_total,used_samples)
            counts_cmp_use.to_csv(os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test)), sep="\t")
            counts_table = os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test))
            no_sig_final_use[ctrl + test] = list(set(exp_total.index) - set(counts_cmp_use.index))
            with open(os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test)), "r") as f:
                header = f.readline().strip('\n').split(sep)
            # script_list = list()
            # for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)+1) for x in self.group_dict[ctrl]])
            else:
                ctrl_ind = str(header.index(test)+1)
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)+1) for x in self.group_dict[test]])
            else:
                test_ind = str(header.index(test)+1)

            script_name = os.path.join(output, 'DEGseq.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(DEGseq)\n')
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            f.write("ctrl <- readGeneExp(file='{}', geneCol=1, valCol=c({}))\n".format(
                counts_table, ctrl_ind))
            f.write('test <- readGeneExp(file="{}", geneCol=1, valCol=c({}))\n'.format(
                counts_table, test_ind))
            tmp_output = os.path.join(output, '{}_vs_{}'.format(ctrl, test)+'.degseq.tmp')
            f.write(
                'DEGexp(geneExpMatrix1=test, geneCol1=1, expCol1=c(2:{}), groupLabel1="{}", '
                'geneExpMatrix2=ctrl, geneCol2=1, expCol2=c(2:{}), groupLabel2="{}", '
                'method="{}", '
                'outputDir="{}")'
                '\n'.format(test_ind.count(',')+2, test,
                            ctrl_ind.count(',')+2, ctrl,
                            method, tmp_output))
            f.close()
        else:
            self.run(script_list)

        # format result
        cmp_result_dirs = [x for x in os.listdir(output) if x.endswith('.degseq.tmp')]
        all_stat_dicts = dict()
        all_stat_details=[]
        final_details=[]
        for each in cmp_result_dirs:
            cmp_result = output + '/' + each + '/output_score.txt'
            final_detail_df = dask.delayed(self.detail_extract_DEGseq)(cmp_result, threshold_kind,each)
            final_details.append(final_detail_df)
        else:
            out_stat = os.path.join(output, 'diff_summary_degseq.xls')
            final_details = dask.compute(*final_details)
            self.summary_stat(final_details, out_stat)

    def edgeR(self, dispersion=0.1, output=None, sep='\t'):
        if output is None:
            output = os.getcwd()

        no_sig_final_use = defaultdict(list)

        script_list = list()
        for ctrl,test in self.cmp_list:
            used_samples = self.group_dict[ctrl] + self.group_dict[test]
            exp_total = pd.read_table(self.exp, index_col="seq_id")
            counts_cmp_use = self.get_filter_counts_df(exp_total, used_samples)
            counts_cmp_use.to_csv(os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test)), sep="\t")
            counts_table = os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test))
            no_sig_final_use[ctrl + test] = list(set(exp_total.index) - set(counts_cmp_use.index))
            with open(os.path.join(output,'{}_vs_{}.counts_filted'.format(ctrl, test)),"r") as f :
                header = f.readline().strip('\n').split(sep)
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)) for x in self.group_dict[ctrl]])
                ctrl_names = ["'{}'".format(ctrl) for x in self.group_dict[ctrl]]
                if self.batch is not None:
                    try:
                        ctrl_batch=["'{}'".format(self.batch_dict[x]) for x in self.group_dict[ctrl]]
                        batch_err = False
                    except:
                        batch_err = True

            else:
                ctrl_ind = str(header.index(test))
                ctrl_names = ["'{}'".format(ctrl)]
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)) for x in self.group_dict[test]])
                test_names = ["'{}'".format(test) for x in self.group_dict[test]]
                if self.batch is not None:
                    try:
                        test_batch = ["'{}'".format(self.batch_dict[x]) for x in self.group_dict[test]]
                        batch_err = False
                    except:
                        batch_err = True
            else:
                test_ind = str(header.index(test))
                test_names = ["'{}'".format(test)]



            if self.batch is not None:
                a = [x for x in ctrl_batch if x in test_batch]
                if a:
                    pass
                else:
                    batch_err = True

            all_sample_ind = []
            all_sample_names = []
            for sample in header[1:]:
                sample_batch = str(header.index(sample))
                all_sample_ind.append(sample_batch)
            all_sample_ind = ",".join(all_sample_ind)
            all_sample_names = ["'{}'".format(self.sample_dict[x]) for x in header[1:]]
            if self.batch is not None:
                all_batch = ["'{}'".format(self.batch_dict[x]) for x in header[1:]]
                if all_batch:
                    pass
                else:
                    batch_err = True

            script_name = os.path.join(output, 'edgeR.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(limma)\n')
            f.write('library(edgeR)\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                counts_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            # f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            # f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))
            f.write('all_counts <- counts[, c({})]\n'.format(all_sample_ind))
            f.write('all_group <- c({})\n'.format(','.join(all_sample_names)))
            if self.batch is not None and batch_err == False:
                # f.write('tmp_batch <- c({})\n'.format(','.join(ctrl_batch + test_batch)))
                # f.write('tmp_batch <- as.factor(tmp_batch)\n')
                f.write('all_batch <- c({})\n'.format(','.join(all_batch)))
                f.write('all_batch <- as.factor(all_batch)\n')
            # f.write('y <- DGEList(counts=tmp_counts, group=tmp_group)\n')
            f.write('y <- DGEList(counts, group=all_group)\n')
            f.write('y <- calcNormFactors(y)\n')
            # 输出normfactors
            tmp_normFactor = os.path.join(output, "{}_vs_{}.edgeR.normFactor.xls".format(ctrl, test))
            f.write('colnames = c("sample", "group","libsize","norm.factors")\n')
            f.write('write.table(t(colnames), file="{}", sep="\\t", quote=F, row.names=F, col.names=F)\n'.format(
                tmp_normFactor))
            f.write(
                'write.table(y$samples, file="{}", sep = "\t", quote = F, row.names = T, col.names = F, append = T)\n'.format(
                    tmp_normFactor))
            ## 输出normalized结果文件
            tmp_normalize = os.path.join(output, "{}_vs_{}.edgeR.normalize.xls".format(ctrl, test))
            f.write('s_normfacotrs=vector()\n')
            f.write('for(i in 1:nrow(y$samples))\n')
            f.write('s_normfacotrs[i]=y$samples[i,3]\n')
            f.write('raw_c=y$counts\n')
            f.write('for(i in 1:ncol(raw_c))\n')
            f.write('for(j in 1:nrow(raw_c))\n')
            f.write('raw_c[j,i]=raw_c[j,i]/s_normfacotrs[i]\n')
            f.write('colnames = c("gene_id",colnames(raw_c))\n')
            f.write('write.table(t(colnames), file="{}", sep="\t", quote=F, col.names=F,row.names=F)\n'.format(
                tmp_normalize))
            f.write('write.table(raw_c, file="{}", sep="\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                tmp_normalize))

            if (',' not in ctrl_ind) and (',' not in test_ind):
                # do sample vs sample if NO replicates
                f.write('result <- exactTest(y, all_group, dispersion={})\n'.format(dispersion))
            else:
                if self.batch is not None:
                    f.write('design <- model.matrix(~0+all_group+all_batch)\n')
                else:
                    f.write('design <- model.matrix(~0+all_group)\n')
                f.write('y <- estimateDisp(y, design, robust=F)\n')
                f.write('fit <- glmQLFit(y, design, robust=F)\n')
                f.write('con = makeContrasts(all_group{}-all_group{}, levels=design)\n'.format(test, ctrl))
                f.write('result <- glmQLFTest(fit, contrast=con)\n')
            f.write('a=topTags(result, n=dim(result$counts)[1], '
                    'adjust.method="BH", sort.by="PValue")\n')
            f.write('write.table(a, "{}/{}_vs_{}.edger.tmp", sep="\\t", row.names=T, '
                    'col.names=NA, quote=FALSE)\n'.format(output, ctrl, test))
            f.close()
        else:
            self.run(script_list)
        # make final report
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.edger.tmp')]
        all_stat_dicts = dict()
        all_stat_details =[]
        final_details=[]
        for each in cmp_result_dirs:
            final_detail_df = dask.delayed(self.detail_extract)(each, "edgeR")
            final_details.append(final_detail_df)

        else:
            out_stat = os.path.join(output, 'diff_summary_edger.xls')
            final_details = dask.compute(*final_details)
            self.summary_stat(final_details, out_stat)


    def DESeq2(self, output=None, sep='\t', padjust_way=None):
        """
        :param output: output directory
        :param padjust_way: the method to pvalue correction in R
        :param sep: the separator used in  count table
        Example of using deseq2:
        cs = read.table("egdeR_input.count.xls", header=T, row.names=1, sep="\t")
        colData <- data.frame(row.names=colnames(cs), group=c('c','c','c','e','e','e'))
        dds <- DESeqDataSetFromMatrix(countData=cs, colData=colData, design= ~group)
        dds <- DESeq(dds)
        res <- results(dds, contrast<-c("group", "e", "c"))
        """


        if output is None:
            output = os.getcwd()
        counts_final_use=defaultdict(str)
        no_sig_final_use=defaultdict(list)

        script_list = list()
        for ctrl,test in self.cmp_list:
            used_samples = self.group_dict[ctrl] + self.group_dict[test]
            exp_total = pd.read_table(self.exp, index_col="seq_id")
            counts_cmp_use = self.get_filter_counts_df(exp_total, used_samples)
            counts_cmp_use.to_csv(os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test)), sep="\t")
            counts_table = os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test))
            no_sig_final_use[ctrl + test] = list(set(exp_total.index) - set(counts_cmp_use.index))
            with open(os.path.join(output,'{}_vs_{}.counts_filted'.format(ctrl, test)),"r") as f :
                header = f.readline().strip('\n').split(sep)
        #script_list = list()
        #for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)) for x in self.group_dict[ctrl]])
                ctrl_names = ["'{}'".format(ctrl) for x in self.group_dict[ctrl]]
                if self.batch is not None:
                    try:
                        ctrl_batch=["'{}'".format(self.batch_dict[x]) for x in self.group_dict[ctrl]]
                        batch_err = False
                    except:
                        batch_err = True

            else:
                ctrl_ind = str(header.index(test))
                ctrl_names = ["'{}'".format(ctrl)]
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)) for x in self.group_dict[test]])
                test_names = ["'{}'".format(test) for x in self.group_dict[test]]
                if self.batch is not None:
                    try:
                        test_batch = ["'{}'".format(self.batch_dict[x]) for x in self.group_dict[test]]
                        batch_err = False
                    except:
                        batch_err = True
            else:
                test_ind = str(header.index(test))
                test_names = ["'{}'".format(test)]
            if self.batch is not None:
                a = [x for x in ctrl_batch if x in test_batch]
                if a:
                    pass
                else:
                    batch_err = True

            all_sample_ind = []
            all_sample_names = []
            for sample in header[1:]:
                sample_batch = str(header.index(sample))
                all_sample_ind.append(sample_batch)
            all_sample_ind = ",".join(all_sample_ind)
            all_sample_names = ["'{}'".format(self.sample_dict[x]) for x in header[1:]]
            if self.batch is not None:
                all_batch = ["'{}'".format(self.batch_dict[x]) for x in header[1:]]
                if all_batch:
                    pass
                else:
                    batch_err = True

            script_name = os.path.join(output, 'DESeq2.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('suppressMessages(library(DESeq2))\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                counts_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            # f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            # f.write('tmp_counts = floor(tmp_counts+0.5)\n')
            # f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))
            f.write('all_counts <- counts[, c({})]\n'.format(all_sample_ind))
            f.write('all_counts = floor(all_counts+0.5)\n')
            f.write('all_group <- c({})\n'.format(','.join(all_sample_names)))
            if self.batch is not None and batch_err is False:
                f.write('all_batch <- c({})\n'.format(','.join(all_batch)))
                f.write('colData <- data.frame(row.names=colnames(all_counts), group=all_group, batch=all_batch)\n')
                f.write('dds <- DESeqDataSetFromMatrix(countData=all_counts, colData=colData, '
                        'design= ~group+batch)\n')
            else:
                f.write('colData <- data.frame(row.names=colnames(all_counts), group=all_group)\n')
                f.write('dds <- DESeqDataSetFromMatrix(countData=all_counts, colData=colData, '
                        'design= ~group)\n')
            f.write('dds <- DESeq(dds)\n')
            ## 输出sizefactor
            tmp_sizeFactor = os.path.join(output, "{}_vs_{}.deseq2.sizeFactor.xls".format(ctrl, test))
            f.write('sizefactor <- data.frame(dds$sizeFactor)\n')
            f.write('colnames = c("sample", "sizefactor")\n')
            f.write('write.table(t(colnames), file="{}", sep="\\t", quote=F, row.names=F, col.names=F)\n'.format(
                tmp_sizeFactor))
            f.write(
                'write.table(sizefactor, file="{}", sep="\\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                    tmp_sizeFactor))
            ## 输出normalized结果文件
            tmp_normalize = os.path.join(output, "{}_vs_{}.deseq2.normalize.xls".format(ctrl, test))
            raw_tmp_normalize = os.path.join(output, "{}_vs_{}.deseq2.rawnormal.xls".format(ctrl, test))
            f.write('normalized_counts <- counts(dds, normalized=TRUE)\n')
            f.write('colnames = c("seq_id",colnames(normalized_counts))\n')
            f.write('write.table(t(colnames), file="{}", sep="\\t", quote=F, col.names=F, row.names=F)\n'.format(
                raw_tmp_normalize))
            f.write(
                'write.table(normalized_counts, file="{}", sep="\\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                    raw_tmp_normalize))
            if padjust_way is None:
                padjust_way = 'BH'
            f.write('res <- results(dds, contrast<-c("group", "{}", "{}"), '
                    'pAdjustMethod="{}")\n'.format(test, ctrl, padjust_way))
            tmp_out_file = os.path.join(output, "{}_vs_{}.deseq2.tmp".format(ctrl, test))
            f.write('write.table(res, file="{}", sep="\\t", quote=F, '
                    'row.names=T, col.names=NA)\n'.format(tmp_out_file))
            #add by fwy in 20200512 输出可以直接计算fc的标准化矩阵
            # f.write('res1 <- results(dds, contrast<-c("group", "{}", "{}"), '
            #         'pAdjustMethod="{}",addMLE = True)\n'.format(test, ctrl, padjust_way))
            f.write('dds1 <- DESeq(dds,betaPrior = TRUE)\n')
            f.write('res1 <- results(dds1, contrast<-c("group", "{}", "{}"), '
                    'pAdjustMethod="{}",addMLE = TRUE)\n'.format(test, ctrl, padjust_way))
            f.write('d = mcols(dds1)\n')
            f.write('d$group1= 2^d$MLE_Intercept\n')
            f.write('co  <- try((d$group2= 2^(d$MLE_Intercept+d$MLE_group_{}_vs_{})), silent=TRUE)\n'.format(test,ctrl))
            f.write('if (\'try-error\' %in% class(co))\n')
            f.write('d$group2= 2^(d$MLE_Intercept+d$MLE_group_{}_vs_{})\n'.format(ctrl,test))
            f.write('lista=rownames(res1)\n')
            f.write('d$seq_id=lista\n')
            f.write('write.table(d,'
                    ' file="{}", sep="\t", quote=F, row.names=T, col.names=NA)'.format(tmp_normalize))
            f.close()
        else:
            self.run(script_list)

        # make final report
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.deseq2.tmp')]
        all_stat_dicts = dict()
        all_stat_details = []
        final_details= []
        for each in cmp_result_dirs:
            final_detail_df = dask.delayed(self.detail_extract)(each,"DESeq2")
            final_details.append(final_detail_df)

        else:

            out_stat = os.path.join(output, 'diff_summary_deseq2.xls')
            final_details = dask.compute(*final_details)
            self.summary_stat(final_details, out_stat)



    def detail_extract(self,tmp_path,method):
        stat_table = pd.read_table(tmp_path, index_col=0)
        if method == "DESeq2":
            pvalues = stat_table['pvalue']
            if self.padjust_way != 3:
                pvalue_correct = self.multtest_correct
                df = pd.DataFrame(dict(pvalue=pvalues,
                                       padjust=pvalue_correct(pvalues, method=self.padjust_way),
                                       log2fc=stat_table['log2FoldChange']), index=pvalues.index)
            else:
                df = pd.DataFrame(dict(pvalue=pvalues,
                                       padjust=stat_table['padj'],
                                       log2fc=stat_table['log2FoldChange']), index=pvalues.index)
            ctrl, test = os.path.basename(tmp_path).split('.deseq2.tmp')[0].split('_vs_')
        elif method == "edgeR":
            pvalues = stat_table['PValue']
            df = pd.DataFrame(dict(pvalue=pvalues,
                                   padjust=self.multtest_correct(pvalues, method=self.padjust_way),
                                   log2fc=stat_table['logFC']), index=pvalues.index)
            ctrl, test = os.path.basename(tmp_path).split('.edger.tmp')[0].split('_vs_')
        elif method == "limma":
            try:
                pvalues = stat_table['P.Value']
            except:
                pvalues = stat_table['PValue']
            df = pd.DataFrame(dict(pvalue=pvalues,
                                   padjust=self.multtest_correct(pvalues, method=self.padjust_way),
                                   log2fc=stat_table['logFC']), index=pvalues.index)
            ctrl, test = os.path.basename(tmp_path).split('.limma.tmp')[0].split('_vs_')
        elif method == "svaseqlimma":
            pvalues = stat_table['P.Value']
            df = pd.DataFrame(dict(pvalue=pvalues,
                                   padjust=self.multtest_correct(pvalues, method=self.padjust_way),
                                   log2fc=stat_table['logFC']), index=pvalues.index)
            ctrl, test = os.path.basename(tmp_path).split('.svaseqlimma.tmp')[0].split('_vs_')

        detail_df, result_table = df, tmp_path.split('.tmp')[0] + '.xls'
        cmp = "{}_vs_{}".format(ctrl, test)
        log_fc_cutoff = math.log(self.fc_cutoff, 2)
        if ctrl in self.group_dict:
            ctrl_samples = self.group_dict[ctrl]
        else:
            ctrl_samples = [ctrl]
        if test in self.group_dict:
            test_samples = self.group_dict[test]
        else:
            test_samples = [test]
        detail_df['regulate'] = detail_df['log2fc'].fillna(0).apply(
            lambda x: 'up' if x > 0 else 'down' if x else 'no change')
        detail_df['significant'] = ['yes' if i and j else 'no' for i, j in
                                    zip((abs(detail_df['log2fc']) >= log_fc_cutoff),
                                        detail_df[self.sig_type].fillna(1) <= self.stat_cutoff)]
        detail_df["fc"] = detail_df["log2fc"].apply(lambda x: math.pow(2, x))
        select_columns = ["fc", "log2fc", "pvalue", "padjust", "significant", "regulate"]
        detail_df = detail_df[select_columns]
        cmp_samples = ctrl_samples + test_samples
        count_df = self.total_count_df[cmp_samples]
        count_df.columns = [i + "_count" for i in count_df.columns]
        exp_df = self.total_exp_df[cmp_samples]
        exp_df[ctrl] = exp_df.loc[:, ctrl_samples].mean(axis=1)
        exp_df[test] = exp_df.loc[:, test_samples].mean(axis=1)
        exp_df.columns = [i + "_" + self.exp_type for i in exp_df.columns]
        cmp_stat_df = pd.concat([count_df, exp_df, detail_df], axis=1)
        cmp_stat_df.index.name = "seq_id"
        cmp_stat_df = cmp_stat_df.fillna(
            {'fc': 1, 'log2fc': 0, "pvalue": 1, "padjust": 1, "significant": "no test", "regulate": "no test"})
        cmp_stat_df = cmp_stat_df.sort_values("pvalue")
        cmp_stat_df.to_csv(result_table, sep="\t")
        return cmp,cmp_stat_df

    def detail_extract_DEGseq(self,tmp_path,threshold_kind,cmpname):
        df = pd.read_table(tmp_path, index_col=0)
        if int(threshold_kind) == 5 or int(threshold_kind) == 4:
            padjust = 'q-value(Storey et al. 2003)'
        elif int(threshold_kind) == 1:
            padjust = 'p-value'
        elif int(threshold_kind) == 3:
            padjust = 'q-value(Benjamini et al. 1995)'
        elif int(threshold_kind) == 2:
            padjust = 'z-score'
        else:
            threshold_kind = 0
        pvalues = df['p-value']
        if threshold_kind:
            stat_df = pd.DataFrame(dict(pvalue=pvalues,
                                        padjust=df[padjust],
                                        log2fc=df['log2(Fold_change) normalized'], ),
                                   index=pvalues.index)
        else:
            correction = self.multtest_correct
            stat_df = pd.DataFrame(dict(pvalue=df['p-value'],
                                        padjust=correction(pvalues, method=self.padjust_way),
                                        log2fc=df['log2(Fold_change) normalized']),
                                   index=pvalues.index)
        ctrl, test = os.path.basename(cmpname).split('.degseq.tmp')[0].split('_vs_')

        detail_df, result_table = stat_df, tmp_path.split('.tmp')[0] + '.xls'
        cmp = "{}_vs_{}".format(ctrl, test)
        log_fc_cutoff = math.log(self.fc_cutoff, 2)
        if ctrl in self.group_dict:
            ctrl_samples = self.group_dict[ctrl]
        else:
            ctrl_samples = [ctrl]
        if test in self.group_dict:
            test_samples = self.group_dict[test]
        else:
            test_samples = [test]
        detail_df['regulate'] = detail_df['log2fc'].fillna(0).apply(
            lambda x: 'up' if x > 0 else 'down' if x else 'no change')
        detail_df['significant'] = ['yes' if i and j else 'no' for i, j in
                                    zip((abs(detail_df['log2fc']) >= log_fc_cutoff),
                                        detail_df[self.sig_type].fillna(1) <= self.stat_cutoff)]
        detail_df["fc"] = detail_df["log2fc"].apply(lambda x: math.pow(2, x))
        select_columns = ["fc", "log2fc", "pvalue", "padjust", "significant", "regulate"]
        detail_df = detail_df[select_columns]
        cmp_samples = ctrl_samples + test_samples
        count_df = self.total_count_df[cmp_samples]
        count_df.columns = [i + "_count" for i in count_df.columns]
        exp_df = self.total_exp_df[cmp_samples]
        exp_df[ctrl] = exp_df.loc[:, ctrl_samples].mean(axis=1)
        exp_df[test] = exp_df.loc[:, test_samples].mean(axis=1)
        exp_df.columns = [i + "_" + self.exp_type for i in exp_df.columns]
        cmp_stat_df = pd.concat([count_df, exp_df, detail_df], axis=1)
        cmp_stat_df.index.name = "seq_id"
        cmp_stat_df = cmp_stat_df.fillna(
            {'fc': 1, 'log2fc': 0, "pvalue": 1, "padjust": 1, "significant": "no test", "regulate": "no test"})
        cmp_stat_df = cmp_stat_df.sort_values("pvalue")
        cmp_stat_df.to_csv(result_table, sep="\t")
        return cmp,cmp_stat_df

    def detail_extract_NOIseq(self,tmp_path):
        ctrl, test = os.path.basename(tmp_path).split('.noiseq.tmp')[0].split('_vs_')
        stat_table = pd.read_table(tmp_path, index_col=0)
        if self.sample == 'True':
            theta = stat_table['theta']
            # df = pd.DataFrame(dict(theta=theta, mean_a=stat_table['{}_mean'.format(ctrl)], mean_b=stat_table['{}_mean'.format(ctrl)],
            #                        log2fc=stat_table['log2FC'], prob=stat_table['prob']), index=theta.index)
            df = pd.DataFrame({'theta': theta, '{}_mean'.format(ctrl): stat_table['{}_mean'.format(ctrl)],
                               '{}_mean'.format(test): stat_table['{}_mean'.format(test)],
                               'log2fc': stat_table['log2FC'], 'prob': stat_table['prob']}, index=theta.index)
        if self.sample == 'False':
            df = pd.DataFrame({'{}_mean'.format(ctrl): stat_table['{}_mean'.format(ctrl)],
                               '{}_mean'.format(test): stat_table['{}_mean'.format(test)],
                               'log2fc': stat_table['M'], 'prob': stat_table['prob']}, index=stat_table.index)
        detail_df, result_table = df, tmp_path.split('.tmp')[0] + '.xls'
        cmp = "{}_vs_{}".format(ctrl, test)
        log_fc_cutoff = math.log(self.fc_cutoff, 2)
        if ctrl in self.group_dict:
            ctrl_samples = self.group_dict[ctrl]
        else:
            ctrl_samples = [ctrl]
        if test in self.group_dict:
            test_samples = self.group_dict[test]
        else:
            test_samples = [test]
        detail_df['regulate'] = detail_df['log2fc'].fillna(0).apply(
            lambda x: 'up' if x > 0 else 'down' if x else 'no change')
        detail_df['significant'] = ['yes' if i and j else 'no' for i, j in
                                    zip((abs(detail_df['log2fc']) >= log_fc_cutoff),
                                        detail_df['prob'].fillna(0) >= self.prob)]
        detail_df["fc"] = detail_df["log2fc"].apply(lambda x: math.pow(2, x))
        detail_df["D"] = detail_df.apply(lambda x: abs(x['{}_mean'.format(ctrl)] - x['{}_mean'.format(test)]), axis=1)
        if self.sample == 'True':
            select_columns = ['{}_mean'.format(ctrl),'{}_mean'.format(test),"fc", "log2fc", "theta", "D", "prob", "significant", "regulate"]
        else:
            select_columns = ['{}_mean'.format(ctrl),'{}_mean'.format(test),"fc", "log2fc", "D", "prob", "significant", "regulate"]
        detail_df = detail_df[select_columns]
        cmp_samples = ctrl_samples + test_samples
        count_df = self.total_count_df[cmp_samples]
        count_df.columns = [i + "_count" for i in count_df.columns]
        exp_df = self.total_exp_df[cmp_samples]
        exp_df[ctrl] = exp_df.loc[:, ctrl_samples].mean(axis=1)
        exp_df[test] = exp_df.loc[:, test_samples].mean(axis=1)
        exp_df.columns = [i + "_" + self.exp_type for i in exp_df.columns]
        cmp_stat_df = pd.concat([count_df, exp_df, detail_df], axis=1)
        cmp_stat_df.index.name = "seq_id"
        cmp_stat_df = cmp_stat_df.fillna(
            {'fc': 1, 'log2fc': 0, "theta": 0, "prob": 0,"D": 0, "significant": "no test", "regulate": "no test"})
        cmp_stat_df = cmp_stat_df.sort_values("prob",ascending=False)
        cmp_stat_df.to_csv(result_table, sep="\t")
        return cmp, cmp_stat_df



    def limma(self, dispersion=0.1, output=None, sep='\t'):
        if output is None:
            output = os.getcwd()

        no_sig_final_use = defaultdict(list)

        script_list = list()
        for ctrl,test in self.cmp_list:
            used_samples = self.group_dict[ctrl] + self.group_dict[test]
            exp_total = pd.read_table(self.exp, index_col="seq_id")
            counts_cmp_use = self.get_filter_counts_df(exp_total, used_samples)
            counts_cmp_use.to_csv(os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test)), sep="\t")
            counts_table = os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test))
            no_sig_final_use[ctrl + test] = list(set(exp_total.index) - set(counts_cmp_use.index))
            with open(os.path.join(output,'{}_vs_{}.counts_filted'.format(ctrl, test)),"r") as f :
                header = f.readline().strip('\n').split(sep)
        #script_list = list()
        #for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)) for x in self.group_dict[ctrl]])
                ctrl_names = ["'{}'".format(ctrl) for x in self.group_dict[ctrl]]
                if self.batch is not None:
                    try:
                        ctrl_batch=["'{}'".format(self.batch_dict[x]) for x in self.group_dict[ctrl]]
                        batch_err = False
                    except:
                        batch_err = True

            else:
                ctrl_ind = str(header.index(test))
                ctrl_names = ["'{}'".format(ctrl)]
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)) for x in self.group_dict[test]])
                test_names = ["'{}'".format(test) for x in self.group_dict[test]]
                if self.batch is not None:
                    try:
                        test_batch = ["'{}'".format(self.batch_dict[x]) for x in self.group_dict[test]]
                        batch_err = False
                    except:
                        batch_err = True
            else:
                test_ind = str(header.index(test))
                test_names = ["'{}'".format(test)]
            if self.batch is not None:
                a = [x for x in ctrl_batch if x in test_batch]
                if a:
                    pass
                else:
                    batch_err = True

            all_sample_ind = []
            all_sample_names = []
            for sample in header[1:]:
                sample_batch = str(header.index(sample))
                all_sample_ind.append(sample_batch)
            all_sample_ind = ",".join(all_sample_ind)
            all_sample_names = ["'{}'".format(self.sample_dict[x]) for x in header[1:]]
            if self.batch is not None:
                all_batch = ["'{}'".format(self.batch_dict[x]) for x in header[1:]]
                if all_batch:
                    pass
                else:
                    batch_err = True

            script_name = os.path.join(output, 'limma.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(limma)\n')
            f.write('library(edgeR)\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                counts_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            # f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            # f.write('tmp_counts <- floor(tmp_counts+0.5)\n')
            # f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))
            f.write('all_counts <- counts[, c({})]\n'.format(all_sample_ind))
            f.write('all_counts <- floor(all_counts+0.5)\n')
            f.write('all_group <- c({})\n'.format(','.join(all_sample_names)))

            if self.batch is not None and batch_err == False:
                # f.write('tmp_batch <- c({})\n'.format(','.join(ctrl_batch + test_batch)))
                # f.write('tmp_batch <- as.factor(tmp_batch)\n')
                f.write('all_batch <- c({})\n'.format(','.join(all_batch)))
                f.write('all_batch <- as.factor(all_batch)\n')
            f.write('y <- DGEList(counts=all_counts, group=all_group)\n')
            f.write('y <- calcNormFactors(y)\n')
            # 输出normfactors
            tmp_normFactor = os.path.join(output, "{}_vs_{}.limma.normFactor.xls".format(ctrl, test))
            f.write('colnames = c("sample", "group","libsize","norm.factors")\n')
            f.write('write.table(t(colnames), file="{}", sep="\\t", quote=F, row.names=F, col.names=F)\n'.format(
                tmp_normFactor))
            f.write(
                'write.table(y$samples, file="{}", sep = "\t", quote = F, row.names = T, col.names = F, append = T)\n'.format(
                    tmp_normFactor))
            ## 输出normalized结果文件
            tmp_normalize = os.path.join(output, "{}_vs_{}.limma.normalize.xls".format(ctrl, test))
            f.write('s_normfacotrs=vector()\n')
            f.write('for(i in 1:nrow(y$samples))\n')
            f.write('s_normfacotrs[i]=y$samples[i,3]\n')
            f.write('raw_c=y$counts\n')
            f.write('for(i in 1:ncol(raw_c))\n')
            f.write('for(j in 1:nrow(raw_c))\n')
            f.write('raw_c[j,i]=raw_c[j,i]/s_normfacotrs[i]\n')
            f.write('colnames = c("gene_id",colnames(raw_c))\n')
            f.write('write.table(t(colnames), file="{}", sep="\t", quote=F, col.names=F,row.names=F)\n'.format(
                tmp_normalize))
            f.write('write.table(raw_c, file="{}", sep="\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                tmp_normalize))

            if (',' not in ctrl_ind) and (',' not in test_ind):
                # do sample vs sample if NO replicates
                f.write('result <- exactTest(y, all_group, dispersion={})\n'.format(dispersion))
                f.write('result = topTags(result, n=dim(result$counts)[1], '
                        'adjust.method="BH", sort.by="PValue")\n')
                # f.write('result <- topTable(result, n = nrow(tmp_counts))\n')
            else:
                if self.batch is not None and batch_err == False:
                    f.write('design <- model.matrix(~0+all_group+all_batch)\n')
                else:
                    f.write('design <- model.matrix(~0+all_group)\n')
                f.write('y_voom <- voom(y, design)\n')
                f.write('fit <- lmFit(y_voom, design)\n')
                f.write('con = makeContrasts(all_group{}-all_group{}, levels=design)\n'.format(
                    test, ctrl))
                f.write('fit2 <- contrasts.fit(fit, con)\n')
                f.write('fit3 <- eBayes(fit2)\n')
                f.write('result <- topTable(fit3, n = nrow(all_counts))\n')


            f.write('write.table(result, "{}/{}_vs_{}.limma.tmp", sep="\\t", row.names=T, '
                    'col.names=NA, quote=FALSE)\n'.format(output, ctrl, test))
            f.close()
        else:
            self.run(script_list)
        # make final report
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.limma.tmp')]
        all_stat_dicts = dict()
        all_stat_details=[]
        final_details=[]
        for each in cmp_result_dirs:
            final_detail_df = dask.delayed(self.detail_extract)(each, "limma")
            final_details.append(final_detail_df)

        else:
            out_stat = os.path.join(output, 'diff_summary_limma.xls')
            final_details = dask.compute(*final_details)
            self.summary_stat(final_details, out_stat)


    def svaseq_limma(self, dispersion=0.1, output=None, sep='\t'):
        if output is None:
            output = os.getcwd()

        no_sig_final_use = defaultdict(list)

        script_list = list()
        for ctrl,test in self.cmp_list:
            used_samples = self.group_dict[ctrl] + self.group_dict[test]
            exp_total = pd.read_table(self.exp, index_col="seq_id")
            counts_cmp_use = self.get_filter_counts_df(exp_total, used_samples)
            counts_cmp_use.to_csv(os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test)), sep="\t")
            counts_table = os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test))
            no_sig_final_use[ctrl + test] = list(set(exp_total.index) - set(counts_cmp_use.index))
            with open(os.path.join(output,'{}_vs_{}.counts_filted'.format(ctrl, test)),"r") as f :
                header = f.readline().strip('\n').split(sep)
        #script_list = list()
        #for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)) for x in self.group_dict[ctrl]])
                ctrl_names = ["'{}'".format(ctrl) for x in self.group_dict[ctrl]]
            else:
                ctrl_ind = str(header.index(test))
                ctrl_names = ["'{}'".format(ctrl)]
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)) for x in self.group_dict[test]])
                test_names = ["'{}'".format(test) for x in self.group_dict[test]]
            else:
                test_ind = str(header.index(test))
                test_names = ["'{}'".format(test)]
            script_name = os.path.join(output, 'svaseqlimma.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(limma)\n')
            f.write('library(edgeR)\n')
            f.write('library(sva)\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                counts_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            f.write('tmp_counts <- floor(tmp_counts+0.5)\n')
            f.write('tmp_counts_filter <- tmp_counts[which(rowSums(tmp_counts)>0),]\n')
            f.write('tmp_counts_filter <- as.matrix(tmp_counts_filter)\n')
            f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))
            f.write('mod <- model.matrix(~0+tmp_group)\n')
            f.write('mod0 <- cbind(rep(1,length(colnames(tmp_counts))))\n')
            f.write('svobj <- svaseq(tmp_counts_filter, mod, mod0)\n')
            f.write('sv <- svobj$sv\n')
            f.write('sv <- as.matrix(sv)\n')
            f.write('colnames(sv) <- colnames(sv,do.NULL=FALSE,prefix = "sv_")\n')
            f.write('modSv <- cbind(mod,sv)\n')
            f.write('y <- DGEList(counts=tmp_counts, group=tmp_group)\n')
            f.write('y <- calcNormFactors(y)\n')
            # 输出normfactors
            tmp_normFactor = os.path.join(output, "{}_vs_{}.svaseqlimma.normFactor.xls".format(ctrl, test))
            f.write('colnames = c("sample", "group","libsize","norm.factors")\n')
            f.write('write.table(t(colnames), file="{}", sep="\\t", quote=F, row.names=F, col.names=F)\n'.format(
                tmp_normFactor))
            f.write(
                'write.table(y$samples, file="{}", sep = "\t", quote = F, row.names = T, col.names = F, append = T)\n'.format(
                    tmp_normFactor))
            ## 输出normalized结果文件
            tmp_normalize = os.path.join(output, "{}_vs_{}.svaseqlimma.normalize.xls".format(ctrl, test))
            f.write('s_normfacotrs=vector()\n')
            f.write('for(i in 1:nrow(y$samples))\n')
            f.write('s_normfacotrs[i]=y$samples[i,3]\n')
            f.write('raw_c=y$counts\n')
            f.write('for(i in 1:ncol(raw_c))\n')
            f.write('for(j in 1:nrow(raw_c))\n')
            f.write('raw_c[j,i]=raw_c[j,i]/s_normfacotrs[i]\n')
            f.write('colnames = c("gene_id",colnames(raw_c))\n')
            f.write('write.table(t(colnames), file="{}", sep="\t", quote=F, col.names=F,row.names=F)\n'.format(
                tmp_normalize))
            f.write('write.table(raw_c, file="{}", sep="\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                tmp_normalize))

            if (',' not in ctrl_ind) and (',' not in test_ind):
                # do sample vs sample if NO replicates
                f.write('result <- exactTest(y, tmp_group, dispersion={})\n'.format(dispersion))
            else:

                f.write('y_voom <- voom(y, modSv)\n')
                f.write('fit <- lmFit(y_voom, modSv)\n')
                f.write('con = makeContrasts(tmp_group{}-tmp_group{}, levels=modSv)\n'.format(
                    test, ctrl))
                f.write('fit2 <- contrasts.fit(fit, con)\n')
                f.write('fit3 <- eBayes(fit2)\n')
            f.write('result <- topTable(fit3, n = nrow(tmp_counts))\n')


            f.write('write.table(result, "{}/{}_vs_{}.svaseqlimma.tmp", sep="\\t", row.names=T, '
                    'col.names=NA, quote=FALSE)\n'.format(output, ctrl, test))
            f.close()
        else:
            self.run(script_list)
        # make final report
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.svaseqlimma.tmp')]
        all_stat_dicts = dict()
        all_stat_details = []
        final_details = []
        for each in cmp_result_dirs:
            final_detail_df = dask.delayed(self.detail_extract)(each, "svaseqlimma")
            final_details.append(final_detail_df)
        else:

            out_stat = os.path.join(output, 'diff_summary_svaseqlimma.xls')
            final_details = dask.compute(*final_details)
            self.summary_stat(final_details, out_stat)



    def NOIseq(self, dispersion=0.1, output=None, sep='\t'):
        if output is None:
            output = os.getcwd()

        no_sig_final_use = defaultdict(list)

        script_list = list()
        for ctrl,test in self.cmp_list:
            used_samples = self.group_dict[ctrl] + self.group_dict[test]
            exp_total = pd.read_table(self.exp, index_col="seq_id")
            counts_cmp_use = self.get_filter_counts_df(exp_total, used_samples)
            counts_cmp_use.to_csv(os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test)), sep="\t")
            counts_table = os.path.join(output, '{}_vs_{}.counts_filted'.format(ctrl, test))
            no_sig_final_use[ctrl + test] = list(set(exp_total.index) - set(counts_cmp_use.index))
            with open(os.path.join(output,'{}_vs_{}.counts_filted'.format(ctrl, test)),"r") as f :
                header = f.readline().strip('\n').split(sep)
        #script_list = list()
        #for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)) for x in self.group_dict[ctrl]])
                ctrl_names = ["'{}'".format(ctrl) for x in self.group_dict[ctrl]]
                if self.batch is not None:
                    ctrl_batch=[self.batch_dict[x] for x in self.group_dict[ctrl]]
            else:
                ctrl_ind = str(header.index(test))
                ctrl_names = ["'{}'".format(ctrl)]
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)) for x in self.group_dict[test]])
                test_names = ["'{}'".format(test) for x in self.group_dict[test]]
                if self.batch is not None:
                    test_batch = [self.batch_dict[x] for x in self.group_dict[test]]
            else:
                test_ind = str(header.index(test))
                test_names = ["'{}'".format(test)]

            all_sample_ind = []
            all_sample_names = []
            for sample in header[1:]:
                sample_batch = str(header.index(sample))
                all_sample_ind.append(sample_batch)
            all_sample_ind = ",".join(all_sample_ind)
            all_sample_names = ["'{}'".format(self.sample_dict[x]) for x in header[1:]]
            if self.batch is not None:
                all_batch = ["'{}'".format(self.batch_dict[x]) for x in header[1:]]
                if all_batch:
                    pass
                else:
                    batch_err = True

            script_name = os.path.join(output, 'NOIseq.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(NOISeq)\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                counts_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            # f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            # # f.write('tmp_counts <- floor(tmp_counts+0.5)\n')
            # f.write('tmp_col <- colnames(tmp_counts)\n')
            # f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))

            f.write('all_counts <- counts[, c({})]\n'.format(all_sample_ind))
            # f.write('tmp_counts <- floor(tmp_counts+0.5)\n')
            f.write('all_col <- colnames(all_counts)\n')
            f.write('all_group <- c({})\n'.format(','.join(all_sample_names)))
            f.write('myfactors <- data.frame(tmp_group, row.names=all_col)\n')
            f.write('myfactors\n')

            f.write("## Normalized \n")
            if self.normalized == 'rpkm':
                f.write('normalized_count <- rpkm(all_counts, long = mylength, lc = 1, k = 0.5)\n')
            if self.normalized == 'tmm':
                f.write('normalized_count <- tmm(all_counts, long = 1000, lc = 0, k = 0.5)\n')
            if self.normalized == 'uqua':
                f.write('normalized_count <- uqua(all_counts, long = 1000, lc = 0, k = 0.5)\n')
            f.write('colnames <- c("seq_id",colnames(normalized_count))\n')
            tmp_normalize = os.path.join(output, "{}_vs_{}.NOIseq.normalize.xls".format(ctrl, test))
            f.write('write.table(t(colnames), file="{}", sep="\t", quote=F, col.names=F,row.names=F)\n'.format(
                tmp_normalize))
            f.write('write.table(normalized_count, file="{}", sep="\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                tmp_normalize))
            f.write("## Normalized \n")
            f.write('mydata <- readData(data = normalized_count, factors = myfactors)\n')
            if self.sample == 'False':
                f.write('mynoiseq <- noiseq(mydata, norm = "n", replicates = "no", factor = "{}", conditions = c("{}", "{}"))\n'.format('all_group',
                                                                                                               ctrl, test))

            if self.sample == 'True':
                f.write('mynoiseq <- noiseqbio(mydata, norm = "n", factor = "{}", conditions = c("{}", "{}"))\n'.format('all_group'
                                                                                                                   , ctrl, test))
            f.write('results <- mynoiseq@results[[1]]\n')
            f.write('colnames <- c("seq_id",colnames(results))\n')
            result = os.path.join(output, "{}_vs_{}.noiseq.tmp".format(ctrl, test))
            f.write('write.table(t(colnames), file="{}", sep="\t", quote=F, col.names=F,row.names=F)\n'.format(
                result))
            f.write('write.table(results, file="{}", sep="\t", quote=F, row.names=T, col.names=F, append=T)\n'.format(
                result))
            f.close()
        else:
            self.run(script_list)
        # make final report
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.noiseq.tmp')]
        all_stat_dicts = dict()
        all_stat_details=[]
        final_details=[]
        for each in cmp_result_dirs:
            prepare_list = locals()
            final_detail_df = dask.delayed(self.detail_extract_NOIseq)(each)
            final_details.append(final_detail_df)

        else:
            out_stat = os.path.join(output, 'diff_summary_noiseq.xls')
            final_details = dask.compute(*final_details)
            self.summary_stat(final_details, out_stat)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-count', type=str, required=True,
                        help='path of read count table. Sample name is used as column name. '
                             'Second column must be gene length if "-exp" is None.')
    parser.add_argument('-exp', type=str, default=None,
                        help="path of expression value table, tab as separator."
                             " If None, the second column of count_matrix must be gene length which"
                             " will be used to calculate fpkm or tpm. NOTE: Expression table "
                             "has nothing to do with differential analysis; Only used in report.")
    parser.add_argument('--exp_type', type=str, default="fpkm", help='fpkm or tpm. Default: fpkm')
    parser.add_argument('-group', type=str, required=True,
                        help="path of group info file with at least two columns. First column must"
                             " consist of sample names. Other columns consist of group names."
                             "if no replicate exist, just use sample name as group name. "
                             "Header line starts with '#'")
    parser.add_argument('-cmp', type=str, required=True,
                        help="path of comparison info file with only two columns(ctrl vs test)."
                             " Header line starts with '#'")
    parser.add_argument('-method', type=str, default="edgeR", help='DEGseq or edgeR or DESeq2')
    parser.add_argument('--no_filter', default=False, action='store_true',
                        help='Do no filtering. This option will be ignored by default.')
    parser.add_argument('-output', type=str, default=None, help='output directory.')
    parser.add_argument('-pool', type=int, default=5, help='process number for computing')
    parser.add_argument('--plot', default=False, action='store_true', help="do plotting")
    parser.add_argument('-pvalue', type=float, default=0.05, help='p(q)value cutoff. Default: 0.05')
    parser.add_argument('-fc', type=float, default=2.0, help='fold change cutoff. Default: 2.0')
    parser.add_argument('--count_cutoff', type=float, default=4.0,
                        help='count number cutoff for filtering before diff analysis. Default: 4.0')
    parser.add_argument('--filter_method', type=str, default=None,
                        help='filter_method for tpm filtering before diff analysis,string in ["num","average","min","max"]'
                             'Default: None')
    parser.add_argument('--tpm_filter_threshold', type=float, default=0,
                        help='tpm_filter_threshold cutoff for filtering before diff analysis,a number >0'
                             'Default: 0.0')
    parser.add_argument('--passed_number_cutoff', type=int, default=None,
                        help='sample( count > count_cutoff ) number cutoff for filtering before '
                             'diff analysis. Let M=passed_number_cutoff, N=total_sample_number, '
                             'the following event must happen for a gene to be tested: '
                             'Each gene_count in M samples out of N must >= "count_cutoff". '
                             'Default: self-determined')
    parser.add_argument('--degseq_method', type=str, default='MARS',
                        help='method of degseq. Default: MARS')
    parser.add_argument('--degseq_padjust_way', type=int, default=0,
                        help="param of degseq. Integer in [0,1,2,3,4,5]. "
                             " Default: 5 for using qValue (Storey et al. 2003). "
                             "If 0, this option will be ignored; If 2, z-score will be used.")
    parser.add_argument('-sig_type', default="padjust", type=str,
                        help="pvalue or padjust, for diff significance judgement. Default: padjust")
    parser.add_argument('--dispersion', type=float, default=0.1,
                        help='Only used for single sample vs single sample with edgeR.Default: 0.1')
    parser.add_argument('--deseq2_padjust_way', type=str, default=None,
                        help='One of ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",)'
                             'Default: None')
    parser.add_argument('-padjust_way', type=int, default=3,
                        help='Integer in [1,2,3,4]. '
                             'This option is not used if other *_padjust_way specified.'
                             ' 1. Bonferroni. ---> bonferroni;'
                             ' 2. Bonferroni Step-down(Holm) ' '---> Holm;'
                             ' 3. Benjamini and Hochberg False Discovery Rate ---> BH;'
                             ' 4. FDR Benjamini-Yekutieli --->BY'
                             ' Default: 3')
    parser.add_argument('--batch', type=str, default=None)
    parser.add_argument('--normailized', type=str, default=None)
    parser.add_argument('--is_duplicate', default=True)
    parser.add_argument('--prob', type=float, default=0.9)

    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    toolbox = DiffExpToolbox(args.count, args.group, args.cmp,
                             exp_matrix=args.exp,
                             exp_type=args.exp_type,
                             sig_type=args.sig_type,
                             fc_cutoff=args.fc,
                             stat_cutoff=args.pvalue,
                             padjust_way=args.padjust_way,
                             pool_size=args.pool,
                             tpm_filter_threshold=args.tpm_filter_threshold,
                             batch_info=args.batch,
                             prob=args.prob,
                             is_duplicate=args.is_duplicate,
                             filter_method=args.filter_method
                             # normalized=args.normalized,
                             # sample=args.sample

                             )
#    if not args.no_filter:
#        toolbox.filter()

    if args.method == 'DEGseq':
        toolbox.DEGseq(method=args.degseq_method, threshold_kind=args.degseq_padjust_way,
                       output=args.output, )
    elif args.method == 'edgeR':
        toolbox.edgeR(dispersion=args.dispersion, output=args.output, )
    elif args.method == 'DESeq2':
        toolbox.DESeq2(output=args.output, padjust_way=args.deseq2_padjust_way, )
    elif args.method.lower() == 'limma':
        toolbox.limma(dispersion=args.dispersion, output=args.output,)
    elif args.method == 'svaseqlimma':
        toolbox.svaseq_limma(dispersion=args.dispersion, output=args.output)
    elif args.method == 'NOIseq':
        toolbox.NOIseq(dispersion=args.dispersion, output=args.output)
    else:
        raise Exception('Method {} is not supported'.format(args.method))

    def diff_plot(table):
        df = pd.read_table(table, index_col=0, header=0)
        colors = pd.DataFrame(['gray']*df.shape[0], index=df.index, columns=['color'])
        colors[df['significant'] == 'yes'] = 'red'
        colors[df['regulate'] == "down"] = 'green'
        colors[df['significant'] == 'no'] = 'gray'
        scatter_df = pd.concat([np.log(df.iloc[:, -7]+1), np.log(df.iloc[:, -6]+1), colors], axis=1)
        scatter_df.sort_values(by='color', axis=0, ascending=True, inplace=True)
        x_label, y_label, _ = scatter_df.columns
        scatter_df.plot.scatter(x=x_label, y=y_label, c=scatter_df['color'])
        plt.savefig(os.path.join(args.output, x_label + '_vs_' + y_label + '.scatter.png'), dpi=300)

        volcano_df = pd.concat([df['log2fc'], -np.log10(df[args.sig_type]), colors], axis=1)
        volcano_df.columns = ['log2fc', '-log10(' + args.sig_type + ')', 'color']
        volcano_df.sort_values(by='color', axis=0, ascending=True, inplace=True)
        x_label2, y_label2, _ = volcano_df.columns
        volcano_df.plot.scatter(x=x_label2, y=y_label2, c=volcano_df['color'])
        plt.savefig(os.path.join(args.output, x_label + '_vs_' + y_label + '.volcano.png'), dpi=300)

    def density_plot():
        exp = pd.read_table(toolbox.exp, index_col=0, header=0)
        exp = exp[exp.mean(axis=1) >= 0.05]
        exp_df = np.log(exp+1).dropna()
        exp_df.plot(kind="density",)
        plt.xlabel('log2(exp)')
        plt.savefig(os.path.join(args.output, 'exp_based.density.png'), dpi=300)

        exp = pd.read_table(toolbox.count, index_col=0, header=0)
        exp = exp[exp.mean(axis=1) >= 0.8]
        exp_df = np.log(exp+1).dropna()
        exp_df.plot(kind="density",)
        plt.xlabel('log2(exp)')
        plt.savefig(os.path.join(args.output, 'count_based.density.png'), dpi=300)

    # plotting
    if args.plot:
        if args.output is None:
            args.output = os.getcwd()

        results = glob.glob(args.output+'/*_vs_*.{}.xls'.format(args.method.lower()))
        # print(results)
        density_plot()
        pool = Pool(args.pool)
        pool.map(diff_plot, results)
        pool.close()
        pool.join()
        # done


