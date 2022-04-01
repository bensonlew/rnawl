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
from collections import defaultdict
from matplotlib import pyplot as plt
import logging

def run_script(codes):
    cmd = "Rscript {}".format(codes)
    retcode = subprocess.call(cmd, shell=True)
    if retcode:
        msg = 'fail to excecute {}'.format(cmd)
        raise Exception(msg)
    else:
        logging.info('succeed in excecuting command ({})'.format(cmd))
    # code = os.system(cmd)
    # if code:
    #     raise Exception("命令{}执行失败！".format(cmd))
    # else:
    #     logging.info("命令{}执行成功！".format(cmd))
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

class Diffgeneset4limma(PvalueCorrect):
    """
    The LIMMA package was used for GSVA different analysis
    """
    def __init__(self, count_matrix, group_info, cmp_info, sig_type='pvalue', stat_cutoff=0.05, fc_cutoff=2,
                 padjust_way=3, pool_size=5,tpm_filter_threshold=0,):
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
        if sig_type not in ['pvalue', 'padjust']:
            raise Exception('sig_type is not pvalue/padjust')
        self.sig_type = sig_type
        self.fc_cutoff = fc_cutoff
        self.stat_cutoff = stat_cutoff
        self.padjust_way = padjust_way
        self.tpm_filter_threshold = tpm_filter_threshold
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
        cmp_list = sorted(list(set(cmp_list)))
        self.cmp_list = cmp_list
        # print sample info
        pprint('group_dict is: ')
        pprint(self.group_dict)
        pprint('comparison list is (ctrl, test): ')
        pprint(self.cmp_list)
        if error_names:
            raise Exception('Each group name of {} is not in {}!'.format(error_names, group_info))
        # check the consistency between group info and count table
        with open(self.count) as f:
            count_samples = f.readline().strip('\n').split('\t')[1:]
        diff = set(self.samples).difference(set(count_samples))
        if diff:
            raise Exception('samples: {} are not contained in count table file'.format(diff))
        df = pd.read_table(self.count, index_col=0)
        self.count_dicts = df.to_dict('index')
    def run(self, script_list):
        pool = Pool(self.pool_size)
        pool.map(run_script, script_list)
        pool.close()
        pool.join()

    def __make_result(self, ctrl, test, target_seqs, stat_dict, out_diff_table, out_deg_list):
        if ctrl in self.group_dict:
            ctrl_samples = self.group_dict[ctrl]
        else:
            ctrl_samples = [ctrl]
        if test in self.group_dict:
            test_samples = self.group_dict[test]
        else:
            test_samples = [test]
        # with open(out_diff_table, 'w') as f, open(out_deg_list, 'w') as f2:
        with open(out_diff_table, 'w') as f:
            count_header = '_count\t'.join(ctrl_samples+test_samples) + '_count'
            # tmp_sep = '_' + self.exp_type + '\t'
            # exp_header_list = ctrl_samples+test_samples
            # if len(ctrl_samples) >= 1:
            #     exp_header_list.append(ctrl)
            # if len(test_samples) >= 1:
            #     exp_header_list.append(test)
            # exp_header = tmp_sep.join(exp_header_list) + '_' + self.exp_type

            f.write('geneset\t{}\tfc\tlog2fc\tpvalue\tpadjust\n'.format(
                count_header))
            cmp_samples = ctrl_samples + test_samples
            for geneset in target_seqs:
                line_list = [geneset]
                tmp_count_dict = self.count_dicts[geneset]
                line_list += [tmp_count_dict[x] for x in cmp_samples]
                # tmp_exp_dict = self.exp_dicts[geneset]
                # tmp_exp_list = [tmp_exp_dict[x] for x in cmp_samples]
                # if len(ctrl_samples) >= 1:
                #     tmp_exp_list.append(sum([tmp_exp_dict[x] for x in ctrl_samples])/len(ctrl_samples))
                # if len(test_samples) >= 1:
                #     tmp_exp_list.append(sum([tmp_exp_dict[x] for x in test_samples])/len(test_samples))
                # line_list += tmp_exp_list
                tmp_stat_dict = stat_dict.get(geneset)
                if tmp_stat_dict:
                    # get fold change
                    tmp_lfc = tmp_stat_dict['log2fc']
                    tmp_lfc = 0 if tmp_lfc != tmp_lfc else tmp_lfc
                    tmp_fc= 2 ** tmp_lfc
                    line_list.append(tmp_fc)
                    # get fold change
                    tmp_lfc = tmp_stat_dict['log2fc']
                    tmp_lfc = 0 if tmp_lfc != tmp_lfc else tmp_lfc
                    line_list.append(tmp_lfc)
                    # get pvalue
                    tmp_pvalue = tmp_stat_dict['pvalue']
                    tmp_pvalue = 1 if tmp_pvalue != tmp_pvalue else tmp_pvalue
                    line_list.append(tmp_pvalue)
                    # get adjusted pvalue
                    tmp_padjust = tmp_stat_dict['padjust']
                    tmp_padjust = 1 if tmp_padjust != tmp_padjust else tmp_padjust
                    line_list.append(tmp_padjust)
                    # judge significant
                    # tmp_stat = tmp_stat_dict[self.sig_type]
                    # if abs(tmp_lfc) >= 0.2 and tmp_stat <= self.stat_cutoff:
                    #     line_list.append('yes')
                    # else:
                    #     line_list.append('no')
                    # # judge regulate
                    # if tmp_lfc == 0:
                    #     reg = 'no change'
                    # elif tmp_lfc > 0:
                    #     reg = 'up'
                    # else:
                    #     reg = 'down'
                    # line_list.append(reg)
                    # save DEG list
                    # if line_list[-2] == 'yes':
                    #     f2.write(geneset + '\t' + reg + '\n')
                # else:
                #     line_list += [1,0, 1, 1, 'no test', 'no test']
                # save
                f.write('\t'.join([str(x) for x in line_list])+'\n')

    # def __diff_stat(self, all_stat_dicts, out_stat):
    #     significant_info = dict()
    #     total_deg = set()
    #     cmp_deg_sum = list()
    #     for each_cmp in all_stat_dicts.keys():
    #         tmp_table = all_stat_dicts[each_cmp]
    #         tmp_table = tmp_table[abs(tmp_table['log2fc']) >= 0.2]
    #         deg_list = tmp_table[tmp_table[self.sig_type] <= self.stat_cutoff].index.values
    #         significant_info[each_cmp] = deg_list
    #         total_deg.update(deg_list)
    #         cmp_deg_sum.append(len(deg_list))
    #     with open(out_stat, 'w') as f:
    #         f.write('geneset\t{}\tsum\n'.format('\t'.join(all_stat_dicts.keys())))
    #         sum_info = '\t'.join([str(x) for x in cmp_deg_sum])
    #         f.write('{}\t{}\t{}\n'.format(len(self.count_dicts.keys()), sum_info, len(total_deg)))
    #         all_seqs = list(total_deg)+list(set(self.count_dicts.keys())-total_deg)
    #         for seq_id in all_seqs:
    #             yes_no_lst = ['yes' if seq_id in significant_info[tmp_cmp] else 'no'
    #                           for tmp_cmp in all_stat_dicts.keys()]
    #             f.write('{}\t{}\t{}\n'.format(seq_id, '\t'.join(yes_no_lst), yes_no_lst.count('yes')))

    def limma(self, dispersion=0.1, output=None, sep='\t'):
        if output is None:
            output = os.getcwd()

        script_list = list()
        for ctrl,test in self.cmp_list:
            used_samples=self.group_dict[ctrl]+self.group_dict[test]
            counts_total=pd.read_table(self.count,index_col="geneset")
            counts_cmp_total=pd.DataFrame(counts_total,columns=used_samples)
            counts_cmp_total.to_csv(os.path.join(output,'{}_vs_{}.counts'.format(ctrl, test)),sep="\t")
            counts_table=os.path.join(output,'{}_vs_{}.counts'.format(ctrl, test))
            with open(os.path.join(output,'{}_vs_{}.counts'.format(ctrl, test)),"r") as f :
                header = f.readline().strip('\n').split(sep)
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
            script_name = os.path.join(output, 'limma.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(limma)\n')
            f.write('library(edgeR)\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                counts_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))
            if (',' not in ctrl_ind) and (',' not in test_ind):
                # do sample vs sample if NO replicates
                f.write('result <- exactTest(tmp_counts, tmp_group, dispersion={})\n'.format(dispersion))
                f.write('result = topTags(result, n=dim(result$counts)[1], '
                        'adjust.method="BH", sort.by="PValue")\n')
                # f.write('result <- topTable(result, n = nrow(tmp_counts))\n')
            else:
                f.write('design <- model.matrix(~0+tmp_group)\n')
                f.write('fit <- lmFit(tmp_counts, design)\n')
                f.write('con = makeContrasts(tmp_group{}-tmp_group{}, levels=design)\n'.format(
                    test, ctrl))
                f.write('fit2 <- contrasts.fit(fit, con)\n')
                f.write('fit3 <- eBayes(fit2)\n')
                f.write('result <- topTable(fit3, n = nrow(tmp_counts))\n')
            f.write('write.table(result, "{}/{}_vs_{}.limma.tmp", sep="\\t", row.names=T, '
                    'col.names=NA, quote=FALSE)\n'.format(output, ctrl, test))
            f.close()
        else:
            self.run(script_list)
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.limma.tmp')]
        all_stat_dicts = dict()
        for each in cmp_result_dirs:
            stat_table = pd.read_table(each, index_col=0)
            try:
                pvalues = stat_table['P.Value']
            except:
                pvalues = stat_table['PValue']
            df = pd.DataFrame(dict(pvalue=pvalues,
                                   padjust=self.multtest_correct(pvalues, method=self.padjust_way),
                                   log2fc=stat_table['logFC']), index=pvalues.index)
            stat_dict = df.to_dict('index')
            ctrl, test = os.path.basename(each).split('.limma.tmp')[0].split('_vs_')
            target_seqs = list(df[self.sig_type].sort_values().index)
            result_table = each.split('.tmp')[0] + '.xls'
            result_delist = each.split('.tmp')[0] + '.DE.list'
            all_stat_dicts['{}_vs_{}'.format(ctrl, test)] = df
            self.__make_result(ctrl, test, target_seqs, stat_dict, result_table, result_delist)
        # else:
        #     out_stat = os.path.join(output, 'diff_summary_limma.xls')
        #     self.__diff_stat(all_stat_dicts, out_stat)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-count', type=str, required=True,
                        help='path of read count table. Sample name is used as column name. '
                             'Second column must be gene length if "-exp" is None.')
    parser.add_argument('-group', type=str, required=True,
                        help="path of group info file with at least two columns. First column must"
                             " consist of sample names. Other columns consist of group names."
                             "if no replicate exist, just use sample name as group name. "
                             "Header line starts with '#'")
    parser.add_argument('-cmp', type=str, required=True,
                        help="path of comparison info file with only two columns(ctrl vs test)."
                             " Header line starts with '#'")
    parser.add_argument('-output', type=str, default=None, help='output directory.')
    parser.add_argument('-pool', type=int, default=5, help='process number for computing')
    parser.add_argument('-pvalue', type=float, default=0.05, help='p(q)value cutoff. Default: 0.05')
    parser.add_argument('-fc', type=float, default=2.0, help='fold change cutoff. Default: 2.0')
    parser.add_argument('--count_cutoff', type=float, default=4.0,
                        help='count number cutoff for filtering before diff analysis. Default: 4.0')
    parser.add_argument('--tpm_filter_threshold', type=float, default=0,
                        help='tpm_filter_threshold cutoff for filtering before diff analysis,a number >0'
                             'Default: 0.0')
    parser.add_argument('--passed_number_cutoff', type=int, default=None,
                        help='sample( count > count_cutoff ) number cutoff for filtering before '
                             'diff analysis. Let M=passed_number_cutoff, N=total_sample_number, '
                             'the following event must happen for a gene to be tested: '
                             'Each gene_count in M samples out of N must >= "count_cutoff". '
                             'Default: self-determined')
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
    args = parser.parse_args()
    diff = Diffgeneset4limma(args.count, args.group, args.cmp,
                             sig_type=args.sig_type,
                             fc_cutoff=args.fc,
                             stat_cutoff=args.pvalue,
                             # padjust_way=args.padjust_way,
                             pool_size=args.pool,
                             tpm_filter_threshold=args.tpm_filter_threshold,
                             )
    diff.limma(dispersion=args.dispersion, output=args.output, )