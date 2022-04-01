
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# go的pvalue矫正改为跟kegg类似的函数--by 封一统 20180713

from biocluster.workflow import Workflow
import glob
import os
import pandas as pd
import numpy as np
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json


class PrKeggEnrichWorkflow(Workflow):
    """
    蛋白集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(PrKeggEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_list", "type": "string"},
            {"name": "protein_list", "type": "string"},
            {"name": "protein_kegg", "type": "string"},
            {"name": "gene_kegg", "type": "string"},
            {"name": "protein_info", "type": "string"},
            {"name": "gene_info", "type": "string"},
            {"name": "rna_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "method", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.protein_tools = list()
        self.rna_tools = list()
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/7_relate/03_relate_anno/05_relate_keggenrich')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(PrKeggEnrichWorkflow, self).send_log(data)

    def run_tools(self):
        for protein_list in self.option('protein_list').split(','):
            tool = self.add_tool("itraq_and_tmt.proteinset.kegg_rich")
            options = {
                "diff_list": protein_list,
                "add_info": self.option("protein_info"),
                "correct": self.option("method"),
                "kegg_table": self.option("protein_kegg"),
            }
            tool.set_options(options)
            self.protein_tools.append(tool)

        for gene_list in self.option('gene_list').split(','):
            tool = self.add_tool("itraq_and_tmt.proteinset.kegg_rich")
            options = {
                "diff_list": gene_list,
                "add_info": self.option("gene_info"),
                "correct": self.option("method"),
                "kegg_table": self.option("gene_kegg"),
            }
            tool.set_options(options)
            self.rna_tools.append(tool)

    def run(self):
        self.run_tools()
        self.on_rely(self.protein_tools + self.rna_tools, self.set_db)
        for tool in self.protein_tools + self.rna_tools:
            tool.run()
        super(PrKeggEnrichWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        def change_qvalue(file):
            pf = pd.read_table(file, header=0, sep="\t")
            pf['p_corrected'] = self.multtest_correct(pf['p_uncorrected'].tolist(), self.option("method").lower())
            pf.drop([9])
            os.remove(file)
            pf.to_csv(file, sep='\t', index=False)
        api_goenrich = self.api.api('protein_transcript.relaset')
        protein_file = glob.glob("{}/*.xls".format(self.protein_tools[0].output_dir))[0]
        # change_qvalue(protein_file)
        # target = os.path.join(self.output_dir, os.path.basename(protein_file).split('.xls')[0] + '_protein.xls')
        target = os.path.join(self.output_dir, os.path.basename(protein_file))
        if os.path.exists(target):
            os.remove(target)
        os.link(protein_file, target)
        rna_file = glob.glob("{}/*.xls".format(self.rna_tools[0].output_dir))[0]
        # change_qvalue(rna_file)
        # target = os.path.join(self.output_dir, os.path.basename(rna_file).split('.xls')[0] + '_rna.xls')
        target = os.path.join(self.output_dir, os.path.basename(rna_file))
        if os.path.exists(target):
            os.remove(target)
        os.link(rna_file, target)
        api_goenrich.add_kegg_enrich_pr(self.option("main_table_id"), protein_file, rna_file)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        pt_keggenrich = self.output_dir
        chart.chart_pt_keggenrich(pt_keggenrich)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = os.path.join(self.work_dir, "pt_keggenrich.shadow_bar_groups.pdf")
        if os.path.exists(pdf_file):
            os.link(pdf_file, os.path.join(self.output_dir, "bar.pdf"))
        pdf_file = os.path.join(self.work_dir, "pt_keggenrich_bubble.bubble_groups.pdf")
        if os.path.exists(pdf_file):
            os.link(pdf_file, os.path.join(self.output_dir, "bubble.pdf"))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析",0],
            ["7_relate/03_relate_anno", "", "功能注释信息", 0],
            ["7_relate/03_relate_anno/05_relate_keggenrich", "", "KEGG富集", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白集与基因集GO富集分析结果文件"],
            ["bar.pdf", "", "直方图"],
            ["bubble.pdf", "", "气泡图"],
        ])
        super(PrKeggEnrichWorkflow, self).end()


    def multtest_correct(self, p_values, methods=3): #复制于package下的kegg_rich.py
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
                return np.arange(1, nobs + 1) / float(nobs)

            pvals = np.asarray(pvals)
            if not is_sorted:
                pvals_sortind = np.argsort(pvals)
                pvals_sorted = np.take(pvals, pvals_sortind)
            else:
                pvals_sorted = pvals  # alias

            if method in ['i', 'indep', 'p', 'poscorr']:
                ecdffactor = _ecdf(pvals_sorted)
            elif method in ['n', 'negcorr']:
                cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))  # corrected this
                ecdffactor = _ecdf(pvals_sorted) / cm
            else:
                raise ValueError('only indep and negcorr implemented')
            reject = pvals_sorted <= ecdffactor * alpha
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

        pvalue_list = list(p_values)
        n = len(pvalue_list)
        fdr = list()
        if methods == 'bonferroni':
            fdr = [eachP * n for eachP in pvalue_list]
        elif methods == 'holm':
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP * (n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
        elif methods == 'bh':
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP * n / (sorted_pvalues.index(eachP) + 1) for eachP in pvalue_list]
        elif methods == 'by':
            _, fdr = fdrcorrection(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
        fdr = np.array(fdr)
        fdr[fdr > 1] = 1.
        return fdr