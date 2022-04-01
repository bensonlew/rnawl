
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
from bson.objectid import ObjectId
import types
import re
import pandas as pd
import numpy as np
import json
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class ProteinsetEnrichWorkflow(Workflow):
    """
    蛋白集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(ProteinsetEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "kegg_table", "type": "string"},
            {"name": "go_list", "type": "string"},
            {"name": "proteinset_list", "type": "string"},
            {"name": "proteinset_id", "type": "string"},
            # {"name": "all_list", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "type", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': "2017"},
            {"name": "go_version", "type": "string", "default": "2018"}, #pir database version
            {"name": "add_info", "type": "string", "default": None},  # 输入两列的列表文件，有head，第一列为pathway，第二列为底图链接
            {"name": "proteinset_kegg", "type": "infile", "format": "itraq_and_tmt.common"},
            {"name": "proteinset_name", 'type': 'string'},  # added for chart
            {'name': 'task_version', 'type': 'string', 'default': "2"},
            
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.enrich_tool = self.add_tool("labelfree.proteinset.go_enrich") if self.option("anno_type") == "go" else self.add_tool("dia_v3.proteinset.kegg_rich")
        self.output_dir = self.enrich_tool.output_dir
        if self.option('anno_type') == 'go':
            self._sheet.output = self._sheet.output.replace('interaction_results',
                                                            'interaction_results/5_Proteinset/04_Enrich/01_GO')
        else:
            self._sheet.output = self._sheet.output.replace('interaction_results',
                                                            'interaction_results/5_Proteinset/04_Enrich/02_KEGG')
        self.inter_dirs = []

    def run(self):
        if self.option("anno_type") == "kegg":
            options = {
                "kegg_table": self.option("kegg_table"),
                # "all_list": background_path,
                "diff_list": self.option("proteinset_list"),
                "proteinset_kegg": self.option("proteinset_kegg").prop['path'],
                "proteinset_id": self.option("proteinset_id"),
                "correct": self.option("method"),
                "add_info": self.option("add_info"),
                "task_id": self.option("task_id"),
                "kegg_version": self.option('kegg_version')
            }
        else:
            options = {
                "diff_list": self.option("proteinset_list"),
                # "all_list": background_path,
                "go_list": self.option("go_list"),
                "go_version": self.option('go_version'),
                # "pval": self.option("pval"),
                "method": self.option("method"),
            }
        self.logger.info(options)
        self.enrich_tool.set_options(options)
        self.enrich_tool.on('end', self.set_db)
        self.enrich_tool.run()
        super(ProteinsetEnrichWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        self.output_dir = self.enrich_tool.output_dir
        api_proteinset = self.api.api('dia.proteinset')
        output_file = glob.glob("{}/*.xls".format(self.output_dir))
        # png_file = glob.glob("{}/*.png".format(self.output_dir))
        # go_png = self.output_dir + "/go_lineage.png"
        # go_pdf = self.output_dir + "/go_lineage.pdf"
        go_adjust_png = self.output_dir + "/adjust_lineage.png"
        go_adjust_pdf = self.output_dir + "/adjust_lineage.pdf"
        if self.option("anno_type") == "kegg":
            kegg_enrich_id =  api_proteinset.add_kegg_enrich_detail(self.option("main_table_id"), output_file[0])
            api_proteinset.add_kegg_enrich_pic(self.option("main_table_id"),  output_file[0], self.enrich_tool.output_dir + '/pathways')
        else:
            pf = pd.read_table(output_file[0], header=0, sep="\t")
            pf['p_corrected'] = self.multtest_correct(pf['p_uncorrected'].tolist(), self.option("method").lower())
            pf = pf.drop(pf.columns.tolist()[9], axis=1)
            os.remove(output_file[0])
            pf.to_csv(output_file[0], sep='\t', index=False)
            api_proteinset.add_go_enrich_detail(self.option("main_table_id"), output_file[0])
            api_proteinset.update_directed_graph(self.option("main_table_id"), go_adjust_png, go_adjust_pdf)
            conn = api_proteinset.db["sg_proteinset_go_enrich"]
            self.workflow_output_tmp = self._sheet.output
            if re.match(r'tsanger:', self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
            else:
                self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
            graph_dir = self.workflow_output
            # xls_dir = glob.glob("{}/*.xls".format(graph_dir))[0]
            # os.remove(xls_dir)
            # os.link(graph_dir, output_file[0])
            record_id = self.option("main_table_id")
            if isinstance(record_id, types.StringTypes):
                record_id = ObjectId(record_id)
            elif isinstance(record_id, ObjectId):
                record_id = record_id
            else:
                raise Exception("main_id参数必须为字符串或者ObjectId类型!")
            conn.update({"_id": record_id}, {"$set": {'result_dir': graph_dir}}, upsert=True)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        if self.option("anno_type") == "go":
            go_enrich_all_file = os.path.join(self.enrich_tool.output_dir, "go_enrich_proteinset_list_protein.xls")
            chart.chart_proteinsetcluster_gorich([''], [go_enrich_all_file], [go_enrich_all_file], [go_enrich_all_file], self.option('proteinset_name'))
            chart.to_pdf()
            for old, new in [["*densebubble*pdf", 'enrichgo.bubble2.pdf'],
                             ["*scatterbubble*pdf", 'enrichgo.bubble.pdf'],
                             ["*shadowbar*pdf", 'enrichgo.bar.pdf']]:
                files = glob.glob(os.path.join(self.work_dir, old))
                for each in files:
                    os.link(each, os.path.join(self.enrich_tool.output_dir, new))
        if self.option("anno_type") == "kegg":
            kegg_enrich_all_file = os.path.join(self.enrich_tool.output_dir, "proteinset_list_protein.list.DE.list.check.kegg_enrichment.xls")
            chart.chart_proteinsetcluster_keggrich_web(kegg_enrich_all_file, self.option('proteinset_name'))
            chart.to_pdf()
            for old, new in [["enrichkegg.densebubble.pdf", 'enrichkegg.bubble2.pdf'],
                             ["enrichkegg.scatterbubble.pdf", 'enrichkegg.bubble.pdf'],
                             ["enrichkegg.shadowbar.pdf", 'enrichkegg.bar.pdf']]:
                if os.path.exists(os.path.join(self.work_dir, old)):
                    os.link(os.path.join(self.work_dir, old), os.path.join(self.enrich_tool.output_dir, new))

    def end(self):
        self.chart()
        rm_files = glob.glob(os.path.join(self.enrich_tool.output_dir, 'pathways/*.png')) + \
                   glob.glob(os.path.join(self.enrich_tool.output_dir, 'pathways/*.pdf'))
        if self.option("task_version") >= "2.1":
            for rm_file in rm_files:
                os.remove(rm_file)
        
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "go":
            self.inter_dirs = [
                ["5_Proteinset", "", "蛋白集分析", 0],
                ["5_Proteinset/04_Enrich/", "", "功能富集", 0],
                ["5_Proteinset/04_Enrich/01_GO", "", "go富集", 0],
            ]
            result_dir.add_relpath_rules([
                [".", "", "蛋白集GO富集分析结果文件", 0, "220076"],
                ["enrichgo.bar.pdf", "", "GO富集分析直方图", 0, "220076"],
                ["enrichgo.bubble.pdf", "", "GO富集分析气泡图", 0, "220076"],
                ["enrichgo.bubble2.pdf", "", "GO富集分析气泡图", 0, "220076"],
            ])
        elif self.option("anno_type") == "kegg":
            self.inter_dirs = [
                ["5_Proteinset", "", "蛋白集分析", 0],
                ["5_Proteinset/04_Enrich/", "", "功能富集", 0],
                ["5_Proteinset/04_Enrich/02_KEGG", "", "kegg富集", 0],
            ]
            result_dir.add_relpath_rules([
                [".", "", "蛋白集KEGG富集分析结果文件", 0, "220077"],
                ["enrichkegg.bar.pdf", "", "KEGG富集分析直方图", 0, "220076"],
                ["enrichkegg.bubble.pdf", "", "KEGG富集分析气泡图", 0, "220076"],
                ["enrichkegg.bubble2.pdf", "", "KEGG富集分析气泡图", 0, "220076"],
            ])
        super(ProteinsetEnrichWorkflow, self).end()

    def multtest_correct(self, p_values, methods=3):  # 复制于package下的kegg_rich.py
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
