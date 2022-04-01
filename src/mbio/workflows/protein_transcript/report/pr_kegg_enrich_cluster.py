
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import glob
import os
import pandas as pd
import numpy as np
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json

class PrKeggEnrichClusterWorkflow(Workflow):
    """
    蛋白集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(PrKeggEnrichClusterWorkflow, self).__init__(wsheet_object)
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
            {"name": "n_clusters", "type": "int"},
            {"name": "gcm", "type": "string"},
            {"name": "gct", "type": "string"},
            {"name": "gcd", "type": "string"},
            {"name": "cluster_col", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.protein_tools = list()
        self.rna_tools = list()
        self.cluster_tool = self.add_tool("protein_transcript.enrich_cluster")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/7_relate/03_relate_anno/06_relate_keggenrich_cluster')
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
        super(PrKeggEnrichClusterWorkflow, self).send_log(data)

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

    def run_cluster(self):
        def change_qvalue(file):
            pf = pd.read_table(file, header=0, sep="\t")
            pf['p_corrected'] = self.multtest_correct(pf['p_uncorrected'].tolist(), self.option("method").lower())
            pf.drop(pf.columns[[9]], axis=1, inplace=True)
            os.remove(file)
            pf.to_csv(file, sep='\t', index=False)
        def cal_col(row):
            a,b = row.split('/')
            return float(a) / float(b)
        col_list = list()
        info_list_tmp = list()
        for tool in self.protein_tools + self.rna_tools:
            enrich_file = glob.glob("{}/*.xls".format(tool.output_dir))[0]
            # change_qvalue(enrich_file)
            target = os.path.join(self.output_dir, os.path.basename(enrich_file))
            if os.path.exists(target):
                os.remove(target)
            os.link(enrich_file, target)
            diff_path = tool.option('diff_list').prop['path']
            # print(diff_path)
            if u'_protein.list' in diff_path:
                name = os.path.basename(diff_path.split('_protein.list')[0])
            else:
                name = os.path.basename(diff_path.split('_gene.list')[0])
            enrich_df = pd.read_table(enrich_file, sep= '\t', header=0,)
            # 由于kegg分map和单物种，所以需要强行转一下
            enrich_df['ID'] = enrich_df['ID'].apply(lambda x: 'map' + filter(str.isdigit, x))
            enrich_df_col = enrich_df.set_index('ID')
            if self.option('cluster_col').lower() == 'ratio_in_study':
                enrich_df_col = enrich_df_col['Ratio_in_study']
                # enrich_df_col['ratio_in_study_'] = enrich_df_col.apply(cal_col)
                # enrich_df_col.drop(enrich_df_col.columns[[0]], axis=1, inplace=True)
                enrich_df_col = enrich_df_col.apply(cal_col)
            if self.option('cluster_col').lower() == 'ratio_in_pop':
                enrich_df_col = enrich_df_col['Ratio_in_pop']
                enrich_df_col = enrich_df_col.apply(cal_col)
            if self.option('cluster_col').lower() == 'p_value':
                enrich_df_col = enrich_df_col['P-Value']
            if self.option('cluster_col').lower() == 'q_value':
                enrich_df_col = enrich_df_col['Corrected P-Value']
            enrich_df_col.name = name
            enrich_df_col.index.name = 'kegg'
            col_list.append(enrich_df_col)
            enrich_df = enrich_df[enrich_df.columns.tolist()[0:4]]
            info_list_tmp += enrich_df.to_dict('record')
        self.info_list = list()
        kegg_list = set()
        for dict_ in info_list_tmp:
            if dict_['ID'] not in kegg_list:
                self.info_list.append(dict_)
                kegg_list.add(dict_['ID'])
        cluster_pd = pd.concat(col_list, axis=1)
        cluster_pd.index.name = 'kegg'
        if self.option('cluster_col').lower() == 'ratio_in_study' or self.option('cluster_col').lower() == 'ratio_in_pop':
            cluster_pd = cluster_pd.fillna(0)
        else:
            cluster_pd = cluster_pd.fillna(1)
        self.cluster_file = os.path.join(self.output_dir, 'col_cluster.xls')
        cluster_pd.to_csv(self.cluster_file, sep="\t", index=True)
        self.cluster_tool.on("end", self.set_db)
        options = dict(
            cluster_file=self.cluster_file,
            n_clusters=int(self.option('n_clusters')),
            sct='no',
            type_='kegg',
            gct=self.option('gct'),
            gcm=self.option('gcm'),
            gcd=self.option('gcd'),
        )
        self.cluster_tool.set_options(options)
        self.cluster_tool.run()

    def run(self):
        self.run_tools()
        self.on_rely(self.protein_tools + self.rna_tools, self.run_cluster)
        for tool in self.protein_tools + self.rna_tools:
            tool.run()
        super(PrKeggEnrichClusterWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_keggenrich = self.api.api('protein_transcript.relaset')
        cluster_file = os.path.join(self.cluster_tool.output_dir, 'cluster_matrix.xls')
        df_cluster = pd.read_table(cluster_file, sep = '\t', header = 0, index_col=0)
        df_info = pd.DataFrame(self.info_list)
        df_info.rename(columns={'#Study_num': "study_num", "Term":"term",
                                "Database": "database", "ID": "kegg"}, inplace=True)
        df_info = df_info.set_index('kegg')
        # df_cluster = pd.concat([df_info, df_cluster], axis=1)
        df_cluster = df_cluster.join(df_info, how='inner')
        df_cluster.index.name = 'kegg'
        os.remove(cluster_file)
        df_cluster.to_csv(cluster_file, sep="\t", index=True)
        target = os.path.join(self.output_dir, os.path.basename(cluster_file))
        if os.path.exists(target):
            os.remove(target)
        os.link(cluster_file, target)
        api_keggenrich.add_enrich_cluster(self.cluster_tool.output_dir, main_id=self.option('main_table_id'), type_ = 'kegg')
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        pt_keggenrich_cluster_matrix = os.path.join(self.cluster_tool.output_dir, "cluster_matrix.xls")
        pt_keggenrich_cluster_tree = os.path.join(self.cluster_tool.output_dir, "seq.cluster_tree.txt")
        if os.path.exists(pt_keggenrich_cluster_matrix) and os.path.exists(pt_keggenrich_cluster_tree):
            chart.chart_pt_keggenrich_cluster(pt_keggenrich_cluster_matrix, pt_keggenrich_cluster_tree)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = os.path.join(self.work_dir, "pt_keggenrich_cluster.heatmap_new.pdf")
            if os.path.exists(pdf_file):
                os.link(pdf_file, os.path.join(self.output_dir, "heat.pdf"))
            pdf_files = glob.glob(self.work_dir + "/pt_keggenrich_cluster__*.pdf")
            for pdf_file in pdf_files:
                os.link(pdf_file, os.path.join(self.output_dir, "sub"+os.path.basename(pdf_file)[23:].split('.')[0].split('_')[1]+'.pdf'))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["7_relate", "", "蛋白组与转录组关联分析",0],
            ["7_relate/03_relate_anno", "", "功能注释信息", 0],
            ["7_relate/03_relate_anno/06_relate_keggenrich_cluster", "", "KEGG富集分类", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "多个蛋白集和基因GO富集及其聚类分析结果文件"],
            ["./heat.pdf", "", "热图"],
            ["./*sub*.pdf", "", "子聚类图"],
        ])
        super(PrKeggEnrichClusterWorkflow, self).end()


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