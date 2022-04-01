
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# go的pvalue矫正改为跟kegg类似的函数--by 封一统 20180713

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import re
from bson.objectid import ObjectId
import types
import pandas as pd
import numpy as np
import tarfile
import shutil
from mbio.packages.prok_rna.chart import Chart


class GenesetEnrichWorkflow(Workflow):
    """
    基因集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "geneset_kegg", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "kegg_table", "type": "string"},
            {"name": "go_list", "type": "string"},
            {'name': 'go_version', 'type': 'string', 'default': '2018'},
            {"name": "geneset_list", "type": "string"},
            {"name": "all_list", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "type", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {'name': 'source', 'type': 'string', 'default': None},
            {"name": "add_info", "type": "string", "default": None},  # 输入两列的列表文件，有head，第一列为pathway，第二列为底图链接
            {"name": "geneset_name", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "go":
            self.enrich_tool = self.add_tool("prok_rna.geneset.go_enrich")
        else:
            self.enrich_tool = self.add_tool("prok_rna.geneset.kegg_rich")
        self.kegg_class = self.add_tool("prok_rna.geneset.kegg_class")
        self.output_dir1 = self.enrich_tool.output_dir
        self.output_dir2 = self.kegg_class.output_dir

    def run(self):
        if self.option("anno_type") == "kegg":
            options = {
                "kegg_table": self.option("kegg_table"),
                # "all_list": background_path,
                "diff_list": self.option("geneset_list"),
                "correct": self.option("method"),
                "add_info": self.option("add_info"),
                "kegg_version": self.option("kegg_version")
            }
        else:
            options = {
                "diff_list": self.option("geneset_list"),
                # "all_list": background_path,
                "go_list": self.option("go_list"),
                'go_version': self.option('go_version'),
                # "pval": self.option("pval"),
                "method": self.option("method"),
            }
        self.logger.info(options)
        self.enrich_tool.set_options(options)
        if self.option("anno_type") == "kegg":
            self.enrich_tool.on('end', self.run_kegg_class)
            self.kegg_class.on('end', self.set_db)
        else:
            self.enrich_tool.on('end', self.set_db)
        self.enrich_tool.run()
        super(GenesetEnrichWorkflow, self).run()

    def replace_kegg_link(self):
        '''
        替换kegg链接
        '''
        enrich_result = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        annot_result = self.kegg_class.output_dir + '/kegg_stat.xls'
        with open(annot_result, 'rb') as f:
            map2link = {line.split("\t")[0]:line.strip().split("\t")[-1] for line in f.readlines()[1:] if line.split("\t")[-1].startswith("http")}
        with open(enrich_result, 'rb') as f, open(enrich_result + "relink", 'w') as fo:
            header = f.readline()
            fo.write(header)
            for line in f:
                cols = line.split("\t")
                if cols[3] in map2link:
                    cols[9] = map2link[cols[3]]
                fo.write("\t".join(cols))
        os.remove(enrich_result)
        os.link(enrich_result + "relink", enrich_result)


    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('prok_rna.geneset')
        output_file = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))
        workflow_output = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))
        workflow_output = self.get_workflow_output_dir() + '/' + workflow_output[0].split('/')[-1]
        # png_file = glob.glob("{}/*.png".format(self.output_dir))
        # go_png = self.output_dir + "/go_lineage.png"
        # go_pdf = self.output_dir + "/go_lineage.pdf"
        go_adjust_png = self.enrich_tool.output_dir + "/adjust_lineage.png"
        go_adjust_pdf = self.enrich_tool.output_dir + "/adjust_lineage.pdf"
        go_adjust_svg = self.enrich_tool.output_dir + "/adjust_lineage.svg"
        if self.option("anno_type") == "kegg":
            if self.option("source") == "diff_exp":
                self.replace_kegg_link()
            api_geneset.add_kegg_enrich_detail(self.option("main_table_id"), output_file[0])
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option("main_table_id"),
                                         result_dir=workflow_output)
            graph_dir = os.path.join(self.get_workflow_output_dir(), 'pathways')
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option("main_table_id"),
                                         graph_dir=graph_dir)
            kegg_stat = self.kegg_class.output_dir + '/kegg_stat.xls'
            # pathway_class = self.option("kegg_table_2").prop['path']
            pathway_file = self.kegg_class.output_dir + '/pathways'

            self.logger.info("开始进行kegg_class的导表")
            # api_geneset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
            # api_geneset.add_kegg_regulate_new(self.option("main_table_id"), self.option("geneset_id"), output_file, self.option("kegg_table_2").prop['path'], self.work_dir)
            # api_geneset.add_kegg_regulate_new2(self.option("main_table_id"), self.option("geneset_kegg").prop['path'],
            #                                    kegg_stat, self.option("kegg_table_2").prop['path'])
            api_geneset.add_kegg_enrich_pic(self.option("main_table_id"), kegg_stat,
                                            pathway_file, source=self.option("source"))
            pngs = os.listdir(pathway_file)
            tar_file = pathway_file + ".tar.gz"
            with tarfile.open(tar_file, mode='w:gz') as f:
                for png in pngs:
                    f.add(pathway_file + "/" + png, arcname=png)
            shutil.rmtree(pathway_file)
        else:
            pf = pd.read_table(output_file[0], header=0, sep="\t")
            pf['p_corrected'] = self.multtest_correct(pf['p_uncorrected'].tolist(), self.option("method").lower())
            if not pf.empty:  # allow empty result
                pf.drop([9])
            os.remove(output_file[0])
            pf.to_csv(output_file[0], sep='\t', index=False)
            api_geneset.add_go_enrich_detail(self.option("main_table_id"), output_file[0])
            api_geneset.update_directed_graph(self.option("main_table_id"), go_adjust_png, go_adjust_pdf)
            conn = api_geneset.db["sg_geneset_go_enrich"]
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
                self.set_error("main_id参数必须为字符串或者ObjectId类型", code = "15000201")
            conn.update({"_id": record_id}, {"$set": {'result_dir': graph_dir}}, upsert=True)

        self.end()

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        if self.option("anno_type") == "go":
            geneset_go_enrich = os.path.join(self.output_dir1, 'go_enrich_geneset_list_gene.xls')
            if os.path.exists(geneset_go_enrich):
                chart.chart_geneset_gorich(geneset_go_enrich, self.option('geneset_name'))
        elif self.option("anno_type") == "kegg":
            geneset_kegg_enrich = os.path.join(self.output_dir1, 'geneset_list_gene.list.DE.list.check.kegg_enrichment.xls')
            if os.path.exists(geneset_kegg_enrich):
                chart.chart_geneset_keggrich(geneset_kegg_enrich, self.option('geneset_name'))
        chart.to_pdf()

        # move pdfs
        if self.option("anno_type") == "go":
            for i in [['*_all_go_enrich_dense_bubble.densebubble.pdf', 'go_enrich_bubble.pdf'],
                      ['*_all_go_enrich_scatter_bubble.scatterbubble.pdf', 'go_enrich_disbubble.pdf'],
                      ['*_go_enrich_bar.shadowbar.pdf', 'go_enrich_bar.pdf']]:
                pdf = glob.glob(os.path.join(self.work_dir, i[0]))[0]
                self.move_pdf(pdf, os.path.join(self.output_dir1, i[1]))

        elif self.option("anno_type") == "kegg":
            for i in [['*_enrichkegg.densebubble.pdf', 'pathway_enrich_bubble.pdf'],
                      ['*_enrichkegg.shadowbar.pdf', 'pathway_enrich_bar.pdf']]:
                pdf = glob.glob(os.path.join(self.work_dir, i[0]))[0]
                self.move_pdf(pdf, os.path.join(self.output_dir1, i[1]))

    def end(self):
        self.chart()

        result_dir = self.add_upload_dir(self.output_dir1)
        if self.option("anno_type") == "go":
            result_dir.add_relpath_rules([
                [".", "", "基因集GO富集分析结果文件"],
                ["adjust_lineage.svg", "", "显著富集GO的层级结构图"],
                ["adjust_lineage.png", "", "显著富集GO的层级结构图"],
                ["go_lineage.svg", "", "GO的层级结构图"],
                ["go_enrich_geneset_list_gene.xls", "", "GO富集分析统计表"],
                ["go_lineage.pdf", "", "GO的层级结构图"],
                ["adjust_lineage.pdf", "", "显著富集GO的层级结构图"],
                ["go_lineage.png", "", "GO的层级结构图"],
                ["go_enrich_bar.pdf", "", "GO富集分析分类统计柱形图"],
                ["go_enrich_bubble.pdf", "", "GO富集分析分类统计气泡图"],
                ["go_enrich_disbubble.pdf", "", "GO富集分析分类统计分散型气泡图"],
            ])
        elif self.option("anno_type") == "kegg":
            if os.path.exists(os.path.join(self.output_dir2, 'pathways.tar.gz')):
                os.link(os.path.join(self.output_dir2, 'pathways.tar.gz'), os.path.join(self.output_dir1, 'pathways.tar.gz'))
            if os.path.exists(os.path.join(self.output_dir2, 'ko')):
                cmd = 'cp -r {} {}'.format(os.path.join(self.output_dir2, 'ko'), os.path.join(self.output_dir1, 'ko'))
                # os.link(os.path.join(self.output_dir2, 'ko'), os.path.join(self.output_dir1, 'ko'))
                os.system(cmd)
            if os.path.exists(os.path.join(self.output_dir2, 'kegg_stat.xls')):
                os.link(os.path.join(self.output_dir2, 'kegg_stat.xls'), os.path.join(self.output_dir1, 'kegg_stat.xls'))
            result_dir.add_relpath_rules([
                [".", "", "基因集KEGG富集分析结果文件"],
                ["ko", "", "Pathway通路图文件"],
                ["kegg_stat.xls", "", "Pathway通路统计表"],
                ["geneset_list_gene.list.DE.list.check.kegg_enrichment.xls", "", "kegg富集结果表"],
                ["pathways.tar.gz", " ", "Pathway通路图压缩文件"],
                ["pathway_enrich_bar.pdf", "", "KEGG富集分析分类统计柱形图"],
                ["pathway_enrich_bubble.pdf", "", "KEGG富集分析分类统计气泡图"],
            ])
            # result_dir2 = self.add_upload_dir(self.output_dir2)
            # result_dir2.add_relpath_rules([
            #     ["pathways.tar.gz", " ", "Pathway通路图压缩文件"],
            # ])
        super(GenesetEnrichWorkflow, self).end()


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

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

    def run_kegg_class(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg"),
            "kegg_table": self.option("kegg_table"),
            "geneset_id": self.option("geneset_id"),
            "background_links": self.option("add_info"),
            "type": self.option("type"),
            "task_id": self.option("task_id"),
            'source': self.option('source'),
            "kegg_version": self.option('kegg_version')
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()
