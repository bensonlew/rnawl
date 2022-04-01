# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import pandas as pd
import os
import glob
from biocluster.file import getsize, exists
from biocluster.file import download
from mbio.packages.ref_rna_v2.functions import tryforgood
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import glob
from mbio.packages.denovo_rna_v2.chart import Chart


class DiffexpBatchNewWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpBatchNewWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="string",),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            dict(name='exp_type', type='string', default="tpm"),
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="filter_method", type="string", default=None),  # filter_method 20190801修改
            dict(name="tpm_filter_threshold", type="float", default="0"), #20190708 添加 by fwy
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),
            # DESeq2 DE test method, Wald|LRT
            dict(name="deseq2_method", type="string", default="Wald"),
            # edger_method DE test method,exactTest|glmLRT|glmQLFTest
            dict(name="edger_method", type="string", default="glmQLFTest"),
            # degseq_method DE test method, LRT|CTR|FET|MARS|MATR|FC
            dict(name="degseq_method", type="string", default="MARS"),
            dict(name="task_id", type="string", default=""),
            dict(name="exp_level", type="string", default=""),
            # batch by zjx 20200917
            dict(name="is_batch", type="bool", default=False),
            dict(name="has_batch", type="bool"),
            dict(name="batch_matrix", type="infile", format="denovo_rna_v2.common"),
            # noiseq
            dict(name='prob', type='float', default=0.8),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("denovo_rna_v3.batch.diffexp_batch")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 Diff_Express')
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
        super(DiffexpBatchNewWorkflow, self).send_log(data)

    def run(self):
        self.tool.on('end', self.run_uniform)
        self.get_run_log()
        self.run_tool()
        super(DiffexpBatchNewWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_diff", main_id=self.option('diff_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.stop_timeout_check()
        self.all_exp = self.api.api("denovo_rna_v2.all_exp")
        # add result info
        self.all_exp.add_diffexp_all(self.uniform.output_dir, self.tool.output_dir,
                                     group_dict=self.option('group_dict'),
                                     main_id=self.option('diff_main_id'),
                                     diff_method=self.option('method'),
                                     create_geneset=False,
                                     pvalue_padjust=self.option('pvalue_padjust')
                                     )

        self.paste_annotation()
        self.set_output()
        self.end()

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["03 Diff_Express", "", "表达量差异分析结果目录",0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达量差异分析文件", 0],
        ])
        result_dir.add_regexp_rules([
            [r'.*_vs_.*\.xls', 'xls', '表达量差异结果表',0],
            [r'.*summary.*\.xls', 'xls', '表达量差异统计表',0],
            [r'.*total_diff_stat.*\.xls', 'xls', '表达量差异详情总表',0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['.*scatter.pdf', 'pdf', '表达量差异散点图', 0],
            ['.*volcano.pdf', 'pdf', '表达量差异火山图', 0],
            ['.*bar_h.*.pdf', 'pdf', '表达量差异统计柱形图', 0],
            ['.*bar_v.*.pdf', 'pdf', '表达量差异统计堆叠图', 0],
        ])
        super(DiffexpBatchNewWorkflow, self).end()

    def set_output(self):
        diff_total=pd.read_table(os.path.join(self.tool.output_dir,"total_diff_stat.{}.xls".format(self.option("method").lower())),index_col="seq_id")
        diff_files_stat = glob.glob(os.path.join(self.tool.output_dir, '*_vs_*.{}.xls'.format(self.option("method").lower())))
        df_anno = pd.read_table(diff_files_stat[0], index_col="seq_id")
        df1_genebase = pd.DataFrame(df_anno, columns=["length"])
        df1_geneanno=pd.DataFrame(df_anno,columns=["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot"])
        dff1 = pd.concat([df1_genebase, diff_total], axis=1)
        dff2 = pd.concat([dff1, df1_geneanno], axis=1)
        dff2.index.set_names("seq_id", inplace=True)
        dff2.to_csv(self.tool.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), sep="\t")
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))

    def run_tool(self):
        opts = dict(
            exp=self.option('exp_matrix'),
            exp_type=self.option('exp_type'),
            count=self.option('count'),
            group=self.option('group'),
            method=self.option('method'),
            cmp=self.option('cmp'),
            fc=self.option('fc'),
            tpm_filter_threshold=self.option("tpm_filter_threshold"),
            is_batch=self.option('is_batch'),
            filter_method=self.option("filter_method"),
            no_filter=True,
            analysis='split'
        )

        if self.option('is_batch') == False:
            if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                opts.update(
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                )
            else:
                opts.update(
                    prob=self.option('prob')
                )
        else:
            if self.option('has_batch') == True:
                opts = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=self.option('count'),
                    group=self.option('group'),
                    method=self.option('method'),
                    cmp=self.option('cmp'),
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                    fc=self.option('fc'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                    batch_matrix=self.option('batch_matrix'),
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
            else:
                opts = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=self.option('count'),
                    group=self.option('group'),
                    method=self.option('method'),
                    cmp=self.option('cmp'),
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                    fc=self.option('fc'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
        self.tool.set_options(opts)
        self.tool.run()

    def run_uniform(self):
        self.uniform = self.add_tool('denovo_rna_v3.batch.uniform')
        opts = {
            'method': self.option('method'),
            'input_dir': self.tool.output_dir,
        }
        if self.option('method').lower() in ['deseq2', 'degseq', 'edger', 'limma']:
            opts.update({'pvalue_padjust': self.option('pvalue_padjust')})
        self.uniform.set_options(opts)
        self.uniform.on('end', self.set_db)
        self.uniform.run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        diff_summary = os.path.join(self.uniform.output_dir, 'json')
        diff_volcano = os.path.join(self.uniform.output_dir, 'all_volcano.txt')
        diff_scatter = os.path.join(self.uniform.output_dir, 'all_scatter.txt')
        chart.denovo_chart_diff_summary(diff_summary, self.option('exp_level'))
        chart.denovo_chart_diff_volcano(diff_volcano, self.option('exp_level'))
        chart.denovo_chart_diff_scatter(diff_scatter, self.option('exp_level'))
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.tool.output_dir + "/" + os.path.basename(p))

    @tryforgood
    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '12003001')
        return to_path


    def paste_annotation(self):
        conn = self.all_exp.db['sg_annotation_stat']
        task_id = self.option('task_id')
        find_result = conn.find_one({"task_id": task_id, "type": "origin"})

        exp_level = self.option('exp_level')
        if exp_level:
            self.logger.info("we find exp_level")
        else:
            self.logger.info("exp_level info miss")
        result_dir=find_result['result_dir']
        if exists(os.path.join(result_dir, 'all_annot.xls')):
            # remote_annot = os.path.join(result_dir, 'anno_stat/gene_anno_detail.xls')
            try:
                remote_annot = os.path.join(result_dir, 'all_annot.xls')
                annot = self.download_s3_file(remote_annot, "all_annot.xls")
            except:
                annot=os.path.join(self.work_dir,"all_annot.xls")
                os.link(os.path.join(result_dir, 'all_annot.xls'),annot)

            all_annot = pd.read_table(annot, header=0)
            if exp_level[0].upper() == 'G':
                annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript', 'is_gene']).set_index(
                    'gene_id')
            else:
                annot_pd = all_annot.reset_index().drop(columns=['is_gene']).set_index('transcript')
        else:
            if exp_level[0].upper() == 'G':
                if exists(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'unigene_anno_detail.xls')):
                    remote_annot = os.path.join(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'unigene_anno_detail.xls'))
                    annot = self.download_s3_file(remote_annot, "gene_anno_detail.xls")
                    annot_pd = pd.read_table(annot, header=0,index_col="gene_id")
                else:
                    self.set_error('annotation result file cannot found', code="12003002")
            else:
                if exists(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'transcript_anno_detail.xls')):
                    remote_annot = os.path.join(os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'transcript_anno_detail.xls'))
                    annot = self.download_s3_file(remote_annot, "trans_anno_detail.xls")
                    all_annot = pd.read_table(annot, header=0)
                    annot_pd = all_annot.reset_index().drop(columns=['gene_id']).set_index('transcript')
                else:
                    self.set_error('annotation result file cannot found', code="12003003")
        self.logger.info("annot dir is {}".format(annot))

        for each in glob.glob(os.path.join(self.tool.output_dir, '*_vs_*.*.xls')):
            if each.endswith('.annot.xls'):
                continue
            if each.endswith('.normalize.xls'):
                continue
            if each.endswith('.sizeFactor.xls'):
                continue
            if each.endswith('.normFactor.xls'):
                continue
            diff_pd = pd.read_table(each, header=0, sep='\t', index_col="seq_id")
            ctrl, test = os.path.basename(each).split('.{}.xls'.format(self.option("method").lower()))[0].split('_vs_')
            if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                b = pd.DataFrame(diff_pd,columns = ["fc","log2fc", "pvalue", "padjust", "significant", "regulate"])
                c=diff_pd.drop(["fc","log2fc", "pvalue", "padjust", "significant", "regulate"], axis=1)
            else:
                try:
                    b = pd.DataFrame(diff_pd,columns = ['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "theta", "D", 'prob', "significant", "regulate"])
                    c = diff_pd.drop(['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "theta", "D", 'prob', "significant", "regulate"], axis=1)
                except:
                    b = pd.DataFrame(diff_pd,columns = ['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "D", 'prob', "significant", "regulate"])
                    c = diff_pd.drop(['{}_mean'.format(ctrl), '{}_mean'.format(test),"fc","log2fc", "D", 'prob', "significant", "regulate"], axis=1)
            result = pd.concat([b, c], axis=1)
            order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot"]
            order2 = ["length"]
            pd1 = pd.DataFrame(annot_pd, columns=order)
            pd2 = pd.DataFrame(annot_pd, columns=order2)
            pd3 = pd.concat([pd2, result], axis=1, join_axes=[result.index])
            pd4 = pd.concat([pd3, pd1], axis=1, join_axes=[pd3.index])
            pd4.to_csv(each, sep='\t', header=True, index=True)

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.denovo_rna_v3.report.diffexp_batch import DiffexpBatchWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'Diff_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'denovo_rna_v3.report.diffexp_batch',
            'options': dict(
                count="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/unigene.count.matrix.xls",
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/exp_matrix",
                method="DESeq2",
                group="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/group",
                cmp="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/cmp_little",
                padjust_way="BH",
                pvalue=0.05,
                fc=1,
                is_batch=False,
                exp_level='G',


            )
        }
        wsheet = Sheet(data=data)
        wf =DiffexpBatchWorkflow(wsheet)
        wf.sheet.id = 'diff'
        wf.sheet.project_sn = 'diff'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


