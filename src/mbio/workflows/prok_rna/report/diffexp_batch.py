# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
import glob
import os
import json
import time
from biocluster.file import getsize, exists
from biocluster.file import download
from mbio.packages.ref_rna_v2.functions import tryforgood
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.prok_rna.chart import Chart


class DiffexpBatchWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpBatchWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="infile", format="prok_rna.common",),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            dict(name='exp_type', type='string', default='tpm'),
            # dict(name="type", type="string"),
            dict(name="task_id", type="string"),
            dict(name="exp_level", type="string"),
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="filter_method", type="string", default=None),  # filter_method 20190801修改
            dict(name="tpm_filter_threshold", type="float", default="0"), #20190708 添加 by fwy
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),
            # batch by zjx 20200917
            dict(name="is_batch", type="bool", default=False),
            dict(name="has_batch", type="bool"),
            dict(name="batch_matrix", type="infile", format="prok_rna.common"),
            # noiseq
            dict(name='prob', type='float', default=0.8),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("prok_rna.diffexp_batch")
        self.all_exp = self.api.api("prok_rna.all_exp")

    def run(self):
        self.tool.on('end', self.run_uniform)
        self.run_tool()
        super(DiffexpBatchWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # add result info
        self.stop_timeout_check()
        self.all_exp = self.api.api("prok_rna.all_exp")
        # add result info
        self.all_exp.add_diffexp_all(self.uniform.output_dir, self.tool.output_dir,
                                     group_dict=self.option('group_dict'),
                                     main_id=self.option('diff_main_id'),
                                     diff_method=self.option('method'),
                                     create_geneset=False,
                                     pvalue_padjust=self.option('pvalue_padjust')
                                     )
        # # self.all_exp.add_diffexp(self.tool.output_dir,
        # #                          main_id=self.option('diff_main_id'),
        # #                          diff_method=self.option('method'),
        # #                          create_geneset=False,
        # #                          pvalue_padjust=self.option('pvalue_padjust')
        # #                     )
        self.paste_annotation()
        self.end()

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        gene_diff_summary = os.path.join(self.tool.output_dir, 'diff_summary_{}.xls'.format(self.option('method').lower()))
        cmp_list = list()
        with open(self.option('cmp'), 'r') as f:
            f.readline()
            for line in f:
                cmp = line.strip().split('\t')
                cmp_list.append([cmp[0], cmp[1]])
                gene_diff_scatter = os.path.join(self.tool.work_dir, '{}_vs_{}.{}.xls'.format(cmp[0], cmp[1], self.option('method').lower()))
                chart.prok_diffexp_scatter(gene_diff_scatter, '_vs_'.join(cmp), soft=self.option('method'))
        chart.chart_diffexp_stat(gene_diff_summary, cmp_list=cmp_list, soft=self.option('method'))
        chart.to_pdf()
        # move pdf
        target_dir = self.tool.output_dir
        volcano_pdfs = glob.glob(os.path.join(self.work_dir, "*.diffexp.*.pdf"))
        for each in volcano_pdfs:
            cmp, _, suffix = os.path.basename(each).split('.', 2)
            new_path = os.path.join(target_dir, cmp + '.' + suffix)
            self.move_pdf(each, new_path)

        bar_pdf = glob.glob(os.path.join(self.work_dir, '*.diffexp_summary.bar.pdf'))
        if bar_pdf:
            new_path = os.path.join(target_dir, self.option('method') + '_diff_bar.pdf')
            self.move_pdf(bar_pdf[0], new_path)

        stacked_pdf = glob.glob(os.path.join(self.work_dir, '*.diffexp_summary.stacked_bar.pdf'))
        if stacked_pdf:
            new_path = os.path.join(target_dir, self.option('method') + '_diff_summary.pdf')
            self.move_pdf(stacked_pdf[0], new_path)



    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_regexp_rules([
            [".", "", "差异分析结果目录"],
            [r".*_vs_.*\.xls", "xls", "差异分析结果表"],
            [r".*\.DE\.list", "xls", "差异基因列表"],
            [r".*summary\.xls", "xls", "差异统计结果表"],
            [r".*_vs_.*\.*annot\.xls", "xls", "差异统计注释结果表"],
            [r".*_vs_.*\..*annot\.xls", "xls", "差异统计注释结果表"],
            [r".*_vs_.*\..*normalize\.xls", "xls", "差异矩阵标准化结果表"],
            [r".*_vs_.*\..*sizeFactor\.xls", "xls", "差异矩阵标准化因子表"],
            [r".*_vs_.*\..*normFactor\.xls", "xls", "差异矩阵标准化因子表"],
            [r".*_vs_.*\.scatter\.pdf", "pdf", "表达量差异散点图"],
            [r".*_vs_.*\.volcano\.pdf", "pdf", "表达量差异火山图"],
            [r".*_diff_bar\.pdf", "pdf", "差异统计统计柱状图"],
            [r".*_diff_summary\.pdf", "pdf", "差异统计统计堆积图"],
        ])
        super(DiffexpBatchWorkflow, self).end()

    def run_tool(self):
        # filter count file
        if self.option("exp_level") == 'mRNA+sRNA':
            count_file = self.option('count').prop['path']
        elif self.option("exp_level") == 'mRNA':
            count_file = self.work_dir + '/known_seqs_count.matrix'
            with open(count_file, 'w') as fw, open(self.option('count').prop['path']) as fr:
                for line in fr:
                    if not line.startswith('sRNA') or line.startswith('seq_id'):
                        fw.write(line)
        else:
            count_file = self.work_dir + '/known_seqs_count.matrix'
            with open(count_file, 'w') as fw, open(self.option('count').prop['path']) as fr:
                for line in fr:
                    if line.startswith('sRNA') or line.startswith('seq_id'):
                        fw.write(line)
        opts = dict(
            exp=self.option('exp_matrix'),
            exp_type=self.option('exp_type'),
            count=count_file,
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
                    count=count_file,
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
                    batch_matrix=self.option('batch_matrix').path,
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
            else:
                opts = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=self.option('count').path,
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
        self.uniform = self.add_tool('prok_rna.batch.uniform')
        opts = {
            'method': self.option('method'),
            'input_dir': self.tool.output_dir,
        }
        if self.option('method').lower() in ['deseq2', 'degseq', 'edger', 'limma']:
            opts.update({'pvalue_padjust': self.option('pvalue_padjust')})
        self.uniform.set_options(opts)
        self.uniform.on('end', self.set_db)
        self.uniform.run()

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
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def paste_annotation(self):
        conn = self.all_exp.db['sg_annotation_stat']
        task_id = self.option('task_id')
        find_result = conn.find_one({"task_id": task_id, "type": "origin", 'status': 'end'})
        annot = os.path.join(find_result['result_dir'], 'summary/all_anno_detail.xls')
        if not os.path.exists(annot):
            annot = self.download_s3_file(annot, "annot.xls")
        all_annot = pd.read_table(annot, header=0, index_col=0)
        all_annot.rename(columns={'gene_id': 'seq_id'}, inplace=True)
        for each in glob.glob(os.path.join(self.tool.output_dir, '*_vs_*.*.xls')):
            if each.endswith('.annot.xls'):
                continue
            if each.endswith('.normalize.xls'):
                continue
            if each.endswith('.sizeFactor.xls'):
                continue
            if each.endswith('.normFactor.xls'):
                continue
            diff_pd = pd.read_table(each, header=0, sep='\t', index_col=0)
            diff_pd = diff_pd.join(all_annot, how='left')
            out_diff = each[:-3] + 'annot.xls'
            diff_pd.to_csv(out_diff, sep='\t', header=True, index=True)
