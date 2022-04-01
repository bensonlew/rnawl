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


class DiffexpWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="infile", format="prok_rna.common",),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            # dict(name="type", type="string"),
            dict(name="task_id", type="string"),
            dict(name="exp_level", type="string"),
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("prok_rna.diffexp")
        self.all_exp = self.api.api("prok_rna.all_exp")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(DiffexpWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # add result info
        self.all_exp.add_diffexp(self.tool.output_dir,
                                 main_id=self.option('diff_main_id'),
                                 diff_method=self.option('method'),
                                 create_geneset=False,
                                 pvalue_padjust=self.option('pvalue_padjust')
                            )
        self.paste_annotation()
        self.end()

    def end(self):
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
        ])
        super(DiffexpWorkflow, self).end()

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
        options = dict(
            exp=self.option('exp_matrix'),
            count=count_file,
            group=self.option('group'),
            method=self.option('method'),
            cmp=self.option('cmp'),
            pvalue=self.option('pvalue'),
            pvalue_padjust=self.option('pvalue_padjust'),
            padjust_way=self.option('padjust_way'),
            fc=self.option('fc'),
        )
        self.tool.set_options(options)
        self.tool.run()

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
