# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'20181011

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.annotation.mg_annotation.pfam_select import pfam_out
from mbio.packages.annotation.mg_annotation.pfam_select import anno_select
import unittest

class PfamAnnoAgent(Agent):
    """
    宏基因组对比对结果文件进行注释
    """
    def __init__(self, parent):
        super(PfamAnnoAgent, self).__init__(parent)
        options = [
            {"name": "hmmscan_result", "type": "infile", "format": "paternity_test.tab"},  # 比对结果table文件
            {"name": "pfam_anno_result", "type": "outfile", "format": "sequence.profile_table"},    #注释详情结果表
            {"name": "result_format", "type": "string", "default": "meta"}, #输出结果类型
            {"name": "best", "type": "bool", "default": True},  # 是否一条序列只挑选一个结果
            {"name": "version", "type": "string", "default": "pfam_v33.1"}  # 是否一条序列只挑选一个结果
        ]
        self.add_option(options)
        self.result_name = ""

    def check_options(self):
        if not self.option("hmmscan_result").is_set:
            raise OptionError("必须提供hmmscan比对结果作为输入文件", code="31205301")
        return True

    def set_resource(self):
        self.cpu = 2
        self._memory = "5G"

    def end(self):
        self.option('pfam_anno_result', os.path.join(self.work_dir, "gene_pfam_anno.xls"))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.',"", "结果输出目录"],
            ['gene_pfam_anno.xls', 'xls', '序列详细注释文件']
        ])
        super(PfamAnnoAgent, self).end()


class PfamAnnoTool(Tool):
    def __init__(self, config):
        super(PfamAnnoTool, self).__init__(config)
        self.python_path =  self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.script_path = self.config.PACKAGE_DIR + "/annotation/scripts/meta_pfam_mongo.py"
        if self.option("version") in ["pfam_v33.1"]:
            self.pfam_table_path = self.config.SOFTWARE_DIR + "/database/pfam_33.1/Pfam-A.database.xls"
        else:
            self.pfam_table_path = self.config.SOFTWARE_DIR + "/database/pfam_31/Pfam-A.database.xls"
        self.result_name = ''

    def run(self):
        """
        运行注释函数
        :return:
        """
        super(PfamAnnoTool, self).run()
        self.get_align_result()
        self.run_pfam_anno()
        self.set_output()
        self.end()

    def get_align_result(self):
        self.result_name = self.option('hmmscan_result').path
        best = self.option('best')
        pfam_out(self.result_name, best)
    """
    def run_pfam_anno(self):
        self.out_name = os.path.join(self.work_dir, "pfam_domain")
        if os.path.exists(self.out_name):
            cmd = '{} {} -i {} -o {}'.format(self.python_path, self.script_path, self.out_name, self.output_dir + "/gene_pfam_anno.xls")
            self.logger.info(cmd)
        else:
            self.set_error('hmmscan结果筛选出错',code="")
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行pfam_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行pfam_anno出错', code="")
    """

    def run_pfam_anno(self):
        table = os.path.join(self.work_dir, "pfam_domain")
        anno_profile_ref = self.pfam_table_path
        anno_select(table, anno_profile_ref)

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + 'gene_pfam_anno.xls'):
            os.remove(self.output_dir + '/' + 'gene_pfam_anno.xls')
        os.link(self.work_dir + '/' + 'gene_pfam_anno.xls', self.output_dir + '/' + 'gene_pfam_anno.xls')

