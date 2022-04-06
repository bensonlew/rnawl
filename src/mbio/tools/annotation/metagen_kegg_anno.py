# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import os
import subprocess


class MetagenKeggAnnoAgent(Agent):
    """
    调用宏基因kegg注释需要的脚本 meta_kegg_mongo.py
    author: zhouxuan
    last_modify: 2017.0607
    """

    def __init__(self, parent):
        super(MetagenKeggAnnoAgent, self).__init__(parent)
        options = [
            {"name": "kegg_xml", "type": "infile", "format": "align.blast.blast_xml"},  # 输入的比对结果xml文件
            {"name": "result_dir", "type": "outfile", "format": "meta_genomic.kegg_dir"}  # 输出结果文件夹
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("kegg_xml").is_set:
            raise OptionError("必须设置输入文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
            ])
        super(MetagenKeggAnnoAgent, self).end()


class MetagenKeggAnnoTool(Tool):
    def __init__(self, config):
        super(MetagenKeggAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/meta_kegg_mongo.py'

    def run(self):
        """
        运行
        :return:
        """
        super(MetagenKeggAnnoTool, self).run()
        self.run_kegg_anno()
        self.set_output()
        self.end()

    def run_kegg_anno(self):
        table = xml2table(self.option('kegg_xml').prop['path'], self.work_dir + '/tmp_kegg_table.xls')
        cmd = '{} {} {} {}'.format(self.python_path, self.python_script, table, self.output_dir)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行kegg_anno完成')
        except subprocess.CalledProcessError:
            self.set_error('运行kegg_anno出错')

    def set_output(self):
        self.logger.info("set_output")
        try:
            self.option("result_dir", self.output_dir)
        except Exception as e:
            raise Exception("SET_OUTFILE FAILED {}".format(e))
        self.logger.info("start set_output")