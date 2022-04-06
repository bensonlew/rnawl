# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from multiprocessing import Process


class FormatTaxonAgent(Agent):
    """
    对传入的ref_taxon 格式进行规范，处理重名情况
    """
    def __init__(self, parent):
        super(FormatTaxonAgent, self).__init__(parent)
        options = [
            {"name": "in_taxon_table", "type": "infile", "format": "taxon.seq_taxon"},  # 输入的taxon文件
            {"name": "out_taxon_table", "type": "outfile", "format": "taxon.seq_taxon"}  # 输出的taxon文件

        ]
        self.add_option(options)



    def check_options(self):
        """
        参数检测
        """
        if not self.option("in_taxon_table").is_set:
            raise OptionError("输入的taxon文件不能为空")


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["taxa.new.xls", "xls", "输出taxon文件"]

        ])
        super(FormatTaxonAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 1
        self._memory = "1G"



class FormatTaxonTool(Tool):
    def __init__(self, config):
        super(FormatTaxonTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.taxon_clean =os.path.join(self.config.PACKAGE_DIR, "meta/scripts/taxon_clean.py")
    def run_taxon(self):
        os.system("dos2unix {}".format(self.option('in_taxon_table').prop['path'])) ## add_by qingchen.zhang @20200813
        self.out =  self.output_dir + "/taxon.new.xls"
        #cmd = "{} {} {} {}".format(self.python_path, self.taxon_clean, self.option('in_taxon_table').prop['path'], self.option('out_taxon_table').prop['path'])
        cmd = "{} {} {} {}".format(self.python_path, self.taxon_clean, self.option('in_taxon_table').prop['path'], self.out)
        self.logger.info(cmd)
        command = self.add_command("format_taxon", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.option('out_taxon_table',self.out)
        else:
            self.set_error("command失败")

    def run(self):
        super(FormatTaxonTool, self).run()
        self.run_taxon()
        self.end()



