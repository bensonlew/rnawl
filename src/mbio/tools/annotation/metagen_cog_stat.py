# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class MetagenCogStatAgent(Agent):
    """
    宏基因cog注释结果统计.py v1.0
    author: zhouxuan
    last_modify: 2017.0605
    """

    def __init__(self, parent):
        super(MetagenCogStatAgent, self).__init__(parent)
        options = [
            {"name": "cog_table", "type": "infile", "format": "annotation.cog.cog_table"},
            # 比对到string库的注释结果文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("cog_table").is_set:
            raise OptionError("必须设置输入文件")
        if not self.option('reads_profile_table').is_set:
            raise OptionError
        return True

    def set_resource(self):
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
            ])
        super(MetagenCogStatAgent, self).end()


class MetagenCogStatTool(Tool):
    def __init__(self, config):
        super(MetagenCogStatTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "miniconda2/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/meta_genomic_cog_stat.py'

    def run(self):
        """
        运行
        :return:
        """
        super(MetagenCogStatTool, self).run()
        self.run_cog_stat()
        self.set_output()
        self.end()

    def run_cog_stat(self):
        self.logger.info("start cog_stat")
        cmd = "{} {} -i {} -r {} -o {}".format(self.python_path, self.python_script, self.option('cog_table').prop['path'],
                                               self.option('reads_profile_table').prop['path'], self.output_dir)
        command = self.add_command('tax_profile', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cog_stat succeed")
        else:
            self.set_error("cog_stat failed")
            raise Exception("cog_stat failed")

    def set_output(self):
        self.logger.info("start set_output")