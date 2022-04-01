# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
from biocluster.core.exceptions import OptionError
import subprocess


class BacSummaryAgent(Agent):
    """
    last_modify: 2018.04.23
    """

    def __init__(self, parent):
        super(BacSummaryAgent, self).__init__(parent)
        options = [
            {"name": "gene_statistics", "type": "infile", "format": "sequence.profile_table"},
            {"name": "assemble", "type": "infile", "format": "sequence.profile_table"},
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "cog", "type": "infile", "format": "sequence.profile_table"},
            {"name": "kegg", "type": "infile", "format": "sequence.profile_table"},
            {"name": "analysis", "type": "string", "default": "uncomplete"},  ###流程分析模式complete，uncomplete
            {'name': 'sample_name', "type": "string"},  # 样本名
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("gene_statistics").is_set:
            raise OptionError("请传入基因的统计文件！", code="31403301")
        if not self.option("assemble").is_set:
            raise OptionError("请传入组装评估统计文件！", code="31403302")
        if not self.option("rrna_gff").is_set:
            raise OptionError("请传入rRNA的gff3文件", code="31403303")
        if not self.option("trna_gff").is_set:
            raise OptionError("请传入tRNA的gff3文件", code="31403304")
        if not self.option("cog").is_set:
            raise OptionError("请传入COG注释表", code="31403305")
        if not self.option("kegg").is_set:
            raise OptionError("请传入KEGG注释表", code="31403306")
        if not self.option("analysis"):
            raise OptionError("请传入分析类型！", code="31403307")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(BacSummaryAgent, self).end()

class BacSummaryTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(BacSummaryTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + '/bacgenome/'

    def run_summary(self):
        cmd = '{} {}Summary_project.pl {} {} {} {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.option('analysis'),
                                                                       self.option('sample_name'),self.option('assemble').prop['path'],self.option('gene_statistics').prop['path'],self.option('rrna_gff').prop['path'],self.option('trna_gff').prop['path'],self.option('kegg').prop['path'],self.option('cog').prop['path'])
        self.logger.info(cmd)
        command = self.add_command("run_summary", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_summary运行完成")
        else:
            self.set_error("run_summary运行出错!", code="31403301")



    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '.project.summary'):
            os.remove(self.output_dir + '/' + self.option('sample_name') + '.project.summary')
        os.link(self.work_dir + '/' + self.option('sample_name') + '.project.summary',self.output_dir + '/' + self.option('sample_name') + '.project.summary')
        self.logger.info("设置结果目录成功")

    def run(self):
        super(BacSummaryTool, self).run()
        self.run_summary()
        self.set_output()
        self.end()
