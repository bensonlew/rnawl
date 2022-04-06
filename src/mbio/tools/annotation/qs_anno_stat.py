# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class QsAnnoStatAgent(Agent):
    """
    qs群体感应统计
    """

    def __init__(self, parent):
        super(QsAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "qs_anno_result", "type": "infile", "format": "sequence.profile_table"},# 比对到QS库的注释结果文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("qs_anno_result").is_set:
            raise OptionError("Please enter an input file of QS!", code="31205601")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("Please enter an input file of abundance!", code="31205602")
        return True

    def set_resource(self):
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(QsAnnoStatAgent, self).end()

class QsAnnoStatTool(Tool):
    def __init__(self, config):
        super(QsAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR  + '/annotation/qs_anno_abundance.pl'
        self.python_path = "/miniconda2/bin/python"
        self.script_path = self.config.PACKAGE_DIR + "/metagenomic/scripts/anno_qs_graph.py"
        self.result_name = ''

    def run(self):
        """
        运行
        :return:
        """
        super(QsAnnoStatTool, self).run()
        self.run_qs_stat()
        self.run_qs_graph()
        self.set_output()
        self.end()


    def run_qs_stat(self):
        self.logger.info("start qs_stat")
        cmd = "{} {} -i {} -geneprofile {} -o {}".format(self.perl_path, self.script, self.option('qs_anno_result').prop['path'],
                                               self.option('reads_profile_table').prop['path'], 'qs')
        command = self.add_command('run_qs_stat', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qs_stat succeed")
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("qs_stat failed", code="31205601")

    def run_qs_graph(self):
        self.logger.info(self.option('group_table').prop['path'])
        cmd = "{} {} {} {}".format(self.python_path, self.script_path, self.work_dir + '/qs_class_profile.xls',
                                               self.option('group_table').prop['path'])
        command = self.add_command('run_qs_graph', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_qs_graph succeed")
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("run_qs_graph failed", code="31205602")

    def set_output(self):
        self.logger.info("start set_output")
        for i in ['qs_class_profile.xls','qs_lowest_profile.xls','anno_qs_graph.xls']:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            os.link(self.work_dir + '/' + i,self.output_dir + '/' + i)
