# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class SoapInfoAgent(Agent):
    """
    version v1.0
    author: zouxuan
    last modified:2017.0228
    """

    def __init__(self, parent):
        super(SoapInfoAgent, self).__init__(parent)
        options = [
            {"name": "map_dir", "type": "string", "default": ""},  # map结果
            {"name": "insertsize", "type": "infile", "format": "sample.insertsize_table"},  # 插入片段文件
            {"name": "soap_info", "type": "outfile", "format": "sequence.profile_table"},  # soap_info文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("map_dir"):
            raise OptionError("必须提供map结果", code="34101201")
        if not self.option("insertsize").is_set:
            raise OptionError("必须插入片段文件", code="34101202")

    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(SoapInfoAgent, self).end()


class SoapInfoTool(Tool):
    def __init__(self, config):
        super(SoapInfoTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/statistical/prepare_profile.pl"

    def mapinfo(self):
        cmd1 = '%s %s %s %s %s' % (self.perl_path, self.script1_path, self.option("map_dir"),
                                      self.option("insertsize").prop['path'],
                                      self.work_dir)
        self.logger.info(cmd1)
        self.logger.info("生成soap_info")
        command1 = self.add_command("saop_info", cmd1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("soap_info succeed")
        else:
            self.set_error("soap_info failed", code="34101201")
            raise Exception("soap_info failed")

    def set_output(self):
        oldfile = os.path.join(self.work_dir, 'soap.info')
        newfile = os.path.join(self.output_dir, 'soap_info.txt')
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile,newfile)
        self.option('soap_info',newfile)
        self.end()

    def run(self):
        super(SoapInfoTool, self).run()
        self.mapinfo()
        self.set_output()
