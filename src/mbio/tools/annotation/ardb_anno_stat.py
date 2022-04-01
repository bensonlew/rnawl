# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re

class ArdbAnnoStatAgent(Agent):
    """
    宏基因ardb注释结果丰度统计表
    author: shaohua.yuan
    last_modify:
    """

    def __init__(self, parent):
        super(ArdbAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "ardb_anno_table", "type": "infile", "format": "sequence.profile_table"},
            # 基因注释具体信息结果表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ardb_anno_table").is_set:
            raise OptionError("必须设置注释文件", code="31200301")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("必须设置基因丰度文件", code="31200302")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'  # 内存2G增加到5G by GHD @20180309 改回 by GHD @ 20180428
        # tmp_mem = 5 * (self._rerun_time + 1)  # 每次因拼接失败而重运行的内存增加5G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('vfdb_anno_stat use memory : ' + self._memory)
"""
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
            ])
        super(ArdbAnnoStatAgent, self).end()
"""

class ArdbAnnoStatTool(Tool):
    def __init__(self, config):
        super(ArdbAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        #self.script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/ardb_anno_abudance.pl'
        # self.script = '/bioinfo/annotation/scripts/ardb_anno_abudance.pl'
        self.perl =  'program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/ardb_anno_abudance.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(ArdbAnnoStatTool, self).run()
        self.run_ardb_stat()
        self.set_output()
        self.end()

    def run_ardb_stat(self):
        self.logger.info("start ardb_stat")
        #cmd = "perl {} -q {} -p {} -o {}".format(self.script, self.option('ardb_anno_table').prop['path'],\
        #                                       self.option('reads_profile_table').prop['path'], self.output_dir)
        cmd = "{} {} -q {} -p {} -o {}".format(self.perl, self.script, self.option('ardb_anno_table').prop['path'],self.option('reads_profile_table').prop['path'], self.output_dir)
        command = self.add_command('ardb_profile', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("ardb_stat succeed")
        elif command.return_code == -9:                        # change return_code by guhaidong @ 20180604
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.set_error("ardb_stat failed", code="31200301")
            raise Exception("ardb_stat failed")

    def set_output(self):
        self.logger.info("start set_output")
