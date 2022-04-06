# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import pandas as pd

class MvirdbAnnoStatAgent(Agent):
    """
    宏基因ardb注释结果丰度统计
    """

    def __init__(self, parent):
        super(MvirdbAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "mvirdb_anno_table", "type": "infile", "format": "sequence.profile_table"},
            # 基因注释具体信息结果表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},
            {"name":"is_modify",'type':'string','default':'F'}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("mvirdb_anno_table").is_set:
            raise OptionError("必须设置注释文件", code="31205401")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("必须设置基因丰度文件", code="31205402")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(MvirdbAnnoStatAgent, self).end()


class MvirdbAnnoStatTool(Tool):
    def __init__(self, config):
        super(MvirdbAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        self.python = '/miniconda2/bin/python'
        self.script = self.config.PACKAGE_DIR + '/annotation/abund_stat_z.py'

    def run(self):
        """
        运行
        :return:
        """
        super(MvirdbAnnoStatTool, self).run()
        if self.option('is_modify') == 'T':
            self.change_table()
        self.run_mvirdb_stat()
        self.set_output()
        self.end()

    def change_table(self):
        anno_file = pd.read_table(self.option('mvirdb_anno_table').prop['path'],sep='\t')
        anno_file['#Query'] = anno_file['Gene ID'].str[:-2]
        del anno_file['Gene ID']
        anno_file.set_index('#Query',inplace=True)
        anno_file.to_csv(self.work_dir + '/mvirdb_anno_rm_1.xls',sep='\t')


    def run_mvirdb_stat(self):
        self.logger.info("start mvirdb_stat")
        if self.option('is_modify') == 'T':
            cmd = "{} {} {} {} mvirdb".format(self.python,self.script, self.option('reads_profile_table').prop['path'], self.work_dir + '/mvirdb_anno_rm_1.xls')
        else:
            cmd = "{} {} {} {} mvirdb".format(self.python,self.script, self.option('reads_profile_table').prop['path'], self.option('mvirdb_anno_table').prop['path'])
        command = self.add_command('mvirdb_profile', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mvirdb_stat succeed")
        elif command.return_code == -9:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("mvirdb_stat failed", code="31205401")
            self.bind_object.set_error("mvirdb_stat failed", code="31203902")

    def set_output(self):
        self.logger.info("start set_output")
        work_files = ['Factor_abund_out.xls','Type_abund_out.xls','Factor_gene_stat.xls','Type_gene_stat.xls']
        for i in work_files:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)

