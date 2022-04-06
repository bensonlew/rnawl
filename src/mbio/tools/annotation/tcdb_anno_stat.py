# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import os

class TcdbAnnoStatAgent(Agent):
    """
    宏基因ardb注释结果丰度统计
    """

    def __init__(self, parent):
        super(TcdbAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "tcdb_anno_table", "type": "infile", "format": "sequence.profile_table"},
            # 基因注释具体信息结果表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": 'is_modify',"type":"string",'default':'F'}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("tcdb_anno_table").is_set:
            raise OptionError("必须设置注释文件", code="31203901")
        if not self.option('reads_profile_table').is_set:
            raise OptionError("必须设置基因丰度文件", code="31203902")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(TcdbAnnoStatAgent, self).end()


class TcdbAnnoStatTool(Tool):
    def __init__(self, config):
        super(TcdbAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        self.python = '/miniconda2/bin/python'
        self.script = self.config.PACKAGE_DIR + '/annotation/abund_stat_z.py'


    def run(self):
        """
        运行
        :return:
        """
        super(TcdbAnnoStatTool, self).run()
        if self.option('is_modify') == 'T':
            self.change_table()
        self.run_tcdb_stat()
        self.set_output()
        self.end()
    def change_table(self):
        anno_file = pd.read_table(self.option('tcdb_anno_table').prop['path'],sep='\t')
        split_id = anno_file['TCDB ID'].str.split('.',expand=True)
        anno_file['class_id'] = split_id[0]
        anno_file['subclass_id'] = split_id[0]+'.'+split_id[1]
        anno_file['family_id'] = split_id[0]+'.'+split_id[1] + '.' + split_id[2]
        anno_file['#Query'] = anno_file['Gene ID'].str[:-2]
        del anno_file['Gene ID']
        anno_file.set_index('#Query',inplace=True)
        anno_file.to_csv(self.work_dir + '/tcdb_anno_split_id.xls',sep='\t')


    def run_tcdb_stat(self):
        self.logger.info("start tcdb_stat")
        if self.option('is_modify') == 'T':
            cmd = "{} {} {} {} tcdb".format(self.python, self.script, self.option('reads_profile_table').prop['path'], self.work_dir + '/tcdb_anno_split_id.xls')
        else:
            cmd = "{} {} {} {} tcdb".format(self.python, self.script, self.option('reads_profile_table').prop['path'],self.option('tcdb_anno_table').prop['path'])
        command = self.add_command('tcdb_profile', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("tcdb_stat succeed")
        elif command.return_code in [-9,1]:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("tcdb_stat failed", code="31203901")
            self.set_error("tcdb_stat failed", code="31203902")

    def set_output(self):
        self.logger.info("start set_output")
        work_files = ['class_abund_out.xls','subclass_abund_out.xls','family_abund_out.xls','tcdb_abund_out.xls',
                      'Class_gene_stat.xls','SubClass_gene_stat.xls','Family_gene_stat.xls','TCDB_gene_stat.xls']
        for i in work_files:
            if os.path.exists(self.output_dir + '/' + i):
                os.remove(self.output_dir + '/' + i)
            os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)

