# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'


from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import os
import re

class CypsAnnoStatAgent(Agent):
    """
    宏基因组cyps注释结果丰度统计表
    """
    def __init__(self, parent):
        super(CypsAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "cyps_anno_table", "type": "infile", "format": "sequence.profile_table"},  #输入文件,基因注释的结果表
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"}   #丰度表
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("cyps_anno_table").is_set:
            raise OptionError('必须设置注释文件', code="31205201")
        if not self.option("reads_profile_table").is_set:
            raise OptionError('必须设置基因丰度文件', code="31205202")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'


class CypsAnnoStatTool(Tool):
    def __init__(self, config):
        super(CypsAnnoStatTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/cyps_anno_abundance.pl'
        #self.script = '/mnt/ilustre/users/sanger-dev/sg-users/zhangqingchen/test/cyps_test/package/cyps_anno_abundance.pl'

    def run(self):
        """
        运行tool
        :return:
        """
        super(CypsAnnoStatTool, self).run()
        self.run_cyps_stat()
        self.run_stat()
        self.set_output()
        self.end()

    def run_cyps_stat(self):
        self.logger.info("start cyps stat")
        cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.script, self.option('cyps_anno_table').prop['path'], self.option('reads_profile_table').prop['path'], self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('cyps_profile', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('cyps_stat succed')
        elif command.return_code == -9:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("cyps_stat failed", code="31205201")
            self.set_error('cyps_stat failed', code="31205202")

    def run_stat(self):
        """
        对统计结果相同的去重复
        :return:
        """
        stat_file_path = self.output_dir + '/cyps_anno_stat.xls'
        new_stat_file_path = self.output_dir + '/new_cyps_anno_stat.xls'
        stat_file = pd.read_table(stat_file_path)
        stat_file = pd.DataFrame(stat_file)
        new_stat_file = stat_file.drop_duplicates().reset_index(drop=True)
        new_stat_file.to_csv(new_stat_file_path, sep='\t',index =False)
        os.system('rm %s'%(stat_file_path))
        os.system('mv %s %s'%(new_stat_file_path, stat_file_path))

    def set_output(self):
        self.logger.info('start set_output')
        if os.path.exists(self.output_dir + '/cyps_anno_stat.xls'):
            pass
        self.logger.info('finish set_output')