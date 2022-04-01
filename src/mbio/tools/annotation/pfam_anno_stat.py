## -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import pandas as pd


class PfamAnnoStatAgent(Agent):
    """
    宏基因组pfam注释结果丰度统计，并生成丰度统计表
    """
    def __init__(self, parent):
        super(PfamAnnoStatAgent, self).__init__(parent)
        options = [
            {"name": "pfam_anno_table", "type": "infile", "format": "sequence.profile_table"},  #pfam注释结果文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},   #基因丰度文件
        ]
        self.add_option(options)
        
    def check_options(self):
        if not self.option("pfam_anno_table").is_set:
            raise OptionError("必须设置注释文件", code="31204001")
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置基因丰度文件", code="31204002")
        return True
    
    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"
        
    def end(self):
        super(PfamAnnoStatAgent, self).end()
        
class PfamAnnoStatTool(Tool):
    def __init__(self, config):
        super(PfamAnnoStatTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + "/annotation/mg_annotation/pfam_anno_abundance.pl"
    
    def run(self):
        """
        运行相关的函数
        :return:
        """
        super(PfamAnnoStatTool, self).run()
        self.run_pfam_stat()
        self.run_stat()
        self.set_output()
        self.end()

    def run_pfam_stat(self):
        self.logger.info('start pfam stat')
        anno_file = self.option('pfam_anno_table').prop["path"]
        new_anno_file = self.work_dir + "/gene_pfam_anno.xls"
        os.system("cp %s %s" % (anno_file, new_anno_file))
        os.system("sed -i '1d' %s" % (new_anno_file))
        cmd = "{} {} -q {} -p {} -o {}".format(self.perl_path, self.script, new_anno_file, self.option("reads_profile_table").prop["path"], self.output_dir)
        command = self.add_command('pfam_profile',cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('pfam_anno_stat succeed')
        elif command.return_code == -9 :
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error('pfam_stat failed', code="31204001")
            self.set_error('pfam_stat failed', code="31204002")

    def run_stat(self):
        """
        对统计结果相同的去重复
        :return:
        """
        stat_file_path = self.output_dir + '/gene_pfam_anno_stat.xls'
        new_stat_file_path = self.output_dir + '/new_pfam_anno_stat.xls'
        stat_file = pd.read_table(stat_file_path)
        stat_file = pd.DataFrame(stat_file)
        new_stat_file = stat_file.drop_duplicates().reset_index(drop=True)
        new_stat_file.to_csv(new_stat_file_path, sep='\t',index =False)
        os.system('rm %s'%(stat_file_path))
        os.system('mv %s %s'%(new_stat_file_path, stat_file_path))

        clan_path = self.output_dir + '/pfam_clan_profile.xls'
        new_clan_path = self.output_dir + '/new_pfam_clan_profile.xls'
        clan_file = pd.read_table(clan_path)
        clan_file = pd.DataFrame(clan_file)
        clan_file = clan_file[clan_file["#Clan ID"].astype('str') != '-']
        clan_file.to_csv(new_clan_path, sep='\t',index =False)
        os.system('rm %s'%(clan_path))
        os.system('mv %s %s'%(new_clan_path, clan_path))

    def set_output(self):
        self.logger.info("start set_output")
        if os.path.exists(self.output_dir + '/gene_pfam_anno_stat.xls'):
            pass
        self.logger.info('finish set_output')