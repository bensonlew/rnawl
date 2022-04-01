# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from mbio.packages.bac_comp_genome.common_function import link_file
import subprocess


class TmhmmAnnoAgent(Agent):
    """
    跨膜蛋白注释结果 v1.0
    author: gaohao
    last_modify: 20190923
    """

    def __init__(self, parent):
        super(TmhmmAnnoAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample", "type": "string"} , # 样品名称
            {"name": "tmhmm", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置输入序列文件", code="31102401")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="31102402")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(TmhmmAnnoAgent, self).end()


class TmhmmAnnoTool(Tool):
    def __init__(self, config):
        super(TmhmmAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.tmhmm_path = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/tmhmm-2.0c/bin/tmhmm'
        self.tmhmm_sh = 'bioinfo//Genomic/Sofware/tmhmm-2.0c/bin/tmhmm.sh'
        self.perl_path = 'program/perl-5.24.0/bin/perl '
        self.script_path = self.config.PACKAGE_DIR + '/annotation/signalp_tmhmmtxt2xls.pl'
        self.python_script_path = self.config.PACKAGE_DIR + '/bacgenome/add_gene_info.py '
        self.python_path = '/program/Python/bin/python'

    def run(self):
        """
        运行
        :return:
        """
        super(TmhmmAnnoTool, self).run()
        self.tmhmm_predict()
        self.set_output()
        self.end()

    def tmhmm_predict(self):
        linkfile = self.work_dir + '/' + os.path.basename(self.option('query').prop['path'])
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(self.option('query').prop['path'], linkfile)
        cmd = '{} {} {} {}'.format(self.tmhmm_sh, self.tmhmm_path, linkfile, self.work_dir + '/tmhmm.txt')
        command = self.add_command('predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("tmhmm预测成功")
        else:
            self.set_error("tmhmm预测失败")
        cmd1 = '{} {} -faa {} -{} {} -outfile {}'.format(self.perl_path, self.script_path,
                                                         self.option('query').prop['path'],
                                                         'tmhmm', self.work_dir + '/tmhmm.txt',
                                                         self.work_dir + '/' + 'tmhmm_stat.xls')
        command1 = self.add_command('stat', cmd1).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("tmhmm注释结果统计成功")
        else:
            self.set_error("tmhmm注释结果统计失败")
        cmd2 = "awk -F \"\t\" \'{if((NR >13 && $7==\"YES\")||($1==\"Gene ID\")){print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$8}}\' " \
               + self.work_dir + "/tmhmm_stat.xls > " + self.work_dir + "/" + "gene_tmhmm_anno.xls"
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info("数据整理成功")
        except subprocess.CalledProcessError:
            self.set_error("数据整理失败")

    def set_output(self):
        link_file(self.work_dir + "/" + "gene_tmhmm_anno.xls", self.output_dir + "/" + self.option("sample") + ".tmhmm_anno.xls")
        self.option("tmhmm", self.output_dir + "/" + self.option("sample") + ".tmhmm_anno.xls")