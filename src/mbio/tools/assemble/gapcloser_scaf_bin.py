# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os, sys
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class GapcloserScafBinAgent(Agent):
    """
    细菌基因组的scaffold的序列根据N50筛选
    version: v1
    author: hao.gao
    last_modify: 2018.01.29
    """

    def __init__(self, parent):
        super(GapcloserScafBinAgent, self).__init__(parent)
        options = [
            {"name": "seq_scaf", "type": "infile", "format": "sequence.fasta"},  # 最佳N50的scaffold文件
            {"name": "config", "type": "infile", "format": "bacgenome.config_file"},  # 配置的config文件
            {"name": "sample_name", "type": "string"},  #基因组名称
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"}  # 输出文件
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('seq_scaf'):
            raise OptionError('必须提供最佳N50的scaffold序列文件！', code="")
        if not self.option('config'):
            raise OptionError('必须提供config文件！', code="")
        if not self.option('sample_name'):
            raise OptionError('必须提供最佳样品名称！', code="")

    def set_resource(self):
        """
        :return:
        """
        self._cpu = 4
        self._memory = "5G"
        table_size = self.option("seq_scaf").get_size() /1024/1024   #guanqing.zou 20180905
        if table_size > 5:
            self._memory ="{}G".format(5+int((table_size -5)/5.0))

    def end(self):
        super(GapcloserScafBinAgent, self).end()


class GapcloserScafBinTool(Tool):
    def __init__(self, config):
        super(GapcloserScafBinTool, self).__init__(config)
        self.seq_scaf = self.option('seq_scaf').prop['path']
        #self.pe_list = self.option('PE_list').prop['path']
        self.sample_name = self.option('sample_name')
        self.perl_path = '/program/perl-5.24.0/bin/perl '
        self.perl_script_path = self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.gapcloser = "/bioinfo/metaGenomic/SOAPdenovo2/v1.12-r6/bin/GapCloser"
        self.remove_seq_scaf = self.work_dir + '/' + 'remove.300.fna'
        self.scaf = self.work_dir + '/' + 'last.gapcloser.fna'
        self.python_path = "/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/bacgenome/"
        #config_file = self.option('config')

    def run(self):
        """
        运行
        :return:
        """
        super(GapcloserScafBinTool, self).run()
        self.run_remove_scaf()
        #self.run_config()
        self.run_gapcloser()
        self.set_output()
        self.end()

    def run_remove_scaf(self):
        if os.path.exists(self.remove_seq_scaf):
            os.remove(self.remove_seq_scaf)
        cmd = '{} {}remove_scaf_fasta.pl {} {}'.format(self.perl_path, self.perl_script_path, self.seq_scaf, self.remove_seq_scaf)
        command = self.add_command('run_remove_scaf', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_remove_scaf运行完成")
        else:
            self.set_error("run_remove_scaf运行出错!", code="31300701")

    def run_gapcloser(self):
        cmd = '{} -a {} -b {} -o {} -t 4'.format(self.gapcloser, self.remove_seq_scaf, self.option('config').prop['path'], self.scaf)
        command = self.add_command('run_gapcloser', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_gapcloser运行完成" )
        else:
            self.set_error("run_gapcloser运行出错!", code="31300703")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + '/' + self.sample_name + '.scaffold.fna'):
            os.remove(self.output_dir + '/' + self.sample_name + '.scaffold.fna')
        self.logger.info("%s" % self.scaf)
        self.logger.info("%s" % (self.output_dir + '/' + self.sample_name + '.scaffold.fna'))
        os.link(self.scaf,self.output_dir + '/' + self.sample_name + '.scaffold.fna')
        self.option('seq').set_path(self.output_dir + '/' + self.sample_name + '.scaffold.fna')
        self.logger.info("设置gapcloser的补洞结果")
