# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os, sys
# import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class ScafAgpContigAgent(Agent):
    """
    细菌基因组的scaffold的序列得到最终scffold和contig序列
    version: v1
    author: hao.gao
    last_modify: 2018.02.06
    """

    def __init__(self, parent):
        super(ScafAgpContigAgent, self).__init__(parent)
        options = [
            {"name": "seq_scaf", "type": "infile", "format": "sequence.fasta"},  # 最佳Gapcloser的scaffold文件
            {"name": "sample_name", "type": "string"}, #样品名
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
            {"name": "contig", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
            {"name": "sort_scaf", "type": "bool", "default": True}
        ]
        self.step.add_steps('mp_cut', 'config_file', 'assemble', 'scaf_select')
        self.add_option(options)


    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('seq_scaf'):
            raise OptionError('必须提供最佳N50的scaffold序列文件！', code="31301501")
        if not self.option('sample_name'):
            raise OptionError('必须提供最佳样品名称！', code="31301502")

    def set_resource(self):
        """
        :return:
        """
        self._cpu = 4
        self._memory = "5G"

    def end(self):
        super(ScafAgpContigAgent, self).end()


class ScafAgpContigTool(Tool):
    def __init__(self, config):
        super(ScafAgpContigTool, self).__init__(config)
        self.seq_scaf = self.option('seq_scaf').prop['path']
        self.sample_name= self.option('sample_name')
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        self.perl_script_path = self.config.PACKAGE_DIR + '/bacgenome/'
        self.gapcloser = "/bioinfo/metaGenomic/SOAPdenovo2/v1.12-r6/bin/GapCloser"
        self.sort_scaf = self.work_dir + '/' + self.sample_name + '.sort.fna'
        self.agp_scaf = self.work_dir + '/' + self.sample_name + '.last.agp.fna'
        self.last_scaf = self.work_dir + '/' + 'result.scaf.fna'
        self.file1= self.work_dir + '/' + self.sample_name + '.scaf2contig'
        self.file2 = self.work_dir + '/' + self.sample_name + '.agp'
        self.sort = "Y" if self.option("sort_scaf") else "N"

    def run(self):
        """
        运行
        :return:
        """
        super(ScafAgpContigTool, self).run()
        self.run_rename1()
        self.run_apg()
        self.run_apg_contig()
        self.run_rename2()
        self.run_last_result()
        self.set_output()
        self.end()

    def run_rename1(self):
        cmd = '{} {}fasta-rename.pl {} scaffold N {}'.format(self.perl_path, self.perl_script_path,self.seq_scaf, self.sort_scaf)
        command = self.add_command('run_rename1', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_rename1运行完成")
        else:
            self.set_error("run_rename1运行出错!", code="31301501")

    def run_apg(self):
        cmd = '{} {}genome_scaffold_split.pl {} {}'.format(self.perl_path, self.perl_script_path,self.sort_scaf,
                                                      self.sample_name)
        command = self.add_command('run_apg', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_apg运行完成")
        else:
            self.set_error("run_apg运行出错!", code="31301502")

    def run_apg_contig(self):
        cmd = '{} {}genome_newFa.pl {} {} {} 200'.format(self.perl_path, self.perl_script_path,self.file1,
                                                  self.file2,self.agp_scaf)
        command = self.add_command('run_apg_contig', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_apg_contig运行完成")
        else:
            self.set_error("rrun_apg_contig运行出错!", code="31301503")

    def run_rename2(self):
        str = 'Scaffold'
        cmd = '{} {}fasta-rename.pl {} {} Y {}'.format(self.perl_path, self.perl_script_path, self.agp_scaf,str,
                                                             self.last_scaf)
        command = self.add_command('run_rename2', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_rename2运行完成")
        else:
            self.set_error("run_rename2运行出错!", code="31301504")

    def run_last_result(self):
        cmd = '{} {}genome_scaffold_split.pl {} {}'.format(self.perl_path, self.perl_script_path,self.last_scaf,
                                                      'result')
        command = self.add_command('run_last_result', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_last_result运行完成")
        else:
            self.set_error("run_last_result运行出错!", code="31301505")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('scaffold').set_path(self.last_scaf)
        self.option('contig').set_path(self.work_dir + '/result.scaf2contig')
        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '_scaf.fna'):
            os.remove(self.output_dir + '/' + self.option('sample_name') + '_scaf.fna')
        os.link(self.last_scaf,self.output_dir + '/' + self.option('sample_name') + '_scaf.fna')
        if os.path.exists( self.output_dir + '/' + self.option('sample_name') + '_ctg.fna'):
            os.remove( self.output_dir + '/' + self.option('sample_name') + '_ctg.fna')
        os.link(self.work_dir + '/result.scaf2contig', self.output_dir + '/' + self.option('sample_name') + '_ctg.fna')
        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '.agp'):
            os.remove(self.output_dir + '/' + self.option('sample_name') + '.agp')
        os.link(self.work_dir + '/result.agp', self.output_dir + '/' + self.option('sample_name') + '.agp')
        if os.path.exists(self.output_dir + '/all.seqid.xls' ):
            os.remove(self.output_dir + '/all.seqid.xls')
        os.link(self.work_dir + '/all.seqid.xls', self.output_dir + '/all.seqid.xls')
        self.logger.info("设置结果完成")