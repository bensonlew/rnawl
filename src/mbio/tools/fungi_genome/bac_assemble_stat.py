# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.3.26
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class BacAssembleStatAgent(Agent):
    """
    真菌基因组扫描图组装结果评估
    """
    def __init__(self, parent):
        super(BacAssembleStatAgent, self).__init__(parent)
        options = [
            {"name": "scaf_seq", "type": "infile","format": "sequence.fasta"},  # 参考序列文件夹
            {"name": "cont_seq", "type": "infile","format": "sequence.fasta"},  # 输入列表文件
            {'name': 'sample_name', "type": "string"},  # 样本名
            {'name': 'seq_type', "type": "string","default":"pacbio"},  #
            ]
        self.add_option(options)
        self.window =[]

    def check_options(self):
        if not self.option("scaf_seq").is_set:
            raise OptionError("必须添加scaf_seq的序列文件！", code="32100101")
        if not self.option("cont_seq").is_set:
            raise OptionError("必须添加cont_seq的序列的文件！", code="32100102")
        if not self.option('sample_name'):
            raise OptionError('必须输入sample_name样品名称！', code="32100103")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(BacAssembleStatAgent, self).end()

class BacAssembleStatTool(Tool):
    def __init__(self, config):
        super(BacAssembleStatTool, self).__init__(config)
        self.scaf_fa = self.option("scaf_seq").prop['path']
        self.cont_fa = self.option("cont_seq").prop['path']
        self.sample_name = self.option("sample_name")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.perl_script2 = self.config.PACKAGE_DIR + "/fungi_genome/"
        self.set_environ(PERL5LIB = self.config.SOFTWARE_DIR + '/program/perl-5.24.0/lib/site_perl/5.24.0')  #add for nb zouguanqing 20181022

    def run_asse_stat(self):
        output = self.work_dir + '/' + self.sample_name + '.summary.xls'
        cmd ='{} {}assemble_summary.pl {} 1000 {}'.format(self.perl_path,self.perl_script,self.scaf_fa,output)
        self.logger.info(cmd)
        command = self.add_command("assemble_summary", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("assemble_summary运行完成")
        else:
            self.set_error("assemble_summary运行出错!", code="32100101")


    def run_asse_seq(self):
        cmd = '{} {}bac_genome_stat.pl {} {}'.format(self.perl_path,self.perl_script,self.scaf_fa,self.cont_fa)
        self.logger.info(cmd)
        command = self.add_command("asse_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("asse_seq运行完成")
        else:
            self.set_error("asse_seq运行出错!", code="32100102")

    def run_asse_graph(self):
            cmd = '{} {}seq_distribution.pl {} {} {}'.format(self.perl_path, self.perl_script2, self.scaf_fa,
                                                                   self.cont_fa,self.sample_name)
            self.logger.info(cmd)
            command = self.add_command('run_asse_graph', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_asse_graph 运行完成")
            else:
                self.set_error("run_asse_graph 运行出错!", code="32100103")

    def set_output(self):
        if not os.path.exists(self.output_dir + '/summary'):
            os.mkdir(self.output_dir + '/summary')
        os.link(self.work_dir + '/scaffolds.stat.xls',self.output_dir + '/summary/' + self.sample_name +'_assembly_scaffold_details.xls')
        os.link(self.work_dir + '/contigs.stat.xls', self.output_dir + '/summary/'  + self.sample_name +'_assembly_contig_details.xls')
        os.link(self.work_dir + '/' + self.sample_name  + '.summary.xls', self.output_dir + '/summary/' + self.sample_name  + '_assembly_summary.xls')
        if not os.path.exists(self.output_dir + '/len'):
            os.mkdir(self.output_dir + '/len')
        path1 = self.output_dir + '/len/' + self.sample_name + '.scaffolds.len.xls'
        path2 = self.output_dir + '/len/' + self.sample_name + '.contigs.len.xls'
        path3 = self.work_dir + '/' + self.sample_name + '.scaffolds.len.xls'
        path4 = self.work_dir + '/' + self.sample_name + '.contigs.len.xls'
        if os.path.exists(path1):
            os.remove(path1)
        else:
            os.link(path3, path1)
        if os.path.exists(path2):
            os.remove(path2)
        else:
            os.link(path4, path2)
        self.end()

    def run(self):
        super(BacAssembleStatTool, self).run()
        self.run_asse_stat()
        self.run_asse_seq()
        self.run_asse_graph()
        self.set_output()
