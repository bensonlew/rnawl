#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil,re
from mbio.packages.metagbin.common_function import link_dir,bin_rename

class ConcoctAgent(Agent):
    """
    用于Concoct生成bin
    version 1.0
    author: gaohao
    last_modify: 2019.01.08
    """

    def __init__(self, parent):
        super(ConcoctAgent, self).__init__(parent)
        options = [
            {"name": "minContig", "type": "string","default":"1000"},  #metabat2最小contigs
            {"name": "contig_fa", "type": "infile", "format": "sequence.fasta"},  # metabat2输入文件contigs.fa
            {"name": "bam_dir", "type": "infile", "format": "metagbin.bam_dir"},
            {"name": "concoct_bin", "type": "outfile", "format": "sequence.fasta_dir"},#生成bin的目录
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('contig_fa').is_set:
            raise OptionError("组装的fasta文件不存在！")
        if not self.option('bam_dir').is_set:
            raise OptionError("bam_dir文件夹不存在！")

    def set_resource(self):
        """
        所需资源
        """
        self.size = os.path.getsize(self.option("contig_fa").prop["path"])
        self._cpu = 16
        num = float(self.size) / float(1000000000)
        if num <= 1:
            num2 = 200
        else:
            num2 = int(num) * 100 + 200
        self._memory = str(num2) + 'G'

    def end(self):
        super(ConcoctAgent, self).end()


class ConcoctTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(ConcoctTool, self).__init__(config)
        self.path=self.config.SOFTWARE_DIR + "/program/Python/bin"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/library/gsl23/lib"
        self.set_environ(PATH=self.path,LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.python = "/miniconda2/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + "/program/Python/bin/"
        self.cut_sh = "../../../../../.." + self.config.PACKAGE_DIR + '/metagbin/concoct_cut.sh'
        self.cov_sh ="../../../../../.." + self.config.PACKAGE_DIR + '/metagbin/concoct_coverage.sh'
        self.merge_sh = "../../../../../.." + self.config.PACKAGE_DIR + '/metagbin/concoct_merge.sh'
        self.fa =self.option('contig_fa').prop['path']
        self.bam = self.option('bam_dir').prop['path'] + "/*.sorted.bam"
        self.mincontig = self.option('minContig')
        self.bed =self.work_dir + '/contigs_10K.bed'
        self.k10_fa = self.work_dir + '/contigs_10K.fa'
        self.cov = self.work_dir + '/coverage_table.tsv'
        if not os.path.exists(self.work_dir + '/concoct'):
            os.mkdir(self.work_dir + '/concoct')
        self.bin_concoct = self.work_dir + '/concoct'
        if not os.path.exists(self.work_dir + '/concoct_bin'):
            os.mkdir(self.work_dir + '/concoct_bin')
        self.concoct = self.work_dir + '/concoct_bin'

    def run_cut_seq(self):
        cmd = "{} {} {} {} {} ".format(self.cut_sh,self.python_script + 'cut_up_fasta.py',self.fa, self.bed,self.k10_fa)
        self.logger.info(cmd)
        self.logger.info("开始运行run_cut_seq")
        command = self.add_command("run_cut_seq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_cut_seq完成")
        else:
            self.set_error("运行run_cut_seq运行出错!")

    def run_seq_coverage(self):
        cmd = "{} {} {} {} {} ".format(self.cov_sh,self.python_script + 'concoct_coverage_table.py',self.bed, self.bam, self.cov)
        self.logger.info(cmd)
        self.logger.info("开始运行run_seq_coverage")
        command = self.add_command("run_seq_coverage", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_seq_coverage完成")
        else:
            self.set_error("运行run_seq_coverage运行出错!")

    def run_concoct(self):
        cmd = "{} {}concoct --composition_file {} --coverage_file {} -t 16 -l {} -r 100 -b {}".format(self.python,self.python_script,self.k10_fa,self.cov,self.mincontig,self.bin_concoct)
        self.logger.info(cmd)
        self.logger.info("开始运行run_concoct")
        command = self.add_command("run_concoct", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_concoct完成")
        else:
            self.set_error("运行run_concoct运行出错!")

    def run_merge_seq(self):
        path_file = ''
        for i in os.listdir(self.work_dir + '/concoct/'):
            if re.search('clustering_gt',i):
                path_file = self.work_dir + '/concoct/' + i
        cmd = "{} {} {} {} ".format(self.merge_sh, self.python_script + 'merge_cutup_clustering.py', path_file, self.work_dir + '/concoct/clustering_merged.csv')
        self.logger.info(cmd)
        self.logger.info("开始运行run_merge_seq")
        command = self.add_command("run_merge_seq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_merge_seq完成")
        else:
            self.set_error("运行run_merge_seq运行出错!")

    def run_tiqu_seq(self):
        cmd = "{} {}extract_fasta_bins.py {} {} --output_path {} ".format(self.python,self.python_script, self.fa, self.work_dir + '/concoct/clustering_merged.csv',self.concoct)
        self.logger.info(cmd)
        self.logger.info("开始运行run_tiqu_seq")
        command = self.add_command("run_tiqu_seq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            bin_rename(self.concoct,'concoct')
            self.logger.info("运行run_tiqu_seq完成")
        else:
            self.set_error("运行run_tiqu_seq运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/concoct_bin'):
            shutil.rmtree(self.output_dir + '/concoct_bin')
        files = os.listdir(self.work_dir + '/concoct_bin')
        os.mkdir(self.output_dir + '/concoct_bin')
        for file in files:
            if os.path.getsize(self.work_dir + '/concoct_bin/' + file) < 500000000:
                os.link(self.work_dir + '/concoct_bin/' + file, self.output_dir + '/concoct_bin/' + file)
        self.option('concoct_bin',self.output_dir + '/concoct_bin')

    def run(self):
        """
        运行
        """
        super(ConcoctTool, self).run()
        self.run_cut_seq()
        self.run_seq_coverage()
        self.run_concoct()
        self.run_merge_seq()
        self.run_tiqu_seq()
        self.set_output()
        self.end()