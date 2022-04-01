#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil,re
from mbio.packages.metagbin.common_function import link_dir,bin_rename

class MaxbinAgent(Agent):
    """
    用于Maxbin2生成bin
    version 1.0
    author: gaohao
    last_modify: 2019.01.08
    """

    def __init__(self, parent):
        super(MaxbinAgent, self).__init__(parent)
        options = [
            {"name": "minContig", "type": "string","default":"1000"},  #metabat2最小contigs
            {"name": "contig_fa", "type": "infile", "format": "sequence.fasta"},  # metabat2输入文件contigs.fa
            {"name": "depth_file", "type": "infile", "format": "sequence.profile_table"},  # metabat2的contigs的覆盖度文件
            {"name": "maxbin_bin", "type": "outfile", "format": "sequence.fasta_dir"},#生成bin的目录
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('contig_fa').is_set:
            raise OptionError("组装的fasta文件不存在！")
        if not self.option('depth_file').is_set:
            raise OptionError("contigs的覆盖度文件不存在！")

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
            num2 = int(num) * 120 + 200
        self._memory = str(num2) + 'G'

    def end(self):
        super(MaxbinAgent, self).end()


class MaxbinTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MaxbinTool, self).__init__(config)
        self.path =self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/MaxBin-2.2.5:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/MaxBin-2.2.5/auxiliary/idba-1.1.1/bin"
        self.LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64"
        self.set_environ(PATH=self.path,LD_LIBRARY_PATH=self.LD_LIBRARY_PATH)
        self.perl = "/program/perl-5.24.0/bin/perl"
        self.maxbin_path =self.config.SOFTWARE_DIR +"/bioinfo/metaGenomic/MaxBin-2.2.5/run_MaxBin.pl"
        self.fa =self.option('contig_fa').prop['path']
        self.depth = self.option('depth_file').prop['path']
        self.mincontig = self.option('minContig')
        if not os.path.exists(self.work_dir + '/maxbin_bin'):
            os.mkdir(self.work_dir + '/maxbin_bin')
        self.out =self.work_dir + '/maxbin_bin/bin'

    def run_maxbin(self):
        cmd = "{} {} -abund {} -contig {} -out {} -min_contig_length {} -thread 16 -prob_threshold 0.8 -preserve_intermediate".format(self.perl,self.maxbin_path, self.depth,self.fa,self.out,self.mincontig)
        self.logger.info(cmd)
        self.logger.info("开始运行run_maxbin")
        command = self.add_command("run_maxbin", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            bin_rename(self.work_dir + '/maxbin_bin', 'maxbin')
            self.logger.info("运行run_maxbin完成")
        else:
            self.set_error("运行run_maxbin运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/maxbin_bin'):
            shutil.rmtree(self.output_dir + '/maxbin_bin')
        os.mkdir(self.output_dir + '/maxbin_bin')
        files =os.listdir(self.work_dir + '/maxbin_bin')
        for file in files:
            if re.search(r'.fa$',file) and os.path.getsize(self.work_dir + '/maxbin_bin/' + file) < 500000000:
                os.link(self.work_dir + '/maxbin_bin/' + file, self.output_dir + '/maxbin_bin/' + file)
        self.option('maxbin_bin',self.output_dir + '/maxbin_bin')

    def run(self):
        """
        运行
        """
        super(MaxbinTool, self).run()
        self.run_maxbin()
        self.set_output()
        self.end()