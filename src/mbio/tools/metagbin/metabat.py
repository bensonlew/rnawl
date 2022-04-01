#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil
from mbio.packages.metagbin.common_function import link_dir,bin_rename

class MetabatAgent(Agent):
    """
    用于metabat2生成bin
    version 1.0
    author: gaohao
    last_modify: 2019.01.08
    """

    def __init__(self, parent):
        super(MetabatAgent, self).__init__(parent)
        options = [
            {"name": "minContig", "type": "string","default":"1500"},  #metabat2最小contigs
            {"name": "contig_fa", "type": "infile", "format": "sequence.fasta"},  # metabat2输入文件contigs.fa
            {"name": "depth_file", "type": "infile", "format": "sequence.profile_table"},  # metabat2的contigs的覆盖度文件
            {"name": "metabat_bin", "type": "outfile", "format": "sequence.fasta_dir"},#生成bin的目录
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
        num =float(self.size)/float(1000000000)
        if num <=0.5:
            num2 = 100
        else:
            num2 = int(round(num))*80
        self._memory =str(num2) + 'G'

    def end(self):
        super(MetabatAgent, self).end()


class MetabatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MetabatTool, self).__init__(config)
        self.metabat_path ="/bioinfo/metaGenomic/metabat/metabat2"
        self.fa =self.option('contig_fa').prop['path']
        self.depth = self.option('depth_file').prop['path']
        self.mincontig = self.option('minContig')
        if not os.path.exists(self.work_dir + '/metabat_bin'):
            os.mkdir(self.work_dir + '/metabat_bin')
        self.out =self.work_dir + '/metabat_bin/bin'

    def run_metabat(self):
        cmd = "{} -a {} -i {} -o {} -m {} -t 16 ".format(self.metabat_path, self.depth,self.fa,self.out,self.mincontig)
        self.logger.info(cmd)
        self.logger.info("开始运行run_metabat")
        command = self.add_command("run_metabat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            bin_rename(self.work_dir + "/metabat_bin", 'metabat')
            self.logger.info("运行run_metabat完成")
        else:
            self.set_error("运行run_metabat运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/metabat_bin'):
            shutil.rmtree(self.output_dir + '/metabat_bin')
        files = os.listdir(self.work_dir + '/metabat_bin')
        os.mkdir(self.output_dir + '/metabat_bin')
        for file in files:
            if os.path.getsize(self.work_dir + '/metabat_bin/' + file) < 500000000:
                os.link(self.work_dir + '/metabat_bin/' + file,self.output_dir + '/metabat_bin/' + file)
        self.option('metabat_bin',self.output_dir + '/metabat_bin')

    def run(self):
        """
        运行
        """
        super(MetabatTool, self).run()
        self.run_metabat()
        self.set_output()
        self.end()