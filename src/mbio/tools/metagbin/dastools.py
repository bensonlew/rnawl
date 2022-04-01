#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil
from mbio.packages.metagbin.common_function import link_dir,bin_rename

class DastoolsAgent(Agent):
    """
    用于DAS_toools的输入文件生成
    version 1.0
    author: gaohao
    last_modify: 2019.01.08
    """
    def __init__(self, parent):
        super(DastoolsAgent, self).__init__(parent)
        options = [
            {"name": "sofware_bin", "type": "string", "default": "meatbat"},
            {"name": "contig_fa", "type": "infile", "format": "sequence.fasta"},  # 输入文件contigs.fa
            {"name": "maxbin_in", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "metabat_in", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "concoct_in", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "out", "type": "outfile", "format": "sequence.fasta_dir"},  # 生成bin的目录
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('contig_fa').is_set:
            raise OptionError("contig_fa不存在！")
        if self.option('sofware_bin') == "":
            raise OptionError("选择软件不存在！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 4
        size = float(os.path.getsize(self.option('contig_fa').path) / (1024 * 1024 * 1024))
        self._memory =str(int(size *15))+"G"

    def end(self):
        super(DastoolsAgent, self).end()


class DastoolsTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(DastoolsTool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/DAS_Tool-1.1.3:" + self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/usearch-v8.1:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/pullseq-1.0.2/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Prodigal-2.6.3:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/ruby/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/align/diamond-0.9.29"
        #self.ld_path = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64"
        self.set_environ(PATH=self.path)
        self.das_tool ="/bioinfo/metaGenomic/DAS_Tool-1.1.3/DAS_Tool"
        self.fa = self.option('contig_fa').path
        self.out = self.work_dir + '/bin'

    def run_dastools(self):
        sof_list = str(self.option('sofware_bin')).split(',')
        self.logger.info(sof_list)
        file = ''
        name = ''
        for sof in sof_list:
            if sof in ["metabat"]:
                file += self.option('metabat_in').prop['path'] + ','
                name += "metabat,"
            elif sof in ["maxbin"]:
                file += self.option('maxbin_in').prop['path'] + ','
                name += "maxbin,"
            elif sof in ["concoct"]:
                file += self.option('concoct_in').prop['path'] + ','
                name += "concoct,"
        cmd = "{} -i {} -l {} -c {} -o {} --write_bins 1 --score_threshold 0.1 -t 3".format(self.das_tool,file,name,self.fa,self.out)
        size = float(os.path.getsize(self.fa)/(1024*1024*1024))
        if size > 3:
            cmd += " --search_engine diamond"
        self.logger.info(cmd)
        self.logger.info("开始运行run_dastools")
        command = self.add_command("run_dastools", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            bin_rename(self.work_dir + '/bin_DASTool_bins',"all")
            self.logger.info("运行run_dastools完成")
        else:
            self.set_error("运行run_dastools运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + '/bin'):
            shutil.rmtree(self.output_dir + '/bin')
        link_dir(self.work_dir + '/bin_DASTool_bins',self.output_dir + '/bin')
        self.option('out',self.output_dir + '/bin')

    def run(self):
        """
        运行
        """
        super(DastoolsTool, self).run()
        self.run_dastools()
        self.set_output()
        self.end()