# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.gene_structure.pfam_domtblout import pfam_out


class LongorfAgent(Agent):
    """
    transdecoder:Longorf,cds预测软件
    hmmscan：Pfam数据库比对工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11
    """

    def __init__(self, parent):
        super(LongorfAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},  # 输入文件
            {"name": "p_length", "type": "int", "default": 50},  # 最小蛋白长度
        ]
        self.add_option(options)
        self.step.add_steps('Longorf')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Longorf.start()
        self.step.update()

    def step_end(self):
        self.step.Longorf.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入fasta序列文件", code = "32005001")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 20
        self._memory = '30G'

    def end(self):
        super(LongorfAgent, self).end()


class LongorfTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(LongorfTool, self).__init__(config)
        self.transdecoder_path = "bioinfo/gene-structure/TransDecoder-3.0.0/"
        self.hmmscan_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/"
        self.pfam_db = self.config.SOFTWARE_DIR + "/database/Pfam/Pfam-A.hmm"# hmm参考库
        self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]# 得到fasta文件名
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # TransDecoder-3.0.0的运行依赖perl的环境，所以程序必须设置一个环境变量，直接在我们的服务器程序里测试单条是不行的，我们没有perl_env
        self.perl_lib = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin'
        self.set_environ(PATH=self.perl_lib)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    # 此处要进行构建训练集，这里的fa文件不能拆分
    def td_longorfs(self):
        self.logger.info(self.option("p_length"))
        cmd = "{}TransDecoder.LongOrfs -t {} -m {}".format(self.transdecoder_path, self.option("fasta").prop["path"],
        self.option("p_length"))
        print(cmd)
        self.logger.info("开始提取长orf")
        self.logger.info(cmd)
        command = self.add_command("transdecoder_longorfs", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("提取结束！")
        else:
            self.set_error("提取orf过程出错", code = "32005002")

    def run(self):
        super(LongorfTool, self).run()
        self.td_longorfs()
        self.end()


        