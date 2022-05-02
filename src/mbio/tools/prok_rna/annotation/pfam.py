# -*- coding: utf-8 -*-
# __author__ = 'litangjian,zoujiaxun'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.gene_structure.pfam_domtblout import pfam_out
from mbio.packages.rna.annot_config import AnnotConfig

class PfamAgent(Agent):
    """
    transdecoder:orf,cds预测软件
    hmmscan：Pfam数据库比对工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11

    now
    hmmscan  Version 3.2.1
    last_modify: 20200110
    """

    def __init__(self, parent):
        super(PfamAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "prok_rna.common"},  # 输入文件
            {"name": "e_value", "type": "float", "default": 1e-3},
            {"name": "pfam_version", "type": "string", "default": "32"},
        ]
        self.add_option(options)
        self.step.add_steps('pfam')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.pfam.start()
        self.step.update()

    def step_end(self):
        self.step.pfam.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入fasta序列文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 20
        self._memory = '30G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_regexp_rules([
            ["./pfam_domain", "", "Pfam比对蛋白域结果信息"]
        ])

        super(PfamAgent, self).end()


class PfamTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(PfamTool, self).__init__(config)
        # self.hmmscan_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/"
        self.hmmscan_path = "miniconda2/bin/"
        # self.pfam_db = self.config.SOFTWARE_DIR + "/database/Pfam/Pfam-A.hmm"# hmm参考库
        self.pfam_db = AnnotConfig().get_file_path(
            file ="Pfam-A",
            db = "pfam",
            version = self.option("pfam_version"))

        self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]# 得到fasta文件名
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def hmmscan(self, pep):
        cmd = "{}hmmscan -E {} --domE {} --cpu 20 --noali --acc --notextw --domtblout {} --tblout {} {} {}".format(
            self.hmmscan_path, self.option("e_value"), self.option("e_value"), "pfam.domtblout", "pfam.tblout", self.pfam_db, pep)
        self.logger.info("开始运行hmmscan")
        self.logger.info(cmd)
        command = self.add_command("hmmscan", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行hmmscan结束")
        else:
            self.set_error("运行hmmscan出错")
        pfam_out("pfam.domtblout")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("将结果文件link到output文件夹下面")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))

        os.link(self.work_dir + "/pfam_domain", self.output_dir + "/pfam_domain")
        self.logger.info("成功将结果文件link到output文件夹下面")

    def run(self):
        super(PfamTool, self).run()
        self.hmmscan(self.option("fasta").prop['path'])
        self.set_output()
        self.end()
