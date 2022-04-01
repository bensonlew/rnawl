# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/23'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_file


class FalconAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(FalconAgent, self).__init__(parent)
        options = [
            {"name": "fq_file", "type": "infile", "format": "sequence.fastq", "required": True},
            {"name": "genome_size", "type": "string", "default": "5000000"},
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "scaf", "type": "outfile", "format": "sequence.fasta", "required": True},
            {"name": "cpu", "type": "int", "default": 8},
            {"name": "mem", "type": "int", "default": 50}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = self.option("cpu")
        self._memory = "%sG" % self.option("mem")


class FalconTool(Tool):
    def __init__(self, config):
        super(FalconTool, self).__init__(config)
        self.python_path = '/program/Python/bin/python'
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.cfg = self.config.PACKAGE_DIR + "/bacgenome/run_falcon_cfg.pl"
        self.falcon = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/falcon/bin/fc_run"
        # falcon环境变量
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin")
        # fq2fa环境变量
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.fq2fa = '/bioinfo/metaGenomic/MaxBin-2.2.5/auxiliary/idba-1.1.1/bin/fq2fa'

    def run_fq2fa(self):
        cmd = self.fq2fa + " --paired %s %s " % (self.option("fq_file").prop["path"], self.work_dir + "/input.fa")
        command = self.add_command("run_fq2fa", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fq2fa运行完成")
        else:
            self.set_error("fq2fa运行出错！")

    def run_cfg(self):
        """
        description
        :return:
        """
        cmd = "{} {} --infile {} --gsize {} --cpu  {} --mem {}".format(self.perl_path, self.cfg,  self.work_dir + "/input.fa", self.option("genome_size"), self.option("cpu"), self.option("mem"))
        command = self.add_command("run_cfg", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_cfg运行完成")
        else:
            self.set_error("run_cfg运行出错！")

    def run_falcon(self):
        cmd = "{} {} falcon.cfg".format(self.python_path, self.falcon)
        command = self.add_command("run_falcon", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_falcon运行完成")
        else:
            self.end()
            # self.set_error("run_falcon运行出错！")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        from_f = os.path.join(self.work_dir, "2-asm-falcon/a_ctg_all.fa")
        to_f = os.path.join(self.output_dir, self.option("sample_name") + ".scaf.fna")
        link_file(from_f, to_f)
        self.option("scaf").set_path(to_f)

    def run(self):
        super(FalconTool, self).run()
        self.run_fq2fa()
        self.run_cfg()
        self.run_falcon()
        self.set_output()
        self.end()