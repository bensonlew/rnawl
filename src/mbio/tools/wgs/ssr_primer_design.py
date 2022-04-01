# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import subprocess


class SsrPrimerDesignAgent(Agent):
    """
    软件：primer3_core，脚本：sample.ssr.p3in.pl、sample.ssr.p3out.pl
    对样本参考基因组的SSR进行引物设计
    """
    def __init__(self, parent):
        super(SsrPrimerDesignAgent, self).__init__(parent)
        options = [
            {"name": "misa_file", "type": "string"},  # scafSeq.misa文件
            {"name": "scafseq_file", "type": "string"},  # scafSeq文件
            {"name": "tm1", "type": "float", "default": 57.0},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 63.0},  # float, Tm2 (℃),要大于tm1
            {"name": "product_size", "type": "string", "default": "300-500"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int", "default": 3},  # Max Pairs Primer Number,范围:[1,5]
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("misa_file"):
            raise OptionError("请设置misa_file文件", code="34506401")

    def set_resource(self):
        self._cpu = 2
        self._memory = "30G"

    def end(self):
        super(SsrPrimerDesignAgent, self).end()


class SsrPrimerDesignTool(Tool):
    def __init__(self, config):
        super(SsrPrimerDesignTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64")
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.ssr_p3in_path = self.config.PACKAGE_DIR + "/wgs/sample.ssr.p3in.pl"
        self.ssr_p3out_path = self.config.PACKAGE_DIR + "/wgs/sample.ssr.p3out.pl"
        self.primer3_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/primer3/src/"
        self.python = "program/Python/bin/python"
        self.primer_download = self.config.PACKAGE_DIR + "/wgs/primer_download.py"

    def run_ssr_primer_in(self):
        """
        sample.ssr.p3in.pl
        """
        cmd = "{} {} -d {} ".format(self.perl_path, self.ssr_p3in_path, self.option("misa_file"))
        cmd += "-r {} -T1 {} -T2 {}".format(self.option("product_size"), self.option("tm1"), self.option("tm2"))
        cmd += " -p {} -ref {}".format(self.option("primer_num"), self.option("scafseq_file"))
        cmd += " -o {}".format(self.work_dir)
        command = self.add_command("ssr_primer_in", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr_primer_in完成")
        else:
            self.set_error("ssr_primer_in失败", code="34506401")

    def run_primer3_core(self):
        """
        ./primer3_core
        """
        self.sample_name = os.path.basename(self.option("misa_file")).split(".")[0]
        p3in_file = os.path.join(self.work_dir, self.sample_name + ".p3in")
        cmd = "cd {} && ./primer3_core --output {} {}".format(self.primer3_path, self.work_dir + "/variation.p3out",
                                                              p3in_file)
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("运行primer3_core完成")
        except subprocess.CalledProcessError:
            self.logger.info("运行primer3_core出错")
            self.set_error("运行primer3_core出错", code="34506402")
        os.system("cd {}".format(self.work_dir))

    def run_ssr_primer_out(self):
        """
        sample.ssr.p3out.pl
        """
        cmd = "{} {} {}".format(self.perl_path, self.ssr_p3out_path, os.path.join(self.work_dir, "variation.p3out"))
        command = self.add_command("ssr_primer_out", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr_primer_out完成")
        else:
            self.set_error("ssr_primer_out失败", code="34506403")

    def run_primer_download(self):
        """
        primer_download.py
        """
        cmd = "{} {} -i {} -o {} -type {}".format(self.python, self.primer_download, self.work_dir + "/variation.result",
                                                  os.path.join(self.output_dir, self.sample_name + "_result.xls"), "ssr")
        command = self.add_command("primer_download", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("primer_download.py运行完成")
        else:
            self.set_error("primer_download.py运行失败", code="34506406")

    def set_output(self):
        os.link(os.path.join(self.work_dir, "variation.result"), os.path.join(self.output_dir,  self.sample_name + ".result"))

    def run(self):
        super(SsrPrimerDesignTool, self).run()
        self.run_ssr_primer_in()
        self.run_primer3_core()
        self.run_ssr_primer_out()
        self.run_primer_download()
        self.set_output()
        self.end()
