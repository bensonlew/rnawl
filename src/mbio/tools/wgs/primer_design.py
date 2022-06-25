# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.23

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import subprocess


class PrimerDesignAgent(Agent):
    """
    引物设计
    """
    def __init__(self, parent):
        super(PrimerDesignAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "string"},
            {"name": "diff_variant", "type": "string"},  # 比较分析的结果diff.variant
            {"name": "tm1", "type": "float", "default": 57.0},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 63.0},  # float, Tm2 (℃),要大于tm1
            {"name": "product_size", "type": "string"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int"},  # Max Pairs Primer Number,范围:[1,5]
            {"name": "project", "type": "string", "default": "wgs"}
            # {"name": "update_info", "type": "string"},
            # {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa文件", code="34504401")
        if not self.option("diff_variant"):
            raise OptionError("请设置比较分析的结果diff_variant文件", code="34504402")
        if self.option("tm1") >= self.option("tm2"):
            raise OptionError("tm1要小于tm2", code="34504403")
        if not self.option("product_size"):
            raise OptionError("请设置product_size", code="34504404")
        if not re.search(r"(\d+)-(\d+).*", self.option("product_size")):
            raise OptionError("%s product_size格式不正确,用-分隔",variables=(self.option("product_size")), code="34504405")
        if not self.option("primer_num"):
            raise OptionError("请设置primer_num", code="34504406")
        if self.option("primer_num") > 5 or self.option("primer_num") < 1:
            raise OptionError("primer_num:%s范围不为[1,5]", variables=(self.option("primer_num")), code="34504407")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(PrimerDesignAgent, self).end()


class PrimerDesignTool(Tool):
    def __init__(self, config):
        super(PrimerDesignTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64")
        self.perl_path = "miniconda2/bin/perl"
        self.python = "miniconda2/bin/python"
        self.p3in_path = self.config.PACKAGE_DIR + "/wgs/1.p3in.pl"
        self.p3out_path = self.config.PACKAGE_DIR + "/wgs/1.p3out.pl"
        self.primer3_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/primer3/src/"
        if self.option("project") == "wgs":
            self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome"  # 参考组配置文件
        else:
            self.json_path = self.config.SOFTWARE_DIR + "/database/dna_geneome"
        self.ref_fa = os.path.join(self.json_path, self.option("ref_fa"))
        self.primer_download = self.config.PACKAGE_DIR + "/wgs/primer_download.py"

    def run_primer_in(self):
        """
        1.p3in.pl
        """
        cmd = "{} {} -d {} ".format(self.perl_path, self.p3in_path, self.option("diff_variant"))
        cmd += "-r {} -T1 {} -T2 {}".format(self.option("product_size"), self.option("tm1"), self.option("tm2"))
        cmd += " -p {} -ref {}".format(self.option("primer_num"), self.ref_fa)
        cmd += " -o {}".format(self.work_dir)
        command = self.add_command("primer_in", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("primer_in完成")
        else:
            self.set_error("primer_in失败", code="34504401")

    def run_primer3_core(self):
        """
        ./primer3_core
        """
        cmd = "cd {} && ./primer3_core --output {} {}".format(self.primer3_path, self.work_dir + "/variation.p3out",
                                                              self.work_dir + "/variation.p3in")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("运行primer3_core完成")
        except subprocess.CalledProcessError:
            self.logger.info("运行primer3_core出错")
            self.set_error("运行primer3_core出错", code="34504402")
        os.system("cd {}".format(self.work_dir))

    def run_primer_out(self):
        """
        1.p3out.pl
        """
        cmd = "{} {} {}".format(self.perl_path, self.p3out_path, self.work_dir + "/variation.p3out")
        command = self.add_command("primer_out", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("primer_out完成")
        else:
            self.set_error("primer_out失败", code="34504403")

    def run_primer_download(self):
        """
        primer_download.py
        """
        cmd = "{} {} -i {} -o {} -type {}".format(self.python, self.primer_download, self.work_dir + "/variation.result",
                                                  self.output_dir + "/variant_result.xls", "primer")
        command = self.add_command("primer_download", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("primer_download.py运行完成")
        else:
            self.set_error("primer_download.py运行失败", code="34504404")

    def set_output(self):
        os.link(os.path.join(self.work_dir, "variation.result"), os.path.join(self.output_dir, "variation.result"))

    def set_db(self):
        """
        将结果导到mongo数据库
        """
        self.logger.info("保存结果到mongo")
        primer_id = self.option("main_id")
        primer_api = self.api.api("wgs.primer_design")
        primer_result = os.path.join(self.work_dir, "variation.result")
        primer_api.add_sg_primer_detail(primer_id, primer_result)

    def run(self):
        super(PrimerDesignTool, self).run()
        self.run_primer_in()
        self.run_primer3_core()
        self.run_primer_out()
        self.run_primer_download()
        # self.set_output()
        # if self.option("main_id"):
        #     self.set_db()
        self.end()
