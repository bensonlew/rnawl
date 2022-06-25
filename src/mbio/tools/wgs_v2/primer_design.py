# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.08.06

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
            raise OptionError("请设置ref_fa文件")
        if not self.option("diff_variant"):
            raise OptionError("请设置比较分析的结果diff_variant文件")
        if self.option("tm1") >= self.option("tm2"):
            raise OptionError("tm1要小于tm2")
        if not self.option("product_size"):
            raise OptionError("请设置product_size")
        if not re.search(r"(\d+)-(\d+).*", self.option("product_size")):
            raise OptionError("%s product_size格式不正确,用-分隔",variables=(self.option("product_size")))
        if not self.option("primer_num"):
            raise OptionError("请设置primer_num", code="34504406")
        if self.option("primer_num") > 5 or self.option("primer_num") < 1:
            raise OptionError("primer_num:%s范围不为[1,5]", variables=(self.option("primer_num")))

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
        self.p3in_path = self.config.PACKAGE_DIR + "/wgs_v2/1.p3in.pl"
        self.p3out_path = self.config.PACKAGE_DIR + "/wgs_v2/11.p3out"
        self.primer3_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/primer3/src/"
        if self.option("project") == "wgs":
            self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome"  # 参考组配置文件
        else:
            self.json_path = self.config.SOFTWARE_DIR + "/database/dna_geneome"
        self.ref_fa = os.path.join(self.json_path, self.option("ref_fa"))
        # self.ref_fa = '/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/ref.fa'
        # self.merge_out = self.config.PACKAGE_DIR + "/wgs_v2/3.merge_p3out.pl"
        self.merge_out = self.config.PACKAGE_DIR + "/wgs_v2/3.merge_p3out.v2.pl"

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
            self.set_error("primer_in失败")

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
            self.set_error("运行primer3_core出错")
        os.system("cd {}".format(self.work_dir))

    def run_merge(self):
        """
        3.merge_p3out.pl
        :return:
        """
        cmd = "{} {} {} {}".format(self.perl_path, self.merge_out, self.work_dir + "/variation.p3out", self.option("primer_num"))
        command = self.add_command("primer_download", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("primer_download.py运行完成")
        else:
            self.set_error("primer_download.py运行失败")

    def write_result(self):
        write_lines = "#CHROM POS Total number Type Ref Alt Marker size(bp) Marker start(bp) Marker end(bp) FORWARD PRIMER1 (5'-3')  Tm(.C)  GC(%)  size  REVERSE PRIMER1 (5'-3')  Tm(.C)  GC(%)  size  PRODUCT1 size (bp)  start (bp)  end (bp)  FORWARD PRIMER2 (5'-3')  Tm(.C)  GC(%)  size  REVERSE PRIMER2 (5'-3')  Tm(.C)  GC(%)  size  PRODUCT2 size (bp)  start (bp)  end (bp)  FORWARD PRIMER3 (5'-3')  Tm(.C)  GC(%)  size  REVERSE PRIMER3 (5'-3')  Tm(.C)  GC(%)  size  PRODUCT3 size (bp)  start (bp)  end (bp)"
        with open(os.path.join(self.output_dir, "variation.result"), "w")as fw:
            fw.write(write_lines)

    def set_output(self):
        os.link(os.path.join(self.work_dir, "variation.result"), os.path.join(self.output_dir, "variation.result"))

    def run(self):
        super(PrimerDesignTool, self).run()
        with open(self.option("diff_variant"), "r")as fr:
            lines = fr.readlines()
            print len(lines)
            if len(lines) > 1:
                self.run_primer_in()
                with open((self.work_dir + "/variation.p3in"), "r")as fr:
                    lines = fr.readlines()
                    if len(lines) > 1:
                        self.run_primer3_core()
                        self.run_merge()
                        self.set_output()
                    else:
                        self.write_result()
            else:
                self.write_result()
        self.end()
