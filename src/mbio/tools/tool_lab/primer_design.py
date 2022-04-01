# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200429

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import math
import subprocess


class PrimerDesignAgent(Agent):
    """
    引物设计
    """
    def __init__(self, parent):
        super(PrimerDesignAgent, self).__init__(parent)
        options = [
            {"name": "seq", "type": "string"},
            {"name": "primer_length", "type": "int", "default": 20},
            {"name": "tm1", "type": "float", "default": 55},  # float, Tm1 (℃)
            {"name": "tm2", "type": "float", "default": 60},  # float, Tm2 (℃),要大于tm1
            {"name": "gc1", "type": "float", "default": 40},  # float, gc1 (℃)
            {"name": "gc2", "type": "float", "default": 60},  # float, gc2 (℃),要大于gc1
            {"name": "product_size", "type": "string", "default": "100-300"},  # Product Size(bp),数值要为int,范围间用-分隔，多个用,分隔,如：100-300
            {"name": "primer_num", "type": "int", "default": 3},  # Max Pairs Primer Number,范围:[1,5]
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq"):
            raise OptionError("请设置seq")
        if not self.option("primer_length"):
            raise OptionError("请输入引物长度")
        if self.option("tm1") >= self.option("tm2"):
            raise OptionError("tm1要小于tm2")
        if self.option("gc1") >= self.option("gc2"):
            raise OptionError("gc1要小于gc2")
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
        self.set_environ(PKG_CONFIG_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64/pkgconfig")
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.python = "program/Python/bin/python"
        self.p3in_path = self.config.PACKAGE_DIR + "/tool_lab/primer_file.py"
        self.primer3_path = self.config.SOFTWARE_DIR + "/bioinfo/WGS/primer3/src/"
        self.merge_out = self.config.PACKAGE_DIR + "/tool_lab/3.merge_p3out.v2.pl"

    def run_primer_in(self):
        """
        primer_file.py
        """
        cmd = "{} {} -s {} -l {} -mt {} -nt {} -mg {} -ng {} -p {} -n {} -o {}".format(self.python, self.p3in_path, self.option("seq")
                                                                 , self.option("primer_length"),
                                                    self.option("tm2"), self.option("tm1"), self.option("gc2")
                                                                 ,self.option("gc1"),
                                                                       self.option("product_size"),
                                                                             self.option("primer_num"), self.work_dir)

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
        # cmd = "cd {} && ./primer3_core --output {} {}".format(self.primer3_path, self.work_dir + "/variation.p3out",
        #                                                       self.work_dir + "/prime_file.txt")
        cmd = " {} -default_version=1 --output {} {}".format(self.primer3_path + "/primer3_core", self.work_dir + "/variation.p3out",
                                                              self.work_dir + "/prime_file.txt")
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

    def run_file(self):
        path1 = os.path.join(self.work_dir, "variation.result")
        path2 = os.path.join(self.work_dir, "primer_results.txt")
        with open(path1) as f, open(path2, "w") as w:
            lines = f.readlines()
            if len(lines) < 2 or lines[1].strip() == "":
                self.logger.info("运行结果为空，请检查序列长度是否太短！")
                self.set_error("运行结果为空，请检查序列长度是否太短！")
            else:
                for line in lines[1:-1]:
                    item = line.strip().split("\t")
                    w.write("OLIGO\tstart\tlen\ttm\tgc%\tany\t3'\tseq\n")
                    w.write("LEFT PRIMER\t" + item[19] + "\t" + item[3] + "\t" + item[1] + "\t" + item[2] + "\t" + item[
                        5] + "\t" + item[4] + "\t" + item[0].upper() + "\n")
                    w.write("RIGHT PRIMER\t" + item[20] + "\t" + item[12] + "\t" + item[10] + "\t" + item[11] + "\t" + item[
                        14] + "\t" + item[13] + "\t" + item[9].upper() + "\n")
                    w.write(
                        "PRODUCT SIZE: " + item[18] + "," + "PAIR ANY COMPL: " + item[22] + "," + "PAIR 3' COMPL: " + item[
                            23] + "\n")
                    sequence = list(item[24].upper())
                    row = int(math.ceil(float(len(sequence)) / 60))
                    sequence_new = ""
                    sequence_bank = ["-" for j in range(len(sequence))]
                    for i in range(int(item[19]), int(item[3]) + int(item[19])):
                        sequence_bank[i] = ">"
                    for i in range(int(item[20]) - int(item[12]) + 1, int(item[20]) + 1):
                        sequence_bank[i] = "<"
                    for i in range(row):
                        sequence_new += ("".join(sequence[60 * i: (60 * (i + 1))]) + "\n")
                        sequence_new += ("".join(sequence_bank[60 * i: (60 * (i + 1))]) + "\n")
                    w.write(sequence_new)
                    w.write("\n\n")
                w.write("注：\n")
                w.write("Forward Primer (5'-3')：正向引物序列；\nTm：正向引物退火温度；\nGC%：正向引物序列GC含量；\nLength：正向引物序列长度；\n")
                w.write(
                    "Reverse Primer (5'-3')：反向引物序列；\nTm：反向引物退火温度；\nGC%：反向引物序列GC含量；\nLength：反向引物序列长度；\nProduct Size：引物扩增产物大小；")

    def set_output(self):
        os.link(os.path.join(self.work_dir, "primer_results.txt"), os.path.join(self.output_dir, "primer_results.txt"))

    def run(self):
        super(PrimerDesignTool, self).run()
        self.run_primer_in()
        self.run_primer3_core()
        self.run_merge()
        self.run_file()
        self.set_output()
        self.end()
