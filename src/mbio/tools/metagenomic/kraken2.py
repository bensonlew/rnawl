# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
from mbio.packages.toolapps.common import link_dir


class Kraken2Agent(Agent):
    """
    Kraken2基于reads进行丰度统计和功能计算
    """

    def __init__(self, parent):
        super(Kraken2Agent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "sample", "type": "string"},
            {"name": "readlen", "type": "int", "default": 100},
            {"name": "confidence", "type": "float", "default": 0.5},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("read1").is_set:
            raise OptionError("请输入read1文件！")
        if not self.option("read2").is_set:
            raise OptionError("请输入read2文件！")
        if not self.option("sample"):
            raise OptionError("请输入sample的样品名称！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 8
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
        ])
        super(Kraken2Agent, self).end()

class Kraken2Tool(Tool):
    def __init__(self, config):
        super(Kraken2Tool, self).__init__(config)
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin:" + self.config.SOFTWARE_DIR + "/gcc/7.2.0/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/nhmmer:"
        self.lib = self.config.SOFTWARE_DIR + "/gcc/7.2.0/lib64:"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.kraken2 = "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin/kraken2"
        self.out = self.output_dir + "/" + self.option("sample") + ".taxon.xls"
        self.db = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/kraken2-2.0.8-beta/database"
        self.brakren = "/bioinfo/metaGenomic/Bracken-master/"
        self.brakren_src = self.config.PACKAGE_DIR + "/metagenomic/est_abundance.py"
        self.kraken = "/bioinfo/metaGenomic/kraken2-2.0.8-beta/bin/"
        self.python = "program/Python/bin/python"

    def run_kraken2(self):
        """
        description
        :return:
        """
        cmd = "{} --threads 8 --db {} --report {} --confidence {} --output {}".format(self.kraken2, self.db, self.out, self.option("confidence"), self.work_dir + "/" + self.option("sample") + ".mapping.stdout")
        cmd +=" --paired {} {}".format(self.option("read1").prop['path'], self.option("read2").prop['path'])
        command = self.add_command("run_kraken2", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            os.system("grep '^C' {} >{}".format(self.work_dir + "/" + self.option("sample") + ".mapping.stdout",
                                                self.work_dir + "/" + self.option("sample") + ".mapping.xls"))
            self.logger.info("run_kraken2运行完成！")
        else:
            self.set_error("run_kraken2运行完成运行出错!")

    def run_brakren_build(self):
        """
        构建brakren的索引
        """
        cmd = "{}bracken-build -d {} -t 8 -x {} -l {} ".format(self.brakren, self.db, self.kraken, self.option("readlen"))
        command = self.add_command("run_brakren_build", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_brakren_build运行完成！")
        else:
            self.set_error("run_brakren_build运行完成运行出错!")

    def run_brakren(self):
        """
        运行brakren的，得到物种各水平分类
        """
        k = "{}/database{}mers.kmer_distrib".format(self.db, self.option("readlen"))
        cmd_f = "{} {} -i {} -k {} -o {}/{}_{}.xls -l {} -t 0"
        for i in ["D", "K", "P", "C", "O", "F", "G", "S"]:
            cmd = cmd_f.format(self.python, self.brakren_src, self.out, k,
                               self.output_dir, self.option("sample"), i, i)
            psd = "run_brakren_" + i.lower()
            command = self.add_command(psd, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("{}运行完成!".format(psd))
            else:
                self.set_error("{}运行完成运行出错!".format(psd))

    def run(self):
        super(Kraken2Tool, self).run()
        self.run_kraken2()
        self.run_brakren_build()
        self.run_brakren()
        self.end()
