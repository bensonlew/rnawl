# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
from mbio.packages.denovo_rna.gene_structure.snp_position import ssr_position


class SsrAgent(Agent):
    """
    misa:SSR分析软件
    primer3：引物设计软件
    version 1.0
    author: qindanhua
    last_modify: 2016.08.03
    """

    def __init__(self, parent):
        super(SsrAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "primer", "type": "bool", "default": True},  # 是否设计SSR引物
            {"name": "bed", "type": "infile", "format": "denovo_rna.gene_structure.bed"}  # bed格式文件
        ]
        self.add_option(options)
        self.step.add_steps('ssr')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.ssr.start()
        self.step.update()

    def step_end(self):
        self.step.ssr.finish()
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
        self._cpu = 10
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./misa_stat.xls", "xls", "ssr类型统计表"]
        ])
        result_dir.add_regexp_rules([
            [r"misa$", "misa", "ssr结果"]
        ])
        # print self.get_upload_files()
        super(SsrAgent, self).end()


class SsrTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(SsrTool, self).__init__(config)
        self.misa_path = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/gene-structure/misa/")
        self.primer3_path = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/gene-structure/primer3-2.3.7/")
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/scripts/"
        self.fasta_name = os.path.basename(self.option("fasta").prop["path"])
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/"
        self.python_path = "program/Python/bin/"

    def misa(self):
        # self.logger.info(self.work_dir + "/" + self.fasta_name)
        fasta_path = self.option("fasta").prop["path"]
        shutil.copy(fasta_path, self.work_dir)
        fasta_copy = self.work_dir + "/" + self.fasta_name
        self.logger.info(fasta_copy + ".misa")
        cmd = "{}perl {}misa.pl {}".format(self.perl_path, self.misa_path, fasta_copy)
        # print(cmd)
        self.logger.info("开始运行misa")
        command = self.add_command("misa", cmd)
        command.run()
        self.wait()
        if command.return_code == 0 or None:
            self.logger.info("运行misa结束！")
            if self.option("bed").is_set:
                self.logger.info("统计ssr位置信息")
                ssr_position(fasta_copy + ".misa", self.option("bed").prop["path"])
                self.logger.info("统计ssr位置信息完成")
        else:
            self.set_error("运行misa过程出错")

    def primer(self):
        cmd = "{}python {}primer3_core.py -p {} -i {} -o {}".format(self.python_path, self.script_path, self.primer3_path, self.fasta_name + ".p3in", self.fasta_name + ".misa.p3out")
        # print(cmd)
        self.logger.info(cmd)
        self.logger.info("开始运行primer")
        command = self.add_command("primer", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行primer结束")
        else:
            self.set_error("运行primer出错")

    def primer_in(self):
        p3_in_cmd = "{}perl {}p3_in.pl {}".format(self.perl_path, self.misa_path, self.fasta_name + ".misa")
        self.logger.info("转换misa结果为primer输入格式")
        print(p3_in_cmd)
        self.logger.info(p3_in_cmd)
        command = self.add_command("p3_in", p3_in_cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行p3_in结束！")
        else:
            self.set_error("运行p3_in过程出错")

    def primer_out(self):
        p3_out_cmd = "{}perl {}p3_out.pl {} {}".format(self.perl_path, self.misa_path, self.fasta_name + ".misa.p3out", self.fasta_name + ".misa")
        self.logger.info("转换misa结果为primer输出格式")
        self.logger.info(p3_out_cmd)
        command = self.add_command("p3_out", p3_out_cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行p3_out结束！")
        else:
            self.set_error("运行p3_out过程出错")

    def ssr_position(self):
        self.logger.info(self.config.SOFTWARE_DIR + self.misa_path)
        ssr_position_cmd = "{}perl {}misa_anno.pl -i {} -orf {}".format(self.config.SOFTWARE_DIR + "/" + self.perl_path, self.misa_path, self.fasta_name + ".misa", self.option("bed").prop["path"])
        self.logger.info("判断ssr位置")
        self.logger.info(ssr_position_cmd)
        try:
            subprocess.check_output(ssr_position_cmd, shell=True)
            self.logger.info("OK")
            return True
        except subprocess.CalledProcessError:
            self.logger.info("过程出错")
            return False

    def set_output(self):
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        misa_stat = self.work_dir+'/' + self.fasta_name + ".statistics"
        if os.path.exists(misa_stat):
            start_n = 0
            with open(misa_stat, "r") as s, open(self.output_dir + "/misa_stat.xls", "w") as w:
                for n, line in enumerate(s):
                    if re.match(r"Frequency of classified repeat types", line):
                        start_n = n + 2
                    if start_n == 0:
                        continue
                    else:
                        if n > start_n:
                            w.write(line)
        os.link(self.work_dir+'/' + self.fasta_name + ".misa", self.output_dir+'/' + self.fasta_name + ".misa")
        self.end()

    def run(self):
        """
        运行
        """
        super(SsrTool, self).run()
        if self.option("primer"):
            self.misa()
            self.primer_in()
            self.primer()
            self.primer_out()
            # self.ssr_position()
        else:
            self.misa()
            # self.primer_in()
            # self.primer_out()
        self.set_output()
