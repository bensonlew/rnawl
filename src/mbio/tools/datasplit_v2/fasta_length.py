# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171206

"""fasta_clipping_histogram_unique 用于统计fasta序列的长度"""
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import glob
import os


class FastaLengthAgent(Agent):
    def __init__(self, parent):
        super(FastaLengthAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # fasta文件
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},  # fasta文件夹
            {"name": "list_file", "type": "infile", "format": "sequence.file_sample"}  #  fasta文件夹对应的list
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("fasta").is_set and not self.option("fasta_dir").is_set:
            raise OptionError("请传入fasta文件或者fasta文件夹")

    def set_resource(self):
        # self._cpu = 2
        # self._memory = "10G"
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(FastaLengthAgent, self).end()


class FastaLengthTool(Tool):
    def __init__(self, config):
        super(FastaLengthTool, self).__init__(config)
        self.perl_path = '/miniconda2/bin/perl '
        self.length_path = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/fasta_clipping_histogram_unique.pl"

    def single_fasta_length_stat(self):
        """统计单个fasta文件的序列长度"""
        fa_path = self.option("fasta").prop["path"]
        out_path = os.path.join(self.work_dir, os.path.basename(fa_path).split(".")[0] + "_length.xls")
        cmd = "{} {} {} {}".format(self.perl_path, self.length_path, fa_path, out_path)
        self.logger.info("开始进行fasta的长度统计")
        command = self.add_command("fasta_length", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fasta的长度统计运行完成")
        else:
            self.set_error("fasta的长度统计运行出错")

    def many_fasta_length_stat(self):
        """统计多个fasta文件的序列长度"""
        self.logger.info("开始进行fasta的长度统计")
        cmds = []
        if self.fa_files:
            for i in range(len(self.fa_files)):
                out_path = os.path.basename(self.fa_files[i]).split(".")[0] + "_length.xls"
                cmd = "{} {} {} {}".format(self.perl_path, self.length_path, self.fa_files[i], out_path)
                command = self.add_command("fasta_length_{}".format(i), cmd).run()
                cmds.append(command)
        self.wait()
        for cmd in cmds:
            if cmd.return_code == 0:
                self.logger.info("运行{}完成".format(cmd.name))
            else:
                self.set_error("运行{}出错".format(cmd.name))
        self.logger.info("fasta的长度统计运行完成")

    def get_many_fasta_files(self):
        """得到多个fasta文件的list"""
        fa_dir = self.option("fasta_dir").prop["path"]
        self.fa_files = []
        list_file = os.path.join(fa_dir, "list.txt")
        if self.option("list_file").is_set:
            with open(self.option("list_file").prop["path"], "r") as f:
                for line in f:
                    item = line.strip().split()
                    fa = os.path.join(fa_dir, item[0])
                    if os.path.isfile(fa):
                        self.fa_files.append(fa)
                    else:
                        self.logger.info(fa)
        elif os.path.exists(list_file):
            with open(list_file, "r") as f:
                for line in f:
                    item = line.strip().split("\t")
                    fa = os.path.join(fa_dir, item[0])
                    if os.path.isfile(fa):
                        self.fa_files.append(fa)
        else:
            for f in os.listdir(fa_dir):
                fa = os.path.join(fa_dir, f)
                if os.path.isfile(fa):
                    self.fa_files.append(fa)
                elif os.path.isdir(fa):
                    for f1 in os.listdir(fa):
                        fa1 = os.path.join(fa, f1)
                        if os.path.isfile(fa1):
                            self.fa_files.append(fa1)
                        else:
                            self.set_error("{}不是一个文件,fasta_dir里只能有一层文件夹,请检查".format(fa1))
                else:
                    self.set_error("{}既不是文件又不是文件夹,请检查".format(fa))

    def set_output(self):
        self.logger.info("set output")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        outflies = glob.glob(r"*_length.xls")
        for f in outflies:
            f1 = os.path.join(self.work_dir, f)
            f2 = os.path.join(self.output_dir, f)
            if os.path.exists(f2):
                os.remove(f2)
            os.link(f1, f2)
        self.logger.info("done")

    def run(self):
        super(FastaLengthTool, self).run()
        if self.option("fasta_dir").is_set:
            self.get_many_fasta_files()
            self.many_fasta_length_stat()
        else:
            self.single_fasta_length_stat()
        self.set_output()
        self.end()
