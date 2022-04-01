# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2018/12/25'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import check_command
from mbio.packages.metagbin.common_function import Fasta


class AniAgent(Agent):
    """
    fastANI is a fast alignment-free implementation for computing whole-genome
    Average Nucleotide Identity (ANI) between genomes
    Example usage:
    fastANI -q genome1.fa -r genome2.fa -o output.txt
    fastANI -q genome1.fa --rl genome_list.txt -o output.txt
    version 1.0

    """

    def __init__(self, parent):
        super(AniAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            # query genome (fasta/fastq)[.gz]
            {"name": "query_list", "type": "string"},
            # a file containing list of query genome files, one genome per line
            {"name": "ref", "type": "infile", "format": "sequence.fasta"},
            # reference genome (fasta/fastq)[.gz]
            {"name": "ref_list", "type": "string"},
            # a file containing list of reference genome files, one genome per line
            {"name": "kmer", "type": "int", "default": 16},
            # kmer size <= 16 [default 16]
            {"name": "fraglen", "type": "int", "default": 3000},
            # fragment length [default : 3,000]
            {"name": "minfrag", "type": "int", "default": 50},
            # minimum fragments for trusting ANI [default : 50]
            {"name": "visualize", "type": "string", "default": "disabled"},
            # output mappings for visualization, can be enabled for single genome to
            # single genome comparison only [disabled by default]
            {"name": "outfile", "type": "outfile", "format": "sequence.profile_table"}
            # output file name
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("query").is_set and not self.option("query_list"):
            raise OptionError("必须输入查询文件")
        if not self.option("ref").is_set and not self.option("ref_list"):
            raise OptionError("必须输入参考文件")
        if self.option("kmer") > 16:
            raise OptionError("kmer不能超过16")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "5G"


class AniTool(Tool):
    def __init__(self, config):
        super(AniTool, self).__init__(config)
        lib_path = os.path.join(self.config.SOFTWARE_DIR, "library/gsl23/lib")
        lib64_path = os.path.join(self.config.SOFTWARE_DIR, "gcc/5.1.0/lib64")
        env_path = lib_path + ":" + lib64_path
        self.set_environ(LD_LIBRARY_PATH=env_path)
        self.logger.debug("export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH" % env_path)
        self.ani = "bioinfo/metaGenomic/FastANI-1.0/fastANI"

    def run_ani(self):
        """
        description
        :return:
        """
        cmd = "%s " % self.ani
        if self.option("ref").is_set:
            ref_path = self.check_file(self.option("ref").path)
            if not ref_path:
                self.set_error("参考序列n50少于10Kbp")
            cmd = cmd + " -r %s" % ref_path
        else:
            cmd = cmd + " --refList %s" % self.check_list(self.option("ref_list"))
        if self.option("query").is_set:
            query_path = self.check_file(self.option("query").path)
            if not query_path:
                self.logger.info("查询序列N50小于10Kbp")
            cmd = cmd + " -q %s" % self.check_file(self.option("query").path)
        else:
            cmd = cmd + " -ql %s" % self.check_list(self.option("query_list"))
        cmd += " -o %s -k %s --fragLen %s --minFrag %s --visualize %s" % (os.path.join(self.work_dir, "output.txt"),
        self.option("kmer"), self.option("fraglen"), self.option("minfrag"), self.option("visualize"))
        command = self.add_command("fastani", cmd, ignore_error=True).run()
        self.wait(command)
        def success():
            self.logger.info("fastani运行完成")
        def fail():
            self.set_error("fastani运行出错")
        check_command(self, command, [0], [-9], success, fail, memory_limit_fun=None)

    def check_file(self, file_path):
        obj = Fasta(file_path)
        obj.parse()
        n50 = obj.calculate_nxx(50)
        self.logger.info("check N50: %s %s %s" % (file_path, n50, obj.max_len))
        if obj.check_n50(10000):
            return file_path
        else:
            return file_path  # 理论上不允许，暂时改成这样，应该提交一个错误信息到页面

    def check_list(self, list_path):
        with open(list_path, "r") as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip()
                if os.path.isfile(line):
                    value = self.check_file(line)
                    if value:
                        pass  # 需要写一个新list文件
                    else:
                        self.logger.info("文件N50低于10Kbp，忽略比较:%s" % line)
                else:
                    self.set_error("没有找到文件%s" % line)
        return list_path

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.option("outfile", os.path.join(self.work_dir, "output.txt"))

    def run(self):
        super(AniTool, self).run()
        self.run_ani()
        self.set_output()
        self.end()