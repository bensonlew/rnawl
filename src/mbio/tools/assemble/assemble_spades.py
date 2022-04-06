# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2020/4/01'#gaohao
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_file


class AssembleSpadesAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(AssembleSpadesAgent, self).__init__(parent)
        options = [
            {"name": "PE_list", "type": "infile", "format": "meta.otu.otu_table", "required": True},  # 二代PE reads
            {"name": "MP_list", "type": "infile", "format": "meta.otu.otu_table"},  # 二代MP reads
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "pacbio", "type": "infile", "format": "sequence.fastq"},  # 三代数据clr reads
            {"name": "nanopore", "type": "infile", "format": "sequence.fastq"},
            {"name": "ccs", "type": "infile", "format": "bacgenome.simple_file"},  # 三代数据ccs reads
            {"name": "kmer", "type": "string"},
            {"name": "scf_seq", "type": "outfile", "format": "sequence.fasta"},
            {"name": "cpu", "type": "int"},
            {"name": "mem", "type": "int"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = self.option("cpu") if self.option("cpu") else 16
        self._memory = "%sG" % self.option("mem") if self.option("mem") else "70G"


class AssembleSpadesTool(Tool):
    def __init__(self, config):
        super(AssembleSpadesTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.spades = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/SPAdes-3.11.0-Linux/bin/spades.py"

    def run_spades(self):
        """
        description
        :return:
        """
        cmd = "{} {}".format(self.python_path, self.spades)
        cmd += self.add_pe()
        cmd += self.add_mp()
        cmd += " --pacbio %s " % self.option("pacbio").prop["path"] if self.option("pacbio").is_set else ""
        cmd += " --nanopore %s " % self.option("nanopore").prop["path"] if self.option("nanopore").is_set else ""
        cmd += " -s %s " % self.option("ccs").prop["path"] if self.option("ccs").is_set else ""
        cmd += " -k %s " % self.option("kmer") if self.option("kmer") else ""
        cmd += " --threads %s" % self.option("cpu") if self.option("cpu") else ''
        cmd += " --memory %s " % self.option("mem") if self.option("mem") else ""
        cmd += " -o " + self.option("sample_name")
        command = self.add_command("spades", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("spades运行完成")
        else:
            self.set_error("spades运行出错！")

    def add_pe(self):
        return_cmd = ""
        with open(self.option("PE_list").prop["path"], "r") as file:
            lines = file.readlines()
            if len(lines) == 1:
                line = lines[0].strip().split("\t")
                return_cmd = " -1 %s -2 %s" % (line[1], line[2])
                if len(line) >3:
                    return_cmd  += " -s " + line[3]
                return return_cmd
            for index,line in enumerate(lines):
                line = line.strip().split("\t")
                return_cmd += " --pe%s-1 %s --pe%s-2 %s " % (index + 1, line[1], index + 1, line[2])
                if line[3]:
                    return_cmd += " --pe%s-s %s " % (index+1, line[3])
        return return_cmd

    def add_mp(self):
        return_cmd = ""
        if not self.option("MP_list").is_set:
            return return_cmd
        with open(self.option("MP_list").prop["path"], "r") as file:
            lines = file.readlines()
            for index,line in enumerate(lines):
                line = line.strip().split("\t")
                return_cmd += " --mp%s-1 %s --mp%s-2 %s " % (index+1, line[1], index+1, line[2])
        return return_cmd

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        scaffold = os.path.join(self.work_dir, self.option("sample_name"), "scaffolds.fasta")
        contig = os.path.join(self.work_dir, self.option("sample_name"), "contigs.fasta")
        output = os.path.join(self.output_dir, self.option("sample_name")) + ".scaffold.fna"
        if os.path.isfile(scaffold):
            self.logger.info("%s: set scaffolds.fasta to output dir" % self.option("sample_name"))
            link_file(scaffold, output)
        else:
            self.logger.info("%s: set contigs.fasta to output dir" % self.option("sample_name"))
            link_file(contig, output)
        self.option("scf_seq").set_path(output)

    def run(self):
        super(AssembleSpadesTool, self).run()
        self.run_spades()
        self.set_output()
        self.end()