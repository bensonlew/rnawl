# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/25'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_file


class AssembleVelvetAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(AssembleVelvetAgent, self).__init__(parent)
        options = [
            # pe_list  pe1 pe2 pes insert read_len
            # mp_list  mp1 mp2 insert read_len
            {"name": "PE_list", "type": "infile", "format": "meta.otu.otu_table", "required": True},  # 二代PE reads
            {"name": "MP_list", "type": "infile", "format": "meta.otu.otu_table"},  # 二代MP reads
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "kmer", "type": "string"},
            {"name": "min_contig_lgth", "type": "int", "default": 200, "min": 100},
            {"name": "min_pair_count", "type": "int", "default": 15, "min": 5},
            {"name": "scf_seq", "type": "outfile", "format": "sequence.fasta"},
            {"name": "cpu", "type": "int"},
            {"name": "mem", "type": "int"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = self.option("cpu") if self.option("cpu") else 16
        self._memory = "%sG" % self.option("mem") if self.option("mem") else "50G"


class AssembleVelvetTool(Tool):
    def __init__(self, config):
        super(AssembleVelvetTool, self).__init__(config)
        self.velvet_path = "bioinfo/Genomic/Sofware/velvet_1.2.10/"
        self.pe_num = 0
        self.mp_num = 0
        self.insert_str = ""

    def run_velveth(self):
        """
        description
        :return:
        """
        cmd = self.velvet_path + "velveth result %s " % self.option("kmer")
        with open(self.option("PE_list").prop["path"], "r") as f1:
            lines = f1.readlines()
            self.pe_num = len(lines)
        if self.option("MP_list").is_set:
            with open(self.option("MP_list").prop["path"], "r") as f2:
                lines += f2.readlines()
                self.mp_num = len(lines) - self.pe_num
        for index,line in enumerate(lines):
            line = line.strip().split("\t")
            if index == 0:
                cmd += " -shortPaired -fastq -separate " + " ".join(line[1:3])
                self.insert_str += " -ins_length %s " % line[-2]
            else:
                cmd += " -shortPaired%s -fastq -separate %s %s" % (index +1, line[1], line[2])
                self.insert_str += " -ins_length%s %s " % (index+1, line[-2])
        command = self.add_command("velveth", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("velveth运行完成")
        else:
            self.set_error("velveth运行出错！")

    def run_velvetg(self):
        cmd = self.velvet_path + "velvetg result -cov_cutoff auto -exp_cov auto -min_contig_lgth %s -min_pair_count %s " \
                                 "-scaffolding yes" % (self.option("min_contig_lgth"), self.option("min_pair_count"))
        cmd += self.insert_str
        for i in range(1, self.pe_num + 1):
            cmd += " -shortMatePaired%s yes" % i
        command = self.add_command("velvetg", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("velvetg运行完成")
        else:
            self.set_error("velvetg运行出错！")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        link_file(os.path.join(self.work_dir, "result/contigs.fa"), os.path.join(self.output_dir, "%s.contig.fa" % self.option("sample_name")))
        self.option("scf_seq").set_path(os.path.join(self.output_dir, "%s.contig.fa" % self.option("sample_name")))

    def run(self):
        super(AssembleVelvetTool, self).run()
        self.run_velveth()
        self.run_velvetg()
        self.set_output()
        self.end()