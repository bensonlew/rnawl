# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.10.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file
import os
from Bio import SeqIO


class GetHousekeepingAgent(Agent):
    """
    获取单个或多个看家基因的序列，并连接
    """

    def __init__(self, parent):
        super(GetHousekeepingAgent, self).__init__(parent)
        options = [
            {"name": "type", "type": "string", "default": "port"},  # 核酸nul or 蛋白port
            {"name": "core_gene_blast", "type": "infile", "format": "sequence.profile_table"}, #同源蛋白聚类cluster表
            {"name": "core_names", "type": "string"},
            {"name": "sample_list", "type": "string"},  # 文件中有样品信息
            {"name": "out_group", "type": "string"},
            {"name": "out", "type": "outfile", "format": "sequence.fasta"}, #比对对齐的序列文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("core_gene_blast").is_set:
            raise OptionError("必须输入homolog文件！")
        if not self.option("sample_list"):
            raise OptionError("必须输入sample_list文件！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GetHousekeepingAgent, self).end()


class GetHousekeepingTool(Tool):
    def __init__(self, config):
        super(GetHousekeepingTool, self).__init__(config)
        self._version = "1.0"
        self.core_gene_blast = self.option("core_gene_blast").prop["path"]
        self.sample_list = self.option("sample_list")
        self.core_names = self.option("core_names")
        self.python_path = "/program/Python/bin/python"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'
        self.out = self.work_dir + "/all.housekeeping.fna"

    def run(self):
        """
        运行
        :return:
        """
        super(GetHousekeepingTool, self).run()
        if self.option("out_group"):
            self.run_get_fasta()
            self.get_outgroup()
            self.set_output()
            self.end()
        else:
            self.run_get_fasta()
            self.set_output()
            self.end()

    def run_get_fasta(self):
        cmd = '{} {}get_housekeeping.py  -i {} -list {} -s {} -o {}'.format(self.python_path, self.package_path, self.core_gene_blast, self.core_names, self.sample_list, self.out)
        command = self.add_command("run_get_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_get_fasta运行完成！")
        else:
            self.set_error("run_get_fasta运行完成运行出错!")

    def get_outgroup(self):
        seqs = ''
        for name in self.option("core_names").split(","):
            if os.path.exists(self.config.SOFTWARE_DIR + "/database/MajorbioDB/" + self.option("out_group") + "/house_keeping/core_gene/" + name + ".fa"):
                os.link(self.config.SOFTWARE_DIR + "/database/MajorbioDB/" + self.option("out_group") + "/house_keeping/core_gene/" + name + ".fa", self.work_dir + "/" + name + ".fa")
                uniprot_iterator = SeqIO.parse(self.work_dir + "/" + name + ".fa", "fasta")
                for i in uniprot_iterator:
                    if i.id == name:
                        seqs += i.seq
        with open(self.work_dir + "/outgroup.fasta", "w") as f:
            f.write(">{}\n{}\n".format(self.option("out_group"), seqs))
        os.system("cat {} >> {}".format(self.work_dir + "/outgroup.fasta", self.out))

    def set_output(self):
        link_file(self.out, self.output_dir + "/all.housekeeping.fna")
        self.option("out", self.output_dir + "/all.housekeeping.fna")