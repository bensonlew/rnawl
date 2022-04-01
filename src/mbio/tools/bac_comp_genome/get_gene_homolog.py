# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.10.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file


class GetGeneHomologAgent(Agent):
    """
    同源基因的进化树分析
    """

    def __init__(self, parent):
        super(GetGeneHomologAgent, self).__init__(parent)
        options = [
            {"name": "type", "type": "string", "default": "port"},  # 核酸nul or 蛋白port
            {"name": "homolog", "type": "infile", "format": "sequence.profile_table"}, #同源蛋白聚类cluster表
            {"name": "clusters", "type": "string"},#逗号隔开的cluster如：Cluster1，cluster2
            {"name": "gene_path", "type": "string"},
            {"name": "sample_list", "type": "string"},  # 文件中有样品信息
            {"name": "out", "type": "outfile", "format": "sequence.fasta"}, #比对对齐的序列文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("homolog").is_set:
            raise OptionError("必须输入homolog文件！")
        if not self.option("sample_list"):
            raise OptionError("必须输入sample_list文件！")
        if not self.option("clusters"):
            raise OptionError("必须输入clusters文件！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GetGeneHomologAgent, self).end()


class GetGeneHomologTool(Tool):
    def __init__(self, config):
        super(GetGeneHomologTool, self).__init__(config)
        self._version = "1.0"
        self.homolog = self.option("homolog").prop["path"]
        self.sample_list = self.option("sample_list")
        self.type = self.option("type")
        self.clusters = self.option("clusters")
        self.gene_dir = self.option("gene_path")
        self.python_path = "/program/Python/bin/python"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'
        self.out = self.work_dir + "/all.homolog.fna"

    def run(self):
        """
        运行
        :return:
        """
        super(GetGeneHomologTool, self).run()
        self.run_get_fasta()
        self.set_output()
        self.end()

    def run_get_fasta(self):
        cmd = '{} {}get_gene_homolog.py -t {} -i {} -list {} -dir {} -s {} -o {}'.format(self.python_path, self.package_path, self.type, self.homolog, self.clusters, self.gene_dir, self.sample_list, self.out)
        command = self.add_command("run_get_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_get_fasta运行完成！")
        else:
            self.set_error("run_get_fasta运行完成运行出错!")

    def set_output(self):
        link_file(self.out, self.output_dir + "/all.homolog.fna")
        self.option("out", self.output_dir + "/all.homolog.fna")