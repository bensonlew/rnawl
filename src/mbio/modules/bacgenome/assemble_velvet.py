# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/25'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file


class AssembleVelvetModule(Module):
    """
    velvet拼接
    """

    def __init__(self, work_id):
        super(AssembleVelvetModule, self).__init__(work_id)
        option = [
            # pe_list  pe1 pe2 pes insert read_len
            {"name": "PE_list", "type": "infile", "format": "meta.otu.otu_table", "required": True},
            # mp_list  mp1 mp2 insert read_len
            {"name": "MP_list", "type": "infile", "format": "meta.otu.otu_table"},  # 可选MPreads
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "kmers", "type": "string", "required": True},
            {"name": "min_contig_lgth", "type": "int", "default": 200, "min": 100},
            {"name": "min_pair_count", "type": "int", "default": 15, "min": 5},
            {"name": "scafSeq", "type": "outfile", "format": "sequence.fasta_dir"}
        ]
        self.add_option(option)
        self.run_tools = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        return True

    def run(self):
        super(AssembleVelvetModule, self).run()
        self.run_velvet()

    def run_velvet(self):
        kmers = self.option("kmers").split(",")
        for k in kmers:
            tool = self.add_tool("assemble.assemble_velvet")
            tool.set_options({
                "PE_list": self.option("PE_list"),
                "MP_list": self.option("MP_list"),
                "sample_name": self.option("sample_name"),
                "kmer": k,
                "min_contig_lgth": self.option("min_contig_lgth"),
                "min_pair_count": self.option("min_pair_count")
            })
            self.run_tools.append(tool)
        self.on_rely(self.run_tools, self.set_output)
        for tool in self.run_tools:
            tool.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.run_tools:
            link_file(tool.option("scf_seq").prop["path"], self.output_dir +
                      "/%s.%skmer.scafSeq" % (self.option("sample_name"), tool.option("kmer")))
        self.option("scafSeq").set_path(self.output_dir)
        self.logger.info("设置结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AssembleVelvetModule, self).end()
