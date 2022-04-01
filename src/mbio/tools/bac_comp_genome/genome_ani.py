# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 20119.10.19

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_file

class GenomeAniAgent(Agent):
    """
    比较基因组的ani计算
    """
    def __init__(self, parent):
        super(GenomeAniAgent, self).__init__(parent)
        options = [
            {"name": "seq_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入参考序列文件夹
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"},  #输出ani计算结果
            {"name": "method", "type": "string", "default": "ANIm"} #ANIm,ANIb,ANIblastall,TETRA
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq_dir").is_set:
            raise OptionError("必须设置参数seq_dir文件夹!")
        if not self.option("method"):
            raise OptionError("必须设置参数method!")
        else:
            if self.option("method") not in ["ANIm", "ANIb", "ANIblastall", "TETRA"]:
                raise OptionError("必须设置正确参数method!")

    def set_resource(self):
        num =len(os.listdir(self.option("seq_dir").prop["path"]))
        self._cpu = 4
        self._memory = str(num*3)+'G'

    def end(self):
        super(GenomeAniAgent, self).end()

class GenomeAniTool(Tool):
    def __init__(self, config):
        super(GenomeAniTool, self).__init__(config)
        lib_path = os.path.join(self.config.SOFTWARE_DIR, "program/Python35/lib")
        path =self.config.SOFTWARE_DIR + "/program/Python35/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/MUMmer3.23:" + self.config.SOFTWARE_DIR + "/bioinfo/align/ncbi-blast-2.3.0+/bin:" + self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/blast-2.2.18/bin"
        self.set_environ(PATH=path, LD_LIBRARY_PATH=lib_path)
        self.ani = "program/Python35/bin/average_nucleotide_identity.py"

    def run_ani(self):
        if os.path.exists(self.work_dir + "/out"):
            shutil.rmtree(self.work_dir + "/out")
        cmd = "{} -i {} -m {} -o {} -g --workers 4".format(self.ani, self.option("seq_dir").prop["path"], self.option("method"), self.work_dir + "/out")
        command = self.add_command("run_ani", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_ani运行完成！")
        else:
            self.set_error("run_ani运行完成运行出错!")

    def set_output(self):
        if self.option("method") in ["ANIm"]:
            link_file(self.work_dir + "/out/ANIm_percentage_identity.tab", self.output_dir + "/all.ANI_summary.xls")
        elif self.option("method") in ["TETRA"]:
            link_file(self.work_dir + "/out/TETRA_correlations.tab", self.output_dir + "/all.ANI_summary.xls")
        elif self.option("method") in ["ANIblastall"]:
            link_file(self.work_dir + "/out/ANIblastall_percentage_identity.tab", self.output_dir + "/all.ANI_summary.xls")
        elif self.option("method") in ["ANIb"]:
            link_file(self.work_dir + "/out/ANIb_percentage_identity.tab", self.output_dir + "/all.ANI_summary.xls")
        self.option("out", self.output_dir + "/all.ANI_summary.xls")

    def run(self):
        super(GenomeAniTool, self).run()
        self.run_ani()
        self.set_output()
        self.end()