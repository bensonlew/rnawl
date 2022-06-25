# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# modified 20190320

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class CrisprRessoAgent(Agent):
    """
    接口“变异位点比较分析Tool, 导表在tool中完成”
    """

    def __init__(self, parent):
        super(CrisprRessoAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile",
                "format": "wgs_v2.bam"},
            {"name": "sgrna_region", "type": "infile", "format": "wgs_v2.sgrna_region"},
            {"name": "name", "type": "string"},
            {"name": "ref_fa", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('crispr_resso')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.crispr_resso.start()
        self.step.update()

    def step_end(self):
        self.step.crispr_resso.finish()
        self.step.update()

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("必须输入bam文件")
        if not self.option("name"):
            raise OptionError("必须输入样本名称")
        if not self.option("sgrna_region"):
            raise OptionError("必须输入sgrna_region文件")

    def set_resource(self):
        self._cpu = 4
        self._memory = "20G"

    def end(self):
        super(CrisprRessoAgent, self).end()


class CrisprRessoTool(Tool):

    def __init__(self, config):
        """
        这里还需要解决路径问题，当测试好后，看看把这个文件夹移动到app文件夹下后还是否能正常使用。
        :param config:
        """
        super(CrisprRessoTool, self).__init__(config)
        self.resso = 'miniconda2/bin/CRISPRessoWGS'
        self.set_environ(
            PATH=self.config.SOFTWARE_DIR +
            "/bioinfo/align/bowtie2-2.2.9/")
        # self.set_environ(
        #     PATH="/mnt/ilustre/users/sanger-dev/CRISPResso_dependencies/bin")
        self.set_environ(
            PATH=os.path.dirname(self.config.SOFTWARE_DIR) + "/CRISPResso_dependencies/bin")
        # self.set_environ(PATH="/mnt/ilustre/users/sanger-dev/app/program/vim/bin")
        # self.set_environ(PATH="/mnt/ilustre/users/sanger-dev/biocluster/bin")
        # self.set_environ(PATH="/mnt/ilustre/users/sanger-dev/biocluster/bin")
        # self.set_environ(PATH="/usr/lib64/qt-3.3/bin")
        # self.set_environ(PATH="/usr/local/bin")
        # self.set_environ(PATH="/usr/local/sbin")
        # self.set_environ(PATH="/usr/local/sbin")
        # self.set_environ(PATH="/opt/ganglia/bin")
        # self.set_environ(PATH="/opt/ganglia/sbin")
        # self.set_environ(PATH="/opt/ganglia/sbin")
        # # self.set_environ(PATH="/opt/pdsh/bin")
        # # self.set_environ(PATH="/opt/rocks/bin")
        # # self.set_environ(PATH="/opt/rocks/sbin")
        # # self.set_environ(PATH="/opt/openmpi/bin")
        self.set_environ(PATH="/usr/java/latest/bin")
        # self.set_environ(PATH="/usr/java/latest/bin")
        # self.set_environ(PATH="/mnt/ilustre/users/sanger-dev/.aspera/connect/bin")

    def run_resso(self):

        cmd = "{} -b {} -f {} -r {} -n {} -o {} --save_also_png --keep_intermediate".format(
            self.resso,
            self.option("bam_file"). prop["path"],
            self.option("sgrna_region"). prop["path"],
            os.path.join(
                self.config. SOFTWARE_DIR +
                "/database/dna_geneome",
                self.option("ref_fa")),
            self.option("name"),
            self.output_dir)
        # cmd = "{} -b {} -f {} -r {} -n {} -o {} --save_also_png --keep_intermediate".format(self.resso,
        #                                                                                     self.option("bam_file").
        #                                                                                     prop["path"],
        #                                                                                     self.option("sgrna_region").
        #                                                                                     prop["path"],
        #                                                                                     self.option("ref_fa"),
        #                                                                                     self.option("name"),
        # self.output_dir)
        command = self.add_command("ressowgs", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ressoWGS计算完成")
        else:
            self.set_error("ressoWGS计算失败")

    def run(self):
        super(CrisprRessoTool, self).run()
        self.run_resso()
        self.end()
