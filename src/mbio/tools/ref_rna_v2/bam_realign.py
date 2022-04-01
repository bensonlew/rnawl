# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.12.13

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime

starttime = datetime.datetime.now()


class BamRealignAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(BamRealignAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "name", "type": "string"}  # 生成文件名字，测试文件中名字为DE1_10.realign.bam，其中.realign.bamf为固定。
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("请设置bam路径", code="33709203")
        if not self.option("fa_file"):
            raise OptionError("请设置ref.fa路径", code="33709204")

    def set_resource(self):
        self._cpu = 10
        self._memory = "{}G".format(int(os.path.getsize(self.option('bam_file').path) / 1024.0 ** 3 * 10 +
                   os.path.getsize(self.option('fa_file').path) / 1024.0 ** 3 * 4) + 60)

    def end(self):
        super(BamRealignAgent, self).end()


class BamRealignTool(Tool):
    def __init__(self, config):
        super(BamRealignTool, self).__init__(config)
        self.set_environ(SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/MajorBio_cluster_201.20.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        self.sentieon = "bioinfo/WGS/sentieon-genomics-201808/bin/sentieon"

    def run_bam_realign(self):
        """
        PopLDdecay

        """
        cmd = "{} driver -t 8 -i {} -r {} --algo Realigner {}".format(self.sentieon,
                                                                      self.option("bam_file").prop["path"],
                                                                      self.option("fa_file").prop["path"],
                                                                      self.output_dir + "/" + self.option("name"))
        command = self.add_command("bam_realign", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam_realign运行完成")
        else:
            self.set_error("bam_realign运行失败", code="33709202")

    def run(self):
        super(BamRealignTool, self).run()
        self.run_bam_realign()
        self.end()
