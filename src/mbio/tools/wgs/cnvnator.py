# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.10

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CnvnatorAgent(Agent):
    """
    软件: cnvnator
    call cnv
    """
    def __init__(self, parent):
        super(CnvnatorAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.sam"},  # bam文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file").is_set:
            raise OptionError("请设置sam文件", code="34501701")

    def set_resource(self):
        self._cpu = 2
        self._memory = "50G"

    def end(self):
        super(CnvnatorAgent, self).end()


class CnvnatorTool(Tool):
    def __init__(self, config):
        super(CnvnatorTool, self).__init__(config)
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin") # 线上环境
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/root-6.12.06/build/bin/")
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/root-6.12.06/build/lib/")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/root/lib")  # nb
        self.cnvnator = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/CNVnator_v0.3/src/cnvnator"
        self.cnvnator_path = "bioinfo/gene-structure/CNVnator_v0.3/src/cnvnator"
        self.cnvnator_call = "bioinfo/WGS/cnvnator_call.sh"

    def run_cnvnator_tree(self):
        """
        cnvnator tree 方法，得到root文件
        """
        self.sample_name = os.path.basename(self.option("bam_file").prop["path"]).split(".")[0]
        self.root_name = os.path.join(self.work_dir, self.sample_name + ".root")
        cmd = "{} -root {} -tree {}".format(self.cnvnator_path, self.root_name, self.option("bam_file").prop["path"])
        command = self.add_command("cnvnator_tree", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cnvnator tree运行完成")
        else:
            self.set_error("cnvnator tree运行失败", code="34501701")

    def run_cnvnator_his(self):
        """
        cnvnator his 方法
        """
        cmd = "{} -root {} -his 300 -d".format(self.cnvnator_path, self.root_name)
        command = self.add_command("cnvnator_his", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cnvnator his运行完成")
        else:
            self.set_error("cnvnator his运行失败", code="34501702")

    def run_cnvnator_stat(self):
        """
        cnvnator stat 方法
        """
        cmd = "{} -root {} -stat 300 -d".format(self.cnvnator_path, self.root_name)
        command = self.add_command("cnvnator_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cnvnator stat运行完成")
        else:
            self.set_error("cnvnator stat运行失败", code="34501703")

    def run_cnvnator_partition(self):
        """
        cnvnator partition 方法
        """
        cmd = "{} -root {} -partition 300 -d".format(self.cnvnator_path, self.root_name)
        command = self.add_command("cnvnator_partition", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cnvnator partition运行完成")
        else:
            self.set_error("cnvnator partition运行失败", code="34501704")

    def run_call_cnv(self):
        """
        call cnv
        """
        call_cnv = os.path.join(self.work_dir, self.sample_name + ".cnv")
        cmd = "{} {} -root {} -call 300 {}".format(self.cnvnator_call, self.cnvnator, self.root_name, call_cnv)
        command = self.add_command("cnvnator_call_cnv", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cnvnator call cnv运行完成")
        else:
            self.set_error("cnvnator call cnv运行失败", code="34501705")
        if os.path.exists(os.path.join(self.output_dir, self.sample_name + ".cnv.xls")):
            os.remove(os.path.join(self.output_dir, self.sample_name + ".cnv.xls"))
        os.link(call_cnv, os.path.join(self.output_dir, self.sample_name + ".cnv.xls"))

    def run(self):
        super(CnvnatorTool, self).run()
        self.run_cnvnator_tree()
        self.run_cnvnator_his()
        self.run_cnvnator_stat()
        self.run_cnvnator_partition()
        self.run_call_cnv()
        self.end()
