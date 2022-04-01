# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/2'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_file


class PacbioConvert2Agent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(PacbioConvert2Agent, self).__init__(parent)
        options = [
            {"name": "files", "type": "string", "required": True},
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "lib_type", "type": "string", "default": "Pacbio"},
            {"name": "bam_result_path", "type": "string", "required": True},
            {"name": "fq_result_path", "type": "string", "required": True}
        ]
        self.add_option(options)
        self.step.add_steps('pacbio_convert')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.pacbio_convert.start()
        self.step.update()

    def step_end(self):
        self.step.pacbio_convert.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        # modified check_option
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class PacbioConvert2Tool(Tool):
    def __init__(self, config):
        super(PacbioConvert2Tool, self).__init__(config)
        self.bax2bam = "/bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/bax2bam"
        self.pbindex = "/bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/pbindex"
        self.bam2fastq = "/bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/bam2fastq"
        self.prefix = self.option("sample_name") + "_" + self.option("lib_type")

    def run_bax2bam(self):
        """
        description
        :return:
        """
        cmd = self.bax2bam + " " + self.option("files") + " -o " + self.prefix + " --output-xml " + self.prefix
        command = self.add_command("bax2bam", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行bax2bam完成")
        else:
            self.set_error("运行bax2bam出错")
        self.run_bam2fastq(self.prefix + ".subreads.bam")

    def run_bam2fastq(self, input):
        if not os.path.isfile(input + ".pbi"):
            cmd = self.pbindex + " " + input
            command = self.add_command("pbindex", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行pbindex完成")
            else:
                self.set_error("运行pbindex出错")
        cmd = self.bam2fastq + " -u -o " + self.prefix + " " + input
        command = self.add_command("bam2fastq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行bam2fastq完成")
        else:
            self.set_error("运行bam2fastq出错")
        link_file(input, os.path.join(self.option("bam_result_path"), self.prefix + ".bam"))
        link_file(input + ".pbi", os.path.join(self.option("bam_result_path"), self.prefix + ".bam.pbi"))

    def run_fastq2fastq(self):
        cmd = "cat " + self.option("files") + " > " + self.prefix + ".fastq"
        self.logger.info("command: %s" % cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("cat done")
        except subprocess.CalledProcessError:
            self.set_error("cat error")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        bam_file = os.path.join(self.work_dir, self.prefix + ".subreads.bam")
        if os.path.isfile(bam_file):
            link_file(bam_file, os.path.join(self.option("bam_result_path"), self.prefix + ".bam"))
            link_file(bam_file + ".pbi", os.path.join(self.option("bam_result_path"), self.prefix + ".bam.pbi"))
        fq_file = os.path.join(self.work_dir, self.prefix + ".fastq")
        if os.path.isfile(fq_file):
            link_file(fq_file, os.path.join(self.option("fq_result_path"), self.prefix + ".fastq"))

    def run(self):
        super(PacbioConvert2Tool, self).run()
        if self.option("files").endswith(".bax.h5"):
            self.run_bax2bam()
        elif self.option("files").endswith("bam"):
            self.run_bam2fastq(self.option("files"))
        else:
            self.run_fastq2fastq()
        self.set_output()
        self.end()