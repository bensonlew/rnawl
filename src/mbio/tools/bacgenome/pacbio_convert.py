# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2020/03/31' #hao.gao
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_file


class PacbioConvertAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(PacbioConvertAgent, self).__init__(parent)
        options = [
            {"name": "files", "type": "string", "required": True},
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "lib_type", "type": "string", "default": "Pacbio"},
            {"name": "fq_result_path", "type": "outfile", "format": "sequence.fastq"}
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


class PacbioConvertTool(Tool):
    def __init__(self, config):
        super(PacbioConvertTool, self).__init__(config)
        self.bax2bam = "/bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/bax2bam"
        self.bam2fastq = "/bioinfo/Genomic/Sofware/bedtools2/bin/bamToFastq"
        self.prefix = self.option("sample_name") + "_" + self.option("lib_type")
        self.path = self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin:"
        self.lib = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64:"
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)

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
        cmd = self.bam2fastq + " -i " + input + " -fq " + self.prefix + ".fastq"
        command = self.add_command("bam2fastq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行bam2fastq完成")
        else:
            self.set_error("运行bam2fastq出错")

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
        fq_file = os.path.join(self.work_dir, self.prefix + ".fastq")
        if os.path.isfile(fq_file):
            self.option("fq_result_path", fq_file)

    def run(self):
        super(PacbioConvertTool, self).run()
        if self.option("files").endswith(".bax.h5"):
            self.run_bax2bam()
        elif self.option("files").endswith("bam"):
            self.run_bam2fastq(self.option("files"))
        else:
            self.run_fastq2fastq()
        self.set_output()
        self.end()