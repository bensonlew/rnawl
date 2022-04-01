# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,shijin'

import glob
import os

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.denovo_rna.qc.fastq_stat import fastq_qual_stat


class DrawFastqInfoAgent(Agent):
    """
    DrawFastqInfo:用于做fastq序列质量统计的工具
    version 1.0
    author: qindanhua
    """

    def __init__(self, parent):
        super(DrawFastqInfoAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq,sequence.fastq_dir"}  # 输入文件fastq序列
        ]
        self.add_option(options)
        self.step.add_steps('quality_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.quality_stat.start()
        self.step.update()

    def step_end(self):
        self.step.quality_stat.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("fastq").is_set:
            raise OptionError("请选择序列文件", code="34000501")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        super(DrawFastqInfoAgent, self).end()


class DrawFastqInfoTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(DrawFastqInfoTool, self).__init__(config)
        self.fastxtoolkit_path = 'bioinfo/seq/fastx_toolkit_0.0.14/'
        self.fastq_name = self.option("fastq").prop['path'].split("/")[-1]

    def run(self):
        super(DrawFastqInfoTool, self).run()
        if self.option("fastq").format == "sequence.fastq":
            fq_path = self.option("fastq").prop['path']
            command = self.fastq_quality_stats(fq_path, "qual_stat")
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行{}完成")
                fastq_qual_stat("qual_stat")
            else:
                if command.return_code is None:
                    command.rerun()
                    self.wait(command)
                    if command.return_code == 0:
                        self.logger.info("运行{}完成".format(command.name))
                        fastq_qual_stat("qual_stat")
                    else:
                        self.set_error("运行%s运行出错!", variables=(command.name,), code="34000502")
                        return False
                else:
                    self.set_error("运行%s运行出错!", variables=(command.name,), code="34000504")
                    return False
        elif self.option("fastq").format == "sequence.fastq_dir":
            commands = self.multi_fastq_quality_stats()
            self.wait()
            for cmd in commands:
                if cmd.return_code == 0:
                    self.logger.info("运行{}完成".format(cmd.name))
                    self.logger.info(os.path.join(self.work_dir, cmd.name[:-14] + "_qual_stat"))
                    fastq_qual_stat(os.path.join(self.work_dir, cmd.name[:-14] + "_qual_stat"))
                else:
                    if cmd.return_code is None:
                        cmd.rerun()
                        self.wait(cmd)
                        if cmd.return_code == 0:
                            self.logger.info("运行{}完成".format(cmd.name))
                            self.logger.info(os.path.join(self.work_dir, cmd.name[:-14] + "_qual_stat"))
                            fastq_qual_stat(os.path.join(self.work_dir, cmd.name[:-14] + "_qual_stat"))
                        else:
                            self.set_error("运行%s运行出错!", variables=(cmd.name,), code="34000506")
                            return False
                    else:
                        self.set_error("运行%s运行出错!", variables=(cmd.name,), code="34000508")
                        return False
        self.set_output()
        self.end()

    def fastq_quality_stats(self, fastq, outfile):
        fastq_name = fastq.split("/")[-1]
        fastq_name = fastq_name.lower()
        cmd = self.fastxtoolkit_path + 'fastx_quality_stats -i {} -Q 33 -o {}'.format(fastq, outfile)
        self.logger.info(cmd)
        self.logger.info("开始运行{}_quality_stats".format(fastq_name))
        command = self.add_command("{}_quality_stats".format(fastq_name), cmd)
        command.run()
        return command

    def multi_fastq_quality_stats(self):
        commands = []
        file_dir = self.option("fastq").prop["path"]
        for f in os.listdir(file_dir):
            if f == "list.txt":
                pass
            else:
                file_path = os.path.join(file_dir, f)
                command = self.fastq_quality_stats(file_path, f.lower() + "_qual_stat")
                commands.append(command)
        return commands

    def set_output(self):
        self.logger.info("set output")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        file_path = glob.glob(r"*qual_stat")
        for f in file_path:
            output_dir = os.path.join(self.output_dir, f)
            os.link(os.path.join(self.work_dir, f), output_dir)
            os.remove(os.path.join(self.work_dir, f))
