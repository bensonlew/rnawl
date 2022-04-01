# -*- coding: utf-8 -*-
# __author__ = 'qingchenzhang'
import os
import re,math
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class Bowtie2AlignAgent(Agent):
    """
     宏基因组用bowtie2软件将PE文库的read1、read2和singleto序列比对到参考序列上，不建索引
    """

    def __init__(self, parent):
        super(Bowtie2AlignAgent, self).__init__(parent)
        options = [
            {"name": "database", "type": "string"},  # 输入文件,sample.contig.fa
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,s 可不传
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,l
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,r
            {"name": "sam_file", "type": "outfile", "format": "metagbin.sam_dir"},  # 输出文件,map结果的sam路径
        ]
        self.add_option(options)
        self.step.add_steps("Bowtie2_align")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 30

    def stepstart(self):
        self.step.Bowtie2_align.start()
        self.step.update()

    def stepfinish(self):
        self.step.Bowtie2_align.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("database"):
            raise OptionError('必须输入建库文件路径')
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        if self.option("fastq1").is_set:
            size = os.path.getsize(self.option("fastq1").prop["path"])
            n = int(size / (1024 * 1024 * 1024))
            if n > 2:
                self._cpu = 4
                self._memory = str(n*10) + "G"
            else:
                self._cpu = 4
                self._memory = "20G"
        elif self.option("fastqs").is_set:
            size = os.path.getsize(self.option("fastqs").prop["path"])
            n = int(size / (1024 * 1024 * 1024))
            if n > 2:
                self._cpu = 4
                self._memory = str(n * 10) + "G"
            else:
                self._cpu = 4
                self._memory = "20G"

    def end(self):
        super(Bowtie2AlignAgent, self).end()


class Bowtie2AlignTool(Tool):
    def __init__(self, config):
        super(Bowtie2AlignTool, self).__init__(config)
        #self.version = "v2.2.9"
        self.bowtie2_path = '/bioinfo/align/bowtie2-2.2.9/'
        self.index_prefix = ''
        if self.option("fastq1").is_set:
            self.samp_name = os.path.basename(self.option('fastq1').prop['path'])
        elif self.option("fastqs").is_set:
            self.samp_name = os.path.basename(self.option('fastqs').prop['path'].split('.')[0])

    def run_bowtie2_map(self):
        """
        运行bowtie2比对
        先进行pair-end reads比对，当有single reads出现时，进行single reads比对，没有则不进行
        :return:
        """
        self.logger.info('运行bowtie2的双端reads开始比对')
        if self.option('fastq1').is_set and self.option('fastq2').is_set:
            cmd = "{}bowtie2 -p 4 -x {} -1 {} -2 {} -S {}/{}.pair.sam ".format(self.bowtie2_path,
                                                                                         self.option('database'),
                                                                                         self.option("fastq1").path,
                                                                                         self.option("fastq2").path,
                                                                                         self.output_dir, self.samp_name)
            command = self.add_command("bowtie2_map_pair", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("bowtie2_map_pair运行完成")
            else:
                self.set_error("bowtie2_map_pair运行出错！", code="31100502")
        elif self.option("fastqs").is_set:
            self.logger.info("运行bowtie2 比对 single reads开始比对")
            cmd = "{}bowtie2  -p 4 -x {} -U {} -S {}/{}.single.sam".format(self.bowtie2_path,
                                                                                     self.option('database'),
                                                                                     self.option("fastqs").path,
                                                                                     self.output_dir, self.samp_name)
            command = self.add_command("bowtie2_map_single", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("bowtie2_map_single运行完成")
            elif command.return_code == -9:
                self.add_state('memory_limit', 'memory is low!')
            else:
                self.logger.info("return code: %s" % command.return_code)
                self.set_error("bowtie2_map_single运行出错！", code="31100503")
        else:
            self.set_error("缺少输入文件，运行出错！")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('sam_file').set_path(self.output_dir)
        self.logger.info('设置组装拼接分析结果目录成功')

    def run(self):
        """
        运行bowtie2比对程序
        :return:
        """
        super(Bowtie2AlignTool, self).run()
        self.run_bowtie2_map()
        self.set_output()
        self.end()
