# -*- coding: utf-8 -*-
# __author__ = 'shijin'

"""fq文件序列拆分"""
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re

class SplitFastqAgent(Agent):

    def __init__(self, parent):
        super(SplitFastqAgent, self).__init__(parent)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq"}
        ]
        self.add_option(options)
        self.step.add_steps('split_fastq')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.split_fastq.start()
        self.step.update()

    def step_end(self):
        self.step.split_fastq.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("in_fastq").is_set:
            raise OptionError("请传入fastq序列文件")
        else:
            return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(SplitFastqAgent, self).end()

class SplitFastqTool(Tool):

    def __init__(self, config):
        super(SplitFastqTool, self).__init__(config)
        self.samples = {}

    def run_fastq_split(self):
        fastq_path = self.option("in_fastq").prop["path"]
        with open(fastq_path,"r") as f:
            for line in f:
                m = re.match("@(.+)_(\d+)", line)
                if not m:
                    raise Exception('fastq文件格式不符合要求')
                sample_name = m.group(1)
                sample = self.return_sample(sample_name)
                sample.add_new_fastq(line, next(f), next(f), next(f))
        for sample in self.samples.values():
            sample.close_all()
        self.logger.info("全部fastq序列处理完毕")

    def return_sample(self,sample_name):
        if sample_name in self.samples:
            return self.samples[sample_name]
        sample = Sample(sample_name,self.output_dir)
        self.samples[sample_name] = sample
        return sample

    def run(self):
        super(SplitFastqTool, self).run()
        self.run_fastq_split()
        self.end()

class Sample(object):
    def __init__(self, name,dir):
        self.name = name
        self._new_fastq_file = open(dir+ "/" + self.name + '.fastq', 'w')

    def add_new_fastq(self, line1,line2,line3,line4):
        self._new_fastq_file.write(line1)
        self._new_fastq_file.write(line2)
        self._new_fastq_file.write(line3)
        self._new_fastq_file.write(line4)


    def close_all(self):
        self._new_fastq_file.close()