# -*- coding: utf-8 -*-
# __author__ = 'shijin'

"""fq文件序列拆分"""

import os
import datetime
import json
import shutil
import re
from biocluster.workflow import Workflow


class FastqSplitWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FastqSplitWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_fastq", "type": "infile", 'format': "sequence.fastq"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.samples = {}

    def run_fastq_split(self):
        fastq_path = self.option("in_fastq").prop["path"]
        with open(fastq_path,"r") as f:
            for line in f:
                m = re.match("@(.+)_(\d+)", line)
                if not m:
                    self.set_error('fastq文件格式不符合要求', code="12701601")
                sample_name = m.group(1)
                sample = self.return_sample(sample_name)
                sample.add_new_fastq(line, next(f), next(f), next(f))
        for sample in self.samples.values():
            sample.close_all()
        self.logger.info("全部fastq序列处理完毕")
        # self.set_db()
        self.end_this()


    def return_sample(self,sample_name):
        if sample_name in self.samples:
            return self.samples[sample_name]
        sample = Sample(sample_name,self.output_dir)
        self.samples[sample_name] = sample
        return sample



    def run(self):
        self.run_fastq_split()
        super(FastqSplitWorkflow, self).run()


    def end_this(self):
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["venn_table.xls", "xls", "Venn表格"]
        ])
        """
        super(FastqSplitWorkflow, self).end()

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