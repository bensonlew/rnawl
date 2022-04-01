#!/usr/bin/env python
# -*- coding: utf-8 -*-
# last modified by shicaiping at 20180507
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import glob
import re
import unittest
from collections import OrderedDict

class LargeFastqStatAgent(Agent):
    """
    用于做fastq序列基本信息统计的工具
    version 1.0
    author: qindanhua
    """

    def __init__(self, parent):
        super(LargeFastqStatAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},
            {"name": "quality", "type": "int", "default": 33}
        ]
        self.add_option(options)
        self._memory_increase_step = 200
        self.step.add_steps('fastq_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fastq_stat.start()
        self.step.update()

    def step_end(self):
        self.step.fastq_stat.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("请说明序列类型，PE or SE?", code = "33705201")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 20
        self._memory = '200G'

    def end(self):
        super(LargeFastqStatAgent, self).end()


class LargeFastqStatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(LargeFastqStatTool, self).__init__(config)
        self.LargeFastqStat_path = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/FastqStat.jar"
        self.fastq_name = self.option("fastq").prop["path"].split("/")[-1]

    def fastq_stat(self):
        file_list = self.get_list_file()
        for n,file in enumerate(file_list):
            fastq_stat = self.add_tool("ref_rna_v3.large.fastq_stat")
            options = {
                "fq_type" : self.option("fq_type"),
                'quality' : self.option("quality"),
                'fastq_info': file
            }
            fastq_stat.set_options(options)
            self.tools.append(fastq_stat)
        for tool in self.tools:
            tool.run()


    def get_list_file(self):
        if self.option("fastq").format == "sequence.fastq":
            with open("fq_list_for_LargeFastqStat", "wb") as w:
                w.write("{}\t{}".format(self.fastq_name, self.option("fastq").prop["path"]))
        elif self.option("fastq").format == "sequence.fastq_dir":
            sample_file = {}
            sample_file =OrderedDict()
            fq_dir = self.option("fastq").prop["path"]
            list_info = os.path.join(fq_dir, "list.txt")
            with open(list_info, "rb") as l:
                for line in l:
                    line = line.strip().split()
                    if line[1] not in sample_file:
                        sample_file[line[1]] = [line[0]]
                    else:
                        sample_file[line[1]].append(line[0])
                self.logger.info(sample_file)

            def yield_sample_sample_list(sample_info_dict,chunksize =50):
                file_info = ""
                n = 0
                for sample in sample_info_dict:
                    if  n == chunksize:
                        yield file_info
                        file_info = ""
                        n=0
                    if len(sample_info_dict[sample]) == 2:
                            file_info += "{}\t{}\t{}\n".format(sample, os.path.join(fq_dir, sample_file[sample][0]),
                                                               os.path.join(fq_dir, sample_file[sample][1]))
                    if len(sample_file[sample]) == 1:
                            file_info += "{}\t{}\n".format(sample, os.path.join(fq_dir, sample_file[sample][0]))
                    n+=1
                if file_info:
                    yield file_info
            for n,file_info in enumerate(yield_sample_sample_list(sample_file)):
                with open("fq_list_for_LargeFastqStat"+"_{}".format(str(n)), "wb") as w:
                    w.write(file_info)
        file_list = glob.glob(self.work_dir +"/fq_list_for_LargeFastqStat_*")
        return file_list

    def set_output(self):
        self.logger.info("set output")
        os.system('rm -rf '+self.output_dir)
        os.system('mkdir '+self.output_dir)
        sample_list = []
        sample_stat_infos ={}
        list_info = os.path.join(self.option("fastq").prop["path"], "list.txt")
        output = self.output_dir+ '/fastq_stat.xls'
        with open (list_info, "r") as f:
            for line in f:
                tmp = line.strip().split()
                if tmp[1] not in sample_list:
                    sample_list.append(tmp[1])
                else:
                    pass
        fatq_stat_files =[]
        for tool in self.tools:
            fatq_stat_file = os.path.join(tool.output_dir,"fastq_stat.xls")
            fatq_stat_files.append(fatq_stat_file)
        with open(fatq_stat_files[0],"r") as  r:
            for line in r:
                if re.search(r'^#Sample_ID', line):
                    header = line
                    break
        with open(output, "w") as w:
            w.write(header)
            for fatq_stat_file in fatq_stat_files:
                with open(fatq_stat_file,"r") as f:
                    for line in r:
                        if re.search(r'^#Sample_ID', line):
                            pass
                        else:
                            line1 = line.strip().split()
                            if line1[0] in sample_list:
                                sample_stat_infos[line1[0]] = line1
            for sample in sample_list:
                w.write(sample_stat_infos[sample])

        self.logger.info("done")
        self.end()

    def run(self):
        super(LargeFastqStatTool, self).run()
        self.on_rely(self.tools,self.set_output)
        self.fastq_stat()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        data = {
            "id": "large_fastq_stat" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "ref_rna_v2.large_fastq_stat",
            "instant": False,
            "options": dict(
                fastq='/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/FastpRna/output/fastq',
                fq_type = "PE"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
