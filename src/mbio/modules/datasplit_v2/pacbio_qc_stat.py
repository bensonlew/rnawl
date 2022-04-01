# -*- coding: utf-8 -*-
# __author__ = "XueQinwen"
# last_modify: 20210908

import os
import re
import json
import time
import unittest
from biocluster.module import Module
from biocluster.core.exceptions import OptionError

class PacbioQcStatModule(Module):
    '''
    整合三代拆分质控和统计部分
    '''
    def __init__(self, work_id):
        super(PacbioQcStatModule, self).__init__(work_id)
        options = [
            {"name":"split_dir",'type':"infile","format":"denovo_rna_v2.common_dir"},
            {"name":"index_path","type":"string"},
            ]
        self.add_option(options)
        self.tools = [] 
        self.qc_tools = []
        self.stat_options = []
    
    def check_option(self):
        """
        参数检查
        """
        if not self.option("split_dir"):
            raise OptionError("没有获取到拆分数据")
        if not self.option("index_path"):
            raise OptionError("没有获取到index文件")
    
    def get_index(self):
        self.index = {}
        with open(self.option("index_path"),'r') as ip:
            while 1:
                line = ip.readline()
                if not line:
                    break 
                fd = line.rstrip().split('\t')
                self.index[fd[0]]={}
                self.index[fd[0]][1]= fd[1]
                self.index[fd[0]][2]= fd[2]
                self.index[fd[0]]['type'] = fd[3]
                self.index[fd[0]]['primer_type'] = fd[4]
                self.index[fd[0]]['sn_name'] = fd[5]
    
    def run(self):
        super(PacbioQcStatModule, self).run()
        self.get_index()
        for i in self.index.keys():
            sample_path = os.path.join(self.option("split_dir").prop['path'],"split.{}.bam".format(i))
            if not os.path.exists(sample_path):
                continue
            if self.index[i]["type"] == "diversity":
                self.run_qc(sample_path,i)
            else:
                options_stat = {
                    "input_bam":sample_path,
                    # "raw_fastq":os.path.join(qc_output,"{}.ccs.fastq".format(i)),
                    "project_type":"non_diversity",
                    # "clean_fastq":os.path.join(qc_output,"{}_value.fastq".format(i))
                }
                self.stat_options.append(options_stat)
                # self.run_stat(options_stat)
        if len(self.qc_tools) > 0:
            if len(self.qc_tools)==1:
                self.qc_tools[0].on('end',self.run_stat)
            else:
                self.on_rely(self.qc_tools,self.run_stat)
            for tool in self.qc_tools:
                tool.run()
        else:
            self.run_stat()
                
    def run_qc(self,sample_path,i):
        qc = self.add_tool('datasplit_v2.pacbio_qc')
        primer_type = self.index[i]['primer_type']
        options = {
            "sn_name":self.index[i]['sn_name'],
            "input_bam":sample_path,
            "primer_type":primer_type,
        }
        qc.set_options(options)
        qc_output = qc.output_dir
        options_stat = {
            "input_bam":sample_path,
            "raw_fastq":os.path.join(qc_output,"{}.ccs.fastq.gz".format(i)),
            "project_type":"diversity",
            "clean_fastq":os.path.join(qc_output,"{}_value.fastq.gz".format(i))
        }
        self.stat_options.append(options_stat)
        self.qc_tools.append(qc)
        # qc.on('end',self.run_stat,options_stat)

    def run_stat(self):
        for i in self.stat_options:
            stat = self.add_tool('datasplit_v2.pacbio_sample_summary')
            stat.set_options(i)
            self.tools.append(stat)
        if len(self.tools)==1:
            self.tools[0].on('end',self.set_output)
        else:
            self.on_rely(self.tools,self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        stat_dir = os.path.join(self.output_dir, "statistic")
        if not os.path.exists(stat_dir):
            os.mkdir(stat_dir)
        for i in self.tools:
            output_dir = i.output_dir
            file_name = ""
            for _,_,files in os.walk(output_dir):
                for name in files:
                    file_name = name
            os.link(os.path.join(output_dir,file_name),
                        os.path.join(stat_dir, file_name))
        data_dir =os.path.join(self.output_dir, "data")
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        for i in self.index.keys():
            if self.index[i]["type"] != "diversity":
                file = os.path.join(self.option("split_dir").prop['path'],"split.{}.bam".format(i))
                file_name = "{}.ccs.bam".format(i)
                try:
                    os.link(file,os.path.join(data_dir,file_name))
                except:
                    continue
        for i in self.qc_tools:
            qc_dir = i.output_dir
            for _,_,files in os.walk(qc_dir):
                for name in files:
                    os.link(os.path.join(qc_dir,name),
                                os.path.join(data_dir, name))
        self.end()

    def end(self):
        super(PacbioQcStatModule,self).end()

class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id':'qc_stat_' + str(random.randint(1,10000)),
            "type": "module",
            "name": "datasplit_v2.pacbio_qc_stat",
            "options": {
                "split_dir":"/mnt/ilustre/users/sanger-dev/workspace/20210915/Single_lima_9691/Lima/output/split",
                "index_path":'/mnt/ilustre/users/sanger-dev/workspace/20210915/PacbioSplit_20181016PE300_20210915_105314/sample_list.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()