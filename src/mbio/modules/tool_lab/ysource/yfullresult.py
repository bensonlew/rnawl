# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest
import random
import datetime
import glob
import subprocess
import re
import os
import sys
import shutil

class YfullresultModule(Module):
    def __init__(self, work_id):
        super(YfullresultModule, self).__init__(work_id)
        self.options = [
            {"name":"yfull_list","type":"string"},
            {"name":"private_list","type":"string"},
            {"name":"source_list","type":"string"},
            {"name":"all_ot",'type':"string"},
            {"name":"all_dep","type":"string"},
            {"name":"size_file","type":"string"}
            # {"name":"pos_list","type":"string"},
            # {"name":"tree","type":"string"},
            # {"name":"snp","type":"string"}
        ]
        self.add_option(self.options)
        self.txttoexcel = self.add_tool("tool_lab.ysource.txttoexcel")
        self.sample_list = []
        self.yfull = {}
        self.private = {}
        self.source_result = {}
        self.tools = []
        # self.get_input()

    def check_option(self):
        """
        参数检查
        """
        for i in self.options:
            op = self.option(i['name'])
            if not os.path.exists(op):
                raise OptionError("找不到{}".format(op))
    
    def get_input(self):
        with open(str(self.option("yfull_list")),"r") as yl:
            while 1:
                line = yl.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                self.yfull[fd[0]] = fd[1]
                self.sample_list.append(fd[0])
        with open(str(self.option("private_list")),"r") as pl:
            while 1:
                line = pl.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                self.private[fd[0]] = fd[1]

        with open(str(self.option("source_list")),"r") as sl:
            while 1:
                line = sl.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                self.source_result[fd[0]] = fd[1]


    def run(self):
        super(YfullresultModule,self).run()
        self.run_tools()

    def run_tools(self):
        self.get_input()
        for sn in self.sample_list:
            self.run_vcf2snp(sn, self.yfull[sn], self.private[sn],self.source_result[sn])
        if len(self.tools) == 1:
            self.tools[0].on("end", self.run_txtToExceltools)
        else:
            self.on_rely(self.tools, self.run_txtToExceltools)
        for tool in self.tools:
            tool.run()

    def run_vcf2snp(self,sn,yfull, private, source_result):
        source = os.path.join(source_result, "database_result/{}.txt".format(sn))
        vcf2snp = self.add_tool("tool_lab.ysource.vcftosnp")
        # vcf2snp = self.add_tool("tool_lab.ysource.add_database")
        options = {
            "yfull_vcf":yfull ,
            "sample_txt": source,
            "private_vcf": private,
            # "snp":self.option("snp"),
            # "pos":self.option("pos_list"),
            # "tree":self.option("tree"),
            "sample_name": sn
        }
        vcf2snp.set_options(options)
        self.tools.append(vcf2snp)

    def run_txtToExceltools(self):
        txtinput = os.path.join(self.work_dir,"all_sample_result.txt")
        with open(txtinput,"w") as asr:
            for sn in self.source_result.keys():
                dir_path = self.source_result[sn]
                with open(os.path.join(dir_path, "all_sample_result.txt"),'r') as asr_in:
                    asr.write(asr_in.readline())
        options = {
            "input_file" : txtinput,
            "output_name" : "all_sample_result.xlsx"
        }
        self.txttoexcel.set_options(options)
        self.txttoexcel.on("end", self.set_output)
        self.txttoexcel.run()
        

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(
            self.__class__.__name__))
        
        for tool in self.tools:

            out_dir = tool.output_dir
            sample_name = tool.option('sample_name')
            try:
                shutil.copytree(out_dir,os.path.join(self.output_dir,sample_name))
            except Exception as e:
                self.logger.info("设置结果目录失败{}".format(e)) 
        all_ot = self.option('all_ot')
        all_dep = self.option('all_dep')
        size_file = self.option('size_file')
        os.link(all_ot, os.path.join(self.output_dir, "all_sample_MAP_OT.xls"))
        os.link(all_dep, os.path.join(self.output_dir, "all_sample_dep.xls"))
        os.link(size_file, os.path.join(self.output_dir, "size.xls"))
        os.link(os.path.join(self.txttoexcel.output_dir,"all_sample_result.xlsx"),os.path.join(self.output_dir,"all_sample_result.xlsx"))
        self.end()
    
    def end(self):
        super(YfullresultModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "yfullresult_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.ysource.yfullresult",
            "options": dict(
                yfull_list="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_bamtovcf_6434/Bamtovcf/output/yfull_list",
                private_list = "/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_bamtovcf_6434/Bamtovcf/output/private_list",                
                source_list="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_yfullsource_4143/Yfullsource/output/yfullsource_list",
                all_ot="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_yfullsource_4143/Yfullsource/output/all_sample_MAP_OT.xls",
                all_dep="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_yfullsource_4143/Yfullsource/output/all_sample_dep.xls",
                size_file="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_Mapping_1564/Mapping/output/size.xls",
                # pos_list="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/file/posdb",
                # tree= "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/file/tree",
                # snp="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/file/SOGG_snp.xls"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
  

