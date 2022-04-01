# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest
import random
import datetime
import subprocess
import re
import glob
import os
import sys
import shutil

class DatapreModule(Module):
    '''
    Yfull流程整合输入数据输出成列表
    '''
    def __init__(self, work_id):
        super(DatapreModule, self).__init__(work_id)
        options = [
            {'name':'bam_list','type':'infile',"format":"denovo_rna_v2.bamlist"},
        ]
        self.add_option(options)
        self.bam_file = {}
        self.bam2fqgz_tools = list()
        self.makelist_tools = list()    
        self.makelist_options = {}   
        # self.bam2fqgz = self.add_tool('tool_lab.ysource.bam2fqgz')
        # self.makelist = self.add_tool('tool_lab.ysource.makelist')

    def check_options(self):
        if not self.option('bam_list').is_set:
            raise OptionError("没有导入bam文件list")
        with open(self.option('bam_list').prop['path'],'r') as bl:
            while 1:
                line = bl.readline()
                if not line:
                    break
                bam_path = line.rstrip()
                bam_name = re.findall(r"YGB\d{,}",bam_path)[0]
                self.logger.info(bam_name)
                self.logger.info(bam_path)
                self.bam_file[bam_name] = bam_path 
        for i in self.bam_file.keys():
            if not os.path.exists(self.bam_file[i]):
                raise OptionError("样本{}的bam文件没有找到".format(i))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(DatapreModule, self).run()
        self.run_tool()

    def run_tool(self):
        self.run_bam2fqgz()
        self.logger.info(len(self.bam2fqgz_tools))
        if len(self.bam2fqgz_tools) == 1 :
            self.bam2fqgz_tools[0].on('end', self.run_makelist)
        else:
            self.on_rely(self.bam2fqgz_tools, self.run_makelist)
        for tool in self.bam2fqgz_tools:
            tool.run()
        

    def run_bam2fqgz(self):
        self.bam_file_keys = sorted(self.bam_file.keys())
        for i in self.bam_file_keys:
            bam_path = self.bam_file[i]
            sample_name = i
            bam2fqgz = self.add_tool('tool_lab.ysource.bam2fqgz')
            options = {
                "bam_file": bam_path,
                "sample_name": sample_name
            }
            bam2fqgz.set_options(options)
            self.bam2fqgz_tools.append(bam2fqgz)
            self.makelist_options[sample_name] = {
                "sample_name" : sample_name ,
                "fastq1" : bam2fqgz.output_dir + "/{}_1.fq.gz".format(sample_name),
                "fastq2" : bam2fqgz.output_dir + "/{}_2.fq.gz".format(sample_name)
            }
    
    def run_makelist(self):
        
        for i in self.bam_file_keys:
            index = self.bam_file_keys.index(i)
            sample_name = i
            makelist = self.add_tool('tool_lab.ysource.makelist')
            # options = {
            #     "sample_name" : i ,
            #     "fastq1" : self.bam2fqgz_tools[index].output_dir + "/{}_1.fq.gz".format(i),
            #     "fastq2" : self.bam2fqgz_tools[index].output_dir + "/{}_2.fq.gz".format(i)
            # }
            options = self.makelist_options[sample_name]
            makelist.set_options(options)
            self.makelist_tools.append(makelist)
            self.logger.info(len(self.makelist_tools))
        if len(self.makelist_tools) == 1 :
            self.makelist_tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.makelist_tools,self.set_output)
        for tool in self.makelist_tools:
            tool.run()

    def set_output(self,event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        with open(os.path.join(self.output_dir, "sample_list"),'w') as spl:
            for tool in self.makelist_tools:
                for source in glob.glob(os.path.join(tool.output_dir, '*')):
                    link_name = os.path.join(self.output_dir,os.path.basename(source))
                    if os.path.exists(link_name):
                        os.remove(link_name)
                    with open(source,'r') as sor:
                        line = sor.readline()
                        spl.write(line)
                        spl.write("\n")
        self.end()

    def end(self):
        super(DatapreModule, self).end()


class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id':'datapre_' + str(random.randint(1,10000)),
            "type": "module",
            "name": "tool_lab.ysource.datapre",
            "options": {
                "bam_list":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/bam1.list"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()