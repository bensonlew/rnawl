# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from collections import namedtuple, defaultdict
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import random
import datetime
import xlsxwriter
import subprocess
import re
import os
import sys
import shutil

class YfullsourceAgent(Agent):
    def __init__(self, parent):
        super(YfullsourceAgent, self).__init__(parent)
        options = [
            {"name":"vcf","type":"string"},
            {"name":"sample_name","type":"string"},
            # {"name":"pos_list","type":"string"},
            # {"name":"tree","type":"string"}
        ]
        self.add_option(options)
    
    def check_options(self):
        """
        参数检查
        """
        if not self.option("vcf"):
            raise OptionError("导入VCF文件失败")
        if not self.option("sample_name"):
            raise OptionError("没有输入样本名")
    
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"
    
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(YfullsourceAgent, self).end()

class YfullsourceTool(Tool):
    def __init__(self, config):
        super(YfullsourceTool, self).__init__(config)

        # self.pos_list = self.option('pos_list')
        # self.tree = self.option('tree')
        self.vcf = self.option("vcf")
        self.sample_name = self.option("sample_name")
        self._version = '1.0'
        self.software_dir = self.config.SOFTWARE_DIR
        self.pos = os.path.join(self.software_dir,"database/Tool_lab/ysource/pos")
        self.tree = os.path.join(self.software_dir,"database/Tool_lab/ysource/tree")
        self.perl = 'program/perl-5.24.0/bin/perl'
        self.script = self.config.PACKAGE_DIR + '/tool_lab/yoogene/yfull/sourse_Y_20201105.pl'

    def run(self):
        """
        运行
        """ 
        super(YfullsourceTool, self).run()
        self.run_yfullsource()
        self.set_output()
        self.end()

    def run_txtToExcel(self,input, output):
        if os.path.exists(input):
            f = open(input)
            wo = xlsxwriter.Workbook(output)
            sheet = wo.add_worksheet()
            x = 0
            while 1:
                line = f.readline()
                if not line:
                    break
                for i in range(len(line.split('\t'))):
                    item = line.split('\t')[i]
                    sheet.write(x, i, item)
                x += 1
            wo.close()
            f.close()
    
    def run_yfullsource(self):
        '''
        perl source_Y.pl -sample_name sample_name -sample_path sample_path
            -pos pos -tree tree -type Y_full
        '''
        cmd = "{} {} ".format(self.perl, self.script)
        cmd = cmd + "-sample_name {} -sample_path {} ".format(self.sample_name,self.vcf)
        cmd = cmd + "-pos {} -tree {} -type Y_full ".format(self.pos, self.tree)
        self.logger.info(cmd)
        self.logger.info("开始计算Yfull")
        command1 = self.add_command("yfullsource_{}".format(self.sample_name.lower()), cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行Yfullsource完成")
        else:
            self.logger.info("运行错误，重新检查参数")
        asr_in = os.path.join(self.work_dir, "1/all_sample_result.txt")
        asr_on = os.path.join(self.work_dir, "1/all_sample_result.xlsx")
        self.run_txtToExcel(asr_in,asr_on)
    
    def set_output(self):
        '''
        将结果文件赋值到output文件夹下面
        :return:
        '''
        if len(self.output_dir) >0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        self.logger.info("设置结果目录")
        try:
            shutil.copytree(os.path.join(self.work_dir, "1"),
                    os.path.join(self.output_dir, "1"))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))
    

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "yfullsource_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.yfullsource",
            "options": dict(
                vcf="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_script/align_final_sort.vcf",
                sample_name="YG202003761",
                # pos_list = "/mnt/ilustre/users/sanger-dev/app/database/Tool_lab/ysource/Y_20200801_list_fix.txt",                
                # tree="/mnt/ilustre/users/sanger-dev/app/database/Tool_lab/ysource/20200801_tree_fix.txt",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

