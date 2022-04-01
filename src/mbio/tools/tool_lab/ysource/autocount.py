# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import datetime
import random
import re
import os
import sys
import shutil

class AutocountAgent(Agent):
    """
    统计Vcf常染色体上的SNP位点
    小工具接口
    """
    def __init__(self, parent):
        super(AutocountAgent, self).__init__(parent)
        options =[
            {"name":"vcf_list","type":"string"},
        ]
        self.add_option(options)

    def check_option(self):
        """
        参数检查
        """
        if not self.option("vcf_list"):
            raise OptionError("没有传入VCFlist")
        with open(self.option("vcf_list"),"r") as vl:
            while 1:
                line = vl.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                if not os.path.exists(fd[1]) or not os.path.exists(fd[2]):
                    raise OptionError("vcf列表有误")
        
    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(AutocountAgent, self).end()

class AutocountTool(Tool):
    def __init__(self, config):
        super(AutocountTool,self).__init__(config)
        self._version = "v1.0"
        self.vcf_path = {}
        self.software_dir = self.config.SOFTWARE_DIR
        with open(self.option("vcf_list"),"r") as vl:
            while 1:
                line = vl.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                self.vcf_path[fd[0]] = fd[1]
        self.script_autoresult = self.config.PACKAGE_DIR + '/tool_lab/yoogene/yfull/autosome_result.pl'
        self.script_combineresult = self.config.PACKAGE_DIR + '/tool_lab/yoogene/yfull/combine_result.py'
        self.snp_info = os.path.join(self.software_dir,"database/Tool_lab/ysource/Y_20200801_list_fix.txt")
        self.perl = "program/perl-5.24.0/bin/perl"
        self.python = 'program/Python/bin/python'

    def run(self):
        '''
        运行
        '''  
        super(AutocountTool, self).run()
        for sn in self.vcf_path.keys():
            self.run_autosome_result(sn, self.vcf_path[sn])
            self.run_combine_result(sn)
        self.set_output()
        self.end()

    def run_autosome_result(self, sn, vcf):
        """
        perl autosome_result.pl -vcf vcf_file 
            -sample_name sample_name -snp_file snp_file
        """
        cmd = "{} {} ".format(self.perl, self.script_autoresult)
        cmd += "-vcf {} -sample_name {} ".format(vcf, sn)
        cmd += "-snp_file {} ".format(self.snp_info)
        self.logger.info(cmd)
        self.logger.info("{} 开始生成常染snp表".format(sn))
        command1 = self.add_command("mkautosome_result_{}".format(sn.lower()),cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("{}常染snp表生成".format(sn))
        else:
            self.logger.info("运行错误")
    
    def run_combine_result(self, sn):
        """
        python combine_result.py -i inputfile -n sn
        """
        cmd = "{} {} ".format(self.python, self.script_combineresult)
        infile = self.work_dir + "/{}/{}all_antosome_snp_result.txt".format(sn,sn)
        cmd += "-i {} -n {}".format(infile, sn)
        self.logger.info(cmd)
        self.logger.info("开始运行{}_combine_result.py".format(sn))
        command1 = self.add_command("{}_antosome_snp_result".format(sn.lower()),cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("{}_antosome_snp_result.txt生成".format(sn))
        else:
            self.logger.info("运行错误")
    
    def set_output(self):
        '''
        将结果文件赋值到output文件夹下面
        :return:
        '''
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        self.logger.info("设置结果目录")
        for i in self.vcf_path.keys():
            try:
                os.link(os.path.join(self.work_dir, "{}_antosome_snp_result.txt".format(i)),
                        os.path.join(self.output_dir, "{}_antosome_snp_result.txt".format(i)))
            except Exception as e:
                self.logger.info("设置结果目录失败{}".format(e))

    

    
    

    

        

            
