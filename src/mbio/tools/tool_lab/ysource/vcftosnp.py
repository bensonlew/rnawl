# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from mbio.packages.whole_transcriptome.utils import runcmd
# import pandas as pd
import unittest
import xlsxwriter
import datetime
import random
import re
import os
import sys
import shutil

class VcftosnpAgent(Agent):
    def __init__(self, parent):
        super(VcftosnpAgent, self).__init__(parent)
        self.options = [
            {"name":"yfull_vcf","type":"string"},
            # {"name":"yfull_xls","type":"string"},
            {"name":"sample_txt","type":"string"},
            {"name":"private_vcf","type":"string"},
            # {"name":"snp","type":"string"},
            # {"name":"tree","type":"string"},
            # {"name":"pos","type":"string"},
            {"name":"sample_name","type":"string"}
        ]
        self.add_option(self.options)

    def check_option(self):
        """
        参数检查
        """
        for i in self.options:
            op = self.option(i['name'])
            if not os.path.exists(op):
                raise OptionError("找不到{}".format(op))

    def set_resource(self):
        self._cpu = 2 
        self._memory = "10G"

    def end(self):
        super(VcftosnpAgent, self).end()

class VcftosnpTool(Tool):
    """
    version 1.0
    """
    def __init__(self,config):
        super(VcftosnpTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.pos = os.path.join(self.software_dir,"database/Tool_lab/ysource/pos")
        self.tree = os.path.join(self.software_dir,"database/Tool_lab/ysource/tree")
        self.snp = os.path.join(self.software_dir,"database/Tool_lab/ysource/SOGG_snp.xls")
        # self.snp = self.option("snp")
        self.perl = 'program/perl-5.24.0/bin/perl'
        self.script_dir = self.config.PACKAGE_DIR + '/tool_lab/yoogene/yfull'

    def run(self):
        super(VcftosnpTool, self).run()
        self.run_vcf2snp_8M()
        self.run_vcf2snp_new()
        self.run_vcf2snp_final()
        self.set_output()
        self.end()

    def run_vcf2snp_final(self):
        """
        perl 2_vcfToSNP_final_20200801fix.pl -in YFULL_SNP_reslt.xls 
            -m sample.txt - tree tree - pos pos
        """
        cmd = "{} {}/2_vcfToSNP_final_20200801fix.pl ".format(self.perl, self.script_dir)
        cmd = cmd +"-in YFULL_SNP_reslt.xls -m {} ".format(self.option("sample_txt"))
        cmd = cmd + "-tree {} -pos {}".format(self.tree,self.pos)
        self.logger.info(cmd)
        command = self.add_command("vcf2snp_final_{}".format(self.option("sample_name")).lower(),cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}vcf2snp_final运行完成".format(self.option("sample_name")))
        else:
            self.logger.info("{}vcf2snp_final运行失败".format(self.option("sample_name")))
        out_file = os.path.join(self.work_dir,"final_list_result.txt")
        self.del_yfull(out_file,"{}_9000_SNP.xlsx".format(self.option("sample_name")))
        
            
    def run_vcf2snp_new(self):
        """
        perl 2_vcfToSNP_new_20200801fix.pl -vcf yfull_sample.vcf 
        -tree tree -pos pos_list -SNP
        """
        cmd = "{} {}/2_vcfToSNP_new_20200801fix.pl ".format(self.perl, self.script_dir)
        cmd = cmd +"-vcf {} ".format(self.option("yfull_vcf"))
        cmd = cmd + "-tree {} -pos {} -SNP {}".format(self.tree,self.pos,self.snp)
        self.logger.info(cmd)
        command = self.add_command("vcf2snp_new_{}".format(self.option("sample_name")).lower(),cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}vcf2snp_new运行完成".format(self.option("sample_name")))
        else:
            self.logger.info("{}vcf2snp_new运行失败".format(self.option("sample_name")))
        out_file = os.path.join(self.work_dir, "YFULL_SNP_reslt.xls")
        self.run_txtToExcel(out_file,"{}_YFULL_SNP.xlsx".format(self.option("sample_name")))


    def run_vcf2snp_8M(self):
        """
        perl 2_vcfToSNP_8M_private.pl -vcf private_vcf 
            -SNP SOGG
        """
        cmd = "{} {}/2_vcfToSNP_8M_private.pl ".format(self.perl, self.script_dir)
        cmd = cmd +"-vcf {} ".format(self.option("private_vcf"))
        cmd = cmd + "-SNP {}".format(self.snp)
        self.logger.info(cmd)
        command = self.add_command("vcf2snp_8m_{}".format(self.option("sample_name")).lower(),cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}vcf2snp_new运行完成".format(self.option("sample_name")))
        else:
            self.logger.info("{}vcf2snp_new运行失败".format(self.option("sample_name")))
        out_file1 = os.path.join(self.work_dir, "sample_SNP.txt")
        out_file2 = os.path.join(self.work_dir, "sample_list.txt")
        self.run_txtToExcel(out_file1, "{}-private.xlsx".format(self.option("sample_name")))
        self.run_txtToExcel(out_file2, "{}-private-list.xlsx".format(self.option("sample_name")))
    
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

    def del_yfull(self, input, output):
        if os.path.exists(input):
            f = open(input)
            wo = xlsxwriter.Workbook(output)
            sheet = wo.add_worksheet()
            x = 0
            while 1:
                line = f.readline()
                line_list = line.split('\t')
                if not line:
                    break
                if len(line_list) >= 3:
                    if line_list[2] != '':
                        for  i in range(len(line.split('\t'))):
                            item = line.split('\t')[i]
                            sheet.write(x, i, item)
                        x += 1
            wo.close()
            f.close()
                        


    def set_output(self):
        """
        设置输出文件夹
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir, "{}_YFULL_SNP.xlsx".format(self.option("sample_name"))),
                    os.path.join(self.output_dir, "{}_YFULL_SNP.xlsx".format(self.option("sample_name"))))
            os.link(os.path.join(self.work_dir, "final_list_result.txt"),
                    os.path.join(self.output_dir, "{}_9000_SNP.xlsx".format(self.option("sample_name"))))
            os.link(os.path.join(self.work_dir, "{}-private.xlsx".format(self.option("sample_name"))),
                    os.path.join(self.output_dir, "{}-private.xlsx".format(self.option("sample_name"))))
            os.link(os.path.join(self.work_dir, "{}-private-list.xlsx".format(self.option("sample_name"))),
                    os.path.join(self.output_dir, "{}-private-list.xlsx".format(self.option("sample_name"))))        
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


