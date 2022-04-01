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

class BamtovcfModule(Module):
    def __init__(self, work_id):
        super(BamtovcfModule, self).__init__(work_id)
        options = [
            {"name":"bam_list","type":"string"},

        ]
        self.add_option(options)
        self.sp_info = {}
        self.sp_list = []
        self.bam2vcf_tools = []
        self.yfull_bed = "database/Tool_lab/ysource/yfull_hg38_20200801.bed"
        self.Y_8M_bed = "database/Tool_lab/ysource/8M_hg38.bed"
        self.pos_db_bed = "database/Tool_lab/ysource/posdb_update_20200801.bed"
    
    def check_options(self):
        if not os.path.exists(self.option("bam_list")):
            raise OptionError("没有找到bam_list文件")
        with open(self.option("bam_list"),'r') as bl:
            while 1:
                line = bl.readline()
                self.logger.info(line)
                if not line:
                    break
                field = line.rstrip().split('\t')
                sn = field[0]
                bam_path = field[1]
                self.sp_list.append(sn)
                if not os.path.exists(bam_path):
                    raise OptionError("缺少bam文件")
                self.sp_info[sn] = bam_path
        return True
    
    def run(self):
        super(BamtovcfModule, self).run()
        self.run_tool()

    def run_tool(self):
        for i in self.sp_list:
            self.run_bam2vcf(self.pos_db_bed,self.sp_info[i],i)
            self.run_bam2vcf(self.yfull_bed,self.sp_info[i],i)
            self.run_bam2vcf(self.Y_8M_bed,self.sp_info[i],i)
        self.on_rely(self.bam2vcf_tools,self.set_output)
        for tool in self.bam2vcf_tools:
            tool.run()

    def run_bam2vcf(self, bed, bam, sn):
        output = ""
        if bed == self.pos_db_bed:
            output = "{}_align_final_sort".format(sn)
        elif bed == self.Y_8M_bed:
            output = "{}_align_final_sort_private".format(sn)
        elif bed == self.yfull_bed:
            output = "{}_align_final_sort_yfull".format(sn)
        bam2vcf = self.add_tool('tool_lab.ysource.bamtovcf')
        bam2vcf.set_options({
            "bam_file" : bam,
            "output_name" : output,
            "bed" : bed
        })
        self.bam2vcf_tools.append(bam2vcf)

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(
            self.__class__.__name__))
        with open(os.path.join(self.output_dir, "vcf_list"), 'w') as bl,open(os.path.join(self.output_dir, "yfull_list"), 'w') as yf, open(os.path.join(self.output_dir, "private_list"), 'w') as pl:
            for tool in self.bam2vcf_tools:
                # sample_name = tool.option('output_name')
                
                for source in glob.glob(os.path.join(tool.output_dir, '*')):
                    link_name = os.path.join(
                        self.output_dir, os.path.basename(source))
                    if os.path.exists(link_name):
                        os.remove(link_name)
                    os.link(source, link_name)
                    self.logger.info(
                        'succeed in linking {} to {}'.format(source, link_name))
                text1 = ""
                if os.path.basename(source)[-21:] == "_align_final_sort.vcf":
                    sample_name = os.path.basename(source)[:-21]
                    text1 = link_name
                    text = '\t' + text1
                    bl.write(sample_name)
                    bl.write(text)
                    bl.write('\n')
                if os.path.basename(source)[-29:] == "_align_final_sort_private.vcf":
                    sample_name = os.path.basename(source)[:-29]
                    text1 = link_name
                    text = '\t' + text1
                    pl.write(sample_name)
                    pl.write(text)
                    pl.write('\n')
                if os.path.basename(source)[-27:] == "_align_final_sort_yfull.vcf":
                    sample_name = os.path.basename(source)[:-27]
                    text1 = link_name
                    text = '\t' + text1
                    yf.write(sample_name)
                    yf.write(text)
                    yf.write('\n')
        self.logger.info('finish set_output at {}'.format(
            self.__class__.__name__))
        self.end()

    def end(self):
        super(BamtovcfModule, self).end()

class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bamtovcf_' + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.ysource.bamtovcf",
            "options": {
                "bam_list": "/mnt/ilustre/users/sanger-dev/workspace/20201226/Single_Mapping_6220/Mapping/output/bam_list"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()