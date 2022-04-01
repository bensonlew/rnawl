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


class YfullsourceModule(Module):
    def __init__(self, work_id):
        super(YfullsourceModule, self).__init__(work_id)
        options = [
            {"name": "vcf_list", "type": "string"},
            # {"name":"pos_list","type":"string"},
            # {"name":"tree","type":"string"},
            {"name": "bam_path", "type": "string"},
            {"name": "autosome", "type": "bool", "default": True},
            {"name": "yfull_list", "type": "string"},
            # {"name":"private_list","type":"string"}
        ]
        self.add_option(options)
        self.sp_info = {}
        self.sp_list = []
        self.yfullsource_tools = list()
        self.adddatabase = self.add_tool("tool_lab.ysource.add_database")
        self.tools = []
        self.mapstat_tools = []
        self.yfull_info = {}
        # self.autocount = self.add_tool("tool_lab.ysource.autocount")

    def check_options(self):
        if not os.path.exists(self.option("vcf_list")):
            raise OptionError("找不到vcflist，请检查一下")
        with open(self.option("vcf_list"), 'r') as vl:
            while 1:
                line = vl.readline()
                if not line:
                    break
                field = line.rstrip().split('\t')
                sn = field[0]
                vcf = field[1]
                self.sp_list.append(sn)
                if not os.path.exists(vcf):
                    raise OptionError("没有vcf文件")
                self.sp_info[sn] = vcf
        with open(self.option("yfull_list"), 'r') as vf:
            while 1:
                line = vf.readline()
                if not line:
                    break
                field = line.rstrip().split('\t')
                sn = field[0]
                yfvcf = field[1]
                # self.sp_list.append(sn)
                if not os.path.exists(yfvcf):
                    raise OptionError("没有yfull_vcf文件")
                self.yfull_info[sn] = yfvcf
        return True

    def run(self):
        super(YfullsourceModule, self).run()
        self.run_tools()

    def run_tools(self):
        if self.option("autosome"):
            self.run_autocount()

        self.run_yfullsource()
        self.tools.extend(self.yfullsource_tools)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.run_add_database)
        else:
            self.on_rely(self.tools, self.run_add_database)
        # if len(self.yfullsource_tools) == 1:
        #     self.yfullsource_tools[0].on('end', self.run_add_database)
        # else:
        #     self.on_rely(self.yfullsource_tools, self.run_add_database)
        # for tool in self.yfullsource_tools:
        #     tool.run()
        for tool in self.tools:
            tool.run()

    def run_yfullsource(self):
        for i in self.sp_list:
            sn = i
            vcf = self.sp_info[sn]
            yfullsource = self.add_tool('tool_lab.ysource.yfullsource')
            options = {
                "sample_name": sn,
                "vcf": vcf,
                # "pos_list" : self.option('pos_list'),
                # "tree" : self.option('tree')
            }
            yfullsource.set_options(options)
            self.yfullsource_tools.append(yfullsource)

    def run_autocount(self):
        autocount = self.add_tool("tool_lab.ysource.autocount")
        options = {
            'vcf_list': self.option('vcf_list')
        }
        autocount.set_options(options)
        self.tools.append(autocount)
        # self.autocount.on('end', self.set_output)
        # self.autocount.run()

    def run_add_database(self):
        text = ""
        for tool in self.yfullsource_tools:
            asr_file = os.path.join(tool.output_dir, "1/all_sample_result.txt")
            with open(asr_file, "r") as sr:
                text = text + sr.readline()
        all_result = os.path.join(self.work_dir, "all_sample_result.txt")
        with open(all_result, "w") as aw:
            aw.write(text)
        options = {
            "source_result": all_result
        }
        self.adddatabase.set_options(options)
        self.adddatabase.on('end', self.run_mapstat)
        self.adddatabase.run()

    def run_mapstat(self):
        bam_path = self.option('bam_path')
        for i in self.sp_info.keys():
            sample_name = i
            bam_file = os.path.join(
                bam_path, "{}_align_final_sort.bam".format(sample_name))
            bam_pre_file = os.path.join(
                bam_path, "{}_aln-pe.sort.bam".format(sample_name))
            vcf_yfull = self.yfull_info[sample_name]
            clean_statistic = os.path.join(
                bam_path, "{}_clean_statistic_reform.xls".format(sample_name))
            mapstat = self.add_tool("tool_lab.ysource.mapstat")
            options = {
                "bam_file": bam_file,
                "bam_prermdup": bam_pre_file,
                "vcf_yfull": vcf_yfull,
                "clean_statistic": clean_statistic,
                "sample_name": sample_name
            }
            mapstat.set_options(options)
            self.mapstat_tools.append(mapstat)
        if len(self.mapstat_tools) == 1:
            self.mapstat_tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.mapstat_tools, self.set_output)
        for tool in self.mapstat_tools:
            tool.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(
            self.__class__.__name__))
        for tool in self.yfullsource_tools:
            out_dir = os.path.join(tool.output_dir, '1')
            sample_name = tool.option('sample_name')
            try:
                shutil.copytree(out_dir,
                                os.path.join(self.output_dir, "{}_1".format(sample_name)))
            except Exception as e:
                self.logger.info("设置结果目录失败{}".format(e))
            with open(os.path.join(self.output_dir, "yfullsource_list"), 'a') as yl:
                text1 = ""
                if os.path.join(self.output_dir, "{}_1".format(sample_name)):
                    text1 = os.path.join(
                        self.output_dir, "{}_1".format(sample_name))
                    text = '\t' + text1
                    yl.write(sample_name)
                    yl.write(text)
                    yl.write('\n')
        with open(os.path.join(self.output_dir, "all_sample_MAP_OT.xls"), "w") as ot, open(os.path.join(self.output_dir, "all_sample_dep.xls"), "w") as ad:
            ot.write(
                "sample\tclean_bp\tmapping_bp\tontarget_bp\tmapping_rate\tOT\tdedup_ratio\n")
            ad.write("sample\tpos\t1X\t2X\t5X\t7X\t10X\t15X\t20X\n")
            for tool in self.mapstat_tools:
                ot_file = os.path.join(
                    tool.output_dir, "{}_sample_MAP_OT.xls".format(tool.option('sample_name')))
                dep_file = os.path.join(
                    tool.output_dir, "{}_sample_dep.xls".format(tool.option('sample_name')))
                with open(ot_file, "r") as sot, open(dep_file, "r") as sdep:
                    sot.readline()
                    sdep.readline()
                    ot.write(sot.readline())
                    ad.write(sdep.readline())
        self.end()

    def end(self):
        super(YfullsourceModule, self).end()


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
            "type": "module",
            "name": "tool_lab.ysource.yfullsource",
            "options": dict(
                vcf_list="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_bamtovcf_6434/Bamtovcf/output/vcf_list",
                bam_path="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_Mapping_1564/Mapping/output",
                yfull_list="/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_bamtovcf_6434/Bamtovcf/output/yfull_list"
                # pos_list = "/mnt/ilustre/users/sanger-dev/app/database/Tool_lab/ysource/Y_20200801_list_fix.txt",
                # tree="/mnt/ilustre/users/sanger-dev/app/database/Tool_lab/ysource/20200801_tree_fix.txt",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
