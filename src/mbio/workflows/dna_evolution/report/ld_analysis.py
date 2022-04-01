# -*- coding: utf-8 -*-
# __author__ = 'Binbin Zhao'
# modified 20181214

import os
import json
import re
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class LdAnalysisWorkflow(Workflow):
    """
    交互分析：连锁不平衡workflow
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LdAnalysisWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_file", "type": 'infile', "format": "dna_gmap.vcf"},
            {"name": "min_dp", "type": "string"},
            {"name": "max_dp", "type": "string"},
            {"name": "max_missing", "type": "string", "default": "0.3"},
            {"name": "min_maf", "type": "string", "default": "0.05"},
            {"name": "max_maf", "type": "string", "default": "1"},
            {"name": "group_dict", "type": "string"},
            {"name": "recode", "type": "bool"},
            {"name": "task_id", "type": "string"},
            {"name": "project_sn", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ld_analysis = self.add_module("dna_evolution.ld_analysis")

    def check_options(self):
        a = vars(self.option("vcf_file"))
        if not self.option("vcf_file"):
            raise OptionError("请输入vcf.file文件")
        if not self.option("max_missing"):
            raise OptionError("请输入缺失率")
        else:
            if float(self.option("max_missing")) < 0 or float(self.option("max_missing")) > 1:
                raise OptionError("缺失率应该为0到1之间的数字")
        if not self.option("min_maf"):
            raise OptionError("请输入次要等位基因最小值")
        else:
            if float(self.option("min_maf")) > 1:
                raise OptionError("次要等位基因频率最大值为1")
        if not self.option("max_maf"):
            raise OptionError("请输入次要等位基因最大值")
        else:
            if float(self.option("max_maf")) > 1:
                raise OptionError("次要等位基因频率最大值为1")
        return True

    def run_ld_analysis(self):
        options = {
            "vcf_file": self.option("vcf_file"),
            "min_dp": self.option("min_dp"),
            "max_dp": self.option("max_dp"),
            "max_missing": self.option("max_missing"),
            "min_maf": self.option("min_maf"),
            "max_maf": self.option("max_maf"),
            "recode": self.option("recode"),
            "task_id": self.option("task_id"),
            "project_sn": self.option("project_sn"),
            "main_id": self.option("main_id"),
            "update_info": self.option("update_info"),
            "graph_path":  self._sheet.output,
            "group_dict": self.option("group_dict")
        }
        self.ld_analysis.set_options(options)
        self.ld_analysis.on("end", self.set_output, "ld_analysis")
        self.ld_analysis.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'ld_analysis':
            self.linkdir(obj.output_dir, 'ld_analysis')
            self.end()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        self.run_ld_analysis()
        super(LdAnalysisWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(LdAnalysisWorkflow, self).end()
