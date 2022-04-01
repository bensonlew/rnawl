# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
import json
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class VcfLddecayWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VcfLddecayWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "group_file","type": "infile","format":"denovo_rna_v2.common"},
            {"name": "max_missing", "type": "float", "default": 30.0},  # 缺失率
            {"name": "maf", "type": "float", "default": 0.05},  # 次要等位基因频率min
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        # self.filter = self.add_tool("tool_lab.vcf_filter")
        # self.plink = self.add_tool("tool_lab.vcf_plink")
        self.lddecay = self.add_module('tool_lab.ld_analysis')

    def check_option(self):
        """
        参数检查
        """
        if not self.option("vcf_file"):
            raise OptionError("必须输入vcf文件")
        if not self.option("group_file"):
            raise OptionError("必须输入分组文件")
        if not self.option("max_missing"):
            raise OptionError("必须设置最大缺失率")
        if not self.option("maf"):
            raise OptionError("必须设置maf")
    
    # def run_filter(self):
    #     self.filter.set_options({
    #         'vcf_file':self.option("vcf_file"),
    #         'max_missing':self.option("max_missing"),
    #         "maf":self.option('maf')
    #     })
    #     self.filter.on('end',self.run_plink)
    #     self.filter.run()

    # def run_plink(self):
    #     self.filter_vcf = os.path.join(self.filter.output_dir,"pop.recode.vcf")
    #     self.plink.set_options({
    #         'vcf_file': self.filter_vcf,
    #     })
    #     self.plink.on('end',self.run_lddecay)
    #     self.plink.run()

    def run_lddecay(self):
        group_file = self.option("group_file").prop["path"]
        group_dict = {}
        with open(group_file,"r") as gf:
            while 1:
                line = gf.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                if group_dict.has_key(fd[1]):
                    group_dict[fd[1]] = group_dict[fd[1]] +","+fd[0]
                else:
                    group_dict[fd[1]] = fd[0]
        self.lddecay.set_options({
            'vcf_file':self.option("vcf_file"),
            'max_missing':str(self.option('max_missing')/100),
            "min_maf":str(self.option("maf")),
            "group_dict":json.dumps(group_dict)
        })
        self.lddecay.on('end',self.set_output)
        self.lddecay.run()
    
    def set_output(self):
        self.logger.info('strat set_output as {}'.format(
            self.__class__.__name__))
        decay_path = os.path.join(self.lddecay.output_dir,"ld_decay_dir")
        draw_path = os.path.join(self.lddecay.output_dir,"ld_draw_dir")
        try:
            shutil.copytree(decay_path,
                                os.path.join(self.output_dir, "ld_decay_dir"))
            shutil.copytree(draw_path,
                                os.path.join(self.output_dir, "ld_draw_dir"))
        except Exception as e:
            self.set_error("设置结果目录失败{}".format(e))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        pictrue_path = self._sheet.output + "/ld_draw_dir/ld.pdf"
        self.logger.info(pictrue_path)
        api_vcf_lddecay = self.api.api("tool_lab.vcf_lddecay")
        api_vcf_lddecay.add_detail(self.option("main_id"),pictrue_path)
        api_vcf_lddecay.vcf_ld_detail(self.option("main_id"),os.path.join(self.lddecay.work_dir,"gro_list"))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        self.run_lddecay()
        super(VcfLddecayWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VcfLddecayWorkflow, self).end()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    data = {
        'name': 'test_vcf_lddecay',
        'id': 'vcf_lddecay_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "vcf_file": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/test_WGS/6.ld-decay/pop.recode.vcf",
        "group_file":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/test_WGS/6.ld-decay/group.list",
        "main_id" : "5e9e6a6017b2bf2049a81be2"
        }
    }
    wsheet = Sheet(data=data)
    wf = VcfLddecayWorkflow(wsheet)
    wf.run()