# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
import random
import glob
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class YfullWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(YfullWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"sample_list","type":"infile","format":"small_rna.common"},
            # {"name":"data_type","type":"string"},
            {"name":"bam_num","type":"int"},
            {"name":"bam1","type":"infile","format":"align.bwa.bam"},
            {"name":"bam2","type":"infile","format":"align.bwa.bam"},
            {"name":"bam3","type":"infile","format":"align.bwa.bam"},
            {"name":"bam4","type":"infile","format":"align.bwa.bam"},
            {"name":"bam5","type":"infile","format":"align.bwa.bam"},
            # {"name":"bam6","type":"infile","format":"align.bwa.bam"},
            # {"name":"bam7","type":"infile","format":"align.bwa.bam"},
            {"name":"main_id","type":"string"},
            {"name":"update_info","type":"string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.autosome = True
        self.datapre = self.add_module("tool_lab.ysource.datapre")
        self.qc = self.add_module("tool_lab.ysource.qc")
        self.mapping = self.add_module("tool_lab.ysource.mapping")
        self.bam2vcf = self.add_module("tool_lab.ysource.bamtovcf")
        self.yfullsource = self.add_module("tool_lab.ysource.yfullsource")
        self.yfullresult = self.add_module("tool_lab.ysource.yfullresult")
        self.zipresult = self.add_tool("tool_lab.ysource.zipresult")
        self.out_api = ""
        self.qc_options = {}
        
    def check_options(self):
        """
        参数检查
        """
        if not self.option('sample_list') and not self.option('bam1'):
            raise OptionError("必须设置输入文件，sample_list或者bam文件")

    def run(self):
        self.logger.info('...options:{}'.format(self._options))
        if self.option("bam_num"):
            # self.autosome = False
            self.run_datapre()
            self.qc_options = {"sample_list":self.datapre.output_dir + "/sample_list"}
            self.datapre.on('end',self.run_qc)
            self.datapre.run()
        else:
            self.qc_options = {"sample_list":self.option('sample_list').prop['path']}
            self.run_qc()
        super(YfullWorkflow,self).run()

    def run_datapre(self):
        with open(self.work_dir + "/bam.list","w") as bl:
            for i in range(1,self.option("bam_num")+1):
                bl.write(self.option("bam{}".format(i)).prop['path'])
                bl.write("\n")
        options = {
            "bam_list": os.path.join(self.work_dir, 'bam.list')
        }
        self.datapre.set_options(options)
        # self.datapre.run()

    def run_qc(self):
        self.qc.set_options(self.qc_options)
        self.qc.on('end',self.run_mapping)
        self.qc.run()        

    def run_bam2vcf(self):
        options = {
            'bam_list' : self.mapping.output_dir + "/bam_list"
        }
        self.bam2vcf.set_options(options)
        self.bam2vcf.on('end',self.run_yfullsource)
        self.bam2vcf.run()

    def run_mapping(self):
        options = {
            "bam_list" : self.qc.output_dir + "/bam_list"
        }
        self.mapping.set_options(options)
        self.mapping.on('end',self.run_bam2vcf)
        self.mapping.run()

    def run_yfullsource(self):
        options = {
            "vcf_list": self.bam2vcf.output_dir + "/vcf_list",
            "bam_path": self.mapping.output_dir,
            # "pos_list": self,
            # "tree" : self,
            "yfull_list": self.bam2vcf.output_dir + "/yfull_list",
            "autosome": self.autosome
        }
        self.yfullsource.set_options(options)
        self.yfullsource.on('end',self.run_yfullresult)
        self.yfullsource.run()

    def run_yfullresult(self):
        options = {
            "yfull_list":self.bam2vcf.output_dir+"/yfull_list",
            "private_list":self.bam2vcf.output_dir + "/private_list",
            "source_list":self.yfullsource.output_dir + "/yfullsource_list",
            "all_ot":self.yfullsource.output_dir + "/all_sample_MAP_OT.xls",
            "all_dep":self.yfullsource.output_dir + "/all_sample_dep.xls",
            "size_file":self.mapping.output_dir + "/size.xls"
        }
        self.yfullresult.set_options(options)
        self.yfullresult.on('end',self.run_zipresult)
        self.yfullresult.run()

    def run_zipresult(self):
        options = {
            "result_file": self.yfullresult.output_dir
        }
        self.zipresult.set_options(options)
        self.zipresult.on('end',self.set_output)
        self.zipresult.run()

    def set_output(self):
        self.logger.info('strat set_output as {}'.format(
            self.__class__.__name__))
        zip_output = self.zipresult.output_dir + "/result.zip"
        yfull_output=self.output_dir + "/result.zip"
        if os.path.exists(yfull_output):
            os.remove(yfull_output)
        os.link(zip_output, yfull_output)
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        for root, dirs, files in os.walk(self.yfullresult.output_dir, topdown=False):
                self.out_api = self.out_api + "\t".join(files)
        api_yfull = self.api.api("tool_lab.ysource.yfull")
        api_yfull.add_yfull_detail(self.option('main_id'),self.out_api)
        self.logger.info("导表结束")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(YfullWorkflow, self).end()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test_yfull',
        'id': 'tsg_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
            "sample_list": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/sample_list",
            "main_id" : "5e9e6a6017b2bf2049a81be3"
        }
    }
    wsheet = Sheet(data=data)
    wf = YfullWorkflow(wsheet)
    wf.run()