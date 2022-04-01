# -*- coding: utf-8 -*-
# __author__ = 'konghualei' 20170420

"""Venn表计算"""

import os
import datetime
import json
import shutil
import re
from biocluster.workflow import Workflow
import pandas as pd

class ExpressVennWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpressVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "express_file", "type":"string"},
            {"name": "group_id", "type":"string"},  #样本的分组信息
            {"name":"group_detail","type":"string"},
            {"name": "update_info", "type": "string"},
            {"name":"type","type":"string"},#对应gene/transcript
            {"name":"express_level","type":"string"}, #对应fpkm/tpm
            {"name":"threshold","type":"float","default":1}, #过滤
            # {"name":"sample_group",'type':"string","default":"sample"},
            {"name":"venn_id","type":"string"},
        ]
        self.logger.info(options)
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.venn = self.add_tool("rna.expressvenn")
        with open(self.option('group_id'), 'r+') as f1:
            f1.readline()
            if not f1.readline():
                self.group_id = 'all'
            else:
                self.group_id = self.option('group_id')
        # self.samples = re.split(',', self.option("specimen"))
        
    def run_venn(self, fpkm_path, specimen):
        """样本间特异性基因venn图"""
        options = {
            "express_matrix": fpkm_path,
            "group_table": self.option("group_id"),
            "threshold":self.option("threshold")
        }
        
        self.logger.info("检查new_fpkm和new_group的路径:")
        self.venn.set_options(options)
        self.venn.on('end', self.set_db)
        self.venn.run()
        
    def set_db(self):
        venn_path = self.venn.output_dir
        # 更新venn_graph.xls文件
        shutil.copy2(self.venn.work_dir+"/new_venn_graph.xls",venn_path+"/venn_graph.xls")
        api_venn = self.api.refrna_corr_express
        venn_id = self.option("venn_id")
        self.logger.info("准备开始向mongo数据库中导入venn图detail表和graph信息！")
        api_venn.add_venn_detail(venn_path + '/venn_table.xls', venn_id, 'ref')
        api_venn.add_venn_graph(venn_path + '/venn_graph.xls', venn_id, 'ref')
        self.logger.info("导入venn图detail表和graph表成功！")
        self.end()
    
    def get_samples(self):
        edger_group_path = self.option("group_id")
        self.logger.info(edger_group_path)
        self.samples=[]
        with open(edger_group_path,'r+') as f1:
            f1.readline()
            for lines in f1:
                line=lines.strip().split("\t")
                self.samples.append(line[0])
        print self.samples
        return self.samples

    def run(self):
        fpkm = self.option("express_file").split(",")[0]
        if self.group_id in ['all','All','ALL']:
            with open(fpkm,'r+') as f1:
                specimen = f1.readline().strip().split("\t")
        else:
            specimen = self.get_samples()
        self.run_venn(fpkm, specimen)
        super(ExpressVennWorkflow, self).run()
    
    def end(self):
        output1_dir = self.venn.output_dir
        result = self.add_upload_dir(output1_dir)
        result.add_relpath_rules([[".", "", "表达量Venn分析结果文件"], ])
        super(ExpressVennWorkflow, self).end()


