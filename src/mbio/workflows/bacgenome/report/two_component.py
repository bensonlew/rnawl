# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifies 20190401

'''双组分调控系统'''

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from biocluster.file import download, exists



class TwoComponentWorkflow(Workflow):
    """
    报告中使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(TwoComponentWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pfam_anno", "type": "string"},   # 基因具体的注释信息表，‘，分割’
            {'name': 'sample_name', "type": "string",'default':'out'} ,  #样本名，‘，分割’
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

        self.tool_list = []
        self.samples = self.option("sample_name").split(',')
        self.pfam_files = self.option("pfam_anno").split(",")
        self.num = len(self.samples)

    def download_anno_files(self):
        for i in range(self.num):
            target_file = self.work_dir+ '/'+ self.samples[i]+'.pfam'
            ori_file = self.pfam_files[i]
            self.logger.info("download:")
            self.logger.info(ori_file)
            self.logger.info(target_file)
            download(ori_file, target_file)

    def run_two_component(self):

        for i in range(self.num):
            options = {
                'pfam_anno': self.work_dir+ '/'+ self.samples[i]+'.pfam',
                'sample_name':  self.samples[i]
            }
            two_component_tool = self.add_tool("annotation.two_component")
            two_component_tool.set_options(options)
            self.tool_list.append(two_component_tool)

        if len(self.tool_list) == 1:
            self.tool_list[0].on('end', self.set_db)
        else:
            self.on_rely(self.tool_list,self.set_db)
        for t in self.tool_list:
            t.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.download_anno_files()
        self.run_two_component()
        super(TwoComponentWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        for i in range(self.num):

            c_tool = self.tool_list[i]
            #sample_name = self.samples[i]
            sample_name = c_tool.option('sample_name')  #tool_list
            result1 = self.output_dir + '/' + sample_name + '.senser_regulator.xls'
            result2 = self.output_dir + '/' + sample_name + '.senser_regulator.stat'
            if os.path.exists(result1):
                os.remove(result1)
            if os.path.exists(result2):
                os.remove(result2)

            #if not os.path.exists(c_tool.output_dir + '/' + sample_name + '.senser_regulator.xls'):
            self.logger.info("file %s" % c_tool.output_dir + '/' + sample_name + '.senser_regulator.xls')

            #if os.path.exists(c_tool.output_dir + '/' + sample_name + '.senser_regulator.stat'):
            self.logger.info("file: %s" % c_tool.output_dir + '/' + sample_name + '.senser_regulator.stat')

            os.link(c_tool.output_dir + '/' + sample_name + '.senser_regulator.xls', result1)
            os.link(c_tool.output_dir + '/' + sample_name + '.senser_regulator.stat', result2)
            api = self.api.api('bacgenome.common_api')
            mongo_key1 = 'gene_id,type,pfam_id,domain,domain_desc,,,,,,location,'  #没有score
            #mongo_key1 = 'gene_id,type,pfam_id,domain,domain_desc,,,,,,,location,' #
            api.add_main_detail(result1,'anno_regulator_detail', self.option('main_id'), mongo_key1, has_head =True,
                            main_name='regulator_id',other_dic={'specimen_id':sample_name})
            mongo_key2 = 'senser,regulator,hybrid'
            api.add_main_detail(result2,'anno_regulator_stat', self.option('main_id'), mongo_key2, has_head =True,
                            main_name='regulator_id', other_dic={'specimen_id':sample_name},
                            main_table='anno_regulator', update_dic={"main_id":ObjectId(self.option('main_id'))})
        self.end()

    def end(self):
        repaths = [
            [".", "", "双组分调控系统目录",0],
        ]
        regexps = [
            [r'.*\.xls$', 'xls', '双组分调控系统分析结果文件',0],
            [r'.*\.stat$', 'stat', '双组分调控系统分析统计结果',0]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(TwoComponentWorkflow, self).end()
