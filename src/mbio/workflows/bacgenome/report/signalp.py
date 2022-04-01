# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modifies 20180320

"""旁系同源分析"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime


class SignalpWorkflow(Workflow):
    """
    报告中使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(SignalpWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "query", "type": "string"},
            {"name": "type", "type": "string", "default": "gram-"},  # 菌株类型gram-，gram+，euk, bac
            {"name": "specimen_id", "type": "string"},  # 样品名称
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "main_name", "type": "string"},
            {"name": "params", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.signalp = self.add_module("bacgenome.signalp")

    def run_signalp(self):
        self.spe_list = self.option('specimen_id').split(",")
        self.sig_dir = {}
        options = {
            'query': self.option('query'),
            #'type': self.option('type'),
            'type': "bac",
            "sample": self.option('specimen_id')
        }
        self.signalp.set_options(options)
        self.signalp.on('end', self.set_db)
        self.signalp.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_signalp()
        super(SignalpWorkflow, self).run()

    def set_db(self):
        self.logger.info("start mongo>>>>>>>>>>>>")
        api_signalp = self.api.api('bacgenome.anno_signalp')
        id = self.option('main_id')
        sig_dir = self.signalp.output_dir
        for i in self.spe_list:
            #each_dir = os.path.join(sig_dir, i)
            for type in  ['Gram-','Gram+']:
                each_dir = '{}/{}_{}_SignalP.txt'.format(sig_dir,i,type)
                self.logger.info(each_dir)
                api_signalp.add_anno_signalp_detail(each_dir, type, main_id=id)
        self.linkdir(sig_dir)
        self.end()

    def linkdir(self, dirpath, dirname=None):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)

        if dirname:
            newdir = os.path.join(self.output_dir, dirname)
        else:
            newdir = self.output_dir
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        self.logger.info(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def end(self):
        repaths = [
            [".", "", "分泌蛋白预测结果目录",0,'130505'],
        ]
        regexps = [
            [r'.*SignalP_.*_Gram+.*', '', '革兰氏阳性菌分泌蛋白预测结果目录',0,'130501'],
            [r'.*SignalP_.*_Gram-.*', '', '革兰氏阴性菌分泌蛋白预测结果目录',0,'130502'],
            [r'.*_Gram\+_SignalP\.txt$', 'txt', '革兰氏阳性菌分泌蛋白预测结果文件',0,'130503'],
            [r'.*_Gram-_SignalP\.txt$', 'txt', '革兰氏阴性菌分泌蛋白预测结果文件',0,'130504']
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(SignalpWorkflow, self).end()
