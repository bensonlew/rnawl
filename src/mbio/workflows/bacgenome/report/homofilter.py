# -*- coding: utf-8 -*-
# __author__ = 'ysh'

import os
import re
import types
from biocluster.workflow import Workflow
from biocluster.config import Config
from mainapp.models.mongo.bacgenome import Bacgenome
from bson import ObjectId
import datetime
from biocluster.file import download, exists
import json


class HomofilterWorkflow(Workflow):
    """
    核心和特有基因分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HomofilterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimens", "type": "string"},  #所有样本，用来下载注释表使用
            {"name": "stat_file", "type": "infile", "format": "sequence.profile_table"},
            {'name': 'filter_type', 'type': 'string'},
            {'name': 'filter', 'type': 'string'},
            {'name': 'select', 'type': 'string'},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("bacgenome.homofilter")
        self.anno_dir = os.path.join(self.work_dir, "anno_dir")

    def run_homofilter(self):
        self.down_load_files()
        all_samples = self.option("specimens").split(",")
        options = {
            'anno_dir': self.anno_dir,
            'stat_file': self.option('stat_file'),
        }
        if not self.option("filter") and not self.option("select"):
            self.logger.info("不能同时没有过滤样本和筛选样本")
        if self.option("filter_type") == "core":
            select_samples = self.option("select").split(",")
            if sorted(all_samples) != sorted(select_samples):
                self.set_error("选择核心基因时样本不正确")
        if self.option("filter_type") == "unique":
            if self.option("filter"):
                filter_sams = self.option("filter").split(",")
                if len(filter_sams) != (len(all_samples) -1):
                    self.set_error("选择特有样本时filter不正确")
            else:
                filter_sams = self.option("specimens").split(",")
                filter_sams.remove(self.option("select"))
                options["filter"] = ",".join(filter_sams)
        if self.option("filter"):
            options["filter"] = self.option("filter")
        if self.option("select"):
            options["select"] =  self.option("select")
        self.tool.set_options(options)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def down_load_files(self):
        bacgenome = Bacgenome()
        bacgenome._config = Config()
        task_id = self.option("task_id")
        if not os.path.exists(self.anno_dir):
            os.mkdir(self.anno_dir)
        anno_paths = bacgenome.get_anno_summary(task_id, self.option("specimens"))
        for each in anno_paths.keys():
            newfile = os.path.join(self.anno_dir, each + ".xls")
            download(anno_paths[each], newfile)

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_homofilter()
        super(HomofilterWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.linkdir(self.tool.output_dir, self.output_dir)
        homofilter_file = self.tool.option("result").prop["path"]
        api_homofilter = self.api.api('bacgenome.homofilter')
        api_homofilter.add_homofilter_detail(self.option('main_id'), homofilter_file)
        self.end()

    def linkdir(self, dirpath, newdir):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
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
            [".", "", "核心和特有基因分析"],
            ["Homofilter_genes.xls", "xls", "同源基因筛选结果表"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(HomofilterWorkflow, self).end()
