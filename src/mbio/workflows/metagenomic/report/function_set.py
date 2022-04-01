# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifiy = modified 20181122

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import os
import types


class FunctionSetWorkflow(Workflow):
    """
    宏基因组功能集
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FunctionSetWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "database", "type": "string"},  # database
            {"name": "level", "type": "string"},  # 创建功能集的功能水平
            {"name": "member", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info","type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

        self.func_tool = self.add_tool('annotation.function_set')


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start FunctionSet!")
        self.run_func()
        super(FunctionSetWorkflow, self).run()

    def run_func(self):
        self.func_tool.set_options({
            'database': self.option('database'),
            'level' : self.option('level'),
            'member' : self.option('member'),
            'main_id' : self.option('main_id')
        })
        self.func_tool.on("end", self.set_db)
        self.func_tool.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        #out = self.output_dir + '/function_detail.xls'
        #if os.path.exists(out):
        #    os.remove(out)
        #os.link(self.func_tool.output_dir + '/function_detail.xls', out)

        #api_func = self.api.api("metagenomic.function_set")
        #self.logger.info("正在写入mongo数据库")
        #main_id = self.option("main_id")
        #if not isinstance(main_id, ObjectId):
        #    if isinstance(main_id, types.StringTypes):
        #        main_id = ObjectId(main_id)
        #    else:
         #       raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        #api_func.add_func_set_detail(out, main_id)
        self.end()

    def end(self):
        #result_dir = self.add_upload_dir(self.output_dir)
        #result_dir.add_relpath_rules([
        #    [".", "", "功能集目录", 0, "120214"],
        #    ["function_detail.xls", "xls", "功能集表", 0, "120215"],
        #])
        super(FunctionSetWorkflow, self).end()