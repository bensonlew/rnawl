# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types



class MultiAlignWorkflow(Workflow):
    """
    代谢样本相关性分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MultiAlignWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "align_method", "type": "string", "default": "muscle"},  #muscle ,mafft, clustalo
            {"name": "out_format", "type": "string", "default": "fasta"},
            {"name": "main_id", "type": "string", "default": ""},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")

        self.run_tool()

        super(MultiAlignWorkflow, self).run()


    def run_tool(self):
        self.this_tool = self.add_tool("tool_lab.multi_align")
        options = {
            "fasta" : self.option("fasta"),
            "align_method" : self.option("align_method"),
            "out_format" : self.option("out_format")
        }
        self.this_tool.set_options(options)
        self.this_tool.on('end',self.set_db)
        self.this_tool.run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        # api_name = self.api.api("tool_lab.multi_align")
        # self.logger.info("正在写入mongo数据库")
        # main_id = self.option("main_id")
        #
        # api_name.add_detail(self.this_tool.output_dir,main_id)
        self.end()


    def end(self):
        result_dir = self.add_upload_dir(self.this_tool.output_dir)
        relpath_rules = [
            [".", "", "比较结果文件夹", 0, ],
        ]
        regexps = [
            [r".*/.*", "", "比对结果", 0],

        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(MultiAlignWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "fasta" : "/mnt/ilustre/users/sanger-dev/sg-users/guhaidong/CS/test_20strains/cut_fa/t.fa",
            "main_id" : "5e4ce5e817b2bf4b326f4b20"
        }
    }

    wsheet = Sheet(data=data)
    wf = MultiAlignWorkflow(wsheet)
    wf.run()