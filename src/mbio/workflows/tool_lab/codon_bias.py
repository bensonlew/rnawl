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



class CodonBiasWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CodonBiasWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "sequence.fasta"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},

        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.logger.info('seq_num:'+ str(self.option('table').prop["seq_number"]))
        if int(str(self.option('table').prop["seq_number"]).strip()) != 1:
            self.logger.error('输入文件必须是1条核苷酸序列')
            self.set_error('输入文件必须是1条核苷酸序列')

        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")
        self.run_tool()
        super(CodonBiasWorkflow, self).run()


    def run_tool(self):
        self.this_tool = self.add_tool("tool_lab.codon_bias")
        options = {
            "fasta" : self.option("table")
        }
        self.this_tool.set_options(options)
        self.this_tool.on('end',self.set_db)
        self.this_tool.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_name = self.api.api("tool_lab.codon_bias")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        api_name.add_detail(self.this_tool.output_dir+'/result.txt',main_id)
        pic_path = self._sheet.output + '/Sample_codon.png'
        api_name.add_stats(self.this_tool.output_dir+'/fasta.result',main_id, pic_path)
        self.end()


    def end(self):
        result_dir = self.add_upload_dir(self.this_tool.output_dir)
        relpath_rules = [
            [".", "", "CodonBias结果文件夹", 0, ],
        ]
        regexps = [
            [r"result.txt", "txt", "明细表", 0],
            [r"fasta.result", "", "统计表", 0],
            ["Sample_codon.png",'','图']
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(CodonBiasWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/codon/input/BJ.circle.fasta",
            "main_id" : "5e4ce5e817b2bf4b326f4b20"
        }
    }

    wsheet = Sheet(data=data)

    wf = CodonBiasWorkflow(wsheet)
    wf.run()