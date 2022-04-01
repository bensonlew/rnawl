# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class RelationshipSankeyWorkflow(Workflow):
    """
    丰度桑基图小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationshipSankeyWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"input_table", "type": "infile", "format": "small_rna.common"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.relationship_sankey = self.add_tool('tool_lab.relationship_sankey')
    
    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_table"):
            raise OptionError("请输入表格文件")

    def run_sankey(self):
        self.relationship_sankey.set_options({
            'input_table':self.option('input_table'),
        })
        self.relationship_sankey.on('end',self.set_db)
        self.relationship_sankey.run()
    
    def set_db(self):
        self.logger.info("开始导表")
        api_sankey = self.api.api("tool_lab.sankey")
        api_sankey.add_detail(self.option("main_id"),
                        os.path.join(self.relationship_sankey.output_dir,"plot.txt"))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_sankey()
        super(RelationshipSankeyWorkflow, self).run()\

    def end(self):
        super(RelationshipSankeyWorkflow, self).end()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    data = {
        'name': 'test_sankey',
        'id': 'sankey_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        "options": {
            "input_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/relationship_sankey/otu.txt",
            "main_id" : "5e9e6a6017b2bf2049a81be6"
                }
    }
    wsheet = Sheet(data=data)
    wf = RelationshipSankeyWorkflow(wsheet)
    wf.run()