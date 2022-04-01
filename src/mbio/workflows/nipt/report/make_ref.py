# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# create_time：20170821
from biocluster.workflow import Workflow
import os
from biocluster.config import Config


class MakeRefWorkflow(Workflow):
    """
    该workflow用于nipt构建参考组的接口
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MakeRefWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample_txt", "type": "infile", "format": "nipt.bed"},  # 输入样本的id
            {"name": "ref_group", "type": "int", "default": 1},
            {"name": "update_info", "type": "string"},
            {"name": "nipt_task_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.MakeRef = self.add_tool('nipt.make_ref')
        self.bed_dir = self.config.SOFTWARE_DIR + "/database/human/hg38_nipt/bed_file"
        self.ref_cor = self.config.SOFTWARE_DIR + "/database/human/hg38_nipt/ref_cor"

    def run_niptanalysis(self):
        options = {
            'sample_txt': self.option("sample_txt").prop['path'],
            'bed_dir': self.bed_dir,
            'ref_cor': self.ref_cor,
            'ref_group': self.option('ref_group')
        }
        self.MakeRef.set_options(options)
        self.MakeRef.on('end', self.end)
        self.MakeRef.run()


    def end(self):
        super(MakeRefWorkflow, self).end()

    def run(self):
        self.run_niptanalysis()
        super(MakeRefWorkflow, self).run()
