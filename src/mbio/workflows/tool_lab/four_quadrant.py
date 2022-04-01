# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from biocluster.core.exceptions import OptionError


class FourQuadrantWorkflow(Workflow):
    """
    Four Quadrant Plot
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FourQuadrantWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'dataframe1', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'dataframe2', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'id_map', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'fc1', 'type': 'string'},
            {'name': 'fc2', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.four_quadrant")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(FourQuadrantWorkflow, self).run()

    def check_options(self):
        if not self.option("dataframe1").is_set:
            raise OptionError("请上传dataframe1数据表")
        if not self.option("dataframe2").is_set:
            raise OptionError("请上传dataframe2数据表")
        return True

    def run_tool(self):
        opts = {
            'dataframe1': self.option('dataframe1'),
            'dataframe2': self.option('dataframe2'),
        }
        if self.option('id_map').is_set:
            opts.update({'id_map': self.option('id_map')})
        if self.option('fc1'):
            opts.update({'fc1': float(self.option('fc1'))})
        if self.option('fc2'):
            opts.update({'fc2': float(self.option('fc2'))})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        four_quadrant = self.api.api("tool_lab.four_quadrant")
        # add result info
        result = glob.glob(os.path.join(self.tool.output_dir, '*png'))[0]
        png = os.path.join(self._sheet.output, os.path.basename(result))
        four_quadrant.add_four_quadrant(self.option('main_id'), result=png)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "四象限图结果目录",0],
            ['./*png', 'PNG', '四象限图PNG', 0],
            ['./*pdf', 'PDF', '四象限图PDF', 0],
            ['./*txt', 'TXT', '四象限图结果文件', 0],
        ])
        super(FourQuadrantWorkflow, self).end()
