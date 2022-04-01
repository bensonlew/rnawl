# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

"""宏基因组交互分析公共类"""

#import os
#import re
#import types
from biocluster.workflow import Workflow


class ReportWorkflow(Workflow):
    """
    宏基因组交互分析公共类
    """

    def __init__(self, wsheet_object):
        self.__sheet = wsheet_object
        super(ReportWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},  # 数据库注释表格
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "level_id", "type": "string"},
            {"name": "level_type_name", "type": "string"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name": "lowest_level", "type": "string", "default": ""},  # 注释表数据库对应的最低分类，eg：KEGG的ko
        ]
        self.add_option(options)
        self.abundance = self.add_tool("meta.create_abund_table")

    def run_abundance(self, *func):
        options = {
            'anno_table': self.option('anno_table'),
            'geneset_table': self.option('geneset_table'),
            'level_type': self.option('level_id'),
            'level_type_name': self.option('level_type_name'),
            'gene_list': self.option('gene_list'),
            'lowest_level': self.option('lowest_level'),
        }
        self.abundance.set_options(options)
        for f in func:
            self.abundance.on('end', f)
        self.abundance.run()