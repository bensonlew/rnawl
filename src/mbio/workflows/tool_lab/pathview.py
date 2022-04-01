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


class PathviewWorkflow(Workflow):
    """
    Pathview
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PathviewWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_file", "type": "infile", 'format': 'ref_rna_v2.common'},
            {"name": "compound_file", "type": "infile", 'format': 'ref_rna_v2.common'},
            {"name": "type", "type": "string", 'default': "gene"},
            {"name": "data_type", "type": "string", 'default': "continuous"},
            {"name": "gene_id_type", "type": "string", 'default': "entrez"},
            {"name": "compound_id_type", "type": "string", 'default': "kegg"},
            {'name': 'species', 'type': 'string', 'default': 'hsa'},
            {'name': 'pathway', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.pathview")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(PathviewWorkflow, self).run()

    def check_options(self):
        if self.option('type') == 'gene' and not self.option("gene_file").is_set:
            raise OptionError("必须设置输入基因文件")
        if self.option('type') == 'compound' and not self.option("compound_file").is_set:
            raise OptionError("必须设置输入化合物文件")
        if self.option('type') == 'both' and not self.option("gene_file").is_set:
            raise OptionError("数据类型与输入文件不匹配，请检查。")
        if self.option('type') == 'both' and not self.option("compound_file").is_set:
            raise OptionError("数据类型与输入文件不匹配，请检查。")
        if self.option('gene_id_type').lower() not in ['entrezid', 'entrez'] and self.option('data_type') == 'continuous' and not self.option('pathway'):
            raise OptionError("暂不支持该基因ID类型自动检测KEGG pathway。请输入指定pathwayID，以逗号分割，或转换成EntrezID。")
        if self.option('compound_id_type').lower() != 'kegg' and self.option('data_type') == 'continuous' and not self.option('pathway'):
            raise OptionError("暂不支持该化合物ID类型自动检测KEGG pathway。请输入指定pathwayID，以逗号分割，或转换成KEGG Compound ID。")
        if self.option('data_type') != 'continuous' and not self.option('pathway'):
            raise OptionError("请输入指定pathwayID，以逗号分割。")
        return True

    def run_tool(self):
        if len(self.option('type').split(',')) > 1:
            d_type = 'both'
        else:
            d_type = self.option('type')
        if self.option('compound_file').is_set:
            data_type = 'continuous'
        else:
            data_type = self.option('data_type')
        opts = {
            'type': d_type,
            'data_type': data_type,
            'species': self.option('species'),
        }
        if self.option('compound_file').is_set:
            opts.update({'compound_file': self.option('compound_file')})
            opts.update({'compound_id_type': self.option('compound_id_type')})
        if self.option('gene_file').is_set:
            opts.update({'gene_file': self.option('gene_file')})
            opts.update({'gene_id_type': self.option('gene_id_type')})
        if self.option('pathway'):
            opts.update({'pathway': self.option('pathway').lower()})
        self.tool.set_options(opts)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        pathview = self.api.api("tool_lab.pathview")
        # add result info
        result = glob.glob(os.path.join(self.tool.output_dir, '*png'))
        png_dict = dict()
        if len(result) > 0:
            for i in result:
                name = os.path.basename(i).split('.', 1)[0]
                png_dict[name] = os.path.join(self._sheet.output, os.path.basename(i))
        pathview.add_pathview(self.option('main_id'), result=png_dict)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "pathview图结果目录",0],
            ['./*png', 'PNG', 'pathview结果图', 0],
            ['./*txt', 'TXT', '富集结果文件', 0],
        ])
        super(PathviewWorkflow, self).end()
