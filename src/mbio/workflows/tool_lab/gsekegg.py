# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class GsekeggWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GsekeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneList", "type": "infile", "format": "ref_rna_v2.common"},
            # A csv file contains two columns, one for gene ID (no duplicated allowed) and another one for fold change.
            {"name": "organism", "type": "string", "default": ""},
            # such as 'Arabidopsis thaliana', 'arabidopsis_thaliana', 'Arabidopsis_thaliana'
            {"name": "id_type", "type": "string", "default": "EntrezGene"},
            # EntrezGene or EnsemblGene
            {"name": "nPerm", "type": "int", "default": 1000},
            # The number of permutations, default 1000.
            {"name": "minGSSize", "type": "int", "default": 10},
            # Gene sets smaller than this number are EXLCUDED from the analysis.
            {"name": "maxGSSize", "type": "int", "default": 500},
            # Gene sets bigger than this number are EXLCUDED from the analysis.
            {"name": "pvalueCutoff", "type": "float", "default": 0.05},
            {"name": "nGenesets", "type": "int", "default": 5},
            # The number of gene sets to be displayed on the same figure.
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.gsekegg")
        self.api = self.api.api("tool_lab.api_base")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_tool()
        super(GsekeggWorkflow, self).run()

    def check_options(self):
        if not self.option('geneList').is_set:
            raise OptionError('基因集文件必须输入')
        if not self.option('organism'):
            raise OptionError('物种信息必须输入')
        if self.option('id_type').lower() not in ['entrezgene', 'ensemblgene']:
            raise OptionError('目前只支持EntrezGene或EnsemblGene两种类型的Gene ID！')
        return True

    def run_tool(self):
        opts = {
            'geneList': self.option('geneList'),
            'organism': self.option('organism'),
            'id_type': self.option('id_type'),
            'nPerm': self.option('nPerm'),
            'minGSSize': self.option('minGSSize'),
            'maxGSSize': self.option('maxGSSize'),
            'pvalueCutoff': self.option('pvalueCutoff'),
            'nGenesets': self.option('nGenesets'),
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_output)
        self.tool.run()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            if os.path.exists(os.path.join(self.output_dir, file)):
                os.remove(os.path.join(self.output_dir, file))
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        query_dict = {"main_id": ObjectId(self.option("main_id"))}
        pdf_list = [os.path.join(self._sheet.output, "gseaplot.pdf"), os.path.join(self._sheet.output, "gseaupset.pdf")]
        update_dict = {"output_dir": pdf_list, "status": "end"}
        self.api.update_db_record('gsekegg', query_dict, update_dict)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "多通路GSEA富集分析结果文件",0],
            [r'gsea.txt', 'txt', '富集分析结果文件', 0],
            [r'gseaplot.pdf', 'pdf', '富集分析plot PDF图', 0],
            [r'gseaupset.pdf', 'pdf', '富集分析upset PDF图', 0],
            [r'gseaplot.png', 'png', '富集分析plot PNG图', 0],
            [r'gseaupset.png', 'png', '富集分析upset PNG图', 0],
        ])
        super(GsekeggWorkflow, self).end()