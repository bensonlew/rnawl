# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modifies 20180319

"""旁系同源分析"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime


class ParalogsWorkflow(Workflow):
    """
    报告中使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(ParalogsWorkflow, self).__init__(wsheet_object)
        options = [
            #{"name": "gene_file", "type": "infile", "format": "sequence.fasta"},
            {"name": "sequence", "type": "infile", "format": "sequence.fasta"},
            {"name": "task_id", "type": "string"},
            {"name": "specimen_id", "type": "string"},
            # {"name": "old_species", "type": "string"},
            #{"name": "search_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "params", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.align = self.add_tool("align.refalign_dnabac")

    def run_align(self):
        # with open(self.option('sequence').prop['path']) as fr,open(self.work_dir+'/new_faa','w') as fw:
        #     for line in fr:
        #         if line[0] == '>':
        #             fw.write(line.split()[0]+'_1\n')
        #         else:
        #             fw.write(line)
        os.link(self.option('sequence').prop['path'], self.work_dir+'/new_faa')
        options = {
            #'query': self.option('gene_file'),
            'query': self.option('sequence'),
            'reference_seq': self.work_dir+'/new_faa',
            #"query_type": "nucl",
            "query_type": "prot",
            #"blast": "blastn",
            "blast": "blastp",
            "num_alignment": 1000,
            #"reference_type": "nucl"
            "reference_type": "prot",
            "evalue": 0.00001,
            "coverage":'T',
            "analysis_type" : "paralog"
        }
        self.align.set_options(options)
        self.align.on('end', self.set_db)
        self.align.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_align()
        super(ParalogsWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        result = self.output_dir + '/' + self.option('specimen_id') + '_paralogs.xls'
        if os.path.exists(result):
            os.remove(result)
        file = os.listdir(self.align.output_dir)
        os.link(self.align.output_dir + '/' + file[0], result)
        #self.option('gene_file').get_info()
        #len = self.option('gene_file').prop['bases']
        api_paralogs = self.api.api('bacgenome.paralogs')
        api_paralogs.add_paralogs(result, main=False, main_id=self.option('main_id'),
                                  len=len, specimen_id=self.option('specimen_id'),
                                  task_id=self.option('task_id'))
        self.end()

    def end(self):
        repaths = [
            [".", "", "旁系同源基因查询目录",0,'130514'],
        ]
        regexps = [
            [r'.*_paralogs.*\.xls$', 'xls', '旁系同源基因查询结果文件',0,'130515']
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(ParalogsWorkflow, self).end()
