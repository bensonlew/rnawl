# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modifies 20180323

import os, glob
import pandas as pd
from biocluster.workflow import Workflow
from biocluster.config import Config


class RefalignDnabacWorkflow(Workflow):
    """
    细菌基因组参考基因组注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(RefalignDnabacWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample_seqpath", "type": "string"},  # "{'a': 'apath', 'b': 'bpath'}"
            {"name": "ref", "type": "string", "default": ""},  # 整理的参考基因组对应Genbank号
            {"name": "ref_path", "type": "infile","format": "gene_structure.gbk"},  # 上传的gbk文件路径
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},  # "{'a':main_id , 'b': main_id}"
            {"name": "task_type", "type": "string"},
            {"name": "main_name", "type": "string"},  # "{'a':main_name , 'b': main_name}"
            {"name": "params", "type": "string"},
            {"name": "evalue", "type":"float","default":10}  # 20190326 zouguanqing
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_list = []
        self.query = eval(self.option('sample_seqpath'))
        self.main_id = eval(self.option('main_id'))
        self.main_name = eval(self.option('main_name'))
        self.ref_tab = {}

    def run_refalign_dnabac(self):
        self.logger.info(self.query)
        self.logger.info(self.main_id)
        self.logger.info(self.main_name)
        for k in sorted(self.query.keys()):
            self.ref = self.add_tool("align.refalign_dnabac")
            options = {
                'query': self.query[k],
                'evalue': self.option('evalue')  #增加evalue
            }
            if self.option('ref') != "":
                options['refalign_database'] = self.option('ref')
            elif self.option('ref_path'):
                options['reference_gbk'] = self.option('ref_path')
            self.logger.info(self.query[k])
            self.ref.set_options(options)
            self.ref_tab[k] = self.ref.output_dir
            self.logger.info(self.ref_tab[k])
            self.tool_list.append(self.ref)
        if len(self.tool_list) > 1:
            self.on_rely(self.tool_list, self.set_db)
            self.logger.info(self.tool_list)
        else:
            self.tool_list[0].on('end', self.set_db)
        for tool in self.tool_list:
            tool.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_refalign_dnabac()
        super(RefalignDnabacWorkflow, self).run()

    def set_db(self):
        api_ref = self.api.api('bacgenome.anno_ref')
        for k in sorted(self.query.keys()):
            self.logger.info(">>>>>>>>>>>>>>>>>>>>>>")
            ref_tab = glob.glob(self.ref_tab[k] + "/*")[0]
            anno = self.output_dir + "/" + self.main_name[k] + '/' + os.path.basename(ref_tab)
            if not os.path.exists(self.ref_tab[k] + "/" + self.main_name[k]):
                os.mkdir(self.output_dir + "/" + self.main_name[k])
            self.logger.info(ref_tab)
            self.logger.info(anno)
            self.anno_table_tidy(ref_tab, anno)
            api_ref.add_anno_ref_detail(self.main_id[k], k, anno, self.option("task_id"))
        self.end()

    def anno_table_tidy(self, anno_table, table_out_path):
        anno_table = pd.read_table(anno_table, sep='\t', header=0)
        anno_table = anno_table.ix[:,
                     ["Query-Name", "Q-Len", "Q-Begin", "Q-End", "HSP-Len", "Hit-Description", "Hit-Name", "Hit-Len",
                      "Hsp-Begin", "Hsp-End", "Identity-%", "E-Value", "Score"]]
        anno_table.rename(
            columns={'Query-Name': 'Gene ID', 'Hit-Name': 'Hit', 'Q-Len': 'Gene Len', 'Q-Begin': 'Gene Start',
                     'Q-End': 'Gene End', 'HSP-Len': 'Hit Len', 'Hsp-Begin': 'Hit Start',
                     'Hsp-End': 'Hit End', 'Identity-%': 'Identity', 'E-Value': 'Evalue', "Hit-Len": "Ref Len"},
            inplace=True)
        table_out = anno_table.ix[:,
                    ["Hit", "Ref Len", "Hit-Description", "Gene Len", "Gene Start", "Gene End", "Hit Len",
                     "Hit Start", "Hit End", "Identity", "Evalue", "Score"]]
        table_out.index = anno_table["Gene ID"]
        table_out.to_csv(table_out_path, sep="\t")
        return table_out

    def end(self):
        repaths = [
            [".", "", "基于参考基因组注释结果目录",0,'130511'],
        ]
        regexps = [
            [r'.*\.xls$', 'xls', '基于参考基因组注释结果文件',0,'130513']
        ]
        dir = self.add_upload_dir(self.output_dir)
        dir.add_relpath_rules(repaths)
        dir.add_regexp_rules(regexps)
        super(RefalignDnabacWorkflow, self).end()
