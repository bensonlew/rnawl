#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/30 10:25
@file    : geneset_gsea.py
"""

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import re
from bson.objectid import ObjectId
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import unittest

class GenesetGseaWorkflow(Workflow):
    """
    基因集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetGseaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            # {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            # {"name": "task_type", "type": "string"},
            # {"name": "task_id", "type": "string"},

            {"name": "preranked", "type": "bool", 'default': False},
            {"name": "geneset_source", "type": "string"},
            # {"name": "geneset_id", "type": "string"},

            {"name": "g1", "type": "string"},
            {"name": "g2", "type": "string"},
            {"name": "group", "type": "infile", 'format': 'ref_rna_v2.common'},

            {"name": "max_num", "type": "string", "default": "500"},
            {"name": "min_num", "type": "string", "default": "15"},

            # {"name": "level", "type": "string", "default": "G"},
            # {"name": "go_genesets", "type": "string", "default": "all"},
            # {"name": "go_list", "type": "infile", "format": "ref_rna_v2.common"},
            # {"name": "go_type", "type": "string", "default": ""},

            {"name": "sort_method", "type": "string"},
            {"name": "permutation_type", "type": "string"},
            {"name": "plot_top_x", "type": "string", "default": "20"},

            # {"name": 'source', 'type': 'string', 'default': 'msigdb'},# [custom, msigdb]
            {"name": 'custom_geneset', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # {"name": 'geneset', 'type': 'string'}, #'species;c1;c2;'
            # {"name": "gmx", "type": "infile", "format": "ref_rna_v2.geneset_gmt"},
            {"name": "rnk", "type": "infile", "format": "ref_rna_v2.geneset_rnk"},
            {"name": "matrix", "type": "infile", "format": "ref_rna_v2.ref_common"},
            {'name': 'species', 'type': 'string'},
            {'name': 'c1', 'type': 'string'},
            {'name': 'c2', 'type': 'string'},
            {'name': 'genesets_1', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'control_group', 'type': 'string', 'default': None}

            # {"name": "genes_detail", "type": "infile", "format": "ref_rna_v2.ref_common"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        if self.option('preranked'):
            self.gsea_tool = self.add_tool('tool_lab.gsea_preranked')
        else:
            self.gsea_tool = self.add_tool('tool_lab.gsea')
        self.group_dict = dict()
        with open(self.option('group').path, 'r') as file:
            for line in file.readlines():
                if line.startswith('#'):
                    continue
                else:
                    sample, group = line.strip().split('\t')
                    if group not in self.group_dict:
                        self.group_dict[group] = [sample]
                    else:
                        self.group_dict[group].append(sample)
        self.group_dict_json = json.dumps(self.group_dict)

        self.g1, self.g2 = self.group_dict.keys()
        if self.option('control_group'):
            if self.g1 == self.option('control_group'):
                self.g1 = self.option('control_group')
                self.g2 = self.g2
            elif self.g2 == self.option('control_group'):
                self.g1 = self.g2
                self.g2 = self.option('control_group')
            else:
                raise Exception('比较组不在分组文件中')

    def gsea_runner(self):
        if self.option('geneset_source') == 'msigdb':
            geneset = ';'.join([self.option('species'), self.option('c1'), self.option('c2'), self.option('genesets_1')])
            gmx = export_msigdb_genesets(geneset, self.work_dir)
        if self.option('geneset_source') == 'custom':
            gmx = self.option('custom_geneset').path
        self.gsea_tool.set_options({
            'g1': self.g1,
            'g2': self.g2,
            'group': self.group_dict_json,
            'set_max': int(self.option('max_num')),
            'set_min': int(self.option('min_num')),
            'metric': self.option('sort_method'),
            'plot_top_x': int(self.option('plot_top_x')),
            'permute': self.option('permutation_type'),
            'geneset_source': self.option('geneset_source'),
            'gmx': gmx,
            'matrix': self.option('matrix'),
            'species': self.option('species')
            # 'gene_detail': self.option('genes_detail')
        })
        self.gsea_tool.run()

    def gsea_preranked_runner(self):
        if self.option('geneset_source') == 'msigdb':
            geneset = ';'.join([self.option('species'), self.option('c1'), self.option('c2'), self.option('genesets_1')])
            gmx = export_msigdb_genesets(geneset, self.work_dir)
        if self.option('geneset_source') == 'custom':
            gmx = self.option('custom_geneset').path
        self.gsea_tool.set_options({
            'set_max': self.option('max_num'),
            'set_min': self.option('min_num'),
            'gmx': gmx,
            'rnk': self.option('rnk'),
            'geneset_source': self.option('geneset_source'),
            'matrix': self.option('matrix'),
            'species': self.option('species'),
            'plot_top_x': int(self.option('plot_top_x')),
        })
        self.gsea_tool.run()

    def run(self):
        self.gsea_tool.on("end", self.set_db)
        if self.option("preranked"):
            self.gsea_preranked_runner()
        else:
            self.gsea_runner()
        super(GenesetGseaWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """

        api_geneset = self.api.api('tool_lab.gsea')
        api_geneset.run_webroot(self.option("main_id"), self.gsea_tool.output_dir, self._sheet.output)

        conn = api_geneset.db["sg_geneset_gsea"]

        conn.update({"_id": ObjectId(self.option("main_id"))}, {"$set": {'result_dir': self._sheet.output}}, upsert=True)

        self.set_output()
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        # self.inter_dirs = [
        #     ["06 Advanced_Analysis", "", "高级分析结果目录",0],
        #     ["06 Advanced_Analysis/03 GSEA_Analysis", "", "GSEA分析结果目录", 0]
        # ]
        # result_dir.add_relpath_rules([
        #     [".", "", "GSEA分析文件",0,"211546"],
        #     ["./*.pdf", "", "该基因集富集结果图pdf",0,"211547"],
        #     ["./*.png", "", "该基因集富集结果图png",0,"211548"],
        #     ["./all_sets.detail", "txt", "所有基因集富集结果详情表",0,"211549"],
        #     ["./gsea_report.xls", "xls", "所有基因集富集结果统计表",0,"211550"],
        #     ["./all_exp.detail", "txt", "所有基因集leading基因表达量信息",0,"211551"]
        # ])
        super(GenesetGseaWorkflow, self).end()

    def set_output(self):
        for file in os.listdir(self.gsea_tool.output_dir):
            os.link(os.path.join(self.gsea_tool.output_dir, file), os.path.join(self.output_dir, file))

def export_msigdb_genesets(data, output_dir):
    project_type = 'ref_rna'
    db = Config().get_mongo_client(mtype=project_type, ref=True)[Config().get_mongo_dbname(project_type, ref=True)]
    # {'c1': data.c1, 'c2': data.c2, 'c3': data.c3}
    # bind_obj.logger.info("data is {}".format(data))
    extr_dic = {}
    species= data.split(";")[0]
    extr_dic['organism'] = species
    c1 = data.split(";")[1]
    extr_dic['c1'] = c1
    c2 = data.split(";")[2]
    if c2.count('|') == 1:
        c2, tmp_c3 = c2.split('|')
        extr_dic['c2'] = c2
        extr_dic['c3'] = tmp_c3
    else:
        extr_dic['c2'] = c2

    results = db['sg_msigdb'].find(extr_dic, {'_id': 0})
    ref_link = re.compile('<[^>]*>')
    res_dict = {dic['name']: '{}\t{}\t{}\n'.format(
        dic["name"],
        re.sub(ref_link, '', dic['brief_description']),
        dic["gene_symbols"].replace(" ", "\t")) for dic in results if 'name' in dic}
    c3 = data.split(";")[3]
    out_file = os.path.join(output_dir, 'gsea.gmt')
    with open(out_file, 'w') as out_handler:
        if c3 == 'all' or c3 =='':
            out_handler.write(''.join(res_dict.values()))
        else:
            keys = {i.strip() for i in c3.split(',')}
            out_handler.write(''.join(v for k, v in res_dict.items() if k in keys))
    return out_file


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.geneset_gsea import GenesetGseaWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "gsea_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.geneset_gsea",
            "options": dict(
                matrix='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/GSEA/gsea.txt',
                group="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/GSEA/group_test.txt",
                species="Mus musculus",
                geneset_source='msigdb',
                c1="C2: curated gene sets",
                c2="CGP: chemical and genetic perturbations",
                genesets_1="LANDIS_BREAST_CANCER_PROGRESSION_UP,LANDIS_BREAST_CANCER_PROGRESSION_DN,GESERICK_TERT_TARGETS_DN",
                sort_method='log2_Ratio_of_Classes',
                max_num='500',
                min_num='15'
            )
        }
        wsheet = Sheet(data=data)
        wf =GenesetGseaWorkflow(wsheet)
        wf.sheet.id = 'extract_gff_fasta'
        wf.sheet.project_sn = 'extract_gff_fasta'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)