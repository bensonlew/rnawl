#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/30 10:25
@file    : geneset_gsva.py
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

class GenesetGsvaWorkflow(Workflow):
    """
    基因集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetGsvaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            # {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            # {"name": "task_type", "type": "string"},
            # {"name": "task_id", "type": "string"},

            {"name": "group", "type": "infile", 'format': 'sample.group_table'},

            {"name": 'source', 'type': 'string', 'default': 'msigdb'},# [custom, msigdb]
            {"name": 'custom_geneset', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {"name": "matrix", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'species', 'type': 'string'},
            {'name': 'c1', 'type': 'string'},
            {'name': 'c2', 'type': 'string'},
            {'name': 'genesets_1', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {"name": 'method', 'type': 'string', 'default': 'gsva'},
            {'name': 'kcdf', 'type': 'string', 'default': 'Gaussian'},
            {'name': 'min_num', 'type': 'float', 'default': 10},
            {'name': 'max_num', 'type': 'float', 'default': 500},
            {'name': 'es_score', 'type': 'string', 'default': 'True'},
            {'name': 'diff', 'type': 'bool', 'default': False},
            {'name': 'cmp', 'type': 'string'},
            {'name': 'fc', 'type': 'float', 'default': 2},
            {'name': 'pvalue_padjust', 'type': 'string', 'default': 'padjust'},
            {'name': 'pvalue', 'type': 'float', 'default': '0.05'}


            # {"name": "genes_detail", "type": "infile", "format": "ref_rna_v2.ref_common"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        # self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/06 Advanced_Analysis/03 GSEA_Analysis')
        # self.inter_dirs = []
        self.step.add_steps('gsva_analysis')
        # if self.option("go_list").is_set:
        #     self.go_gene_gmx = self.add_tool('ref_rna_v2.geneset.go_gene_gmx')

        self.gsva_tool = self.add_tool('tool_lab.gsva')
        self.cluster_heatmap = self.add_tool('tool_lab.exp_cluster2gsva')
        self.gsva_tool.on('end', self.run_cluster_heatmap)
        if self.option('diff') == True:
            self.cluster_heatmap.on('end', self.run_diff)
        else:
            self.cluster_heatmap.on('end', self.set_db)
        self.group_dict = dict()
        # with open(self.option('group').path, 'r') as file:
        #     for line in file.readlines():
        #         if line.startswith('#'):
        #             continue
        #         else:
        #             sample, group = line.strip().split('\t')
        #             if group not in self.group_dict:
        #                 self.group_dict[group] = [sample]
        #             else:
        #                 self.group_dict[group].append(sample)
        # self.group_dict_json = json.dumps(self.group_dict)
        # self.g1, self.g2 = self.group_dict.keys()

    def gsva_runner(self):

        if self.option('source') == 'msigdb':
            geneset = ';'.join([self.option('species'), self.option('c1'), self.option('c2'), self.option('genesets_1')])
            self.gmx = export_msigdb_genesets(geneset, self.work_dir)
        if self.option('source') == 'custom':
            self.gmx = self.option('custom_geneset')
        self.gsva_tool.set_options({
            'gmx': self.gmx,
            'matrix': self.option('matrix'),
            'method': self.option('method'),
            'kcdf': self.option('kcdf'),
            'min_num': self.option('min_num'),
            'max_num': self.option('max_num'),
            'es_score': self.option('es_score'),
            'species': self.option('species')
        })
        self.gsva_tool.on('start', self.set_step, {'start': self.step.gsva_analysis})
        self.gsva_tool.on('end', self.set_step, {'end': self.step.gsva_analysis})
        self.gsva_tool.run()

    def run_cluster_heatmap(self):
        options = dict(
            exp=self.gsva_tool.option('es_exp_table').path,
            group=self.option('group'),
            # use_group=self.option('use_group'),
        )
        self.cluster_heatmap.set_options(options)
        self.cluster_heatmap.on('start', self.set_step, {'start': self.step.gsva_analysis})
        self.cluster_heatmap.on('end', self.set_step, {'end': self.step.gsva_analysis})
        self.cluster_heatmap.run()

    def run_diff(self):
        self.diff = self.add_tool('tool_lab.diff_geneset')
        options = dict(
            count=self.gsva_tool.option('es_exp_table').path,
            group=self.option('group'),
            cmp=self.option('cmp'),
            fc=self.option('fc'),
            pvalue_padjust=self.option('pvalue_padjust'),
            pvalue=self.option('pvalue')
        )
        self.diff.set_options(options)
        self.diff.on('start', self.set_step, {'start': self.step.gsva_analysis})
        self.diff.on('end', self.set_step, {'end': self.step.gsva_analysis})
        self.diff.on('end', self.set_db)
        self.diff.run()

    def run(self):
        self.gsva_runner()
        super(GenesetGsvaWorkflow, self).run()


    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_gsva = self.api.api('tool_lab.gsva')
        main_id = api_gsva.add_gsva(self.option("main_id"), self.gsva_tool.option('es_exp_table').path, self.gmx)
        api_gsva.add_gsva_cluster(self.cluster_heatmap.output_dir, main_id)
        if self.option('diff') == True:
            api_gsva.add_gsva_diff(self.diff.output_dir, main_id, self.option('group').path)
        self.set_output()
        self.end()

    def set_output(self):
        for file in os.listdir(self.gsva_tool.output_dir):
            os.link(os.path.join(self.gsva_tool.output_dir, file), os.path.join(self.output_dir, file))

    def end(self):
        result_dir = self.add_upload_dir(self.gsva_tool.output_dir)
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
        super(GenesetGsvaWorkflow, self).end()

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
    out_file = os.path.join(output_dir, 'gsva.gmt')
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
        from mbio.workflows.tool_lab.geneset_gsva import GenesetGsvaWorkflow
        from biocluster.wsheet import Sheet
        import random
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "gsva_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.geneset_gsva",
            "options": dict(
                matrix='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/GSEA/gsea.txt',
                group="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/GSEA/group_test.txt",
                species="Mus musculus",
                source='msigdb',
                c1="C2: curated gene sets",
                c2="CGP: chemical and genetic perturbations",
                genesets_1="LANDIS_BREAST_CANCER_PROGRESSION_UP,LANDIS_BREAST_CANCER_PROGRESSION_DN,GESERICK_TERT_TARGETS_DN",
                max_num='500',
                min_num='15',
                cmp='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/GSVA/compare_1.txt'
            )
        }
        wsheet = Sheet(data=data)
        wf =GenesetGsvaWorkflow(wsheet)
        wf.sheet.id = 'extract_gff_fasta'
        wf.sheet.project_sn = 'extract_gff_fasta'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)