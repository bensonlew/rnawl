# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import pandas as pd
import json
from collections import OrderedDict
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class AsprofileSearchWorkflow(Workflow):
    """
    可变剪切事件

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AsprofileSearchWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "as_result_merge", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'target_genes', 'type': 'string'},
            {'name': 'geneset_asprofile', 'type': 'string',"default":None},
            {'name': 'type', 'type': 'string'},
            #id_type, [gene_name,gene_id]
            {'name': 'id_type', 'type': 'string',"default":"gene_id"},
            {'name': 'anno', 'type': "infile", 'format': "ref_rna_v2.common"},
            {'name': 'sample', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 Gene_structure_analysis/01 AS/07 search_for_ASprofile')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(AsprofileSearchWorkflow, self).send_log(data)

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))


    def run(self):
        self.get_run_log()
        self.prepare_all_genes()
        self.run_asprofile_search()
        super(AsprofileSearchWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_asprofile_search", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def prepare_all_genes(self):
        target_genes = set()
        final_target_genes = set()
        geneset_genes =set()
        if self.option("id_type") == "gene_id":
            if self.option("target_genes") != "":
                for i in self.option("target_genes").strip().split(","):
                        target_genes.add(i)
                with open(self.option("geneset_asprofile")) as g:
                    setgenes = g.readlines()
                    if len(setgenes)> 0:
                        for i in setgenes:
                            geneset_genes = geneset_genes | set(i.strip().split(","))
                        for i in target_genes:
                            if i in geneset_genes:
                                final_target_genes.add(i)
                            else:
                                pass
                        if len(final_target_genes) == 0:
                            self.set_error('no target gene in geneset')
                    else:
                            final_target_genes =  target_genes
            else:
                with open(self.option("geneset_asprofile")) as g:
                    setgenes = g.readlines()
                    if len(setgenes)> 0:
                        for i in setgenes:
                            geneset_genes = geneset_genes | set(i.strip().split(","))
                        final_target_genes = geneset_genes
                    else:
                        pass
        else:
            name2geneid=self.get_name2id(self.option("anno").prop["path"])
            if self.option("target_genes") != "":
                for i in self.option("target_genes").strip().split(","):
                        target_genes.add(name2geneid[i])
                with open(self.option("geneset_asprofile")) as g:
                    setgenes = g.readlines()
                    if len(setgenes)> 0:
                        for i in setgenes:
                            geneset_genes = geneset_genes |set(i.strip().split(","))
                        for i in target_genes:
                            if i in geneset_genes:
                                final_target_genes.add(i)
                            else:
                                pass
                        if len(final_target_genes) == 0:
                            self.set_error('no target gene in geneset')
                    else:
                          final_target_genes =  target_genes
            else:
                with open(self.option("geneset_asprofile")) as g:
                    setgenes = g.readlines()
                    if len(setgenes)> 0:
                        for i in setgenes:
                            geneset_genes = geneset_genes | set(i.strip().split(","))
                        final_target_genes = geneset_genes
                    else:
                        pass
        if len(final_target_genes) != 0:
            with open(os.path.join(self.work_dir, "final_target_genes"), "w") as ft:
                for i in final_target_genes:
                    ft.write(i + "\n")
            self.final_target_genes = os.path.join(self.work_dir, "final_target_genes")
        else:
            self.final_target_genes = "all"

    def get_name2id(self,anno_path):
        gene_annot = pd.read_table(anno_path)
        gene_annot = gene_annot.loc[:, ['gene_id', 'gene_name']]
        gene_annot = gene_annot.dropna(subset=["gene_name"], axis=0)
        df_gene_tmp = gene_annot.drop_duplicates('gene_name', keep='first', inplace=False)
        name2iddict = OrderedDict(zip(df_gene_tmp.gene_name, df_gene_tmp.gene_id))
        return name2iddict

    def run_asprofile_search(self):
        self.as_search = self.add_tool('medical_transcriptome.ASprofile.asprofile_search')
        self.as_search.set_options({
                'as_result_merge': self.option('as_result_merge'),
                # 'target_genes': os.path.join(self.work_dir,"final_target_genes"),
                'target_genes': self.final_target_genes,
                'event_type':self.option("type"),
                'sample': self.option("sample")
            })
        self.as_search.on('end', self.set_db)
        self.as_search.run()

    def end(self):
        if os.path.exists(os.path.join(self.as_search.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.as_search.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.as_search.output_dir, os.path.basename(self.run_log)))
        os.rename(os.path.join(self.as_search.output_dir, 'target_detail.txt'),
                  os.path.join(self.as_search.output_dir, 'ASprofile_detail.xls'))
        result_dir = self.add_upload_dir(self.as_search.output_dir)
        self.inter_dirs = [
            ["03 Gene_structure_analysis", "", "基因结构分析数据挖掘结果目录", 0],
            ["03 Gene_structure_analysis/01 AS", "", "可变剪切分析结果目录", 0],
            ["03 Gene_structure_analysis/01 AS/07 search_for_ASprofile", "", "目标可变剪切时间查询结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "", 0],
            ['ASprofile_detail.xls', 'xls', '可变剪切ASprofile搜索结果文件'],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(AsprofileSearchWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        asprofile = self.api.api("medical_transcriptome.asprofile")
        # add result info
        as_search = os.path.join(self.as_search.output_dir, 'target_detail.txt')
        asprofile.add_asprofile_search_result(as_search, self.option('main_id'))
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.asprofile_search import AsprofileSearchWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "Diff_ASprofile" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.report.asprofile_search",
            "options": dict(
                as_result_merge='/mnt/ilustre/users/sanger-dev/workspace/20200813/Asprofile_tsg_38314_8419_7589/Asprofile/output/AS_result_merge.txt',
                target_genes='AK2,E2F2',
                # target_genes='ENSG00000004455,ENSG00000007968',
                sample = "SNU16_7",
                type ="TSS",
                anno = "/mnt/ilustre/users/sanger-dev/workspace/20200813/GenesetPpi_tsg_38314_2767_2975/seq_annot.xls",
                id_type = "gene_name"
            )
        }

        wsheet = Sheet(data=data)
        wf =AsprofileSearchWorkflow(wsheet)
        wf.sheet.id = 'asprofile_search'
        wf.sheet.project_sn = 'asprofile_search'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
