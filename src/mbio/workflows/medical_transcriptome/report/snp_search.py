# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os
from biocluster.workflow import Workflow
from biocluster.file import getsize, exists
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import pandas as pd
import json
from collections import OrderedDict
from biocluster.file import download
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class SnpSearchWorkflow(Workflow):
    """
    Snp搜索

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpSearchWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "result_path", "type": "string"},
            {'name': 'target_genes', 'type': 'string'},
            {'name': 'geneset_snp', 'type': 'string',"default":None},
            {'name': 'type', 'type': 'string'},
            {'name': 'region', 'type': 'string'},
            #{"大于":"greater","小于":"less","大于等于":greateroreq,"小于等于":"lessoreq","等于":"equal"}
            {'name': 'depth_comare', 'type': 'string'},
            {'name': 'depth', 'type': 'string'},
            #id_type, [gene_name,gene_id]
            {'name': 'id_type', 'type': 'string',"default":"gene_id"},
            {'name': 'anno', 'type': "infile", 'format': "ref_rna_v2.common"},
            {'name': 'sample', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 Gene_structure_analysis/02_SNP_InDel_Analysis')
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
        super(SnpSearchWorkflow, self).send_log(data)

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def prepare_all_genes(self):
        if self.option("type").lower() == "snp":
            result_file = os.path.join(self.option("result_path"),"snp_anno.xls")
            self.result_file=os.path.join(self.work_dir,"snp_anno.xls")
            self.download_s3_file(result_file,self.result_file)
        else:
            result_file = os.path.join(self.option("result_path"), "indel_anno.xls")
            self.result_file = os.path.join(self.work_dir, "indel_anno.xls")
            self.download_s3_file(result_file, self.result_file)
        target_genes = set()
        final_target_genes = set()
        geneset_genes =set()
        if self.option("id_type") == "gene_id":
            if self.option("target_genes") != "":
                for i in self.option("target_genes").strip().split(","):
                        target_genes.add(i)
                with open(self.option("geneset_snp")) as g:
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
                with open(self.option("geneset_snp")) as g:
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
                with open(self.option("geneset_snp")) as g:
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
                with open(self.option("geneset_snp")) as g:
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

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_snp_search", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run(self):
        self.get_run_log()
        self.prepare_all_genes()
        self.run_snp_search()
        super(SnpSearchWorkflow, self).run()

    # def prepare_all_genes(self):
    #     target_genes = set()
    #     final_target_genes = set()
    #     geneset_genes =set()
    #     if self.option("id_type") == "gene_id":
    #         for i in self.option("target_genes").strip().split(","):
    #                 target_genes.add(i)
    #         if self.option("geneset_asprofile"):
    #             with open(self.option("geneset_asprofile")) as g:
    #                 for i in g.readlines():
    #                     geneset_genes.add(i.strip())
    #             for i in target_genes:
    #                 if i in geneset_genes:
    #                     final_target_genes.add(i)
    #                 else:
    #                     pass
    #         else:
    #             final_target_genes =  target_genes
    #     else:
    #         name2geneid=self.get_name2id(self.option("anno").prop["path"])
    #         for i in self.option("target_genes").strip().split(","):
    #                 target_genes.add(name2geneid[i])
    #         if self.option("geneset_asprofile"):
    #             with open(self.option("geneset_asprofile")) as g:
    #                 for i in g.readlines():
    #                     geneset_genes.add(i.strip())
    #             for i in target_genes:
    #                 if i in geneset_genes:
    #                     final_target_genes.add(i)
    #                 else:
    #                     pass
    #         else:
    #             final_target_genes =  target_genes
    #     with open(os.path.join(self.work_dir, "final_target_genes"), "w") as ft:
    #         for i in final_target_genes:
    #             ft.write(i + "\n")

    def get_name2id(self,anno_path):
        gene_annot = pd.read_table(anno_path)
        gene_annot = gene_annot.loc[:, ['gene_id', 'gene_name']]
        gene_annot = gene_annot.dropna(subset=["gene_name"], axis=0)
        df_gene_tmp = gene_annot.drop_duplicates('gene_name', keep='first', inplace=False)
        name2iddict = OrderedDict(zip(df_gene_tmp.gene_name, df_gene_tmp.gene_id))
        return name2iddict


    def run_snp_search(self):
        self.logger.info("参数如下")
        self.logger.info("snp_result_file:{}".format(self.result_file))
        self.logger.info("target_genes:{}".format(self.final_target_genes))
        self.logger.info("type:{}".format(self.option("type")))
        self.logger.info("sample:{}".format(self.option("sample")))
        self.logger.info("region:{}".format(self.option("depth")))
        self.logger.info("depth:{}".format(self.result_file))
        self.logger.info("depth_comare:{}".format(self.option("depth_comare")))
        self.snp_search = self.add_tool('medical_transcriptome.snp.snp_search')
        self.snp_search.set_options({
                'snp_result_file': self.result_file,
                'target_genes': self.final_target_genes,
                'type':self.option("type"),
                'sample': self.option("sample"),
                'region': self.option("region"),
                'depth': self.option("depth"),
                'depth_comare': self.option("depth_comare"),
            })

        self.snp_search.on('end', self.set_db)
        self.snp_search.run()

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path


    def end(self):
        if self.option("type").lower() == "snp":
            os.rename(os.path.join(self.snp_search.output_dir, 'target_detail.txt'), os.path.join(self.snp_search.output_dir, 'snp_detail_anno.xls'))
        else:
            os.rename(os.path.join(self.snp_search.output_dir, 'target_detail.txt'), os.path.join(self.snp_search.output_dir, 'indel_detail_anno.xls'))
        if os.path.exists(os.path.join(self.snp_search.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.snp_search.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.snp_search.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.snp_search.output_dir)
        self.inter_dirs = [
            ["03 Gene_structure_analysis", "", "基因结构分析数据挖掘结果目录", 0],
            ["03 Gene_structure_analysis/02_SNP_InDel_Analysis", "", "SNP/InDel分析结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "", 0],
            ['indel_detail_anno.xls', 'xls', 'InDel分析结果注释搜索详情表'],
            ['snp_ detail_anno.xls', 'xls', 'SNP分析结果注释搜索详情表'],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(SnpSearchWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        snp = self.api.api("medical_transcriptome.snp")
        # add result info
        snp_search = os.path.join(self.snp_search.output_dir, 'target_detail.txt')
        snp.add_snp_search_result(snp_search, self.option('main_id'))
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.snp_search import SnpSearchWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "SNP_search" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.report.snp_search",
            "options": dict(
                # result_path='s3://refrnav2/files/m_188/188_5d01dede4f911/tsg_38314/workflow_results/08SNP/',
                result_path="s3://refrnav2/files/m_188/188_5d01dede4f911/tsg_38314/workflow_results/08SNP",
                # target_genes='',
                target_genes='ENSG00000004455,ENSG00000007968',
                geneset_snp="/mnt/ilustre/users/sanger-dev/workspace/20200813/GenesetPpi_tsg_38314_2767_2975/genes",
                sample = "H1581_1",
                type ="snp",
                anno = "/mnt/ilustre/users/sanger-dev/workspace/20200813/GenesetPpi_tsg_38314_2767_2975/seq_annot.xls",
                id_type = "gene_id",
                region='exonic',
                depth_comare="greateroreq",
                depth="1"

            )
        }

        wsheet = Sheet(data=data)
        wf =SnpSearchWorkflow(wsheet)
        wf.sheet.id = 'snp_search'
        wf.sheet.project_sn = 'snp_search'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
