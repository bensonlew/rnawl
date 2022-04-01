#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/30 10:25
@file    : geneset_gsva.py
"""
import unittest
from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import re
from bson.objectid import ObjectId
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import shutil
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset


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
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {'name': 'sample_list_str', 'type': 'string'},

            {"name": "geneset_source", "type": "string"},
            {"name": "geneset_id", "type": "string"},

            {"name": "group", "type": "string"},
            {'name': 'cmp', 'type': 'string'},
            {'name': 'id2name', 'type': 'string'},

            {"name": "max_num", "type": "string", "default": "500"},
            {"name": "min_num", "type": "string", "default": "15"},

            {"name": "level", "type": "string", "default": "G"},
            {"name": "go_genesets", "type": "string", "default": "all"},
            {"name": "go_list", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "go_type", "type": "string", "default": ""},

            {'name': 'stat_type', 'type': 'string', 'default': 'padjust'},
            {'name': 'fc', 'type': 'float', 'default': 2},
            {'name': 'pvalue', 'type': 'float', 'default': 0.05},

            {"name": "gmx", "type": "infile", "format": "ref_rna_v2.geneset_gmt"},
            {"name": "matrix", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "genes_detail", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'es_method', 'type': 'string', 'default': 'diff'},

            dict(name="scm", type="string", default='complete'),
            dict(name="scd", type="string", default='euclidean'),
            dict(name="sct", type="string", default='hierarchy'),
            dict(name="gct", type="string", default='hierarchy'),
            dict(name="gcm", type="string", default='complete'),
            dict(name="gcd", type="string", default='euclidean'),

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 Advanced_Analysis/04 GSVA_Analysis')
        # self.inter_dirs = []
        self.step.add_steps('gsva_analysis')
        if self.option("go_list").is_set:
            self.go_gene_gmx = self.add_tool('medical_transcriptome.geneset.go_gene_gmx')
        self.gsva_tool = self.add_tool('medical_transcriptome.geneset.gsva')
        self.gsva_tool.on('end', self.run_cluster_heatmap)

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
        super(GenesetGsvaWorkflow, self).send_log(data)

    def get_go_gmx(self):
        opts = {
            "go_list": self.option("go_list").prop['path'],
            "min_num": self.option("min_num"),
            "max_num": self.option("max_num"),
            "go_type": self.option("go_type"),
            "go_sets": self.option("go_genesets")
        }
        self.go_gene_gmx.set_options(opts)
        self.go_gene_gmx.run()


    def gsva_runner(self):
        # 输出: biomart.xls, biomart.json
        if self.option("go_list").is_set:
            self.gmx = self.go_gene_gmx.output_dir + '/go.gmt'
        else:
            self.gmx =self.option('gmx').path
        self.gsva_tool.set_options({
            'max_num': int(self.option('max_num')),
            'min_num': int(self.option('min_num')),
            'gmx': self.gmx,
            'matrix': self.option('matrix'),
            'genes_detail': self.option('genes_detail'),
            'es_method': self.option('es_method'),
            'geneset_source': self.option('geneset_source')
        })
        self.gsva_tool.on('start', self.set_step, {'start': self.step.gsva_analysis})
        self.gsva_tool.on('end', self.set_step, {'end': self.step.gsva_analysis})
        self.gsva_tool.run()

    def run_cluster_heatmap(self):
        self.cluster_heatmap = self.add_tool('medical_transcriptome.geneset.exp_cluster2gsva')
        options = dict(
            exp=self.gsva_tool.option('es_exp_table').path,
            group=self.option('group'),
            sct=self.option('sct'),
            gct=self.option('gct'),
            scm=self.option('scm'),
            gcm=self.option('gcm'),
            scd=self.option('scd'),
            gcd=self.option('gcd'),
            # use_group=self.option('use_group'),
        )
        self.cluster_heatmap.set_options(options)
        self.cluster_heatmap.on('start', self.set_step, {'start': self.step.gsva_analysis})
        self.cluster_heatmap.on('end', self.set_step, {'end': self.step.gsva_analysis})
        self.cluster_heatmap.on('end', self.run_diff)
        self.cluster_heatmap.run()


    def run_diff(self):
        self.diff = self.add_tool('medical_transcriptome.geneset.diff_geneset')
        options = dict(
            count=self.gsva_tool.option('es_exp_table').path,
            group=self.option('group'),
            cmp=self.option('cmp'),
            fc=self.option('fc'),
            pvalue_padjust=self.option('stat_type'),
            pvalue=self.option('pvalue')
        )
        self.diff.set_options(options)
        self.diff.on('start', self.set_step, {'start': self.step.gsva_analysis})
        self.diff.on('end', self.set_step, {'end': self.step.gsva_analysis})
        self.diff.on('end', self.set_output)
        self.diff.run()

    def run(self):
        self.get_run_log()
        if self.option("go_list").is_set:
            self.go_gene_gmx.on('end', self.gsva_runner)
            self.get_go_gmx()
        else:
            self.gsva_runner()
        super(GenesetGsvaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_geneset_gsva", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

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
        geneset_description = dict()
        with open(self.gmx, 'r') as g:
            lines = g.readlines()
            for line in lines:
                geneset, description = line.strip().split('\t')[0:2]
                geneset_description[geneset] = description
        api_gsva = self.api.api('medical_transcriptome.gsva')
        # main_id = api_gsva.add_gsva(self.option("main_table_id"), self.gsva_tool.option('es_exp_table').path, self.option('gmx').path, self.gsva_tool.option('gsva').path, self.option('id2name'), geneset_description)
        main_id = api_gsva.add_gsva(self.option("main_table_id"), self.gsva_tool.option('es_exp_table').path,
                                    self.gmx, self.gsva_tool.option('gsva').path, self.option('id2name'),
                                    geneset_description)
        if self.option("go_list").is_set:
            self.modify_go_cluster_dir(self.cluster_heatmap.output_dir)
        api_gsva.add_gsva_cluster(self.cluster_heatmap.output_dir, main_id)
        api_gsva.add_gsva_diff(self.diff.output_dir, main_id, self.option('group'))
        self.end()

    def modify_go_cluster_dir(self,cluster_heatmap_output):
        results = os.listdir(cluster_heatmap_output)
        for result in results:
            file_path = os.path.join(cluster_heatmap_output,result)
            with open(file_path,"r") as r:
                info = r.read()
                final_info = info.replace("-",":")
            with open(file_path,"w") as w:
                w.write(final_info)


    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        group_dict = dict()
        with open(self.option('group'), 'r') as g:
            for line in g.readlines():
                if line.startswith('#'):
                    continue
                else:
                    sample, group = line.strip().split('\t')
                    if group not in group_dict:
                        group_dict[group] = [sample]
                    else:
                        group_dict[group].append(sample)
        cluster_exp = self.gsva_tool.output_dir + "/gsva_table.txt"
        cluster_tree = self.cluster_heatmap.output_dir + "/seq.cluster_tree.txt"
        sample_tree = self.cluster_heatmap.output_dir + "/sample.cluster_tree.txt"
        chart.chart_geneset_gsva_cluster(cluster_exp, cluster_tree, sample_tree, group_dict=group_dict)
        chart.to_pdf()

        limma_file = glob.glob(self.diff.output_dir + '/*.xls')
        for i in limma_file:
            chart.chart_geneset_gsva_diff(i)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            print "copy", p, self.cluster_heatmap.output_dir + "/" + os.path.basename(p)
            shutil.copyfile(p, self.output_dir + "/" + os.path.basename(p))


    def set_output(self):
        gsva_output_file = glob.glob(os.path.join(self.gsva_tool.output_dir, '*'))
        for i in gsva_output_file:
            os.link(i, os.path.join(self.output_dir, os.path.basename(i)))
        diff_output_dir = glob.glob(os.path.join(self.diff.output_dir, '*'))
        for d in diff_output_dir:
            os.link(d, os.path.join(self.output_dir, os.path.basename(d)))
        cluster_output_dir = glob.glob(os.path.join(self.cluster_heatmap.output_dir, '*'))
        for c in cluster_output_dir:
            os.link(c, os.path.join(self.output_dir, os.path.basename(c)))
        # shutil.rmtree(self.output_dir)
        # shutil.copytree(self.gsva_tool.output_dir, self.output_dir)
        # shutil.copytree(self.diff.output_dir, self.output_dir)
        # shutil.copytree(self.cluster_heatmap.output_dir, self.output_dir)
        # # os.rename(os.path.join(self.output_dir, 'result.tsv'), os.path.join(self.output_dir, 'result.xls'))
        # # self.option('result').set_path(os.path.join(self.output_dir, 'result.xls'))
        self.set_db()


    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 Advanced_Analysis", "", "高级分析结果目录",0],
            ["05 Advanced_Analysis/04 GSVA_Analysis", "", "GSVA分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "GSVA分析结果目录",0],
            ["./gsva_table.txt", "", "所有基因集gsva结果表",0],
            ["./sample.cluster_tree.txt", "", "样本聚类树文件",0],
            ["./seq.cluster_tree.txt", "", "基因集聚类树文件",0],
            ["./*.limma.xls", "xls", "基因集差异分析结果文件",0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ['./*.heat_corr.pdf', "pdf", "基因集聚类热图", 0],
            ['./geneset_gsva_diff*.pdf', "pdf", "基因集差异柱状图", 0]
        ])
        super(GenesetGsvaWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.geneset_gsva import GenesetGsvaWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "GSVA" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.report.geneset_gsva",
            "options": dict(
                gmx='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/medical_transcriptome/GSVA/homo/gsea.gmt',
                matrix='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/medical_transcriptome/GSVA/homo/gene.tpm.matrix.xls',
                genes_detail='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/medical_transcriptome/GSVA/homo/gsea.chip',
                geneset_source='msigdb',
                fc=2,
                stat_type='padjust',
                pvalue=0.05,
                group='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/medical_transcriptome/GSVA/homo/example_group_1528169151.txt',
                cmp='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/medical_transcriptome/GSVA/homo/example_control_1528169151.txt',
                sct="hierarchy",
                gct="hierarchy",
                scm="complete",
                gcm="complete",
                scd="euclidean",
                gcd="euclidean",
                # sct="no",
                # gct="no",

            )
        }

        wsheet = Sheet(data=data)
        wf =GenesetGsvaWorkflow(wsheet)
        wf.sheet.id = 'Gsva'
        wf.sheet.project_sn = 'Gsva'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)