# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

from biocluster.workflow import Workflow
import os
import shutil
import glob
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir

class GenesetDoEnrichWorkflow(Workflow):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetDoEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'task_id', 'type': 'string'},
            {'name': 'geneset_kegg', 'type': 'string'},
            {'name': 'kegg_table', 'type': 'string'},
            {'name': 'do_list', 'type': 'string'},
            {'name': 'geneset_list', 'type': 'string'},
            {'name': 'all_list', 'type': 'string'},
            {'name': 'anno_type', 'type': 'string'},
            {'name': 'geneset_type', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_table_id', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'submit_location', 'type': 'string'},
            {'name': 'task_type', 'type': 'string'},
            {'name': 'method', 'type': 'string'},
            {'name': 'type', 'type': 'string'},
            {'name': 'add_info', 'type': 'string', 'default': None},  # 输入两列的列表文件，有head，第一列为pathway，第二列为底图链接
            {'name': 'geneset_id', 'type': 'string'},
            {'name': 'kegg_table_2', 'type': 'string'},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {'name': 'source', 'type': 'string'},
            {'name': 'result', 'type': 'infile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/03 Enrich/04 DO')
        self.inter_dirs = []
        self.has_results = False
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="medical_transcriptome",
                                                  main_id=self.option('main_table_id'))
            interactiondelete.delete_interactions_records()

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
        super(GenesetDoEnrichWorkflow, self).send_log(data)

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        self.get_run_log()
        self.run_tool()
        super(GenesetDoEnrichWorkflow, self).run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        geneset_do_table = self.enrich_tool.option('result').path
        if self.has_results:
            chart.chart_geneset_enrich_do(do_enrich_table=geneset_do_table)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = glob.glob(self.work_dir + "/*.pdf")
            for p in pdf_file:
                linkfile(p, self.output_dir + "/" + os.path.basename(p))

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_geneset_do_enrich", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        self.run_do_enrich()

    def run_do_enrich(self):
        self.step.add_steps('go_enrich')
        self.enrich_tool = self.add_tool('medical_transcriptome.geneset.do_enrich')
        self.enrich_tool.set_options({
            'diff_list': self.option('geneset_list'),
            'do_list': self.option('do_list'),
            'method': self.option('method')
        })
        self.enrich_tool.on('start', self.set_step, {'start': getattr(self.step, 'go_enrich')})
        self.enrich_tool.on('end', self.set_step, {'end': getattr(self.step, 'go_enrich')})
        self.enrich_tool.on('end', self.set_output)
        self.enrich_tool.run()


    def set_output(self):
        add_name = GeneInfoSupple(id2namedict = None, task_id = self.option("task_id"), level ="G")
        add_columns = ["Genes"]
        result = os.path.join(self.output_dir, "do_enrich_stat.xls")
        results_num = len(open(self.enrich_tool.option('result').path, 'r').readlines())
        if results_num > 1:
            self.has_results = True
            try:
                add_name.add_gene_name(self.enrich_tool.option('result').path, split="|", annot_type="do", add_columns=add_columns, outfile_path=result)
            except Exception as e:
                self.logger.info("添加基因name 失败 {}".format(e))
                # CopyFile().linkfile(output_file, result)
            self.option('result').set_path(result)
        self.set_db()


    def set_db(self):



        api = self.api.api('medical_transcriptome.geneset_annot')
        if self.option("source") == "diff_exp":
            api = self.api.api('medical_transcriptome.diff_geneset_annot')
        geneset_do_table = self.enrich_tool.option('result').path
        if self.has_results :
            api.add_geneset_do_enrich(geneset_do_table, self.option("main_id"))
        else:
            api.update_db_record('sg_geneset_do_enrich', self.option('main_id'),
                                         has_results=False)

        self.end()

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/03 Enrich", "", "基因集功能富集", 0],
            ["04 GeneSet/03 Enrich/04 DO", "", "基因集DO功能富集", 0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '基因集DO富集分析文件',0],
            ["*bar.pdf", "pdf", "DO富集分析柱形图", 0],
            ["*bar_line.pdf", "pdf", "DO富集分析柱形图(带折线)", 0],
            ["*buble.pdf", "pdf", "DO富集分析气泡图", 0],
            ["*buble2.pdf", "pdf", "DO富集分析气泡图(分散性)", 0],
            ['*.xls', 'xls', '基因集DO富集分析统计表 ',0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0]
        ])

        super(GenesetDoEnrichWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
