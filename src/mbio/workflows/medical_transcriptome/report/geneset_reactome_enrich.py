# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

from biocluster.workflow import Workflow
import os
import shutil
import glob
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import pandas as pd
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir

class GenesetReactomeEnrichWorkflow(Workflow):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetReactomeEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'task_id', 'type': 'string'},
            {'name': 'geneset_reactome', 'type': 'string'},
            {'name': 'reactome_annot', 'type': 'string'},
            {'name': 'geneset_list', 'type': 'string'},
            {'name': 'geneset_names', 'type': 'string'},
            {'name': 'multi_geneset_list', 'type': 'string'},
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
            {'name': 'reactome_version', 'type': 'string', 'default': None},
            {'name': 'source', 'type': 'string'},
            {'name': 'result', 'type': 'infile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/03 Enrich/03 Reactome')

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
        super(GenesetReactomeEnrichWorkflow, self).send_log(data)

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
        super(GenesetReactomeEnrichWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_geneset_reactome_enrich", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        self.run_reactome_rich()

    def run_reactome_rich(self):
        self.step.add_steps('reactome_rich')
        self.enrich_tool = self.add_tool('medical_transcriptome.geneset.reactome_enrich')
        self.enrich_tool.set_options({
            'reactome_list': self.option('reactome_annot'),
            'diff_list': self.option('geneset_list'),
            'method': self.option('method'),
        })
        self.enrich_tool.on('start', self.set_step, {'start': getattr(self.step, 'reactome_rich')})
        self.enrich_tool.on('end', self.set_step, {'start': getattr(self.step, 'reactome_rich')})
        self.enrich_tool.on('end', self.run_reactome_class)
        self.enrich_tool.run()

    def run_reactome_class(self):
        self.step.add_steps('reactome_class')
        self.reactome_class = self.add_tool('medical_transcriptome.geneset.reactome_class')
        self.reactome_class.set_options({
            'geneset_ids': self.option('multi_geneset_list'),
            'reactome_annot': self.option('reactome_annot'),
            'reactome_version': self.option('reactome_version')
        })
        self.reactome_class.on('start', self.set_step, {'start': getattr(self.step, 'reactome_class')})
        self.reactome_class.on('end', self.set_step, {'start': getattr(self.step, 'reactome_class')})
        self.reactome_class.on('end', self.set_output)
        self.reactome_class.run()

    def set_output(self):
        self.replace_reactome_link()
        self.set_db()

    def replace_reactome_link(self):
        '''
        替换reactome链接
        '''
        enrich_result = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        enrich_result_out = self.output_dir + '/Reactome_enrich_stat.xls'
        annot_result = self.reactome_class.output_dir + '/reactome_path.xls'
        annot_df = pd.read_table(annot_result, sep='\t',header=0)
        map2link = dict(zip(annot_df['Pathway ID'], annot_df['link']))
        map2des = dict(zip(annot_df['Pathway ID'], annot_df['Description']))
        map2category = dict(zip(annot_df['Pathway ID'], annot_df['category']))

        enrich_df = pd.read_table(enrich_result, sep='\t',header=0)
        enrich_df['Description'] = [map2des.get(x, "") for x in enrich_df['Pathway ID']]
        enrich_df['link'] = [map2link.get(x, "") for x in enrich_df['Pathway ID']]
        enrich_df['category'] = [map2category.get(x, "") for x in enrich_df['Pathway ID']]
        enrich_df.to_csv(enrich_result_out, sep = '\t', index=False)




    def set_db(self):
        api_geneset = self.api.api('medical_transcriptome.geneset_annot')
        if self.option("source") == "diff_exp":
            api_geneset = self.api.api('medical_transcriptome.diff_geneset_annot')

        enrich_result = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        enrich_result_out = self.output_dir + '/Reactome_enrich_stat.xls'
        results_num = len(open(enrich_result_out, 'r').readlines())
        if results_num>1:
            self.has_results = True
        if self.has_results:
            svg_path = self.reactome_class.output_dir + '/svg'
            api_geneset.add_geneset_reactome_enrich(self.option("main_table_id"),
                                            self.option("geneset_names"),
                                            enrich_result_out,
                                            svg_path)

            add_name = GeneInfoSupple(id2namedict = None, task_id = self.option("task_id"), level ="G")
            add_columns = ["Genes"]
            result = os.path.join(self.work_dir, "Reactome_enrich_stat.xls")
            try:
                add_name.add_gene_name(enrich_result_out, split="|", annot_type="reactome", add_columns=add_columns, outfile_path=result)
                CopyFile().linkfile(result, enrich_result_out)
            except Exception as e:
                self.logger.info("添加基因name 失败 {}".format(e))
            if os.path.exists(os.path.join(self.output_dir, "svg.tar.gz")):
                os.remove(os.path.join(self.output_dir, "svg.tar.gz"))
            os.link(os.path.join(self.reactome_class.output_dir, "svg.tar.gz"),
                    os.path.join(self.output_dir, "Reactome_pathways.tar.gz"))
        else:
            api_geneset.update_db_record('sg_geneset_reactome_enrich', self.option('main_table_id'),
                                         has_results=False)



        # CopyFile().linkdir(svg_path, self.output_dir + "/Reactome_pathways")
        # mark_file = glob.glob(self.output_dir + "/Reactome_pathways/*.mark")
        # for file in mark_file:
        #     os.remove(file)
        self.end()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        enrich_result_out = self.output_dir + '/Reactome_enrich_stat.xls'
        if self.has_results:
            chart.chart_geneset_enrich_reactome(reactome_enrich_table=enrich_result_out)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = glob.glob(self.work_dir + "/*.pdf")
            for p in pdf_file:
                linkfile(p, self.output_dir + "/" + os.path.basename(p))
        else:
            pass

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/03 Enrich", "", "基因集功能富集", 0],
            ["04 GeneSet/03 Enrich/03 Reactome", "", "基因集Reactome富集分析分析", 0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '基因集Reactome富集分析文件',0],
            ['Reactome_enrich_stat.xls', 'xls', 'Reactome富集分析结果表 ',0],
            ["Reactome_pathways.tar.gz", "", "基因集Reactome通路图",0],
            ["*bar.pdf", "pdf", "Reactome富集分析柱形图", 0],
            ["*bar_line.pdf", "pdf", "Reactome富集分析柱形图(带折线)", 0],
            ["*buble.pdf", "pdf", "Reactome富集分析气泡图", 0],
            ["*buble2.pdf", "pdf", "Reactome富集分析气泡图(分散性)", 0],
            ["Reactome_pathways", "", "基因集Reactome通路图",0],
            ['Reactome_pathways/*.svg', '', 'Reactome通路图片svg',0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])


        super(GenesetReactomeEnrichWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
