# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

from biocluster.workflow import Workflow
import os
import shutil
import glob
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset
from biocluster.config import Config
from bson import ObjectId
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir


class GenesetEnrichWorkflow(Workflow):
    '''
    last_modify: 2019.06.12
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'task_id', 'type': 'string'},
            {'name': 'geneset_kegg', 'type': 'string'},
            {'name': 'kegg_table', 'type': 'string'},
            {'name': 'go_list', 'type': 'string'},
            {'name': 'go_version', 'type': 'string', 'default': '20200628'},
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
            # {'name': 'type', 'type': 'string'},
            {'name': 'add_info', 'type': 'string', 'default': None},  # 输入两列的列表文件，有head，第一列为pathway，第二列为底图链接
            {'name': 'geneset_id', 'type': 'string'},
            {'name': 'kegg_table_2', 'type': 'string'},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {'name': 'source', 'type': 'string'},
            {'name': 'result', 'type': 'infile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "kegg":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/03 Enrich/02 KEGG')
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/03 Enrich/01 GO')
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
        super(GenesetEnrichWorkflow, self).send_log(data)

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
        # if self.option('anno_type') == 'kegg':
        #     options = {
        #         'kegg_table': self.option('kegg_table'),
        #         # 'all_list': background_path,
        #         'diff_list': self.option('geneset_list'),
        #         'correct': self.option('method'),
        #         'add_info': self.option('add_info')
        #     }
        # else:
        #     options = {
        #         'diff_list': self.option('geneset_list'),
        #         # 'all_list': background_path,
        #         'go_list': self.option('go_list'),
        #         # 'pval': self.option('pval'),
        #         'method': self.option('method'),
        #     }
        # self.logger.info(options)
        # self.enrich_tool.set_options(options)
        # if self.option('anno_type') == 'kegg':
        #     self.enrich_tool.on('end', self.run_kegg_class)
        #     self.kegg_class.on('end', self.set_db)
        # else:
        #     self.enrich_tool.on('end', self.set_db)
        # self.enrich_tool.run()
        self.get_run_log()
        self.run_tool()
        super(GenesetEnrichWorkflow, self).run()

    def get_run_log(self):
        if self.option('anno_type') == 'go':
            table = "sg_geneset_go_enrich"
        elif self.option('anno_type') == 'kegg':
            table = "sg_geneset_kegg_enrich"
        get_run_log = GetRunLog("medical_transcriptome", table=table, main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        if self.option('anno_type') == 'go':
            self.run_go_enrich()
        elif self.option('anno_type') == 'kegg':
            self.run_kegg_rich()

    def run_go_enrich(self):
        self.step.add_steps('go_enrich')
        self.enrich_tool = self.add_tool('medical_transcriptome.geneset.go_enrich')
        self.enrich_tool.set_options({
            'diff_list': self.option('geneset_list'),
            'go_list': self.option('go_list'),
            'go_version': self.option('go_version'),
            'method': self.option('method')
        })
        self.enrich_tool.on('start', self.set_step, {'start': getattr(self.step, 'go_enrich')})
        self.enrich_tool.on('end', self.set_step, {'end': getattr(self.step, 'go_enrich')})
        self.enrich_tool.on('end', self.set_output)
        self.enrich_tool.run()

    def run_kegg_rich(self):
        self.step.add_steps('kegg_rich')
        self.enrich_tool = self.add_tool('medical_transcriptome.geneset.kegg_rich')
        self.enrich_tool.set_options({
            'kegg_table': self.option('kegg_table'),
            'diff_list': self.option('geneset_list'),
            'kegg_table2': self.option('kegg_table_2'),
            'correct': self.option('method'),
            'add_info': self.option('add_info'),
            'kegg_version': self.option('kegg_version'),
        })
        self.enrich_tool.on('start', self.set_step, {'start': getattr(self.step, 'kegg_rich')})
        self.enrich_tool.on('end', self.set_step, {'start': getattr(self.step, 'kegg_rich')})
        self.enrich_tool.on('end', self.run_kegg_class)
        self.enrich_tool.run()

    def run_kegg_class(self):
        self.step.add_steps('kegg_class')
        self.kegg_class = self.add_tool('medical_transcriptome.geneset.kegg_class')
        self.kegg_class.set_options({
            'geneset_kegg': self.option('geneset_kegg'),
            'kegg_table': self.option('kegg_table'),
            'kegg_table2': self.option('kegg_table_2'),
            'geneset_id': self.option('geneset_id'),
            'background_links': self.option('add_info'),
            # 'type': self.option('type'),
            'task_id': self.option('task_id'),
            "kegg_version": self.option('kegg_version'),
            'source': self.option('source')
        })
        self.kegg_class.on('start', self.set_step, {'start': getattr(self.step, 'kegg_class')})
        self.kegg_class.on('end', self.set_step, {'start': getattr(self.step, 'kegg_class')})
        self.kegg_class.on('end', self.set_output)
        self.kegg_class.run()

    def set_output(self):
        if self.option('anno_type') == 'go':
            go_enrich_result_file = os.path.join(self.enrich_tool.output_dir, 'go_enrich_geneset_list_gene.xls')
            if os.path.exists(go_enrich_result_file):
                self.has_results = True
            if self.has_results:
                result = os.path.join(self.output_dir, os.path.basename(self.enrich_tool.option('result').path))
                shutil.copy(self.enrich_tool.option('result').path, result)
                self.option('result').set_path(result)
            else:
                pass
        self.set_db()

    def replace_kegg_link(self):
        '''
        替换kegg链接
        '''
        enrich_result = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        annot_result = self.kegg_class.output_dir + '/kegg_stat.xls'
        with open(annot_result, 'rb') as f:
            map2link = {line.split("\t")[0]:line.strip().split("\t")[-1] for line in f.readlines()[1:] if line.split("\t")[-1].startswith("http")}
        if os.path.exists(enrich_result + "relink"):
            os.remove(enrich_result + "relink")
        with open(enrich_result, 'rb') as f, open(enrich_result + "relink", 'w') as fo:
            header = f.readline()
            fo.write(header)
            for line in f:
                cols = line.split("\t")
                if cols[3] in map2link:
                    cols[9] = map2link[cols[3]]
                fo.write("\t".join(cols))
        os.remove(enrich_result)
        os.link(enrich_result + "relink", enrich_result)



    def set_db(self):
        if self.option('anno_type') == 'go':
            api = self.api.api('medical_transcriptome.go_enrich')
            if self.has_results:
                api.add_go_enrich(self.enrich_tool.option('result').path, self.option('main_id'), self._sheet.output)
                split_str = ';'
                column_name_judge_str = '_list'
                needed_add_genename_file = self.output_dir + "/go_enrich_geneset_list_gene.xls"
                add_name = GeneInfoSupple(task_id=self.option("task_id"), level=self.option("geneset_type"))
                try:
                    geneid_list = [i for i in open(needed_add_genename_file).readline().strip().split('\t') if
                                   column_name_judge_str in i]
                    add_name.add_gene_name(needed_add_genename_file, split_str, self.option('anno_type'), geneid_list,
                                           needed_add_genename_file)
                    print
                    "geneid_list:    " + str(geneid_list)
                except Exception as e:
                    self.logger.info("添加基因name 失败 {}".format(e))
            else:
                # conn = api.db["sg_geneset_kegg_class"]
                # conn.update({"_id": ObjectId(self.option('main_table_id'))}, {"$set": {'has_results': False, 'status': 'end'}}, upsert=True)
                api.update_db_record('sg_geneset_go_enrich', self.option('main_table_id'),
                                     has_results=False)
        elif self.option('anno_type') == 'kegg':
            if self.option("source") == "diff_exp":
                self.replace_kegg_link()
            self.output_dir1 = self.enrich_tool.output_dir
            self.output_dir2 = self.kegg_class.output_dir

            api_geneset = self.api.api('medical_transcriptome.medical_transcriptome_geneset')
            output_file = glob.glob("{}/*.xls".format(self.output_dir1))
            workflow_output = glob.glob("{}/*.xls".format(self.output_dir1))
            workflow_output = self.get_workflow_output_dir() + '/' + workflow_output[0].split('/')[-1]
            results_num = len(open(output_file[0], 'r').readlines())
            if results_num > 1:
                self.has_results = True
                api_geneset.add_kegg_enrich_detail(self.option('main_table_id'), output_file[0],
                                                   self.option('geneset_list'), self.option('all_list'))
                api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'),
                                             result_dir=workflow_output)
                api_geneset.add_kegg_enrich_pic(self.option('main_table_id'), output_file[0],
                                                self.output_dir2 + '/pathways', source=self.option("source"))
                graph_dir = os.path.join(self.get_workflow_output_dir(), 'pathways')
                api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'),
                                             graph_dir=graph_dir)

                ## 文件加入对应于基因id列的【基因名】这一新列
                split_str = '|'
                column_name_judge_str = 'Genes'
                needed_add_genename_file = os.path.join(self.output_dir1,
                                                        "geneset_list_gene.list.DE.list.check.kegg_enrichment.xls")
                add_name = GeneInfoSupple(task_id=self.option("task_id"), level=self.option("geneset_type"))
                try:
                    geneid_list = [i for i in open(needed_add_genename_file).readline().strip().split('\t') if
                                   column_name_judge_str in i]
                    add_name.add_gene_name(needed_add_genename_file, split_str, self.option('anno_type'), geneid_list,
                                           needed_add_genename_file)
                    print
                    "geneid_list:    " + str(geneid_list)
                except Exception as e:
                    self.logger.info("添加基因name 失败 {}".format(e))
            else:
                api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'),
                                             has_results=False)

        # if self.option('anno_type') == 'go':
        #     self.output_dir1 = self.enrich_tool.output_dir
        # elif self.option('anno_type') == 'kegg':
        #     self.output_dir2 = self.kegg_class.output_dir
        # api_geneset = self.api.api('medical_transcriptome.medical_transcriptome_geneset')
        # output_file = glob.glob('{}/*.xls'.format(self.output_dir1))
        # workflow_output = glob.glob('{}/*.xls'.format(self.output_dir1))
        # workflow_output = self.get_workflow_output_dir() + '/' + workflow_output[0].split('/')[-1]
        # png_file = glob.glob('{}/*.png'.format(self.output_dir))
        # go_png = self.output_dir + '/go_lineage.png'
        # go_pdf = self.output_dir + '/go_lineage.pdf'
        # go_adjust_png = self.output_dir1 + '/adjust_lineage.png'
        # go_adjust_pdf = self.output_dir1 + '/adjust_lineage.pdf'
        # if self.option('anno_type') == 'kegg':
        #     api_geneset.add_kegg_enrich_detail(self.option('main_table_id'), output_file[0], self.option('geneset_list'), self.option('all_list'))
        #     api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'), result_dir= workflow_output)
        #     api_geneset.add_kegg_enrich_pic(self.option('main_table_id'), output_file[0], self.output_dir2 + '/pathways')
        #     graph_dir = os.path.join(self.get_workflow_output_dir(), 'pathways')
        #     api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'), graph_dir=graph_dir)
        # else:
        #     api_geneset.add_go_enrich_detail(self.option('main_table_id'), output_file[0])
            # api_geneset.update_directed_graph(self.option('main_table_id'), go_adjust_png, go_adjust_pdf)
            # api_geneset.update_db_record('sg_geneset_go_enrich', self.option('main_table_id'), result_dir=workflow_output)
        self.end()

    def get_geneset_name(self):
        name = "geneset"
        try:
            db = Config().get_mongo_client(mtype="medical_transcriptome")[Config().get_mongo_dbname("medical_transcriptome")]
            geneset_info = db["sg_geneset"].find_one({"_id": ObjectId(self.option("geneset_id"))})
            name = geneset_info["name"]
        except:
            pass
        return name

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"

        if self.option('anno_type') == 'go':
            if self.has_results:
                go_enrich_table = self.enrich_tool.option('result').path
                with open(go_enrich_table, "r") as f:
                    lines = f.readlines()
                    if len(lines) >= 2:
                        chart.chart_geneset_enrich_go(go_enrich_table, geneset_name=self.get_geneset_name())
                        chart.to_pdf()
                        # move pdf to result dir
                        pdf_file = glob.glob(self.work_dir + "/*.pdf")
                        for p in pdf_file:
                            linkfile(p, self.output_dir + "/" + os.path.basename(p))

        elif self.option('anno_type') == 'kegg':
            if self.has_results:
                kegg_enrich_table = os.path.join(self.output_dir1,
                                                 "geneset_list_gene.list.DE.list.check.kegg_enrichment.xls")
                with open(kegg_enrich_table, "r") as f:
                    lines = f.readlines()
                    if len(lines) >= 2:
                        try:
                            chart.chart_geneset_enrich_kegg(kegg_enrich_table, geneset_name=self.get_geneset_name())
                            chart.to_pdf()
                            # move pdf to result dir
                            pdf_file = glob.glob(self.work_dir + "/*.pdf")
                            for p in pdf_file:
                                linkfile(p, self.output_dir + "/" + os.path.basename(p))
                        except:
                            pass
    def end(self):
        self.chart()
        if self.option('anno_type') == 'go':
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录", 0],
                ["04 GeneSet/03 Enrich", "", "基因集功能富集", 0],
                ["04 GeneSet/03 Enrich/01 GO", "", "基因集GO功能富集", 0]
            ]
            result_dir.add_relpath_rules([
                ['.', '', '基因集GO富集分析文件',0],
                ['go_enrich_geneset_list_gene.xls', 'xls','基因集GO富集分析统计表 ',0],
                ["*bar.pdf", "pdf", "GO富集分析柱形图",0],
                ["*bar_line.pdf", "pdf", "GO富集分析柱形图(带折线)", 0],
                ["*buble.pdf", "pdf", "GO富集分析气泡图", 0],
                ["*buble2.pdf", "pdf", "GO富集分析气泡图(分散型)", 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0]
            ])
        elif self.option('anno_type') == 'kegg':
            # result_dir1 = self.add_upload_dir(self.output_dir1)
            # result_dir1.add_relpath_rules([
            #     ['.', '', '基因集KEGG功能富集', 0],
            #     ['kegg_enrich_geneset_list_gene.xls', '', 'KEGG富集分析统计表', 0],
            # ])
            linkfile(os.path.join(self.output_dir1, "geneset_list_gene.list.DE.list.check.kegg_enrichment.xls"),
                     os.path.join(self.output_dir, "kegg_enrich_geneset_list_gene.xls"))
            linkdir(os.path.join(self.output_dir2, "pathways"), os.path.join(self.output_dir, "pathways"))
            # shutil.rmtree(self.output_dir2 + "/ko")
            rm_file = glob.glob(self.output_dir + "/pathways/*.html.mark")
            for file in rm_file:
                os.remove(file)
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir2 = self.add_upload_dir(self.output_dir)

            # os.link(os.path.join(self.output_dir1, "geneset_list_gene.list.DE.list.check.kegg_enrichment.xls"),
            #         os.path.join(self.output_dir2, "kegg_enrich_geneset_list_gene.xls"))
            # shutil.rmtree(self.output_dir2 + "/ko")
            # if os.path.exists(self.output_dir2 + "/pathways"):
            #     shutil.rmtree(self.output_dir2 + "/pathways")
            # rm_file = glob.glob(self.output_dir2 + "/pathways/*.html.mark") + glob.glob(self.output_dir2 + "/kegg_stat.xls")
            # for file in rm_file:
            #     os.remove(file)
            # if os.path.exists(os.path.join(self.output_dir2, os.path.basename(self.run_log))):
            #     os.remove(os.path.join(self.output_dir2, os.path.basename(self.run_log)))
            # os.link(self.run_log, os.path.join(self.output_dir2, os.path.basename(self.run_log)))
            # result_dir2 = self.add_upload_dir(self.output_dir2)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录", 0],
                ["04 GeneSet/03 Enrich", "", "基因集功能富集", 0],
                ["04 GeneSet/03 Enrich/02 KEGG", "", "基因集KEGG功能富集", 0]
            ]
            result_dir2.add_relpath_rules([
                ['.', '', '基因集KEGG富集分析文件',0],
                ['kegg_enrich_geneset_list_gene.xls', '', '基因集KEGG富集分析统计表',0],
                ['pathways', ' ', 'KEGG富集通路图',0,"211539"],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
                ["*bar.pdf", "pdf", "KEGG富集分析柱形图", 0],
                ["*bar_line.pdf", "pdf", "KEGG富集分析柱形图(带折线)", 0],
                ["*buble.pdf", "pdf", "KEGG富集分析气泡图", 0],
                ["*buble2.pdf", "pdf", "KEGG富集分析气泡图(分散型)", 0],
            ])
            result_dir2.add_regexp_rules([
                [r'pathways/.*\.html', '', 'KEGG通路html文件',0,"211540"],
                [r'pathways/.*\.png', '', 'KEGG通路图片png',0,"211541"],
            ])
        super(GenesetEnrichWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
