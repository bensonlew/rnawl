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
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
from biocluster.config import Config
from bson import ObjectId

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
            {'name': 'go_version', 'type': 'string', 'default': '2019'},
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
        if self.option("anno_type") == "kegg":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 GeneSet/06 KEGG_Enrich')
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 GeneSet/05 GO_Enrich')
        self.inter_dirs = []
        # self.enrich_tool = self.add_tool('ref_rna_v2.geneset.go_enrich') if self.option('anno_type') == 'go' else self.add_tool('ref_rna_v2.geneset.kegg_rich')
        # self.kegg_class = self.add_tool('ref_rna_v2.geneset.kegg_class')
        # self.output_dir1 = self.enrich_tool.output_dir
        # self.output_dir2 = self.kegg_class.output_dir
        # modify for 重运行 by fwy 20210624
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="ref_rna_v2",
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
        else:
            table = "sg_geneset_kegg_enrich"
        get_run_log = GetRunLog("ref_rna_v2", table=table, main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        if self.option('anno_type') == 'go':
            self.run_go_enrich()
        elif self.option('anno_type') == 'kegg':
            self.run_kegg_rich()


    def run_go_enrich(self):
        self.step.add_steps('go_enrich')
        self.enrich_tool = self.add_tool('ref_rna_v2.geneset.go_enrich')
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
        self.enrich_tool = self.add_tool('ref_rna_v2.geneset.kegg_rich')
        self.enrich_tool.set_options({
            'kegg_table': self.option('kegg_table'),
            'diff_list': self.option('geneset_list'),
            'kegg_table2': self.option('kegg_table_2'),
            "kegg_version": self.option('kegg_version'),
            'correct': self.option('method'),
            'add_info': self.option('add_info')
        })
        self.enrich_tool.on('start', self.set_step, {'start': getattr(self.step, 'kegg_rich')})
        self.enrich_tool.on('end', self.set_step, {'start': getattr(self.step, 'kegg_rich')})
        self.enrich_tool.on('end', self.run_kegg_class)
        self.enrich_tool.run()

    def run_kegg_class(self):
        self.step.add_steps('kegg_class')
        self.kegg_class = self.add_tool('ref_rna_v2.geneset.kegg_class')
        self.kegg_class.set_options({
            'geneset_kegg': self.option('geneset_kegg'),
            'kegg_table': self.option('kegg_table'),
            'kegg_table2': self.option('kegg_table_2'),
            'geneset_id': self.option('geneset_id'),
            'background_links': self.option('add_info'),
            'type': self.option('type'),
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
            result = os.path.join(self.output_dir, os.path.basename(self.enrich_tool.option('result').path))
            shutil.copy(self.enrich_tool.option('result').path, result)
            self.option('result').set_path(result)

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
        os.rename(enrich_result + "relink", enrich_result)
        # os.link(enrich_result + "relink", enrich_result)



    def set_db(self):
        if self.option('anno_type') == 'go':
            api = self.api.api('ref_rna_v2.go_enrich')
            api.add_go_enrich(self.option('result').path, self.option('main_id'), self._sheet.output)
            api.update_db_record('sg_geneset_go_enrich', self.option('main_table_id'),
                                 has_enrich_file="yes")
        elif self.option('anno_type') == 'kegg':
            if self.option("source") == "diff_exp" or not self.option("kegg_version"):
                self.replace_kegg_link()
            self.output_dir1 = self.enrich_tool.output_dir
            self.output_dir2 = self.kegg_class.output_dir

            api_geneset = self.api.api('ref_rna_v2.ref_rna_v2_geneset')
            output_file = glob.glob("{}/*.xls".format(self.output_dir1))
            workflow_output = glob.glob("{}/*.xls".format(self.output_dir1))
            # workflow_output = self.get_workflow_output_dir() + '/' + workflow_output[0].split('/')[-1]
            workflow_output = self.get_workflow_output_dir() + '/' + "kegg_enrich_geneset_list_gene.xls"
            api_geneset.add_kegg_enrich_detail(self.option('main_table_id'), output_file[0], self.option('geneset_list'), self.option('all_list'))
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'), result_dir= workflow_output, has_enrich_file = "yes")
            api_geneset.add_kegg_enrich_pic(self.option('main_table_id'), output_file[0], self.output_dir2 + '/pathways', source=self.option("source"))
            graph_dir = os.path.join(self.get_workflow_output_dir(), 'pathways')
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'), graph_dir=graph_dir)

        # if self.option('anno_type') == 'go':
        #     self.output_dir1 = self.enrich_tool.output_dir
        # elif self.option('anno_type') == 'kegg':
        #     self.output_dir2 = self.kegg_class.output_dir
        # api_geneset = self.api.api('ref_rna_v2.ref_rna_v2_geneset')
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
            db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
            geneset_info = db["sg_geneset"].find_one({"_id": ObjectId(self.option("geneset_id"))})
            name = geneset_info["name"]
        except:
            pass
        return name

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"

        if self.option('anno_type') == 'go':
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
            kegg_enrich_table = os.path.join(self.output_dir1, "geneset_list_gene.list.DE.list.check.kegg_enrichment.xls")
            with open(kegg_enrich_table, "r") as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    chart.chart_geneset_enrich_kegg(kegg_enrich_table, geneset_name=self.get_geneset_name())
                    chart.to_pdf()
                    # move pdf to result dir
                    pdf_file = glob.glob(self.work_dir + "/*.pdf")
                    for p in pdf_file:
                        linkfile(p, self.output_dir + "/" + os.path.basename(p))

    def end(self):
        self.chart()
        if self.option('anno_type') == 'go':
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["03 GeneSet", "", "基因集分析结果目录", 0],
                ["03 GeneSet/05 GO_Enrich", "", "GO功能富集", 0]
            ]
            result_dir.add_relpath_rules([
                ['.', '', 'GO富集分析文件',0,"211535"],
                ['go_enrich_geneset_list_gene.xls', 'xls', 'GO富集分析统计表',0,"211536"],
                ["*bar.pdf", "pdf", "GO富集分析柱形图",0],
                ["*bar_line.pdf", "pdf", "GO富集分析柱形图(带折线)", 0],
                ["*buble.pdf", "pdf", "GO富集分析气泡图(单基因集)", 0],
                ["*buble2.pdf", "pdf", "GO富集分析气泡图(分散型)", 0],
                ['./run_parameter.txt', 'txt', '运行参数日志', 0],
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
            self.inter_dirs = [
                ["03 GeneSet", "", "基因集分析结果目录", 0],
                ["03 GeneSet/06 KEGG_Enrich", "", "KEGG功能富集", 0]
            ]
            result_dir2.add_relpath_rules([
                ['.', '', 'KEGG富集分析文件',0,"211537"],
                ['kegg_enrich_geneset_list_gene.xls', '', 'KEGG富集分析统计表',0,"211538"],
                ['pathways', ' ', 'KEGG富集通路图',0,"211539"],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
                ["*bar.pdf", "pdf", "KEGG富集分析柱形图", 0],
                ["*bar_line.pdf", "pdf", "KEGG富集分析柱形图(带折线)", 0],
                ["*buble.pdf", "pdf", "KEGG富集分析气泡图(单基因集)", 0],
                ["*buble2.pdf", "pdf", "KEGG富集分析气泡图(分散型)", 0],
            ])
            result_dir2.add_regexp_rules([
                [r'pathways/map.*\.html', '', 'KEGG通路html文件',0,"211540"],
                [r'pathways/map.*\.png', '', 'KEGG通路图片png',0,"211541"],
            ])
        super(GenesetEnrichWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
