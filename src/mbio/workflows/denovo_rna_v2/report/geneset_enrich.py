
# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import os
import json
from mbio.packages.denovo_rna_v2.chart_geneset import ChartGeneset
import glob
from bson import ObjectId


class GenesetEnrichWorkflow(Workflow):
    """
    基因集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "kegg_table", "type": "string"},
            {"name": "go_list", "type": "string"},
            {'name': 'go_version', 'type': 'string', 'default': '2019'},
            {"name": "geneset_list", "type": "string"},
            # {"name": "all_list", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "type", "type": "string"},
            {"name": "add_info", "type": "string", "default": None}  # 输入两列的列表文件，有head，第一列为pathway，第二列为底图链接
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "go":
            self.enrich_tool = self.add_tool("denovo_rna_v2.go_enrich")
        else:
            self.enrich_tool = self.add_tool("denovo_rna_v2.kegg_rich")
        # self.output_dir = self.enrich_tool.output_dir
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/05 GO_Enrich')
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
        super(GenesetEnrichWorkflow, self).send_log(data)

    def run(self):
        if self.option("anno_type") == "kegg":
            options = {
                "kegg_table": self.option("kegg_table"),
                # "all_list": background_path,
                "diff_list": self.option("geneset_list"),
                "correct": self.option("method"),
                "add_info": self.option("add_info")
            }
        else:
            options = {
                "diff_list": self.option("geneset_list"),
                # "all_list": background_path,
                "go_list": self.option("go_list"),
                'go_version': self.option('go_version'),
                # "pval": self.option("pval"),
                "method": self.option("method"),
            }
        self.logger.info(options)
        self.enrich_tool.set_options(options)
        self.enrich_tool.on('end', self.set_db)
        self.get_run_log()
        self.enrich_tool.run()
        super(GenesetEnrichWorkflow, self).run()

    def get_run_log(self):
        if self.option('anno_type') == 'go':
            table = "sg_geneset_go_enrich"
        else:
            table = "sg_geneset_kegg_enrich"
        get_run_log = GetRunLog("denovo_rna_v2", table=table, main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('denovo_rna_v2.geneset_new')
        output_file = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))
        # png_file = glob.glob("{}/*.png".format(self.output_dir))
        # go_png = self.output_dir + "/go_lineage.png"
        # go_pdf = self.output_dir + "/go_lineage.pdf"
        #go_adjust_png = self.output_dir + "/adjust_lineage.png"
        #go_adjust_pdf = self.output_dir + "/adjust_lineage.pdf"
        if self.option("anno_type") == "kegg":
            api_geneset.add_kegg_enrich_detail(self.option("main_table_id"), output_file[0])

        else:
            api_geneset = self.api.api('denovo_rna_v2.geneset')
            result = os.path.join(self.enrich_tool.output_dir,"go_enrich_geneset_list_gene.xls")
            api_geneset.add_go_enrich_detail(self.option("main_table_id"),result)
            workflow_output = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))
            workflow_output = self.get_workflow_output_dir() + '/' + workflow_output[0].split('/')[-1]
            api_geneset.update_db_record('sg_geneset_go_enrich', self.option('main_table_id'),
                                         result_dir=workflow_output)
            # api_geneset.update_directed_graph(self.option("main_table_id"), go_adjust_png, go_adjust_pdf)
        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

    def get_geneset_name(self):
        name = "geneset"
        try:
            db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("denovo_rna_v2")]
            geneset_info = db["sg_geneset"].find_one({"_id": ObjectId(self.option("geneset_id"))})
            name = geneset_info["name"]
        except:
            pass
        return name

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"

        if self.option('anno_type') == 'go':
            go_enrich_table = os.path.join(self.enrich_tool.output_dir,"go_enrich_geneset_list_gene.xls")
            with open(go_enrich_table, "r") as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    chart.denovo_chart_geneset_enrich_go(go_enrich_table, geneset_name=self.get_geneset_name())
                    chart.to_pdf()
                    # move pdf to result dir
                    pdf_file = glob.glob(self.work_dir + "/*.pdf")
                    for p in pdf_file:
                        os.link(p, self.enrich_tool.output_dir + "/" + os.path.basename(p))


    def end(self):
        rm_file = glob.glob(self.enrich_tool.output_dir+ "/*.pdf") + glob.glob(self.enrich_tool.output_dir + "/*.png")
        for file in rm_file:
            os.remove(file)
        if os.path.exists(os.path.join(self.enrich_tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.enrich_tool.output_dir, os.path.basename(self.run_log)))
        self.chart()
        os.link(self.run_log, os.path.join(self.enrich_tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.enrich_tool.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/05 GO_Enrich", "", "GO功能富集", 0]
        ]
        if self.option("anno_type") == "go":
            result_dir.add_relpath_rules([
                [".", "", "GO富集分析文件",0,],
                ['go_enrich_geneset_list_gene.xls', 'xls', 'GO富集分析统计表',0],
                ["*bar.pdf", "pdf", "GO富集分析柱形图", 0],
                ["*bar_line.pdf", "pdf", "GO富集分析柱形图(带折线)", 0],
                ["*buble1.pdf", "pdf", "GO富集分析气泡图(单基因集)", 0],
                ["*buble2.pdf", "pdf", "GO富集分析气泡图(分散型)", 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
        elif self.option("anno_type") == "kegg":
            result_dir.add_relpath_rules([
                [".", "", "KEGG富集分析文件", 0],
                ['kegg_enrich_geneset_list_gene.xls', 'xls', 'KEGG富集分析统计表',0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
        super(GenesetEnrichWorkflow, self).end()


