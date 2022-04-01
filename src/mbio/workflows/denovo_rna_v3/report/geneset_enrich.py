# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import json

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import re
import shutil
from bson.objectid import ObjectId
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
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
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "geneset_list", "type": "string"},
            {"name": "type", "type": "string"},
            #{"name": "all_list", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "source", "type": "string"},
            {"name": "geneset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "string"},
            {'name': 'kegg_table_2', 'type': 'string'},
            {"name": "go_list", "type": "string"},
            {'name': 'result', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {"name": "add_info", "type": "string", "default": None},  # 输入两列的列表文件，有head，第一列为pathway，第二列为底图链接
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

        if self.option("anno_type") == "go":
            self.enrich_tool = self.add_tool("denovo_rna_v2.geneset.go_enrich")
        else:
            self.enrich_tool = self.add_tool("denovo_rna_v3.geneset.kegg_enrich")
            self.kegg_class = self.add_tool("denovo_rna_v3.geneset.kegg_class")
            self.output_dir2 = self.kegg_class.output_dir
        self.output_dir1 = self.enrich_tool.output_dir
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/06 KEGG_Enrich')
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
                "kegg_version": self.option('kegg_version'),
                # "all_list": background_path,
                "diff_list": self.option("geneset_list"),
                'kegg_table2': self.option('kegg_table_2'),
                "correct": self.option("method"),
                "add_info": self.option("add_info")
            }
        else:
            options = {
                "diff_list": self.option("geneset_list"),
                # "all_list": background_path,
                "go_list": self.option("go_list"),
                # "pval": self.option("pval"),
                "method": self.option("method"),
            }
        self.logger.info(options)
        self.enrich_tool.set_options(options)
        if self.option("anno_type") == "kegg":
            self.enrich_tool.on('end', self.run_kegg_class)
            self.kegg_class.on('end', self.set_db)
        else:
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
        api_geneset = self.api.api('denovo_rna_v3.geneset_new')
        output_file = glob.glob("{}/*.xls".format(self.output_dir1))
        if self.option('anno_type') == 'kegg':
            if self.option("source") == "diff_exp":
                self.replace_kegg_link()
            self.output_dir1 = self.enrich_tool.output_dir
            self.output_dir2 = self.kegg_class.output_dir
        workflow_output = glob.glob("{}/*.xls".format(self.output_dir1))
        workflow_output = self.get_workflow_output_dir() + '/' + workflow_output[0].split('/')[-1]
        # png_file = glob.glob("{}/*.png".format(self.output_dir))
        # go_png = self.output_dir + "/go_lineage.png"
        # go_pdf = self.output_dir + "/go_lineage.pdf"
        # go_adjust_png = self.output_dir1 + "/adjust_lineage.png"
        # go_adjust_pdf = self.output_dir1 + "/adjust_lineage.pdf"
        if self.option("anno_type") == "kegg":
            api_geneset.add_kegg_enrich_detail(self.option("main_table_id"), output_file[0], self.option("geneset_list"))
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_table_id'),result_dir=workflow_output)
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option("main_table_id"), result_dir= workflow_output)
            api_geneset.add_kegg_enrich_pic(self.option("main_table_id"), output_file[0], self.output_dir2 + "/pathways",source=self.option("source"))
            graph_dir = os.path.join(self.get_workflow_output_dir(), 'pathways')
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option("main_table_id"),graph_dir=graph_dir)
        else:
            api_geneset = self.api.api('denovo_rna_v2.geneset')
            api_geneset.add_go_enrich_detail(self.option("main_table_id"), output_file[0])
            # api_geneset.update_directed_graph(self.option("main_table_id"), go_adjust_png, go_adjust_pdf)
            api_geneset.update_db_record('sg_geneset_go_enrich', self.option("main_table_id"),
                                         result_dir=workflow_output)
        self.end()

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
            go_enrich_table = self.enrich_tool.option('result').path
            with open(go_enrich_table, "r") as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    chart.denovo_chart_geneset_enrich_go(go_enrich_table, geneset_name=self.get_geneset_name())
                    chart.to_pdf()
                    # move pdf to result dir
                    pdf_file = glob.glob(self.work_dir + "/*.pdf")
                    for p in pdf_file:
                        os.link(p, self.output_dir + "/" + os.path.basename(p))
        elif self.option('anno_type') == 'kegg':
            kegg_enrich_table = os.path.join(self.output_dir1, "geneset_list_gene.list.DE.list.check.kegg_enrichment.xls")
            with open(kegg_enrich_table, "r") as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    chart.denovo_chart_geneset_enrich_kegg(kegg_enrich_table, geneset_name=self.get_geneset_name())
                    chart.to_pdf()
                    # move pdf to result dir
                    pdf_file = glob.glob(self.work_dir + "/*.pdf")
                    for p in pdf_file:
                        os.link(p, self.output_dir2 + "/" + os.path.basename(p))

    def end(self):
        if self.option("anno_type") == "go":
            if os.path.exists(os.path.join(self.output_dir1, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir1, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir1, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir1)
            result_dir.add_relpath_rules([
                [".", "", "基因集GO富集分析结果文件", 0, "201478"],
                ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
        elif self.option("anno_type") == "kegg":
            os.link(os.path.join(self.output_dir1, "geneset_list_gene.list.DE.list.check.kegg_enrichment.xls"),
                    os.path.join(self.output_dir2, "kegg_enrich_geneset_list_gene.xls"))
            shutil.rmtree(self.output_dir2 + "/ko")
            rm_file = glob.glob(self.output_dir2 + "/pathways/*.html.mark") + glob.glob(
                self.output_dir2 + "/kegg_stat.xls")
            for file in rm_file:
                os.remove(file)
            if os.path.exists(os.path.join(self.output_dir2, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir2, os.path.basename(self.run_log)))
            self.chart()
            os.link(self.run_log, os.path.join(self.output_dir2, os.path.basename(self.run_log)))
            result_dir2 = self.add_upload_dir(self.output_dir2)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录",0],
                ["04 GeneSet/06 KEGG_Enrich", "", "KEGG富集分析", 0]
            ]
            result_dir2.add_relpath_rules([
                ['.', '', 'KEGG富集分析文件',0],
                ['kegg_enrich_geneset_list_gene.xls', 'xls', 'KEGG富集分析统计表',0,"201480"],
                ['pathways', '', 'KEGG富集通路图',0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
                ["*bar.pdf", "pdf", "KEGG富集分析柱形图", 0],
                ["*bar_line.pdf", "pdf", "KEGG富集分析柱形图(带折线)", 0],
                ["*buble1.pdf", "pdf", "KEGG富集分析气泡图(单基因集)", 0],
                ["*buble2.pdf", "pdf", "KEGG富集分析气泡图(分散型)", 0],
            ])
            result_dir2.add_regexp_rules([
                [r'pathways/map.*\.html', '', 'KEGG通路html文件',0,"201482"],
                [r'pathways/map.*\.png', '', 'KEGG通路png图片',0,"201483"],
            ])
        super(GenesetEnrichWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

    def run_kegg_class(self):
        opts = {
            "geneset_kegg": self.option('geneset_kegg'),
            "kegg_table": self.option("kegg_table"),
            'kegg_table2': self.option('kegg_table_2'),
            "geneset_id": self.option("geneset_id"),
            "background_links": self.option("add_info"),
            "type": self.option("type"),
            "task_id": self.option("task_id"),
            "kegg_version": self.option('kegg_version'),
            "source": self.option('source')
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()

    def replace_kegg_link(self):
        '''
        替换kegg链接
        '''
        enrich_result = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        annot_result = self.kegg_class.output_dir + '/kegg_stat.xls'
        with open(annot_result, 'rb') as f:
            map2link = {line.split("\t")[0]:line.strip().split("\t")[-1] for line in f.readlines()[1:] if line.split("\t")[-1].startswith("http")}
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
