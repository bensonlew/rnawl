# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import types
from mainapp.models.mongo.metabolome import Metabolome
from bson.objectid import ObjectId
from mbio.packages.metabolome.get_data import dump_trans_data
from biocluster.config import Config

class RelationKeggpEnrichWorkflow(Workflow):
    """
    基因集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationKeggpEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 by ghd @20191015
            {"name": "metab_set_table", "type": "infile", "format": "metabolome.metabset"},
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "trans_geneset_main_id", "type": "string"},
            {"name": "trans_kegg_main_id", "type": "string"},
            {"name": "correct", "type": "string"},  # 多重检验方法
            {"name": "species", "type": "string", "default": "all"},  # 代谢富集背景为物种
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.metab_enrich_tool = self.add_tool("metabolome.metabset.enrich")
        self.trans_enrich_tool = self.add_tool("metabolome.relation.trans_kegg_enrich")
        self.task_id = "_".join(self._sheet.id.split("_")[0:2])
        self.metablome = Metabolome()
        self.metablome._config = Config()
        self.version_value = self.metablome.find_version_from_task(self.task_id)

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run Relation_keggp_enrich workflow")
        self.select_tools = [self.metab_enrich_tool, self.trans_enrich_tool]
        self.on_rely(self.select_tools, self.set_db)
        self.run_metab_enrich()
        self.run_trans_enrich()
        super(RelationKeggpEnrichWorkflow, self).run()

    def run_metab_enrich(self):
        opts = {
            'anno_overview': self.option("anno_overview"),
            'ko_overview': self.option("ko_overview"),
            'metabset': self.option("metab_set_table"),
            'correct': self.option('correct'),
            'bg': "species",
            'species' : self.option("species"),
            "version": self.version_value
        }
        self.metab_enrich_tool.set_options(opts)
        self.metab_enrich_tool.run()

    def run_trans_enrich(self):
        metab_client = Config().get_mongo_client(mtype="metabolome")
        relation_db = metab_client[Config().get_mongo_dbname("metabolome")]
        relation_info = relation_db['sg_relation_analysis'].find_one({"task_id": self.task_id, "delete_by": ""})
        relate_task_id = relation_info["relate_task_id"]
        relate_project_type = relation_info["relate_project_type"]
        self.logger.info("relate_task_id为{}".format(relate_task_id))
        self.logger.info("relate_project_type为{}".format(relate_project_type))
        trans_client = Config().get_mongo_client(mtype=relate_project_type)
        trans_db = trans_client[Config().get_mongo_dbname(relate_project_type)]
        if relate_project_type == "whole_transcriptome":
            sg_task = trans_db["task"].find_one({"task_id": relate_task_id})
            self.logger.info("sg_task为{}".format(sg_task))
            try:  # 转录一些旧项目存在没有kegg版本库的问题，没有的话默认为202007版本
                version = sg_task["long_task"]["database_version"]["kegg"]
            except:
                self.logger.info("没有找到转录kegg数据库版本")
                version = "202007"
        else:
            sg_task = trans_db["sg_task"].find_one({"task_id": relate_task_id})
            try:
                version = sg_task["database_version"]["kegg"]
            except:
                self.logger.info("没有找到转录kegg数据库版本")
                version = "202007"
        self.logger.info("version为{}".format(version))
        opts = {
            "trans_geneset_main_id": self.option("trans_geneset_main_id"),
            "trans_kegg_main_id": self.option("trans_kegg_main_id"),
            'correct': self.option('correct'),
            "task_id": self.task_id,
            "version": version
        }
        self.trans_enrich_tool.set_options(opts)
        self.trans_enrich_tool.run()

    def set_db(self):
        api_name = self.api.api("metabolome.relation_kegg_enrich")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700501")
        metab_enrich_result = self.metab_enrich_tool.output_dir + "/DE.list.check.kegg_enrichment.xls"
        os.link(metab_enrich_result, self.output_dir + "/metabolome.DE.list.check.kegg_enrichment.xls")
        trans_enrich_result = self.trans_enrich_tool.work_dir + "/DE.list.check.kegg_enrichment.xls"
        os.link(trans_enrich_result, self.output_dir + "/transcript.DE.list.check.kegg_enrichment.xls")
        metab_desc = self.option("metab_desc").prop['path']
        gene_name2id = None
        if os.path.exists(os.path.join(self.trans_enrich_tool.work_dir, "gene_id2name.xls")):
            gene_name2id = os.path.join(self.trans_enrich_tool.work_dir, "gene_id2name.xls")
        api_name.add_relation_keggp_enrich_detail(main_id, metab_enrich_result, trans_enrich_result, metab_desc=metab_desc,
                                                  gene_name2id=gene_name2id, species_name=self.option("species"),
                                                  version=self.version_value)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "relationcorr",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()
            
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "关联分析kegg富集分析结果文件夹", 0, "150069"],
            ["metabolome.DE.list.check.kegg_enrichment.xls", "xls", "代谢富集结果文件", 0, "150070"],
            ["transcript.DE.list.check.kegg_enrichment.xls", "xls", "转录富集结果文件", 0, "150028"],
        ]
        regexps = [
            [r"metab\.subcluster_.*\.xls", "xls", "代谢物相关性kmeans各子类结果表", 0, "150068"],
            [r"gene\.subcluster_.*\.xls", "xls", "基因相关性kmeans各子类结果表", 0, "150074"],
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(RelationKeggpEnrichWorkflow, self).end()