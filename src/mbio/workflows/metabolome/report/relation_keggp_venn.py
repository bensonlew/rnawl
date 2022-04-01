# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last_modifiy = modified 2022.1.10

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import os
import types
from mbio.packages.metabolome.get_data import dump_trans_data
from mainapp.models.mongo.metabolome import Metabolome
from biocluster.config import Config


class RelationKeggpVennWorkflow(Workflow):
    """
    转录代谢关联分析--kegg通路注释venn图、统计图
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationKeggpVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_keggp_table_id", "type": "string"},
            {"name": "metab_set_table_id",  "type": "string"},
            {"name": "trans_keggp_main_id", "type": "string"},  # 转录表达量表id
            {"name": "trans_geneset_main_id", "type": "string"},  # 转录基因集表id
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},  # 代谢物详情表用于根据metabid取metabname
            {"name": "species", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.keggp_tool = self.add_tool("metabolome.relation.relation_keggp")
        self.task_id = "_".join(self._sheet.id.split("_")[0:2])
        self.metablome = Metabolome()
        self.metablome._config = Config()
        self.version_value = self.metablome.find_version_from_task(self.task_id)

    def run(self):
        self.keggp_tool.on('end', self.set_db)
        self.run_trans_keggp()
        super(RelationKeggpVennWorkflow, self).run()

    def run_trans_keggp(self):
        self.logger.info("start run_trans_keggp")
        options = {
            "metab_set_table_id":self.option("metab_set_table_id"),
            "trans_keggp_main_id" : self.option("trans_keggp_main_id"),
            "trans_geneset_main_id": self.option("trans_geneset_main_id"),
            "task_id" : self.task_id
        }
        self.keggp_tool.set_options(options)
        self.keggp_tool.run()

    def set_db(self):
        api_name = self.api.api("metabolome.relation_keggp_venn")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="14700501")
        metab_desc = self.option("metab_desc").prop['path']
        trans_keggp_table = self.link_file("trans_keggp_table.xls", "trans_keggp_table.xls")
        metab_keggp_table = self.link_file("metab_keggp_table.xls", "metab_keggp_table.xls")
        gene_id2name = None
        if os.path.exists(os.path.join(self.keggp_tool.work_dir, "gene_id2name.xls")):
            gene_id2name = os.path.join(self.keggp_tool.work_dir, "gene_id2name.xls")
        api_name.add_relation_keggp_venn_detail(main_id, metab_desc = metab_desc, metab_keggp_table=metab_keggp_table,
                                                gene_name2id=gene_id2name, trans_keggp_table=trans_keggp_table,
                                                species_name=self.option("species"), version=self.version_value)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "relationkeggpvenn",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.keggp_tool.output_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)
        return newfile

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        relpath_rules = [
            [".", "", "关联分析相关性结果文件夹", 0, "150069"],
            ["metab_keggp_table.xls", "xls", "代谢物通路注释表", 0, "150070"],
            ["trans_keggp_table.xls", "xls", "基因通路注释表", 0, "150028"]
        ]
        regexps = [
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(RelationKeggpVennWorkflow, self).end()