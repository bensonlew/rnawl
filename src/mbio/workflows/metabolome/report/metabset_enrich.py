
# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
from bson.objectid import ObjectId


class MetabsetEnrichWorkflow(Workflow):
    """
    基因集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        # self.rpc = False
        super(MetabsetEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 by ghd @20191015
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset"},
            {"name": "correct", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "bg", "type": "string", "default": "project"},
             # 背景，project:本项目鉴定到的代谢物合集; species:本物种全部代谢物合集; kegg:KEGG数据库全部代谢物合集
            {"name": "species", "type": "string", "default": "all"},
            {"name": "method", "type": "string", "default": "none"},  ##topo方法
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.enrich_module = self.add_module("metabolome.enrich_topo")
        self.output_dir = self.enrich_module.output_dir
        if self.option("method") == 'rc':
            self.method = 'rbc'
        elif self.option("method") == 'oc':
            self.method = 'rod'
        else:
            self.method = 'none'

    def run(self):
        # super(MetabsetEnrichWorkflow, self).run()
        options = {
            "anno_overview": self.option('anno_overview'),
            "ko_overview": self.option("ko_overview"), # add by ghd @20191015
            "metabset": self.option("metabset"),
            "correct": self.option("correct"),
            "bg": self.option("bg"),
            "species":  self.option("species"),
            "method" :   self.method, #20190620
            "task_id": "_".join(self._sheet.id.split("_")[0:2])
        }
        self.enrich_module.set_options(options)
        self.enrich_module.on('end', self.set_db)
        self.enrich_module.run()
        super(MetabsetEnrichWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        enrich_api = self.api.api('metabolome.enrich')
        output_file = glob.glob("{}/*.xls".format(self.output_dir))

        overview_table = self.option('anno_overview').path
        enrich_table = output_file[0]

        enrich_api.add_enrich_detail(self.option('main_table_id'), overview_table, enrich_table)
        topo_table = self.output_dir + '/topology_png/kegg_topology.xls'
        remote_topo = self._sheet.output + '/topology_png/'  ##self._sheet.output
        enrich_api.add_topology(self.option('main_table_id'),topo_table,remote_topo)
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetkeggenrich",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "代谢集KEGG富集分析结果", 0, "1500042"],
            ["./topology_png", "", "代谢集KEGG富集分析拓扑图", 0],
            ["DE.list.check.kegg_enrichment.xls", "xls", "代谢集KEGG富集分析结果表", 0, "150043"]
        ])
        result_dir.add_regexp_rules([
            [r"topology_png/*png", "png", "结果图", 0],  #20190620
            [r"topology_png/kegg_topology.xls",'xls','结果表']
        ])
        super(MetabsetEnrichWorkflow, self).end()

