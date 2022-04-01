# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_file,link_dir
import os,time
import gevent


class VolcanoWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VolcanoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "p_method1", "type": "string", "default": "P"},
            {"name": "p_method2", "type": "string", "default": "less"},
            {"name": "p_method3", "type": "float", "default": 0.05},
            {"name": "vip_method1", "type": "string", "default": "vip_plsda"},
            {"name": "vip_method2", "type": "string", "default": "more"},
            {"name": "vip_method3", "type": "float", "default": 1},
            {"name": "up_method1", "type": "string", "default": "more"},
            {"name": "up_method2", "type": "float", "default": 1},
            {"name": "table", "type": "infile", 'format': "bacgenome.simple_file"},
            {"name": "x_method", "type": "string","default": "log2"},
            {"name": "y_data", "type": "string", "default": "P_value"},
            {"name": "y_method", "type": "string", "default": "-log10"},
            {"name": "vip_method", "type": "string", "default": "true"},
            {"name": "vip", "type": "string", "default": "None"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "diff_gene", "type": "string", "default": "false"},
            {"name": "diff_gene_num", "type": "int", "default": 10},
            {"name": "diff_gene_name", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("table").is_set:
            raise OptionError("必须输入文件！")

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")
        gevent.spawn_later(5, self.set_db)
        super(VolcanoWorkflow, self).run()

    def set_db(self):
        self.logger.info("开始导表")
        api_main = self.api.api("tool_lab.volcano")
        api_main.add_main(file=self.option("table").prop["path"], x_method=self.option("x_method"),y_data=self.option("y_data"),y_method=self.option("y_method"),
                          vip_method=self.option("vip_method"),vip=self.option("vip"), main_id = self.option('main_id'),p_method1=self.option("p_method1"),
                          p_method2=self.option("p_method2"),p_method3=self.option("p_method3"), vip_method1=self.option("vip_method1"),
                          vip_method2=self.option("vip_method2"),vip_method3=self.option("vip_method3"),up_method1=self.option("up_method1"),
                          up_method2=self.option("up_method2"),diff_gene=self.option("diff_gene"),diff_gene_num=self.option("diff_gene_num"),diff_gene_name=self.option("diff_gene_name"))
        self.end()


    def end(self):
        super(VolcanoWorkflow, self).end()