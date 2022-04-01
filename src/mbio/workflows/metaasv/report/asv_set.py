# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class AsvSetWorkflow(Workflow):
    """
    metaasv lefse分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AsvSetWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table_name", "type": "string"},
            {"name": "table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "asv_name", "type": "string", "default": ""},
            {"name": "desc", "type": "string", "default": ""},
            {"name": "pvalue", "type": "float", "default": 0.0},
            {"name": "qvalue", "type": "float", "default": 0.0},
            {"name": "species_name", "type": "string", "default": ""},
            {"name": "label", "type": "string"},
            {"name": "lda", "type": "float", "default": 0.0},
            {"name": "top", "type": "int", "default": 0},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False

    def check_options(self):
        """
        参数二次检查
        """
        if self.option("table_name") in ["multiple_group", "two_group", "two_sample"]:
            if (not self.option("pvalue")) and (not self.option("qvalue")) and (not self.option("species_name")):
                raise OptionError("{},{},{}不存在请检查".format(self.option("pvalue"), self.option("qvalue"), self.option("species_name")))
        if self.option("table_name") in ["lefse"]:
            if (not self.option("pvalue")) and (not self.option("lda")):
                raise OptionError("{},{}不存在请检查".format(self.option("pvalue"), self.option("lda")))
        if self.option("table_name") in ["randomforest"]:
            if (not self.option("top")):
                raise OptionError("{}不存在请检查".format(self.option("top")))
        if self.option("table_name") in ["venn"]:
            if (not self.option("species_name")) and (not self.option("label")):
                raise OptionError("{}不存在请检查".format(self.option("species_name"), self.option("label")))

    def end(self):
        self.logger.info("结束啦")
        super(AsvSetWorkflow, self).end()

    def set_db(self):
        """
        直接运行导表函数
        """
        api_set = self.api.api("metaasv.asv_set")
        if self.option("table_name") in ["multiple_group", "two_group", "two_sample"]:
            if self.option("species_name") != "":
                if self.option("pvalue") != 0.0:
                    if self.option("qvalue") != 0.0:
                        api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"),species_name=self.option("species_name"), pvalue=self.option("pvalue"), qvalue=self.option("qvalue"))
                    else:
                        api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"),species_name=self.option("species_name"), pvalue=self.option("pvalue"))
                else:
                    if self.option("qvalue") != 0.0:
                        api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"),species_name=self.option("species_name"), qvalue=self.option("qvalue"))
                    else:
                        api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"),species_name=self.option("species_name"))
            else:
                if self.option("pvalue") != 0.0:
                    if self.option("qvalue") != 0.0:
                        api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"), pvalue=self.option("pvalue"), qvalue=self.option("qvalue"))
                    else:
                        api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"), pvalue=self.option("pvalue"))
                else:
                    if self.option("qvalue") != 0.0:
                        api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"), qvalue=self.option("qvalue"))
        if self.option("table_name") in ["lefse"]:
            if self.option("pvalue") != 0.0:
                if self.option("lda") != 0.0:
                    api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"),pvalue=self.option("pvalue"), lda=self.option("lda"))
                else:
                    api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"),pvalue=self.option("pvalue"))
            else:
                if self.option("lda") != 0.0:
                    api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"), lda=self.option("lda"))

        if self.option("table_name") in ["randomforest"]:
            api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"), top=self.option("top"))
        if self.option("table_name") in ["venn"]:
            api_set.add_main_table(self.option("main_id"), self.option("table_id"),self.option("table_name"), species_name=self.option("species_name"), label=self.option("label"))
        gevent.spawn_later(5, self.end)

    def run(self):
        """
        运行
        """
        self.logger.info("开始啦")
        self.set_db()
        super(AsvSetWorkflow, self).run()
