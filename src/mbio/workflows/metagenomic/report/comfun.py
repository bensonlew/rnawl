# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

"""宏基因组交互分析功能筛选公共类"""


from biocluster.workflow import Workflow


class ComfunWorkflow(Workflow):
    """
    宏基因组交互分析功能筛选公共类
    """

    def __init__(self, wsheet_object):
        self.__sheet = wsheet_object
        super(ComfunWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "database", "type": "string", "required":True},
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_anno", "type": "infile", "format": "sequence.profile_table, annotation.mg_anno_dir"},
            {"name": "lowest_level_profile", "type": "infile", "format": "sequence.profile_table"},
            #{"name": "personal_funset_id", "type": "string"},
            {"name": "personal_funset_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "samples", "type": "string"},
            {"name": "identity", "type": "float", "default": 0.0},
            {"name": "align_length", "type": "int", "default": 0},
            {"name": "level_select", "type": "string", "default": "all"},
            {"name": "abu_num", "type": "string", "default": "all"},
            {"name": "abu_proportion", "type": "string", "default": "all"},
            {"name": "gene_type", "type": "string"},
            {"name": "group", "type": "int", "default": 1},  # 1 所有样本，2 剔除样本
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            # {"name": "xml_file", "type": "infile", "format": "sequence.profile_table"},  # kegg使用
            {"name": "align_database", "type": "string", "default": "all"},  # vfdb使用

        ]
        self.add_option(options)
        self.filter = self.add_module("annotation.mg_fun_select")

    def run_filter(self, *func):
        self.profile = "T"
        self.abu_act = "F"
        self.fun_select = "F"
        if self.option("group") == 1 and self.option("gene_type") == "Origin":
            self.profile = "F"
        if self.option("abu_num") != "all" or self.option("abu_proportion") != "all":
            self.abu_act = "T"
        if self.option("identity") != 0 or self.option("align_length") != 0 or self.option("level_select") != "all":
            self.fun_select = "T"
        if self.option("personal_funset_file").is_set:
            self.fun_select = "T"
        self.logger.info("common fun workflow")
        self.logger.info("profile;abu;fun")
        self.logger.info(self.option("database"))
        self.logger.info(self.option("level_select"))
        self.logger.info(self.profile)
        self.logger.info(self.abu_act)
        self.logger.info(self.fun_select)
        self.logger.info(self.option("group"))
        self.logger.info(self.option("gene_type"))
        options = {
            "gene_anno": self.option("gene_anno"),
            "database": self.option("database"),
            "geneset_table": self.option("geneset_table"),
            "lowest_level_profile": self.option("lowest_level_profile"),
            "gene_profile": self.option("gene_profile"),
            "samples": self.option("samples"),
            "identity": self.option("identity"),
            "align_length": self.option("align_length"),
            "level_select": self.option("level_select"),
            "abu_num": self.option("abu_num"),
            "abu_proportion": self.option("abu_proportion"),
            "abu_filter": self.abu_act ,
            "fun_filter": self.fun_select,
            "gene_filter": self.profile,
            "task_id": "_".join(self._sheet.id.split("_")[0:2]) ## add by qingchen.zhang @20201120
        }
        if self.option("align_database"):
            options["align_database"] = self.option("align_database")
        if self.option("database") == "kegg":
            # options["xml_file"] = self.option("xml_file")
            options["gene_anno"] = self.option("gene_anno").path + "/gene_kegg_anno.xls"
        if self.option("group_table").is_set:
            options["group_table"] = self.option("group_table")
        self.logger.info("test>>>>>>>>>>>>>>>>>>>>>>>>>")
        self.logger.info(self.get_option_object().keys())
        if "personal_funset_file" in self.get_option_object().keys() and self.option("personal_funset_file").is_set:
            options["personal_fun"] = self.option("personal_funset_file")
        self.filter.set_options(options)
        for f in func:
            self.filter.on('end', f)
        self.filter.run()
