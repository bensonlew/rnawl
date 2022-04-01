# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

"""宏基因组交互分析功能筛选公共类"""


from biocluster.workflow import Workflow


class CommTableWorkflow(Workflow):
    """
    宏基因组交互分析功能筛选公共类
    """

    def __init__(self, wsheet_object):
        self.__sheet = wsheet_object
        super(CommTableWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_id", "type": "string"},
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},  # 各数据库的注释表格
            {"name": "func_anno_table", "type": "infile", "format": "meta.profile"},  # 各数据库的注释表格
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name": "method", "type": "string", "default": "rpkm"},
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},  # 同geneset_table，兼容不同程序
            {"name": "level_id", "type": "string"},
            {"name": "level", "type": "string"}, # 同level_id,兼容不同程序
            {"name": "func_level", "type": "string", "default": ""}, # 同level_id,兼容不同程序
            {"name": "level_type", "type": "string"}, # 同level_id,兼容不同程序
            {"name": "second_level", "type": "string"},
            {"name": "level_type_name", "type": "string", "default": ""}, #同second_level,兼容不同程序
            {"name": "lowest_level", "type": "string", "default": ""},  # 注释表数据库对应的最低分类，eg：KEGG的ko
            {"name": "anno_type", "type": "string"},
            {"name": "level_color", "type": "string", "default": ""},   #add by qingchen.zhang@20181202 用于bubble和heatmap组成分析的颜色水平的筛选
            {"name": "graphic_type", "type": "string", "default": ""},  # bar,heatmap,circos，bubble
            {"name": "clean_stat", "type": "infile", "format": "meta.profile"},
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def run_abundance(self, *func):
        self.other_step = func
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            self.run_ppm_abun_one()
        else:
            self.run_com_abun()

    def run_ppm_abun_one(self):
        self.logger.info("start cal ppm abundance table>>>>>>>>>>>")
        self.abundance = self.add_tool("meta.cal_ppm")
        options = {
            'clean_stat': self.option('clean_stat').prop["path"],
            'geneset_table': self.option('geneset_table').prop["path"],
        }
        self.abundance.set_options(options)
        self.abundance.on('end', self.run_ppm_abun_two)
        #self.abundance.run()

    def run_ppm_abun_two(self):
        final_geneset_table = self.abundance.option("out_table").prop["path"]
        if self.option("anno_type") == "go":
            self.logger.info("now start step two go database")
            self.abundance_two = self.add_module("annotation.anno_go_stat")
            options = {
                #'gene_anno': self.option('anno_table').prop["path"],
                'go1234level_out': self.option('go1234level_out'),
                'reads_profile_table': final_geneset_table,
                #'level': self.option('level_id'),
            }
            if self.option("anno_table").is_set:
                options["gene_anno"] = self.option('anno_table').prop["path"]
            elif self.option("func_anno_table").is_set:
                options["gene_anno"] = self.option('func_anno_table').prop["path"]
            if self.option('level_id'):
                options["level"] = self.option('level_id')
            elif self.option('level'):
                options["level"] = self.option('level')
            elif self.option('func_level'):
                options["level"] = self.option('func_level')
            elif self.option('level_type'):
                options["level"] = self.option('level_type')
            self.abundance_two.set_options(options)
        else:
            self.abundance_two = self.add_tool("meta.create_abund_table")
            options = {
                'anno_table': self.option('anno_table'),
                'geneset_table': final_geneset_table,
                'gene_list': self.option('gene_list'),
                'level_type': self.option('level_id'),
                #'level_type_name': self.option('second_level'),
                'lowest_level': self.option('lowest_level'),
            }
            if self.option("second_level"):
                options["level_type_name"] = self.option("second_level")
            elif self.option("level_type_name"):
                options["level_type_name"] = self.option("level_type_name")
            if self.option("level_color") != "":
                options["level_color"] = self.option("level_color")
            if self.option("graphic_type") != "":
                options["graphic_type"] = self.option("graphic_type")
            self.abundance_two.set_options(options)
        for f in self.other_step:
            self.abundance_two.on('end',self.set_table)
            self.abundance_two.on('end', f)
            self.abundance_two.run()

    def set_table(self):
        self.abundance.option("out_table", self.abundance_two.option("out_table").prop["path"])

    def run_com_abun(self):
        self.logger.info("start creat abundance table>>>>>>>>>>>")
        self.final_geneset_table = self.option("geneset_table")
        if self.option("anno_type") == "go":
            self.run_go_abu()
        else:
            self.abundance = self.add_tool("meta.create_abund_table")
            options = {
                'anno_table': self.option('anno_table'),
                'geneset_table': self.final_geneset_table,
                'gene_list': self.option('gene_list'),
                'level_type': self.option('level_id'),
                #'level_type_name': self.option('second_level'),
                'lowest_level': self.option('lowest_level'),
                'anno_type': self.option("anno_type")
            }
            if self.option("second_level"):
                options["level_type_name"] = self.option("second_level")
            elif self.option("level_type_name"):
                options["level_type_name"] = self.option("level_type_name")
            if self.option("level_color") != "":
                options["level_color"] = self.option("level_color")
            if self.option("graphic_type") != "":
                options["graphic_type"] = self.option("graphic_type")
            if self.option("anno_type") == "card" and self.option("level_id") == "Class":
                if "level_type_name" in options:
                    options["level_type_name"] = options["level_type_name"].split(' ')[0]
                options["anno_type"] = "card"
            self.abundance.set_options(options)
        for f in self.other_step:
            self.abundance.on('end', f)

    def run_cal_ppm(self, **kargs):
        self.logger.info("start calculate ppm>>>>>>>>>>>>>>")
        if kargs.has_key("tool_name"):
            tool_name = kargs["tool_name"]
        else:
            self.cal_ppm = self.add_tool("meta.cal_ppm")
            tool_name = self.cal_ppm
        options = {
            'clean_stat': self.option('clean_stat').prop["path"],
            #'geneset_table': self.option('geneset_table').prop["path"],
        }
        if self.option("geneset_table").is_set:
            options["geneset_table"] = self.option('geneset_table').prop["path"]
        elif self.option("gene_profile").is_set:
            options["geneset_table"] = self.option('gene_profile').prop["path"]
        tool_name.set_options(options)
        if kargs.has_key("rely") and kargs["rely"]:
            tool_name.on('end', kargs["rely"])
        else:
            tool_name.on('end', self.run_com_abun)
        if kargs.has_key("is_run") and kargs["is_run"]:
            tool_name.run()

    def run_go_abu(self, rely=None, is_run=False, level=None):
        self.abundance = self.add_module("annotation.anno_go_stat")
        options = {
            #'gene_anno': self.option('anno_table').prop["path"],
            'go1234level_out': self.option('go1234level_out'),
            #'reads_profile_table': self.option('geneset_table').prop["path"],
            #'level': self.option('level_id'),
        }
        if self.option("anno_table").is_set:
            options["gene_anno"] = self.option('anno_table').prop["path"]
        elif self.option("func_anno_table").is_set:
            options["gene_anno"] = self.option('func_anno_table').prop["path"]
        if self.option('level_id'):
            options["level"] = self.option('level_id')
        elif self.option('level'):
            options["level"] = self.option('level')
        elif self.option('func_level'):
            options["level"] = self.option('func_level')
        if self.option("geneset_table").is_set:
            options["reads_profile_table"] = self.option('geneset_table').prop["path"]
        elif self.option("gene_profile").is_set:
            options["reads_profile_table"] = self.option('gene_profile').prop["path"]
        if level:
            options["level"] = level
        self.abundance.set_options(options)
        if rely:
            self.abundance.on('end', rely)
        if is_run:
            self.abundance.run()



