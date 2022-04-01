# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:2018.2.6
# last_modify by: shaohua.yuan

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class MgFunSelectModule(Module):
    def __init__(self, work_id):
        super(MgFunSelectModule, self).__init__(work_id)
        options = [
            {"name": "geneset_table", "type": "infile", "format": "sequence.profile_table"},  # 基因list file，以GeneID为开头
            {"name": "gene_anno", "type": "infile", "format": "sequence.profile_table"},  # 基因注释file
            {"name": "lowest_level_profile", "type": "infile", "format": "sequence.profile_table"},  #最低层级丰度表
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "personal_fun", "type": "infile", "format": "sequence.profile_table"}, # 个性化功能集使用
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},  # 基因丰度文件
            {"name": "samples", "type": "string"},  # 选择样本
            {"name": "identity", "type": "float", "default": 0.0},
            {"name": "align_length", "type": "int", "default": 0},
            {"name": "level_select", "type": "string", "default": "all"},  # 功能选择水平
            {"name": "database", "type": "string"},
            {"name": "abu_num", "type": "string", "default": "all"},  # 在n个样本中最低层级丰度大于等于
            {"name": "abu_proportion", "type": "string", "default": "all"},  # 最低层级丰度占比
            {"name": "abu_filter", "type": "string", "default": "F"},
            {"name": "fun_filter", "type": "string", "default": "F"},
            {"name": "gene_filter", "type": "string", "default": "T"},
            {"name": "final_gene_anno", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "xml_file", "type": "infile", "format": "sequence.profile_table"},  # kegg使用
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # vfdb使用
            {"name": "align_database", "type": "string"},  # vfdb使用
            {"name": "filter_zero", "type": "string", "default": "True"},#对total为0的丰度进行过滤add by qingchen.zhang@20190315
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老任务
        ]
        self.add_option(options)
        self.logger.info("test>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        self.gene_profile_select = self.add_tool("sequence.select_table")  # 基因profile筛选
        self.fun_filter_tool = self.add_tool("meta.annotation.anno_fun_select")  # 功能筛选tool
        self.abu_filter_tool = self.add_tool("meta.annotation.anno_abu_select")  # 丰度筛选tool
        self.creat_lowest_tool = self.add_tool('meta.association_model.creat_level_table')  # 新最低层级file生成
        #self.recal_abu_tool = self.add_tool("meta.annotation.recal_anno_abu")  # 重新计算新的各层级丰度
        self.recal_abu_tool = self.add_module("annotation.recal_anno_abu")  # 重新计算新的各层级丰度
        self.lowest_select_tool = self.add_tool("meta.annotation.anno_fun_select")  # 最低层级作为功能筛选
        self.lowest_names = {"nr": "Species", "cog": "NOG", "kegg": "Gene", "cazy": "Family", "ardb": "ARG",
                             "card": "ARO", "vfdb": "VFs", "probio": "Probiotic_name", "go": "GO Term (Lev4)",
                             "phi": "protein", "mvirdb": "Virulence Factor ID", "qs": "QS_id", "pfam": "Pfam ID",
                            "tcdb": "TCDB ID", "p450": "Sid"}
        self.final_gene_anno = ""
        self.final_gene_profile = ""
        self.new_lowest_file = ""
        self.vfdb_only_database = "F"

    def check_options(self):
        if not self.option("database"):
            raise OptionError("必须设置数据库名称database!", code="21201201")
        if not self.option("gene_anno").is_set:
            raise OptionError("必须设置基因注释文件!", code="21201202")
        if not self.option("gene_profile").is_set:
            raise OptionError("必须设置基因丰度文件!", code="21201203")
        if self.option("abu_filter") == "T":
            if not self.option("lowest_level_profile").is_set:
                raise OptionError("进行丰度筛选必须设置最低层级丰度文件!", code="21201204")
        if self.option("gene_filter") == "T":
            if not self.option("geneset_table").is_set:
                raise OptionError("基因profile筛选时必须输入基因集list文件!", code="21201205")
        if self.option("fun_filter") == "T":
            if self.option("identity") == 0 and self.option("align_length") == 0:
                if self.option("level_select") == "all" and not self.option("personal_fun").is_set:
                    raise OptionError("进行功能筛选必须设置identity或align_length或功能水平level_select", code="21201206")
        # if self.option("database") == "kegg":
            # if not self.option("xml_file").is_set:
                # raise OptionError("database为kegg时必须输入xml_file文件！", code="21201207")
        if self.option("database") == "vfdb":
            if not self.option("align_database"):
                raise OptionError("database为vfdb时必须输入align_database参数！", code="21201208")
        return True

    def set_output(self, event):
        all_files = os.listdir(self.recal_abu_tool.output_dir)
        if self.option("database") == "kegg":
            all_files = ["kegg_enzyme_profile.xls", "kegg_gene_profile.xls", "kegg_KO_profile.xls",
                         "kegg_level1_profile.xls", "kegg_level2_profile.xls", "kegg_level3_profile.xls",
                         "kegg_module_profile.xls", "kegg_pathway_eachmap.xls", "kegg_pathway_profile.xls"]
        for i in all_files:
            old = os.path.join(self.recal_abu_tool.output_dir, i)
            link = os.path.join(self.output_dir, i)
            if os.path.exists(link):
                os.remove(link)
            os.link(old, link)
        self.option("final_gene_anno", self.final_gene_anno)
        self.end()

    def run(self):
        super(MgFunSelectModule, self).run()
        self.logger.info("start mg_select_fun module")
        if self.option("gene_filter") == "F":
            if self.option("fun_filter") == "T":
                self.gene_anno_sel = "T"
            elif self.option("abu_filter") == "T":
                self.gene_anno_sel = "F"
            else:
                if self.option("database") in ["nr","cog","kegg","ardb","card","cazy","vfdb","go", "qs", "probio", "pfam", "p450", "tcdb", "mvirdb", "phi"]:
                    name = self.option("database").upper()
                    self.logger.info("分组方案仅用于筛选样本，当前参数没做任何筛选退出运行,请查看{}_Origin结果！".format(name))
                    self.set_error("分组方案仅用于筛选样本，当前参数没做任何筛选退出运行,请查看%s_Origin结果！", variables=(name), code="21201201")
                else:
                    self.gene_anno_sel = "F"
            self.final_gene_profile = self.option("gene_profile")
        else:
            self.gene_anno_sel = "T"
        # 定义依赖关系
        if self.gene_anno_sel == "T" and self.option("abu_filter") == "T":
            self.fun_filter_tool.on('end', self.creat_lowest_file)
            self.creat_lowest_tool.on('end', self.abundance_select)
            self.abu_filter_tool.on('end', self.anno_select_lowest)
            self.lowest_select_tool.on('end', self.run_profile)
        elif self.gene_anno_sel == "F" and self.option("abu_filter") == "T":
            self.abu_filter_tool.on('end', self.anno_select_lowest)
            self.lowest_select_tool.on('end', self.run_profile)
        elif self.gene_anno_sel == "T" and self.option("abu_filter") == "F":
            self.fun_filter_tool.on('end', self.run_profile)
        self.recal_abu_tool.on('end', self.set_output)
        ### 确定开始函数
        if self.option("gene_filter") == "F":
            if self.option("fun_filter") == "T":
                self.run_anno_select()
            elif self.option("abu_filter") == "T":
                self.abundance_select()
            elif self.option("database") == "vfdb" and self.option("align_database"):
                self.vfdb_only_database = "T"
                self.run_profile()
        else:
            self.gene_profile_select.on('end', self.run_anno_select)
            self.gene_profile()

    def gene_profile(self):
        self.logger.info("gene_profile_select!")
        options = {
            "select_genes": self.option("geneset_table"),
            "origin_table": self.option("gene_profile"),
            "filter_zero": self.option("filter_zero") #用于过滤total丰度为0 add by qingchen.zhang@20190315
        }
        if self.option("samples"):
            options["samples"] = self.option("samples")
        self.gene_profile_select.set_options(options)
        self.gene_profile_select.run()

    def creat_lowest_file(self):
        self.logger.info("start run creat lowest_file!")
        level = self.lowest_names[self.option("database")]
        gene_anno = self.fun_filter_tool.option("sel_gene_anno")
        if self.option("gene_filter") == "T":
            self.final_gene_profile = self.gene_profile_select.option("select_table")
        options = {
            "anno_file": gene_anno,
            "gene_list": self.option("geneset_table"),
            "gene_profile": self.final_gene_profile,
            "level": level
        }
        if self.option("samples"):
            options["samples"] = self.option("samples")
        self.creat_lowest_tool.set_options(options)
        self.creat_lowest_tool.run()

    def abundance_select(self):
        self.logger.info("start run abundance_select!")
        abu_num = self.option("abu_num")
        abu_proportion = self.option("abu_proportion")
        if self.gene_anno_sel == "T":
            self.lowest_file = self.creat_lowest_tool.option("outprofile")
        else:
            self.lowest_file = self.option("lowest_level_profile")
        options = {
            "lowest_level_profile": self.lowest_file,
            "abu_num": abu_num,
            "abu_proportion": abu_proportion,
        }
        self.abu_filter_tool.set_options(options)
        self.abu_filter_tool.run()

    def run_anno_select(self):
        '''
        丰度筛选之外的其他筛选
        '''
        self.logger.info("srart function select!")
        options = {
            "gene_anno_table": self.option("gene_anno"),
            "identity": self.option("identity"),
            "align_length": self.option("align_length"),
            "level_select": self.option("level_select"),
            "database": self.option("database"),
        }
        if self.option("gene_filter") == "T":
            options["geneset"] = self.gene_profile_select.option("select_table")
        else:
            options["geneset"] = self.option("geneset_table")
        if "personal_fun" in self.get_option_object().keys() and self.option("personal_fun").is_set:
            options["personal_fun"] = self.option("personal_fun")
        self.fun_filter_tool.set_options(options)
        self.fun_filter_tool.run()

    def anno_select_lowest(self):
        '''
        gene_anno丰度筛选
        '''
        self.logger.info("start run lowest gene_anno!")
        if self.gene_anno_sel == "T" and self.option("abu_filter") == "T":
            self.final_gene_anno = self.fun_filter_tool.option("sel_gene_anno")
            self.new_lowest_file = self.abu_filter_tool.option("abu_select_level")
            self.logger.info("final lowest ok!")
        elif self.gene_anno_sel == "F" and self.option("abu_filter") == "T":
            self.final_gene_anno = self.option("gene_anno")
            self.new_lowest_file = self.abu_filter_tool.option("abu_select_level")
        options = {
            "gene_anno_table": self.final_gene_anno,
            "sel_lowest_level": self.new_lowest_file,
            "database": self.option("database"),
        }
        if self.gene_anno_sel == "F":
            options["geneset"] = self.option("geneset_table")
        self.lowest_select_tool.set_options(options)
        self.lowest_select_tool.run()

    def run_profile(self):
        self.logger.info("start run recalculate profile!")
        if self.option("abu_filter") == "T":
            self.final_gene_anno = self.lowest_select_tool.option("sel_gene_anno")
        elif self.option("abu_filter") == "F":
            self.final_gene_anno = self.fun_filter_tool.option("sel_gene_anno")
        if self.option("gene_filter") == "T":
            self.final_gene_profile = self.gene_profile_select.option("select_table")
        if self.vfdb_only_database == "T":
            self.final_gene_anno = self.option("gene_anno")
        options = {
            "gene_anno_table": self.final_gene_anno,
            "gene_profile": self.final_gene_profile,
            "database": self.option("database"),
            "task_id": self.option("task_id") ## add by qingchen.zhang
        }
        if self.option("database") == "kegg":
            # options["xml_file"] = self.option("xml_file")
            options["origin_anno"] = self.option("gene_anno")
        if self.option("database") == "vfdb":
            options["vfdb_type"] = self.option("align_database")
        if self.option("group_table").is_set:
            options["group"] = self.option("group_table")
        if self.option("database") == "go":
            options["go1234level_out"] = self.option("go1234level_out")
        self.recal_abu_tool.set_options(options)
        self.recal_abu_tool.run()

    def end(self):
        super(MgFunSelectModule, self).end()
