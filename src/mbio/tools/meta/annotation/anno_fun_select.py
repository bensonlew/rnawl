# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180718
# last modified : guhaidong

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
import pandas as pd

class AnnoFunSelectAgent(Agent):
    """
    宏基因组注释交互分析功能筛选
    """

    def __init__(self, parent):
        super(AnnoFunSelectAgent, self).__init__(parent)
        options = [
            {"name": "gene_anno_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "identity", "type": "float","default": 0.0 },
            {"name": "align_length", "type": "int","default": 0 },
            {"name": "database", "type": "string" },
            {"name": "geneset", "type": "infile", "format": "sequence.profile_table"},
            {"name": "level_select", "type": "string","default":"all"},
            {"name": "personal_fun", "type": "infile", "format": "sequence.profile_table"}, # 个性化功能集使用
            {"name": "sel_lowest_level", "type": "infile", "format": "sequence.profile_table"},
            {"name": "sel_gene_anno", "type": "outfile", "format": "sequence.profile_table"},

        ]
        self.add_option(options)
        self._memory_increase_step = 40  # modified by GHD @ 20180718

    def check_options(self):
        if not self.option("gene_anno_table").is_set:
            raise OptionError("必须设置基因丰度文件", code="32700501")
        if self.option("identity") != 0.0:
            if not self.option("identity") >= 0 and self.option("identity") <= 1:
                raise OptionError("identity必须在0-1之间", code="32700502")
        if self.option("align_length") != 0:
            if not isinstance(self.option("align_length"), int):
                raise OptionError("align_length必须为整数", code="32700503")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '8G' # 改回 by guhaidong @ 20180427
        # memory = 8 + 10 * self._rerun_time  # 每次重运行增加5G内存 by guhaidong @ 20180417
        # self._memory = "%sG" % memory

    def end(self):
        super(AnnoFunSelectAgent, self).end()


class AnnoFunSelectTool(Tool):
    def __init__(self, config):
        super(AnnoFunSelectTool, self).__init__(config)
        self.python_path =  "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/function_select.py'

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoFunSelectTool, self).run()
        if self.option("level_select") != "all" or \
                ("personal_fun" in self.get_option_object().keys() and self.option("personal_fun").is_set):
            self.to_level_file()
        self.run_select_abu()
        self.set_output()
        self.end()

    def run_select_abu(self):
        self.logger.info("start annotation function select")
        if self.option("level_select") != "all" or self.option("personal_fun").is_set:
            level_select_file = self.work_dir + "/select_level.xls"
        else:
            level_select_file = "all"
        cmd = self.python_path  + ' {} -i {} -l {} -database {} -o {}'. \
            format(self.script, self.option("gene_anno_table").prop['path'],
                   level_select_file, self.option("database"), self.output_dir + "/function_select_anno.xls")
        if self.option("identity") != 0.0:
            cmd +=  " -iden {}".format( self.option("identity"))
        if self.option("align_length") != 0:
            cmd += " -len {}".format(self.option("align_length"))
        if self.option("sel_lowest_level").is_set:
            sel_lowest_level_file = self.option("sel_lowest_level").prop['path']
            sel_lowest_level_table = pd.read_table(sel_lowest_level_file, sep='\t', header=0)
            if len(sel_lowest_level_table) > 0:
                cmd += " -s {}".format(sel_lowest_level_file)
            else:
                self.set_error("丰度筛选结果为无！", code="32700501")
        if self.option("geneset").is_set:
            geneset = self.option("geneset").prop['path']
            cmd += " -g {}".format(geneset)
        command = self.add_command('function_select', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("annotation function select succeed")
        elif command.return_code in [1, -9]:  # add memory limit by guhaidong @ 20180417 modified @ 20180711
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("annotation function select failed", code="32700502")
            self.set_error("annotation function select failed", code="32700504")

    def set_output(self):
        outputfile = self.output_dir + "/function_select_anno.xls"
        table = pd.read_table(outputfile, sep='\t', header=0)
        if len(table) > 0:
            self.option("sel_gene_anno", outputfile)
        else:
            self.set_error("您选择的数据结果为空！", code="32700503")

    def to_level_file(self):
        if self.option("level_select") != "all":
            level_select = eval(self.option("level_select"))
            self.logger.info(level_select)
            metagenomic = MetagenomicController()
            origin_table = pd.DataFrame()
            for each in level_select:
                mytype = each['type']
                level_id = each['level_id']
                if int(level_id) == 14:
                    level = "Pathway"   ## KEGG特殊处理，注释筛选pathway id，其他分析显示level3名称
                else:
                    level = metagenomic.level_id(level_id)
                name_id = each['name_id']
                if int(level_id) in [24, 71]:
                    name = name_id
                else:
                    name = metagenomic.level_convert(name_id, level_id)
                if int(level_id) == 10:
                    name = name.split(':')[0]
                tmp_list = [str(mytype), level, name]
                tmp_dataframe = pd.DataFrame(tmp_list).T
                origin_table = pd.concat([origin_table,tmp_dataframe])
        if "personal_fun" in self.get_option_object().keys() and self.option("personal_fun").is_set:
            print "add personal"
            personal_fun_path = self.option("personal_fun").prop["path"]
            personal_table = pd.read_table(personal_fun_path,sep="\t",header=None)
            if self.option("level_select") != "all":
                origin_table = pd.concat([origin_table,personal_table])
            else:
                origin_table = personal_table
            print origin_table
        origin_table = origin_table.drop_duplicates()
        origin_table.to_csv(self.work_dir + "/select_level.xls", sep = "\t", index=False, header=False)




