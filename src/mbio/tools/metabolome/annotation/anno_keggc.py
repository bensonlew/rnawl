# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.metabolome.anno_keggc import AnnoKeggc
from mbio.packages.metabolome.common import Relation
from  mainapp.models.mongo.metabolome import  Metabolome
from biocluster.config import Config


class AnnoKeggcAgent(Agent):
    """
    代谢组化合物注释
    version: v1.0
    author: guhaidong
    last_modify: 2018.06.15
    """

    def __init__(self, parent):
        super(AnnoKeggcAgent, self).__init__(parent)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_desc"},  # 预处理的metab_desc结果
            {"name": "level_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释层级结果
            {"name": "stat_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释统计结果
            {"name": "database_name", "type": "string", "default" : "CBR"},
            {"name": "database_version", "type": "string", "default" : ""},
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
        ]
        self.add_option(options)
        self.step.add_steps("keggc")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.keggc.start()
        self.step.update()

    def stepfinish(self):
        self.step.keggc.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('metab_table'):
            raise OptionError('必须输入代谢物详情', code="34700501")
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AnnoKeggcAgent, self).end()


class AnnoKeggcTool(Tool):
    def __init__(self, config):
        super(AnnoKeggcTool, self).__init__(config)
        self.level_out = os.path.join(self.work_dir, 'level.xls')
        self.stat_out = os.path.join(self.work_dir, 'stat.xls')
        self.metablome = Metabolome()
        self.metablome._config = Config()
        if self.option("database_version"):
            self.version_value = self.option("database_version")
        elif self.option("task_id") in ['']:## 工作流
            self.version_value = 'v2021.09.18'
        else:
            self.version_value = self.metablome.find_version_from_task(self.option("task_id"))

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoKeggcTool, self).run()
        obj = AnnoKeggc()
        if self.version_value in ['kegg']:
            map_dic ={
                "CBR" : "kegg_compound",
                "BP" : "kegg_compound_br08005",
                "EDC": "kegg_compound_br08006",
                "P" : "kegg_compound_br08007" ,
                "PC" : "kegg_compound_br08003",
                "L" : "kegg_compound_br08002"
            }
        elif self.version_value in ['v94.2', '94.2']:
            map_dic ={
                "CBR" : "kegg_v202007_br08001",
                "BP" : "kegg_v202007_br08005",
                "EDC": "kegg_v202007_br08006",
                "P" : "kegg_v202007_br08007" ,
                "PC" : "kegg_v202007_br08003",
                "L" : "kegg_v202007_br08002"
            }
        else:
            map_dic ={
                "CBR" : "kegg_{}_br08001",
                "BP" : "kegg_{}_br08005",
                "EDC": "kegg_{}_br08006",
                "P" : "kegg_{}_br08007" ,
                "PC" : "kegg_{}_br08003",
                "L" : "kegg_{}_br08002"
            }
            map_dic = {k: v.format(self.version_value) for k, v in map_dic.items()}

        self.logger.info(self.option("metab_table").prop["path"])
        self.logger.info(map_dic)
        try:
            database_name = map_dic[self.option('database_name')]
            obj.run(self.option('metab_table').path, self.level_out, self.stat_out,database_name=database_name)
        except Exception, e:
            self.logger.error(e)
            self.set_error(e, code="34700501")
        self.id_to_name()
        self.set_output()
        self.end()

    def id_to_name(self):
        # 增加一列代谢物名称
        self.level_out_final = os.path.join(self.output_dir, 'level.xls')
        self.stat_out_final = os.path.join(self.output_dir, 'stat.xls')
        table_path = self.option("metab_table").prop["path"]
        self.metab_trans = Relation()
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
        self.metab_trans.add_metabolites_column(id_name_dict, self.level_out, self.level_out_final)
        self.metab_trans.add_metabolites_column(id_name_dict, self.stat_out, self.stat_out_final)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('level_out').set_path(self.level_out_final)
        self.option('stat_out').set_path(self.stat_out_final)
        self.logger.info("设置anno_keggc分析结果目录成功")
