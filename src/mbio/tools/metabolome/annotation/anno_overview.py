# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.metabolome.common import Relation
from  mainapp.models.mongo.metabolome import  Metabolome
from biocluster.config import Config


class AnnoOverviewAgent(Agent):
    """
    代谢组化合物注释
    version: v1.0
    author: guhaidong
    last_modify: 2018.06.19
    """

    def __init__(self, parent):
        super(AnnoOverviewAgent, self).__init__(parent)
        options = [
            {"name": "metab_pos", "type": "infile", "format": "metabolome.metab_desc"},  # 原始数据的metab_desc结果
            {"name": "metab_neg", "type": "infile", "format": "metabolome.metab_desc"},
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset"},  # 代谢集
            {"name": "work_type", "type": "string", "default": "GC"},  # 工作流类型，"GC"或"LC"
            {"name": "organism", "type": "string", "default": "All"},  # 参考物种
            {"name": "anno_out", "type": "outfile", "format": "metabolome.overview_table"},  # 注释总览结果
            {"name": "ko_out", "type": "outfile", "format": "sequence.profile_table"},  # 通路统计结果
            #{"name": "metab_mix", "type": "infile", "format": "metabolome.metab_desc"}  #LC
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
            {"name": "database_version", "type": "string", "default": ""}, # 用于区分新老工作流
        ]
        self.add_option(options)
        self.step.add_steps("overview")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.overview.start()
        self.step.update()

    def stepfinish(self):
        self.step.overview.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('metab_pos'):
            raise OptionError('必须输入代谢物详情', code="34700401")
        if self.option('work_type') == "LC" and not self.option('metab_neg'):
            raise OptionError('LC流程必须输入阴离子代谢物详情', code="34700402")
        if not self.option('metabset'):
            raise OptionError("必须输入代谢集", code="34700403")
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
        super(AnnoOverviewAgent, self).end()


class AnnoOverviewTool(Tool):
    def __init__(self, config):
        super(AnnoOverviewTool, self).__init__(config)
        self.anno_out = os.path.join(self.output_dir, 'anno.xls')
        self.ko_out = os.path.join(self.output_dir, 'ko.xls')
        self.metablome = Metabolome()
        self.metablome._config = Config()
        self.version_value = ''
        if self.option("database_version"):## 工作流
            self.version_value = self.option("database_version")
        else:# 交互分析
            self.version_value = self.metablome.find_version_from_task(self.option("task_id"))


    def run(self):
        """
        运行
        :return:
        """
        super(AnnoOverviewTool, self).run()
        if self.version_value not in ['kegg']:
            from mbio.packages.metabolome.anno_overview import AnnoOverview
            self.orginisms = self.option('organism').strip("_")
            obj = AnnoOverview(self.option('work_type'), self.orginisms, db_v=self.version_value)
        else:
            from mbio.packages.metabolome.anno_overview2 import AnnoOverview
            self.orginisms = self.option('organism')
            obj = AnnoOverview(self.option('work_type'), self.orginisms)
        metab_table = self.option('metab_pos').path
        if self.option('work_type') == "LC":
            metab_table += ',' + self.option('metab_neg').path
        self.logger.info("DEBUG metab_table %s" % metab_table)
        self.logger.info(self.option('work_type'))
        self.logger.info(self.option('organism'))
        self.logger.info(metab_table)
        self.logger.info(self.option('metabset').path)
        self.logger.info(self.anno_out)
        self.logger.info(self.ko_out)
        obj.run(metab_table, self.option('metabset').path, self.anno_out, self.ko_out, from_mongo=False)
        #try:
        #    obj.run(metab_table, self.option('metabset').path, self.anno_out, self.ko_out, from_mongo=False)
        #except Exception,e:
        #    self.logger.error(e)
        #    self.set_error(e, code="34700401")
        # self.metab_trans = Relation()
        # self.metab_trans.add_oid_fun(self.anno_out, table_path,link_k='metab_id')  ##20190617
        self.logger.info("success")

        self.set_output()
        self.end()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option('anno_out').set_path(self.anno_out)
        self.option('ko_out').set_path(self.ko_out)
        self.logger.info("设置anno_overview分析结果目录成功")
