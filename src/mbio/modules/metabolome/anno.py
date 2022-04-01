# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class AnnoModule(Module):
    """
    代谢集注释模块
    author: guhaidong
    last_modify: 2018.06.20
    """

    def __init__(self, work_id):
        super(AnnoModule, self).__init__(work_id)
        option = [
            {"name": "anno_list", "type": "string", "default": "keggc,keggp,overveiw"},  # 要做什么注释
            {"name": "metab_table", "type": "infile", "format": "metabolome.metab_desc"},  # 预处理的metab_desc结果
            {"name": "org_metab_pos", "type": "infile", "format": "metabolome.metab_desc"}, # 用于总览表的原始metab_desc
            {"name": "org_metab_neg", "type": "infile", "format": "metabolome.metab_desc"},
            {"name": "organism", "type": "string", "default": "False"},  # 参考物种
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset"},  # 代谢集
            {"name": "work_type", "type": "string", "default": "GC"},  # 工作流类型，"GC"或"LC"
            {"name": "keggc_level", "type": "outfile", "format": "sequence.profile_table"}, # 化合物注释结果
            {"name": "keggc_stat", "type": "outfile", "format": "sequence.profile_table"},  # 化合物注释结果
            {"name": "keggp_level", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "keggp_stat", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "overview_out", "type": "outfile", "format": "metabolome.overview_table"},
            {"name": "overview_ko", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "task_id", "type": "string"}, # 用于区分新老工作流
            {"name": "database_version", "type": "string"},
        ]
        self.add_option(option)
        self.run_tools = []  # run_metagene并行运行tools列表
        self.anno_list = []
        self.keggc_tool = self.add_tool('metabolome.annotation.anno_keggc')
        self.keggp_tool = self.add_tool('metabolome.annotation.anno_keggp')
        self.overview_tool = self.add_tool('metabolome.annotation.anno_overview')
        self.step.add_steps("keggp", "keggc", "overview")

    def check_options(self):
        """
        检查参数
        :return:
        """

        self.anno_list = self.option('anno_list').split(',')
        for anno in self.anno_list:
            if anno not in ["keggc", "keggp", "overview"]:
                raise OptionError("注释list参数不正确:%s" , variables=(anno), code="24700101")
        if "keggc" in self.anno_list or "keggp" in self.anno_list:
            if not self.option('metab_table'):
                raise OptionError('必须输入数据表', code="24700102")
        if "overview" in self.anno_list:
            if not self.option('metabset'):
                raise OptionError("注释总览必须输入代谢集", code="24700103")
            if not self.option('org_metab_pos'):
                raise OptionError("必须输入质控前的阳离子数据表", code="24700104")
            if self.option('work_type') == "GC" and not self.option("org_metab_neg"):
                raise OptionError("必须输入质控前的阴离子数据表", code="24700105")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_keggc(self):
        opts = {
            'metab_table': self.option('metab_table'),
            'database_version': self.option("database_version")
        }
        self.keggc_tool.set_options(opts)
        self.keggc_tool.run()

    def run_keggp(self):
        self.logger.info("database_version: {}".format(self.option("database_version")))
        opts = {
            'metab_table': self.option('metab_table'),
            'organism': self.option('organism'),
            'database_version': self.option("database_version")
        }
        self.keggp_tool.set_options(opts)
        self.keggp_tool.run()

    def run_overview(self):
        opts = {
            'metab_pos': self.option('org_metab_pos'),
            'metabset': self.option('metabset'),
            'work_type': self.option('work_type'),
            'organism': self.option('organism'),
            'database_version': self.option("database_version")
        }
        if self.option('work_type') == "LC":
            opts['metab_neg'] = self.option('org_metab_neg')
        self.overview_tool.set_options(opts)
        self.overview_tool.run()

    def run(self):
        """
        运行module
        :return:
        """
        super(AnnoModule, self).run()
        for anno in self.anno_list:
            if anno == "keggc":
                self.run_tools.append(self.keggc_tool)
            elif anno == "keggp":
                self.run_tools.append(self.keggp_tool)
            elif anno == "overview":
                self.run_tools.append(self.overview_tool)
        self.on_rely(self.run_tools, self.set_output)
        for anno in self.anno_list:
            if anno == "keggc":
                self.run_keggc()
            elif anno == "keggp":
                self.run_keggp()
            elif anno == "overview":
                self.run_overview()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for anno in self.anno_list:
            if anno == "keggc":
                self.option('keggc_level', self.keggc_tool.option('level_out'))
                self.option('keggc_stat', self.keggc_tool.option('stat_out'))
            elif anno == "keggp":
                self.option('keggp_level', self.keggp_tool.option('level_out'))
                self.option('keggp_stat', self.keggp_tool.option('stat_out'))
            elif anno == "overview":
                self.option('overview_out', self.overview_tool.option('anno_out'))
                self.option('overview_ko', self.overview_tool.option('ko_out'))
        self.logger.info("设置注释结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AnnoModule, self).end()
