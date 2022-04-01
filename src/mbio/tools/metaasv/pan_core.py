# -*- coding: utf-8 -*-
# __author__ = qingchen.zhang'
import os
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.meta.otu.pan_core_otu import pan_core


class PanCoreAgent(Agent):
    """
    需要pan_core_otu.py的package包
    需要R软件包
    计算pan、core
    """
    def __init__(self, parent):
        super(PanCoreAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir"},  # 输入的OTU文件
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 输入的group表格
            {"name": "level", "type": "string", "default": "asv"},  # 物种水平
            {"name": "pan_otu_table", "type": "outfile", "format": "meta.otu.pan_core_table"},  # 输出的pan_otu表格
            {"name": "core_otu_table", "type": "outfile", "format": "meta.otu.pan_core_table"}]  # 输出的core_otu表格
        self.add_option(options)
        self.step.add_steps("create_pan_core")
        self.on('start', self.start_pan_core)
        self.on('end', self.end_pan_core)

    def start_pan_core(self):
        self.step.create_pan_core.start()
        self.step.update()

    def end_pan_core(self):
        self.step.create_pan_core.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("in_otu_table").is_set:
            raise OptionError("参数otu_table不能为空")
        if self.option("level") not in ['asv', 'domain', 'kindom', 'phylum', 'class',
                                        'order', 'family', 'genus', 'species']:
            raise OptionError("请选择正确的分类水平")
        return True

    def end(self):
        super(PanCoreAgent, self).end()

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 3
        self._memory = '5G'


class PanCoreTool(Tool):
    def __init__(self, config):
        super(PanCoreTool, self).__init__(config)
        self.R_path = os.path.join(Config().SOFTWARE_DIR, "program/R-3.3.1/bin/R")
        self._version = 1.0

    def _create_pan_core(self):
        """
        用脚本pan_core_otu.py,输出pan_otu表格
        """
        self.logger.info("开始生成R脚本")
        otu_path = ""
        if self.option("in_otu_table").format == "meta.otu.tax_summary_dir":
            otu_path = self.option("in_otu_table").get_table(self.option("level"))
        else:
            otu_path = self.option("in_otu_table").prop['path']

        # 检测otu表行数
        if self.option("group_table").is_set:
            group_table = self.option("group_table").prop['path']
            pan_otu = pan_core(otu_path, "pan", group_table, self.work_dir)
            core_otu = pan_core(otu_path, "core", group_table, self.work_dir)
        else:
            pan_otu = pan_core(otu_path, "pan", "none", self.work_dir)
            core_otu = pan_core(otu_path, "core", "none", self.work_dir)

        self.logger.info("R脚本生成完毕")
        self.logger.info("运行R,生成表格文件")
        try:
            pan_cmd = self.R_path + " --restore --no-save < " + pan_otu
            core_cmd = self.R_path + " --restore --no-save < " + core_otu
            subprocess.check_call(pan_cmd, shell=True)
            subprocess.check_call(core_cmd, shell=True)
            self.logger.info("表格生成完毕")
        except subprocess.CalledProcessError:
            self.set_error("表格生成失败")
            raise Exception("运行R出错")
        tmp_pan = os.path.join(os.path.dirname(pan_otu), "pan.richness.xls")
        tmp_core = os.path.join(os.path.dirname(core_otu), "core.richness.xls")
        pan_dir = os.path.join(self.work_dir, "output", "Pan.richness.xls")
        core_dir = os.path.join(self.work_dir, "output", "Core.richness.xls")
        self.del_na(tmp_pan, pan_dir)
        self.del_na(tmp_core, core_dir)
        self.option("pan_otu_table").set_path(pan_dir)
        self.option("core_otu_table").set_path(core_dir)

    def del_na(self, in_path, out_path):
        """
        输出的pan和core表格存在多余的NA列的情况，删除这些列
        """
        none_na = list()
        with open(in_path, "rb") as r:
            r.next()
            for line in r:
                line = line.rstrip().split("\t")
                for i in range(1, len(line)):
                    if line[i] != "NA":
                        none_na.append(i)
                    else:
                        break
        max_len = max(none_na) + 1
        with open(in_path, "rb") as r, open(out_path, "wb") as w:
            for line in r:
                line = line.rstrip().split("\t")
                tmp_list = list()
                for i in range(max_len):
                    tmp_list.append(line[i])
                w.write("\t".join(tmp_list) + "\n")

    def run(self):
        """
        运行
        """
        super(PanCoreTool, self).run()
        self._create_pan_core()
        self.logger.info("程序完成，即将退出")
        self.end()
