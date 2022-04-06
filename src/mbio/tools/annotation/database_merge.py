# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# modified 20180107
# last modified by guhaidong 20180408 增加内存
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess


class DatabaseMergeAgent(Agent):
    """
    宏基因组注释总览各数据库gene_*_anno结果合并
    """
    def __init__(self, parent):
        super(DatabaseMergeAgent, self).__init__(parent)
        options = [
            {"name": "gene_length_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_nr_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_cog_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_kegg_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_cazy_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_ardb_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_card_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_vfdb_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_overview", "type": "outfile", "format": "sequence.profile_table"},
            ###
            {"name": "other", "type":"string", 'default':''},   # 其他分析类型（比如个性化注释）的输入文件,逗号分割
            {"name": "otherf", "type": "string", 'default':''},   # 其他分析类型（比如个性化注释）的输入文件,逗号分割
            {"name": "pre_sum", "type": "infile", "format": "sequence.profile_table"}

        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by guhaidong @ 20100214

    def check_options(self):
        if not self.option("gene_length_table").is_set:
            raise OptionError("必须提供gene_length_table作为输入文件", code="31201301")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"  # 10G增加内存至20G by GHD @20180309
        # tmp_mem = 20 * (self._rerun_time + 1)  # 每次因拼接失败而重运行的内存增加10G by GHD @ 20180408
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('database_merge use memory : ' + self._memory)

    def end(self):
        super(DatabaseMergeAgent, self).end()


class DatabaseMergeTool(Tool):
    def __init__(self, config):
        super(DatabaseMergeTool, self).__init__(config)
        self.python_path = "miniconda2/bin/python"
        self.script_path = self.config.PACKAGE_DIR + "/annotation/mg_annotation/overview_merge.py"

    def run(self):
        super(DatabaseMergeTool, self).run()
        self.run_merge()
        self.set_output()
        self.end()

    def run_merge(self):
        gene_length_table = self.option('gene_length_table').prop['path']
        outfile = os.path.join(self.output_dir, "gene_overview_anno.xls")
        cmd = '{} {} -gl {} -o {}'.format(self.python_path, self.script_path, gene_length_table, outfile)
        if self.option("gene_nr_anno").is_set:
            cmd += " -nr " + self.option("gene_nr_anno").prop['path']
        if self.option("gene_cog_anno").is_set:
            cmd += " -cog " + self.option("gene_cog_anno").prop['path']
        if self.option("gene_kegg_anno").is_set:
            cmd += " -kegg " + self.option("gene_kegg_anno").prop['path']
        if self.option("gene_cazy_anno").is_set:
            cmd += " -cazy " + self.option("gene_cazy_anno").prop['path']
        if self.option("gene_ardb_anno").is_set:
            cmd += " -ardb " + self.option("gene_ardb_anno").prop['path']
        if self.option("gene_card_anno").is_set:
            cmd += " -card " + self.option("gene_card_anno").prop['path']
        if self.option("gene_vfdb_anno").is_set:
            cmd += " -vfdb " + self.option("gene_vfdb_anno").prop['path']
        ##zouguanqing
        if self.option('other') !='' and not self.option('pre_sum').is_set:
            cmd += ' -other {} -otherf {} '.format(self.option('other'), self.option('otherf'))
        elif self.option('other') !='' and self.option('pre_sum').is_set:
            pre_sum = self.option('pre_sum').prop['path']
            #all_sum = self.output_dir + '/gene_overview_anno_all.xls'
            cmd += ' -other {} -otherf {} -add {} '.format(self.option('other'), self.option('otherf'), pre_sum)
        ###
        self.logger.info("start merge overview anno")
        command = self.add_command("overview_anno", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("merge overview anno successed")
        elif command.return_code in [1,-9]:  # 内存错误返回值-9 @ 20180408
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("merge overview anno error", code="31201301")
            self.set_error("merge overview anno error", code="31201302")



    def set_output(self):
        try:
            if os.path.exists(os.path.join(self.output_dir, "gene_overview_anno.xls")):
                self.option("gene_overview",os.path.join(self.output_dir, "gene_overview_anno.xls"))
        except Exception as e:
            self.set_error("设置gene_overview失败--%s", variables=(e), code="31201303")
