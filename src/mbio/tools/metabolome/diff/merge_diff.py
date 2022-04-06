# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.06.14

import os, re, subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.common import Relation

class MergeDiffAgent(Agent):
    """
    合并差异检验和pls分析结果并筛选条件创建代谢集
    """

    def __init__(self, parent):
        super(MergeDiffAgent, self).__init__(parent)
        options = [
            {'name': 'test_dir', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},  # 差异检验结果文件夹
            {"name": "pls_dir", "type": "infile", "format": "annotation.mg_anno_dir"},  # pls分析差异结果文件夹
            {'name': 'metab_trans', 'type': 'infile', "format": "sequence.profile_table"},# 转化id用文件。 提取用代谢物名称的meta id用。
            {"name": "sub_group", "type": "string"},  # 所对应的分组
            #{"name": "metab_des", "type": "string"},  #
            {'name': 'metabset', 'type': 'bool', 'default': False},  # 是否创建代谢集
            {'name': 'mul_metabset', 'type': 'outfile', 'format': 'metabolome.mul_metabset'},  # venn图用metab_list
            {'name': 'merge_dir', 'type': 'outfile', 'format': 'annotation.mg_anno_dir'},  #
            {'name': 'filter_k', 'type':'string','default':''},
            {"name": "filter_t","type":"string"},
            {"name": "filter_v", "type": "string"}

        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("test_dir").is_set:
            raise OptionError("请传入差异检验结果文件夹！", code="34700201")
        if not self.option("pls_dir").is_set:
            raise OptionError("请传入pls分析结果文件夹！", code="34700202")
        if not self.option("sub_group"):
            raise OptionError("请传入需要合并的分组名！", code="34700203")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        super(MergeDiffAgent, self).end()


class MergeDiffTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(MergeDiffTool, self).__init__(config)
        self.python_path =  "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/metabolome/scripts/merge_table.py'

    def run(self):
        """
        运行
        """
        super(MergeDiffTool, self).run()
        self.logger.info("开始运行命令！")
        self.run_merge()
        self.set_output()

    def run_merge(self):
        groups = str(self.option("sub_group"))
        test_dir = self.option("test_dir").prop["path"]
        pls_dir = self.option("pls_dir").prop["path"]
        cmd = self.python_path + ' {} -d1 {} -d2 {} -name {} -o {}'.format(self.script, test_dir, pls_dir, groups,
                                                                        self.output_dir)
        if self.option("metabset"):
            cmd += " --ms "
            cmd += ' -des {}'.format(self.option('metab_trans').path)
            if self.option('filter_k'):
                cmd += ' -filter_k %s'%(self.option('filter_k'))
                cmd += ' -filter_t %s'%(self.option('filter_t'))
                cmd += ' -filter_v %s'%(self.option('filter_v'))

        self.logger.info(cmd)
        command = self.add_command('merge_table', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("merge succeed")
        else:
            self.set_error("merge failed", code="34700201")
            raise Exception("merge failed")

    def set_output(self):
        self.logger.info("set output")
        metabset_dir = os.path.join(self.output_dir,"Metabset/mul.metabset.list.xls")
        if os.path.exists(metabset_dir):
            self.option('mul_metabset',os.path.join(self.output_dir,"Metabset/mul.metabset.list.xls"))
        self.logger.info("进行id转化...")
        diffexp_dir = os.path.join(self.output_dir,"tmp_DiffStat")
        new_dir = os.path.join(self.output_dir,"DiffStat")
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)
        table_path = self.option("metab_trans").prop["path"]
        #table_path = table_path.replace("metab_abund.txt","metab_desc.txt")
        self.metab_trans = Relation()
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
        diffexp_files = os.listdir(diffexp_dir)
        for eachfile in diffexp_files:
            self.trans_data(eachfile, eachfile, map_table, id_name_dict)


        self.end()

    def trans_data(self, oldfile, newfile, map_table, id_name_dict):
        diffexp_dir = os.path.join(self.output_dir,"tmp_DiffStat")
        new_dir = os.path.join(self.output_dir,"DiffStat")
        oldfile = os.path.join(diffexp_dir,oldfile)
        newfile = os.path.join(new_dir,newfile)
        self.metab_trans.get_trans_file(map_table, id_name_dict, oldfile, newfile)
