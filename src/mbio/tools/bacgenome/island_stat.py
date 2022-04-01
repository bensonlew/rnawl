#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.bacgenome.common import sum_stat


class IslandStatAgent(Agent):
    """
    生成基因组岛文件
    version 1.0
    author: gaohao
    last_modify: 2018.04.13
    """

    def __init__(self, parent):
        super(IslandStatAgent, self).__init__(parent)
        options = [
            {"name": "diomb", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "islander", "type": "infile", "format": "bacgenome.island"},
            {"name": "anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "sample_name", "type": "string"},
            {"name": "analysis", "type": "string", "default": "uncomplete"}  ###流程分析模式complete，uncomplete
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('diomb').is_set:
            raise OptionError("请设置diomb软件的结果文件！", code="31402401")
        if not self.option('islander').is_set:
            raise OptionError("请设置islander软件的结果文件！", code="31402402")
        if not self.option('anno').is_set:
            raise OptionError("请设置基因组注释总览表！", code="31402403")
        if not self.option('sample_name'):
            raise OptionError("请设置样品名！", code="31402404")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["", "", ""],
            ["", "", ""]
        ])
        super(IslandStatAgent, self).end()


class IslandStatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(IslandStatTool, self).__init__(config)
        self.diomb =self.option('diomb').prop['path']
        self.islander = self.option('islander').prop['path']
        self.anno = self.option('anno').prop['path']
        self.sample = self.option('sample_name')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"



    def run_island(self):
        cmd = '{} {}combin_island.pl {} {} {}'.format(self.perl_path, self.perl_script,self.option('analysis'),self.diomb,self.islander)
        self.logger.info(cmd)
        command = self.add_command("run_island", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_island运行完成")
        else:
            self.set_error("run_island运行出错!", code="31402401")

    def run_island_stat(self):
        out = self.work_dir + '/all.island.xls'
        cmd = '{} {}get_island_gene.pl {} {} {} {}'.format(self.perl_path, self.perl_script, self.option('analysis'),
                                                           out, self.anno,
                                                           self.sample)
        self.logger.info(cmd)
        command = self.add_command("run_island_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_island_stat运行完成")
        else:
            self.set_error("run_island_stat运行出错!", code="31402402")

    def set_output(self):
        file =self.output_dir + '/' + self.sample + '_GI_summary.xls'
        file2 = self.output_dir + '/' + self.sample + '_GI_detail.xls'
        for fil in [file,file2]:
            if os.path.exists(fil):
               os.remove(fil)
        if os.path.exists(self.work_dir + '/' + self.sample + '.GI_summary.xls'):
            num =self.get_num(self.work_dir + '/' + self.sample + '.GI_summary.xls')
            if num >1:
                os.link(self.work_dir + '/' + self.sample + '.GI_detail.xls',self.output_dir + '/' + self.sample + '_GI_detail.xls')
                os.link(self.work_dir + '/' + self.sample + '.GI_summary.xls', self.output_dir + '/' + self.sample + '_GI_summary.xls')
                stat_num_out = os.path.join(self.output_dir, "sample_stat.xls")
                sum_stat(file,"Location","GI No.",stat_num_out,stat_method="count",ocname="GI No.")

    def run(self):
        """
        运行
        """
        super(IslandStatTool, self).run()
        self.run_island()
        self.run_island_stat()
        self.set_output()
        self.end()

    def get_num(self,file):
        with open(file,'r') as f:
            lines =f.readlines()
            num =len(lines)
        return num
