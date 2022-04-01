# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# __last_modify__ = '2019/4'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class DiffKeggAgent(Agent):
    """
    version 1.0
    """

    def __init__(self, parent):
        super(DiffKeggAgent, self).__init__(parent)
        options = [
            {"name": "pathway_id", "type": "string"},
            {"name": "mark_dir", "type":"string"},
            {"name": "k_list","type":"string"},

        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option('pathway_id'):
            raise OptionError('缺失pathway_id参数')


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'


class DiffKeggTool(Tool):
    def __init__(self, config):
        super(DiffKeggTool, self).__init__(config)
        self.r_path = '/program/R-3.3.1/bin/Rscript'
        self.pic_path = self.config.PACKAGE_DIR + "/bacgenome/pathway_mark.py"
        self.python_path = "/program/Python/bin/python"
        self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"

    def w_pathway_file(self):
        with open('pathway_id.list', 'w') as fw:
            ko = self.option('pathway_id').replace('map','ko')
            fw.write('#pathway_id K_list\n')
            fw.write(' '.join([ko, self.option('k_list')])+'\n')


    def run_pic(self):
        """
        description
        :return:
        """
        self.logger.info("start output kegg img")
        pathway_file = self.work_dir + "/pathway_id.list"


        cmd = "{} {}  -o {} -p {}   -html {} ".format(
            self.python_path, self.pic_path, self.output_dir, pathway_file, self.html_path )

        command = self.add_command("output_kegg_pathway_img", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("output_kegg_pathway_img succeed")
        else:
            self.set_error("output kegg pathway img failed")


    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        self.option('mark_dir', self.output_dir + '/pathway_img')

    def run(self):
        super(DiffKeggTool, self).run()
        self.w_pathway_file()
        self.run_pic()
        self.set_output()
        self.end()