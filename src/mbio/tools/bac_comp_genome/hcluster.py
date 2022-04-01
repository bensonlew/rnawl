# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import subprocess
from biocluster.core.exceptions import OptionError


class HclusterAgent(Agent):
    """
    脚本plot-hcluster_tree.pl
    version v1.0
    author: shenghe
    last_modified:2016.3.24
    """

    def __init__(self, parent):
        super(HclusterAgent, self).__init__(parent)
        options = [
            {"name": "dis_matrix", "type": "infile", "format": "sequence.profile_table"},
            {"name": "newicktree", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},
            {"name": "linkage", "type": "string", "default": "average"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('dis_matrix').is_set:
            raise OptionError('必须提供输入距离矩阵表')
        else:
            self.option('dis_matrix').check()
        if self.option('linkage') not in ['average', 'single', 'complete']:
            raise OptionError('错误的层级聚类方式：%s', variables=(self.option('linkage')))

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "层次聚类结果目录"],
            ["./hcluster.tre", "tre", "层次聚类树"]
        ])
        super(HclusterAgent, self).end()


class HclusterTool(Tool):
    def __init__(self, config):
        super(HclusterTool, self).__init__(config)
        self._version = 'v2.1-20140214'  # plot-hcluster_tree.pl版本
        ld_path = self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64:" + self.config.SOFTWARE_DIR + "/library/lib"
        path = self.config.SOFTWARE_DIR + "gcc/5.1.0/bin"
        self.set_environ(PATH= path, LD_LIBRARY_PATH=ld_path)
        self.cmd_path = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/plot-hcluster_tree.pl')

    def run(self):
        """
        运行
        """
        super(HclusterTool, self).run()
        self.run_hcluster()
        self.end()

    def run_hcluster(self):
        """
        运行plot-hcluster_tree.pl
        """
        real_dis_matrix = self.option("dis_matrix").prop['path']
        cmd = self.cmd_path
        cmd += ' -i %s -o %s -m %s' % (
            real_dis_matrix, self.work_dir, self.option('linkage'))
        self.logger.info('运行plot-hcluster_tree.pl程序计算Hcluster')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 hc.cmd.r 文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 hc.cmd.r 文件失败')
            self.set_error('无法生成 hc.cmd.r 文件')
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1_gcc5.1/bin/R --restore --no-save < %s/hc.cmd.r' % self.work_dir, shell=True)
            self.logger.info('生成树文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成树文件失败')
            self.set_error('无法生成树文件')
        filename = self.work_dir + '/hcluster_tree_' + \
            os.path.basename(real_dis_matrix) + '_' + self.option('linkage') + '.tre'
        linkfile = self.output_dir + '/hcluster.tre'
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename , linkfile)
        self.option('newicktree', linkfile)