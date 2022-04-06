# coding=utf-8
#__author__ = 'qingchen.zhang' @20191128
import re,os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import pandas as pd
import subprocess


class CorrelationAgent(Agent):
    """
    Expression clustering analysis
    """
    def __init__(self, parent):
        super(CorrelationAgent, self).__init__(parent)
        options = [
            {'name': 'dis_matrix', 'type': 'infile', 'format': 'meta.beta_diversity.distance_matrix'},##用于计算矩阵
            {'name': 'corr_method', 'type': 'string', 'default': 'pearsonr'},
            {"name": "hcluster_method", "type": "string", "default": "average"},
            {"name": "cor_table", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "pvalue_table", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "newicktree", "type": "outfile","format": "meta.beta_diversity.newick_tree"},
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"}, ##用于计算相关性heatmap图
            # 做样本聚类不需要传递group_dict
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('dis_matrix').is_set:
            raise OptionError('必须提供输入距离矩阵表')
        if not self.option('corr_method') in ['pearsonr', 'kendalltau', 'spearmanr']:
            raise OptionError('必须提供计算相关性方法')
        if not self.option('hcluster_method') in ['average', 'single', 'complete']:
            raise OptionError('请提供正确的层级聚类方法')

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('20')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性分析"]
            ])
        super(CorrelationAgent, self).end()


class CorrelationTool(Tool):
    """
    Expression clustering analysis
    """
    def __init__(self, config):
        super(CorrelationTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.cmd_path = os.path.join(
            self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/plot-hcluster_tree.pl')
        self.correlation_path = '/miniconda2/bin/python {}/bac_comp_genome/correlation.py' \
            .format(self.config.PACKAGE_DIR)

    def run_hcluster(self):
        """
        运行plot-hcluster_tree.pl得到层级聚类树
        """
        real_dis_matrix = self.option("dis_matrix").prop['path']
        cmd = self.cmd_path
        cmd += ' -i %s -o %s -m %s' % (
            real_dis_matrix, self.work_dir, self.option('hcluster_method'))
        self.logger.info('运行plot-hcluster_tree.pl程序计算Hcluster')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 hc.cmd.r 文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 hc.cmd.r 文件失败')
            self.set_error('无法生成 hc.cmd.r 文件', code="32702201")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/hc.cmd.r' % self.work_dir, shell=True)
            self.logger.info('生成树文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成树文件失败')
            self.set_error('无法生成树文件', code="32702202")
        filename = self.work_dir + '/hcluster_tree_' + \
            os.path.basename(real_dis_matrix) + '_' + self.option('hcluster_method') + '.tre'
        linkfile = self.output_dir + '/hcluster.tre'
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename, linkfile)
        self.option("newicktree", linkfile)


    def run_correlation(self):
        """
        用丰度表计算相关性系数表和pvalue
        """
        cmd = self.correlation_path
        cmd += " %s %s" % (self.option("otu_table").prop['path'], self.option("corr_method"))
        self.logger.info('运行correlation.py计算correlation')
        self.logger.info(cmd)
        command = self.add_command("correlation", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format("correlation"))
        else:
            self.set_error("%s Failed. >>>%s" %("correlation", cmd))


    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        tree_file = self.output_dir + '/hcluster.tre'
        if os.path.exists(tree_file):
            self.logger.info("已经生成了正确的树文件")
        else:
            self.set_error("未能正确的生成树文件")
        corr_file = os.path.join(self.work_dir, 'correlation_{}.xls'.format(self.option("corr_method")))
        linkfile = os.path.join(self.output_dir, 'correlation_{}.xls'.format(self.option("corr_method")))
        if os.path.exists(linkfile):
            os.remove(linkfile)
        if os.path.exists(corr_file):
            os.link(corr_file, linkfile)
        pvalue_file = os.path.join(self.work_dir, "{}_pvalue.xls".format(self.option("corr_method")))
        linkfile2 = os.path.join(self.output_dir, "{}_pvalue.xls".format(self.option("corr_method")))
        if os.path.exists(linkfile2):
            os.remove(linkfile2)
        if os.path.exists(pvalue_file):
            os.link(pvalue_file, linkfile2)
        self.logger.info("设置结果文件完成")


    def run(self):
        super(CorrelationTool, self).run()
        self.run_hcluster()
        self.run_correlation()
        self.set_output()
        self.end()
