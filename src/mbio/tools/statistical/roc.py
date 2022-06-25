# -*- coding: utf-8 -*-
# __author__ = "gaohao"

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import string
import types
import subprocess
from biocluster.core.exceptions import OptionError
from mbio.packages.statistical.normalization import Normalization
import pandas as pd


class RocAgent(Agent):
    """
    计算roc，需要roc_plus_gh.pl
    author: gaohao
    last_modifued:2018.09.29
    """

    def __init__(self, parent):
        super(RocAgent, self).__init__(parent)
        options = [
            {"name": "mode", "type": "int", "default": 1},# 1.丰度表的前n物种计算模式;2.根据物种名提取特殊物种丰度计算模式（目前页面不存在2模式）
            {"name": "abu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method", "type": "string", "default": "sum"},
            {"name": "confidence_interval", "type": "string", "default":"0.95"},
            {"name": "top_n", "type": "int", "default": 100},
            {"name": "norm_method", "type": "string", "default": ""},
            {"name": "env_labs", "type": "string", "default": ''},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('abu_table') and self.option('mode') in [1]:
            raise OptionError('abu_table file must be provided !', code="34101501")
        if self.option('mode') in [1]:
            self.option('abu_table').get_info()
            if self.option('abu_table').prop['sample_num'] < 2:
                raise OptionError('The number of samples in the abu_table table is less than 2, and ROC calculation is not possible', code="34101502")
            if not self.option('group_table').is_set:
                raise OptionError("group_table file must be provided !", code="34101503")
            if self.option('method') != '':
                if self.option('method') not in ['sum', 'average', 'median']:
                    raise OptionError("Abundance calculation method can only choose one of sum, average, median", code="34101504")

    def set_resource(self):
        """
        设置内存和CPU
        """
        self._cpu = 10
        self._memory = '15G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ROC分析结果目录"],
            ["./roc_curve.xls", "xls", "ROC受试者工作特征曲线数据"],
            ["./roc_auc.xls", "xls", "ROC受试者工作特征曲线-AUC VALUE"],
            ["./roc_auc_smooth.xls", "xls", "ROC受试者工作特征曲线-AUC VALUE"],
            ["./roc_curve_smooth.xls", "xls", "ROC受试者工作特征曲线平滑处理后数据"],
            ["./roc_interval.xls", "xls", "ROC受试者工作特征曲线的可信区间"],
            ["./best_loc.xls", "xls", "ROC受试者工作特征曲线最佳cut off点"]
        ])
        super(RocAgent, self).end()


class RocTool(Tool):
    def __init__(self, config):
        super(RocTool, self).__init__(config)
        self._version = '1.0.1'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        """
        super(RocTool, self).run()
        self.run_roc()
        self.set_output()
        self.end()

    def data_process(self):
        '''
        当没有top_n==-1及不设置top时 top_n = 丰度表的所有行
        当配置设置了top_n且参数中有 env_labs时，先去除包含env_labs的行；
           再取top_n，然后和env_labs的行合并，最后top_n = top_n + len(env_labs)
        '''
        import commands
        abu_table = self.option('abu_table').path
        if self.option('top_n') == -1:
            cmd = '/usr/bin/wc -l ' + self.option('abu_table').path
            status, output = commands.getstatusoutput(cmd)
            lines = int(output.split()[0]) - 1
            self.option('top_n', lines)
        elif self.option('env_labs'):
            top_n = self.option('top_n')
            env_labs = self.option('env_labs').split(',')
            tb = pd.read_csv(abu_table, sep='\t')
            first = tb.columns[0]
            envs = tb[tb[first].isin(env_labs)]
            otu = tb[tb[first].isin(env_labs).isin([False])]
            if self.option('method') == 'sum':
                filters = otu.sum(axis=1).sort_values(ascending=False)
            elif self.option('method') == 'average':
                filters = otu.mean(axis=1).sort_values(ascending=False)
            elif self.option('method') == 'median':
                filters = otu.median(axis=1).sort_values(ascending=False)
            otu = otu.head(top_n)
            tb = pd.concat([otu, envs])
            print(top_n)
            print(len(otu.index))
            print(len(envs.index))
            abu_table = os.path.join(self.work_dir, 'new_abu.xls')
            tb.to_csv(abu_table, sep='\t', index=False)
            self.option('top_n', top_n + len(env_labs))

        if self.option('norm_method'):
            new_abu = os.path.join(self.work_dir, 'normalized.xls')
            Normalization(abu_table, norm_meth=self.option('norm_method'), out=new_abu).run()
            abu_table = new_abu
        return abu_table

    def run_roc(self):
        """
        运行roc_plus_gh.pl
        """
        abu_table = self.data_process()
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.config.PACKAGE_DIR + '/metagenomic/scripts/roc_plus_gh.pl '
        cmd += '-o %s ' % (self.work_dir + '/ROC/')
        if not os.path.exists(self.work_dir + '/ROC/'):
            os.mkdir(self.work_dir + '/ROC/')
        if self.option('mode') in [1]:
            cmd += '-i %s ' % (abu_table)
        cmd += '-mode %d ' % (self.option('mode'))
        cmd += '-group %s ' % (self.option('group_table').prop['path'])
        if self.option('method'):
            cmd += '-method %s ' % (self.option('method'))
        if self.option('mode') == 1:
            cmd += '-n %d ' % (self.option('top_n'))
        if self.option('confidence_interval'):
            cmd += '-conf_level %s ' % (self.option('confidence_interval'))
        self.logger.info('开始运行roc_plus_gh.pl计算ROC相关数据')
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 roc.cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 roc.cmd.r 失败')
            self.set_error('Unable to generate roc.cmd.r file', code="34101501")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/roc.cmd.r' % (self.work_dir + '/ROC'), shell=True)
            self.logger.info('ROC计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('ROC计算失败')
            self.set_error('R running calculation ROC failed', code="34101502")
        self.logger.info('运行calc_roc.pl程序进行ROC计算完成')

    def set_output(self):
        files =os.listdir(self.work_dir + '/ROC')
        for file in files:
            if re.search(r'.xls',file):
                if os.path.exists(self.output_dir + '/' +file):
                    os.remove(self.output_dir + '/' +file)
                os.link(self.work_dir + '/ROC/' + file,self.output_dir + '/' +file)

