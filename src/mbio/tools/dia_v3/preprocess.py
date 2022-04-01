# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import pandas as pd
from mbio.packages.dia_v3.preprocess import preprocess
import subprocess
import glob
import unittest
from collections import OrderedDict


class PreprocessAgent(Agent):
    def __init__(self, parent):
        super(PreprocessAgent, self).__init__(parent)
        options = [
            {"name": "raw_path", "type": "infile", "format": "labelfree.common"},   # exp file path
            {"name": "group_table", "type": "infile", "format": "labelfree.group_table"},   # group file path
            {"name": "all_eliminate", "type": "string", "default": "all"},  #
            {"name": "all_percent", "type": "float", "default": 90},
            {"name": "if_group", "type": "string", "default": "yes"},   # if perform specific
            {"name": "group_specific", "type": "string", "default": "any"},
            {"name": "group_percent", "type": "float", "default": 50},
            {"name": "fillna", "type": "string", "default": "seqknn"},  # method for fill na
            {"name": "fill_type", "type": "string", "default": "group"},    # based on all/group
            {"name": "interactive", "type": "string", "default": "no"},
            {"name": "preprocess_exp", "type": "outfile", "format": "labelfree.common"},    # post-preprocess exp file
        ]
        self.add_option(options)
        self.step.add_steps("preprocess")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.preprocess.start()
        self.step.update()

    def step_end(self):
        self.step.preprocess.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        if not self.option("raw_path"):
            raise OptionError("必须设置输入文件:蛋白表达量文件")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(PreprocessAgent, self).end()


class PreprocessTool(Tool):
    def __init__(self, config):
        super(PreprocessTool, self).__init__(config)
        # self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        # self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = self.config.SOFTWARE_DIR + '/bioinfo/dia_v3/miniconda2/bin/Rscript'

    def run_preprocess(self):
        preprocess(exp_file=self.option("raw_path").prop['path'], group_file=self.option("group_table").prop['path'],
                   na_ratio=float(self.option("all_percent")/100), method=self.option("fillna"),
                   specific=self.option("if_group"), group_threshold=float(self.option("group_percent")/100),
                   fill_type=self.option("fill_type"), out_prefix=self.option("fillna"),
                   interactive=self.option('interactive'))
        self.cmd = self.r_path + " preprocess_%s" % self.option('fillna') + ".r" #这里有一个空格，作为运行命令
        return self.cmd

    def run_cmd(self):
        """
        运行命令行
        """
        self.logger.info("开始运行命令！")
        cmd = self.run_preprocess()
        # self.logger.info(cmd)
        # try:
        #     subprocess.check_call(cmd, shell=True)
        #     self.logger.info(cmd + " 运行完成！")
        # except subprocess.CalledProcessError:
        #     import traceback
        #     print(traceback.format_exc())
        #     self.logger.info('CMD:{}'.format(cmd))

        cmd_name = 'run_preprocess'
        command = self.add_command(cmd_name, cmd, shell=True, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))
        else:
            out_file = glob.glob(os.path.join(self.work_dir, '*.o'))[0]
            with open(out_file, 'r') as o:
                o.readline()
                f = o.next().strip().split('"')[1]
            print(f)
            if f == 'OOPS':
                self.set_error('请检查参数，筛选后的数据集为空。')
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))

    def cv_stats(self):
        cv_paths = glob.glob(os.path.join(self.work_dir, '*cv_stats*'))
        for cv_matrix in cv_paths:
            stat_dict_list = list()
            box_file = cv_matrix.replace('{}_cv_stats'.format(self.option('fillna')), 'boxdata')
            group_cv_pd = pd.read_table(cv_matrix, index_col=0, header=0)
            group_cv_pd.index.name = 'accession_id'
            group_cv_pd = group_cv_pd[group_cv_pd.sum(axis=1) > 0.001]
            target_columns = group_cv_pd.columns
            for each in target_columns:
                cv_pd = group_cv_pd[each]
                cv_pd = cv_pd[cv_pd != 0]
                summary = cv_pd.describe()
                summary.index = [u'count', u'mean', u'std', u'min', u'q1', u'median', u'q3', u'max']
                tmp_dict = summary.to_dict(OrderedDict)
                lt25 = cv_pd[cv_pd <= tmp_dict['q1']].shape[0]
                lt50 = cv_pd[cv_pd <= tmp_dict['median']].shape[0]
                lt75 = cv_pd[cv_pd <= tmp_dict['q3']].shape[0]
                upper_whisker = tmp_dict['q3'] + 1.5*(tmp_dict['q3'] - tmp_dict['q1'])
                lower_whisker = tmp_dict['q1'] - 1.5*(tmp_dict['q3'] - tmp_dict['q1'])
                if lower_whisker < 0:
                    lower_whisker = 0
                upper_outliers = map(str, list(cv_pd[cv_pd > upper_whisker]))
                lower_outliers = map(str, list(cv_pd[cv_pd < lower_whisker]))
                tmp_list = [each.split('cv_')[1]] + tmp_dict.values() + [lt25, lt50-lt25, lt75-lt50,
                                                                         cv_pd.shape[0]-lt75, upper_whisker,
                                                                         lower_whisker, ';'.join(upper_outliers),
                                                                         ';'.join(lower_outliers)]
                stat_dict_list.append(tmp_list)
            with open(box_file, 'w') as o:
                o.write('group\tcount\tmean\tstd\tmin\tq1\tmedian\tq3\tmax\tmin-q1\tq1-median\tmedian-q3\tq3-max\tupper_whisker\tlower_whisker\tupper_outliers\tlower_outliers\n')
                for each in stat_dict_list:
                    o.write('\t'.join(map(str, each)) + '\n')

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        for f in os.listdir(self.work_dir):
            if f.startswith(self.option("fillna")):
                os.link(self.work_dir + "/" + f, os.path.join(self.output_dir, f))
        self.logger.info("设置搜库结果目录")
        fillna_path = os.path.join(self.work_dir, "{}_fillna.txt".format(self.option("fillna")))
        self.option("preprocess_exp", fillna_path)

    def run(self):
        super(PreprocessTool, self).run()
        self.run_cmd()
        self.cv_stats()
        self.set_output()
        self.end()

