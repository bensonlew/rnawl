# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.otu_table import OtuTableFile
import os
import unittest
import re


class DistanceCalcAgent(Agent):
    """
    PCoA
    """
    def __init__(self, parent):
        super(DistanceCalcAgent, self).__init__(parent)
        options = [
            {"name": "method", "type": "string", "default": "euclidean"},
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "dis_matrix", "type": "outfile", "format": "meta.beta_diversity.distance_matrix"},
        ]
        self.add_option(options)
        self.step.add_steps('distance_calc')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.distance_calc.start()
        self.step.update()

    def step_end(self):
        self.step.distance_calc.finish()
        self.step.update()

    def check_options(self):
        if not self.option('otutable').is_set:
            raise OptionError('必须提供输入文件', code="32701901")
        if self.option('method') not in ['euclidean', 'bray_curtis', 'manhattan']:
            raise OptionError('错误或者不支持的距离矩阵计算方法', code="32701903")
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 5
        size = os.path.getsize(self.option('otutable').prop['path'])
        num = float(size) / float(1000000)
        if num <= 10:
            num2 = 20
        else:
            num2 = int(num) + 10
        self._memory = str(num2) + 'G'

    def end(self):
        super(DistanceCalcAgent, self).end()


class DistanceCalcTool(Tool):
    def __init__(self, config):
        super(DistanceCalcTool, self).__init__(config)
        self._version = "1.0"
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.program = {
            'python': "program/Python/bin/python",
        }
        self.script = {
            'beta_div':  os.path.join(self.config.SOFTWARE_DIR, "program/Python/bin/beta_diversity.py")
        }

    def run_beta_diversity(self):
        """
        运行qiime:beta_diversity.py
        """
        biom = self.biom_otu_table()
        cmd = self.program['python'] + ' ' + self.script['beta_div']
        cmd += ' -m %s -i %s -o %s' % (self.option('method'), biom,
                                       self.work_dir)
        self.logger.info('运行qiime:beta_diversity.py程序')
        self.logger.info(cmd)
        dist_matrix_command = self.add_command('distance_matrix', cmd)
        dist_matrix_command.run()
        self.wait()
        if dist_matrix_command.return_code == 0:
            self.command_successful()
        elif dist_matrix_command.return_code is None:
            self.logger.warn("运行命令出错，返回值为None，尝试重新运行")
            dist_matrix_command.rerun()
            self.wait()
            if dist_matrix_command.return_code is 0:
                self.command_successful()
            else:
                self.set_error("运行qiime:beta_diversity.py出错", code="32701901")

        else:
            self.set_error('运行qiime:beta_diversity.py出错', code="32701902")

    def command_successful(self):
        self.logger.info('运行qiime:beta_diversity.py完成')
        filename = self.work_dir + '/' + \
                   self.option('method') + '_temp.txt'
        basename = os.path.splitext(os.path.basename(self.option('otutable').prop['path']))[0]
        linkfile = self.output_dir + '/' + \
                   self.option('method') + '_' + basename + '.xls'
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename, linkfile)
        self.dis_matrix_check(linkfile)
        self.option('dis_matrix').set_path(linkfile)
        self.end()

    def dis_matrix_check(self, linkfile):  # 20171031 by zengjing 增加距离矩阵结果检查，便于将具体报错带到页面
        """聚类矩阵结果检查"""
        self.logger.info("开始进行距离矩阵检查")
        dist_dict = dict()
        all_values = []
        with open(linkfile, 'r') as f:
            head = f.readline().rstrip().split('\t')
            head_len = len(head)
            head = head[1:]
            for line in f:
                all_nums = line.rstrip().split('\t')
                if len(all_nums) != head_len:
                    self.set_error('距离矩阵每行数据量格式不正确', code="32701903")
                values = dict(zip(head, all_nums[1:]))
                all_values.extend(all_nums[1:])
                dist_dict[all_nums[0]] = values
            for samp1 in head:
                for samp2 in head:
                    if dist_dict[samp1][samp2] != dist_dict[samp2][samp1]:
                        self.set_error('距离矩阵数据不对称', code="32701904")
            all_values = [float(i) for i in all_values]
            all_plus = sum(all_values)
            if all_plus == 0 and len(all_values) > 1:  # 只有一个样本时距离就为零，不做处理，但是后续不能做任何分析
                self.set_error('所有距离矩阵值全部为零', code="32701905")
            if len(head) != len(set(head)):
                self.set_error('距离矩阵存在重复的样本名', code="32701906")
        self.logger.info("距离矩阵检查正确")

    def biom_otu_table(self):
        """
        将otutable转化成biom格式
        :return biom_path:返回生成的biom文件路径
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            # 如果提供的是otu分类文件夹，需要重新创建类，再使用类的convert_to_biom方法
            newtable = OtuTableFile()
            newtable.set_path(self.option('otutable').prop['path'])
            newtable.check()
        else:
            newtable = self.option('otutable')
        newtable.get_info()
        biom_path = os.path.join(self.work_dir, 'temp.biom')
        if os.path.isfile(biom_path):
            os.remove(biom_path)
        newtable.convert_to_biom(biom_path)
        return biom_path

    def run(self):
        super(DistanceCalcTool, self).run()
        self.run_beta_diversity()

