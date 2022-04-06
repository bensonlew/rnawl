# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import re
import numpy as np

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from collections import Counter


class DiffStatisticalAgent(Agent):
    def __init__(self, parent):
        super(DiffStatisticalAgent, self).__init__(parent)
        options = [
            {'name': 'test', 'type': 'string'},  # 检验类型
            {'name': 'testtype', 'type': 'string', 'default': 'two.side'},
            {'name': 'multitest', 'type': 'string', 'default': 'fdr'},
            {'name': 'ci', 'type': 'float', 'default': 0.95},
            {'name': 'intable', 'type': 'string', 'default': ''},  # 输入table文件, 有行、列名
            {'name': 'samp1', 'type': 'string', 'default': ''},  # 两样本检验所选样本
            {'name': 'samp2', 'type': 'string', 'default': ''},  # 两样本检验所选样本
            {'name': 'gfile', 'type': 'string', 'default': ''},  # 分组文件
            {'name': 'coverage', 'type': 'float', 'default': 0.95},
            {'name': 'meth', 'type': 'string', 'default': ''},
            {'name': 'norm', 'type': 'string', 'default': 'T'},
        ]
        self.add_option(options)
        self.step.add_steps('diff_stats')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.diff_stats.start()
        self.step.update()

    def stepfinish(self):
        self.step.diff_stats.finish()
        self.step.update()

    def check_options(self):
        if not self.option('intable'):
            raise OptionError('比较给出 输入 文件')
        if not self.option('test'):
            raise OptionError('请提供检验类型')
        else:
            test = self.option('test')
            samples = []
            if test in ['chi', 'fisher']:
                if not self.option('samp1') or not self.option('samp2'):
                    raise OptionError('两样本比较给出要比较的两样本的名字')
                self.check_samps()
                samples = [self.option('samp1'), self.option('samp2')]
            elif test in ['kru_H', 'anova', 'student', 'welch', 'mann']:
                if not self.option('gfile'):
                    raise OptionError('两组或多组比较必须给出样本分组文件')
                self.check_gfile()
            else:
                raise OptionError('不存在此检验名称%s' % test)
        return True

    def check_samps(self):
        with open(self.option('intable'), 'r') as f:
            l = f.readline().strip().split('\t')[1:]
            if self.option('samp1') not in l:
                raise OptionError('结果中样本没有%s' % self.option('samp1'))
            if self.option('samp2') not in l:
                raise OptionError('结果中样本没有%s' % self.option('samp2'))
            if self.option('samp1') == self.option('samp2'):
                raise OptionError('需要比较选择两个不同的样本')

    def check_gfile(self):
        with open(self.option('intable'), 'r') as f1,\
                open(self.option('gfile'), 'r') as f2:
            allsamps = f1.readline().strip().split('\t')[1:]
            gsamps = []
            for l in f2.readlines():
                if l.startswith('#'):
                    continue
                gsamps.append(l.split('\t'))
            gsamps = np.array(gsamps)
            if len(gsamps[:, 0]) != len(set(gsamps[:, 0])):
                raise OptionError('分组中不能有重复样本，请检查并重新分组！')
            if not set(gsamps[:, 0]) <= set(allsamps):
                raise OptionError('分组中存在不属于结果文件中的样本名')
            c = Counter(gsamps[:, 1])
            for m, n in c.items():
                if n < 5:
                    raise OptionError('每个分组中至少5个样本，请重新分组！')

    def set_resource(self):
        self._cpu = 2
        self._memory = '4G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     ['.', '', '结果输出目录']
        # ])
        # result_dir.add_regexp_rules([
        #     []
        # ])
        super(DiffStatisticalAgent, self).end()


class DiffStatisticalTool(Tool):
    def __init__(self, config):
        super(DiffStatisticalTool, self).__init__(config)
        self.python = '/miniconda2/bin/python'
        self.package_path = self.config.PACKAGE_DIR +\
            '/bac_comp_genome/diff_statistical.py'

    def run(self):
        super(DiffStatisticalTool, self).run()
        self.run_stat()
        self.set_output()
        self.end()

    def run_stat(self):
        cmd = self.python + ' ' + self.package_path +\
            ' -t {} -tp {} -mt {} -ci {} -i {} -c {} -n {}'.format(
                self.option('test'), self.option('testtype'),
                self.option('multitest'), self.option('ci'),
                self.option('intable'), self.option('coverage'),
                self.option('norm')
            )
        if self.option('meth'):
            cmd += ' -m {}'.format(self.option('meth'))
        if self.option('samp1'):
            cmd += ' -s1 {} -s2 {}'.format(
                self.option('samp1'), self.option('samp2')
            )
        else:
            cmd += ' -g {}'.format(self.option('gfile'))
        command = self.add_command('diff_stats', cmd, ignore_error=True).run()
        self.wait(command)

        if command.return_code == 0:
            pass
        else:
            self.set_error('package diff_statistical.py运行出错')

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        test = self.option('test')
        for r, d, f in os.walk(self.work_dir, topdown=False):
            files = f
        try:
            for f in files:
                if re.match(r'.*_result.xls|.*CI.*xls', f):
                    o = '/' + f
                    os.link(self.work_dir + o, self.output_dir + o)
            self.logger.info('设置{}检验结果目录成功！'.format(test))
        except:
            self.logger.info('设置{}检验结果目录失败！'.format(test))
