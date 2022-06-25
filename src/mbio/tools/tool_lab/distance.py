# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os, re
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.otu_table import OtuTableFile
import unittest
import pandas as pd
from mbio.packages.whole_transcriptome.utils import runcmd

class DistanceAgent(Agent):
    """
    qiime
    version 1.0
    author shenghe
    last_modified:2018.02.02 增加环境因子距离计算功能 add by zhujuan
    """
    METHOD = ['abund_jaccard', 'binary_chisq', 'binary_chord',
              'binary_euclidean', 'binary_hamming', 'binary_jaccard',
              'binary_lennon', 'binary_ochiai',
              'binary_pearson', 'binary_sorensen_dice',
              'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
              'canberra', 'chisq', 'chord', 'euclidean', 'gower',
              'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
              'pearson', 'soergel', 'spearman_approx', 'specprof',]


    def __init__(self, parent):
        super(DistanceAgent, self).__init__(parent)
        options = [
            {"name": "method", "type": "string", "default": "bray_curtis"},
            {"name": "otutable", "type": "infile","format": "ref_rna_v2.common"},
            {"name": "dis_matrix", "type": "outfile","format": "ref_rna_v2.common"},
            {'name': 'sep', 'type':'string'}

        ]
        self.add_option(options)
        self._memory_increase_step = 30     # added by zhangyitong on 20210805
        self.step.add_steps('distance')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.distance.start()
        self.step.update()

    def step_end(self):
        self.step.distance.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        pass
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "距离矩阵计算结果输出目录"],
        # ])
        # result_dir.add_regexp_rules([
        #     [r'%s.*\.xls' % self.option('method'), 'xls', '样本距离矩阵文件']
        # ])
        # # print self.get_upload_files()
        super(DistanceAgent, self).end()


class DistanceTool(Tool):
    def __init__(self, config):
        super(DistanceTool, self).__init__(config)
        # self._version = '1.9.1'  # qiime版本
        # self.cmd_path = 'miniconda2/bin/beta_diversity.py'
        # # 设置运行环境变量
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        # self.real_otu = self.gettable()  # 获取真实的丰度表路劲
        # self.biom = self.biom_otu_table()  # 传入丰度表需要转化为biom格式
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'distance': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/distance_kit.py')
        }
        self.file = {
            'dis_matrix': os.path.join(self.output_dir, '{}'.format(os.path.basename(self.option("otutable").prop['path']))),
        }

    def run(self):
        """
        运行
        """
        super(DistanceTool, self).run()
        self.run_distance()
        self.set_output()
        self.end()

    def run_distance(self):
        cmd = '{} {} '.format(self.program['python'], self.script['distance'])
        cmd += '--input {} '.format(self.option('otutable').path)
        cmd += '--sep {} '.format(self.option('sep'))
        cmd += '--metrics {} '.format(self.option('method').lower())
        cmd += '--output {} '.format(self.file['dis_matrix'])
        cmd_name = 'run_distance'
        # runcmd(self,cmd_name,cmd)
        self.logger.info("start calculating distance matrix")
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("succeed in calculating distance matrix")
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in calculating distance matrix")
        else:
            self.set_error("fail to calculate distance matrix")

    def set_output(self):
        self.option('dis_matrix').set_path(self.file['dis_matrix'])

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "distance_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.distance_calc",
            "options": dict(
                method="bray_curtis",
                otutable="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/cluster/otu_table.xls",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()