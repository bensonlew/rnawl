# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class StatisticsAgent(Agent):
    '''
    last_modify: 2019.04.10
    '''
    def __init__(self, parent):
        super(StatisticsAgent, self).__init__(parent)
        options = [
            {'name': 'exp_type', 'type': 'string', 'default': None},  # T G
            {'name': 'method', 'type': 'string', 'default': None},  # rsem salmon kallisto
            {'name': 'loc2name', 'type': 'infile', 'format': 'ref_rna_v2.loc2name'},
            {'name': 't2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps('statistics')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.statistics.start()
        self.step.update()

    def stepfinish(self):
        self.step.statistics.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        super(StatisticsAgent, self).end()

class StatisticsTool(Tool):
    def __init__(self, config):
        super(StatisticsTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.statistics_py = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/expression/statistics.py')

    def run(self):
        super(StatisticsTool, self).run()
        self.run_statistics()
        self.set_output()
        self.end()

    def run_statistics(self):
        cmd = '{} {}'.format(self.python, self.statistics_py)
        cmd += ' -i {}'.format(self.option('loc2name').path)
        cmd += ' -m {}'.format(self.option('method'))
        cmd += ' -t {}'.format(self.option('exp_type'))
        cmd += ' -c {}'.format(self.option('t2g').path)
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_statistics'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_kallisto_g(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'statistics_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.statistics',
            'instant': False,
            'options': {
                'exp_type': 'G',
                'method': 'kallisto',
                'loc2name': '/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_expression_6347_1439/Expression/location2name.kallisto.G.tsv',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_expression_6347_1439/Expression/Prepare/output/t2g.pairs'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_salmon_t(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'statistics_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.statistics',
            'instant': False,
            'options': {
                'exp_type': 'T',
                'method': 'salmon',
                'loc2name': '/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_expression_1264_1613/Expression/location2name.salmon.T.tsv',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_expression_1264_1613/Expression/Prepare/output/t2g.pairs'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_kallisto_g')])
    unittest.TextTestRunner(verbosity=2).run(suite)
