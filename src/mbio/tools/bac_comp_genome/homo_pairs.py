# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class HomoPairsAgent(Agent):
    def __init__(self, parent):
        super(HomoPairsAgent, self).__init__(parent)
        options = [
            {'name': 'gene_dir', 'type': 'string'},
            {'name': 'samples', 'type': 'string'},
        ]
        self.add_option(options)
        self.step.add_steps('homo_pairs')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.homo_pairs.start()
        self.step.update()

    def stepfinish(self):
        self.step.homo_pairs.finish()
        self.step.update()

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = len(self.option('samples').split(',')) + 1
        self._memory = str(self._cpu) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', '']
        ])
        super(HomoPairsAgent, self).end()


class HomoPairsTool(Tool):
    def __init__(self, config):
        super(HomoPairsTool, self).__init__(config)
        self.python = '/miniconda2/bin/python'
        self.diamond = self.config.SOFTWARE_DIR + '/bioinfo/align/diamond-0.9.11/diamond'
        self.package = self.config.PACKAGE_DIR + '/bac_comp_genome/homo_pairs.py'

    def run(self):
        super(HomoPairsTool, self).run()
        self.homo_pairs_run()
        self.set_output()
        self.end()

    def homo_pairs_run(self):
        samples = self.get_faa()
        cmd = self.python + ' ' + self.package + ' -d {} -s {}'.format(
            self.diamond, samples
        )
        command = self.add_command('homo_pairs', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('package ' + self.package + ' 运行成功')
        else:
            self.set_error(self.package + ' 运行出错')

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))
        pass

    def end(self):
        super(HomoPairsTool, self).end()

    def get_faa(self):
        sp_list = self.option('samples').split(',')
        samples = []
        for f in sp_list:
            faa_dir = os.path.join(self.option('gene_dir'), f)
            samples.append(faa_dir + '/' + f + '_CDS.faa')
        return ','.join(samples)
