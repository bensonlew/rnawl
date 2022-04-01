# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os

from biocluster.agent import Agent
from biocluster.tool import Tool


class GcInfosAgent(Agent):
    def __init__(self, parent):
        super(GcInfosAgent, self).__init__(parent)
        options = [
            {'name': 'ref', 'type': 'string'},
            {'name': 'seq_dir', 'type': 'string'},
            {'name': 'window', 'type': 'int', 'default': 10000},
            {'name': 'step', 'type': 'int', 'default': 10000}
        ]
        self.add_option(options)
        self.step.add_steps('gc_infos')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gc_infos.start()
        self.step.update()

    def stepfinish(self):
        self.step.gc_infos.finish()
        self.step.update()

    def check_options(self):
        pass

    def set_resource(self):
        self._memory = '4G'
        self._cpu = 1

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', "", "结果输出目录"],
            ])
        super(GcInfosAgent, self).end()


class GcInfosTool(Tool):
    def __init__(self, config):
        super(GcInfosTool, self).__init__(config)
        self.python = '/program/Python/bin/python'
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/gc_infos.py'

    def run(self):
        super(GcInfosTool, self).run()
        self.gc_run()
        self.set_output()
        self.end()

    def gc_run(self):
        self.get_genome()
        cmd = self.python + ' ' + self.package_path +\
            ' {} {} {}'.format(
                self.option('ref') + '.fna',
                self.option('window'), self.option('step'),
            )
        command = self.add_command('gc_run', cmd).run()
        self.wait(command)

        if command.return_code == 0:
            self.logger.info('package gc_infos.py运行成功')
        else:
            self.set_error('package gc_infos.py运行出错')

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))

    def end(self):

        super(GcInfosTool, self).end()

    def get_genome(self):
        fasta_dir = os.path.join(self.option('seq_dir'), self.option('ref'))
        cmd = '/bin/cat {}/* > {}.fna'.format(fasta_dir, self.option('ref'))
        self.add_command('gc_infos__get_genome', cmd, shell=True).run()
        self.wait()
