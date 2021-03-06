# -*- coding: utf-8 -*-
# __author__ = 'zengjing, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class EggnogAnnotAgent(Agent):
    '''
    last_modify: 2019.02.14
    '''
    def __init__(self, parent):
        super(EggnogAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'blast_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'blast_table', 'type': 'outfile', 'format': 'lnc_rna.blast_table'},
            {'name': 'cog_table', 'type': 'outfile', 'format': 'lnc_rna.cog_table'},
            {"name": "eggnog_version", "type": "string", "default": ""},

        ]
        self.add_option(options)
        self.step.add_steps('eggnog_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.eggnog_annot.start()
        self.step.update()

    def step_end(self):
        self.step.eggnog_annot.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('blast_xml').is_set:
            self.logger.debug('{} = {}'.format('blast_xml', self.option('blast_xml').prop['path']))
            self.infile_size = os.path.getsize(self.option('blast_xml').prop['path'])
        else:
            raise OptionError('input BLAST XML must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 4 + 4))

    def end(self):
        super(EggnogAnnotAgent, self).end()

class EggnogAnnotTool(Tool):
    def __init__(self, config):
        super(EggnogAnnotTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.xml2table_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/xml2table.py')
        self.cog_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/cog_annotation.py')
        self.xml = self.option('blast_xml').prop['path']
        if self.xml[-4:] == '.xml':
            self.xls = os.path.join(self.work_dir, '{}.xls'.format(os.path.basename(self.xml)[:-4]))
        else:
            self.xls = os.path.join(self.work_dir, '{}.xls'.format(os.path.basename(self.xml)))
        self.cog_xls = os.path.join(self.work_dir, 'cog.xls')

    def run(self):
        super(EggnogAnnotTool, self).run()
        self.run_xml2table()
        self.run_cog_annotation()
        self.set_output()
        self.end()

    def run_xml2table(self):
        cmd = '{} {} -a'.format(self.python, self.xml2table_py)
        cmd += ' -i {}'.format(self.xml)
        cmd += ' -o {}'.format(self.xls)
        cmd_name = 'run_xml2table'
        self.run_code(cmd_name, cmd)

    def run_cog_annotation(self):
        cmd = '{} {}'.format(self.python, self.cog_annotation_py)
        cmd += ' -i {}'.format(self.xls)
        cmd += ' -o {}'.format(self.cog_xls)
        cmd_name = 'run_cog_annotation'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
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
        blast_table = os.path.join(self.output_dir, os.path.basename(self.xls))
        if os.path.exists(blast_table):
            os.remove(blast_table)
        os.link(self.xls, blast_table)
        self.option('blast_table', blast_table)
        self.logger.info('succeed in linking {} to {}'.format(self.xls, blast_table))
        cog_table = os.path.join(self.output_dir, os.path.basename(self.cog_xls))
        if os.path.exists(cog_table):
            os.remove(cog_table)
        os.link(self.cog_xls, cog_table)
        self.option('cog_table', cog_table)
        self.logger.info('succeed in linking {} to {}'.format(self.cog_xls, cog_table))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'eggnog_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.eggnog_annot',
            'instant': False,
            'options': {
                'blast_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/eggnog/blast.filter.xml',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
