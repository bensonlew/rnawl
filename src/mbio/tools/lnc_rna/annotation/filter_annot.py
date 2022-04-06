# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
from biocluster.tool import Tool
import unittest

class FilterAnnotAgent(Agent):
    '''
    last_modify: 2019.02.12
    '''
    def __init__(self, parent):
        super(FilterAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'types', 'type': 'string', 'default': ''},
            {'name': 'xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'hmm', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'identity', 'type': 'float', 'default': 0},
            {'name': 'similarity', 'type': 'float', 'default': 0},
            {'name': 'keep_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'outxml', 'type': 'outfile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'outtable', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('types', self.option('types')))
        if self.option('types') == 'xml':
            if self.option('xml').is_set:
                self.logger.debug('{} - {}'.format('xml', self.option('xml').prop['path']))
                self.infile_size = os.path.getsize(self.option('xml').prop['path'])
            else:
                raise OptionError('BLAST XML must be provided')
        elif self.option('types') == 'hmm':
            if self.option('hmm').is_set:
                self.logger.debug('{} - {}'.format('hmm', self.option('hmm').prop['path']))
                self.infile_size = os.path.getsize(self.option('hmm').prop['path'])
            else:
                raise OptionError('PFAM DOMAIN must be provided')
        elif self.option('types') == 'go':
            if self.option('blast2go_annot').is_set:
                self.logger.debug('{} - {}'.format('blast2go_annot', self.option('blast2go_annot').prop['path']))
                self.infile_size = os.path.getsize(self.option('blast2go_annot').prop['path'])
            else:
                raise OptionError('BLAST2GO XLS must be provided')
        else:
            raise OptionError('option types must be xml or hmm or go')
        self.logger.debug('{} - {}'.format('evalue', self.option('evalue')))
        if not 0 <= self.option('evalue') <= 0.001:
            raise OptionError('e-value must be set in [0-0.001]')
        if self.option('types') != 'hmm':
            self.logger.debug('{} - {}'.format('identity', self.option('identity')))
            if not 0 <= self.option('identity') <= 100:
                raise OptionError('identity must be set in [0-100]')
            self.logger.debug('{} - {}'.format('similarity', self.option('similarity')))
            if not 0 <= self.option('identity') <= 100:
                raise OptionError('similarity must be set in [0-100]')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 8 + 16))

    def end(self):
        super(FilterAnnotAgent, self).end()


class FilterAnnotTool(Tool):
    def __init__(self, config):
        super(FilterAnnotTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.filter_blast_xml_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_blast_xml.py')
        self.filter_pfam_domain_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_pfam_domain.py')
        self.filter_blast2go_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_blast2go.py')
        self.blast_filter_xml = os.path.join(self.work_dir, 'blast.filter.xml')
        self.pfam_filter_xls = os.path.join(self.work_dir, 'pfam.filter.xls')
        self.blast2go_filter_xls = os.path.join(self.work_dir, 'blast2go.filter.xls')

    def run(self):
        super(FilterAnnotTool, self).run()
        if self.option('types') == 'xml':
            self.run_filter_xml()
        elif self.option('types') == 'hmm':
            self.run_filter_hmm()
        elif self.option('types') == 'go':
            self.run_filter_go()
        self.set_output()
        self.end()

    def run_filter_xml(self):
        cmd = '{} {}'.format(self.python, self.filter_blast_xml_py)
        cmd += ' -i {}'.format(self.option('xml').prop['path'])
        cmd += ' -e {}'.format(self.option('evalue'))
        cmd += ' -c {}'.format(self.option('identity'))
        cmd += ' -s {}'.format(self.option('similarity'))
        cmd += ' -o {}'.format(self.blast_filter_xml)
        if self.option('keep_list').is_set:
            cmd += ' -k {}'.format(self.option('keep_list').prop['path'])
        cmd_name = 'run_filter_xml'
        self.run_code(cmd_name, cmd)

    def run_filter_hmm(self):
        cmd = '{} {}'.format(self.python, self.filter_pfam_domain_py)
        cmd += ' -i {}'.format(self.option('hmm').prop['path'])
        cmd += ' -e {}'.format(self.option('evalue'))
        cmd += ' -o {}'.format(self.pfam_filter_xls)
        if self.option('keep_list').is_set:
            cmd += ' -k {}'.format(self.option('keep_list').prop['path'])
        cmd_name = 'run_filter_hmm'
        self.run_code(cmd_name, cmd)

    def run_filter_go(self):
        cmd = '{} {}'.format(self.python, self.filter_blast2go_py)
        cmd += ' -i {}'.format(self.option('blast2go_annot').prop['path'])
        cmd += ' -e {}'.format(self.option('evalue'))
        cmd += ' -c {}'.format(self.option('identity'))
        cmd += ' -s {}'.format(self.option('similarity'))
        cmd += ' -o {}'.format(self.blast2go_filter_xls)
        if self.option('keep_list').is_set:
            cmd += ' -k {}'.format(self.option('keep_list').prop['path'])
        cmd_name = 'run_filter_go'
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
        if self.option('types') == 'xml':
            source = self.blast_filter_xml
            link_name = os.path.join(self.output_dir, 'blast.filter.xml')
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.option('outxml', link_name)
        elif self.option('types') == 'hmm':
            source = self.pfam_filter_xls
            link_name = os.path.join(self.output_dir, 'pfam.filter.xls')
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.option('outtable', link_name)
        elif self.option('types') == 'go':
            source = self.blast2go_filter_xls
            link_name = os.path.join(self.output_dir, 'blast2go.filter.xls')
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.option('outtable', link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_nr(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'filter_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.filter_annot',
            'instant': False,
            'options': {
                'types': 'xml',
                'xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/nr/blast.xml',
                'evalue': '1e-6',
                'identity': '90',
                'similarity': '90',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_swissprot(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'filter_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.filter_annot',
            'instant': False,
            'options': {
                'types': 'xml',
                'xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/swissprot/blast.xml',
                'evalue': '1e-6',
                'identity': '90',
                'similarity': '90',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_eggnog(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'filter_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.filter_annot',
            'instant': False,
            'options': {
                'types': 'xml',
                'xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/eggnog/blast.xml',
                'evalue': '1e-6',
                'identity': '90',
                'similarity': '90',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_kegg(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'filter_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.filter_annot',
            'instant': False,
            'options': {
                'types': 'xml',
                'xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/kegg/blast.xml',
                'evalue': '1e-6',
                'identity': '90',
                'similarity': '90',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_pfam(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'filter_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.filter_annot',
            'instant': False,
            'options': {
                'types': 'hmm',
                'hmm': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_orfpfam/pfam_domain',
                'evalue': '1e-6',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_go(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'filter_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.filter_annot',
            'instant': False,
            'options': {
                'types': 'go',
                'blast2go_annot': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/GO/blast2go_merge.xls',
                'evalue': '1e-6',
                'identity': '90',
                'similarity': '90',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
