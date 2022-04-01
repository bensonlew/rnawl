# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import glob
import os
import unittest

class MirnaEditModule(Module):
    def __init__(self, work_id):
        super(MirnaEditModule, self).__init__(work_id)
        options = [
            {'name': 'list_file', 'type': 'infile', 'format': 'datasplit.list_file'},
            {'name': 'species', 'type': 'string', 'default': ''},
            {'name': 'index', 'type': 'string', 'default': ''},
            {'name': 'hairpin_fa', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'mature_fa', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'result_dir', 'type': 'outfile', 'format': 'small_rna.common_dir'},
        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))

        self.logger.debug('{} - {}'.format('list_file', self.option('list_file')))
        if not self.option('list_file').is_set:
            raise OptionError('list_file must be provided')

        self.logger.debug('{} - {}'.format('species', self.option('species')))
        if self.option('species') == '':
            raise OptionError('specie abbreviation must be provided')

        self.logger.debug('{} - {}'.format('index', self.option('index')))
        if self.option('index') == '':
            raise OptionError('bowtie index location must be provided')

        self.logger.debug('{} - {}'.format('hairpin_fa', self.option('hairpin_fa').prop['path']))
        if not self.option('hairpin_fa').is_set:
            raise OptionError('pre-miRNA fasta from mirbase must be provided')

        self.logger.debug('{} - {}'.format('mature_fa', self.option('mature_fa').prop['path']))
        if not self.option('mature_fa').is_set:
            raise OptionError('mature-miRNA fasta from mirbase must be provided')

        self.logger.debug('{} - {}'.format('result_dir', self.option('result_dir')))

        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        super(MirnaEditModule, self).run()
        self.fastqs = self.option('list_file').prop['samples']
        self.logger.debug('{}'.format(type(self.fastqs)))
        self.logger.debug('{}'.format(self.fastqs))
        self.run_tools()

    def run_tools(self):
        options = {
            'species': self.option('species'),
            'index': self.option('index'),
            'hairpin_fa': self.option('hairpin_fa').prop['path'],
            'mature_fa': self.option('mature_fa').prop['path'],
        }
        tools = list()
        for key in self.fastqs.keys():
            options['sample'] = key
            options['filtered_fq'] = self.fastqs[key][0]
            tool = self.add_tool('small_rna.mirna_edit')
            tool.set_options(options)
            self.tools.append(tool)
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def set_output(self):
        '''
        link result in output of tool to output_dir of module
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        with open(os.path.join(self.output_dir, 'result.info.txt'), 'w') as f:
            for tool in self.tools:
                for source in glob.glob(os.path.join(tool.output_dir, '*')):
                    link_name = os.path.join(self.output_dir, os.path.basename(source))
                    if os.path.exists(link_name):
                        os.remove(link_name)
                    os.link(source, link_name)
                    self.logger.info('succeed in linking {} to {}'.format(source, link_name))
                    f.write('{}\t{}\n'.format(link_name, tool.option('sample')))
        self.option('result_dir', self.output_dir)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', 'mirna_edit_module_output_dir']
        ])
        super(MirnaEditModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''
    def test_hsa(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mirna_edit_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'small_rna.mirna_edit',
            'instant': False,
            'options': {
                'list_file': '/mnt/ilustre/users/sanger-dev/workspace/20181214/Smallrna_tsg_33050/MirnaQc/output/clean_data/list.txt',
                'species': 'hsa',
                'index': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean_index',
                'hairpin_fa': '/mnt/ilustre/users/sanger-dev/app/database/mirbase/hairpin.fa',
                'mature_fa': '/mnt/ilustre/users/sanger-dev/app/database/mirbase/mature.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()