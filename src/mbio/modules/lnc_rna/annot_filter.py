# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import os
import shutil
import unittest

class AnnotFilterModule(Module):
    '''
    last_modify: 2019.03.01
    '''
    def __init__(self, work_id):
        super(AnnotFilterModule, self).__init__(work_id)
        options = [
            {'name': 'db', 'type': 'string', 'default': 'nr,swissprot,kegg,eggnog,pfam,go'},
            # input files
            {'name': 'blast_nr_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'blast_swissprot_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'blast_eggnog_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'blast_kegg_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'lnc_rna.common'},
            # input parameters
            {'name': 'nr_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'nr_similarity', 'type': 'float', 'default': 0},
            {'name': 'nr_identity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'swissprot_similarity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_identity', 'type': 'float', 'default': 0},
            {'name': 'eggnog_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'eggnog_similarity', 'type': 'float', 'default': 0},
            {'name': 'eggnog_identity', 'type': 'float', 'default': 0},
            {'name': 'kegg_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'kegg_similarity', 'type': 'float', 'default': 0},
            {'name': 'kegg_identity', 'type': 'float', 'default': 0},
            {'name': 'pfam_evalue', 'type': 'float', 'default': 1e-5},
            # input file containing seq_ids that needs retaining
            {'name': 'keep_list', 'type': 'infile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.tools = list()
        self.step.add_steps('nr_filter', 'swissprot_filter', 'kegg_filter', 'eggnog_filter', 'pfam_filter', 'go_filter')

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('db', self.option('db')))
        if self.option('db') == '':
            raise OptionError('required database must be specified')
        if 'nr' in self.option('db'):
            if self.option('blast_nr_xml').is_set:
                self.logger.debug('{} - {}'.format('blast_nr_xml', self.option('blast_nr_xml').prop['path']))
            else:
                raise OptionError('NR BLAST XML must be provided')
        if 'swissprot' in self.option('db'):
            if self.option('blast_swissprot_xml').is_set:
                self.logger.debug('{} - {}'.format('blast_swissprot_xml', self.option('blast_swissprot_xml').prop['path']))
            else:
                raise OptionError('SWISSPROT BLAST XML must be provided')
        if 'eggnog' in self.option('db'):
            if self.option('blast_eggnog_xml').is_set:
                self.logger.debug('{} - {}'.format('blast_eggnog_xml', self.option('blast_eggnog_xml').prop['path']))
            else:
                raise OptionError('EGGNOG BLAST XML must be provided')
        if 'kegg' in self.option('db'):
            if self.option('blast_kegg_xml').is_set:
                self.logger.debug('{} - {}'.format('blast_kegg_xml', self.option('blast_kegg_xml').prop['path']))
            else:
                raise OptionError('KEGG BLAST XML must be provided')
        if 'pfam' in self.option('db'):
            if self.option('pfam_domain').is_set:
                self.logger.debug('{} - {}'.format('pfam_domain', self.option('pfam_domain').prop['path']))
            else:
                raise OptionError('PFAM DOMAIN must be provided')
        if 'go' in self.option('db'):
            if self.option('blast2go_annot').is_set:
                self.logger.debug('{} - {}'.format('blast2go_annot', self.option('blast2go_annot').prop['path']))
            else:
                raise OptionError('BLAST2GO XLS must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        if 'nr' in self.option('db'):
            self.nr_filter = self.add_tool('lnc_rna.annotation.filter_annot')
            self.set_nr_filter()
            self.tools.append(self.nr_filter)
        if 'swissprot' in self.option('db'):
            self.swissprot_filter = self.add_tool('lnc_rna.annotation.filter_annot')
            self.set_swissprot_filter()
            self.tools.append(self.swissprot_filter)
        if 'eggnog' in self.option('db'):
            self.eggnog_filter = self.add_tool('lnc_rna.annotation.filter_annot')
            self.set_eggnog_filter()
            self.tools.append(self.eggnog_filter)
        if 'kegg' in self.option('db'):
            self.kegg_filter = self.add_tool('lnc_rna.annotation.filter_annot')
            self.set_kegg_filter()
            self.tools.append(self.kegg_filter)
        if 'pfam' in self.option('db'):
            self.pfam_filter = self.add_tool('lnc_rna.annotation.filter_annot')
            self.set_pfam_filter()
            self.tools.append(self.pfam_filter)
        if 'go' in self.option('db'):
            self.go_filter = self.add_tool('lnc_rna.annotation.filter_annot')
            self.set_go_filter()
            self.tools.append(self.go_filter)
        super(AnnotFilterModule, self).run()
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def set_nr_filter(self):
        options = {
            'types': 'xml',
            'xml': self.option('blast_nr_xml'),
            'keep_list': self.option('keep_list'),
            'evalue': self.option('nr_evalue'),
            'identity': self.option('nr_identity'),
            'similarity': self.option('nr_similarity')
        }
        self.nr_filter.on('start', self.set_step, {'start': self.step.nr_filter})
        self.nr_filter.on('end', self.set_step, {'end': self.step.nr_filter})
        self.nr_filter.set_options(options)

    def set_swissprot_filter(self):
        options = {
            'types': 'xml',
            'xml': self.option('blast_swissprot_xml'),
            'keep_list': self.option('keep_list'),
            'evalue': self.option('swissprot_evalue'),
            'identity': self.option('swissprot_identity'),
            'similarity': self.option('swissprot_similarity')
        }
        self.swissprot_filter.on('start', self.set_step, {'start': self.step.swissprot_filter})
        self.swissprot_filter.on('end', self.set_step, {'end': self.step.swissprot_filter})
        self.swissprot_filter.set_options(options)

    def set_eggnog_filter(self):
        options = {
            'types': 'xml',
            'xml': self.option('blast_eggnog_xml'),
            'keep_list': self.option('keep_list'),
            'evalue': self.option('eggnog_evalue'),
            'identity': self.option('eggnog_identity'),
            'similarity': self.option('eggnog_similarity')
        }
        self.eggnog_filter.on('start', self.set_step, {'start': self.step.eggnog_filter})
        self.eggnog_filter.on('end', self.set_step, {'end': self.step.eggnog_filter})
        self.eggnog_filter.set_options(options)

    def set_kegg_filter(self):
        options = {
            'types': 'xml',
            'xml': self.option('blast_kegg_xml'),
            'keep_list': self.option('keep_list'),
            'evalue': self.option('kegg_evalue'),
            'identity': self.option('kegg_identity'),
            'similarity': self.option('kegg_similarity')
        }
        self.kegg_filter.on('start', self.set_step, {'start': self.step.kegg_filter})
        self.kegg_filter.on('end', self.set_step, {'end': self.step.kegg_filter})
        self.kegg_filter.set_options(options)

    def set_pfam_filter(self):
        options = {
            'types': 'hmm',
            'hmm': self.option('pfam_domain'),
            'keep_list': self.option('keep_list'),
            'evalue': self.option('pfam_evalue'),
        }
        self.pfam_filter.on('start', self.set_step, {'start': self.step.pfam_filter})
        self.pfam_filter.on('end', self.set_step, {'end': self.step.pfam_filter})
        self.pfam_filter.set_options(options)

    def set_go_filter(self):
        options = {
            'types': 'go',
            'blast2go_annot': self.option('blast2go_annot'),
            'keep_list': self.option('keep_list'),
            'evalue': self.option('nr_evalue'),
            'identity': self.option('nr_identity'),
            'similarity': self.option('nr_similarity')
        }
        self.go_filter.on('start', self.set_step, {'start': self.step.go_filter})
        self.go_filter.on('end', self.set_step, {'end': self.step.go_filter})
        self.go_filter.set_options(options)

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for db in self.option('db').split(','):
            if db == 'nr':
                src = self.nr_filter.output_dir
                dst = os.path.join(self.output_dir, 'nr')
            elif db == 'swissprot':
                src = self.swissprot_filter.output_dir
                dst = os.path.join(self.output_dir, 'swissprot')
            elif db == 'eggnog':
                src = self.eggnog_filter.output_dir
                dst = os.path.join(self.output_dir, 'eggnog')
            elif db == 'kegg':
                src = self.kegg_filter.output_dir
                dst = os.path.join(self.output_dir, 'kegg')
            elif db == 'pfam':
                src = self.pfam_filter.output_dir
                dst = os.path.join(self.output_dir, 'pfam')
            elif db == 'go':
                src = self.go_filter.output_dir
                dst = os.path.join(self.output_dir, 'go')
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            self.logger.info('succeed in copying {} to {}'.format(src, dst))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        super(AnnotFilterModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_filter_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'lnc_rna.annot_filter',
            'instant': False,
            'options': {
                'db': 'nr,swissprot,kegg,eggnog,pfam,go',
                'blast_nr_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/nr/blast.xml',
                'blast_swissprot_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/swissprot/blast.xml',
                'blast_eggnog_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/eggnog/blast.xml',
                'blast_kegg_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/kegg/blast.xml',
                'pfam_domain': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_orfpfam/pfam_domain',
                'blast2go_annot': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annotation/filter_annot/annot_mapdb/GO/blast2go_merge.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()