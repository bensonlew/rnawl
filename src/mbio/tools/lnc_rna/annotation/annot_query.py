# -*- coding: utf-8 -*-
# __author__ = 'zengjing, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class AnnotQueryAgent(Agent):
    '''
    last_modify: 2019.02.19
    '''
    def __init__(self, parent):
        super(AnnotQueryAgent, self).__init__(parent)
        options = [
            {'name': 'gene2trans', 'type': 'infile', 'format': 'lnc_rna.common'}, # yes
            {'name': 'new_gtf_path', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'ref_gtf_path', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'length_path', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'gene_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'cog_list', 'type': 'infile', 'format': 'lnc_rna.cog_table'}, # yes
            {'name': 'kegg_table', 'type': 'infile', 'format': 'lnc_rna.kegg_table'}, # yes
            {'name': 'gos_list', 'type': 'infile', 'format': 'lnc_rna.go_list'}, # yes
            {'name': 'blast_nr_table', 'type': 'infile', 'format': 'lnc_rna.blast_table'}, # yes
            {'name': 'blast_swissprot_table', 'type': 'infile', 'format': 'lnc_rna.blast_table'}, # yes
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'lnc_rna.common'}, # yes
            {'name': 'des', 'type': 'infile', 'format': 'lnc_rna.common'}, # yes
            {'name': 'subloc', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'enterz', 'type': 'infile', 'format': 'lnc_rna.common'}, # yes
            {'name': 'des_type', 'type': 'string', 'default': 'type3'}, # yes
            {"name": "kegg_version", "type": "string", "default": "202003"},
        ]
        self.add_option(options)
        self.step.add_steps('annot_query')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.annot_query.start()
        self.step.update()

    def stepfinish(self):
        self.step.annot_query.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('gene2trans').is_set:
            self.logger.debug('{} = {}'.format('gene2trans', self.option('gene2trans').prop['path']))
        else:
            raise OptionError('input isoform2unigene file must be provided')
        if self.option('cog_list').is_set:
            self.logger.debug('{} = {}'.format('cog_list', self.option('cog_list').prop['path']))
        if self.option('kegg_table').is_set:
            self.logger.debug('{} = {}'.format('kegg_table', self.option('kegg_table').prop['path']))
        if self.option('gos_list').is_set:
            self.logger.debug('{} = {}'.format('gos_list', self.option('gos_list').prop['path']))
        if self.option('blast_nr_table').is_set:
            self.logger.debug('{} = {}'.format('blast_nr_table', self.option('blast_nr_table').prop['path']))
        if self.option('blast_swissprot_table').is_set:
            self.logger.debug('{} = {}'.format('blast_swissprot_table', self.option('blast_swissprot_table').prop['path']))
        if self.option('pfam_domain').is_set:
            self.logger.debug('{} = {}'.format('pfam_domain', self.option('pfam_domain').prop['path']))
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('des').is_set:
            self.logger.debug('{} = {}'.format('des', self.option('des').prop['path']))
        else:
            raise OptionError('input des file must be provided')
        if self.option('enterz').is_set:
            self.logger.debug('{} = {}'.format('enterz', self.option('enterz').prop['path']))
        else:
            raise OptionError('input enterz file must be provided')
        if self.option('des_type') != '':
            self.logger.debug('{} = {}'.format('des_type', self.option('des_type')))
        else:
            raise OptionError('des_type must be specified')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '16G'

    def end(self):
        super(AnnotQueryAgent, self).end()

class AnnotQueryTool(Tool):
    def __init__(self, config):
        super(AnnotQueryTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.annotation_query_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/annotation_query.py')

    def run(self):
        super(AnnotQueryTool, self).run()
        self.run_query()
        self.end()

    def run_query(self):
        gene2trans = self.option('gene2trans').prop['path']
        tran_outpath = os.path.join(self.output_dir, 'all_anno_detail.xls')
        gene_outpath = os.path.join(self.output_dir, 'gene_anno_detail.xls')
        if self.option('new_gtf_path').is_set:
            new_gtf_path = self.option('new_gtf_path').prop['path']
        else:
            new_gtf_path = 'None'
        if self.option('ref_gtf_path').is_set:
            ref_gtf_path = self.option('ref_gtf_path').prop['path']
        else:
            ref_gtf_path = 'None'
        if self.option('length_path').is_set:
            length_path = self.option('length_path').prop['path']
        else:
            length_path = 'None'
        if self.option('gene_file').is_set:
            gene_file = self.option('gene_file').prop['path']
        else:
            gene_file = 'None'
        if self.option('cog_list').is_set:
            cog_list = self.option('cog_list').prop['path']
        else:
            cog_list = 'None'
        if self.option('kegg_table').is_set:
            kegg_table = self.option('kegg_table').prop['path']
        else:
            kegg_table = 'None'

        if self.option('gos_list').is_set:
            gos_list = self.option('gos_list').prop['path']
        else:
            gos_list = 'None'
        if self.option('blast_nr_table').is_set:
            blast_nr_table = self.option('blast_nr_table').prop['path']
        else:
            blast_nr_table = 'None'
        if self.option('blast_swissprot_table').is_set:
            blast_swissprot_table = self.option('blast_swissprot_table').prop['path']
        else:
            blast_swissprot_table = 'None'
        if self.option('pfam_domain').is_set:
            pfam_domain = self.option('pfam_domain').prop['path']
        else:
            pfam_domain = 'None'
        des = self.option('des').prop['path']
        if self.option('subloc').is_set:
            subloc = self.option('subloc').prop['path']
        else:
            subloc = 'None'
        enterz = self.option('enterz').prop['path']
        des_type = self.option('des_type')
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            self.python,
            self.annotation_query_py,
            gene2trans,
            tran_outpath,
            gene_outpath,
            new_gtf_path,
            ref_gtf_path,
            length_path,
            gene_file,
            cog_list,
            kegg_table,
            gos_list,
            blast_nr_table,
            blast_swissprot_table,
            pfam_domain,
            des,
            subloc,
            enterz,
            des_type,
            self.option("kegg_version")
        )
        cmd_name = 'run_query'
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

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_query_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.annot_query',
            'instant': False,
            'options': {
                'gene2trans': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/annot_file/i2u',
                'blast_nr_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/nr_annot/blast.filter.xls',  #
                'blast_swissprot_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/swissprot_annot/blast.filter.xls', #
                'pfam_domain': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/pfam/pfam.filter.xls',
                'gos_list': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/go_annot/query_gos.list',
                'kegg_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/kegg_annot/kegg_table.xls',
                'cog_list': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_class/eggnog_annot/cog.xls',
                'des': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                'des_type': 'type1',
                'enterz': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/NCBI/Homo_sapiens.GRCh38.biomart_enterz.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
