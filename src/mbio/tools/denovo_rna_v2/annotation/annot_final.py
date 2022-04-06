# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class AnnotFinalAgent(Agent):
    '''
    last_modify: 2019.04.22
    '''
    def __init__(self, parent):
        super(AnnotFinalAgent, self).__init__(parent)
        options = [
            {'name': 't2g2r2l2p', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'loc2db2type', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'nr_table', 'type': 'infile', 'format': 'ref_rna_v2.blast_table'}, # yes
            {'name': 'swissprot_table', 'type': 'infile', 'format': 'ref_rna_v2.blast_table'}, # yes
            {'name': 'cog_table', 'type': 'infile', 'format': 'ref_rna_v2.cog_table'}, # yes
            {'name': 'kegg_table', 'type': 'infile', 'format': 'ref_rna_v2.kegg_table'}, # yes
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'go_list', 'type': 'infile', 'format': 'ref_rna_v2.go_list'}, # yes
            {'name': 'des', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'enterz', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'des_type', 'type': 'string', 'default': None}, # yes
            {'name': 'statistics', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'query', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {"name": "kegg_version", "type": "string", "default": "202003"}
        ]
        self.add_option(options)
        self.step.add_steps('annot_final')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.annot_final.start()
        self.step.update()

    def stepfinish(self):
        self.step.annot_final.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '16G'

    def end(self):
        super(AnnotFinalAgent, self).end()

class AnnotFinalTool(Tool):
    def __init__(self, config):
        super(AnnotFinalTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.statistics_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/statistics.py')
        self.query_py = os.path.join(self.config.PACKAGE_DIR, 'denovo_rna_v2/new_annotation_query2.py')
        self.statistics_tsv = os.path.join(self.output_dir, 'statistics.tsv')
        self.query_tsv = os.path.join(self.output_dir, 'query.tsv')

    def run(self):
        super(AnnotFinalTool, self).run()
        self.run_query()
        self.run_statistics()
        self.set_output()
        self.end()

    def run_query(self):
        gene2trans = self.option('t2g2r2l2p').prop['path']
        tran_outpath = os.path.join(self.output_dir, 'query.tsv')
        gene_outpath = os.path.join(self.output_dir, 'none.txt')
        new_gtf_path = 'None'
        ref_gtf_path = 'None'
        length_path = 'None'
        gene_file = 'None'
        cog_table = self.option('cog_table').prop['path']
        kegg_table = self.option('kegg_table').prop['path']
        go_list = self.option('go_list').prop['path']
        nr_table = self.option('nr_table').prop['path']
        swissprot_table = self.option('swissprot_table').prop['path']
        pfam_domain = self.option('pfam_domain').prop['path']
        '''
        des = self.option('des').prop['path']
        subloc = 'None'
        enterz = self.option('enterz').prop['path']
        des_type = self.option('des_type')
        '''
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            self.python,
            self.query_py,
            gene2trans,
            tran_outpath,
            gene_outpath,
            new_gtf_path,
            ref_gtf_path,
            length_path,
            gene_file,
            cog_table,
            kegg_table,
            go_list,
            nr_table,
            swissprot_table,
            pfam_domain,
            self.option("kegg_version")
        )
        cmd_name = 'run_query'
        self.run_code(cmd_name, cmd, block=False)

    def run_statistics(self):
        cmd = '{} {}'.format(self.python, self.statistics_py)
        cmd += ' --input {}'.format(self.option('loc2db2type').path)
        cmd += ' --map {}'.format(self.option('t2g2r2l2p').path)
        cmd += ' --output {}'.format(self.statistics_tsv)
        cmd_name = 'run_statistics'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run {}, abord'.format(n))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('query').set_path(self.query_tsv)
        self.option('statistics').set_path(self.statistics_tsv)
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
            'id': 'annot_final_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.annotation.annot_final',
            'instant': False,
            'options': {
                't2g2r2l2p': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/output/all_tran2gene.txt',
                'loc2db2type': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/loc2db2type.tsv',
                'nr_table': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/output/nr/blast.filter.T.tsv',
                'swissprot_table': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/output/swissprot/blast.filter.T.tsv',
                'cog_table': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/output/cog/cog.tsv',
                'kegg_table': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/output/kegg/T/kegg_table.tsv',
                'pfam_domain': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/output/pfam/pfam.filter.T.tsv',
                'go_list': '/mnt/ilustre/users/sanger-dev/workspace/20190418/Single_annot_class_beta_2665_9026/AnnotClassBeta/output/go/id2terms.T.tsv',
                'des': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Test_ref/Ensemble_release_89/biomart/Oreochromis_niloticus.Orenil1.0.biomart_gene.txt',
                'des_type': 'type1',
                'enterz': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Test_ref/Ensemble_release_89/NCBI/Oreochromis_niloticus.Orenil1.0.biomart_enterz.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
