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
            {'name': 'nr_table', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'uniprot_table', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'cog_table', 'type': 'infile', 'format': 'ref_rna_v2.cog_table'}, # yes
            {'name': 'kegg_table', 'type': 'infile', 'format': 'ref_rna_v2.kegg_table'}, # yes
            {'name': 'kegg_table_spe', 'type': 'infile', 'format': 'ref_rna_v2.kegg_table'}, # yes
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'go_list', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'des', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'enterz', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'des_type', 'type': 'string', 'default': None}, # yes
            {'name': 'statistics', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'nr_table_gene', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'uniprot_table_gene', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'cog_table_gene', 'type': 'infile', 'format': 'ref_rna_v2.cog_table'}, # yes
            {'name': 'kegg_table_gene', 'type': 'infile', 'format': 'ref_rna_v2.kegg_table'}, # yes
            {'name': 'kegg_table_gene_spe', 'type': 'infile', 'format': 'ref_rna_v2.kegg_table'}, # yes
            {'name': 'pfam_domain_gene', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'go_list_gene', 'type': 'infile', 'format': 'ref_rna_v2.go_list'}, # yes
            {'name': 'do_gene', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'reactome_gene', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes
            {'name': 'disgenet_gene', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # yes

            {'name': 'query', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'query_gene', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {"name": "kegg_version", "type": "string", "default": "202003"},
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
        self.python = 'program/Python/bin/python'
        self.statistics_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/statistics.py')
        self.query_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/query.py')
        self.query_gene_py = os.path.join(self.config.PACKAGE_DIR, 'ref_genome_db_medical/annotation/query_gene.py')
        self.statistics_tsv = os.path.join(self.output_dir, 'statistics.tsv')
        self.query_tsv = os.path.join(self.output_dir, 'query.tsv')

    def run(self):
        super(AnnotFinalTool, self).run()
        self.run_query()
        self.run_query_gene()
        self.run_statistics()
        self.run_statistics_gene()
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
        kegg_table += "," + self.option('kegg_table_spe').prop['path']
        go_list = self.option('go_list').prop['path']
        nr_table = self.option('nr_table').prop['path']
        uniprot_table = self.option('uniprot_table').prop['path']
        pfam_domain = self.option('pfam_domain').prop['path']
        des = self.option('des').prop['path']
        subloc = 'None'
        enterz = self.option('enterz').prop['path']
        des_type = self.option('des_type')
        kegg_version = self.option('kegg_version')
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
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
            uniprot_table,
            pfam_domain,
            des,
            subloc,
            enterz,
            des_type,
            kegg_version
        )
        cmd_name = 'run_query'
        self.run_code(cmd_name, cmd, block=False)

    def run_query_gene(self):
        gene2trans = self.option('t2g2r2l2p').prop['path']
        tran_outpath = os.path.join(self.output_dir, 'query_gene.tsv')
        gene_outpath = os.path.join(self.output_dir, 'none.txt')
        new_gtf_path = 'None'
        ref_gtf_path = 'None'
        length_path = 'None'
        gene_file = 'None'
        cog_table = self.option('cog_table_gene').prop['path']
        kegg_table = self.option('kegg_table_gene').prop['path']
        kegg_table += "," + self.option('kegg_table_gene_spe').prop['path']
        go_list = self.option('go_list_gene').prop['path']
        nr_table = self.option('nr_table_gene').prop['path']
        uniprot_table = self.option('uniprot_table_gene').prop['path']
        pfam_domain = self.option('pfam_domain_gene').prop['path']
        des = self.option('des').prop['path']
        subloc = 'None'
        enterz = self.option('enterz').prop['path']
        des_type = self.option('des_type')
        kegg_version = self.option('kegg_version')
        do_list = self.option('do_gene').prop['path']
        reactome_list = self.option('reactome_gene').prop['path']
        disgenets_list = self.option('disgenet_gene').prop['path']
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            self.python,
            self.query_gene_py,
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
            uniprot_table,
            pfam_domain,
            des,
            subloc,
            enterz,
            do_list,
            reactome_list,
            disgenets_list

        )
        cmd_name = 'run_query_gene'
        self.run_code(cmd_name, cmd, block=False)

    def run_statistics(self):
        cmd = '{} {}'.format(self.python, self.statistics_py)
        cmd += ' --input {}'.format(self.option('loc2db2type').path)
        cmd += ' --map {}'.format(self.option('t2g2r2l2p').path)
        cmd += ' --output {}'.format(self.statistics_tsv)
        cmd_name = 'run_statistics'
        self.run_code(cmd_name, cmd)

    def run_statistics_gene(self):
        cmd = '{} {}'.format(self.python, self.statistics_py)
        cmd += ' --input {}'.format(self.option('loc2db2type').path)
        cmd += ' --map {}'.format(self.option('t2g2r2l2p').path)
        cmd += ' --output {}'.format(self.statistics_tsv)
        cmd_name = 'run_statistics_gene'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        # print "command {}".format(dir(command))
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run %s, abord', variables=(n), code="33709402")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('query').set_path(self.query_tsv)
        self.option('query_gene').set_path(os.path.join(self.output_dir, 'query_gene.tsv'))
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
            'name': 'ref_genome_db_medical.annotation.annot_final',
            'instant': False,
            'options': {
                't2g2r2l2p': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/begin/t2g2r2l2p.tsv',
                'loc2db2type': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/loc.x.list',
                'nr_table': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/nr_annot/annot_nr.T.ids',
                'uniprot_table': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/uniprot_annot/annot_uniprot.T.ids',
                'cog_table': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/egg_nog_annot/cog.tsv',
                'kegg_table': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/kegg_annot/T/kegg_table.tsv',
                'pfam_domain': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/pfam_annot/pfam.T.tsv',
                'go_list': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/go_annot/id2terms.T.tsv',
                'des': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/biomart/biomart1.txt',
                'des_type': 'type1',
                'enterz': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/NCBI/entrenz.txt',
                'nr_table_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/nr_annot/annot_nr.G.ids',
                'uniprot_table_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/uniprot_annot/annot_uniprot.G.ids',
                'cog_table_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/egg_nog_annot/cog.G.tsv',
                'kegg_table_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/kegg_annot/G/kegg_table.tsv',
                'pfam_domain_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/pfam_annot/pfam.G.tsv',
                'go_list_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/go_annot/id2terms.G.tsv',
                'do_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/do_annot/G/do_detail.xls',
                'reactome_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/reactome_annot/reactome.G.detail.tsv',
                'disgenet_gene': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/disgenet_annot/disgenet.G.tsv',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
