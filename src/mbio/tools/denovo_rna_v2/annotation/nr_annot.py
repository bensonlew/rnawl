# -*- coding: utf-8 -*-
# __author__ = 'zengjing, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest
from mbio.packages.align.blast.blastout_statistics import *
from mbio.packages.denovo_rna_v2.nr_stat import nr_stat
from mbio.packages.rna.annot_config import AnnotConfig


class NrAnnotAgent(Agent):
    '''
    last_modify: 2019.04.22
    '''
    def __init__(self, parent):
        super(NrAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'blast_xml', 'type': 'infile', 'format': 'ref_rna_v2.blast_xml'},
            {'name': 'longest_t2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'database', 'type': 'string', 'default': 'nr'},
            {'name': 'txpt_table', 'type': 'outfile', 'format': 'ref_rna_v2.blast_table'},
            {'name': 'gene_table', 'type': 'outfile', 'format': 'ref_rna_v2.blast_table'},
            {'name': 'txpt_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {"name": "ncbi_taxonmy_version", "type": "string", "default":"2019"},

        ]
        self.add_option(options)
        self.step.add_steps('nr_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.nr_annot.start()
        self.step.update()

    def step_end(self):
        self.step.nr_annot.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        infile_size = os.path.getsize(self.option('blast_xml').prop['path'])
        self._memory = '{}G'.format(int(float(infile_size) / 1024 ** 3 * 16 + 6))

    def end(self):
        super(NrAnnotAgent, self).end()

class NrAnnotTool(Tool):
    def __init__(self, config):
        super(NrAnnotTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.xml2table_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/xml2table.py')
        self.txml2gxml_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/txml2gxml.py')
        self.get_venn_list_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/get_venn_list.py')
        self.blast2taxon_py = os.path.join(self.config.PACKAGE_DIR, 'denovo_rna_v2/nr_taxon_stat.py')
        self.taxon_tree_py = self.config.PACKAGE_DIR + "/denovo_rna_v2/taxon_treeview.py"
        self.prefix = os.path.basename(self.option('blast_xml').path)[:-4]
        self.txpt_xml = os.path.join(self.work_dir, '{}.T.xml'.format(self.prefix))
        self.gene_xml = os.path.join(self.work_dir, '{}.G.xml'.format(self.prefix))
        self.txpt_tsv = os.path.join(self.output_dir, '{}.T.tsv'.format(self.option('database')))
        self.gene_tsv = os.path.join(self.output_dir, '{}.G.tsv'.format(self.option('database')))
        self.txpt_list = os.path.join(self.output_dir, '{}.T.list'.format(self.option('database')))
        self.gene_list = os.path.join(self.output_dir, '{}.G.list'.format(self.option('database')))
        self.ncbi_taxon_db = AnnotConfig().get_file_path(version=self.option("ncbi_taxonmy_version"), db="ncbi_tax", file="taxonomy.db")
        self.acc2taxon_db = AnnotConfig().get_file_path(version=self.option("ncbi_taxonmy_version"), db="ncbi_tax", file="accession2taxid.db")


    def run(self):
        super(NrAnnotTool, self).run()
        self.run_txml2gxml()
        self.run_xml2table_t()
        self.run_xml2table_g()
        self.run_get_venn_list_t()
        self.run_get_venn_list_g()
        self.wait()
        self.run_nr2taxon_t()
        self.run_nr2taxon_g()
        self.wait()

        self.run_taxon2tree_t()
        self.run_taxon2tree_g()
        self.wait()
        self.run_nr_stat()
        self.set_output()
        self.end()

    def run_nr2taxon_g(self):
        cmd = '{} {} {} {} {} {}'.format(self.python, self.blast2taxon_py, self.txpt_tsv, self.txpt_tsv + "species.txt", self.ncbi_taxon_db, self.acc2taxon_db)
        cmd_name = 'run_nr2taxon_g'
        self.run_code(cmd_name, cmd)

    def run_nr2taxon_t(self):
        cmd = '{} {} {} {} {} {}'.format(self.python, self.blast2taxon_py, self.gene_tsv, self.gene_tsv + "species.txt", self.ncbi_taxon_db, self.acc2taxon_db)
        cmd_name = 'run_nr2taxon_t'
        self.run_code(cmd_name, cmd)

    def run_taxon2tree_t(self):
        cmd = "{} {} {}".format(self.python, self.taxon_tree_py, self.txpt_tsv + "species.txt")
        cmd_name = 'run_taxon2tree_t'
        self.run_code(cmd_name, cmd)

    def run_taxon2tree_g(self):
        cmd = "{} {} {}".format(self.python, self.taxon_tree_py, self.gene_tsv + "species.txt")
        cmd_name = 'run_taxon2tree_g'
        self.run_code(cmd_name, cmd)


    def run_nr_stat(self):
        try:
            blastout_statistics(blast_table=self.txpt_tsv, evalue_path=self.output_dir + '/trans_nr_evalue.xls', similarity_path=self.output_dir + '/trans_nr_similar.xls')
            blastout_statistics(blast_table=self.gene_tsv, evalue_path=self.output_dir + '/gene_nr_evalue.xls', similarity_path=self.output_dir + '/gene_nr_similar.xls')
            self.logger.info("End: evalue,similar for gene nr blast table ")
        except Exception as e:
            self.set_error("运行nr evalue,similar for gene nr blast table出错:%s", variables = (e), code = "32001612")
            self.logger.info("Error: evalue,similar for gene nr blast table")
        try:
            nr_stat().nr_stat_info2(
                tran_taxonfile = self.txpt_tsv + "species.txt",
                outpath = self.output_dir + '/tran_nr_species_stat.xls'
            )
            nr_stat().nr_stat_info2(
                tran_taxonfile = self.gene_tsv + "species.txt",
                outpath = self.output_dir + '/gene_nr_species_stat.xls'
            )
            self.logger.info("nr taxomy End ")
        except Exception as e:
            self.set_error("运行nr taxomy出错:%s", variables = (e), code = "32001613")
            self.logger.info("Error: taxomy blast table")



    def run_txml2gxml(self):
        shutil.copy(self.option('blast_xml').path, self.txpt_xml)
        cmd = '{} {}'.format(self.python, self.txml2gxml_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --xml {}'.format(self.txpt_xml)
        cmd += ' --output {}'.format(self.gene_xml)
        cmd_name = 'run_txml2gxml'
        self.run_code(cmd_name, cmd)

    def run_xml2table_t(self):
        cmd = '{} {} -a -t top'.format(self.python, self.xml2table_py)
        cmd += ' -i {}'.format(self.txpt_xml)
        cmd += ' -o {}'.format(self.txpt_tsv)
        cmd_name = 'run_xml2table_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_xml2table_g(self):
        cmd = '{} {} -a -t top'.format(self.python, self.xml2table_py)
        cmd += ' -i {}'.format(self.gene_xml)
        cmd += ' -o {}'.format(self.gene_tsv)
        cmd_name = 'run_xml2table_g'
        self.run_code(cmd_name, cmd)

    def run_get_venn_list_t(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.txpt_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.txpt_list)
        cmd_name = 'run_get_venn_list_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_get_venn_list_g(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.gene_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.gene_list)
        cmd_name = 'run_get_venn_list_g'
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
        self.option('txpt_table').set_path(self.txpt_tsv)
        self.option('gene_table').set_path(self.gene_tsv)
        self.option('txpt_list').set_path(self.txpt_list)
        self.option('gene_list').set_path(self.gene_list)
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
            'id': 'nr_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.annotation.nr_annot',
            'instant': False,
            'options': {
                'blast_xml': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotClass/AnnotFile/result/blast.filter.xml',
                'longest_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotClass/AnnotFile/result/longest.t2g.tsv',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
