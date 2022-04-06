# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest
from mbio.packages.rna.annot_config import AnnotConfig
from biocluster.option import Option

class AnnotMergeAgent(Agent):
    '''
    last_modify: 2020.08.20
    '''
    def __init__(self, parent):
        super(AnnotMergeAgent, self).__init__(parent)
        options = [
            {'name': 'ref_class_dir', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},
            {'name': 'new_class_dir', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},
            {'name': 'database', 'type': 'string', 'default': 'nr,uniprot,cog,kegg,pfam,go'},
            {'name': 'new_mapdb_dir', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},
            {'name': 'is_assemble', 'type': 'bool', 'default': True},
            {'name': 'kegg_version', 'type': 'string', 'default':"202007"},
            {"name": "kegg_species", "type": "string", "default": None},
        ]
        self.add_option(options)
        self._memory_increase_step = 40
        self.queue = 'BLAST2GO'
        self.step.add_steps('annot_merge')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.annot_merge.start()
        self.step.update()

    def step_end(self):
        self.step.annot_merge.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '50G'

    def end(self):
        super(AnnotMergeAgent, self).end()

class AnnotMergeTool(Tool):
    def __init__(self, config):
        super(AnnotMergeTool, self).__init__(config)
        # 临时测试
        if "kegg_species" not in self._options:
            self._options["kegg_species"] = Option({"name": "kegg_species", "type": "string", "default": None})
        # software and script
        self.python = 'miniconda2/bin/python'
        self.merge_py = os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/annotation/merge.py')
        self.merge_query_gene_py = os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/annotation/merge_gene_query.py')
        self.go_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/go_annotation.py')
        self.go_annotation2_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/go_annotation2.py')
        self.kegg_merge_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/kegg_merge.py')
        if self.option("kegg_species"):
            self.kegg_merge_py = os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/annotation/kegg_merge_spe.py')

        self.statistics_py = os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/annotation/statistics.py')
        self.map4_r = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/map4.r')
        self.rscript = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.3/bin/Rscript')
        self.image_magick_convert = os.path.join(self.config.SOFTWARE_DIR, 'program/ImageMagick/bin/convert')

        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))
        self.map_html = self.kegg_files_dict['html'] + '/'

        # output
        self.ref_dir = os.path.join(self.output_dir, 'refannot_class')
        self.new_dir = os.path.join(self.output_dir, 'newannot_class')
        self.all_dir = os.path.join(self.output_dir, 'allannot_class')
        self.ref_map = os.path.join(self.ref_dir, 'all_tran2gene.txt')
        self.new_map = os.path.join(self.new_dir, 'all_tran2gene.txt')
        self.all_map = os.path.join(self.new_dir, 'all_tran2gene.txt')
        # kegg
        self.ref_kegg_pathway_tran_xls = os.path.join(self.ref_dir, 'kegg/kegg_pathway_tran.xls')
        self.new_kegg_pathway_tran_xls = os.path.join(self.new_dir, 'kegg/kegg_pathway_tran.xls')
        self.all_kegg_pathway_tran_xls = os.path.join(self.all_dir, 'kegg/kegg_pathway_tran.xls')
        self.all_kegg_pathway_tran_dir = os.path.join(self.all_dir, 'kegg/kegg_pathway_tran_dir')
        self.all_kegg_gene_tran_xls = os.path.join(self.all_dir, 'kegg/kegg_gene_tran.xls')
        self.ref_kegg_pathway_gene_xls = os.path.join(self.ref_dir, 'kegg/kegg_pathway_gene.xls')
        self.new_kegg_pathway_gene_xls = os.path.join(self.new_dir, 'kegg/kegg_pathway_gene.xls')
        self.all_kegg_pathway_gene_xls = os.path.join(self.all_dir, 'kegg/kegg_pathway_gene.xls')
        self.all_kegg_pathway_gene_dir = os.path.join(self.all_dir, 'kegg/kegg_pathway_gene_dir')
        self.all_kegg_gene_gene_xls = os.path.join(self.all_dir, 'kegg/kegg_gene_gene.xls')

        if self.option("kegg_species"):
            self.ref_kegg_pathway_tran_xls = os.path.join(self.ref_dir, 'kegg/spe_T/spe_pathway_table.tsv')
            self.new_kegg_pathway_tran_xls = os.path.join(self.new_dir, 'kegg/spe_T/spe_pathway_table.tsv')
            self.all_kegg_pathway_tran_xls = os.path.join(self.all_dir, 'kegg/kegg_pathway_tran.xls')
            self.all_kegg_pathway_tran_dir = os.path.join(self.all_dir, 'kegg/kegg_pathway_tran_dir')
            self.all_kegg_gene_tran_xls = os.path.join(self.all_dir, 'kegg/kegg_gene_tran.xls')
            self.ref_kegg_pathway_gene_xls = os.path.join(self.ref_dir, 'kegg/spe_G/spe_pathway_table.tsv')
            self.new_kegg_pathway_gene_xls = os.path.join(self.new_dir, 'kegg/spe_G/spe_pathway_table.tsv')
            self.all_kegg_pathway_gene_xls = os.path.join(self.all_dir, 'kegg/kegg_pathway_gene.xls')
            self.all_kegg_pathway_gene_dir = os.path.join(self.all_dir, 'kegg/kegg_pathway_gene_dir')
            self.all_kegg_gene_gene_xls = os.path.join(self.all_dir, 'kegg/kegg_gene_gene.xls')

        # go
        self.b2g_host = 'localhost'
        self.b2g_user = 'biocluster102'
        self.b2g_passwd = 'sanger-dev-123'
        self.b2g_db = 'b2gdb'
        self.go_txpt_dir = os.path.join(self.work_dir, 'go_txpt')
        self.go_gene_dir = os.path.join(self.work_dir, 'go_gene')
        self.id2terms_txpt_tsv = os.path.join(self.all_dir, 'go/go_list_tran.xls')
        self.id2terms_gene_tsv = os.path.join(self.all_dir, 'go/go_list_gene.xls')
        # useless here, for modify_output in workflow
        self.new_mapdb_dir = os.path.join(self.output_dir, 'newannot_mapdb')
        # reset directory
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        shutil.copytree(self.option('ref_class_dir').path, self.ref_dir)
        if self.option('is_assemble'):
            shutil.copytree(self.option('new_class_dir').path, self.new_dir)
            shutil.copytree(self.option('new_mapdb_dir').path, self.new_mapdb_dir)
            for database in self.option('database').split(','):
                os.makedirs(os.path.join(self.all_dir, database))
        else:
            shutil.copytree(self.ref_dir, self.all_dir)

    def run(self):
        super(AnnotMergeTool, self).run()
        if self.option('is_assemble'):
            self.merge_begin()
            self.merge_final()
            self.set_output()
        self.end()

    def merge_begin(self):
        self.merge_nr()
        self.merge_uniprot()
        self.merge_cog()
        self.merge_kegg()
        self.merge_pfam()
        self.merge_go()
        self.merge_map()
        do_file = os.path.join(self.ref_dir, "do/G/do_detail.xls")
        if os.path.exists(do_file) and os.path.getsize(do_file) != 0 :
            self.merge_do()
        self.merge_reactome()
        disgenet_file = os.path.join(self.ref_dir, "disgenet/disgenet.G.tsv")
        if os.path.exists(disgenet_file) and os.path.getsize(disgenet_file) != 0 :
            self.merge_disgenet()

    def merge_final(self):
        self.run_statistics()
        self.run_kegg_merge_t()
        self.run_kegg_merge_g()
        self.run_go_annotation_t()
        self.run_go_annotation_g()
        self.merge_query()

    def merge_nr(self):
        args_list = [
            ('nr/nr_tran.xls', 'nr/nr_blast_tran.xls', 'nr/nr_blast_tran.xls', 'yes'),
            ('nr/nr_gene.xls', 'nr/nr_blast_gene.xls', 'nr/nr_blast_gene.xls', 'yes'),
            ('nr/nr_venn_tran.txt', 'nr/nr_venn_tran.txt', 'nr/nr_venn_tran.txt', 'no'),
            ('nr/nr_venn_gene.txt', 'nr/nr_venn_gene.txt', 'nr/nr_venn_gene.txt', 'no')
        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_nr_{}'.format(n)
            self.run_code(cmd_name, cmd, block=False)

    def merge_uniprot(self):
        args_list = [
            ('uniprot/uniprot_tran.xls', 'uniprot/uniprot_blast_tran.xls', 'uniprot/uniprot_blast_tran.xls', 'yes'),
            ('uniprot/uniprot_gene.xls', 'uniprot/uniprot_blast_gene.xls', 'uniprot/uniprot_blast_gene.xls', 'yes'),
            ('uniprot/uniprot_venn_tran.txt', 'uniprot/uniprot_venn_tran.txt', 'uniprot/uniprot_venn_tran.txt', 'no'),
            ('uniprot/uniprot_venn_gene.txt', 'uniprot/uniprot_venn_gene.txt', 'uniprot/uniprot_venn_gene.txt', 'no')
        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_uniprot_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_cog(self):
        args_list = [
            ('cog/cog_list_tran.xls', 'cog/cog_list_tran.xls', 'cog/cog_list_tran.xls', 'yes'),
            ('cog/cog_venn_tran.txt', 'cog/cog_venn_tran.txt', 'cog/cog_venn_tran.txt', 'no'),
            ('cog/cog_venn_gene.txt', 'cog/cog_venn_gene.txt', 'cog/cog_venn_gene.txt', 'no')
        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_cog_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_kegg(self):
        args_list = [
            ('kegg/kegg_gene_tran.xls', 'kegg/kegg_gene_tran.xls', 'kegg/kegg_gene_tran.xls', 'yes'),
            ('kegg/kegg_gene_gene.xls', 'kegg/kegg_gene_gene.xls', 'kegg/kegg_gene_gene.xls', 'yes'),
            ('kegg/kegg_venn_tran.txt', 'kegg/kegg_venn_tran.txt', 'kegg/kegg_venn_tran.txt', 'no'),
            ('kegg/kegg_venn_gene.txt', 'kegg/kegg_venn_gene.txt', 'kegg/kegg_venn_gene.txt', 'no')
        ]

        if self.option("kegg_species"):
            args_list = args_list[2:] + [
                ('kegg/spe_T/spe_kegg_table.tsv', 'kegg/spe_T/spe_kegg_table.tsv', 'kegg/kegg_gene_tran.xls', 'yes'),
                ('kegg/spe_G/spe_kegg_table.tsv', 'kegg/spe_G/spe_kegg_table.tsv', 'kegg/kegg_gene_gene.xls', 'yes')
            ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_kegg_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_pfam(self):
        args_list = [
            ('pfam/pfam_domain_tran.xls', 'pfam/pfam_domain_tran.xls', 'pfam/pfam_domain_tran.xls', 'yes'),
            ('pfam/pfam_domain_gene.xls', 'pfam/pfam_domain_gene.xls', 'pfam/pfam_domain_gene.xls', 'yes'),
            ('pfam/pfam_venn_tran.txt', 'pfam/pfam_venn_tran.txt', 'pfam/pfam_venn_tran.txt', 'no'),
            ('pfam/pfam_venn_gene.txt', 'pfam/pfam_venn_gene.txt', 'pfam/pfam_venn_gene.txt', 'no')
        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_pfam_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_go(self):
        args_list = [
            ('go/go_list_tran.xls', 'go/go_list_tran.xls', 'go/go_list_tran.xls', 'no'),
            ('go/go_list_gene.xls', 'go/go_list_gene.xls', 'go/go_list_gene.xls', 'no'),
            ('go/go_venn_tran.txt', 'go/go_venn_tran.txt', 'go/go_venn_tran.txt', 'no'),
            ('go/go_venn_gene.txt', 'go/go_venn_gene.txt', 'go/go_venn_gene.txt', 'no')
        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_go_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_reactome(self):
        args_list = [
            ('reactome/reactome.G.detail.tsv', '', 'reactome/reactome.G.detail.tsv', 'yes'),
            ('reactome/reactome.G.list', '', 'reactome/reactome.G.list', 'no'),
            ('reactome/reactome.T.list', '', 'reactome/reactome.T.list', 'no'),
        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_reactome_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_do(self):
        args_list = [
            ('do/id2terms.G.tsv', '', 'do/id2terms.G.tsv', 'no'),
            ('do/id2terms.G.tsv', '', 'do/id2terms.G.tsv', 'no'),
            ('do/do.G.list', '', 'do/do.G.list', 'no'),
            ('do/do.T.list', '', 'do/do.T.list', 'no'),
            ('do/G/do_detail.xls', '', 'do/G/do_detail.xls', 'yes')
        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_do_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_disgenet(self):
        args_list = [
            ('disgenet/disgenet.G.tsv', '', 'disgenet/disgenet.G.tsv', 'yes'),
            ('disgenet/disgenet.G.list', '', 'disgenet/disgenet.G.list', 'no'),
            ('disgenet/disgenet.T.list', '', 'disgenet/disgenet.T.list', 'no'),

        ]
        for n, args in enumerate(args_list):
            cmd = '{} {}'.format(self.python, self.merge_py)
            cmd += ' --ref {} --new {} --output {} --head {}'.format(*self.proc_args(args))
            cmd_name = 'merge_disgenet_{}'.format(n)
            self.run_code(cmd_name, cmd, block=True)

    def merge_map(self):
        cmd = '{} {} --head no'.format(self.python, self.merge_py)
        cmd += ' --ref {}'.format(os.path.join(self.ref_dir, 'all_tran2gene.txt'))
        cmd += ' --new {}'.format(os.path.join(self.new_dir, 'all_tran2gene.txt'))
        cmd += ' --output {}'.format(os.path.join(self.all_dir, 'all_tran2gene.txt'))
        cmd_name = 'merge_map'
        self.run_code(cmd_name, cmd)

    def run_statistics(self):
        lines = list()
        for database in self.option('database').split(','):
            lines.extend([
                '{}\t{}\t{}\n'.format(os.path.join(self.all_dir, '{}/{}_venn_tran.txt'.format(database, database)), database, 'transcript'),
                '{}\t{}\t{}\n'.format(os.path.join(self.all_dir, '{}/{}_venn_gene.txt'.format(database, database)), database, 'gene')
            ])
        else:
            loc2db2type = os.path.join(self.work_dir, 'loc2db2type.tsv')
            open(loc2db2type, 'w').writelines(lines)
        cmd = '{} {}'.format(self.python, self.statistics_py)
        cmd += ' --input {}'.format(loc2db2type)
        cmd += ' --map {}'.format(os.path.join(self.all_dir, 'all_tran2gene.txt'))
        cmd += ' --output {}'.format(os.path.join(self.all_dir, 'all_stat.xls'))
        cmd_name = 'run_statistics'
        self.run_code(cmd_name, cmd, block=False)

    def run_kegg_merge_t(self):
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            self.python,
            self.kegg_merge_py,
            self.map4_r,
            self.rscript,
            self.image_magick_convert,
            self.ref_kegg_pathway_tran_xls,
            self.new_kegg_pathway_tran_xls,
            self.all_kegg_pathway_tran_xls,
            self.all_kegg_pathway_tran_dir,
            self.all_kegg_gene_tran_xls,
            self.map_html,
            self.option("kegg_version"),
            self.option("kegg_species")

        )
        cmd_name = 'run_kegg_merge_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_kegg_merge_g(self):
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            self.python,
            self.kegg_merge_py,
            self.map4_r,
            self.rscript,
            self.image_magick_convert,
            self.ref_kegg_pathway_gene_xls,
            self.new_kegg_pathway_gene_xls,
            self.all_kegg_pathway_gene_xls,
            self.all_kegg_pathway_gene_dir,
            self.all_kegg_gene_gene_xls,
            self.map_html,
            self.option("kegg_version"),
            self.option("kegg_species")
        )
        cmd_name = 'run_kegg_merge_g'
        self.run_code(cmd_name, cmd, block=False)

    def run_go_annotation_t(self):
        if os.path.isdir(self.go_txpt_dir):
            shutil.rmtree(self.go_txpt_dir)
        os.mkdir(self.go_txpt_dir)

        '''
        cmd = '{} {} {} {} {} {} {} {}'.format(
            self.python,
            self.go_annotation_py,
            self.id2terms_txpt_tsv,
            self.b2g_host,
            self.b2g_user,
            self.b2g_passwd,
            self.b2g_db,
            self.go_txpt_dir
        )
        '''

        cmd = '{} {} {} {}'.format(
            self.python,
            self.go_annotation2_py,
            self.id2terms_txpt_tsv,
            self.go_txpt_dir
        )

        cmd_name = 'run_go_annotation_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_go_annotation_g(self):
        if os.path.isdir(self.go_gene_dir):
            shutil.rmtree(self.go_gene_dir)
        os.mkdir(self.go_gene_dir)
        '''
        cmd = '{} {} {} {} {} {} {} {}'.format(
            self.python,
            self.go_annotation_py,
            self.id2terms_gene_tsv,
            self.b2g_host,
            self.b2g_user,
            self.b2g_passwd,
            self.b2g_db,
            self.go_gene_dir
        )
        '''

        cmd = '{} {} {} {}'.format(
            self.python,
            self.go_annotation2_py,
            self.id2terms_gene_tsv,
            self.go_gene_dir
        )

        cmd_name = 'run_go_annotation_g'
        self.run_code(cmd_name, cmd, block=False)

    def merge_query(self):
        cmd = '{} {} --head yes'.format(self.python, self.merge_py)
        cmd += ' --ref {}'.format(os.path.join(self.ref_dir, 'all_annot_tran.xls'))
        cmd += ' --new {}'.format(os.path.join(self.new_dir, 'all_annot_tran.xls'))
        cmd += ' --output {}'.format(os.path.join(self.all_dir, 'all_annot_tran.xls'))
        cmd_name = 'merge_query'
        self.run_code(cmd_name, cmd)

        cmd = '{} {} --head yes'.format(self.python, self.merge_query_gene_py)
        cmd += ' --ref {}'.format(os.path.join(self.ref_dir, 'all_annot_gene.xls'))
        cmd += ' --new {}'.format(os.path.join(self.new_dir, 'all_annot_gene.xls'))
        cmd += ' --output {}'.format(os.path.join(self.all_dir, 'all_annot_gene.xls'))
        cmd_name = 'merge_query_gene'
        self.run_code(cmd_name, cmd)

    def proc_args(self, args):
        ret = (
            os.path.join(self.ref_dir, args[0]), os.path.join(self.new_dir, args[1]),
            os.path.join(self.all_dir, args[2]), args[3]
        )
        return ret

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
                        self.set_error('fail to run %s, abord', variables=(n), code="33710302")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        shutil.copy(os.path.join(self.go_txpt_dir, 'level2.stat.tsv'), os.path.join(self.all_dir, 'go/go_lev2_tran.stat.xls'))
        shutil.copy(os.path.join(self.go_txpt_dir, 'level3.stat.tsv'), os.path.join(self.all_dir, 'go/go_lev3_tran.stat.xls'))
        shutil.copy(os.path.join(self.go_txpt_dir, 'level4.stat.tsv'), os.path.join(self.all_dir, 'go/go_lev4_tran.stat.xls'))
        shutil.copy(os.path.join(self.go_gene_dir, 'level2.stat.tsv'), os.path.join(self.all_dir, 'go/go_lev2_gene.stat.xls'))
        shutil.copy(os.path.join(self.go_gene_dir, 'level3.stat.tsv'), os.path.join(self.all_dir, 'go/go_lev3_gene.stat.xls'))
        shutil.copy(os.path.join(self.go_gene_dir, 'level4.stat.tsv'), os.path.join(self.all_dir, 'go/go_lev4_gene.stat.xls'))
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
            'id': 'annot_merge_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.annotation.annot_merge',
            'instant': False,
            'options': {
                'ref_class_dir': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class2',
                'new_class_dir': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/ref_annot_class_new',
                'new_mapdb_dir': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/annot_mapdb_new',
                'is_assemble': True,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
