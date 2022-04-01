# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.module import Module
import os
import shutil
import unittest

class AnnotClassBetaModule(Module):
    '''
    last_modify: 2019.05.07
    '''
    def __init__(self, work_id):
        super(AnnotClassBetaModule, self).__init__(work_id)
        options = [
            # essential options for determining the process
            {'name': 'db', 'type': 'string', 'default': 'nr,swissprot,eggnog,kegg,pfam,go'},
            # ['ref', 'new']
            {'name': 'type', 'type': 'string', 'default': ''},
            # essential options for doing statistics on annotation
            {'name': 'gtf', 'type': 'infile', 'format': 'ref_rna_v2.gtf'}, # ref
            {'name': 'fasta', 'type': 'infile', 'format': 'ref_rna_v2.fasta'}, # ref
            {'name': 'g2t2p', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # ref
            {'name': 'gene2trans', 'type': 'infile', 'format': 'ref_rna_v2.common'}, # new
            {'name': 'des', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # ['type1', 'type2', 'type3']
            {'name': 'des_type', 'type': 'string', 'default': None},
            {'name': 'enterz', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # input files and related options
            {'name': 'blast_nr_xml', 'type': 'infile', 'format': 'ref_rna_v2.blast_xml'}, # nr
            {'name': 'blast_swissprot_xml', 'type': 'infile', 'format': 'ref_rna_v2.blast_xml'}, # swissprot
            {'name': 'blast_eggnog_xml', 'type': 'infile', 'format': 'ref_rna_v2.blast_xml'}, # eggnog
            {'name': 'blast_kegg_xml', 'type': 'infile', 'format': 'ref_rna_v2.blast_xml'}, # kegg
            {'name': 'known_ko', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # ['All', 'Animals', 'Plants', 'Protists']
            {'name': 'taxonomy', 'type': 'string', 'default': ''},
            # {'ref': 'yellow', 'new': 'green'}
            {'name': 'link_bgcolor', 'type': 'string', 'default': ''},
            # {'ref': 'FFFF00', 'new': '00CD00'}
            {'name': 'png_bgcolor', 'type': 'string', 'default': ''},
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # go
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # pfam
            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "2019"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "string_version", "type": "string", "default": ""},
            {"name": "cog_version", "type": "string", "default": "2019"},
            {"name": "kegg_version", "type": "string", "default": "2019"},
            {"name": "kegg_subtax1", "type": "string", "default": None},
            {"name": "kegg_subtax2", "type": "string", "default": None},
            {"name": "kegg_species", "type": "string", "default": None},
            {"name": "pir_version", "type": "string", "default":"2019"},
            {"name": "ncbi_taxonmy_version", "type": "string", "default":"2019"},
            {"name": "go_version", "type": "string", "default": "2019"},
            {"name": "pfam_version", "type": "string", "default": "32"},
            {"name": "merge_type", "type": "string", "default": "partial"},

        ]
        self.add_option(options)
        self.tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data']:
            event['data']['start'].start()
        if 'end' in event['data']:
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(AnnotClassBetaModule, self).run()
        self.run_annot_begin()

    def run_annot_begin(self):
        self.step.add_steps('annot_begin')
        self.annot_begin = self.add_tool('ref_rna_v2.annotation.annot_begin')
        if self.option('type') == 'ref':
            options = {
                'gtf': self.option('gtf'),
                'fasta': self.option('fasta'),
                'g2t2p': self.option('g2t2p')
            }
        elif self.option('type') == 'new':
            options = {
                'gene2trans': self.option('gene2trans')
            }
        self.annot_begin.set_options(options)
        self.annot_begin.on('start', self.set_step, {'start': self.step.annot_begin})
        self.annot_begin.on('end', self.set_step, {'end': self.step.annot_begin})
        self.annot_begin.on('end', self.set_output, 'annot_begin')
        self.annot_begin.on('end', self.run_annot_main)
        self.annot_begin.run()

    def run_annot_main(self):
        self.run_nr_annot()
        self.run_swissprot_annot()
        self.run_eggnog_annot()
        self.run_kegg_annot()
        self.run_go_annot()
        self.run_pfam_annot()
        self.on_rely(self.tools, self.run_annot_final)

    def run_nr_annot(self):
        self.step.add_steps('nr_annot')
        self.nr_annot = self.add_tool('denovo_rna_v2.annotation.nr_annot')
        options = {
            'blast_xml': self.option('blast_nr_xml'),
            'longest_t2g': self.annot_begin.option('longest_t2g'),
            'ncbi_taxonmy_version': self.option("ncbi_taxonmy_version")
        }
        self.nr_annot.set_options(options)
        self.nr_annot.on('start', self.set_step, {'start': self.step.nr_annot})
        self.nr_annot.on('end', self.set_step, {'end': self.step.nr_annot})
        self.nr_annot.on('end', self.set_output, 'nr_annot')
        self.nr_annot.run()
        self.tools.append(self.nr_annot)

    def run_swissprot_annot(self):
        self.step.add_steps('swissprot_annot')
        self.swissprot_annot = self.add_tool('denovo_rna_v2.annotation.swissprot_annot')
        options = {
            'blast_xml': self.option('blast_swissprot_xml'),
            'longest_t2g': self.annot_begin.option('longest_t2g')
        }
        self.swissprot_annot.set_options(options)
        self.swissprot_annot.on('start', self.set_step, {'start': self.step.swissprot_annot})
        self.swissprot_annot.on('end', self.set_step, {'end': self.step.swissprot_annot})
        self.swissprot_annot.on('end', self.set_output, 'swissprot_annot')
        self.swissprot_annot.run()
        self.tools.append(self.swissprot_annot)

    def run_eggnog_annot(self):
        self.step.add_steps('eggnog_annot')
        self.eggnog_annot = self.add_tool('ref_rna_v2.annotation.eggnog_annot')
        options = {
            'blast_xml': self.option('blast_eggnog_xml'),
            'eggnog_version': self.option('eggnog_version'),
            'longest_t2g': self.annot_begin.option('longest_t2g'),
        }
        self.eggnog_annot.set_options(options)
        self.eggnog_annot.on('start', self.set_step, {'start': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_step, {'end': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_output, 'eggnog_annot')
        self.eggnog_annot.run()
        self.tools.append(self.eggnog_annot)

    def run_kegg_annot(self):
        self.step.add_steps('kegg_annot')
        self.kegg_annot = self.add_tool('ref_rna_v2.annotation.kegg_annot')
        options = {
            'blast_xml': self.option('blast_kegg_xml'),
            'longest_t2g': self.annot_begin.option('longest_t2g'),
            'known_ko': self.option('known_ko'),
            'taxonomy': self.option('taxonomy'),
            'link_bgcolor': self.option('link_bgcolor'),
            'png_bgcolor': self.option('png_bgcolor'),
            'kegg_version': self.option('kegg_version'),
            'kegg_subtax1': self.option('kegg_subtax1'),
            'kegg_subtax2': self.option('kegg_subtax2'),
            'kegg_species': self.option('kegg_species')
        }
        self.kegg_annot.set_options(options)
        self.kegg_annot.on('start', self.set_step, {'start': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_step, {'end': self.step.kegg_annot})
        self.kegg_annot.on('end', self.set_output, 'kegg_annot')
        self.kegg_annot.run()
        self.tools.append(self.kegg_annot)

    def run_go_annot(self):
        self.step.add_steps('go_annot')
        self.go_annot = self.add_tool('ref_rna_v2.annotation.go_annot')
        options = {
            'blast2go_annot': self.option('blast2go_annot'),
            'longest_t2g': self.annot_begin.option('longest_t2g'),
            'go_version': self.option('go_version')
        }
        self.go_annot.set_options(options)
        self.go_annot.on('start', self.set_step, {'start': self.step.go_annot})
        self.go_annot.on('end', self.set_step, {'end': self.step.go_annot})
        self.go_annot.on('end', self.set_output, 'go_annot')
        self.go_annot.run()
        self.tools.append(self.go_annot)

    def run_pfam_annot(self):
        self.step.add_steps('pfam_annot')
        self.pfam_annot = self.add_tool('ref_rna_v2.annotation.pfam_annot')
        options = {
            'pfam_domain': self.option('pfam_domain'),
            'longest_t2g': self.annot_begin.option('longest_t2g')
        }
        self.pfam_annot.set_options(options)
        self.pfam_annot.on('start', self.set_step, {'start': self.step.pfam_annot})
        self.pfam_annot.on('end', self.set_step, {'end': self.step.pfam_annot})
        self.pfam_annot.on('end', self.set_output, 'pfam_annot')
        self.pfam_annot.run()
        self.tools.append(self.pfam_annot)

    def run_annot_final(self):
        lines = list()
        for tool in self.tools:
            lines.extend([
                '{}\t{}\t{}\n'.format(tool.option('txpt_list').path, tool.option('database'), 'transcript'),
                '{}\t{}\t{}\n'.format(tool.option('gene_list').path, tool.option('database'), 'gene')
            ])
        else:
            loc2db2type = os.path.join(self.work_dir, 'loc2db2type.tsv')
            open(loc2db2type, 'w').writelines(lines)
        self.step.add_steps('annot_final')
        self.annot_final = self.add_tool('denovo_rna_v2.annotation.annot_final')
        options = {
            't2g2r2l2p': self.annot_begin.option('t2g2r2l2p'),
            'loc2db2type': loc2db2type,
            'nr_table': self.nr_annot.option('txpt_table'),
            'swissprot_table': self.swissprot_annot.option('txpt_table'),
            'cog_table': self.eggnog_annot.option('cog_table'),
            'kegg_table': self.kegg_annot.option('txpt_kegg_table'),
            'kegg_version': self.option('kegg_version'),
            'pfam_domain': self.pfam_annot.option('txpt_pfam_domain'),
            'go_list': self.go_annot.option('txpt_id2terms'),

        }
        self.annot_final.set_options(options)
        self.annot_final.on('start', self.set_step, {'start': self.step.annot_final})
        self.annot_final.on('end', self.set_step, {'end': self.step.annot_final})
        self.annot_final.on('end', self.set_output, 'annot_final')
        self.annot_final.run()

    def set_output(self, event):
        if event['data'] == 'annot_begin':
            shutil.copy(event['bind_object'].option('t2g2r2l2p').path, os.path.join(self.output_dir, 'all_tran2gene.txt'))
        if event['data'] == 'nr_annot':
            nr_dir = os.path.join(self.output_dir, 'nr')
            if os.path.isdir(nr_dir):
                shutil.rmtree(nr_dir)
            shutil.copytree(event['bind_object'].output_dir, nr_dir)
            os.rename(os.path.join(nr_dir, 'nr.G.tsv'), os.path.join(nr_dir, 'nr_blast_gene.xls'))
            os.rename(os.path.join(nr_dir, 'nr.T.tsv'), os.path.join(nr_dir, 'nr_blast_tran.xls'))
            os.rename(os.path.join(nr_dir, 'nr.G.list'), os.path.join(nr_dir, 'nr_venn_gene.txt'))
            os.rename(os.path.join(nr_dir, 'nr.T.list'), os.path.join(nr_dir, 'nr_venn_tran.txt'))
        if event['data'] == 'swissprot_annot':
            swissprot_dir = os.path.join(self.output_dir, 'swissprot')
            if os.path.isdir(swissprot_dir):
                shutil.rmtree(swissprot_dir)
            shutil.copytree(event['bind_object'].output_dir, swissprot_dir)
            os.rename(os.path.join(swissprot_dir, 'swissprot.G.tsv'), os.path.join(swissprot_dir, 'swissprot_blast_gene.xls'))
            os.rename(os.path.join(swissprot_dir, 'swissprot.T.tsv'), os.path.join(swissprot_dir, 'swissprot_blast_tran.xls'))
            os.rename(os.path.join(swissprot_dir, 'swissprot.G.list'), os.path.join(swissprot_dir, 'swissprot_venn_gene.txt'))
            os.rename(os.path.join(swissprot_dir, 'swissprot.T.list'), os.path.join(swissprot_dir, 'swissprot_venn_tran.txt'))
        if event['data'] == 'eggnog_annot':
            cog_dir = os.path.join(self.output_dir, 'cog')
            if os.path.isdir(cog_dir):
                shutil.rmtree(cog_dir)
            shutil.copytree(event['bind_object'].output_dir, cog_dir)
            os.rename(os.path.join(cog_dir, 'cog.tsv'), os.path.join(cog_dir, 'cog_list_tran.xls'))
            os.rename(os.path.join(cog_dir, 'cog.G.list'), os.path.join(cog_dir, 'cog_venn_gene.txt'))
            os.rename(os.path.join(cog_dir, 'cog.T.list'), os.path.join(cog_dir, 'cog_venn_tran.txt'))
        if event['data'] == 'kegg_annot':
            kegg_dir = os.path.join(self.output_dir, 'kegg')
            if os.path.isdir(kegg_dir):
                shutil.rmtree(kegg_dir)
            shutil.copytree(event['bind_object'].output_dir, kegg_dir)
            os.rename(os.path.join(kegg_dir, 'G/kegg_layer.tsv'), os.path.join(kegg_dir, 'kegg_layer_gene.xls'))
            os.rename(os.path.join(kegg_dir, 'T/kegg_layer.tsv'), os.path.join(kegg_dir, 'kegg_layer_tran.xls'))
            os.rename(os.path.join(kegg_dir, 'G/kegg_table.tsv'), os.path.join(kegg_dir, 'kegg_gene_gene.xls'))
            os.rename(os.path.join(kegg_dir, 'T/kegg_table.tsv'), os.path.join(kegg_dir, 'kegg_gene_tran.xls'))
            os.rename(os.path.join(kegg_dir, 'G/pathways'), os.path.join(kegg_dir, 'kegg_pathway_gene_dir'))
            os.rename(os.path.join(kegg_dir, 'G/pathway_table.tsv'), os.path.join(kegg_dir, 'kegg_pathway_gene.xls'))
            os.rename(os.path.join(kegg_dir, 'T/pathways'), os.path.join(kegg_dir, 'kegg_pathway_tran_dir'))
            os.rename(os.path.join(kegg_dir, 'T/pathway_table.tsv'), os.path.join(kegg_dir, 'kegg_pathway_tran.xls'))
            os.rename(os.path.join(kegg_dir, 'kegg.G.list'), os.path.join(kegg_dir, 'kegg_venn_gene.txt'))
            os.rename(os.path.join(kegg_dir, 'kegg.T.list'), os.path.join(kegg_dir, 'kegg_venn_tran.txt'))
            shutil.rmtree(os.path.join(kegg_dir, 'G'))
            shutil.rmtree(os.path.join(kegg_dir, 'T'))
        if event['data'] == 'go_annot':
            go_dir = os.path.join(self.output_dir, 'go')
            if os.path.isdir(go_dir):
                shutil.rmtree(go_dir)
            shutil.copytree(event['bind_object'].output_dir, go_dir)
            os.rename(os.path.join(go_dir, 'G/level2.stat.tsv'), os.path.join(go_dir, 'go_lev2_gene.stat.xls'))
            os.rename(os.path.join(go_dir, 'G/level3.stat.tsv'), os.path.join(go_dir, 'go_lev3_gene.stat.xls'))
            os.rename(os.path.join(go_dir, 'G/level4.stat.tsv'), os.path.join(go_dir, 'go_lev4_gene.stat.xls'))
            os.rename(os.path.join(go_dir, 'T/level2.stat.tsv'), os.path.join(go_dir, 'go_lev2_tran.stat.xls'))
            os.rename(os.path.join(go_dir, 'T/level3.stat.tsv'), os.path.join(go_dir, 'go_lev3_tran.stat.xls'))
            os.rename(os.path.join(go_dir, 'T/level4.stat.tsv'), os.path.join(go_dir, 'go_lev4_tran.stat.xls'))
            os.rename(os.path.join(go_dir, 'id2terms.G.tsv'), os.path.join(go_dir, 'go_list_gene.xls'))
            os.rename(os.path.join(go_dir, 'id2terms.T.tsv'), os.path.join(go_dir, 'go_list_tran.xls'))
            os.rename(os.path.join(go_dir, 'go.G.list'), os.path.join(go_dir, 'go_venn_gene.txt'))
            os.rename(os.path.join(go_dir, 'go.T.list'), os.path.join(go_dir, 'go_venn_tran.txt'))
            shutil.rmtree(os.path.join(go_dir, 'G'))
            shutil.rmtree(os.path.join(go_dir, 'T'))
        if event['data'] == 'pfam_annot':
            pfam_dir = os.path.join(self.output_dir, 'pfam')
            if os.path.isdir(pfam_dir):
                shutil.rmtree(pfam_dir)
            shutil.copytree(event['bind_object'].output_dir, pfam_dir)
            os.rename(os.path.join(pfam_dir, 'pfam.G.tsv'), os.path.join(pfam_dir, 'pfam_domain_gene.xls'))
            os.rename(os.path.join(pfam_dir, 'pfam.T.tsv'), os.path.join(pfam_dir, 'pfam_domain_tran.xls'))
            os.rename(os.path.join(pfam_dir, 'pfam.G.list'), os.path.join(pfam_dir, 'pfam_venn_gene.txt'))
            os.rename(os.path.join(pfam_dir, 'pfam.T.list'), os.path.join(pfam_dir, 'pfam_venn_tran.txt'))
        if event['data'] == 'annot_final':
            shutil.copy(event['bind_object'].option('query').path, os.path.join(self.output_dir, 'all_annot.xls'))
            shutil.copy(event['bind_object'].option('statistics').path, os.path.join(self.output_dir, 'all_stat.xls'))
            self.end()

    def end(self):
        super(AnnotClassBetaModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run script to do test.
    '''
    def test_ref(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_class_beta_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'ref_rna_v2.annot_class_beta',
            'instant': False,
            'options': {
                'db': 'nr,swissprot,kegg,eggnog,pfam,go',
                'type': 'ref',
                'gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'fasta': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
                'g2t2p': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/Annotation_v2/g2t2p',
                'des': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                'des_type': 'type1',
                'enterz': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/NCBI/Homo_sapiens.GRCh38.biomart_enterz.txt',
                'blast_nr_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterRef/nr/blast.xml.filter.xml',
                'blast_swissprot_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterRef/swissprot/blast.xml.filter.xml',
                'blast_eggnog_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterRef/eggnog/blast.xml.filter.xml',
                'blast_kegg_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterRef/kegg/blast.xml.filter.xml',
                'known_ko': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/KEGG/Homo_sapiens.GRCh38.pathway',
                'taxonomy': 'Animals',
                'link_bgcolor': 'yellow',
                'png_bgcolor': '#FFFF00',
                'blast2go_annot': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterRef/go/blast2go_merge.xls.filter.xls',
                'pfam_domain': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterRef/pfam/pfam_domain.filter.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_new(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annot_class_beta_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'ref_rna_v2.annot_class_beta',
            'instant': False,
            'options': {
                'db': 'nr,swissprot,kegg,eggnog,pfam,go',
                'type': 'new',
                'gene2trans': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotOrfpfamNew/all_tran2gen.txt',
                'des': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                'des_type': 'type1',
                'enterz': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/NCBI/Homo_sapiens.GRCh38.biomart_enterz.txt',
                'blast_nr_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterNew/nr/blast.xml.filter.xml',
                'blast_swissprot_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterNew/swissprot/blast.xml.filter.xml',
                'blast_eggnog_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterNew/eggnog/blast.xml.filter.xml',
                'blast_kegg_xml': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterNew/kegg/blast.xml.filter.xml',
                'taxonomy': 'Animals',
                'link_bgcolor': 'green',
                'png_bgcolor': '00CD00',
                'blast2go_annot': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterNew/go/blast2go_merge.xls.filter.xls',
                'pfam_domain': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/data/AnnotFilterNew/pfam/pfam_domain.filter.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ref')])
    unittest.TextTestRunner(verbosity=2).run(suite)
