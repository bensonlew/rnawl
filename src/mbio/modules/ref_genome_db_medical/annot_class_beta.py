# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

import os
import shutil
import unittest

from biocluster.module import Module
from mbio.packages.lnc_rna.copy_file import CopyFile


class AnnotClassBetaModule(Module):
    '''
    last_modify: 2019.05.07
    '''

    def __init__(self, work_id):
        super(AnnotClassBetaModule, self).__init__(work_id)
        options = [
            # essential options for determining the process
            {'name': 'db', 'type': 'string', 'default': 'nr,uniprot,eggnog,kegg,pfam,go'},
            # ['ref', 'new']
            {'name': 'type', 'type': 'string', 'default': ''},
            # essential options for doing statistics on annotation
            {'name': 'gtf', 'type': 'infile', 'format': 'ref_rna_v2.gtf'},  # ref
            {'name': 'fasta', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},  # ref
            {'name': 'g2t2p', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # ref
            {'name': 'gene2trans', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # new
            {'name': 'des', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # ['type1', 'type2', 'type3']
            {'name': 'des_type', 'type': 'string', 'default': None},
            {'name': 'enterz', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # input files and related options
            {'name': 'nr_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # nr
            {'name': 'uniprot_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # uniprot
            {'name': 'eggnog_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # eggnog
            {'name': 'kegg_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'go_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'do_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'reactome_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'disgenet_ids', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # ['All', 'Animals', 'Plants', 'Protists']
            {'name': 'taxonomy', 'type': 'string', 'default': ''},
            # {'ref': 'yellow', 'new': 'green'}
            {'name': 'link_bgcolor', 'type': 'string', 'default': ''},
            # {'ref': 'FFFF00', 'new': '00CD00'}
            {'name': 'png_bgcolor', 'type': 'string', 'default': ''},
            {'name': 'blast2go_annot', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # go
            {'name': 'pfam_domain', 'type': 'infile', 'format': 'ref_rna_v2.common'},  # pfam
            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "uniprot_version", "type": "string", "default": "2019"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "string_version", "type": "string", "default": "2019"},
            {"name": "cog_version", "type": "string", "default": "2019"},
            {"name": "kegg_version", "type": "string", "default": "202007"},
            {"name": "kegg_subtax1", "type": "string", "default": None},
            {"name": "kegg_subtax2", "type": "string", "default": None},
            {"name": "kegg_species", "type": "string", "default": None},
            {"name": "pir_version", "type": "string", "default":"2019"},
            {"name": "ncbi_taxonmy_version", "type": "string", "default":"2019"},
            {"name": "go_version", "type": "string", "default": "2019"},
            {"name": "do_version", "type": "string", "default": "202008"},
            {"name": "reactome_version", "type": "string", "default": "202007"},
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
        self.annot_begin = self.add_tool('ref_genome_db_medical.annotation.annot_begin')
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
        self.run_uniprot_annot()
        self.run_eggnog_annot()
        self.run_kegg_annot()
        self.run_go_annot()
        self.run_pfam_annot()
        self.run_do_annot()
        self.run_reactome_annot()
        self.run_disgenet_annot()
        self.on_rely(self.tools, self.run_annot_final)

    def run_nr_annot(self):
        self.step.add_steps('nr_annot')
        self.nr_annot = self.add_tool('ref_genome_db_medical.annotation.nr_annot')
        options = {
            'nr_ids': self.option('nr_ids'),
            'longest_t2g': self.annot_begin.option('longest_t2g')
        }
        self.nr_annot.set_options(options)
        self.nr_annot.on('start', self.set_step, {'start': self.step.nr_annot})
        self.nr_annot.on('end', self.set_step, {'end': self.step.nr_annot})
        self.nr_annot.on('end', self.set_output, 'nr_annot')
        self.nr_annot.run()
        self.tools.append(self.nr_annot)

    def run_uniprot_annot(self):
        self.step.add_steps('uniprot_annot')
        self.uniprot_annot = self.add_tool('ref_genome_db_medical.annotation.uniprot_annot')
        options = {
            'uniprot_ids': self.option('uniprot_ids'),
            'longest_t2g': self.annot_begin.option('longest_t2g')
        }
        self.uniprot_annot.set_options(options)
        self.uniprot_annot.on('start', self.set_step, {'start': self.step.uniprot_annot})
        self.uniprot_annot.on('end', self.set_step, {'end': self.step.uniprot_annot})
        self.uniprot_annot.on('end', self.set_output, 'uniprot_annot')
        self.uniprot_annot.run()
        self.tools.append(self.uniprot_annot)

    def run_eggnog_annot(self):
        self.step.add_steps('eggnog_annot')
        self.eggnog_annot = self.add_tool('ref_genome_db_medical.annotation.eggnog_annot')
        options = {
            'eggnog_ids': self.option('eggnog_ids'),
            'longest_t2g': self.annot_begin.option('longest_t2g'),
            'eggnog_version': self.option('eggnog_version')
        }
        self.eggnog_annot.set_options(options)
        self.eggnog_annot.on('start', self.set_step, {'start': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_step, {'end': self.step.eggnog_annot})
        self.eggnog_annot.on('end', self.set_output, 'eggnog_annot')
        self.eggnog_annot.run()
        self.tools.append(self.eggnog_annot)

    def run_kegg_annot(self):
        self.step.add_steps('kegg_annot')
        self.kegg_annot = self.add_tool('ref_genome_db_medical.annotation.kegg_annot')
        options = {
            'longest_t2g': self.annot_begin.option('longest_t2g'),
            'known_ko': self.option('kegg_ids'),
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
        self.go_annot = self.add_tool('ref_genome_db_medical.annotation.go_annot')
        options = {
            'go_ids': self.option('go_ids'),
            'longest_t2g': self.annot_begin.option('longest_t2g'),
            'go_version': self.option('go_version')
        }
        self.go_annot.set_options(options)
        self.go_annot.on('start', self.set_step, {'start': self.step.go_annot})
        self.go_annot.on('end', self.set_step, {'end': self.step.go_annot})
        self.go_annot.on('end', self.set_output, 'go_annot')
        self.go_annot.run()
        self.tools.append(self.go_annot)

    def run_do_annot(self):
        self.step.add_steps('do_annot')
        self.do_annot = self.add_tool('ref_genome_db_medical.annotation.do_annot')
        options = {
            'do_ids': self.option('do_ids'),
            'longest_t2g': self.annot_begin.option('longest_t2g')
        }
        self.do_annot.set_options(options)
        self.do_annot.on('start', self.set_step, {'start': self.step.do_annot})
        self.do_annot.on('end', self.set_step, {'end': self.step.do_annot})
        self.do_annot.on('end', self.set_output, 'do_annot')
        self.do_annot.run()
        self.tools.append(self.do_annot)

    def run_reactome_annot(self):
        self.step.add_steps('reactome_annot')
        self.reactome_annot = self.add_tool('ref_genome_db_medical.annotation.reactome_annot')
        options = {
            'reactome_ids': self.option('reactome_ids'),
            'longest_t2g': self.annot_begin.option('longest_t2g'),
            'reactome_version': self.option('reactome_version')
        }
        self.reactome_annot.set_options(options)
        self.reactome_annot.on('start', self.set_step, {'start': self.step.reactome_annot})
        self.reactome_annot.on('end', self.set_step, {'end': self.step.reactome_annot})
        self.reactome_annot.on('end', self.set_output, 'reactome_annot')
        self.reactome_annot.run()
        self.tools.append(self.reactome_annot)

    def run_disgenet_annot(self):
        self.step.add_steps('disgenet_annot')
        self.disgenet_annot = self.add_tool('ref_genome_db_medical.annotation.disgenet_annot')
        options = {
            'disgenet_ids': self.option('disgenet_ids'),
            'longest_t2g': self.annot_begin.option('longest_t2g')
        }
        self.disgenet_annot.set_options(options)
        self.disgenet_annot.on('start', self.set_step, {'start': self.step.disgenet_annot})
        self.disgenet_annot.on('end', self.set_step, {'end': self.step.disgenet_annot})
        self.disgenet_annot.on('end', self.set_output, 'disgenet_annot')
        self.disgenet_annot.run()
        self.tools.append(self.disgenet_annot)

    def run_pfam_annot(self):
        self.step.add_steps('pfam_annot')
        self.pfam_annot = self.add_tool('ref_genome_db_medical.annotation.pfam_annot')
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
        self.annot_final = self.add_tool('ref_genome_db_medical.annotation.annot_final')
        options = {
            't2g2r2l2p': self.annot_begin.option('t2g2r2l2p'),
            'loc2db2type': loc2db2type,
            'nr_table': self.nr_annot.option('txpt_table').path,
            'uniprot_table': self.uniprot_annot.option('txpt_table').path,
            'cog_table': self.eggnog_annot.option('cog_table').path,
            'kegg_table': self.kegg_annot.option('txpt_kegg_table').path,
            'kegg_table_spe': self.kegg_annot.output_dir + '/spe_T/spe_kegg_table.tsv',
            'pfam_domain': self.pfam_annot.option('txpt_pfam_domain').path,
            'go_list': self.go_annot.option('txpt_id2terms').path,
            'des': self.option('des'),
            'des_type': self.option('des_type'),
            'kegg_version': self.option('kegg_version'),
            'enterz': self.option('enterz'),
            'nr_table_gene': self.nr_annot.option('gene_table').path,
            'uniprot_table_gene': self.uniprot_annot.option('gene_table').path,
            'cog_table_gene': self.eggnog_annot.option('cog_gene_table').path,
            'kegg_table_gene': self.kegg_annot.option('gene_kegg_table').path,
            'kegg_table_gene_spe': self.kegg_annot.output_dir + '/spe_G/spe_kegg_table.tsv',
            'pfam_domain_gene': self.pfam_annot.option('gene_pfam_domain').path,
            'go_list_gene': self.go_annot.option('gene_id2terms').path,
            'do_gene': self.do_annot.option('gene_id2terms').path,
            'reactome_gene': self.reactome_annot.option('reactome_table').path,
            'disgenet_gene': self.disgenet_annot.option('disgenet_table').path,
        }
        self.annot_final.set_options(options)
        self.annot_final.on('start', self.set_step, {'start': self.step.annot_final})
        self.annot_final.on('end', self.set_step, {'end': self.step.annot_final})
        self.annot_final.on('end', self.set_output, 'annot_final')
        self.annot_final.run()

    def set_output(self, event):
        if event['data'] == 'annot_begin':
            shutil.copy(event['bind_object'].option('t2g2r2l2p').path,
                        os.path.join(self.output_dir, 'all_tran2gene.txt'))
        if event['data'] == 'nr_annot':
            nr_dir = os.path.join(self.output_dir, 'nr')
            if os.path.isdir(nr_dir):
                shutil.rmtree(nr_dir)
            shutil.copytree(event['bind_object'].output_dir, nr_dir)
            os.link(event['bind_object'].option("gene_table").path, os.path.join(nr_dir, 'nr_gene.xls'))
            os.link(event['bind_object'].option("txpt_table").path, os.path.join(nr_dir, 'nr_tran.xls'))
            # os.link(os.path.join(nr_dir, 'nr.T.ids'), os.path.join(nr_dir, 'nr_tran.xls'))
            os.link(os.path.join(nr_dir, 'nr.G.list'), os.path.join(nr_dir, 'nr_venn_gene.txt'))
            os.link(os.path.join(nr_dir, 'nr.T.list'), os.path.join(nr_dir, 'nr_venn_tran.txt'))
        if event['data'] == 'uniprot_annot':
            uniprot_dir = os.path.join(self.output_dir, 'uniprot')
            if os.path.isdir(uniprot_dir):
                shutil.rmtree(uniprot_dir)
            shutil.copytree(event['bind_object'].output_dir, uniprot_dir)
            os.link(event['bind_object'].option("gene_table").path,
                      os.path.join(uniprot_dir, 'uniprot_gene.xls'))
            os.link(event['bind_object'].option("txpt_table").path,
                      os.path.join(uniprot_dir, 'uniprot_tran.xls'))
            os.rename(os.path.join(uniprot_dir, 'uniprot.G.list'),
                      os.path.join(uniprot_dir, 'uniprot_venn_gene.txt'))
            os.rename(os.path.join(uniprot_dir, 'uniprot.T.list'),
                      os.path.join(uniprot_dir, 'uniprot_venn_tran.txt'))

        if event['data'] == 'eggnog_annot':
            cog_dir = os.path.join(self.output_dir, 'cog')
            if os.path.isdir(cog_dir):
                shutil.rmtree(cog_dir)
            shutil.copytree(event['bind_object'].output_dir, cog_dir)
            os.rename(os.path.join(cog_dir, 'cog.T.tsv'), os.path.join(cog_dir, 'cog_list_tran.xls'))
            os.rename(os.path.join(cog_dir, 'cog.G.tsv'), os.path.join(cog_dir, 'cog_list_gene.xls'))
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
        if event['data'] == 'reactome_annot':
            reactome_dir = os.path.join(self.output_dir, 'reactome')
            if os.path.isdir(reactome_dir):
                shutil.rmtree(reactome_dir)
            shutil.copytree(event['bind_object'].output_dir, reactome_dir)
        if event['data'] == 'do_annot':
            do_dir = os.path.join(self.output_dir, 'do')
            if os.path.isdir(do_dir):
                shutil.rmtree(do_dir)
            shutil.copytree(event['bind_object'].output_dir, do_dir)
        if event['data'] == 'disgenet_annot':
            disgenet_dir = os.path.join(self.output_dir, 'disgenet')
            if os.path.isdir(disgenet_dir):
                shutil.rmtree(disgenet_dir)
            shutil.copytree(event['bind_object'].output_dir, disgenet_dir)
        if event['data'] == 'annot_final':
            shutil.copy(event['bind_object'].option('query').path, os.path.join(self.output_dir, 'all_annot_tran.xls'))
            shutil.copy(event['bind_object'].option('query_gene').path, os.path.join(self.output_dir, 'all_annot_gene.xls'))
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
        import json
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        map_db_dir = "/mnt/lustre/users/sanger-dev/sg-users/liubinxu/medical_transcriptome/annot_mapdb2"
        data = {
            'id': 'annot_class_beta_db_rerun2',
            'type': 'module',
            "rerun": True,
            "WPM": True,
            "wfm_port": "7321",
            "project_sn": "test",
            "client": "client03",
            "CLUSTER": "sanger-dev",

            'name': 'ref_genome_db_medical.annot_class_beta',
            'instant': False,
            "skip_all_success": True,
            'options': {
                'db': 'nr,uniprot,kegg,eggnog,pfam,go',
                'type': 'ref',
                'gtf': '/mnt/lustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/gtf/Homo_sapiens.GRCh38.98.gtf',
                'fasta': '/mnt/lustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa',
                'g2t2p': '/mnt/lustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/Annotation_v2/g2t2p.txt',
                'des': '/mnt/lustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/biomart/biomart1.txt',
                'des_type': 'type1',
                'enterz': '/mnt/lustre/users/sanger-dev/sg-users/liubinxu/genome_annotation/ids.tsv',
                'nr_ids': map_db_dir + "/" + 'nr/nr_annot.tsv',
                'uniprot_ids': map_db_dir + "/" + 'uniprot/uniprot_annot.tsv',
                'eggnog_ids': map_db_dir + "/" + 'eggnog/eggnog_annot.tsv',
                'kegg_ids': map_db_dir + "/" + 'kegg/kegg_annot.tsv',
                'taxonomy': 'Animals',
                'kegg_species': 'hsa',
                'link_bgcolor': 'yellow',
                'png_bgcolor': 'FFFF00',
                'go_ids': map_db_dir + "/" + 'go/all2go_annot.xls',
                'do_ids': map_db_dir + "/" + 'do/do_annot.tsv',
                'reactome_ids': map_db_dir + "/" + 'reactome/reactome_annot.tsv',
                'disgenet_ids': map_db_dir + "/" + 'disgenet/disgenet_annot.tsv',
                'pfam_domain': map_db_dir + "/" + 'pfam/pfam_annot.tsv',
            }
        }
        with open("data.json", "w") as f:
            f.write(json.dumps(data, indent=4))

        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.is_skip=True
        # wf.run()




if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ref')])
    unittest.TextTestRunner(verbosity=2).run(suite)
