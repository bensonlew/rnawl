# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.config import Config
from biocluster.module import Module


class AnnotationModule(Module):
    def __init__(self, work_id):
        super(AnnotationModule, self).__init__(work_id)
        options = [
            {'name': 'genome_id', 'type': 'string', 'default': None},

            {'name': 'nr_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'swissprot_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'kegg_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'cog_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'pfam_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'nr_identity', 'type': 'float', 'default': 0},
            {'name': 'nr_similarity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_identity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_similarity', 'type': 'float', 'default': 0},
            {'name': 'kegg_identity', 'type': 'float', 'default': 0},
            {'name': 'kegg_similarity', 'type': 'float', 'default': 0},
            {'name': 'cog_identity', 'type': 'float', 'default': 0},
            {'name': 'cog_similarity', 'type': 'float', 'default': 0},

            {'name': 'ref_mrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            # ['Animals', 'Plants', 'Protists', 'Fungi', 'All']
            {'name': 'kegg_database', 'type': 'string', 'default': None},
            {'name': 'new_mrna_fasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            # ['metazoa', 'viridiplantae', 'protist', 'fungi', 'nr']
            {'name': 'nr_database', 'type': 'string', 'default': None},
            {'name': 'new_mrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},

            {'name': 'is_assemble', 'type': 'bool', 'default': True},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
        collection = database['sg_genome_db']
        self.genome_doc = collection.find_one({'genome_id': self.option('genome_id')})
        self.genome_doc_anno_path_v2 = self.genome_doc['anno_path_v2']
        self.genome_doc_fa = self.genome_doc['dna_fa']
        self.genome_doc_g2t2p = self.genome_doc['g2t2p']
        self.genome_doc_bio_mart_annot = self.genome_doc['bio_mart_annot']
        self.genome_doc_biomart_gene_annotype = self.genome_doc['biomart_gene_annotype']
        self.genome_doc_ensemble2entrez = self.genome_doc['ensemble2entrez']
        self.genome_doc_kegg = self.genome_doc['kegg']
        self.db_path = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish')
        super(AnnotationModule, self).run()
        self.annot_filter_ref = self.add_module('ref_rna_v2.annot_filter')
        self.annot_class_beta_ref = self.add_module('ref_rna_v2.annot_class_beta')
        self.annot_mapdb = self.add_module('ref_rna_v2.annot_mapdb')
        self.annot_orfpfam = self.add_module('ref_rna_v2.annot_orfpfam')
        self.annot_filter_new = self.add_module('ref_rna_v2.annot_filter')
        self.annot_class_beta_new = self.add_module('ref_rna_v2.annot_class_beta')
        self.annot_merge = self.add_tool('ref_rna_v2.annotation.annot_merge')
        self.run_annot_filter_ref()
        self.run_annot_mapdb()
        self.run_annot_orfpfam()
        self.on_rely([self.annot_mapdb, self.annot_orfpfam], self.run_annot_filter_new)
        self.on_rely([self.annot_class_beta_ref, self.annot_class_beta_new], self.run_annot_merge)

    def run_annot_filter_ref(self):
        annot_path_ref = os.path.join(self.db_path, self.genome_doc_anno_path_v2)
        opts = {
            'blast_nr_xml': os.path.join(annot_path_ref, 'annot_mapdb/nr/blast.xml'),
            'blast_eggnog_xml': os.path.join(annot_path_ref, 'annot_mapdb/eggnog/blast.xml'),
            'blast_kegg_xml': os.path.join(annot_path_ref, 'annot_mapdb/kegg/blast.xml'),
            'blast_swissprot_xml': os.path.join(annot_path_ref, 'annot_mapdb/swissprot/blast.xml'),
            'pfam_domain': os.path.join(annot_path_ref, 'annot_orfpfam/pfam_domain'),
            'blast2go_annot': os.path.join(annot_path_ref, 'annot_mapdb/GO/blast2go_merge.xls'),
            'nr_evalue': self.option('nr_evalue'),
            'nr_identity': self.option('nr_identity'),
            'nr_similarity': self.option('nr_similarity'),
            'swissprot_evalue': self.option('swissprot_evalue'),
            'swissprot_identity': self.option('swissprot_identity'),
            'swissprot_similarity': self.option('swissprot_similarity'),
            'eggnog_evalue': self.option('cog_evalue'),
            'eggnog_identity': self.option('cog_identity'),
            'eggnog_similarity': self.option('cog_similarity'),
            'kegg_evalue': self.option('kegg_evalue'),
            'kegg_identity': self.option('kegg_identity'),
            'kegg_similarity': self.option('kegg_similarity'),
            'pfam_evalue': self.option('pfam_evalue')
        }
        self.annot_filter_ref.set_options(opts)
        self.annot_filter_ref.on('end', self.run_annot_class_beta_ref)
        self.annot_filter_ref.run()

    def run_annot_class_beta_ref(self):
        gtf = self.option('ref_mrna_gtf').path
        fasta = os.path.join(self.db_path, self.genome_doc_fa)
        g2t2p = os.path.join(self.db_path, self.genome_doc_g2t2p)
        des = os.path.join(self.db_path, self.genome_doc_bio_mart_annot)
        des_type = self.genome_doc_biomart_gene_annotype
        entrez = os.path.join(self.db_path, self.genome_doc_ensemble2entrez)
        blast_nr_xml = os.path.join(self.annot_filter_ref.output_dir, 'nr/blast.xml.filter.xml')
        blast_kegg_xml = os.path.join(self.annot_filter_ref.output_dir, 'kegg/blast.xml.filter.xml')
        known_ko = os.path.join(self.db_path, self.genome_doc_kegg)
        taxonomy = self.option('kegg_database')
        blast_eggnog_xml = os.path.join(self.annot_filter_ref.output_dir, 'eggnog/blast.xml.filter.xml')
        blast_swissprot_xml = os.path.join(self.annot_filter_ref.output_dir, 'swissprot/blast.xml.filter.xml')
        pfam_domain = os.path.join(self.annot_filter_ref.output_dir, 'pfam/pfam_domain.filter.xls')
        blast2go_annot = os.path.join(self.annot_filter_ref.output_dir, 'go/blast2go_merge.xls.filter.xls')
        opts = {
            'type': 'ref',
            'gtf': gtf,
            'fasta': fasta,
            'g2t2p': g2t2p,
            'des': des,
            'des_type': des_type,
            'enterz': entrez,
            'blast_nr_xml': blast_nr_xml,
            'blast_kegg_xml': blast_kegg_xml,
            'known_ko': known_ko,
            'taxonomy': taxonomy,
            'link_bgcolor': 'yellow',
            'png_bgcolor': 'FFFF00',
            'blast_eggnog_xml': blast_eggnog_xml,
            'blast_swissprot_xml': blast_swissprot_xml,
            'pfam_domain': pfam_domain,
            'blast2go_annot': blast2go_annot
        }
        self.annot_class_beta_ref.set_options(opts)
        self.annot_class_beta_ref.run()

    def run_annot_mapdb(self):
        opts = {
            'query': self.option('new_mrna_fasta').path,
            'nr_db': self.option('nr_database')
        }
        self.annot_mapdb.set_options(opts)
        self.annot_mapdb.run()

    def run_annot_orfpfam(self):
        opts = {
            'fasta': self.option('new_mrna_fasta').path,
            'gtf': self.option('new_mrna_gtf').path
        }
        self.annot_orfpfam.set_options(opts)
        self.annot_orfpfam.run()

    def run_annot_filter_new(self):
        blast_nr_xml = os.path.join(self.annot_mapdb.output_dir, 'nr/blast.xml')
        blast_eggnog_xml = os.path.join(self.annot_mapdb.output_dir, 'eggnog/blast.xml')
        blast_kegg_xml = os.path.join(self.annot_mapdb.output_dir, 'kegg/blast.xml')
        blast_swissprot_xml = os.path.join(self.annot_mapdb.output_dir, 'swissprot/blast.xml')
        pfam_domain = os.path.join(self.annot_orfpfam.output_dir, 'pfam_domain')
        blast2go_annot = os.path.join(self.annot_mapdb.output_dir, 'GO/blast2go_merge.xls')
        opts = {
            'blast_nr_xml': blast_nr_xml,
            'blast_eggnog_xml': blast_eggnog_xml,
            'blast_kegg_xml': blast_kegg_xml,
            'blast_swissprot_xml': blast_swissprot_xml,
            'pfam_domain': pfam_domain,
            'blast2go_annot': blast2go_annot,
            'nr_evalue': self.option('nr_evalue'),
            'nr_identity': self.option('nr_identity'),
            'nr_similarity': self.option('nr_similarity'),
            'swissprot_evalue': self.option('swissprot_evalue'),
            'swissprot_identity': self.option('swissprot_identity'),
            'swissprot_similarity': self.option('swissprot_similarity'),
            'eggnog_evalue': self.option('cog_evalue'),
            'eggnog_identity': self.option('cog_identity'),
            'eggnog_similarity': self.option('cog_similarity'),
            'kegg_evalue': self.option('kegg_evalue'),
            'kegg_identity': self.option('kegg_identity'),
            'kegg_similarity': self.option('kegg_similarity'),
            'pfam_evalue': self.option('pfam_evalue')
        }
        self.annot_filter_new.set_options(opts)
        self.annot_filter_new.on('end', self.run_annot_class_beta_new)
        self.annot_filter_new.run()

    def run_annot_class_beta_new(self):
        gene2trans = os.path.join(self.annot_orfpfam.output_dir, 'all_tran2gen.txt')
        des = os.path.join(self.db_path, self.genome_doc_bio_mart_annot)
        des_type = self.genome_doc_biomart_gene_annotype
        entrez = os.path.join(self.db_path, self.genome_doc_ensemble2entrez)
        blast_nr_xml = os.path.join(self.annot_filter_new.output_dir, 'nr/blast.xml.filter.xml')
        blast_kegg_xml = os.path.join(self.annot_filter_new.output_dir, 'kegg/blast.xml.filter.xml')
        taxonomy = self.option('kegg_database')
        blast_eggnog_xml = os.path.join(self.annot_filter_new.output_dir, 'eggnog/blast.xml.filter.xml')
        blast_swissprot_xml = os.path.join(self.annot_filter_new.output_dir, 'swissprot/blast.xml.filter.xml')
        pfam_domain = os.path.join(self.annot_filter_new.output_dir, 'pfam/pfam_domain.filter.xls')
        blast2go_annot = os.path.join(self.annot_filter_new.output_dir, 'go/blast2go_merge.xls.filter.xls')
        opts = {
            'type': 'new',
            'gene2trans': gene2trans,
            'des': des,
            'des_type': des_type,
            'enterz': entrez,
            'blast_nr_xml': blast_nr_xml,
            'blast_kegg_xml': blast_kegg_xml,
            'taxonomy': taxonomy,
            'link_bgcolor': 'green',
            'png_bgcolor': '00CD00',
            'blast_eggnog_xml': blast_eggnog_xml,
            'blast_swissprot_xml': blast_swissprot_xml,
            'pfam_domain': pfam_domain,
            'blast2go_annot': blast2go_annot
        }
        self.annot_class_beta_new.set_options(opts)
        self.annot_class_beta_new.run()

    def run_annot_merge(self):
        opts = {
            'ref_class_dir': self.annot_class_beta_ref.output_dir,
            'is_assemble': self.option('is_assemble'),
            'new_class_dir': self.annot_class_beta_new.output_dir,
            'new_mapdb_dir': self.annot_mapdb.output_dir
        }
        self.annot_merge.set_options(opts)
        self.annot_merge.on('end', self.set_output)
        self.annot_merge.run()

    def set_output(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.annot_merge.output_dir, self.output_dir)
        self.end()

    def end(self):
        super(AnnotationModule, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.annotation',
            'instant': False,
            'options': {
                'genome_id': 'GM0259',
                'ref_mrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20190911/WholeTranscriptome_workflow_7956_6842/output/large_gush/filter_by_express/filtered_file/known_mrna.gtf',
                'kegg_database': 'Animals',
                'new_mrna_fasta': '/mnt/ilustre/users/sanger-dev/workspace/20190911/WholeTranscriptome_workflow_7956_6842/output/large_gush/filter_by_express/filtered_file/novel_mrna.fa',
                'nr_database': 'metazoa',
                'new_mrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20190911/WholeTranscriptome_workflow_7956_6842/output/large_gush/filter_by_express/filtered_file/novel_mrna.gtf',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
