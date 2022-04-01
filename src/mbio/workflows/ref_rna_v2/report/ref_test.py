# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.workflow import Workflow
from biocluster.config import Config
import os
import shutil
import unittest

class RefTestWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RefTestWorkflow, self).__init__(wsheet_object)
        options = [
            # basic
            {'name': 'sample_num', 'type': 'string', 'default': None}, # ['multiple', 'single']
            {'name': 'fq_type', 'type': 'string', 'default': None}, # ['PE', 'SE']
            {'name': 'quality_score_system', 'type': 'string', 'default': None},  # ['phred+33', 'phred+64']
            {'name': 'strand_specific', 'type': 'bool', 'default': False}, # [False, True]
            {'name': 'strand_dir', 'type': 'string', 'default': None}, # {'forward': 'RF', 'reverse': 'FR'}
            {'name': 'is_duplicate', 'type': 'bool', 'default': True}, # [True, False]
            {'name': 'ref_genome', 'type': 'string', 'default': None}, # sg_genome_db.name
            {'name': 'genome_version', 'type': 'string', 'default': None},  # sg_genome_db.assembly
            {'name': 'genome_annot_version', 'type': 'string', 'default': None}, # sg_genome_db.genome_annot_version
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'control_file', 'type': 'infile', 'format': 'sample.control_table'},
            # pretreatment
            {'name': 'level', 'type': 'string', 'default': None}, # ['transcript', 'gene']
            {'name': 'align_method', 'type': 'string', 'default': 'hisat'}, # ['hisat', 'tophat']
            {'name': 'map_assess_method', 'type': 'string', 'default': 'saturation,distribution,coverage,chr_stat'},
            {'name': 'is_assemble', 'type': 'bool', 'default': True}, # [True, False]
            {'name': 'assemble_method', 'type': 'string', 'default': None}, # ['stringtie', 'cufflinks']
            # annotation
            {'name': 'nr_database', 'type': 'string', 'default': None},  # ['Animal', 'Plant', 'Protist', 'Fungi', 'All']
            {'name': 'kegg_database', 'type': 'string', 'default': None}, # ['Animal', 'Plant', 'Protist', 'All']
            {'name': 'nr_blast_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'nr_blast_identity', 'type': 'float', 'default': 0},
            {'name': 'nr_blast_similarity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_blast_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'swissprot_blast_identity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_blast_similarity', 'type': 'float', 'default': 0},
            {'name': 'cog_blast_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'cog_blast_identity', 'type': 'float', 'default': 0},
            {'name': 'cog_blast_similarity', 'type': 'float', 'default': 0},
            {'name': 'kegg_blast_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'kegg_blast_identity', 'type': 'float', 'default': 0},
            {'name': 'kegg_blast_similarity', 'type': 'float', 'default': 0},
            {'name': 'pfam_blast_evalue', 'type': 'float', 'default': 1e-5},
        ]
        
        self.task_id = self.sheet.id
        self.project_sn = self.sheet.project_sn
        self.add_option(options)
        self.set_options(self._sheet.options())

        db_collection = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]['sg_genome_db']
        db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish')
        genome_info = db_collection.find_one({'name' : self.option('ref_genome'), 'assembly' : self.option('genome_version'), 'annot_version' : self.option('genome_annot_version')})
        self.ref_gtf = os.path.join(db_path, genome_info['gtf'])
        self.ref_annot_dir = os.path.join(db_path, genome_info["anno_path_v2"])
        self.g2t2p = os.path.join(db_path, genome_info['g2t2p'])
        self.des = os.path.join(db_path, genome_info["bio_mart_annot"])
        self.des_type = genome_info["biomart_gene_annotype"]
        self.entrez = os.path.join(db_path, genome_info["ensemble2entrez"])
        self.known_ko = os.path.join(db_path, genome_info['kegg'])
        self.ref_genome = os.path.join(db_path, genome_info["dna_fa"])
        self.ref_fasta = os.path.join(db_path, genome_info["dna_fa"])

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        if self.option('kegg_database') == 'All':
            self.option('kegg_database','All')
        elif self.option('kegg_database') == 'Animal':
            self.option('kegg_database','Animals')
        elif self.option('kegg_database') == 'Plant':
            self.option('kegg_database','Plants')
        elif self.option('kegg_database') == 'Protist':
            self.option('kegg_database','Protists')
        if self.option('nr_database') == 'All':
            self.option('nr_database', 'nr')
        elif self.option('nr_database') == 'Animal':
            self.option('nr_database', 'metazoa')
        elif self.option('nr_database') == 'Plant':
            self.option('nr_database', 'viridiplantae')
        else:
            self.option('nr_database', self.option('nr_database').lower())
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        # add analysis items
        self.filecheck = self.add_tool('ref_rna_v2.file_check')
        self.qc = self.add_module('ref_rna_v2.hiseq_qc')
        self.qc_stat_before = self.add_module('ref_rna_v2.hiseq_reads_stat')
        self.qc_stat_after = self.add_module('ref_rna_v2.hiseq_reads_stat')
        self.annot_filter_ref = self.add_module('ref_rna_v2.annot_filter')
        self.mapping = self.add_module('ref_rna_v2.rnaseq_mapping')
        self.assembly = self.add_module('ref_rna_v2.refrna_assemble')
        self.annot_mapdb = self.add_module('ref_rna_v2.annot_mapdb')
        self.annot_orfpfam = self.add_module('ref_rna_v2.annot_orfpfam')
        self.annot_filter_new = self.add_module('ref_rna_v2.annot_filter')
        self.annot_class_ref = self.add_module('ref_rna_v2.annot_class_beta')
        self.annot_class_new = self.add_module('ref_rna_v2.annot_class_beta')
        self.annot_merge = self.add_tool('ref_rna_v2.annotation.annot_merge')
        # set operation logic
        self.filecheck.on('end', self.run_qc)
        self.filecheck.on('end', self.run_qc_stat, False)
        self.filecheck.on('end', self.run_annot_filter_ref)
        self.annot_filter_ref.on('end', self.run_annot_class_ref)
        self.qc.on('end', self.run_qc_stat, True)
        self.qc.on('end', self.run_mapping)
        self.mapping.on('end', self.run_assembly)
        self.assembly.on('end', self.run_annot_mapdb)
        self.assembly.on('end', self.run_annot_orfpfam)
        self.on_rely([self.annot_mapdb, self.annot_orfpfam], self.run_annot_filter_new)
        self.annot_filter_new.on('end', self.run_annot_class_new)
        self.on_rely([self.annot_class_new, self.annot_class_ref], self.run_annot_merge)
        self.run_filecheck()
        super(RefTestWorkflow, self).run()

    def run_filecheck(self):
        self.step.add_steps('filecheck')
        options = {
            'sample_num': self.option('sample_num'),
            'fq_type': self.option('fq_type'),
            'fastq_dir': self.option('fastq_dir'),
            'in_gtf': self.ref_gtf,
            'group_table': self.option('group_table'),
            'control_file': self.option('control_file')
        }
        self.filecheck.set_options(options)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_qc(self):
        self.step.add_steps('rna_qc')

        self.qc.set_options({
            'fq_type': self.option('fq_type'),
            'fastq_dir': self.option('fastq_dir'),
            'quality_score_system': self.option('quality_score_system')
        })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.rna_qc})
        self.qc.on('end', self.set_step, {'end': self.step.rna_qc})
        self.qc.run()

    def run_qc_stat(self, event):
        if self.option('quality_score_system') == "phred+33":
            quality = 33
        elif self.option('quality_score_system') == "phred+64":
            quality = 64
        if event['data']:
            self.step.add_steps('qc_stat_after')

            self.qc_stat_after.set_options({
                'fastq_dir': self.qc.option('sickle_dir'),
                'fq_type': self.option('fq_type'),
                'quality': quality,
                'dup': True
            })
        else:
            self.step.add_steps('rna_qc')

            self.qc_stat_before.set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type'),
                'quality': quality,
            })
        if event['data']:
            self.step.add_steps('qc_stat_after')
            self.qc_stat_after.on('end', self.set_output, 'qc_stat_after')
            self.qc_stat_after.on("start", self.set_step, {"start": self.step.qc_stat_after})
            self.qc_stat_after.on("end", self.set_step, {"end": self.step.qc_stat_after})
            self.qc_stat_after.run()
        else:
            self.step.add_steps('qc_stat_before')
            self.qc_stat_before.on('end', self.set_output, 'qc_stat_before')
            self.qc_stat_before.on("start", self.set_step, {"start": self.step.qc_stat_before})
            self.qc_stat_before.on("end", self.set_step, {"end": self.step.qc_stat_before})
            self.qc_stat_before.run()

    def run_annot_filter_ref(self):
        self.step.add_steps('annot_filter_ref')
        options = {
            "blast_nr_xml" : os.path.join(self.ref_annot_dir, 'annot_mapdb/nr/blast.xml'),
            "blast_eggnog_xml" : os.path.join(self.ref_annot_dir, 'annot_mapdb/eggnog/blast.xml'),
            "blast_kegg_xml": os.path.join(self.ref_annot_dir, 'annot_mapdb/kegg/blast.xml'),
            "blast_swissprot_xml" : os.path.join(self.ref_annot_dir, 'annot_mapdb/swissprot/blast.xml'),
            "pfam_domain" : os.path.join(self.ref_annot_dir, 'annot_orfpfam/pfam_domain'),
            "blast2go_annot" : os.path.join(self.ref_annot_dir, 'annot_mapdb/GO/blast2go_merge.xls'),
            'nr_evalue': self.option('nr_blast_evalue'),
            'nr_identity': self.option('nr_blast_identity'),
            'nr_similarity': self.option('nr_blast_similarity'),
            'swissprot_evalue': self.option('swissprot_blast_evalue'),
            'swissprot_identity': self.option('swissprot_blast_identity'),
            'swissprot_similarity': self.option('swissprot_blast_similarity'),
            'eggnog_evalue': self.option('cog_blast_evalue'),
            'eggnog_identity': self.option('cog_blast_identity'),
            'eggnog_similarity': self.option('cog_blast_similarity'),
            'kegg_evalue': self.option('kegg_blast_evalue'),
            'kegg_identity': self.option('kegg_blast_identity'),
            'kegg_similarity': self.option('kegg_blast_similarity'),
            'pfam_evalue': self.option('pfam_blast_evalue'),
        }
        self.annot_filter_ref.set_options(options)
        self.annot_filter_ref.on('start', self.set_step, {'start': self.step.annot_filter_ref})
        self.annot_filter_ref.on('end', self.set_step, {'end': self.step.annot_filter_ref})
        self.annot_filter_ref.on('end', self.set_output, "annot_filter_ref")
        self.annot_filter_ref.run()

    def run_annot_class_ref(self):
        self.step.add_steps('annot_class_ref')
        self.ref_filter_dir = self.annot_filter_ref.output_dir
        options = {
            'type': 'ref',
            'gtf': self.ref_gtf,
            'fasta': self.ref_fasta,
            'g2t2p': self.g2t2p,
            'des': self.des,
            'des_type': self.des_type,
            'enterz': self.entrez,
            'blast_nr_xml': os.path.join(self.ref_filter_dir, 'nr/blast.xml.filter.xml'),
            'blast_swissprot_xml': os.path.join(self.ref_filter_dir, 'swissprot/blast.xml.filter.xml'),
            'blast_eggnog_xml': os.path.join(self.ref_filter_dir, 'eggnog/blast.xml.filter.xml'),
            'blast_kegg_xml': os.path.join(self.ref_filter_dir, 'kegg/blast.xml.filter.xml'),
            'known_ko': self.known_ko,
            'taxonomy': self.option('kegg_database'),
            'link_bgcolor': 'yellow',
            'png_bgcolor': 'FFFF00',
            'blast2go_annot': os.path.join(self.ref_filter_dir, 'go/blast2go_merge.xls.filter.xls'),
            'pfam_domain': os.path.join(self.ref_filter_dir, 'pfam/pfam_domain.filter.xls'),
        }
        self.annot_class_ref.set_options(options)
        self.annot_class_ref.on('start', self.set_step, {'start': self.step.annot_class_ref})
        self.annot_class_ref.on('end', self.set_step, {'end': self.step.annot_class_ref})
        self.annot_class_ref.on('end', self.set_output, "annot_class_ref")
        self.annot_class_ref.run()

    def run_mapping(self):
        self.step.add_steps('mapping')
        options = {
            "ref_genome": self.option("ref_genome"),
            "genome_version": self.option("genome_version"),
            "genome_annot_version": self.option("genome_annot_version"),
            "mapping_method": self.option("align_method"),
            "seq_method": self.option("fq_type"),
            "fastq_dir": self.qc.option("sickle_dir"),
            "assemble_method": self.option("assemble_method"),
        }
        if self.option("strand_specific"):
            options.update({"strand_specific":True})
        self.mapping.set_options(options)
        self.mapping.on("end", self.set_output, "mapping")
        self.mapping.on("start", self.set_step, {"start": self.step.mapping})
        self.mapping.on("end", self.set_step, {"end": self.step.mapping})
        self.mapping.run()

    def run_assembly(self):
        self.step.add_steps('assembly')

        options = {
            "sample_bam_dir": self.mapping.option("bam_output"),
            "assemble_method": self.option("assemble_method").lower(),
            "ref_gtf": self.filecheck.option("gtf"),
            "ref_fa": self.ref_genome,
        }
        if self.option("strand_specific"):
            if self.option("strand_dir") == "forward":
                strand_dir = "firststrand"
            else:
                strand_dir = "secondstrand"
            options.update({
                "strand_direct": strand_dir,
                "fr_stranded": "fr-stranded"
                })
        else:
            options.update({
                "fr_stranded": "fr-unstranded"
                })
        self.assembly.set_options(options)
        self.assembly.on("end", self.set_output, "assembly")
        self.assembly.on('start', self.set_step, {'start': self.step.assembly})
        self.assembly.on('end', self.set_step, {'end': self.step.assembly})
        self.assembly.run()

    def run_annot_mapdb(self, event):
        self.step.add_steps('annot_mapdb')
        options = {
            "query" : self.assembly.option("new_transcripts_fa"),
            "nr_db" : self.option("nr_database")
        }
        self.annot_mapdb.set_options(options)
        self.annot_mapdb.on("end", self.set_output, "annot_mapdb")
        self.annot_mapdb.on('start', self.set_step, {'start': self.step.annot_mapdb})
        self.annot_mapdb.on('end', self.set_step, {'end': self.step.annot_mapdb})
        self.annot_mapdb.run()

    def run_annot_orfpfam(self):
        self.step.add_steps('annot_orfpfam')
        options = {
            "fasta" : self.assembly.option("new_transcripts_fa"),
            "gtf" : self.assembly.option("new_transcripts_gtf")
        }
        self.annot_orfpfam.set_options(options)
        self.annot_orfpfam.on("end", self.set_output, "annot_orfpfam")
        self.annot_orfpfam.on('start', self.set_step, {'start': self.step.annot_orfpfam})
        self.annot_orfpfam.on('end', self.set_step, {'end': self.step.annot_orfpfam})
        self.annot_orfpfam.run()

    def run_annot_filter_new(self):
        self.step.add_steps('annot_filter_new')
        options = {
            "blast_nr_xml" : os.path.join(self.annot_mapdb.output_dir, 'nr/blast.xml'),
            "blast_eggnog_xml" : os.path.join(self.annot_mapdb.output_dir, 'eggnog/blast.xml'),
            "blast_kegg_xml": os.path.join(self.annot_mapdb.output_dir, 'kegg/blast.xml'),
            "blast_swissprot_xml" : os.path.join(self.annot_mapdb.output_dir, 'swissprot/blast.xml'),
            "pfam_domain" : os.path.join(self.annot_orfpfam.output_dir, 'pfam_domain'),
            "blast2go_annot" : os.path.join(self.annot_mapdb.output_dir, 'GO/blast2go_merge.xls'),
            'nr_evalue': self.option('nr_blast_evalue'),
            'nr_identity': self.option('nr_blast_identity'),
            'nr_similarity': self.option('nr_blast_similarity'),
            'swissprot_evalue': self.option('swissprot_blast_evalue'),
            'swissprot_identity': self.option('swissprot_blast_identity'),
            'swissprot_similarity': self.option('swissprot_blast_similarity'),
            'eggnog_evalue': self.option('cog_blast_evalue'),
            'eggnog_identity': self.option('cog_blast_identity'),
            'eggnog_similarity': self.option('cog_blast_similarity'),
            'kegg_evalue': self.option('kegg_blast_evalue'),
            'kegg_identity': self.option('kegg_blast_identity'),
            'kegg_similarity': self.option('kegg_blast_similarity'),
            'pfam_evalue': self.option('pfam_blast_evalue'),
        }
        self.annot_filter_new.set_options(options)
        self.annot_filter_new.on('start', self.set_step, {'start': self.step.annot_filter_new})
        self.annot_filter_new.on('end', self.set_step, {'end': self.step.annot_filter_new})
        self.annot_filter_new.on('end', self.set_output, "annot_filter_new")
        self.annot_filter_new.run()

    def run_annot_class_new(self):
        self.step.add_steps('annot_class_new')
        self.new_filter_dir = self.annot_filter_new.output_dir
        options = {
            'type': 'new',
            'gene2trans': os.path.join(self.annot_orfpfam.output_dir, 'all_tran2gen.txt'),
            'des': self.des,
            'des_type': self.des_type,
            'enterz': self.entrez,
            'blast_nr_xml': os.path.join(self.new_filter_dir, 'nr/blast.xml.filter.xml'),
            'blast_swissprot_xml': os.path.join(self.new_filter_dir, 'swissprot/blast.xml.filter.xml'),
            'blast_eggnog_xml': os.path.join(self.new_filter_dir, 'eggnog/blast.xml.filter.xml'),
            'blast_kegg_xml': os.path.join(self.new_filter_dir, 'kegg/blast.xml.filter.xml'),
            'taxonomy': self.option("kegg_database"),
            'link_bgcolor': 'green',
            'png_bgcolor': '00CD00',
            'pfam_domain': os.path.join(self.new_filter_dir, 'pfam/pfam_domain.filter.xls'),
            'blast2go_annot': os.path.join(self.new_filter_dir, 'go/blast2go_merge.xls.filter.xls'),
        }
        self.annot_class_new.set_options(options)
        self.annot_class_new.on('start', self.set_step, {'start': self.step.annot_class_new})
        self.annot_class_new.on('end', self.set_step, {'end': self.step.annot_class_new})
        self.annot_class_new.on('end', self.set_output, "annot_class_new")
        self.annot_class_new.run()

    def run_annot_merge(self):
        self.step.add_steps('annot_merge')
        options = {
            'ref_class_dir': self.annot_class_ref.output_dir,
            'new_class_dir': self.annot_class_new.output_dir,
            'new_mapdb_dir': self.annot_mapdb.output_dir,
        }
        self.annot_merge.set_options(options)
        self.annot_merge.on('start', self.set_step, {'start': self.step.annot_merge})
        self.annot_merge.on('end', self.set_step, {'end': self.step.annot_merge})
        self.annot_merge.on('end', self.end, "annot_merge")
        self.annot_merge.run()

    def set_output(self, event):
        self.logger.info('start set_output for {}'.format(event))
        self.logger.info('finish set_output for {}'.format(event))

    def set_db(self):
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))

    def end(self):
        super(RefTestWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test_wf(self):
        from mbio.workflows.ref_rna_v2.report.ref_test import RefTestWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.ref_test',
            'options': {
                'sample_num': 'multiple',
                'fq_type': 'PE',
                'quality_score_system': 'phred+33',
                'strand_specific': False,
                'ref_genome': 'Homo_sapiens',
                'genome_version': 'GRCh38.p10',
                'genome_annot_version': 'Ensemble_release_89',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/remote_input/fastq_dir/raw_data',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/remote_input/group_table/example_group.txt',
                'control_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/remote_input/control_file/example_control.txt',
                'level': 'transcript',
                'align_method': 'hisat',
                'assemble_method': 'stringtie',
                'nr_database': 'Animal',
                'kegg_database': 'Animala',
            }
        }
        wsheet = Sheet(data=data)
        wf = RefTestWorkflow(wsheet)
        wf.sheet.id = 'ref_rna_v2_upgrade'
        wf.sheet.project_sn = 'ref_rna_v2_upgrade'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

    def test_wpm(self):
        from biocluster.wpm.client import worker_client
        import random
        worker = worker_client()
        number = random.randint(1000, 9999)
        data = {
            'id': 'workflow_{}'.format(number),
            'project_sn': 'workflow_{}'.format(number),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.ref_test',
            'options': {
                'sample_num': 'multiple',
                'fq_type': 'PE',
                'quality_score_system': 'phred+33',
                'strand_specific': False,
                'ref_genome': 'Homo_sapiens',
                'genome_version': 'GRCh38.p10',
                'genome_annot_version': 'Ensemble_release_89',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/remote_input/fastq_dir/raw_data',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/remote_input/group_table/example_group.txt',
                'control_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/remote_input/control_file/example_control.txt',
                'level': 'transcript',
                'align_method': 'hisat',
                'assemble_method': 'stringtie',
                'nr_database': 'Animal',
                'kegg_database': 'Animal',
            }
        }
        info = worker.add_task(data)
        print info

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_wf')])
    unittest.TextTestRunner(verbosity=2).run(suite)
