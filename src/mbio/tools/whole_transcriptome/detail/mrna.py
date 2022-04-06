# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os, re
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class MrnaAgent(Agent):
    '''
    last_modify: 2019.10.31
    '''

    def __init__(self, parent):
        super(MrnaAgent, self).__init__(parent)
        options = [
            {'name': 'gene_type', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'trans_type', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'all_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'mrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'dna_fa', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'relation_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_type', 'type': 'string', 'default': None},
            {'name': 'ref_cds', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'ref_pep', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_cds', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_pep', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'genome_id', 'type': 'string', 'default': None},
            {'name': 'transcript_bed', 'type': 'outfile', 'format': 'gene_structure.bed'},
            {'name': 'transcript_fa', 'type': 'outfile', 'format': 'sequence.fasta'},
            {'name': 'gene_bed', 'type': 'outfile', 'format': 'gene_structure.bed'},
            {'name': 'gene_fa', 'type': 'outfile', 'format': 'sequence.fasta'}
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('dna_fa').path) / 1024 ** 3 + 16)

    def end(self):
        super(MrnaAgent, self).end()


class MrnaTool(Tool):
    def __init__(self, config):
        super(MrnaTool, self).__init__(config)
        self.program = {
            'gffread': 'bioinfo/rna/cufflinks-2.2.1/gffread',
            'gtftogenepred': 'bioinfo/align/ucsc_tools/gtfToGenePred',
            'genepredtobed': 'bioinfo/align/ucsc_tools/genePredToBed',
            'python': 'miniconda2/bin/python',
            'bedtools': 'bioinfo/ref_rna_v2/miniconda2/bin/bedtools'
        }
        self.script = {
            'get_gene_bed': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/get_gene_bed.py'),
            'fasta_clean': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/detail/fasta_clean.py'),
            'base_build': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/detail/base_build.py')
        }
        self.file = {
            'transcript_fasta': os.path.join(self.output_dir, 'transcript.fasta'),
            'transcript_genepred': os.path.join(self.work_dir, 'transcript.genepred'),
            'transcript_bed': os.path.join(self.output_dir, 'transcript.bed'),
            'gene_bed': os.path.join(self.output_dir, 'gene.bed'),
            'gene_fasta': os.path.join(self.output_dir, 'gene.fasta'),
            'input_json': os.path.join(self.work_dir, 'input.json'),
            'database': os.path.join(self.output_dir, 'mrna.db'),
            'gene_detail': os.path.join(self.output_dir, 'gene_detail.pk'),
            'transcript_detail': os.path.join(self.output_dir, 'transcript_detail.pk')
        }

    def run(self):
        super(MrnaTool, self).run()
        if self.option("all_gtf").is_set:
            self.filter_gtf()
        self.run_gffread()
        self.run_gtftogenepred()
        self.run_genepredtobed()
        self.run_get_gene_bed()
        self.run_bedtools_getfasta()
        self.run_fasta_clean()
        self.pre_base_build()
        self.run_base_build()
        self.set_output()
        self.end()

    def filter_gtf(self):
        mrna_gene = dict()
        with open(self.option("gene_type").prop['path'], "r") as f:
            for line in f:
                items = line.strip().split("\t")
                if items[2].lower() == "mrna":
                    mrna_gene[items[0]] = items[1]
        filtered_gtf = os.path.join(self.work_dir, "filter.gtf")
        with open(self.option("all_gtf").path, "r") as f1, open(filtered_gtf, "w") as w:
            pg = re.compile(r'gene_id "(\S+)";')
            for line in f1:
                if line.strip() and line[0] != '#':
                    eles = line.strip().split('\t')
                    if len(eles) >= 8 and 'gene_id' in eles[8]:
                        mg = re.search(pg, eles[8])
                        if mg:
                            g = mg.group(1)
                            if g in mrna_gene:
                                w.write(line)
        self.option('mrna_gtf', filtered_gtf)

    def run_gffread(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.option('mrna_gtf').path)
        cmd += ' -g {}'.format(self.option('dna_fa').path)
        cmd += ' -w {}'.format(self.file['transcript_fasta'])
        runcmd(self, 'run_gffread', cmd, block=False)

    def run_gtftogenepred(self):
        cmd = '{} {} {}'.format(
            self.program['gtftogenepred'], self.option('mrna_gtf').path, self.file['transcript_genepred']
        )
        runcmd(self, 'run_gtftogenepred', cmd)

    def run_genepredtobed(self):
        cmd = '{} {} {}'.format(self.program['genepredtobed'], self.file['transcript_genepred'],
                                self.file['transcript_bed'])
        runcmd(self, 'run_genepredtobed', cmd)

    def run_get_gene_bed(self):
        cmd = '{} {}'.format(self.program['python'], self.script['get_gene_bed'])
        cmd += ' --gtf {}'.format(self.option('mrna_gtf').path)
        cmd += ' --bed {}'.format(self.file['transcript_bed'])
        cmd += ' --output {}'.format(self.file['gene_bed'])
        runcmd(self, 'run_get_gene_bed', cmd)

    def run_bedtools_getfasta(self):
        cmd = '{} getfasta -name -s'.format(self.program['bedtools'])
        cmd += ' -fi {}'.format(self.option('dna_fa').path)
        cmd += ' -fo {}'.format(self.file['gene_fasta'])
        cmd += ' -bed {}'.format(self.file['gene_bed'])
        runcmd(self, 'run_bedtools_getfasta', cmd)

    def run_fasta_clean(self):
        cmd = '{} {}'.format(self.program['python'], self.script['fasta_clean'])
        cmd += ' -i {}'.format(self.file['gene_fasta'])
        cmd += ' -o {}'.format(self.file['gene_fasta'])
        runcmd(self, 'run_fasta_clean', cmd)

    def pre_base_build(self):
        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
        collection = database['sg_genome_db']
        genome_doc = collection.find_one({'genome_id': self.option('genome_id')})
        organism_name = genome_doc['organism_name']
        if 'ensembl' in genome_doc['ensemble_web'].lower():
            source = 'ensembl'
        elif 'ncbi' in genome_doc['ensemble_web'].lower():
            source = 'ncbi'
        else:
            source = None
        species_urls = genome_doc['ensemble_web']
        entrez = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/" +  genome_doc['ensemble2entrez']

        json.dump({'mrna_gtf': self.option('mrna_gtf').path,
                   'biomart_file': self.option('biomart_file').path,
                   'biomart_type': self.option('biomart_type'),
                   'entrez': entrez,
                   'relation_file': self.option('relation_file').path,
                   'ref_cds_fa': self.option('ref_cds').path,
                   'new_cds_fa': self.option('new_cds').path,
                   'ref_pep_fa': self.option('ref_pep').path,
                   'new_pep_fa': self.option('new_pep').path,
                   'gene_fa': self.file['gene_fasta'],
                   'transcript_fa': self.file['transcript_fasta'],
                   'gene_bed': self.file['gene_bed'],
                   'transcript_bed': self.file['transcript_bed'],
                   'trans_type': self.option('trans_type').path,
                   'organism_name': organism_name,
                   'source': source,
                   'species_urls': species_urls,
                   'gene_detail_pk': self.file['gene_detail'],
                   'transcript_detail_pk': self.file['transcript_detail']},
                  open(self.file['input_json'], 'w'), indent=4)

    def run_base_build(self):
        cmd = '{} {}'.format(self.program['python'], self.script['base_build'])
        cmd += ' --json {}'.format(self.file['input_json'])
        cmd += ' --database {}'.format(self.file['database'])
        runcmd(self, 'run_base_build', cmd)

    def set_output(self):
        self.option('transcript_bed').set_path(self.file['gene_fasta'])
        self.option('transcript_fa').set_path(self.file['transcript_fasta'])
        self.option('gene_bed').set_path(self.file['gene_bed'])
        self.option('gene_fa').set_path(self.file['gene_fasta'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.detail.mrna',
            'instant': False,
            'options': {
                'mrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/large_gush/filter_by_express/filtered_file/all_mrna.gtf',
                'dna_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'relation_file': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/annotation/allannot_class/all_tran2gene.txt',
                'biomart_file': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt',
                'biomart_type': 'type1',
                'ref_cds': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/cds/Homo_sapiens.GRCh38.cds.all.fa',
                'ref_pep': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/cds/Homo_sapiens.GRCh38.pep.all.fa',
                'new_cds': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/annotation/newannot_orfpfam/novel_mrna.fa.transdecoder.cds',
                'new_pep': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/annotation/newannot_orfpfam/novel_mrna.fa.transdecoder.pep',
                'genome_id': 'GM0259'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
