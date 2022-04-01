# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from Bio import SeqIO
from biocluster.config import Config
from biocluster.module import Module


class LargeGushModule(Module):
    def __init__(self, work_id):
        super(LargeGushModule, self).__init__(work_id)
        options = [
            {'name': 'organism_name', 'type': 'string', 'default': None},
            {'name': 'cpc', 'type': 'string', 'default': 'True'},
            {'name': 'cnci', 'type': 'string', 'default': 'True'},
            {'name': 'cpat', 'type': 'string', 'default': 'True'},
            {'name': 'pfamscan', 'type': 'string', 'default': 'True'},

            {'name': 'identify_num', 'type': 'int', 'default': 2},
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},

            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'taxonmy', 'type': 'string', 'default': None},
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},

            {'name': 'des', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'des_type', 'type': 'string', 'default': 'type1'},

            {'name': 'new_assembly_fasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_assembly_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'ref_mrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'ref_lncrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},

            {'name': 'ref_mrna_fasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'ref_lncrna_fasta', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},

            {'name': 'program', 'type': 'string', 'default': 'salmon'},
            {'name': 'strand_specific', 'type': 'bool', 'default': True},
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'filter', 'type': 'bool', 'default': True},
            {'name': 'threshold', 'type': 'float', 'default': 0.0},

            {'name': 'ids_mapping', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'database', 'type': 'string', 'default': 'ensembl'},
            {'name': 'known_ko', 'type': 'infile', 'format': 'whole_transcriptome.common'},

            {'name': 'new_gene_list', 'type': 'infile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)
        self.config = Config()

        global WORK_DIR
        WORK_DIR = self.work_dir

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(LargeGushModule, self).run()
        self.run_new_lncrna_predict()

    def run_new_lncrna_predict(self):
        self.new_lncrna_predict = self.add_module('lnc_rna.new_lncrna_predict')
        biomart = self.option('des').path
        biomart_type = self.option('des_type')
        mrna_gtf = self.option('ref_mrna_gtf').path
        lnc_db_gtf = self.option('ref_lncrna_gtf').path
        new_fasta = self.option('new_assembly_fasta').path
        new_gtf = self.option('new_assembly_gtf').path
        opts = {
            'biomart': biomart,
            'biomart_type': biomart_type,
            'mrna_gtf': mrna_gtf,
            'lnc_db_gtf': lnc_db_gtf,
            'cpc': self.option('cpc') == 'True',
            'cnci': self.option('cnci') == 'True',
            'cpat': self.option('cpat') == 'True',
            'pfamscan': self.option('pfamscan') == 'True',
            'new_fasta': new_fasta,
            'new_gtf': new_gtf,
            'identify_num': self.option('identify_num'),
            'transcript_len': self.option('transcript_len'),
            'exon_num': self.option('exon_num'),
            'orf_len': self.option('orf_len'),
            'cpc_score': self.option('cpc_score'),
            'cnci_score': self.option('cnci_score'),
            'taxonmy': self.option('taxonmy'),
            'cpat_score': self.option('cpat_score'),
        }
        if self.option('organism_name') in ['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster']:
            common_name = {
                'Homo_sapiens': 'Human',
                'Mus_musculus': 'Mouse',
                'Danio_rerio': 'Zebrafish',
                'Drosophila_melanogaster': 'Fly'
            }[self.option('organism_name')]
            hexamer_dat = os.path.join(self.config.SOFTWARE_DIR,
                                       'bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_Hexamer.tsv'.format(common_name))
            logit_model = os.path.join(self.config.SOFTWARE_DIR,
                                       'bioinfo/lnc_rna/CPAT-1.2.4/dat/{}_logitModel.RData'.format(common_name))
            opts.update({'hexamer_dat': hexamer_dat, 'logit_model': logit_model})
        self.new_lncrna_predict.set_options(opts)
        self.new_lncrna_predict.on('end', self.run_merge_known_new)
        self.new_lncrna_predict.on('end', self.set_output, 'new_lncrna_predict')
        self.new_lncrna_predict.run()

    def run_merge_known_new(self):
        self.merge_known_new = self.add_tool('lnc_rna.merge_known_new')
        all_known_gtf = combine_gtf([self.option('ref_mrna_gtf').path, self.option('ref_lncrna_gtf').path])
        all_known_fa = combine_fasta([self.option('ref_mrna_fasta').path, self.option('ref_lncrna_fasta').path])
        new_mrna_gtf = os.path.join(self.new_lncrna_predict.output_dir, 'novel_mrna.gtf')
        new_mrna_fa = os.path.join(self.new_lncrna_predict.output_dir, 'novel_mrna.fa')
        new_lncrna_gtf = os.path.join(self.new_lncrna_predict.output_dir, 'novel_lncrna.gtf')
        new_lncrna_fa = os.path.join(self.new_lncrna_predict.output_dir, 'novel_lncrna.fa')
        opts = {
            'all_known_gtf': all_known_gtf,
            'all_known_fa': all_known_fa,
            'new_mrna_gtf': new_mrna_gtf,
            'new_mrna_fa': new_mrna_fa,
            'new_lncrna_gtf': new_lncrna_gtf,
            'new_lncrna_fa': new_lncrna_fa,
        }
        self.merge_known_new.set_options(opts)
        self.merge_known_new.on('end', self.run_expression)
        self.merge_known_new.on('end', self.set_output, 'merge_known_new')
        self.merge_known_new.run()

    def run_expression(self):
        self.expression = self.add_module('whole_transcriptome.expression')
        transcripts = os.path.join(self.merge_known_new.output_dir, 'known_and_new.fa')
        t2g = os.path.join(self.merge_known_new.output_dir, 't2g.txt')
        opts = {
            'program': self.option('program'),
            'transcripts': transcripts,
            'strand_specific': self.option('strand_specific'),
            'strand_dir': self.option('strand_dir'),
            'fastq_dir': self.option('fastq_dir'),
            't2g': t2g
        }
        self.expression.set_options(opts)
        self.known_lnc_identify = self.add_module('lnc_rna.known_lnc_identify')
        self.filter_by_express = self.add_tool('lnc_rna.filter_by_express')
        self.on_rely([self.known_lnc_identify, self.filter_by_express], self.run_lncrna_stat)
        self.expression.on('end', self.run_known_lnc_identify)
        self.expression.on('end', self.run_filter_by_express)
        self.expression.on('end', self.set_output, 'expression')
        self.expression.run()

    def run_known_lnc_identify(self):
        biomart = self.option('des').path
        biomart_type = self.option('des_type')
        lnc_db_gtf = self.option('ref_lncrna_gtf').path
        lnc_db_fa = self.option('ref_lncrna_fasta').path
        mrna_gtf = self.option('ref_mrna_gtf').path
        mrna_fa = self.option('ref_mrna_fasta').path
        exp_matrix = os.path.join(self.expression.output_dir, 'T.tpm.txt')
        ids_mapping = self.option('ids_mapping').path
        opts = {
            'biomart': biomart,
            'biomart_type': biomart_type,
            'lnc_db_gtf': lnc_db_gtf,
            'lnc_db_fa': lnc_db_fa,
            'mrna_gtf': mrna_gtf,
            'mrna_fa': mrna_fa,
            'exp_matrix': exp_matrix,
            'ids_mapping': ids_mapping,
            'database': self.option('database'),
        }
        self.known_lnc_identify.set_options(opts)
        self.known_lnc_identify.on('end', self.set_output, 'known_lnc_identify')
        self.known_lnc_identify.run()

    def run_filter_by_express(self):
        known_mrna_fa = self.option('ref_mrna_fasta').path
        known_lncrna_fa = self.option('ref_lncrna_fasta').path
        known_mrna_gtf = self.option('ref_mrna_gtf').path
        known_lncrna_gtf = self.option('ref_lncrna_gtf').path
        novel_mrna_fa = os.path.join(self.new_lncrna_predict.output_dir, 'novel_mrna.fa')
        novel_lncrna_fa = os.path.join(self.new_lncrna_predict.output_dir, 'novel_lncrna.fa')
        novel_mrna_gtf = os.path.join(self.new_lncrna_predict.output_dir, 'novel_mrna.gtf')
        novel_lncrna_gtf = os.path.join(self.new_lncrna_predict.output_dir, 'novel_lncrna.gtf')
        tpm_matrix = os.path.join(self.expression.output_dir, 'T.tpm.txt')
        fpkm_matrix = os.path.join(self.expression.output_dir, 'T.fpkm.txt')
        count_matrix = os.path.join(self.expression.output_dir, 'T.count.txt')
        tpm_matrix_g = os.path.join(self.expression.output_dir, 'G.tpm.txt')
        fpkm_matrix_g = os.path.join(self.expression.output_dir, 'G.fpkm.txt')
        count_matrix_g = os.path.join(self.expression.output_dir, 'G.count.txt')
        lnc_new_dir = self.new_lncrna_predict.output_dir
        opts = {
            'known_mrna_fa': known_mrna_fa,
            'known_lncrna_fa': known_lncrna_fa,
            'known_mrna_gtf': known_mrna_gtf,
            'known_lncrna_gtf': known_lncrna_gtf,
            'novel_mrna_fa': novel_mrna_fa,
            'novel_lncrna_fa': novel_lncrna_fa,
            'novel_mrna_gtf': novel_mrna_gtf,
            'novel_lncrna_gtf': novel_lncrna_gtf,
            'tpm_matrix': tpm_matrix,
            'fpkm_matrix': fpkm_matrix,
            'count_matrix': count_matrix,
            'tpm_matrix_g': tpm_matrix_g,
            'fpkm_matrix_g': fpkm_matrix_g,
            'count_matrix_g': count_matrix_g,
            'lnc_new_dir': lnc_new_dir
        }
        if self.option('known_ko').is_set:
            opts.update({'known_ko': self.option('known_ko')})
        self.filter_by_express.set_options(opts)
        self.filter_by_express.on('end', self.set_output, 'filter_by_express')
        self.filter_by_express.run()

    def run_lncrna_stat(self):
        self.lncrna_stat = self.add_tool('lnc_rna.lncrna_identification.lncrna_stat')
        novel_lncrna_detail = os.path.join(self.filter_by_express.output_dir, 'filtered_lncnovel/novel_lncrna_predict_detail.xls')
        known_lncrna_detail = os.path.join(self.known_lnc_identify.output_dir, 'known_lncrna_detail.xls')
        exp_matrix = self.expression.output_dir + "/T.tpm.txt",
        gene_type = os.path.join(self.filter_by_express.output_dir, 'filtered_file/gene_type.xls')
        new_gene_list = self.option('new_gene_list').path
        opts = {
            'novel_lncrna_detail': novel_lncrna_detail,
            'known_lncrna_detail': known_lncrna_detail,
            'exp_matrix': exp_matrix,
            'gene_type': gene_type,
            'new_gene_list': new_gene_list
        }
        self.lncrna_stat.set_options(opts)
        self.lncrna_stat.on('end', self.set_output, 'lncrna_stat')
        self.lncrna_stat.run()

    def set_output(self, event):
        obj = event['bind_object']
        shutil.copytree(obj.output_dir, os.path.join(self.output_dir, event['data']))
        if event['data'] == 'lncrna_stat':
            self.end()

    def end(self):
        super(LargeGushModule, self).end()


def combine_gtf(gtf_list, drop_comment=True):
    ret_file = os.path.join(WORK_DIR, 'combined.gtf')
    with open(ret_file, 'w') as handle:
        for gtf in gtf_list:
            for line in open(gtf):
                if line.strip():
                    if drop_comment and line[0] == '#':
                        continue
                    handle.write(line)
    return ret_file


def combine_fasta(fasta_list):
    ret_file = os.path.join(WORK_DIR, 'combined.fasta')
    sequences = list()
    for fasta in fasta_list:
        sequences.extend(SeqIO.parse(fasta, 'fasta'))
    else:
        SeqIO.write(sequences, ret_file, 'fasta')
        return ret_file


class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'large_gush_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.large_gush',
            'instant': False,
            'options': {
                'organism_name': 'Homo_sapiens',
                'taxonmy': 'Animal',
                'des': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/biomart/biomart.txt',
                'des_type': 'type1',
                'new_assembly_fasta': '/mnt/ilustre/users/sanger-dev/workspace/20190917/WholeTranscriptome_workflow_6312_7216/Assembly/output/new.fasta',
                'new_assembly_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20190917/WholeTranscriptome_workflow_6312_7216/Assembly/output/new.gtf',
                'ref_mrna_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/mrna.gtf',
                'ref_lncrna_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/lncrna.gtf',
                'ref_mrna_fasta': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/mrna.fa',
                'ref_lncrna_fasta': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/lncrna.fa',

                'program': 'salmon',
                'strand_specific': True,
                'strand_dir': 'RF',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/workspace/20190917/WholeTranscriptome_workflow_6312_7216/FastpRna/output/fastq',
                'filter': True,
                'threshold': 0.0,

                'ids_mapping': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/ids_matrix.xls',

                'new_gene_list': '/mnt/ilustre/users/sanger-dev/workspace/20190917/WholeTranscriptome_workflow_6312_7216/Assembly/output/new.gene_id.list',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
