# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.11.13

import os
import unittest
import glob
from biocluster.module import Module
import shutil
from mbio.packages.lnc_rna.copy_file import CopyFile


class GeneDetailModule(Module):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GeneDetailModule, self).__init__(wsheet_object)
        options = [
            {"name": "rna_type", "type": "string", "default": "mrna,lncrna,circrna,mirna"},
            {'name': 'mrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'all_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'dna_fa', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'gene_type', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'trans_type', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'relation_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'biomart_type', 'type': 'string', 'default': None},
            {'name': 'ref_cds', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'ref_pep', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_cds', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'new_pep', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'genome_id', 'type': 'string', 'default': None},

            {'name': 'lncrna_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            {'name': 'lnc_relation_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'known_lncrna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'novel_lncrna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'known_mirna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'novel_mirna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'circrna_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},

        ]
        self.add_option(options)

        self.mrna = self.add_tool("whole_transcriptome.detail.mrna")
        self.lncrna = self.add_tool("whole_transcriptome.detail.lncrna")
        self.mirna = self.add_tool("whole_transcriptome.detail.mirna")
        self.circrna = self.add_tool("whole_transcriptome.detail.circrna")

        # self.target_predict.on('end', self.set_db)

    def run_mrna(self):
        options = {
            'mrna_gtf': self.option('mrna_gtf'),
            'trans_type': self.option('trans_type'),
            'gene_type': self.option('gene_type'),
            'all_gtf': self.option('all_gtf'),
            'dna_fa': self.option('dna_fa'),
            'relation_file': self.option('relation_file'),
            'biomart_file': self.option('biomart_file'),
            'biomart_type': self.option('biomart_type'),
            'ref_cds': self.option('ref_cds'),
            'ref_pep': self.option('ref_pep'),
            'new_cds': self.option('new_cds'),
            'new_pep': self.option('new_pep'),
            'genome_id': self.option('genome_id')
        }

        self.mrna.on("end", self.set_output, "mrna")
        self.mrna.set_options(options)
        self.mrna.run()

    def run_lncrna(self):
        options = {
            'lncrna_gtf': self.option('lncrna_gtf'),
            'dna_fa': self.option('dna_fa'),
            'relation_file': self.option('lnc_relation_file'),
            'biomart_file': self.option('biomart_file'),
            'biomart_type': self.option('biomart_type'),
            'known_lncrna_detail': self.option('known_lncrna_detail'),
            'novel_lncrna_detail': self.option('novel_lncrna_detail'),
            'genome_id': self.option('genome_id')
        }

        self.lncrna.on("end", self.set_output, "lncrna")
        self.lncrna.set_options(options)
        self.lncrna.run()

    def run_circrna(self):
        options = {
            'biomart_file': self.option('biomart_file'),
            'biomart_type': self.option('biomart_type'),
            'circrna_detail': self.option('circrna_detail'),
        }

        self.circrna.on("end", self.set_output, "circrna")
        self.circrna.set_options(options)
        self.circrna.run()

    def run_mirna(self):
        options = {
            'known_mirna_detail': self.option('known_mirna_detail'),
            'novel_mirna_detail': self.option('novel_mirna_detail'),
        }

        self.mirna.on("end", self.set_output, "mirna")
        self.mirna.set_options(options)
        self.mirna.run()

    def run(self):
        super(GeneDetailModule, self).run()
        tool_list = list()
        for rna_type in self.option("rna_type").split(","):
            tool_list.append(getattr(self, rna_type))
        self.on_rely(tool_list, self.set_db)
        for rna_type in self.option("rna_type").split(","):
            run_tool = getattr(self, "run_" + rna_type)
            run_tool()

    def set_output(self, event):
        obj = event["bind_object"]
        name = event['data']
        CopyFile().linkdir(obj.output_dir, os.path.join(self.output_dir, name))

    def set_db(self):
        seqdetail_path = os.path.join(self.output_dir, "seqdetail/")
        if os.path.exists(seqdetail_path):
            shutil.rmtree(seqdetail_path)
        os.mkdir(seqdetail_path)
        seq_download_file_list = glob.glob(os.path.join(self.output_dir, '*/seqdownload*'))
        for i in seq_download_file_list:
            newfile_path = os.path.join(seqdetail_path, os.path.basename(i))
            CopyFile().linkfile(i, newfile_path)
        self.end()

    def end(self):
        super(GeneDetailModule, self).end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet

            data = {
                'id': 'whole_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': "module",
                'name': 'whole_transcriptome.gene_detail',
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
                    'genome_id': 'GM0259',
                    'known_mirna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer1/output/srna/known_mirna/known_mirna_detail.xls',
                    'novel_mirna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer1/output/srna/novel_mirna/novel_mirna_detail.xls',

                    'lncrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/all_lncrna.gtf',
                    'known_lncrna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/known_lnc_identify/known_lncrna_detail.xls',
                    'novel_lncrna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/new_lncrna_predict/novel_lncrna_predict_detail.xls',
                    'lnc_relation_file': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/trans_type.xls',
                    'circrna_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088/Transfer/output/circ_brush/detail.txt'

                }
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
