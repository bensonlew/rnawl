# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-12-11"

from biocluster.workflow import Workflow
import unittest
from biocluster.wpm.client import *
import datetime
import re


class ExpressPrCorrWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpressPrCorrWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='protein_matrix', type='string'),
            dict(name='rna_matrix', type='string'),
            dict(name='geneset_list', type='string'),
            dict(name='proteinset_list', type='string'),
            dict(name='rna_group_dict', type='string'),
            dict(name='protein_group_dict', type='string'),
            dict(name="group_rna", type="string"),
            dict(name="group_protein", type="string"),
            dict(name="rna_type", type='string', default="ref_rna_v2"),
            dict(name="corr_method", type='string', default="spearman"),
            dict(name="cor_cutoff", type='float', default=0.8),
            dict(name="pvalue_cutoff", type='float', default=0.05),
            dict(name="qvalue_cutoff", type='float', default=0.05),
            dict(name="padjust_way", type='string', default="fdr_bh"),
            dict(name="sig_type", type='int', default=1),
            dict(name="main_id", type='string'),
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run_tool(self):
        self.tool = self.add_tool("protein_transcript_dia_v3.express_pr_corr")
        self.tool.on("end", self.set_db)
        options = dict(
            group_rna=self.option('group_rna'),
            group_protein=self.option('group_protein'),
            protein_matrix=self.option('protein_matrix'),
            rna_matrix=self.option('rna_matrix'),
            proteinset_list=self.option('proteinset_list'),
            geneset_list=self.option('geneset_list'),
            cor_cutoff=self.option('cor_cutoff'),
            corr_method=self.option('corr_method'),
            padjust_way=self.option('padjust_way'),
            sig_type=self.option('sig_type'),
            pvalue_cutoff=self.option('pvalue_cutoff'),
            qvalue_cutoff=self.option('qvalue_cutoff'),
        )
        self.tool.set_options(options)
        self.tool.run()

    def run(self):
        self.run_tool()
        super(ExpressPrCorrWorkflow, self).run()

    def set_db(self):
        exp_cor = self.api.api("protein_transcript_dia_v3.express_pr_corr")
        # add result info
        exp_cor.add_expcorr(self.get_workflow_output_dir(), self.tool.output_dir,
                                 main_id=self.option('main_id'), )
        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'^\w+://\S+/.+$', workflow_output):
            workflow_output = workflow_output
        else:
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        return workflow_output

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达相关性分析结果目录"],
        ])
        super(ExpressPrCorrWorkflow, self).end()

class TestFunction(unittest.TestCase):
    def test(self):
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        data = {
            "id": "protein_transcript_labelfree_" + id,
            "type": "workflow",
            "name": "protein_transcript_labelfree.protein_transcript_labelfree",
            "options": dict(
                rna_type="ref_rna_v1",
                gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
                pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
                protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
                transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
                protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
                )
        }

        info = worker.add_task(data)
        print(info)


    # def test(self):
    #     import random
    #     from mbio.workflows.single import SingleWorkflow
    #     from biocluster.wsheet import Sheet
    #     data = {
    #         "id": "protein_transcript_labelfree",
    #         "type": "workflow",
    #         "name": "protein_transcript_labelfree.protein_transcript_labelfree",
    #         "options": dict(
    #             rna_type="ref_rna_v1",
    #             gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/biomart/Homo_sapiens.GRCh37.biomart",
    #             pep="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/ensemble_release_87/cds/Homo_sapiens.GRCh37.pep.fa",
    #             protein_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/protein.list",
    #             transcript_list="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/gene.list",
    #             protein_faa="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript_labelfree/blast_testdata/DB.fasta",
    #             )
    #     }
    #     # data['options']['method'] = 'rsem'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     # #
    #     # data['id'] += '1'
    #     # data['options']['method'] = 'salmon'
    #     # wsheet = Sheet(data=data)
    #     # wf = SingleWorkflow(wsheet)
    #     # wf.run()
    #     #
    #     data['id'] += '_fyt'
    #     wsheet = Sheet(data=data)
    #     wf = SingleWorkflow(wsheet)
    #     wf.run()


if __name__ == '__main__':
    unittest.main()
