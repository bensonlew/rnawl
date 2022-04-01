# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import pandas as pd
import re


class TfbsPredictWorkflow(Workflow):
    """
    TfbsPredict description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TfbsPredictWorkflow, self).__init__(wsheet_object)
        options = list()
        for each in self._sheet.options():
            if each in ["seqdb", "motifs_user"]:
                if re.match(r'^\w+://\S+/.+$', each):
                    options.append(dict(name=each, type="infile", format="ref_rna.common"))
                else:
                    options.append(dict(name=each, type="string"))
            else:
                options.append(dict(name=each, type="string"))

        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tfbs_predict = self.add_tool("rna.tfbs_predict")
        self.dump_tool = self.api.api("ref_rna.tfbs_predict")

    def run(self):
        # self.start_listener(); self.fire("start") # if you have no tools, you should use this line
        self.tfbs_predict.on("end", self.set_db)
        self.run_tfbs_predict()
        super(TfbsPredictWorkflow, self).run()

    def run_tfbs_predict(self):
        # 获得new_gtf
        if "refall" not in self.option("geneset_id"):
            conn = self.dump_tool.db['sg_specimen']
            bam_path = conn.find_one({"about_qc" : "after", "task_id" : self.option('task_id'),})['bam_path']
            workflow_result_path = bam_path.split('Align/')[0]
            new_gene_gtf = workflow_result_path + 'Assemble/AssembleResults/new_trans.gtf'
        else:
            new_gene_gtf = ''
        if re.match(r'^\w+://\S+/.+$', self.option('motifs_user')):
            motifs_user = self.option('motifs_user').prop['path']
        else:
            motifs_user = self.option('motifs_user')
        if re.match(r'^\w+://\S+/.+$', self.option('seqdb')):
            seqdb = self.option('seqdb').prop['path']
        else:
            seqdb = self.option('seqdb')
        options = dict(
            motifs_user=motifs_user,
            tf_selected=self.option('tf_selected'),
            genes=self.option('geneset_list'),
            seqdb=seqdb,
            thresh=self.option('thresh'),
            # qv_thresh=self.option('qv_thresh'),
            exp_matrix=self.option('exp_matrix'),
            geneid2tfid=self.option('tf_predict_id'),
            corr_cutoff=self.option('corr_cutoff'),
            corr_pvalue=self.option('corr_pvalue'),
            # corr_use_padjust=self.option('corr_use_padjust'),
            backward=self.option('backward'),
            forward=self.option('forward'),
            organism=self.option("organism"),  # 用于获gtf和genome
            motif_species=self.option("motif_species"),
            motif_species_type=self.option("motif_species_type"),
            gtf=new_gene_gtf,
        )
        # set option
        self.tfbs_predict.set_options(options)
        # run tool
        self.tfbs_predict.run()

    def set_db(self):
        """
        dump data to db
        """
        workflow_output = self.get_workflow_output_dir()
        # add result info
        self.dump_tool.add_tfbs_predict_detail(self.tfbs_predict.output_dir, self.option("main_id"))
        self.dump_tool.update_db_record('sg_tfbs_predict', self.option('main_id'), output_dir=workflow_output, status="end")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tfbs_predict.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TfbsPredict"],
        ])
        super(TfbsPredictWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

