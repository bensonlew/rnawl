# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import pandas as pd
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class TfbsPredictWorkflow(Workflow):
    """
    TfbsPredict description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TfbsPredictWorkflow, self).__init__(wsheet_object)
        options = list()
        for each in self._sheet.options():
            if each in ["seqdb", "motifs_user", 'gtf']:
                pass
            else:
                options.append(dict(name=each, type="string"))
        for each in ["seqdb", "motifs_user", 'gtf']:
            options.append(dict(name=each, type="infile", format="whole_transcriptome.common"))

        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tfbs_predict = self.add_tool("ref_rna_v2.tfbs_predict")
        self.dump_tool = self.api.api("whole_transcriptome.tfbs_predict")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/03 TF_Analysis')

    def run(self):
        # self.start_listener(); self.fire("start") # if you have no tools, you should use this line
        self.tfbs_predict.on("end", self.set_db)
        self.get_run_log()
        self.run_tfbs_predict()
        super(TfbsPredictWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="tfbs_predict", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tfbs_predict(self):
        # 获得new_gtf
        if self.option("gtf").is_set:
            new_gene_gtf = self.option("gtf").prop['path']
        elif "refall" not in self.option("geneset_id").lower():
            conn = self.dump_tool.db['task']
            find_dict = conn.find_one({"task_id": self.option('task_id')})
            if find_dict:
                if 'fastq' in find_dict:
                    fastq_path = find_dict['fastq']
                else:
                    transcript_path = find_dict['assemble_fa']
            else:
                self.set_error('Failed to find task_id: %s in task table', variables = (self.option('task_id')), code = "13703101")
            if 'fastq' in find_dict:
                workflow_result_path = fastq_path.split('QC/')[0]
            else:
                workflow_result_path = transcript_path.split('Assemble/')[0]
            new_gene_gtf = workflow_result_path + 'Assemble/new_transcripts.gtf'
        else:
            new_gene_gtf = ''

        '''
        if re.match(r'^\w+://\S+/.+$', self.option('motifs_user')):
            motifs_user = self.option('motifs_user').prop['path']
        else:
            motifs_user = self.option('motifs_user')
        '''
        '''
        if re.match(r'^\w+://\S+/.+$', self.option('seqdb')):
            seqdb = self.option('seqdb').prop['path']
        else:
            seqdb = self.option('seqdb')
        '''
        options = dict(
            tf_selected=self.option('tf_selected'),
            genes=self.option('geneset_list'),
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
            genome_id=self.option('genome_id'),
        )
        if self.option('motifs_user').is_set:
            options.update({
                "motifs_user" : self.option('motifs_user').prop['path']
            })
        if self.option('seqdb').is_set:
            options.update({
                "seqdb" : self.option('se_db').prop['path']
            })

        
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
        self.dump_tool.update_db_record('tfbs_predict', self.option('main_id'), output_dir=workflow_output, status="end")
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tfbs_predict.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tfbs_predict.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tfbs_predict.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tfbs_predict.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "转录因子靶基因预测文件",0,"211585"],
            ["./tf_target_number.xls", "", "转录因子靶基因预测统计表", 0],
            ["./tfbs_predict.xls", "", "转录因子靶基因预测详情表", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(TfbsPredictWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

