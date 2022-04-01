# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd


class TfPredictWorkflow(Workflow):
    """
    TfPredict description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TfPredictWorkflow, self).__init__(wsheet_object)
        options = list()
        for each in self._sheet.options():
            options.append(dict(name=each, type="string"))
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.tf_predict = self.add_tool("rna.tf_predict")
        self.tf_predict = self.add_module("rna.tf_predict")
        self.dump_tool = self.api.api("ref_rna.tf_predict")

    def run(self):
        self.tf_predict.on("end", self.set_db)
        self.run_tf_predict()
        super(TfPredictWorkflow, self).run()

    def run_tf_predict(self):
        options = dict(
            s=self.option('s'),
            organism=self.option('organism'),
            blast_all=self.option('blast_all'),
            seqfile=self.option('seq_db'),
            E=self.option('E'),
            # domE=self.option('domE'),
            evalue=self.option('evalue'),
        )
        self.tf_predict.set_options(options)
        self.tf_predict.run()

    def set_db(self):
        """
        dump data to db
        """
        workflow_output = self.get_workflow_output_dir()
        # modify result info
        seq_annot_pd = pd.read_table(self.work_dir + '/seq_annot.xls', header=0, index_col=0)
        tf_result = pd.read_table(self.tf_predict.output_dir + '/final_tf_predict.xls', header=0, index_col=0)
        seq_annot_pd = tf_result.join(seq_annot_pd)
        seq_annot_pd.index.name = "transcript_id"
        seq_annot_pd.to_csv(self.tf_predict.output_dir + '/final_tf_predict.xls', header=True, index=True, sep='\t')
        if self.option('s').lower() == 'plant':
            known_meme_tf_list_file = self.config.SOFTWARE_DIR + "/database/TFDB/plant_tf_motif/plant.meme.gene.list"
            with open(known_meme_tf_list_file) as f:
                known_meme_tf_list = [x.strip() for x in f]
                known_meme_tf_list = tuple(known_meme_tf_list)
        else:
            known_meme_tf_list = list()
        seq_annot_pd = seq_annot_pd.fillna("")
        tmp_tuple = zip(seq_annot_pd['gene_id'], seq_annot_pd['blast_hit'], seq_annot_pd['family'])
        tf_select_values = list()
        tf_select_shows = list()
        for x in tmp_tuple:
            if x[1] not in known_meme_tf_list:
                if '.' in x[1]:
                    tmp_list = x[1].split(".")
                    for i in range(0,len(tmp_list)):
                        tmp_x = '.'.join(tmp_list[0:len(tmp_list)-i-1])
                        if tmp_x in known_meme_tf_list:
                            show = x[0]+'('+x[1]+','+x[2]+')'
                            if show not in tf_select_shows:
                                tf_select_shows.append(show)
                                tf_select_values.append(x[0]+"|"+x[1])
                            break
                        else:
                            continue
            else:
                show = x[0]+'('+x[1]+','+x[2]+')'
                if show not in tf_select_shows:
                    tf_select_shows.append(show)
                    tf_select_values.append(x[0]+"|"+x[1])
        # dump
        self.dump_tool.add_tf_predict_detail(self.tf_predict.output_dir, self.option("main_id"))
        self.dump_tool.update_db_record('sg_tf_predict', self.option('main_id'),
                                        output_dir=workflow_output, status="end",
                                        tf_select_values=tf_select_values,
                                        tf_select_shows=tf_select_shows,)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tf_predict.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TfPredict"],
        ])
        super(TfPredictWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
