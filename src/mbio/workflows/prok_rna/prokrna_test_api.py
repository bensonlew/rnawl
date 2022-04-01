# -*- coding: utf-8 -*-

"""空workflow用于测试导表"""

from biocluster.workflow import Workflow
import os
import json

class ProkrnaTestApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProkrnaTestApiWorkflow, self).__init__(wsheet_object)
        options = []
        self.task_id = self.sheet.id
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.config.DBVersion=1

    def check_options(self):
        """
        检查参数设置
        """

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()


    def export_gensets_analysis(self):
        self.export_temporary = os.path.join(self.work_dir, "temporary")
        api = self.api.api('prok_rna.diff_geneset_work_pipline')
        diff_geneset_pipline_result = "/mnt/lustre/users/sanger-dev/wpm2/workspace/20210916/prok_pipline_8/DiffGenesetAllPipline/output"
        diff_id = None
        task_id = "f0kc_qn8lr0i6f23ej4hglqbnib"
        analysis_names = ["kegg", "go", "cog"]
        file_json_path = "/mnt/lustre/users/sanger-dev/wpm2/workspace/20210916/prok_pipline_8/DiffGenesetAllPipline/DiffGenesetAll/output/prepare_json"
        with open(file_json_path, "r") as j:
            file_dict = json.load(j)
        kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
        # self.exp_id_T = ""
        # self.long_group_id = '6105eb2c4f670a777ba7c681'
        exp_id = "613a4722dab19089c4541471"
        group_id = '613a46fddab19089c451d835'
        geneset_file = file_dict["genesets"]["All_Diff_mRNA"]['file_path']["kegg_class"]['multi_gene_list_path']
        # all_annot_path = os.path.join(self.annot_merge.output_dir, "allannot_class", "all_annot_gene.xls")
        # all_annot_df = pd.read_table(all_annot_path)
        # annot_df = all_annot_df[["transcript_id", "gene_name", "description"]]
        # annot_df.to_csv(os.path.join(self.export_temporary, "gene_detail"), sep="\t", index=False)
        api.add_diff_genest_pipline_table(diff_geneset_pipline_result, diff_id=diff_id, task_id=task_id,
                                          analysis_names=analysis_names,
                                          kegg_level_path=kegg_level_path, inter_path=self.export_temporary,
                                          group_id=group_id,geneset_file=geneset_file,exp_id=exp_id)
        self.end()


    def set_output(self, event):
        obj = event["bind_object"]

    def run(self):
        self.export_gensets_analysis()
        super(ProkrnaTestApiWorkflow, self).run()

    def end(self):
        super(ProkrnaTestApiWorkflow, self).end()
