# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import pandas as pd
import math
import time


class TfStatWorkflow(Workflow):
    """
    TfStat description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TfStatWorkflow, self).__init__(wsheet_object)
        options = list()
        for each in self._sheet.options():
            options.append(dict(name=each, type="string"))
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.dump_tool = self.api.api("ref_rna.tf_stat")

    def run(self):
        self.start_listener()
        self.fire("start")
        self.run_tf_stat()
        self.set_db()

    def run_tf_stat(self):
        time.sleep(7)
        options = dict(
            geneset_id=self.option('geneset_id'),
            families=self.option('families'),
            tf_predict_id=self.option('tf_predict_id'),
            e_value=float(self.option('e_value')),
        )
        # get geneset
        if "all" not in options['geneset_id'].lower():
            conn = self.dump_tool.db['sg_geneset_detail']
            geneset_records = conn.find_one({"geneset_id": ObjectId(options['geneset_id'])})
            geneset = geneset_records['gene_list']
            with open(self.work_dir + "/geneset.list", 'w') as f:
                for each in geneset:
                    f.write(each+'\n')
            conn = self.dump_tool.db['sg_geneset']
            result = conn.find_one({"_id": ObjectId(options['geneset_id'])})
            geneset_type = result['type']
            print("geneset type: ", geneset_type)
        else:
            geneset = 'all'
            geneset_type = "all"
        # filter tf_predict_result
        conn = self.dump_tool.db['sg_tf_predict_detail']
        predict_result = conn.find({"tf_predict_id": ObjectId(options['tf_predict_id'])}, {'_id': 0, 'tf_predict_id': 0})
        new_result = list()
        if options['families'] == "all" or (not options['families'].strip()):
            families = None
        else:
            families = options['families'].strip().split(",")
        if geneset == "all":
            for each in predict_result:
                each['e_value'] = float(each['e_value'])
                if families is None:
                    if each['e_value'] <= options['e_value']:
                        new_result.append(each)
                else:
                    if each['e_value'] <= options['e_value'] and each['family'] in families:
                        new_result.append(each)
        else:
            for each in predict_result:
                each['e_value'] = float(each['e_value'])
                if families is None:
                    if geneset_type == "gene":
                        if each['gene_id'] in geneset and each['e_value'] <= options['e_value']:
                            new_result.append(each)
                    else:
                        if each['transcript_id'] in geneset and each['e_value'] <= options['e_value']:
                            new_result.append(each)
                else:
                    if geneset_type == "gene":
                        if each['gene_id'] in geneset and each['e_value'] <= options['e_value'] and each['family'] in families:
                            new_result.append(each)
                    else:
                        if each['transcript_id'] in geneset and each['e_value'] <= options['e_value'] and each['family'] in families:
                            new_result.append(each)
        # stat
        if not new_result:
            print("基因集中的基因均不是转录因子")
            return
        new_result_pd = pd.DataFrame(new_result)
        geneid2evalue = dict(zip(new_result_pd['gene_id'], new_result_pd['e_value']))
        transid2evalue = dict(zip(new_result_pd['transcript_id'], new_result_pd['e_value']))
        gene_result = new_result_pd.loc[:, ['gene_id', 'family']].drop_duplicates()
        trans_result = new_result_pd.loc[:, ['transcript_id', "family"]].drop_duplicates()
        num_stat_pd = gene_result.groupby('family').count().join(trans_result.groupby('family').count())
        num_stat_pd.columns = ['gene_num', 'transcript_num']
        num_stat_pd.to_csv(self.output_dir + "/all_bar_stat.txt", header=True, index=True, sep='\t')
        column_order = ['transcript_id', 'gene_id', 'gene_name', 'gene_desc', 'family', 'domain', 'domain_link', 'description', 'e_value', 'score']
        column_order += ['blast_hit', 'hit_link', 'hit_family', 'hit_pident', 'hit_evalue']
        new_result_pd = new_result_pd.loc[:, column_order]
        new_result_pd.to_csv(self.output_dir + '/filtered_tf_predict.xls', header=True, index=False, sep='\t')
        # stat2
        if families is None:
            families = set(gene_result['family'])
        tmp_dict_list = list()
        for gene, fam in zip(gene_result['gene_id'], gene_result['family']):
            if float(geneid2evalue[gene]) > 0:
                tmp_dict = dict(gene_id=gene, type='gene', e_value=-math.log10(float(geneid2evalue[gene])))
            else:
                tmp_dict = dict(gene_id=gene, type='gene', e_value=0)
            for each in families:
                if each == fam:
                    tmp_dict[each] = 1
                else:
                    tmp_dict[each] = 0
            tmp_dict_list.append(tmp_dict)
        pd.DataFrame(tmp_dict_list).set_index('gene_id').to_csv(self.output_dir + '/gene_circos_stat.txt', header=True, index=True, sep='\t')
        # stat3
        if families is None:
            families = set(trans_result['family'])
        tmp_dict_list = list()
        for gene, fam in zip(trans_result['transcript_id'], trans_result['family']):
            if float(transid2evalue[gene]) > 0:
                tmp_dict = dict(transcript_id=gene, type="transcript", e_value=-math.log10(float(transid2evalue[gene])))
            else:
                tmp_dict = dict(transcript_id=gene, type="transcript", e_value=0)
            for each in families:
                if each == fam:
                    tmp_dict[each] = 1
                else:
                    tmp_dict[each] = 0
            tmp_dict_list.append(tmp_dict)
        pd.DataFrame(tmp_dict_list).set_index('transcript_id').to_csv(self.output_dir + '/transcript_circos_stat.txt', header=True, index=True, sep='\t')

    def set_db(self):
        """
        dump data to db
        """
        workflow_output = self.get_workflow_output_dir()
        # add result info
        self.dump_tool.add_tf_stat_detail(self.output_dir, self.option("main_id"))
        self.dump_tool.update_db_record('sg_tf_stat', self.option('main_id'), output_dir=workflow_output, status="end")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TfStat"],
        ])
        super(TfStatWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
