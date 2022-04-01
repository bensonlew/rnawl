# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import pandas as pd
import math
import time
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import os,re
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_advance import ChartAdvance
import glob

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
        # options = [
        #     {"name": "geneset_id", "type": "string"},
        #     {'name': 'families', 'type': 'string'},
        #     {'name': 'tf_predict_id', 'type': 'string'},
        #     {'name': 'e_value', 'type': 'float'},
        #     {'name': 'update_info', 'type': 'string'},
        #     {'name': "main_id", 'type': 'string'},
        # ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 Advanced_Analysis/05 TF_Analysis')
        self.inter_dirs = []
        self.dump_tool = self.api.api("medical_transcriptome.tf_stat")

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(TfStatWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        self.run_tf_stat()
        self.set_db()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_tf_stat", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

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
            geneset = geneset_records['seq_list']
            with open(self.work_dir + "/geneset.list", 'w') as f:
                for each in geneset:
                    f.write(each+'\n')
            conn = self.dump_tool.db['sg_geneset']
            result = conn.find_one({"main_id": ObjectId(options['geneset_id'])})
            geneset_type = result['level']
            print("geneset type: ", geneset_type)
        else:
            geneset = 'all'
            geneset_type = "all"
        # filter tf_predict_result
        conn = self.dump_tool.db['sg_tf_predict_detail']
        predict_result = conn.find({"tf_predict_id": ObjectId(options['tf_predict_id'])}, {'_id': 0, 'tf_predict_id': 0})
        print('tf_predict_id: ', options['tf_predict_id'])
        new_result = list()
        if options['families'] == "all" or (not options['families'].strip()):
            families = None
        else:
            families = options['families'].strip().split(",")
        if geneset == "all":
            for each in predict_result:
                if 'e_value' not in each:
                    each['e_value'] = 0
                each['e_value'] = float(each['e_value'])
                if families is None:
                    if each['e_value'] <= options['e_value']:
                        new_result.append(each)
                else:
                    try:
                        if each['e_value'] <= options['e_value'] and each['family'] in families:
                            new_result.append(each)
                    except:
                        continue
        else:
            for each in predict_result:
                if 'e_value' not in each:
                    each['e_value'] = 0
                    continue
                each['e_value'] = float(each['e_value'])
                # self.logger.info("e_value is {}, gene_id is{}, geneset_type is {}, families is {}, geneset is {}".format(each['e_value'], each['gene_id'], geneset_type, families, geneset))
                if families is None:
                    if geneset_type == "G":
                        if each['gene_id'] in geneset and each['e_value'] <= options['e_value']:
                            new_result.append(each)
                    else:
                        if each['transcript_id'] in geneset and each['e_value'] <= options['e_value']:
                            new_result.append(each)
                else:
                    if geneset_type == "G":
                        try:
                            if each['gene_id'] in geneset and each['e_value'] <= options['e_value'] and each['family'] in families:
                                new_result.append(each)
                        except:
                            print each
                    else:
                        try:
                            if each['transcript_id'] in geneset and each['e_value'] <= options['e_value'] and each['family'] in families:
                                new_result.append(each)
                        except:
                            print each
        # stat
        if not new_result:
            print("基因集中的基因均不是转录因子")
            return
        new_result_pd = pd.DataFrame(new_result)
        geneid2evalue = dict(zip(new_result_pd['gene_id'], new_result_pd['e_value']))
        transid2evalue = dict(zip(new_result_pd['transcript_id'], new_result_pd['e_value']))
        gene_result = new_result_pd.loc[:, ['gene_id', 'family']].drop_duplicates()
        trans_result = new_result_pd.loc[:, ['transcript_id', "family"]].drop_duplicates()
        gene_dict = dict()
        trans_dict = dict()
        for i in gene_result.index:
            gene_id = gene_result.loc[i]['gene_id']
            family = gene_result.loc[i]['family']
            if family not in gene_dict:
                gene_dict[family] = [gene_id]
            else:
                gene_dict[family].append(gene_id)
        for i in trans_result.index:
            trans_id = trans_result.loc[i]['transcript_id']
            family = trans_result.loc[i]['family']
            if family not in trans_dict:
                trans_dict[family] = [trans_id]
            else:
                trans_dict[family].append(trans_id)
        num_stat_pd = gene_result.groupby('family').count().join(trans_result.groupby('family').count())
        gene_list = list()
        trans_list = list()
        for i in num_stat_pd.index:
            gene_list.append(';'.join(gene_dict[i]))
            trans_list.append(';'.join(trans_dict[i]))
        num_stat_pd.columns = ['gene_num', 'transcript_num']
        num_stat_pd['genes_str'] = gene_list
        num_stat_pd['trans_str'] = trans_list
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
                tmp_dict = dict(gene_id=gene, type='G', e_value=-math.log10(float(geneid2evalue[gene])))
            else:
                tmp_dict = dict(gene_id=gene, type='G', e_value=0)
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
                tmp_dict = dict(transcript_id=gene, type="T", e_value=-math.log10(float(transid2evalue[gene])))
            else:
                tmp_dict = dict(transcript_id=gene, type="T", e_value=0)
            for each in families:
                if each == fam:
                    tmp_dict[each] = 1
                else:
                    tmp_dict[each] = 0
            tmp_dict_list.append(tmp_dict)
        pd.DataFrame(tmp_dict_list).set_index('transcript_id').to_csv(self.output_dir + '/transcript_circos_stat.txt', header=True, index=True, sep='\t')

    def chart(self):
        chart = ChartAdvance()
        chart.work_dir = self.work_dir + "/"

        chart.chart_tf_stat(self.output_dir + "/all_bar_stat.txt", self.output_dir + "/gene_circos_stat.txt", self.output_dir + "/transcript_circos_stat.txt")
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))


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
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 Advanced_Analysis", "", "高级分析结果目录",0],
            ["05 Advanced_Analysis/05 TF_Analysis", "", "转录因子分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "转录因子家族统计文件", 0],
            ["./all_bar_stat.txt", "", "转录因子家族统计表", 0],
            ["./transcript_circos_stat.txt", "", "转录因子家族circos转录本统计表", 0],
            ["./gene_circos_stat.txt", "", "转录因子家族circos基因统计表", 0],
            ["./filtered_tf_predict.xls", "", "转录因子家族详情表", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*.pdf', 'pdf', "转录因子家族统计图", 0]
        ])
        super(TfStatWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.tf_stat import TfStatWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "Diff_ASprofile" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.report.tf_stat",
            "options": dict(
                geneset_id="all",
                families="all",
                tf_predict_id="5f4e145217b2bf2095bdd6a4",
                e_value="0.00001",
            )
        }

        wsheet = Sheet(data=data)
        wf =TfStatWorkflow(wsheet)
        wf.sheet.id = 'tf_stat'
        wf.sheet.project_sn = 'tf_stat'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
