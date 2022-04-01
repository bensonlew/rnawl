# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import os
import pandas as pd
import unittest
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete
from mbio.packages.ref_rna_v2.chart import Chart
import glob

class ExpCorrWorkflow(Workflow):
    """
    表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name='scm', type='string', default='complete'),
            dict(name='scd', type='string', default='correlation'),
            dict(name='corr_method', type='string', default='pearson'),
            dict(name="corr_main_id", type='string'),
            dict(name="log_base", type='int', default=None),
            dict(name="type", type='string'),
            dict(name="Draw_in_groups", type="string", default="no"),  # 20190513新增 是否以组别数据作图 by fwy
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("ref_rna_v2.exp_corr")
        if self.option("Draw_in_groups")=="yes":
            group_dict = json.loads(self.option("group_dict"))
            ori_exp=pd.read_table(self.option('exp_matrix'),index_col="seq_id",header=0)
            for i in group_dict.keys():
                t = i + "=" + "(" + group_dict[i][0]
                if len(group_dict[i]) > 1:
                    for j in group_dict[i][1:]:
                        t = t + "+" + j
                t = t + ")" + "/" + str(len(group_dict[i]))
                ori_exp.eval(t, inplace=True)
            final_exp= ori_exp.loc[:, group_dict.keys()]
            final_exp.to_csv(self.work_dir+"/work_exp_matrix", sep='\t')
            self.final_exp=self.work_dir+"/work_exp_matrix"
        else:
            self.final_exp=self.option('exp_matrix')
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Express/02 Exp_Corr')
        self.inter_dirs = []
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            interactiondelete =  InteractionDelete(bind_object=self,project_type="ref_rna_v2",main_id=self.option('corr_main_id'))
            interactiondelete.delete_interactions_records()


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
        super(ExpCorrWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(ExpCorrWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("ref_rna_v2.all_exp")
        # add result info
        all_exp.add_exp_corr2(self.tool.work_dir, main_id=self.option('corr_main_id'))
        self.end()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_exp_corr", main_id=self.option('corr_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()


    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        if self.option("Draw_in_groups")=="yes":
            group_dict = json.loads(self.option("group_dict"))
            group_dict2 = dict()
            for group in group_dict:
                group_dict2[group] = [group]
            group_dict = group_dict2
        else:
            group_dict = json.loads(self.option("group_dict"))
        exp_corr_file = self.tool.work_dir + "/sample_correlation.xls"
        exp_corr_tree_file = self.tool.work_dir + "/sample.cluster_tree.txt"

        chart.chart_exp_corr(exp_corr_tree_file, exp_corr_file, group_dict)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
        if os.path.exists(os.path.join(self.tool.output_dir,"sample_correlation.pdf")):
            os.remove(os.path.join(self.tool.output_dir,"sample_correlation.pdf"))
        os.link(pdf_file, self.tool.output_dir + "/sample_correlation.pdf")

    def end(self):
        self.chart()
        if os.path.exists(self.tool.output_dir + "/expression_matrix.xls"):
            os.remove(self.tool.output_dir + "/expression_matrix.xls")
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        # self.set_error("我就报错看看有没有问题")
        self.inter_dirs = [
            ["01 Express", "", "表达量分析结果目录",0],
            ["01 Express/02 Exp_Corr", "", "样本间相关性分析文件", 0]
        ]

        result_dir.add_relpath_rules([
            [".", "", "相关性分析结果目录", 0, "211064"],
            ['sample_correlation.xls', 'xls', '样本间相关性系数表',0,"211521"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*.pdf', 'pdf', '样本间相关性系数图', 0],
        ])
        super(ExpCorrWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.final_exp,
            scm=self.option('scm'),
            scd=self.option('scd'),
            corr_method=self.option('corr_method'),
        )
        if self.option('log_base'):
            options.update(log_base=self.option('log_base'))
        self.tool.set_options(options)
        self.tool.run()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.ref_rna_v2.report.exp_corr import ExpCorrWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'exp_corr_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.exp_corr',
            'options': {
                'exp_matrix' : '/mnt/ilustre/users/sanger-dev/workspace/20200519/Single_batch_limma6161_1693/Batch/output/count_batch.txt',
                # 'group_dict' : json.dumps('{"A1":["A1_1", "A1_2", "A1_3"],"A2": [ "A2_1", "A2_2", "A2_3"], "B1": [ "B1_1", "B1_2", "B1_3"], "B2": [ "B2_1", "B2_2", "B2_3"]}'.replace('"', '\\"')),
                'scd': 'euclidean'
            }
        }
        wsheet = Sheet(data=data)
        wf =ExpCorrWorkflow(wsheet)
        wf.sheet.id = 'exp_corr'
        wf.sheet.project_sn = 'exp_corr'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
