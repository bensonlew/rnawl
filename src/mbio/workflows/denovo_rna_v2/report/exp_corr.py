# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import os
import pandas as pd
import unittest
from collections import OrderedDict
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.denovo_rna_v2.chart import Chart
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
            dict(name="log_base", type='int', default=None),                    #20190705新增 是否以取log进行分析 by fwy
            dict(name="draw_in_groups", type="string", default="no"),       # 20190705新增 是否以组别数据作图 by fwy
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("denovo_rna_v2.exp_corr")
        if self.option("draw_in_groups")=="yes":
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
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/02 Express/02 Exp_Corr')
        self.inter_dirs = []

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

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_exp_corr", main_id=self.option('corr_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        if self.option("draw_in_groups")=="yes":
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
        os.link(pdf_file, self.tool.output_dir + "/sample_correlation.pdf")

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        # add result info
        all_exp.add_exp_corr2(self.tool.work_dir, main_id=self.option('corr_main_id'), )
        self.end()

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.system('cp {} {}'.format(os.path.join(self.tool.output_dir,'sample_correlation.xls'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.tool.output_dir,'run_parameter.txt'),self.output_dir))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["02 Express", "", "表达量结果目录",0],
            ["02 Express/02 Exp_Corr", "", "样本间相关性分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "样本间相关性分析文件",0,"201369"],
            ['sample_correlation.xls', 'xls', '样本间相关性系数表',0,"201370"],
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
            log_base=self.option('log_base'),
        )
        self.tool.set_options(options)
        self.tool.run()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.denovo_rna_v2.report.exp_corr import ExpCorrWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "testExpcorr_denovo" + str(random.randint(1, 10000))+"yyyy",
            "type": "workflow",
            "name": "denovo_rna_v2.exp_corr",
            "instant": False,
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/denovo_rna_v2/test/Quant/gene.tpm.matrix",
                #group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                #corr_main_id="test_0711",
                scm="complete",
                scd="correlation",
                corr_method="pearson",
                draw_in_groups="yes",
                log_base="10",
            )
        }
        wsheet = Sheet(data=data)
        wf = ExpCorrWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    unittest.main()

