# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import pandas as pd
import unittest
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import glob
from mbio.packages.denovo_rna_v2.chart import Chart


class ExpPcaWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpPcaWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="pca_main_id", type='string'),
            #dict(name="type", type='string'),
            dict(name="group_table",type="string"),  #20190705新增 group_table参数  by fwy
            dict(name="analysis_type",type="string",default="pca"), #20190705新增 group_tablec参数  by fwy
            dict(name="draw_in_groups",type="string",default="no"), #20190513新增 是否以组别数据作图 by fwy
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("denovo_rna_v2.exp_pca")
        self.ellipse = None
        if self.option("draw_in_groups") == "no":
            group_dict = json.loads(self.option("group_dict"))
            i = 0
            for sample_infos in group_dict.values():
                if len(sample_infos) < 3:
                    i += 1
                else:
                    continue
            if i == 0:
                self.ellipse = self.add_tool("graph.ellipse")
        if self.option("draw_in_groups") == "yes":
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
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/02 Express/03 Exp_PCA')
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
        super(ExpPcaWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        if not self.ellipse is None:
            self.tool.on("end", self.run_ellipse)
            self.ellipse.on("end", self.set_db)
            self.run_tool()
        else:
            self.tool.on("end", self.set_db)
            self.run_tool()
        super(ExpPcaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_exp_pca", main_id=self.option('pca_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_ellipse(self):
        options = {}
        if self.option("group_table"):
            options['group_table'] = self.option("group_table")
        pc_map = {'pca':"/PCA.xls",'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda':'/Dbrda/db_rda_sites.xls','nmds':'/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda'
                  ##'rda_cca'
                  }
        options['analysis'] = self.option('analysis_type')
        options['pc_table'] = self.tool.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(options)
        self.ellipse.run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        if self.option("draw_in_groups") == "yes":
            group_dict = json.loads(self.option("group_dict"))
            group_dict2 = dict()
            for group in group_dict:
                group_dict2[group] = [group]
            group_dict = group_dict2
        else:
            group_dict = json.loads(self.option("group_dict"))
        exp_pca_file = self.tool.work_dir + '/PCA.xls'
        exp_pca_var_file = self.tool.work_dir + '/Explained_variance_ratio.xls'
        if self.ellipse is not None:
            exp_pca_ellipse = self.ellipse.work_dir + '/ellipse_out.xls'
        else:
            exp_pca_ellipse = None
        chart.chart_exp_pca(exp_pca_file, exp_pca_var_file, group_dict=group_dict, exp_pca_ellipse=exp_pca_ellipse,
                            pcs=["PC1", "PC2"])
        chart.to_pdf()

        # move pdf to result dir
        pdf_files = glob.glob(self.work_dir + "/*.pdf")
        for file in pdf_files:
            if os.path.basename(file).endswith("ell.scatter.pdf"):
                os.link(file, self.tool.output_dir + "/sample_pca_ellipse.pdf")
            else:
                os.link(file, self.tool.output_dir + "/sample_pca.pdf")

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        # add result info
        all_exp.add_exp_pca2(self.tool.work_dir, main_id=self.option('pca_main_id'), )
        if not self.ellipse is None:
           all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(self.option('pca_main_id')))
        self.end()

    def end(self):
        self.chart()
        if os.path.exists(self.tool.output_dir + "/PCA.xls"):
            os.remove(self.tool.output_dir + "/PCA.xls")
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.system('cp {} {}'.format(os.path.join(self.tool.output_dir,'Explained_variance_ratio.xls'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.tool.output_dir,'run_parameter.txt'),self.output_dir))
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["02 Express", "", "表达量结果目录",0],
            ["02 Express/03 Exp_PCA", "", "样本间PCA分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "样本间PCA分析文件",0,"201374"],
            ['Explained_variance_ratio.xls', 'xls', '主成分解释表',0,"201375"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*pca.pdf', 'pdf', '样本间PCA图', 0],
            ['*pca_ellipse.pdf', 'pdf', '样本间PCA带置信圈图', 0],
        ])
        super(ExpPcaWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.final_exp,
        )
        self.tool.set_options(options)
        self.tool.run()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.denovo_rna_v2.report.exp_pca import ExpPcaWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "testExppca_denovo" + str(random.randint(1, 10000))+"yyyy",
            "type": "workflow",
            "name": "denovo_rna_v2.exp_pca",
            "instant": False,
            "options": dict(
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/denovo_rna_v2/test/Quant/gene.tpm.matrix",
                #group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                #corr_main_id="test_0711",
                #scm="complete",
                #scd="correlation",
                #corr_method="pearson",
                Draw_in_groups="no",
                #log_base="10",
            )
        }
        wsheet = Sheet(data=data)
        wf = ExpPcaWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    unittest.main()
