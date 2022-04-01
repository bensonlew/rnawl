# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import glob
import pandas as pd
from biocluster.file import getsize, exists
from biocluster.file import download
import os
from biocluster.config import Config
from bson.objectid import ObjectId
import unittest
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import glob
import re

class DiffexpWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name='express_id', type='string'),     # Obtain fill_type
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            dict(name='sig_type', type='string', default='pvalue'),
            dict(name="pvalue", type="float", default=0.05),
            # dict(name='padjust_way', type='int', default=3),
            dict(name="fc_down", type="float", default=0.83),
            dict(name="fc_up", type="float", default=1.2),
            dict(name="correct_method", type='string', default="two.sided"),
            dict(name="method_type", type="string"),
            dict(name="control_id", type="string"),
            dict(name="result_dir", type="string"),
            dict(name="protein_sliced", type="infile", format="labelfree.common"),
            dict(name="log", type='string', default="none"),
            dict(name="padjust_method", type='string', default='fdr'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("dia_v3.diff")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/4_ExpDiff')
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
        super(DiffexpWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(DiffexpWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        num_summary = os.path.join(self.tool.work_dir, "num_summary.xls")
        allsummary = os.path.join(self.tool.work_dir, "allsummary.xls")
        diff_files = glob.glob(os.path.join(self.tool.output_dir, "*_diff.xls"))
        diff_group = [os.path.basename(i).split('_diff')[0] for i in diff_files]
        if os.path.exists(num_summary) and os.path.exists(allsummary):
            chart.chart_diff(diff_group, num_summary, allsummary, diff_files)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(os.path.join(self.work_dir, "*.pdf"))
        for p in pdf_file:
            cmp_raw = os.path.splitext(os.path.basename(p))[0].split('_volcano.')
            if len(cmp_raw) > 1:
                p_type = 'volcano'
                cmp = cmp_raw[0]
            else:
                p_type = 'plot'
                cmp = cmp_raw[0].split('_scatter.')[0]
            if not os.path.exists(os.path.join(self.tool.output_dir, cmp)):
                os.makedirs(os.path.join(self.tool.output_dir, cmp))
            newfile = os.path.join(self.tool.output_dir, cmp, "Diff_" + self.option('method_type') + "_{}.pdf".format(p_type))
            if os.path.exists(newfile):
                os.remove(newfile)
            os.link(p, newfile)

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        diff = self.api.api("dia.diff")
        # add result info
        compare_dict_xls = self.tool.work_dir + "/compare_dict.xls"
        num_summary_xls = self.tool.work_dir + "/num_summary.xls"
        allsummary_xls = self.tool.work_dir + "/all_summary.xls"
        diff.add_diff(self.tool.work_dir, compare_dict_xls=compare_dict_xls,
        num_summary_xls=num_summary_xls, allsummary_xls=allsummary_xls, group_dict=self.option('group_dict'),
        main_id=self.option('diff_main_id'), method_type=self.option('method_type'), fill_type=self.fill_type)
        self.end()

    def end(self):
        # df_protein_sliced = pd.read_table(self.option("protein_sliced").prop['path'], header=0, sep="\t").round(6)
        # df_protein_sliced.rename(columns={df_protein_sliced.columns[0]: "accession_id"}, inplace=True)
        # df_protein_sliced.drop(["Description"], axis=1,inplace=True)
        # anno_file = self.download_s3_file(self.option("result_dir") + "anno_stat/proteins_anno_detail.xls", 'anno_file')
        # df_result_dir = pd.read_table(anno_file, header=0, sep="\t")
        # df_result_dir.rename(columns={df_result_dir.columns[0]: "accession_id"}, inplace=True)
        # files = glob.glob(self.tool.output_dir + "/" + "*_diff.xls*")
        # i = 0
        # for file in files:
        #     df = pd.read_table(file, header=0, sep="\t")
        #     df['compare'].replace(r"\|", "_vs_", inplace=True, regex=True)
        #     result = df.round(6)
        #     df_sum = pd.merge(pd.merge(result, df_result_dir, on='accession_id'), df_protein_sliced, on='accession_id')
        #     df_sum.to_csv(str(file).split('.xls')[0] + '.detail.xls', sep='\t', index=False)
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["4_ExpDiff", "", "差异蛋白数据挖掘结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "差异分析结果目录"],
        ])
        result_dir.add_regexp_rules([
            [r".*_vs_.*/.*plot.pdf", '', '差异散点图'],
            [r".*_vs_.*/.*volcano.pdf", '', '差异火山图']
        ])
        super(DiffexpWorkflow, self).end()

    def run_tool(self):
        options = dict(
            ratio_exp=self.option('exp_matrix'),
            group=self.option('group'),
            method_type=self.option('method_type'),
            cmp=self.option('cmp'),
            sig_type=self.option('sig_type'),
            pvalue=self.option('pvalue'),
            # padjust_way=self.option('padjust_way'),
            fc_down=self.option('fc_down'),
            fc_up=self.option('fc_up'),
            correct_method=self.option('correct_method'),
            group_dict=self.option('group_dict'),
            padjust_method=self.option('padjust_method'),
            log=self.option('log'),
        )
        db = Config().get_mongo_client(mtype='dia')[Config().get_mongo_dbname("dia")]
        exp_main = db["sg_express"]
        params = exp_main.find_one({'main_id': ObjectId(self.option('express_id'))})['params']
        if str(type(params)) == 'str':
            params_dict = json.loads(params)
            if 'fill_type' in params_dict.keys():
                self.fill_type = params_dict['fill_type']
                options['fill_type'] = self.fill_type
            else:
                self.fill_type = 'none'
        else:
            self.fill_type = 'none'
        # self.fill_type = 'group'
        # options['fill_type'] = self.fill_type
        self.tool.set_options(options)
        self.tool.run()

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.dia_v3.report.diffexp import DiffexpWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'dia_v3.report.diffexp',
            'options': {
                "exp_matrix": "/mnt/ilustre/users/sanger-dev/workspace/20201117/Dia_tsg_248820/treat_ref",
                "group": "/mnt/ilustre/users/sanger-dev/workspace/20201117/Dia_tsg_248820/remote_input/protein_group/group.txt",
                "cmp": "/mnt/ilustre/users/sanger-dev/workspace/20201117/Dia_tsg_248820/remote_input/protein_control/control.txt",
                "method_type": "student",
            }
        }
        wsheet_object = Sheet(data=data)
        wf = DiffexpWorkflow(wsheet_object)
        wf.sheet.id = 'dia'
        wf.sheet.project_sn = 'dia'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()



if __name__ == '__main__':
    unittest.main()
