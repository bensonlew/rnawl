# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.workflow import Workflow
import unittest
import os
from mbio.packages.dia_v3.chart import Chart
import glob
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json


class PreprocessWorkflow(Workflow):
    """
    dia v2.0 预处理
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PreprocessWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'raw_path', 'type': 'infile', 'format': 'labelfree.common'},
            {"name": "group_table", "type": "infile", "format": "labelfree.group_table"},
            {'name': 'searchdb', 'type': 'infile', 'format': 'labelfree.common'},
            {"name": "fillna", "type": "string", "default": "seqknn"},
            {"name": "update_info", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "fill_type", "type": "string", "default": "group"},
            {"name": "if_group", "type": "string", "default": "yes"},
            {"name": 'group_specific', 'type': "string", "default": "any"},
            {'name': 'group_percent', 'type': 'float', 'default': 50},
            {'name': 'all_eliminate', 'type': 'string', 'default': 'all'},
            # {'name': 'group_id', 'type': 'string'},
            # {'name': 'group_dict', 'type': 'string'},
            {'name': 'all_percent', 'type': 'float', 'default': 90},
            {"name": "raw_exp_id", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.preprocess = self.add_tool("dia_v3.preprocess")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/1_DataInfo/02_QcInfo')
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
        super(PreprocessWorkflow, self).send_log(data)

    def run(self):
        options = {
            "raw_path": self.option("raw_path"),
            "fillna": self.option("fillna"),
            "group_table": self.option("group_table"),
            "all_percent": self.option("all_percent"),
            "group_percent": self.option("group_percent"),
            "group_specific": self.option("group_specific"),
            "all_eliminate": self.option("all_eliminate"),
            "if_group": self.option("if_group"),
            "fill_type": self.option("fill_type"),
            'interactive': 'yes',
        }
        self.preprocess.set_options(options)
        self.preprocess.on('end', self.set_db)
        self.preprocess.run()
        super(PreprocessWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        preprocess_assess = os.path.join(self.preprocess.work_dir, '{}_nrmse.txt'.format(self.option('fillna')))
        preprocess_box = os.path.join(self.preprocess.work_dir, 'boxdata.xls')
        preprocess_summary = os.path.join(self.preprocess.work_dir, '{}_cv_summary.xls'.format(self.option('fillna')))
        if os.path.exists(preprocess_assess):
            chart.chart_preprocess_assessment(preprocess_assess)
        if os.path.exists(preprocess_box):
            chart.chart_preprocess_inter_cv(preprocess_box)
        if os.path.exists(preprocess_summary):
            chart.chart_preprocess_intra_cv(preprocess_summary)
        chart.to_pdf()

    def set_db(self):
        """
        导表程序
        """
        preprocess_api = self.api.api('dia.preprocess')
        exp_path = os.path.join(self.preprocess.work_dir, "{}_fillna.txt".format(self.option("fillna")))
        cv_path = os.path.join(self.preprocess.work_dir, "{}_cv_stats.xls".format(self.option("fillna")))
        nrmse_path = os.path.join(self.preprocess.work_dir, "{}_nrmse.txt".format(self.option("fillna")))
        cv_summary = os.path.join(self.preprocess.work_dir, "{}_cv_summary.xls".format(self.option("fillna")))
        preprocess_api.add_preprocess_detail(self.option("main_table_id"), exp_path, cv_path, nrmse_path=nrmse_path,
                                             cv_summary=cv_summary, searchdb=self.option('searchdb').prop['path'],
                                             interactive='yes')
        # preprocess_api.add_cv_box(cv_path, self.option("main_table_id"))
        # na_path = os.path.join(self.preprocess.work_dir, "{}_na_stats.xls".format(self.option("fillna")))
        # preprocess_api.add_preprocess_exp(exp_path, cv_path, nrmse_path, na_path=na_path, cv_summary=cv_summary,
        #                                   raw_exp_id=self.option("main_table_id"))

        self.end()

    def end(self):
        self.chart()

        # move pdf to result dir
        pdf_file = glob.glob(os.path.join(self.work_dir, "*.pdf"))
        for p in pdf_file:
            os.link(p, os.path.join(self.preprocess.output_dir, os.path.basename(p)))

        result_dir = self.add_upload_dir(self.preprocess.output_dir)
        self.inter_dirs = [
            ["1_DataInfo", "", "搜库数据结果目录", 0],
            ["1_DataInfo/02_QcInfo", "", "质控信息结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "预处理结果文件夹", 0],
            ["./*cv_stats.xls", "XLS", "蛋白表达量表", 0],
            ["./*cv_summary.xls", "XLS", "蛋白信息表", 0],
            ["./*fillna.txt", "TXT", "蛋白信息表", 0],
            ['./*nrmse.txt', 'TXT', "蛋白信息表", 0],
        ])
        result_dir.add_regexp_rules([
            [r".*_assessment.*pdf", '', '预处理效果图'],
            [r".*inter-cv.*pdf", '', '组间CV分布图'],
            [r".*intra-cv.*pdf", '', '组内CV分布图'],
        ])
        super(PreprocessWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.dia_v3.report.preprocess import PreprocessWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'exp_preprocess_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'dia_v3.report.preprocess',
            'options': {
                # "raw_path": "/mnt/ilustre/users/sanger-dev/workspace/20201126/Diav3_202011261345/raw_treat_ref",
                # "group_table": "/mnt/ilustre/users/sanger-dev/sg-users/xuxi/dia_v3/test_main_workflow/remote_input_from_majorbio_297008/protein_group/group.txt",
                "raw_path": "/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/exp.txt",
                "group_table": "/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/group(1).txt",
                "all_eliminate": "all",
                "all_percent": 90,
                "group_specific": "any",
                "fillna": "min",
                "if_group": "yes",
                "fill_type": "group",
                "main_table_id": "5f7763e317b2bf57d7ed373b",
                # "group_dict": r' ({"C_3": ["C_3_1", "C_3_2", "C_3_3", "C_3_4", "C_3_5"], '
                #               r'"R_3": ["R_3_1", "R_3_2", "R_3_3", "R_3_4", "R_3_5"]}'.replace('"', '\\"'),
            }
        }
        wsheet_object = Sheet(data=data)
        wf = PreprocessWorkflow(wsheet_object)
        wf.sheet.id = 'dia_test'
        wf.sheet.project_sn = 'dia_test'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    unittest.main()
