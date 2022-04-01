# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import unittest
import glob
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.medical_transcriptome.chart.chart import Chart


class ExpVennWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpVennWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="venn_main_id", type='string'),
            dict(name="threshold", type='string'),
            dict(name="group", type="string"),
            dict(name="type", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("medical_transcriptome.exp_venn")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/02 Annotation_Express/05 Exp_Venn')
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
        super(ExpVennWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(ExpVennWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        chart.chart_exp_venn(os.path.join(self.tool.output_dir, 'venn_graph.xls'))
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.tool.output_dir + "/" + os.path.basename(p))

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("medical_transcriptome.all_exp")
        # add result info
        graph_table = os.path.join(self.tool.output_dir, 'venn_graph.xls')
        all_exp.add_exp_venn(graph_table, main_id=self.option('venn_main_id'), )
        self.end()

    def end(self):
        self.chart()
        venn_pdf_file = glob.glob(self.work_dir + "/*venn.pdf")[0]
        upset_pdf_file = glob.glob(self.work_dir + "/*upset.pdf")[0]
        if os.path.exists(self.output_dir + "/sample_venn.pdf"):
            os.remove(self.output_dir + "/sample_venn.pdf")
        os.link(venn_pdf_file, self.output_dir + "/sample_venn.pdf")
        if os.path.exists(self.output_dir + "/sample_upset.pdf"):
            os.remove(self.output_dir + "/sample_upset.pdf")
        os.link(upset_pdf_file, self.output_dir + "/sample_upset.pdf")

        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["02 Annotation_Express", "", "所有基因数据挖掘结果目录", 0],
            ["02 Annotation_Express/05 Exp_Venn", "", "样本间Venn分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达量venn分析结果目录"],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*venn.pdf', 'txt', '样本间Venn图', 0],
            ['*upset.pdf', 'txt', '样本间Upset图', 0],
        ])

        super(ExpVennWorkflow, self).end()

    def run_tool(self):
        options = dict(
            express_matrix=self.option('exp_matrix'),
            group_table=self.option('group'),
            threshold=self.option('threshold'),
        )
        self.tool.set_options(options)
        self.tool.run()




class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.exp_venn import ExpVennWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'ExpVenn_workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'medical_transcriptome.report.exp_venn',
            'options': {
                'exp_matrix':'/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38312/Quant/gene.count.matrix',
                'threshold':"1",
                'group':'/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38312/remote_input/group_table/group.txt',
            }
        }
        wsheet_object = Sheet(data=data)
        wf = ExpVennWorkflow(wsheet_object)
        wf.sheet.id = 'ExpVennWorkflow'
        wf.sheet.project_sn = 'ExpVennWorkflow'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    unittest.main()
