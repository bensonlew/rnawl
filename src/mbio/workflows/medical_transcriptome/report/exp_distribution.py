# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import unittest
import glob
import os
import re
from biocluster.core.function import CJsonEncoder
from shutil import copyfile
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart import Chart
from biocluster.config import Config
from bson import ObjectId


class ExpDistributionWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpDistributionWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name="graph_main_id", type="string"),
            dict(name="group_dict", type="string"),
            dict(name="exp_matrix", type="string"),
            dict(name="kind", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/06 Basic_Analysis/01 Exp_Distribution')
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
        super(ExpDistributionWorkflow, self).send_log(data)

    def run(self):
        # 没有tool时，意味着不需要self.run, 进一步意味着需要发起监听。
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        self.chart()
        self.set_db()
        self.end()
        # super(ExpDistributionWorkflow, self).run() 不需要

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        group_dict = json.loads(self.option("group_dict"))
        exp = self.option("exp_matrix")

        # 获取exp_type
        exp_type = "TPM"
        print "exp_type", exp_type
        try:
            db = Config().get_mongo_client(mtype="medical_transcriptome")[Config().get_mongo_dbname("medical_transcriptome")]
            # print "db is ", db
            params = db["sg_exp_graph"].find_one({"_id": ObjectId(self.option("graph_main_id"))})["params"]
            params = json.loads(params)
            # print params
            exp_info = db["sg_exp"].find_one({"_id": ObjectId(params["exp_id"])})

            # print exp_info
            exp_type = exp_info["exp_type"]
        except Exception as e:
            print e
            pass

        # 绘制分组均值分布图
        chart.chart_exp_dis_one(exp, group_dict, exp_type)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            if os.path.exists(self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "group_distribution")):
                os.remove(self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "group_distribution"))
            copyfile(p, self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "group_distribution"))
        # 绘制样本分布图
        
        chart.chart_exp_dis_one(exp, group_dict=None, exp_type=exp_type)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            if os.path.exists(self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution")):
                os.remove(self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))
            copyfile(p, self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_exp_graph", main_id=self.option('graph_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # time.sleep(1)
        all_exp = self.api.api("medical_transcriptome.all_exp")
        all_exp.add_distribution(
            self.option('exp_matrix'),
            main_id=self.option('graph_main_id'),
            group_dict=json.loads(self.option("group_dict")),
        )

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["06 Basic_Analysis", "", "表达量分析结果目录", 0],
            ["06 Basic_Analysis/01 Exp_Distribution", "", "表达量分布", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "表达量分布分析结果目录"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*box.pdf', 'pdf', '表达量分布盒形图', 0],
            ['*density.pdf', 'pdf', '表达量分布密度图', 0],
            ['*violin.pdf', 'pdf', '表达量分布小提琴图', 0],
        ])
        super(ExpDistributionWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.exp_distribution import ExpDistributionWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'exp_distribution_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'medical_transcriptome.report.exp_distribution',
            'options': {
                'exp_matrix' : '/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/Quant/output/gene.tpm.matrix',
                'group_dict' : json.dumps({"H1581":["H1581_1","H1581_2","H1581_3","H1581_4","H1581_5","H1581_6","H1581_7","H1581_8","H1581_9"],"SNU16":["SNU16_1","SNU16_2","SNU16_3","SNU16_4","SNU16_5","SNU16_6","SNU16_7","SNU16_8","SNU16_9"]}),
                'graph_main_id' : '5f3a584617b2bf32602a4b1f',
            }
        }
        wsheet = Sheet(data=data)
        wf =ExpDistributionWorkflow(wsheet)
        wf.sheet.id = 'exp_distribution'
        wf.sheet.project_sn = 'exp_distribution'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)