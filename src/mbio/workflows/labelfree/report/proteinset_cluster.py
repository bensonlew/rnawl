# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
from mbio.packages.labelfree.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import os
import re
import glob
from biocluster.config import Config

class ProteinsetClusterWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetClusterWorkflow, self).__init__(wsheet_object)

        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="cluster_main_id", type="string"),
            dict(name="n_clusters", type='int'),
            dict(name="use_group", type="string"),
            dict(name="group", type="string"),
            dict(name="group_id", type="string"),
            dict(name="scm", type="string"),
            dict(name="scd", type="string"),
            dict(name="sct", type="string"),
            dict(name="gct", type="string"),
            dict(name="gcm", type="string"),
            dict(name="gcd", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.logger.info("参数为{}".format(self._sheet.options()))
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("labelfree.exp_cluster")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/01_Cluster')
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
        super(ProteinsetClusterWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(ProteinsetClusterWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("labelfree.all_exp")
        # add result info
        all_exp.add_geneset_cluster(self.tool.output_dir, main_id=self.option('cluster_main_id') )
        self.end()

    def chart(self):
        for (key, value) in [["LD_LIBRARY_PATH","$LD_LIBRARY_PATH:"+Config().SOFTWARE_DIR + "/bioinfo/sg_chart/miniconda2/lib"],["NODE_PATH",Config().SOFTWARE_DIR + "/bioinfo/sg_chart/node-v14.16.0-linux-x64/lib/node_modules"]]:
            if key not in os.environ.keys():
                os.environ[key] = value
            else:
                os.environ[key] = value + ":" + os.environ[key]

        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        expression_matrix = os.path.join(self.tool.output_dir, "expression_matrix.xls")
        seq_tree = os.path.join(self.tool.output_dir, "seq.cluster_tree.txt")
        sample_tree = os.path.join(self.tool.output_dir, "sample.cluster_tree.txt")
        chart.chart_proteinsetcluster(['cluster_webpage_'+str(self.option('n_clusters'))], [expression_matrix], [seq_tree], [sample_tree])
        chart.to_pdf()

        # move pdf to result dir
        for file in glob.glob(self.work_dir + "/*.pdf"):
            one_pdf = os.path.basename(file)
            if "____" in one_pdf:
                if "____cluster_webpage_" in one_pdf:
                    os.link(file, os.path.join(self.tool.output_dir, 'heat.pdf'))
                else:
                    os.link(file, os.path.join(self.tool.output_dir, one_pdf.split('.')[0].split('____')[1]+'.pdf'))

    def end(self):
        self.chart()
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析",0],
            ["5_Proteinset/01_Cluster", "", "表达模式聚类", 0]
        ]
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "蛋白集聚类分析结果目录"],
            ["heat.pdf", "", "热图", 0],
        ])
        result_dir.add_regexp_rules([
            [r"sub.*pdf", '', '子聚类趋势图'],
        ])
        super(ProteinsetClusterWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            group=self.option('group'),
            n_clusters=int(self.option('n_clusters')),
            sct=self.option('sct'),
            gct=self.option('gct'),
            scm=self.option('scm'),
            gcm=self.option('gcm'),
            scd=self.option('scd'),
            gcd=self.option('gcd'),
            use_group=self.option('use_group'),
        )
        self.tool.set_options(options)
        self.tool.run()
