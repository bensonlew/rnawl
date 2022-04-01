# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
import glob

class GenesetClusterWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetClusterWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="cluster_main_id", type="string"),
            dict(name="n_clusters", type='int'),
            dict(name="use_group", type="string"),
            dict(name="group", type="string"),
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
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("denovo_rna_v2.exp_cluster")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/01 Cluster_Analysis')
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
        super(GenesetClusterWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(GenesetClusterWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_geneset_cluster", main_id=self.option('cluster_main_id'), dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        group_dict = json.loads(self.option("group_dict"))
        cluster_exp = self.tool.output_dir + "/expression_matrix.xls"
        cluster_tree = self.tool.output_dir + "/seq.cluster_tree.txt"
        sample_tree = self.tool.output_dir + "/sample.cluster_tree.txt"
        subcluster_list = glob.glob(self.tool.output_dir + "/*subcluster_*.xls")
        chart.chart_geneset_cluster(cluster_exp, cluster_tree, sample_tree, subcluster_list, group_dict=group_dict)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        # add result info
        all_exp.add_geneset_cluster(self.tool.output_dir, main_id=self.option('cluster_main_id'), )
        self.end()

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.system('cp {} {}'.format(os.path.join(self.tool.output_dir,'run_parameter.txt'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.tool.output_dir,'expression_matrix.xls'),self.output_dir))
        os.system('find {} -maxdepth 1 -regextype posix-egrep -regex \'.*subcluster_.*.xls\' -exec cp -t {} {{}} +'.format(self.tool.output_dir, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/01 Cluster_Analysis", "", "聚类分析", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "聚类分析文件", 0],
            ["expression_matrix.xls", "xls", "聚类热图分析表",0,"201384"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ["./heatmap.pdf", 'pdf', "聚类热图", 0],
            ["./*heat_corr.pdf", 'pdf', "聚类热图", 0],
            ["./*.line.pdf", 'pdf', "子聚类折线图", 0],
        ])
        result_dir.add_regexp_rules([
            [r'.*subcluster_.*\.xls', 'xls', '子聚类分析表',0,"201385"],
        ])
        super(GenesetClusterWorkflow, self).end()

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
