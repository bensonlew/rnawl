# -*- coding: utf-8 -*-
# __author__ = "chenyanyan, 2016.10.12"
# last_modify by khl 20170504
import os
from biocluster.workflow import Workflow
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import re
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
            dict(name="gene_detail", type="string"),
            dict(name="group_id", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("lnc_rna.geneset.exp_cluster")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/04 GeneSet/01 Cluster_Analysis')
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
        get_run_log = GetRunLog("lnc_rna", table="sg_geneset_cluster", main_id=self.option('cluster_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("lnc_rna.lnc_rna_geneset")
        seq_id_marker = os.path.join(os.path.dirname(self.option('exp_matrix')), 'seq_id_marker.json')
        # add result info
        # self, cluster_output_dir, gene_detail_file, main_id=None, project_sn='lnc_rna', task_id='lnc_rna', params=None
        all_exp.add_geneset_cluster(self.tool.output_dir, self.option("gene_detail"), marker=seq_id_marker,
                                    main_id=self.option('cluster_main_id'))
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tool.output_dir, os.path.basename(self.run_log)))

        rm_files = glob.glob(self.tool.output_dir + '/*cluster*tree*')
        rm_files += glob.glob(self.tool.output_dir + '/seq.*cluster.txt')
        for each in rm_files:
            os.remove(each)

        result_dir = self.add_upload_dir(self.tool.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/01 Cluster_Analysis", "", "聚类分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "基因集聚类分析结果目录", 0, "211072"],
            ['./expression_matrix.xls', 'xls', '聚类热图分析表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'seq.subcluster_.*\.xls', 'xls', '子聚类分析表', 0],
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
