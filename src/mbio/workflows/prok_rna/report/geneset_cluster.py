# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
from mbio.packages.prok_rna.chart import Chart
import os
import glob
from collections import OrderedDict


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
            dict(name="group_id", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.logger.info("参数为{}".format(self._sheet.options()))
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("prok_rna.exp_cluster")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(GenesetClusterWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("prok_rna.all_exp")
        # add result info
        all_exp.add_geneset_cluster(self.tool.output_dir, main_id=self.option('cluster_main_id') )
        self.end()

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        subcluster_list = glob.glob(os.path.join(self.tool.work_dir, 'seq.subcluster*xls'))
        cluster_exp = os.path.join(self.tool.work_dir, 'expression_matrix.xls')
        cluster_tree = os.path.join(self.tool.work_dir, 'seq.cluster_tree.txt')
        sample_tree = os.path.join(self.tool.work_dir, 'sample.cluster_tree.txt')
        if self.option('group_dict'):
            group_dict = json.loads(self.option('group_dict'), object_pairs_hook=OrderedDict)
        elif self.option('group'):
            group_dict = json.loads(self.option('group'), object_pairs_hook=OrderedDict)
        samples = list()
        for each in group_dict.keys():
            samples.extend(group_dict[each])
        chart.chart_geneset_cluster(cluster_exp, cluster_tree, sample_tree, subcluster_list, group_dict=group_dict, samples_order=samples)
        chart.to_pdf()

        # move pdf
        target_dir = self.tool.output_dir
        cluster_pdfs = glob.glob(os.path.join(self.work_dir, 'subcluster*pdf'))
        for each in cluster_pdfs:
            num = os.path.basename(each).split('.')[1]
            new_path = os.path.join(target_dir, 'seq.subcluster_{}.pdf'.format(num))
            self.move_pdf(each, new_path)
        pdf = os.path.join(self.work_dir, 'geneset.cluster.heat_corr.pdf')
        self.move_pdf(pdf, os.path.join(target_dir, 'expression_matrix.pdf'))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_regexp_rules([
            [".", "", "基因集聚类分析结果目录"],
            ["expression_matrix.xls", "", "聚类分析结果表"],
            ["sample.cluster_tree.txt", "", "样本间聚类树"],
            ["expression_matrix.pdf", "", "基因聚类热图"],

        ])
        result_dir.add_regexp_rules([
            [r"seq\.subcluster_.*\.xls", "", "子聚类分析表"],
            [r"seq\.subcluster_.*\.pdf", "", "子聚类趋势图"],
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
