# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""绘制群落组成的热图时需要hcluster"""

import os
import re
import datetime
import json
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.otu_table import OtuTableFile


class HeatClusterWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HeatClusterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "samples", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "linkage", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "newick_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.matrix = self.add_tool('meta.beta_diversity.distance_calc')
        self.hcluster = self.add_tool('meta.beta_diversity.hcluster')
        self.samples = re.split(',', self.option("samples"))

    def check_options(self):
        if self.option('linkage') not in ['average', 'single', 'complete']:
            raise OptionError('错误的层级聚类方式：%s', variables=(self.option('linkage')), code="12701901")

    def run_matrix(self, trans_otu):
        options = {
            "method": "bray_curtis",
            "otutable": trans_otu
        }
        self.matrix.set_options(options)
        self.matrix.run()

    def run_cluster(self):
        options = {
            "dis_matrix": self.matrix.option('dis_matrix'),
            "linkage": self.option('linkage')
        }
        self.logger.debug("cluster")
        self.hcluster.set_options(options)
        self.hcluster.on('end', self.set_db)
        self.hcluster.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        api_heat_cluster = self.api.heat_cluster
        myParams = json.loads(self.sheet.params)
        name = "heat_cluster_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
        newick_id = api_heat_cluster.create_newick_table(self.sheet.params, self.option("linkage"), myParams["otu_id"], name)
        self.hcluster.option("newicktree").get_info()
        tree_path = self.hcluster.option("newicktree").prop['path']
        api_heat_cluster.update_newick(tree_path, newick_id)
        self.add_return_mongo_id("sg_newick_tree", newick_id)
        self.end()

    def run(self):
        self.matrix.on('end', self.run_cluster)
        no_zero_otu = os.path.join(self.work_dir, "otu.nozero")
        my_sps = self.samples
        self.option("in_otu_table").sub_otu_sample(my_sps, no_zero_otu)
        no_zero_file = OtuTableFile()
        no_zero_file.set_path(no_zero_otu)
        trans_otu = os.path.join(self.work_dir, "otu.trans")
        no_zero_file.transposition(trans_otu)
        self.run_matrix(trans_otu)
        super(HeatClusterWorkflow, self).run()
