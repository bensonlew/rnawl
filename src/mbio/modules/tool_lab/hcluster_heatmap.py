# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200411

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import subprocess
from biocluster.config import Config
import pandas as pd

class HclusterHeatmapModule(Module):
    """
    该module为了基因组坐标转化功能，针对没有已经chain文件的分析项目开发的module,其目的功能是通过给定的
    已知fasta和目标fasta，通过互相比对生成chain文件，后通过liftover进行分析，来生成最终的结果文件
    """
    def __init__(self, work_id):
        super(HclusterHeatmapModule, self).__init__(work_id)
        options = [
            {"name":'otutable', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {"name": "sep", "type": 'string', "default": "tab"},
            {"name": 'log_change', "type": 'string'},
            {'name': 'scale', 'type': 'string'},
            {"name": 'cluster', 'type': 'bool'},
            {'name': 'gcd', 'type': 'string'},   #行距离算法
            {"name": 'gcm', 'type': 'string'},   #行聚类方式
            {"name": 'scd', 'type': 'string'},   #列距离算法
            {"name": 'scm', 'type': 'string'},   #列聚类方式
            {"name": 'title', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.distance_col = self.add_tool('tool_lab.distance')
        self.cluster_col = self.add_tool('tool_lab.hcluster')
        self.distance_row = self.add_tool('tool_lab.distance')
        self.cluster_row = self.add_tool('tool_lab.hcluster')
        self.heatmap = self.add_tool('tool_lab.heatmap')

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(HclusterHeatmapModule, self).run()
        # if self.option('cluster') == False:
        #     self.run_heatmap()
        #     self.on_rely([self.heatmap], self.set_output)
        # else:
        #     self.run_heatmap()
        #
        #     self.on_rely([self.heatmap], self.run_distance_col)
        #     self.on_rely([self.heatmap], self.run_distance_row)
        self.run_heatmap()
        self.on_rely([self.cluster_col, self.cluster_row, self.heatmap], self.set_output)
    def run_distance_col(self):
        sep_in = 'tab'
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[sep_in]
        otutable = pd.read_table(self.heatmap.option('scale_table').path,sep=sep)
        otutable.to_csv(os.path.join(self.distance_col.work_dir, 'sample_otutable.txt'), sep='\t',index=False,header=True)
        otutable_path = os.path.join(self.distance_col.work_dir, 'sample_otutable.txt')
        a = otutable < 0
        c = a.values.flatten().tolist()
        # a = otutable.values.flatten().tolist()
        # c = [b < 0 for b in a]
        if any(c) and self.option('scd').lower() in ['bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
                                                          'canberra', 'hellinger', 'kulczynski', 'abund_jaccard', 'morisita_horn',
                                                          'soergel', 'specprof', 'binary_sorensen_dice', 'binary_hamming',
                                                          'binary_jaccard', 'binary_lennon', 'binary_ochiai', 'chisq']:
            self.set_error("原始表中或者标准化之后的表中含有负值，无法使用{}距离算法，请重新选择参数".format(self.option('scd')))
            # raise ValueError("原始表中或者标准化之后的表中含有负值，无法使用{}距离算法，请重新选择参数".format(self.option('scd')))
        self.distance_col.set_options({
            "method": self.option("scd"),
            "otutable": otutable_path,
            'sep': sep_in

        })
        self.distance_col.on('end', self.run_cluster_col)
        self.distance_col.run()

    def run_cluster_col(self):

        self.cluster_col.set_options({
            "dis_matrix": self.distance_col.option('dis_matrix'),
            "linkage": self.option('scm').lower()
        })
        # self.cluster_col.on('end', self.set_output)
        self.cluster_col.run()

    def run_distance_row(self):
        sep_in = 'tab'
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[sep_in]
        otutable = pd.read_table(self.heatmap.option('scale_table').path, sep=sep)
        otutable = otutable.T
        otutable.to_csv(os.path.join(self.distance_row.work_dir, 'gene_otutable.txt'), sep='\t',header=False)
        otutable_path = os.path.join(self.distance_row.work_dir, 'gene_otutable.txt')
        self.distance_row.set_options({
            "method": self.option("gcd"),
            "otutable": otutable_path,
            'sep': sep_in

        })
        self.distance_row.on('end', self.run_cluster_row)
        self.distance_row.run()

    def run_cluster_row(self):

        self.cluster_row.set_options({
            "dis_matrix": self.distance_row.option('dis_matrix'),
            "linkage": self.option('gcm').lower()
        })
        # self.cluster_row.on('end', self.set_output)
        self.cluster_row.run()

    def run_heatmap(self):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[self.option('sep')]
        otutable = pd.read_table(self.option('otutable').path,sep=sep)
        otutable.to_csv(os.path.join(self.heatmap.work_dir, 'otutable.txt'), sep='\t',index=False,header=True)
        otutable_path = os.path.join(self.heatmap.work_dir, 'otutable.txt')
        self.heatmap.set_options({
            'otutable': otutable_path,
            'log_change': self.option('log_change'),
            'scale': self.option('scale')
        })
        if self.option('cluster') == False:
            self.heatmap.on('end', self.set_output)
        else:
            self.heatmap.on('end', self.run_distance_col)
            self.heatmap.on('end', self.run_distance_row)
        self.heatmap.run()
    def set_output(self):
        if self.option('cluster') == True:
            for file_name in os.listdir(self.cluster_col.output_dir):
                source = os.path.join(self.cluster_col.output_dir, file_name)
                link_name = os.path.join(self.output_dir, '{}_{}'.format('sample',file_name))
                if os.path.isfile(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
            for file_name in os.listdir(self.cluster_row.output_dir):
                source = os.path.join(self.cluster_row.output_dir, file_name)
                link_name = os.path.join(self.output_dir, '{}_{}'.format('feature', file_name))
                if os.path.isfile(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
            for file_name in os.listdir(self.heatmap.output_dir):
                source = os.path.join(self.heatmap.output_dir, file_name)
                link_name = os.path.join(self.output_dir, file_name)
                if os.path.isfile(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
        else:
            for file_name in os.listdir(self.heatmap.output_dir):
                source = os.path.join(self.heatmap.output_dir, file_name)
                link_name = os.path.join(self.output_dir, file_name)
                if os.path.isfile(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
        self.end()


    def end(self):
        super(HclusterHeatmapModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "cluster_heatmap" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.hcluster_heatmap",
            "instant": False,
            "options": dict(
                gcd="euclidean",
                otutable="/mnt/ilustre/users/sanger-dev/workspace/20200508/GenesetCluster_tsg_37259_1072_2943/exp_matrix",
                gcm="complete",
                scd="euclidean",
                scm="single",
                log_change='log10',
                scale='zscore',
                cluster=True
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()