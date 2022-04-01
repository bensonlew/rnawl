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

class HclusterModule(Module):
    """
    该module为了基因组坐标转化功能，针对没有已经chain文件的分析项目开发的module,其目的功能是通过给定的
    已知fasta和目标fasta，通过互相比对生成chain文件，后通过liftover进行分析，来生成最终的结果文件
    """
    def __init__(self, work_id):
        super(HclusterModule, self).__init__(work_id)
        options = [
            {"name": "method", "type": "string", "default": "bray_curtis"},
            {"name": "otutable", "type": "infile",
             "format": "meta.otu.otu_table"},
            {"name": "linkage", "type": 'string', "default": "average"},
            {"name": "sep", "type": 'string', "default": "tab"},

        ]
        self.add_option(options)



    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def run(self):
        super(HclusterModule, self).run()
        self.run_distance()

    def run_distance(self):
        # self.distance = self.add_tool('tool_lab.distance_calc')
        self.distance = self.add_tool('tool_lab.distance')
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[self.option('sep')]
        otutable = pd.read_table(self.option('otutable').path,sep=sep)
        otutable.to_csv(os.path.join(self.distance.work_dir, 'otutable.txt'), sep='\t',index=False,header=True)
        otutable_path = os.path.join(self.distance.work_dir, 'otutable.txt')
        self.distance.set_options({
            "method": self.option("method"),
            "otutable": otutable_path,
            "sep": self.option('sep')

        })
        self.distance.on('end', self.run_cluster)
        self.distance.run()

    def run_cluster(self):
        self.cluster = self.add_tool('tool_lab.hcluster')
        self.cluster.set_options({
            "dis_matrix": self.distance.option('dis_matrix'),
            "linkage": self.option('linkage')
        })
        self.cluster.on('end', self.set_output)
        self.cluster.run()

    def set_output(self):
        for file_name in os.listdir(self.cluster.output_dir):
            source = os.path.join(self.cluster.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        for table in os.listdir(self.distance.output_dir):
            tab_source = os.path.join(self.distance.output_dir, table)
            link_names = os.path.join(self.output_dir, table)
            if os.path.isfile(link_names):
                os.remove(link_names)
            os.link(tab_source, link_names)
        self.end()


    def end(self):
        super(HclusterModule, self).end()



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "cluster" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.hcluster",
            "instant": False,
            "options": dict(
                method="bray_curtis",
                otutable="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/cluster/otu_table.xls",
                linkage="average"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()