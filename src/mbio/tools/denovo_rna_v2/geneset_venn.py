#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
__author__ = 'gdq'


class GenesetVennAgent(Agent):
    """
    输入gene_sets示例,一般由数据库导出：
    ------------------------------------
    set1 AT5G57625,AT1G05250,AT4G05110
    set1 AT5G57625,other,AT4G05110
    set1 AT5G57625,AT1G05250
    ........
    ------------------------------------
    """
    def __init__(self, parent):
        super(GenesetVennAgent, self).__init__(parent)
        options = [
            dict(name="gene_sets", type="string"),
            dict(name="graph", type="string"),
            dict(name="table", type="string"),
        ]
        self.add_option(options)

    def check_options(self):
        if not os.path.exists(self.option('gene_sets')):
            raise OptionError("Gene set file: %s not exist", variables=(self.option('gene_sets')), code="32004501")

    def set_resource(self):
        file_size = os.path.getsize(self.option('gene_sets'))
        self._cpu = 1
        self._memory = '{}G'.format(round(float(file_size)/1024**3)+0.5)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["./venn_graph.xls", "xls", "venn图输出目录1"],
            ["./venn_table.xls", "xls", "venn图输出目录2"]
        ])
        super(GenesetVennAgent, self).end()


class GenesetVennTool(Tool):
    def __init__(self, config):
        super(GenesetVennTool, self).__init__(config)

    def geneset_venn(self):
        geneset_info = pd.read_table(self.option("gene_sets"), index_col=0, header=None)
        geneset_info.index.name = "#set_name"
        geneset_info.columns = ['set_detail']
        geneset_info.to_csv(os.path.join(self.output_dir, 'venn_graph.xls'), sep='\t')
        with open(os.path.join(self.output_dir, 'venn_table.xls'), 'w') as f:
            combination_num = len(geneset_info.index)
            all_sets = {x: geneset_info.loc[x][0].split(',') for x in geneset_info.index}
            while combination_num >= 2:
                combinations = itertools.combinations(geneset_info.index, combination_num)
                for cmb in combinations:
                    cmb_detail = [all_sets[x] for x in cmb]
                    common = set(cmb_detail[0]).intersection(*cmb_detail)
                    cmb_name = ' & '.join(cmb)
                    f.write(cmb_name+'\t'+str(len(common))+'\t'+','.join(common)+'\n')
                combination_num -= 1
            #
            for each in all_sets:
                each_detail = all_sets.pop(each)
                other = all_sets.values()
                each_only = set(each_detail).difference(*other)
                all_sets[each] = each_detail
                f.write(each+' only'+'\t'+str(len(each_only))+'\t'+','.join(each_only)+'\n')

    def set_output(self):
        # all ready write results to output
        self.option('graph', os.path.join(self.output_dir, 'venn_graph.xls'))
        self.option('table', os.path.join(self.output_dir, 'venn_table.xls'))

    def run(self):
        super(GenesetVennTool, self).run()
        self.geneset_venn()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'
        data = {
            "id": "GenesetVenn" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.geneset_venn",
            "instant": True,
            "options": dict(
                gene_sets=test_dir + '/geneset_file_geneset_venn',
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
