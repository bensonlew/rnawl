# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from decimal import *
__author__ = 'gdq'


class ExpVennAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(ExpVennAgent, self).__init__(parent)
        options = [
            # count 和 exp 均为定量结果，其中count用于差异分析。
            dict(name="express_matrix", type="infile", format="labelfree.express_matrix"),
            # 没有重复实验时，可以用样本名作为组名。
            dict(name="group_table", type="infile", format="labelfree.group_table"),
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        file_size = os.path.getsize(self.option('express_matrix').prop['path'])
        tmp = int(float(file_size)/1024**3)+5
        self._memory = '{}G'.format(tmp)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量venn图结果目录"],
        ])
        super(ExpVennAgent, self).end()


class ExpVennTool(Tool):
    def __init__(self, config):
        super(ExpVennTool, self).__init__(config)

    def express_venn(self):
        group_dict = self.option("group_table").prop['group_dict']
        exp_table = pd.read_table(self.option("express_matrix").prop['path'], index_col=0, header=0)
        all_sets = dict()
        f = open(os.path.join(self.output_dir, 'venn_graph.xls'), 'w')
        f.write('#set_name\tset_detail\n')
        def determine_venn(row):
            count = 0
            for i in row[::]:
                # if i == "-":
                #     continue
                if Decimal(float(i)) == 0:
                    count += 1
            return count
        for each_group in group_dict:
            samples = group_dict[each_group]
            expressed = exp_table.index[exp_table.loc[:, samples].apply(determine_venn, axis=1) < len(samples)/2.0]
            all_sets[each_group] = list(expressed)
            f.write(each_group+'\t'+','.join(expressed)+'\n')
        f.close()

    def set_output(self):
        # all ready write results to output
        pass

    def run(self):
        super(ExpVennTool, self).run()
        self.express_venn()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/litangjian/labelfree_dev'
        data = {
            "id": "ExpVenn" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "labelfree.exp_venn",
            "instant": True,
            "options": {
                "express_matrix": test_dir + '/' + 'exp1.txt',
                "group_table": test_dir + '/' + 'group.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

