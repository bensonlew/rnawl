# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
__author__ = 'gdq'


class ExpVennAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(ExpVennAgent, self).__init__(parent)
        options = [
            # count 和 exp 均为定量结果，其中count用于差异分析。
            dict(name="express_matrix", type="infile", format="denovo_rna_v2.express_matrix"),
            # 没有重复实验时，可以用样本名作为组名。
            dict(name="group_table", type="infile", format="denovo_rna_v2.group_table"),
            dict(name="threshold", type="float", default=1),
        ]
        self.add_option(options)

    def check_options(self):
        try:
            threshold = float(self.option('threshold'))
        except:
            raise OptionError('无法将threshold的转换为数值，请检查输入', code = "32003101")
        else:
            if threshold <= 0:
                raise OptionError('threshold必须大于0', code = "32003102")

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

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
        threshold = float(self.option("threshold"))
        group_dict = self.option("group_table").prop['group_dict']
        exp_table = pd.read_table(self.option("express_matrix").prop['path'], index_col=0, header=0)
        all_sets = dict()
        f = open(os.path.join(self.output_dir, 'venn_graph.xls'), 'w')
        f.write('#set_name\tset_detail\n')
        for each_group in group_dict:
            samples = group_dict[each_group]
            expressed = exp_table.index[exp_table.loc[:, samples].mean(axis=1) >= threshold]
            all_sets[each_group] = list(expressed)
            f.write(each_group+'\t'+','.join(expressed)+'\n')
        f.close()
        # 计算交集
        # with open(os.path.join(self.output_dir, 'venn_table.xls'), 'w') as f:
        #     combination_num = len(all_sets)
        #     set_names = all_sets.keys()
        #     while combination_num >= 2:
        #         combinations = itertools.combinations(set_names, combination_num)
        #         for cmb in combinations:
        #             cmb_detail = [all_sets[x] for x in cmb]
        #             common = set(cmb_detail[0]).intersection(*cmb_detail)
        #             cmb_name = ' & '.join(cmb)
        #             f.write(cmb_name+'\t'+str(len(common))+'\t'+','.join(common)+'\n')
        #         combination_num -= 1
        #     #
        #     for each in set_names:
        #         each_detail = all_sets.pop(each)
        #         other = all_sets.values()
        #         each_only = set(each_detail).difference(*other)
        #         all_sets[each] = each_detail
        #         f.write(each+' only'+'\t'+str(len(each_only))+'\t'+','.join(each_only)+'\n')

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
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'
        data = {
            "id": "ExpVenn" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.exp_venn",
            "instant": True,
            "options": {
                "express_matrix": test_dir + '/' + 'transcript.tpm.matrix',
                "group_table": test_dir + '/' + 'default_group.txt',
                "threshold": 1,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

