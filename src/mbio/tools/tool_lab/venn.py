# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
__author__ = 'zoujiaxun'


class VennAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(VennAgent, self).__init__(parent)
        options = [
            # count 和 exp 均为定量结果，其中count用于差异分析。
            dict(name="venn_file", type="infile", format="ref_rna_v2.common"),
            dict(name='sep', type='string')
            # 没有重复实验时，可以用样本名作为组名。
            # dict(name="group_table", type="infile", format="denovo_rna_v2.group_table"),
            # dict(name="threshold", type="float", default=1),
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('venn_file').is_set:
            raise OptionError('必须输入venn分析文件')
        # if not self.option('sep').is_set:
        #     raise OptionError('必须输入文件分隔符类型')
        return True

    def set_resource(self):
        self._cpu = 1
        file_size = os.path.getsize(self.option('venn_file').prop['path'])
        tmp = int(float(file_size)/1024**3)+5
        self._memory = '{}G'.format(tmp)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "venn图结果目录"],
        ])
        super(VennAgent, self).end()


class VennTool(Tool):
    def __init__(self, config):
        super(VennTool, self).__init__(config)

    def express_venn(self):
        # threshold = float(self.option("threshold"))
        # group_dict = self.option("group_table").prop['group_dict']
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[self.option('sep')]
        exp_table = pd.read_table(self.option("venn_file").prop['path'], index_col=0, header=0, sep=sep)
        all_sets = dict()
        f = open(os.path.join(self.output_dir, 'venn_graph.xls'), 'w')
        f.write('#set_name\tset_detail\n')
        for i in exp_table:
            expressed = exp_table.index[exp_table.loc[:, i] > 0]
            all_sets[i] = list(expressed)
            f.write(i + '\t' + ','.join(expressed) + '\n')
        # for each_group in group_dict:
        #     samples = group_dict[each_group]
        #     expressed = exp_table.index[exp_table.loc[:, samples].mean(axis=1) >= threshold]
        #     all_sets[each_group] = list(expressed)
        #     f.write(each_group+'\t'+','.join(expressed)+'\n')
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
        super(VennTool, self).run()
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
        data = {
            "id": "Venn" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.venn",
            "instant": True,
            "options": {
                # "express_matrix": '/mnt/ilustre/users/sanger-dev/workspace/20190514/ExpVenn_tsg_33538_3815_1478/exp_matrix',
                # "group_table": '/mnt/ilustre/users/sanger-dev/workspace/20190514/ExpVenn_tsg_33538_3815_1478/group',
                # "threshold": 1,
                'venn_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/transcript.count.matrix',
                'sep': '\t'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

