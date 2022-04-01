# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
__author__ = 'gdq'


class FusionVennAgent(Agent):
    """
    需输入表达量矩阵及样本分组信息
    """
    def __init__(self, parent):
        super(FusionVennAgent, self).__init__(parent)
        options = [
            # count 和 exp 均为定量结果，其中count用于差异分析。
            dict(name="fusion_matrix", type="infile", format="ref_rna_v2.common"),
            # 没有重复实验时，可以用样本名作为组名。
            dict(name="group_table", type="infile", format="denovo_rna_v2.group_table"),
            dict(name="filter_threshold", type="float", default=50),
            dict(name="use_group", type="string", default="no"),
        ]
        self.add_option(options)

    def check_options(self):
        try:
            use_group = self.option("use_group")
        except:
            raise OptionError('请选择是否以分组进行计算')
        if use_group == "yes":
            try:
                threshold = float(self.option('filter_threshold'))
            except:
                raise OptionError('无法将threshold的转换为数值，请检查输入', code="33704903")
            else:
                if threshold <= 0 or threshold > 100:
                    raise OptionError('threshold必须大于0且小于100')

    def set_resource(self):
        self._cpu = 1
        file_size = os.path.getsize(self.option('fusion_matrix').prop['path'])
        tmp = int(float(file_size)/1024**3)+5
        self._memory = '{}G'.format(tmp)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量venn图结果目录"],
        ])
        super(FusionVennAgent, self).end()


class FusionVennTool(Tool):
    def __init__(self, config):
        super(FusionVennTool, self).__init__(config)

    def get_group_info(self,group_dict):
        group_dict = group_dict
        group_num = {}
        for key in group_dict:
            group_num[key] = len(group_dict[key])
        sample_dict = {}
        for key in group_dict:
            for sample in group_dict[key]:
                sample_dict[sample] = key
        return group_dict,group_num,sample_dict

    def fusion_venn(self):
        group_dict,group_num,sample_dict = self.get_group_info(self.option("group_table").prop['group_dict'])
        if self.option("use_group") == "yes":
            threshold = float(self.option("filter_threshold")/100)
            fusion_df = pd.read_table(self.option("fusion_matrix").prop["path"])
            fusion_df["group"] = fusion_df["sample"].map(sample_dict)
            gp_fusion_df = fusion_df.groupby(["leftbreakpoint", "rightbreakpoint","group"])
            final_fusion_df = gp_fusion_df.filter(lambda x: x.shape[0] >= group_num[x["group"].values[0]]*threshold)
            group_details = {}
            for group in  group_dict:
                group_details[group] = final_fusion_df[final_fusion_df["group"] == group]["fusion_unique_id"].drop_duplicates().tolist()
            f = open(os.path.join(self.output_dir, 'fusion_venn_graph.xls'), 'w')
            f.write('name\tids\n')
            for each_group in group_details:
                ids = group_details[each_group]
                f.write(each_group + '\t' + ','.join(ids) + '\n')
        else:
            fusion_df = pd.read_table(self.option("fusion_matrix").prop["path"])
            sample_details = {}
            for sample in self.option("group_table").prop["sample"]:
                sample_details[sample] =  fusion_df[fusion_df["sample"] == sample]["fusion_unique_id"].drop_duplicates().tolist()
            f = open(os.path.join(self.output_dir, 'fusion_venn_graph.xls'), 'w')
            f.write('name\tids\n')
            for sample in self.option("group_table").prop["sample"]:
                ids = sample_details[sample]
                f.write(sample + '\t' + ','.join(ids) + '\n')

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
        super(FusionVennTool, self).run()
        self.fusion_venn()
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
            "id": "FusionVenn" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.exp_venn",
            "instant": True,
            "options": {
                "express_matrix": '/mnt/ilustre/users/sanger-dev/workspace/20190514/FusionVenn_tsg_33538_3815_1478/exp_matrix',
                "group_table": '/mnt/ilustre/users/sanger-dev/workspace/20190514/FusionVenn_tsg_33538_3815_1478/group',
                "threshold": 1,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

