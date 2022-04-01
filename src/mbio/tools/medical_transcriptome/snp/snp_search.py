# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
__author__ = 'fwy'


class SnpSearchAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(SnpSearchAgent, self).__init__(parent)
        options = [

            dict(name="target_genes", type="string"),
            dict(name='snp_result_file', type='infile', format='ref_rna_v2.common'),
            dict(name='type', type='string'),
            dict(name='sample', type='string'),
            dict(name='region', type='string'),
            dict(name='depth', type='int'),
            dict(name='depth_comare', type='string'),
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        pass

    def set_resource(self):
        self._cpu = 1

        self._memory = '20G'

    def end(self):
        super(SnpSearchAgent, self).end()


class SnpSearchTool(Tool):
    def __init__(self, config):
        super(SnpSearchTool, self).__init__(config)

    def run(self):
        super(SnpSearchTool, self).run()
        self.run_snp_search()
        self.end()

    def run_snp_search(self):
        snp_detail = pd.read_table(self.option("snp_result_file").prop["path"])
        sample = self.option("sample")
        depth = self.option("depth")
        depth_comare = self.option("depth_comare")
        target_genes =list()
        if not self.option("target_genes") == "all":
            with open(self.option("target_genes")) as t:
                for i in t.readlines():
                    target_genes.append(i.strip())
        else:
            target_genes = "all"
        #过滤条件较多,分步过滤
        #step1:按照区域和种类分类
        if self.option("region") == "all":
            search_df = snp_detail
        else:
            search_df = snp_detail[(snp_detail["type"] == self.option("type")) & (snp_detail["Anno"] == self.option("region") )]
        #step2 按照基因过滤,如果客户没选择基因过滤,则会传入一个all，否则将是一个文件。
        if target_genes ==  "all":
            pass
        else:
            search_df = search_df[search_df["GENE(in or nearby)"].isin(target_genes)]
        #step3 按照样本过滤
        #step4 按照depth过滤,因为选择all和单样本时逻辑原则不同因此合并两步分别讨论计算
        def get_all(x, depth, columns):
            numset = set()
            for i in columns:
                numset.add(x[i])
            if len(numset) == 1:
                for num in numset:
                    if num == depth:
                        return 1
                    else:
                        return 0
            else:
                return 0
        if self.option("sample") == "all":
            search_df = search_df
            sample_depths,samples,sample_t_depths =self.get_depth_columns(search_df)
            if depth_comare == "greater":
                search_df["min_depth"] = search_df.apply(lambda x: x[sample_t_depths].min(), axis=1)
                search_df = search_df[search_df["min_depth"] > depth]
                search_df = search_df.drop("min_depth", axis=1)
            elif depth_comare == "greateroreq":
                search_df["min_depth"] = search_df.apply(lambda x: x[sample_t_depths].min(), axis=1)
                search_df = search_df[search_df["min_depth"] >= depth]
                search_df = search_df.drop("min_depth", axis=1)
            elif depth_comare == "less":
                search_df["max_depth"] = search_df.apply(lambda x: x[sample_t_depths].max(), axis=1)
                search_df = search_df[search_df["max_depth"] < depth]
                search_df = search_df.drop("max_depth", axis=1)
            elif depth_comare == "lessoreq":
                search_df["max_depth"] = search_df.apply(lambda x: x[sample_t_depths].max(), axis=1)
                search_df = search_df[search_df["max_depth"] <= depth]
                search_df = search_df.drop("max_depth", axis=1)
            elif depth_comare == "equal":
                search_df["reamin"] = search_df.apply(get_all,args=(depth,sample_t_depths,), axis=1)
                search_df = search_df[search_df["reamin"] == 1]
                search_df = search_df.drop("reamin", axis=1)
            search_df = search_df.drop(sample_t_depths, axis=1)
        else:
            search_df = search_df[
                ["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt",
                 "Total depth", "QUAL", "Anno", "MUT type", "MUT info", "type", sample + "_sampledepth",
                 sample + "genotype"]]
            selected_column = sample+"_sampledepth"
            search_df[sample + "_depth"] = search_df[selected_column].apply(
                lambda x: sum([int(i) if int(i) else 0 for i in x.replace(".","0").split(",")]))
            new_selected_column = sample + "_depth"
            if depth_comare == "greater":
                search_df = search_df[search_df[new_selected_column] > depth]
            elif depth_comare == "greateroreq":
                search_df = search_df[search_df[new_selected_column] >= depth]
            elif depth_comare == "less":
                search_df = search_df[search_df[new_selected_column] < depth]
            elif depth_comare == "lessoreq":
                search_df = search_df[search_df[new_selected_column] <= depth]
            elif depth_comare == "equal":
                search_df = search_df[search_df[new_selected_column] == depth]
            search_df = search_df.drop([new_selected_column], axis=1)
        if search_df.shape[0] == 0:
            self.set_error(u"没有符合要求的snp事件")
        else:
            search_df.to_csv(os.path.join(self.output_dir, "target_detail.txt"), index=False, sep="\t")

        if self.option("sample") == "all":
            use_col = ["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt",
                      "Total depth", "QUAL", "Anno", "MUT type", "MUT info", "type"]
            for sample in samples:
                use_col.extend([sample + "_sampledepth",sample + "genotype"])
            search_df = search_df[["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt",
                      "Total depth", "QUAL", "Anno", "MUT type", "MUT info", "type", sample + "_sampledepth",
                      sample + "genotype"]]
        else:
            search_df = search_df[
                ["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt",
                 "Total depth", "QUAL", "Anno", "MUT type", "MUT info", "type", sample + "_sampledepth",
                 sample + "genotype"]]
        # search_df = snp_detail[(snp_detail["event_type"] == self.option("event_type"))
        #                              & (snp_detail["gene_id"].isin(target_genes)) & (snp_detail["sample"] ==  self.option("sample"))]
        # search_df.to_csv(os.path.join(self.output_dir,"target_deati.txt"),index=False,sep="\t")

    def get_depth_columns(self,search_df):
        all_columns = search_df.columns
        sample_depths = []
        sample_t_depths = []
        sample_names = []
        for column in all_columns:
            if column.endswith("_sampledepth"):
                sample = column.split("_sampledepth")[0]
                sample_names.append(sample)
                sample_depths.append(column)
                search_df[sample+"_depth"] = search_df[column].apply(lambda x: sum([int(i) for i in x.replace(".","0").split(",")]))
                sample_t_depths.append(sample+"_depth")
        return sample_depths,sample_names,sample_t_depths


    def set_output(self):
        # all ready write results to output
        pass




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "ASprofile" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ASprofile.asprofile",
            "instant": True,
            "options": {
                "transcripts": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/merged.gtf',
                "hdrs": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/ref.fa.hdrs',
                'sample': 'Merge'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

