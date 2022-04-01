# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import defaultdict
import pandas as pd
import os
from bson.objectid import ObjectId

class KeggRichstatAgent(Agent):
    """
    Kegg富集分析
    version v1.0.1
    author: qiuping
    last_modify: 2016.11.23
    """
    def __init__(self, parent):
        super(KeggRichstatAgent, self).__init__(parent)
        options = [
            {"name": "geneset_kegg_enrich_info", "type": "string"},
            # 筛选过后的基因集kegg详情
            {"name": "geneset_kegg_enrich_stat_id", "type": "string"},
            # 客户输入的gene名字文件
            #{"name": "geneset_name", "type": "string"}  # 基因集名称

        ]
        self.add_option(options)
        self.step.add_steps("kegg_enrich_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kegg_enrich_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.kegg_enrich_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '4G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r"kegg_enrichment.xls$", "xls", "kegg富集分析结果"]
        ])
        super(KeggRichstatAgent, self).end()


class KeggRichstatTool(Tool):
    def __init__(self, config):
        super(KeggRichstatTool, self).__init__(config)


    def run(self):
        """
        运行
        :return:
        """
        super(KeggRichstatTool, self).run()
        list_custom = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
                       'Cellular Processes', 'Organismal Systems', 'Human Diseases',
                       'Drug Development']
        new_data = pd.read_table(self.option("geneset_kegg_enrich_info"), sep='\t', header=0)
        new_data.groupby("second_category")
        group_obj = new_data.groupby("second_category")
        groups = group_obj.groups.keys()
        genesets = ["seq_list"]
        result = defaultdict(dict)
        for each in groups:
            first = new_data.loc[new_data["second_category"] == each]['first_category']
            first = first.to_dict().values()[0]
            for geneset in genesets:
                group_detail = group_obj.get_group(each)
                genes = list()
                for g in group_detail[geneset]:
                    if not pd.isnull(g):
                        tmp = g.split('|')
                        genes += [x for x in tmp]
                        genes = list(set(genes))
                    else:
                        genes = []
                result[geneset][each] = [len(genes), first]
        a = pd.DataFrame(result)
        a.reset_index(inplace=True)
        a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
        a.to_csv(self.work_dir + "/" + "k", sep='\t', index=False)
        with open(self.work_dir + "/" + "k") as f1, open(self.work_dir + "/" + "kegg_statistic", "w") as fw:
            header = f1.readline()
            fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + "num" + "\n")
            for line in f1:
                line_split = line.strip().split("\t")
                sec = line_split[0]
                num1 = line_split[1].strip("[]").split(",")[0]
                first_cate = line_split[1].strip("[]").split(",")[1].strip().strip("'")
                fw.write(first_cate + "\t" + sec + "\t" + num1 + "\n")
        df_a = pd.read_table(self.work_dir + "/" + "kegg_statistic", header=0, sep="\t")
        appended_data_new_1 = []
        for i in list_custom:
            if i in list(df_a.first_category):
                data = df_a.loc[df_a['first_category'] == i]
                appended_data_new_1.append(data)

        appended_data_new_1 = pd.concat(appended_data_new_1)
        data_new = appended_data_new_1.to_dict('records')
        appended_data_new_1.to_csv(self.work_dir + "/" + "kegg_statistic", sep='\t', index=False)
        self.end()


