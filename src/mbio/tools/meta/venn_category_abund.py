# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modifiy:2017.11.01

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
import math


class VennCategoryAbundAgent(Agent):
    """
    生成venn每个Category对应的百分比文件，用于venn分析中的物种/功能/基因分布饼图
    """
    def __init__(self, parent):
        super(VennCategoryAbundAgent, self).__init__(parent)
        options = [
            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},  # 用于计算venn的丰度表格
            {"name": "venn_table", "type": "infile", "format": "graph.venn_table"},  # venn分析得到的“venn_table.xls”
            {"name": "other", "type": "float", "default": "0.0001,0.0005,0.001,0.01,0.02,0.05"},
            {"name": "out_table", "type": "outfile", "format": "graph.venn_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("abund_file").is_set:
            raise OptionError("请传入venn分析的丰度表格文件！", code="32706701")
        if not self.option("venn_table").is_set:
            raise OptionError("请传入venn分析结果的venn_table.xls文件！", code="32706702")

    def set_resource(self):
        self._cpu = 1
        table_path = self.option("abund_file").prop['path']
        table_size = os.path.getsize(table_path) / float(1024*1024)
        memory_number = int(math.ceil(table_size) + 3)
        self._memory = '{}G'.format(str(memory_number))

    def end(self):
        super(VennCategoryAbundAgent, self).end()


class VennCategoryAbundTool(Tool):
    def __init__(self, config):
        super(VennCategoryAbundTool, self).__init__(config)

    def get_venn_pie_abund(self):
        others_set = self.option('other').split(',')
        abund_file = pd.read_table(self.option("abund_file").prop["path"], sep='\t', header=0, index_col=0)
        venn_table = pd.read_table(self.option("venn_table").prop["path"], sep='\t', header=None)
        venn_pie_data = os.path.join(self.output_dir, "venn_pie_abund.xls")
        with open(venn_pie_data, 'wb') as w:
            for n in range(len(venn_table)):
                category_name = venn_table.iat[n, 0]
                #abund = abund_file.ix[venn_table.iat[n, 2].split(',')]
                if venn_table.iat[n, 1]:
                    abund = abund_file.ix[venn_table.iat[n, 2].split(';')]  # 与venn中物种/功能间的连接符相关
                    abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
                    abund_tmp = abund[["Col_sum"]]
                    col_sum = sum(abund['Col_sum'])
                    abund_tmp = abund_tmp.apply(lambda x: x/col_sum, axis=1)
                    for i in others_set:
                        o = float(i)
                        df = abund_tmp.ix[list((abund_tmp > o).any(axis=1))]
                        # new_df2 = df.copy
                        others = abund_tmp.ix[list((abund_tmp < o).all(axis=1))]
                        if len(others) > 0:
                            df.loc["others"] = others.apply(lambda x: x.sum(), axis=0)
                        data_list = []
                        for k in range(len(df)):
                            data = str(df.index[k]) + ";" + str(df.iat[k, 0])
                            data_list.append(data)
                        data = str(data_list)
                        w.write(category_name + "\t" + i + "\t" + data + "\n")

    def run(self):
        super(VennCategoryAbundTool, self).run()
        self.get_venn_pie_abund()
        self.option("out_table").set_path(os.path.join(self.output_dir, "venn_pie_abund.xls"))
        self.end()
