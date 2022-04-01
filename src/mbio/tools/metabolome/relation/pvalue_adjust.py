# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last modify date: 20220125
# last modified: zhaoyuzhuo

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd
from statsmodels.stats.multitest import multipletests


class PvalueAdjustAgent(Agent):
    """
    代谢转录关联分析--筛选转录丰度表用于后续相关性分析
    """
    def __init__(self, parent):
        super(PvalueAdjustAgent, self).__init__(parent)
        options = [
            {"name": "pvalue_table", "type": "infile", "format": "sequence.profile_table"},  # 测试用的输入转录丰度文件参数
            {"name": "corr_table", "type": "infile", "format": "sequence.profile_table"},  # 所选择的转录基因集的gene_id文件
            {"name": "padjust_method", "type": "string", "default": "fdr_bh"},  # 多重校验方法,
            {"name": "sort", "type": "string", "default": "pvalue"},
            {"name": "top", "type": "int", "default": 200}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(PvalueAdjustAgent, self).end()


class PvalueAdjustTool(Tool):
    def __init__(self, config):
        super(PvalueAdjustTool, self).__init__(config)
        self.output_file = self.output_dir + "/qvalue.xls"

    def run(self):
        """
        运行
        :return:
        """
        super(PvalueAdjustTool, self).run()
        self.p_adjust()
        self.sort_result()
        self.set_output()
        self.end()

    def p_adjust(self):
        df_corr_table = pd.read_table(self.option("corr_table").path, '\t')
        df_pvalue_table = pd.read_table(self.option("pvalue_table").path, '\t')
        df_metab = df_pvalue_table["Metabolite"]
        dfpvals = df_pvalue_table.ix[:,df_pvalue_table.columns[1:].tolist()]
        flat_pvalue = [item for sublist in dfpvals.values.tolist() for item in sublist]
        flat_adjust_pvalue = multipletests(flat_pvalue, alpha=0.05, method=self.option("padjust_method"), is_sorted=False, returnsorted=False)
        adjust_p_list = flat_adjust_pvalue[1].tolist()
        size = dfpvals.shape[1]
        nested_adjust_pvalue = [adjust_p_list[i:i + size] for i in range(0, len(adjust_p_list), size)]
        df_qvalue = pd.DataFrame(nested_adjust_pvalue, columns=dfpvals.columns)
        df_qvalue_output = pd.concat([df_metab, df_qvalue], axis=1)
        df_qvalue_output.to_csv(self.output_file, '\t', index=False)
        self.out_dir_file = os.path.join(self.output_dir, "express_correlation_info.txt")
        with open(self.out_dir_file, "w") as info_file:
            seq = ("metab_name", "gene_id", "corr", "p_value", "q_value")
            info_file.write("\t".join(seq) + "\n")
            for i in df_corr_table.index:
                for j in df_corr_table.columns[1:].tolist():
                    corr = [df_corr_table.ix[i,"Metabolite"], j, str(df_corr_table.ix[i, j])]
                    pvalue = str(df_pvalue_table.ix[i,j])
                    qvalue = str(df_qvalue_output.ix[i,j])
                    info_file.write("\t".join(corr) + "\t" + pvalue + "\t" + qvalue + '\n')

    def sort_result(self):
        df_info = pd.read_table(self.out_dir_file, "\t")
        if self.option("sort") == "pvalue":  # 以pvalue排序
            pvalue_list = df_info["p_value"].tolist()
            self.logger.info('pvalue_list为{}'.format(pvalue_list))
            pvalue_list.sort(reverse=False)  # pvalue倒序排列
            self.logger.info('pvalue_list排序后为{}'.format(pvalue_list))
            threshold = pvalue_list[self.option("top")]
            df_info = df_info[(df_info["p_value"] <= threshold)]
            self.logger.info('df_info为{}'.format(df_info))
        elif self.option("sort") == "corr":  # 以相关性系数值排序
            corr_list = df_info["corr"].tolist()
            abs_corr_list = []
            for i in corr_list:
                abs_corr_list.append(abs(i))
            abs_corr_list.sort(reverse=True)
            threshold = abs_corr_list[self.option("top")]  # 根据设置参数取第{}值作为筛选条件
            df_info = df_info[((df_info["corr"] >= threshold) | (df_info["corr"] <= -threshold))]
            self.logger.info('df_info为{}'.format(df_info))
        sort_out_file = self.output_dir + "/sorted_express_correlation_info.txt"
        df_info.to_csv(sort_out_file, "\t", index=False)

    def set_output(self):
        self.logger.info('开始设置输出结果文件')