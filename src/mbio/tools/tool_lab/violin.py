# -*- coding: utf-8 -*-
# __author__ = 'wuqin'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import math
import pandas as pd


class ViolinAgent(Agent):
    def __init__(self, parent):
        super(ViolinAgent, self).__init__(parent)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "grouptable", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "log", "type": "string", "default": "none"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        if self.option('tooltable').is_set == 'False':
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('log') not in ['log2', 'log10', 'none']:
            raise OptionError('数据处理方式不存在', code="32702909")
        if self.option('tooltable').prop['sample_num'] < 2:
            raise OptionError('列数少于2，不可进行分析', code="32702904")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(ViolinAgent, self).end()


class ViolinTool(Tool):
    def __init__(self, config):
        super(ViolinTool, self).__init__(config)
        # self.table = self.get_table()

    def run(self):
        """
        运行
        """
        super(ViolinTool, self).run()
        self.get_table()
        self.run_violin()
        self.linkfile()
        # # self.set_db()
        self.end()

    def get_table(self):
        """
        取对数计算，并检查分组信息是否正确
        """
        new_log_table = self.option('tooltable').path
        if self.option("log") in ["log2", "log10"]:
            with open(new_log_table, 'r') as f, open(self.work_dir + '/new_table.txt', 'w'):
                n1 = 0
                df = pd.DataFrame()
                for i in f.readlines():
                    n1 += 1
                    new_list = []
                    if n1 == 1:
                        l_list = i.rstrip().split('\t')
                        df.insert((n1 - 1), "name", l_list)
                    else:
                        l_list = i.split('\t')
                        num = 0
                        for n in range(len(l_list)):
                            num += 1
                            if num == 1:
                                new_list.append(l_list[0])
                            else:
                                b = float(l_list[num - 1])
                                try:
                                    if self.option('log') == 'log2':
                                        a = math.log(b, 2)
                                    else:
                                        a = math.log(b, 10)
                                except:
                                    self.set_error("{}, 数值不能为0或者负数".format(b))
                                else:
                                    new_list.append(a)
                        df.insert((n1 - 1), (n1 - 1), new_list)
                    new_df = df.T
                    new_df.to_csv(self.work_dir + '/new_table.txt', sep='\t', header=0, index=0)
        else:
            with open(self.work_dir + '/new_table.txt', 'w'):
                path0 = self.work_dir + '/new_table.txt'
                if os.path.exists(path0):
                    os.remove(path0)
                os.link(new_log_table, self.work_dir + '/new_table.txt')
        result1 = pd.read_table(self.work_dir + '/new_table.txt')
        if self.option('grouptable').is_set:
            group0 = pd.read_table(self.option('grouptable').path)
            g_sample_list = list(group0[group0.columns[0]])
            r_sample_list1 = list(result1[result1.columns[0]])
            same_list1 = [x for x in g_sample_list if x in r_sample_list1]
            if len(same_list1) == 0:
                with open(self.work_dir + '/new_table.txt') as f, open('new_test.txt', 'w') as w:
                    table_list = [i.rstrip().split('\t') for i in f.readlines()]
                    table_list = map(lambda * t: '\t'.join(t) + '\n', *table_list)
                    w.writelines(table_list)
                result2 = pd.read_table('new_test.txt')
                r_sample_list2 = list(result2[result2.columns[0]])
                same_list2 = [y for y in g_sample_list if y in r_sample_list2]
                if len(same_list2) != 0:
                    result2.to_csv(self.work_dir + '/final_table.txt', sep='\t', index=0)
                else:
                    self.set_error("{}, 样本信息不在分组中")
            else:
                with open(self.work_dir + '/final_table.txt', 'w'):
                    path1 = self.work_dir + '/final_table.txt'
                    if os.path.exists(path1):
                        os.remove(path1)
                    os.link(self.work_dir + '/new_table.txt', self.work_dir + '/final_table.txt')
        else:
            with open(self.work_dir + '/final_table.txt', 'w'):
                path2 = self.work_dir + '/final_table.txt'
                if os.path.exists(path2):
                    os.remove(path2)
                os.link(self.work_dir + '/new_table.txt', self.work_dir + '/final_table.txt')

    def run_violin(self):
        result = pd.read_table(self.work_dir + '/final_table.txt')
        if self.option('grouptable').is_set:
            group = pd.read_table(self.option('grouptable').path)
            col_list = list(result)
            g_col_list = list(group)
            col_list[0] = g_col_list[0]
            result.columns = col_list
            new_table = pd.merge(result, group, on=g_col_list[0])
            new_table.rename(columns={g_col_list[1]: 'group'}, inplace=True)
            result_table = new_table.groupby('group').sum()
            result_table1 = new_table.groupby('group').mean()
            result_table2 = new_table.groupby('group').median()
            result_table.to_csv(self.work_dir + '/violin_sum.xls', sep='\t')
            result_table1.to_csv(self.work_dir + '/violin_mean.xls', sep='\t')
            result_table2.to_csv(self.work_dir + '/violin_median.xls', sep='\t')
            new_table.to_csv(self.work_dir + '/violin_none.xls', sep='\t', index=0)
        else:
            result_table = result.T
            result_table.to_csv(self.work_dir + '/violin.xls', sep='\t', header=0)

    def linkfile(self):
        if self.option('grouptable').is_set:
            s_file = self.work_dir + '/violin_sum.xls'
            m_file = self.work_dir + '/violin_mean.xls'
            z_file = self.work_dir + '/violin_median.xls'
            n_file = self.work_dir + '/violin_none.xls'
            s_link = self.output_dir + '/violin_sum.xls'
            m_link = self.output_dir + '/violin_mean.xls'
            z_link = self.output_dir + '/violin_median.xls'
            n_link = self.output_dir + '/violin.xls'
            if os.path.exists(s_link):
                os.remove(s_link)
            os.link(s_file, s_link)
            if os.path.exists(m_link):
                os.remove(m_link)
            os.link(m_file, m_link)
            if os.path.exists(z_link):
                os.remove(z_link)
            os.link(z_file, z_link)
            if os.path.exists(n_link):
                os.remove(n_link)
            os.link(n_file, n_link)
        else:
            t_file = self.work_dir + '/violin.xls'
            link = self.output_dir + '/violin.xls'
            if os.path.exists(link):
                os.remove(link)
            os.link(t_file, link)

    def set_db(self):
        self.logger.info("开始导表")
        api_violin = self.api.api("tool_lab.violin")
        api_violin.add_violin_detail(self.option('main_id'), self.output_dir)
        self.logger.info("导表结束")
