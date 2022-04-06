# !usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import numpy as np
import subprocess
from mbio.packages.tool_lab.correlation_new import correlation
import os
import math
import re


class CorrelationAgent(Agent):
    """
    CorHeatmapAgent:用于生成之间的correlation
    """
    def __init__(self, parent):
        super(CorrelationAgent, self).__init__(parent)
        options = [
            {"name": "corr_table", "type": "infile", "format": "toolapps.table"},
            {"name": "strategy", "type": "string", "default": "sample"},
            {"name": "sample", "type": "string", "default": "col"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},
            {"name": "group_method", "type": "string", "default": "mean"},  # mean, sum, median
            {"name": "method", "type": "string", "default": "pearson"},
            {"name": "distance_method", "type": "string", "default": "euclidean"},
            {"name": "cluster", "type": "string", "default": "complete"},
            {"name": "cluster_tree", "type": "string", "default": "True"},
            {"name": "log", "type": "string", "default": "True"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if self.option("strategy") not in ["group", "sample"]:
            raise OptionError("不支持该分析策略")
        if self.option("method") not in ["pearson", "spearman"]:
            raise OptionError("不支持该相关系数方法")
        if self.option("distance_method") not in ["euclidean", "maximum", "manhattan",
                                                  "canberra", "binary", "minkowski"]:
            raise OptionError("不支持该距离算法")
        if self.option("cluster") not in ["complete", "single", "average", "no"]:
            raise OptionError("不支持该聚类方式")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        """
        计算结束
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Correlation计算结果输出目录"],
            ["./correlation_matrix.xls", "xls", "Correlation矩阵"],
            ["./cluster.tre", "tre", "聚类结果"]
        ])
        super(CorrelationAgent, self).end()


class CorrelationTool(Tool):
    """
    用于相关性分析
    """
    def __init__(self, config):
        super(CorrelationTool, self).__init__(config)
        self._version = 'v2.1-20140214'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        # self.fpkm_path = self.option("fpkm").prop["path"]
        self.r_path = '/program/R-3.3.1/bin/'
        self.perl_path = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin/'
        self.Rscript_path = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/'
        self.hcluster_script_path = self.config.SOFTWARE_DIR + "/bioinfo/statistical/scripts/"
        self.cluster_method = self.option('cluster') if self.option('cluster') != 'no' else "complete"  # 如果不聚类，
        # 后面不导表
        # self.cmd_path = os.path.join(self.config.SOFTWARE_DIR, '/corr_heatmap.r')
        # self.cmd_path = '{}/miniconda2/bin/python {}/src/mbio/packages/metabolome/scripts/cluster.py'\
        #     .format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)

    def get_corr_table(self):
        """
        选择要进行分析的数据
        """
        corr_table = self.option('corr_table').path
        with open(corr_table, 'r') as r:
            for line in r:
                if line.split('\t')[0].strip().startswith('#'):
                    self.set_error("上传的表格中第一行不能以#号开始！")
                    break
        group_table = self.option('group_table').path
        if self.option("strategy") in ["sample"]:
            if self.option("sample") in ["col"]:
                with open(self.work_dir + '/new_corr_table.csv', 'w'):
                    path1 = self.work_dir + '/new_corr_table.csv'
                    if os.path.exists(path1):
                        os.remove(path1)
                    os.link(corr_table, self.work_dir + '/new_corr_table.csv')
            else:
                with open(corr_table, 'r') as f:
                    corr_list = [i.rstrip().split('\t') for i in f.readlines()]
                    new_corr_list = map(lambda *a: '\t'.join(a) + '\n', *corr_list)
                    with open(self.work_dir + '/new_corr_table1.csv', 'w') as fw1:
                        fw1.writelines(new_corr_list)
                    with open(self.work_dir + '/new_corr_table.csv', 'w'):
                        path2 = self.work_dir + '/new_corr_table.csv'
                        if os.path.exists(path2):
                            os.remove(path2)
                        os.link(self.work_dir + '/new_corr_table1.csv', self.work_dir + '/new_corr_table.csv')
        else:
            if self.option("sample") in ["col"]:
                with open(corr_table, 'r') as f:
                    corr_list2 = [i.rstrip().split('\t') for i in f.readlines()]
                    new_corr_list2 = map(lambda *a: '\t'.join(a) + '\n', *corr_list2)
                    with open(self.work_dir + "/new_corr_table1.csv", 'w') as fw1:
                        fw1.writelines(new_corr_list2)
            else:
                with open(self.work_dir + '/new_corr_table1.csv', 'w'):
                    path3 = self.work_dir + '/new_corr_table1.csv'
                    if os.path.exists(path3):
                        os.remove(path3)
                    os.link(corr_table, self.work_dir + '/new_corr_table1.csv')
            table = self.work_dir + '/new_corr_table1.csv'
            df = pd.read_table(table)
            n1 = 0
            list0 = []
            g_list2 = []
            with open(table, 'r')as fw1:
                for i in fw1.readlines():
                    corr_list3 = i.rstrip().split('\t')
                    if n1 == 0:
                        n1 += 1
                    else:
                        with open(group_table, 'r') as g:
                            for j in g.readlines():
                                group_list = j.strip().split('\t')
                                if corr_list3[0] == group_list[0]:
                                    list0.append(group_list[1])
                                    n1 += 1
            df.insert(1, 'group', list0)
            if self.option("group_method") in ["sum"]:
                g_df = df.groupby('group').sum()
            elif self.option("group_method") in ["mean"]:
                g_df = df.groupby('group').mean()
            else:
                g_df = df.groupby('group').median()
            g_df.to_csv(self.work_dir + '/new_corr_table2.csv', sep='\t', header=1, index=0)
            new_g_df = g_df.T
            g_list2.append("name")
            g_list1 = list(new_g_df.keys())
            g_list2.extend(g_list1)
            g_list = [i.replace("\r\n", "") for i in g_list2]
            g_list = [i.strip() for i in g_list]  # 过滤掉样本后面有\n
            # new_g_df.to_csv(self.work_dir + '/new_corr_table2.csv', sep='\t', header=0)
            df0 = pd.read_csv(self.work_dir + '/new_corr_table2.csv', sep='\t', header=None, index_col=False)
            df0.insert(0, "name", g_list)
            df0_t = df0.T
            df0_t.to_csv(self.work_dir + '/new_corr_table.csv', sep='\t', header=0, index=0)

    def run(self):
        """
        运行
        """
        super(CorrelationTool, self).run()
        self.get_corr_table() # 对输入数据进行处理，得到new_corr_table.csv
        self.run_log() # 取对数，得到new_table.csv
        self.run_correlation() # 相关性计算并聚类，得到correlation.xls， corr_col.tre, corr_row.tre
        # self.run_cluster() # 聚类，得到cluster.xls
        self.set_output()
        # self.set_table()
        self.end()

    def run_log(self):
        """
        取对数计算
        """
        new_log_table = self.work_dir + '/new_corr_table.csv'
        if self.option("log") in ["True"]:
            with open(new_log_table, 'r') as f, open(self.work_dir + '/new_table.csv', 'w'):
                n1 = 0
                df = pd.DataFrame()
                for i in f.readlines():
                    n1 += 1
                    new_list = []
                    if n1 == 1:
                        list = i.rstrip().split('\t')
                        df.insert((n1-1), "name", list)
                    else:
                        list = i.split('\t')
                        num = 0
                        for n in range(len(list)):
                            num += 1
                            if num == 1:
                                new_list.append(list[0])
                            else:
                                b = float(list[num-1])
                                try:
                                    a = math.log(b, 10)
                                except:
                                    self.set_error("{}, 数值不能为0或者负数".format(b))
                                else:
                                    new_list.append(a)
                        df.insert((n1-1), (n1-1), new_list)
                    new_df = df.T
                    new_df.to_csv(self.work_dir + '/new_table.csv', sep='\t', header=0, index=0)
        else:
            with open(self.work_dir + '/new_table.csv', 'w'):
                path7 = self.work_dir + '/new_table.csv'
                if os.path.exists(path7):
                    os.remove(path7)
                os.link(new_log_table, self.work_dir + '/new_table.csv')

    def run_correlation(self):
        """
        相关性计算
        """
        self.logger.info(self.cluster_method)
        new_table = self.work_dir + '/new_table.csv'
        self.logger.info(self.work_dir + 'correlation_matrix.xls')
        correlation(new_table, self.work_dir + '/correlation_matrix.xls', self.work_dir + '/pvalue_matrix.xls',
                    self.work_dir + '/tvalue_matrix.xls', self.work_dir + '/correlation_heatmap.pdf',
                    self.work_dir + '/corr_col.tre', self.work_dir + '/corr_row.tre', self.option('method'),
                    self.option('distance_method'), self.option('distance_method'), self.cluster_method)
        cmd = self.r_path + "Rscript run_correlation_new.r"
        self.logger.info("开始运行correlation检验")
        command = self.add_command("correlation_new", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("correlation运行完成")
        else:
            self.set_error("correlation运行出错！")

    def set_output(self):
        """
        设置输出文件路径
        """
        corr_path = self.output_dir + "/correlation_matrix.xls"
        if os.path.exists(corr_path):
            os.remove(corr_path)
        os.link(self.work_dir + "/correlation_matrix.xls", self.output_dir + "/correlation_matrix.xls")
        pvalue_path = self.output_dir + "/pvalue_matrix.xls"
        if os.path.exists(pvalue_path):
            os.remove(pvalue_path)
        os.link(self.work_dir + "/pvalue_matrix.xls", self.output_dir + "/pvalue_matrix.xls")
        if self.option("cluster_tree") in ["True"]:
            if self.cluster_method != "no":
                cluster_tree_row_path = self.output_dir + "/corr_row.tre"
                cluster_tree_col_path = self.output_dir + "/corr_col.tre"
                if os.path.exists(cluster_tree_row_path):
                    os.remove(cluster_tree_row_path)
                if os.path.exists(cluster_tree_col_path):
                    os.remove(cluster_tree_col_path)
                os.link(self.work_dir + "/corr_col.tre", self.output_dir + "/corr_col.tre")
                os.link(self.work_dir + "/corr_row.tre", self.output_dir + "/corr_row.tre")
        else:
            pass

    def set_table(self):
        """
        保存结果到Mongo库
        """
        api_correlation = self.api.api("tool_lab.correlation_api")
        api_correlation.insert_detail_table(self.option("main_id"), self.output_dir + "/correlation_matrix.xls",
                                            self.output_dir + "/cluster.tre")
