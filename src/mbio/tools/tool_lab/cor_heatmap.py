# -*- coding: utf-8 -*-
# author: zzg 20201123
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.statistical.correlation import corr_heatmap
import subprocess
import os
import re


class CorHeatmapAgent(Agent):
    """
    CorHeatmapAgent:相关性热图
    version: 2
    """
    def __init__(self, parent):
        super(CorHeatmapAgent, self).__init__(parent)
        options = [
            {"name": "table1", "type": "infile", "format": "tool_lab.table"},
            {"name": "table2", "type": "infile", "format": "tool_lab.table"},
            {"name": "method", "type": "string", "default": "spearman"},
            {"name": "row_cluster", "type": "string", "default": "none"},
            {"name": "column_cluster", "type": "string", "default": "none"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("table1").is_set:
            raise OptionError('必须提供矩阵1文件')
        if self.option("method") not in ["pearson", "spearman", "kendalltau"]:
            raise OptionError('不支持该相关系数方法')
        if not self.option('table2').is_set:
            raise OptionError('必须提供矩阵2文件')
        table1 = self.option("table1").prop['path']
        table2 = self.option("table2").prop['path']
        with open(table1) as f:
            a = f.readlines()
            table1_sample = a[0]
        with open(table2) as f:
            b = f.readlines()
            table2_sample = b[0]
        table1_list = table1_sample.strip("\r\n").split("\t")[1:]
        table2_list = table2_sample.strip("\r\n").split("\t")[1:]
        if len(table1_list) != len(table2_list):
            raise OptionError('两个矩阵文件样本数量不一致')
        different_list = []
        for i in table2_list:
            if i not in table1_list:
                different_list.append(i)
        for i in table1_list:
            if i not in table2_list:
                different_list.append(i)
        if different_list:
            raise OptionError(('样本名检测不一致:{}').format(",".join(different_list)))


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 5
        self._memory = '5G'

    def end(self):
        super(CorHeatmapAgent, self).end()


class CorHeatmapTool(Tool):
    def __init__(self, config):
        super(CorHeatmapTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.perl_path = self.config.SOFTWARE_DIR + "/program/perl/perls/perl-5.24.0/bin/"
        self.r_path = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/Rscript'
        self.hcluster_script_path = self.config.SOFTWARE_DIR + "/bioinfo/statistical/scripts/"
        self.Rscript_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/"
        self.cmd_path = '{}/miniconda2/bin/python {}/statistical/pearsonsCorrelation.py'\
            .format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        # self.cmd_path=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/statistical/scripts/pearsonsCorrelation.py')
        self.table1 = self.option('table1').prop['path']
        self.table2 = self.t_table(self.option('table2').prop['path'], self.work_dir + "/new_table2.xls")
        self.name_to_name = {}
        self.env_name = {}

    def run(self):
        """
        运行
        """
        super(CorHeatmapTool, self).run()
        self.run_pearsonsCorrelation()
        self.run_heatmap()
        self.set_output()
        self.end()

    def run_pearsonsCorrelation(self):
        """
        run pearsonsCorrelation.py
        """
        global method
        if self.option("method") in ['pearson']:
            method = 'pearsonr'
        elif self.option("method") in ['spearman']:
            method = 'spearmanr'
        else:
            method = self.option("method")
        cmd = self.cmd_path
        cmd += " %s %s %s %s %s" % (self.table1, self.table2, "./correlation_"+ method + ".xls"
                                    , "./pvalue_"+ method + ".xls",method)
        self.logger.info('运行pearsonsCorrelation.py计算correlation')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('Pearsons Correlation 计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('Pearsons Correlation 计算失败')
            self.set_error('pearsonsCorrelation.py 计算失败')
        self.logger.info('运行pearsonsCorrelation.py计算correlation完成')


    def get_name(self, table):
        with open(table, "r") as f, open(self.work_dir + "/tem.collection.xls", "w") as w, \
                open("name_to_name.xls", "w") as nf, open("env_name.xls", "w") as ew:
            first_line = f.readline()
            # w.write(first_line)
            col_names = first_line.strip().split("\t")
            w.write(col_names[0])
            e = 1
            for c in col_names[1:]:
                env_new_name = "colnew"+str(e)
                w.write("\t" + env_new_name)
                self.env_name[env_new_name] = c
                e += 1
            n = 1
            w.write("\n")
            self.logger.info(self.env_name)
            ew.write(str(self.env_name))
            for line in f:
                line = line.split("\t")
                name = line[0]
                new_name = "name"+str(n)
                nf.write(new_name + "\t" + name + "\n")
                self.name_to_name[new_name] = name
                n += 1
                new_line = new_name+"\t"+"\t".join(line[1:])
                w.write(new_line)
        return n

    def run_heatmap(self):
        line_num = self.get_name(self.work_dir + "/correlation_"+ method + ".xls")
        if line_num < 2:
            # raise Exception('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行') #modified by hongdongxuan 20170406
            self.set_error('相关系数矩阵行数/物种数小于2，请尝试切换水平重新运行')
        if self.option("row_cluster") == "" and self.option("column_cluster") == "": # modified by qingchen.zhang @20200219
            pass
        else:
            corr_heatmap(self.work_dir + "/tem.collection.xls", "row_tree_tmp.tre", "column_tree_tmp.tre",
                         self.option("row_cluster"), self.option("column_cluster"))
            cmd = self.r_path + " run_corr_heatmap.r"
            try:
                subprocess.check_output(cmd, shell=True)
                self.logger.info('heatmap计算成功')
            except subprocess.CalledProcessError:
                self.logger.info('heatmap计算失败')
                self.set_error('heatmap计算失败')
            self.logger.info('生成树文件成功')

    def dashrepl(self, matchobj):
        return self.name_to_name[matchobj.groups()[0]]

    def dashrepl_env(self, matchobj):
        return self.env_name[matchobj.groups()[0]]

    def set_output(self):
        newpath = self.output_dir + "/correlation_"+ method + ".xls"
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(self.work_dir + "/correlation_"+ method + ".xls",
                self.output_dir + "/correlation_"+ method + ".xls")

        newpath2 = self.output_dir + "/pvalue_"+ method + ".xls"
        if os.path.exists(newpath2):
            os.remove(newpath2)
        os.link(self.work_dir + "/pvalue_"+ method + ".xls",
                self.output_dir + "/pvalue_"+ method + ".xls")

        row_tree_path = self.work_dir + "/column_tree_tmp.tre"
        if os.path.exists(row_tree_path):
            with open(row_tree_path, "r") as f, open(self.work_dir + "/column_tree.tre", "w") as w:
                species_tree = f.readline().strip()
                new_species_tree = re.sub(r"(name\d+)", self.dashrepl, species_tree)
                w.write(new_species_tree)
            os.link(self.work_dir + "/column_tree.tre",self.output_dir + "/column_tree.tre")
        column_tree_path = self.work_dir + "/row_tree_tmp.tre"
        if os.path.exists(column_tree_path):
            with open(column_tree_path, "r") as f, open(self.work_dir + "/row_tree.tre", "w") as w:
                env_tree = f.readline().strip()
                new_species_tree = re.sub(r"(colnew\d+)", self.dashrepl_env, env_tree)
                w.write(new_species_tree)
            os.link(self.work_dir + "/row_tree.tre", self.output_dir + "/row_tree.tre")

    def t_table(self, table_file, new_table):  # 表格转置
        """
    	转换颠倒表格内容
    	"""
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)
        return new_table