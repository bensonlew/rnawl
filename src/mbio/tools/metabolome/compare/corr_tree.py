# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
import numpy as np
from biocluster.core.exceptions import OptionError
import os, re, glob
import subprocess
from mbio.packages.metabolome.common import Relation
import pandas as pd


class CorrTreeAgent(Agent):
    """
    计算相关系数和及相关性聚类树的工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11
    """

    def __init__(self, parent):
        super(CorrTreeAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'infile', 'format': 'metabolome.express'},  # 表达矩阵文件
            {'name': 'sct', 'type': 'string', 'default': ''},  # 聚类算法，hierarchy、kmeans
            {'name': 'scd', 'type': 'string', 'default': ''},  # 距离计算方式
            {'name': 'corr_method', 'type': 'string', 'default': 'pearson'},  # 相关性计算方法,"spearman", "pearson", "kendall"
            {'name': 'n_cluster', 'type': 'int', 'default': 0},  # 聚类数目，kmeans时使用
            {'name': 'scm', 'type': 'string', 'default': ''},  # 聚类方式, hierarchy时使用，"complete","average","single",""
            {'name': 'metab_trans', 'type': 'infile', "format": "sequence.profile_table"},
            # 结果文件转化id使用，第一列为old_name,第二列为最后结果要使用的id, 可选参数
            {'name': 'file_tran', 'type': 'bool', 'default': False},
            {'name': 'transform', 'type':'string', 'default': 'none'} #默认值none，绝对不能改动。
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp").is_set:
            raise OptionError("请传入丰度矩阵！", code="34700701")
        if self.option("sct") == "kmeans":
            if not self.option("n_cluster") > 0:
                raise OptionError("使用kmeans时必须聚类数目必须大于1！", code="34700702")
        if self.option("sct") == "hierarchy":
            if not self.option("scm") or not self.option("scd"):
                raise OptionError("使用hierarchy时必须传入层次聚类方式和距离计算方式！", code="34700703")
        elif self.option("sct") == "kmeans":
            if not self.option("scd") or not self.option("n_cluster"):
                raise OptionError("使用kmeans时必须传入子聚类数目和距离计算方式！", code="34700704")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        super(CorrTreeAgent, self).end()


class CorrTreeTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CorrTreeTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.script = self.config.PACKAGE_DIR + "/metabolome/scripts/corr_cluster.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def scale_data(self, ori_table, outDir, method='UV'):   #add method by zouguanqiing 20190605
        table = pd.read_table(ori_table, sep='\t', index_col=0)
        if method=='Par':
            scaled_data = table.apply(lambda x: (x - np.mean(x)) / np.sqrt(np.std(x, ddof=1)), axis=1)
        elif method == 'Ctr':
            scaled_data = table.apply(lambda x: (x-np.mean(x)), axis=1)
        elif method == 'UV':
            scaled_data = table.apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)
        else:
            return ori_table

        exp_profile = os.path.join(outDir, "scale_data.xls")
        scaled_data = scaled_data.fillna(0)
        scaled_data.to_csv(exp_profile, index=True, header=True, sep="\t")
        return exp_profile

    def correlation(self,table):
        #table = self.option("exp").prop["path"]
        cmd = '{} {} -exp {} -out {} --corr --ngc '.format(self.python_path, self.script, table, self.work_dir)
        cmd += " -corr_method " + self.option("corr_method")
        '''
        if self.option("row_col") == "c":
            cmd += " --ngc "
        '''
        if self.option("file_tran"):
            cmd += " --T "
        if self.option("sct"):
            cmd += " -sct " + self.option("sct")
        else:
            cmd += " --nsc "
        if self.option("sct") == "kmeans":
            cmd += " -n_cluster " + str(self.option("n_cluster"))
            cmd += " -scd " + self.option("scd")
        if self.option("sct") == "hierarchy":
            cmd += " -scm " + self.option("scm")
            cmd += " -scd " + self.option("scd")
        command = self.add_command("correlation", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("correlation运行完成")
        else:
            self.set_error("correlation运行出错!", code="34700701")

    def set_output(self):
        self.logger.info("set output")
        if self.option("metab_trans").is_set:
            self.logger.info("进行id转化...")
            table_path = self.option("metab_trans").prop["path"]
            self.metab_trans = Relation()
            map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
            if self.option("sct") == "hierarchy":
                old_tree = os.path.join(self.work_dir, "corr.cluster_tree.txt")
                if os.path.exists(old_tree):
                    self.trans_tree_name(id_name_dict, "corr.cluster_tree.txt", "corr.cluster_tree.xls")
                    self.link_file("corr.cluster_tree.txt", "metab_id.corr_tree.xls")
            elif self.option("sct") == "kmeans":
                subclusters = glob.glob(self.work_dir + '/*subcluster*')
                for each in subclusters:
                    if os.path.exists(each):
                        each = each.split("/")[-1]
                        self.trans_data(each, each, map_table, id_name_dict)
                self.metab_trans.trans_col_data(id_name_dict, self.work_dir + "/corr.kmeans_cluster.txt",
                                            self.output_dir + "/corr.kmeans_cluster.xls")
            self.trans_data("corr.xls", "corr.xls", map_table, id_name_dict)
            self.trans_data("pvalue.xls", "pvalue.xls", map_table, id_name_dict)
            self.metab_trans.add_oid_fun(self.output_dir+'/corr.xls',table_path,oid='o_id',link_k='metab_id',oldfile_link_id=0)   #20190617
            self.metab_trans.add_oid_fun(self.output_dir+'/pvalue.xls',table_path,oid='o_id',link_k='metab_id',oldfile_link_id=0) #20190617
        else:
            if self.option("sct"):
                old_tree = os.path.join(self.work_dir, "corr.cluster_tree.txt")
                if os.path.exists(old_tree):
                    self.link_file("corr.cluster_tree.txt", "corr.cluster_tree.xls")
            self.link_file("corr.xls", "corr.xls")
            self.link_file("pvalue.xls", "pvalue.xls")

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)

    def trans_data(self, oldfile, newfile, map_table, id_name_dict):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        self.metab_trans.get_trans_file(map_table, id_name_dict, oldfile, newfile, mycol=True)
        target_data = pd.read_table(newfile,sep='\t',header=0)
        map_data = pd.read_table(self.option('metab_trans').prop['path'],sep='\t',header=0)
        if 'ID' in map_data.columns:
            new_map = map_data[['ID','Metabolite']]
            new_data = pd.merge(new_map,target_data,how="right", on="Metabolite")
            new_data.to_csv(newfile,sep='\t',index=False, quoting=3)


    def trans_tree_name(self, namedict, treefile, newfile):
        """
        用正则的方式替换复原树文件中的名称
        """
        treefile = os.path.join(self.work_dir, treefile)
        newfile = os.path.join(self.output_dir, newfile)
        if not isinstance(namedict, dict):
            self.set_error('复原树的枝名称需要旧名称和当前名称的字典', code="34700702")
        #namedict = {item[1]: item[0] for item in namedict.iteritems()}
        try:
            new_names = []
            with open(treefile, 'rb') as f, open(newfile, 'wb') as w:
                tree = f.readline().strip()
                for item in namedict.iteritems():
                    tree = re.sub(item[0] + ':', item[1] + ':', tree)
                w.write(tree + "\n")
                name = f.readline().strip().split(";")
                for each in name:
                    new_name = namedict[each]
                    new_names.append(new_name)
                w.write(";".join(new_names) + "\n")
        except IOError, e:
            self.set_error('聚类树文件无法找到或者无法打开：%s' , variables=(e), code="34700703")

    def run(self):
        """
        运行
        """
        super(CorrTreeTool, self).run()
        self.exp_table = self.scale_data(self.option('exp').prop['path'],self.work_dir,method=self.option('transform'))
        self.correlation(self.exp_table)
        # self.correlation()
        # self.plot_hcluster()
        self.set_output()
        self.end()
