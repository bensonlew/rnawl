# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.06.05

import os, re, subprocess, glob
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.scripts.mul_diff_stat import MulDiffStat
import traceback
from mbio.packages.metabolome.common import Relation


class MetabClusterAgent(Agent):
    """
    代谢集代谢物聚类分析
    """

    def __init__(self, parent):
        super(MetabClusterAgent, self).__init__(parent)
        options = [
            {'name': 'exp', 'type': 'infile', 'format': 'metabolome.express,sequence.profile_table'},  # 表达矩阵文件
            {'name': 'sct', 'type': 'string', 'default': ''},  # 样本聚类算法，hierarchy、无
            {'name': 'scd', 'type': 'string', 'default': ''},  # 样本距离计算方式
            {'name': 'scm', 'type': 'string', 'default': ''},  # 样本聚类方式,"complete","average","single"
            {'name': 'mct', 'type': 'string', 'default': ''},  # 代谢物聚类算法，hierarchy、kmeans、无
            {'name': 'mcd', 'type': 'string', 'default': ''},  # 代谢物距离计算方式
            {'name': 'n_cluster', 'type': 'int', 'default': 0},  # 代谢物聚类数目，kmeans时使用
            {'name': 'mcm', 'type': 'string', 'default': ''},  # 代谢物聚类方式, hierarchy时使用，"complete","average","single"
            {'name': 'metab_trans', 'type': 'infile', "format": "sequence.profile_table"},  # 转化metab_id
            {'name': 'before_scale', 'type': 'infile', "format": "sequence.profile_table"},  # 使用scale时，生成的原始丰度文件
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp").is_set:
            raise OptionError("请传入丰度矩阵！", code="34701301")
        if self.option("mct") == "kmeans":
            if not self.option("n_cluster") > 0:
                raise OptionError("使用kmeans时必须聚类数目必须大于1！", code="34701302")
        if self.option("sct") == "hierarchy":
            if not self.option("scm"):
                raise OptionError("样本使用hierarchy时必须传入层次聚类方式！", code="34701303")
        if self.option("mct") == "hierarchy":
            if not self.option("mcm"):
                raise OptionError("代谢物使用hierarchy时必须传入层次聚类方式！", code="34701304")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '8G'

    def end(self):
        super(MetabClusterAgent, self).end()


class MetabClusterTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MetabClusterTool, self).__init__(config)
        env_path = ':'.join([
            os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'),
            os.path.join(self.config.SOFTWARE_DIR, 'gcc/5.1.0/bin'),
            os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin')
        ])
        ld_lib_path = ':'.join([
            os.path.join(self.config.SOFTWARE_DIR, 'gcc/5.1.0/lib64'),
            os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/lib64/R/lib')
        ])
        self._r_home = os.path.join(self.config.SOFTWARE_DIR, "program/R-3.3.1/lib64/R")
        self.set_environ(PATH=env_path, R_HOME=self._r_home, LD_LIBRARY_PATH=ld_lib_path)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'cluster': os.path.join(self.config.PACKAGE_DIR, 'metabolome/scripts/corr_cluster.py')
        }

    def run(self):
        """
        运行
        """
        super(MetabClusterTool, self).run()
        self.logger.info("开始运行命令！")
        self.run_cluster()
        self.set_output()

    def run_cluster(self):
        """
        """
        table = self.option("exp").prop["path"]
        cmd = '{} {} -exp {} -out {}'.format(self.program['python'], self.script['cluster'], table, self.work_dir)
        if self.option("sct"):
            cmd += " -sct " + self.option("sct")
        else:
            cmd += " --nsc"
        if self.option("mct"):
            cmd += " -gct " + self.option("mct")
        else:
            cmd += " --ngc"
        if self.option("scd"):
            cmd += " -scd " + self.option("scd")
        if self.option("mcd"):
            cmd += " -gcd " + self.option("mcd")
        #if self.option("mct") == "kmeans":
        cmd += " -n_clusters " + str(self.option("n_cluster"))
        if self.option("sct") == "hierarchy":
            cmd += " -scm " + self.option("scm")
        if self.option("mct") == "hierarchy":
            cmd += " -gcm " + self.option("mcm")
        command = self.add_command("cluster", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("metab cluster succeed")
        else:
            self.set_error("metab cluster failed!", code="34701301")
            raise Exception("metab cluster failed!")

    def set_output(self):
        self.logger.info("set output")
        table_path = self.option("metab_trans").prop["path"]
        self.metab_trans = Relation()
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
        self.logger.info(map_table.head())
        self.logger.info(table_path)
        if self.option("sct"):
            self.link_file("col.cluster_tree.txt", "col.cluster_tree.xls")
        self.logger.info("进行id转化...")
        if self.option("mct") == "hierarchy":
            self.trans_tree_name(id_name_dict, "row.cluster_tree.txt", "row.cluster_tree.xls")
            # self.link_file("row.cluster_tree.txt", "row_id.cluster_tree.xls")
            #20190618 zouguanqing
            #row.subcluster_*.xls 更名，移位， 转id
            # subclusters = glob.glob(self.work_dir + '/*subcluster*')
            # for each in subclusters:
            #     each = each.split("/")[-1]
            #     newname = each.replace("row", "metab")
            #     self.trans_data(each, newname, map_table, id_name_dict)
            self.metab_trans.trans_col_data(id_name_dict, self.work_dir + "/row.cluster.txt",     #zouguanqing 20190618
                                            self.output_dir + "/row.hierarchy_cluster.xls")

        elif self.option("mct") == "kmeans":
            # subclusters = glob.glob(self.work_dir + '/*subcluster*')
            # for each in subclusters:
            #     each = each.split("/")[-1]
            #     self.trans_data(each, each, map_table, id_name_dict)
            self.metab_trans.trans_col_data(id_name_dict, self.work_dir + "/row.kmeans_cluster.txt",
                                            self.output_dir + "/row.kmeans_cluster.xls")
        if "before_scale" in self.get_option_object().keys() and self.option("before_scale").is_set:
            origin_file = self.option("before_scale").prop["path"]
            self.order_origin(origin_file, "cluster_exp.xls", "cluster_exp_ori.xls")
            self.trans_data("cluster_exp_ori.xls", "cluster_exp.xls", map_table, id_name_dict)
            self.trans_data("cluster_exp.xls", "cluster_scale_exp.xls", map_table, id_name_dict,drop_metab_id=False)
            self.metab_trans.add_oid_fun(self.output_dir+ '/cluster_scale_exp.xls', table_path, link_k='Metabolite',oldfile_link_id=0)
        else:
            self.trans_data("cluster_exp.xls", "cluster_exp.xls", map_table, id_name_dict,drop_metab_id=False)
            self.metab_trans.add_oid_fun(self.output_dir+ '/cluster_exp.xls', table_path, link_k='Metabolite',oldfile_link_id=0)
        self.wait()
        self.end()

    def link_file(self, oldfile, newfile):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        if os.path.exists(oldfile):
            os.link(oldfile, newfile)

    def trans_data(self, oldfile, newfile, map_table, id_name_dict,drop_metab_id=True):
        oldfile = os.path.join(self.work_dir, oldfile)
        newfile = os.path.join(self.output_dir, newfile)
        self.metab_trans.get_trans_file(map_table, id_name_dict, oldfile, newfile,drop_metab_id=drop_metab_id)

    def order_origin(self, origin_file, order_file, outfile):
        order_file = os.path.join(self.work_dir, order_file)
        outfile = os.path.join(self.work_dir, outfile)
        table = pd.read_table(order_file, sep="\t", index_col=0)
        select_names = table.index
        origin_table = pd.read_table(origin_file, sep="\t", index_col=0)
        selecl_origin = origin_table.loc[select_names,]
        selecl_origin.to_csv(outfile, index=True, header=True,sep="\t")

    def trans_tree_name(self, namedict, treefile, newfile):
        """
        用正则的方式替换复原树文件中的名称
        """
        treefile = os.path.join(self.work_dir, treefile)
        newfile = os.path.join(self.output_dir, newfile)
        if not isinstance(namedict, dict):
            self.set_error('复原树的枝名称需要旧名称和当前名称的字典', code="34701303")
        if os.path.exists(treefile):
            try:
                new_names = []
                with open(treefile, 'rb') as f, open(newfile, 'wb') as w:
                    tree = f.readline().strip()
                    for item in namedict.iteritems():
                        tree = re.sub(str(item[0]) + ':', str(item[1]) + ':', tree)
                    w.write(tree + "\n")
                    name = f.readline().strip().split(";")
                    for each in name:
                        new_name = namedict[each]
                        new_names.append(new_name)
                    w.write(";".join(new_names) + "\n")
            except IOError, e:
                self.set_error('聚类树文件无法找到或者无法打开：%s', variables=(e), code="34701304")