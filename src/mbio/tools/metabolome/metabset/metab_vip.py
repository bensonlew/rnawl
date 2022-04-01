# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.06.05

import os, re, subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.scripts.mul_diff_stat import MulDiffStat
import traceback
import unittest
from mbio.packages.metabolome.common import Relation
import glob


class MetabVipAgent(Agent):
    """
    代谢集代谢物聚类分析
    """

    def __init__(self, parent):
        super(MetabVipAgent, self).__init__(parent)
        options = [
            {'name': 'diff_dir', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},  # 差异分析结果目录
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_trans", "type": "infile", "format": "sequence.profile_table"},  # id转化使用
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            #{'name': 'group_name', 'type': 'string', 'default': ''},  # 差异分组 “P_vs_C;C_vs_X”
            {'name': 'group_method', 'type': 'string', 'default': 'none'},  # 分组计算方法
            {'name': 'vip_type', 'type': 'string', 'default': 'oplsda'},  # vip来源 plsda，oplsda
            {'name': 'vip_cut', 'type': 'float', 'default': 1.0},  # vip 阈值
            {'name': 'vip_top', 'type': 'int', 'default': 30},  # Top vip 数目
            {'name': 'mct', 'type': 'string', 'default': ''},  # 代谢物聚类算法，hierarchy、kmeans、无
            {'name': 'mcd', 'type': 'string', 'default': ''},  # 代谢物距离计算方式
            {'name': 'n_cluster', 'type': 'int', 'default': 0},  # 代谢物聚类数目，kmeans时使用
            {'name': 'mcm', 'type': 'string', 'default': ''},  # 代谢物聚类方式, hierarchy时使用，"complete","average","single"
            ## 样本相关参数暂备用
            {'name': 'sct', 'type': 'string', 'default': ''},  # 样本聚类算法，hierarchy、无
            {'name': 'scd', 'type': 'string', 'default': ''},  # 样本距离计算方式
            {'name': 'scm', 'type': 'string', 'default': ''},  # 样本聚类方式,"complete","average","single"
            {"name": "scale", "type": "bool", "default": False},  # 标准化使用
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"}, # 标准化时用，用原始丰度表
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("diff_dir").is_set:
            raise OptionError("请传入差异分析结果目录！", code="34701101")
        if self.option("mct") == "kmeans":
            if not self.option("n_cluster") > 0:
                raise OptionError("使用kmeans时必须聚类数目必须大于1！", code="34701102")
        if self.option("mct") == "hierarchy":
            if not self.option("mcm"):
                raise OptionError("代谢物使用hierarchy时必须传入层次聚类方式！", code="34701103")
        if self.option("mct") == "hierarchy" or self.option("mct") == "kmeans":
            if not self.option("mcd"):
                raise OptionError("必须输入代谢物距离计算方式！", code="34701104")
        if self.option("group_method") != "none":
            if not self.option("group").is_set:
                raise OptionError("有分组计算方法时必须输入分组文件！", code="34701105")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '8G'

    def end(self):
        super(MetabVipAgent, self).end()


class MetabVipTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MetabVipTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.script = self.config.PACKAGE_DIR + "/metabolome/scripts/vip.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.group = []

    def run(self):
        """
        运行
        """
        super(MetabVipTool, self).run()
        self.logger.info("开始运行命令！")
        self.run_vip()
        self.set_output()

    def run_vip(self):
        """
        """
        diff_dir = self.option("diff_dir").prop["path"]
        self.logger.info(diff_dir)
        vip_type = self.option("vip_type")
        vip_cut = self.option("vip_cut")
        vip_top = self.option("vip_top")
        diff_files = os.listdir(diff_dir)
        number = 0
        for each_diff in diff_files:
            if ".diff" not in each_diff:
                continue
            number += 1
            group = each_diff.split(".diff")[0]
            name = "vip_" + str(number)
            each_diff = os.path.join(diff_dir, each_diff)
            ### 兼容新旧结果
            infile = each_diff
            if not "tmp" in diff_dir:
                self.metab_to_id_trans = Relation()
                table_path = self.option("metab_trans").prop["path"]
                map_table, id_name_dict = self.metab_to_id_trans.get_dataframe_and_dict(table_path, mode=2)
                infile = os.path.join(self.work_dir, "tmp_diff_" + group)
                self.trans_metab(each_diff, infile, map_table, id_name_dict)
            output = os.path.join(self.work_dir, group)
            self.group.append(group)
            if not os.path.exists(output):
                os.mkdir(output)
            cmd = '{} {} -exp {} -out {} --nsc -vty {} -vcut {} -vtop {}'.format(self.python_path, self.script, infile,
                                                                                 output, vip_type, vip_cut, vip_top)
            cmd += " -dn " + group.replace("_vs_", ",")
            if self.option("mct") != "none":
                cmd += " -gct " + self.option("mct")
            else:
                cmd += " --ngc"
            if self.option("mcd"):
                cmd += " -gcd " + self.option("mcd")
            if self.option("mct") == "kmeans":
                cmd += " -n_cluster " + str(self.option("n_cluster"))
            if self.option("mct") == "hierarchy":
                cmd += " -gcm " + self.option("mcm")
            if self.option("metab_set_table").is_set:
                select_file = self.option("metab_set_table").prop["path"]
                cmd += " -metaf " + select_file
            if self.option("group_method") != "none":
                cmd += " -group " + self.option("group").prop["path"] + " -gm " + self.option("group_method")
            if "scale" in self.get_option_object().keys() and self.option("scale"):
                cmd += " --scale -t_abu " + self.option("metab_table").prop["path"]
            command = self.add_command(name, cmd)
            command.run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("metab vip succeed")
            else:
                self.set_error("metab vip failed!", code="34701101")
                raise Exception("metab vip failed!")



    def set_output(self):
        self.logger.info("进行id转化...")
        table_path = self.option("metab_trans").prop["path"]
        diff_dir = self.option("diff_dir").prop["path"]
        self.metab_trans = Relation()
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
        for group in self.group:
            od_dir = self.work_dir + "/" + group
            new_dir = self.output_dir + "/" + group
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            if self.option("mct") == "hierarchy":
                old_tree = od_dir + "/row.cluster_tree.txt"
                new_tree = new_dir + "/metab.cluster_tree.xls"
                if os.path.exists(old_tree):
                    self.trans_tree_name(id_name_dict, old_tree, new_tree)
                    id_tree = new_dir + "/metab_id.cluster_tree.xls"
                    self.link_file(old_tree, id_tree)
            elif self.option("mct") == "kmeans":
                subclusters = glob.glob(od_dir + '/*subcluster*')
                for each in subclusters:
                    #each = eachfile.split("/")[-1]
                    newname = each.replace("row", "metab").replace(od_dir, new_dir)
                    self.trans_data(each, newname, map_table, id_name_dict)
                if os.path.exists(od_dir + "/row.kmeans_cluster.txt"):
                    #self.link_file(od_dir + "/row.kmeans_cluster.txt", new_dir + "/metab.kmeans_cluster.xls")
                    self.metab_trans.trans_col_data(id_name_dict, od_dir + "/row.kmeans_cluster.txt",
                                                    new_dir + "/metab.kmeans_cluster.xls")
            old_file = od_dir + "/Vip.xls"
            self.logger.info(old_file)
            new_file = new_dir + "/Vip_exp.xls"
            if "scale" in self.get_option_object().keys() and self.option("scale"):
                new_file_scale = new_dir + "/Vip_scale_exp.xls"
                select_outfile = od_dir + "/Vip_before_scale.xls"
                self.trans_data(select_outfile, new_file, map_table, id_name_dict)
                self.trans_data(old_file, new_file_scale, map_table, id_name_dict)
                self.metab_trans.add_oid_fun(new_file, table_path,link_k='Metabolite') ##20190617
                self.metab_trans.add_oid_fun(new_file_scale, table_path,link_k='Metabolite')  #20190617
            else:
                self.trans_data(old_file, new_file, map_table, id_name_dict)
                self.metab_trans.add_oid_fun(new_file, table_path,link_k='Metabolite') ##20190617
        self.end()

    def link_file(self, oldfile, newfile):
        #oldfile = os.path.join(self.work_dir, oldfile)
        #newfile = os.path.join(self.output_dir, newfile)
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)

    def trans_data(self, oldfile, newfile, map_table, id_name_dict):
        try:
            self.metab_trans.get_trans_file(map_table, id_name_dict, oldfile, newfile)
        except Exception as e:
            self.set_error("id转化异常——%s", variables=(e), code="34701103")
            raise Exception("id转化异常——{}".format(e))

    def trans_metab(self, oldfile, newfile, map_table, id_name_dict):
        try:
            self.metab_to_id_trans.get_metab_to_idfile(map_table, id_name_dict, oldfile, newfile)
        except Exception as e:
            self.set_error("代谢物转化异常——%s", variables=(e), code="34701105")
            raise Exception("代谢物转化异常——%s", e)

    def trans_tree_name(self, namedict, treefile, newfile):
        """
        用正则的方式替换复原树文件中的名称
        """
        if not isinstance(namedict, dict):
            self.set_error('复原树的枝名称需要旧名称和当前名称的字典', code="34701107")
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
            self.set_error('聚类树文件无法找到或者无法打开：%s', variables=(e), code="34701108")


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "MetabVip",
            "type": "tool",
            "name": "metabolome.metabset.metab_vip",
            "instant": True,
            "options": dict(
                diff_dir="/mnt/ilustre/users/sanger-dev/workspace/20180627/Metabolome_metabolome/DiffPls/MergeDiff/output/DiffStat/",
                vip_type="oplsda",
                vip_cut=1.5,
                vip_top=30,
                mct="hierarchy",
                mcd="euclidean",
                mcm="complete"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
