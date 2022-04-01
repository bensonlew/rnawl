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

class DiffMulStatAgent(Agent):
    """
    两组差异检验，student T检验，welch T检验，wilcox秩和检验
    """

    def __init__(self, parent):
        super(DiffMulStatAgent, self).__init__(parent)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'metabolome.express,metabolome.metab_abun'},  # 表达矩阵文件
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {'name': 'group_name', 'type': 'string', 'default': ''},  # 差异组名，eg："A|B";"A|C" ,两两比较
            {'name': 'mul_type', 'type': 'string', 'default': 'pca'},  # 多元统计类型，pca，plsda, oplsda ，可以是“pca;plsda”形式
            {'name': 'confidence', 'type': 'string', 'default': '0.95'},  # 置信度，与mul_type对应的“0.95;0.95”形式
            {'name': 'perm', 'type': 'string', 'default': ''},  # 置换次数，与mul_type对应的“200;100”形式
            {'name': 'data_trans', 'type': 'string', 'default': 'UV'},  # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数
            {'name': 'id_trans', 'type': 'bool', 'default': False},  # 代谢id是否需要转化为代谢物名称
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp_file").is_set:
            raise OptionError("请传入丰度矩阵！", code="34700101")
        if not self.option("group_file").is_set:
            raise OptionError("请传入group文件！", code="34700102")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        super(DiffMulStatAgent, self).end()


class DiffMulStatTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(DiffMulStatTool, self).__init__(config)
        #self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.set_environ(R_LIBS=self.config.SOFTWARE_DIR + '/program/R-3.3.1/lib64/R/library')
        self.r_path = "/program/R-3.3.1/bin/Rscript"
        self.cmds = []
        self.results = []

    def run(self):
        """
        运行
        """
        super(DiffMulStatTool, self).run()
        self.logger.info("开始运行命令！")
        self.logger.info(self.option("exp_file").prop["path"])
        self.logger.info(self.option("group_file").prop["path"])
        self.group_file_sorted = self.sort_group_file_fun()  #20200318 对同组内样本做排序，oplsda和plsda 对组内样本顺序有要求，排序以保证结果一致性
        self.check_exp_file()
        if self.option("group_name") != "":
            self.diff_g_paths = self.make_group(self.group_file_sorted, self.option("group_name"))
        else:
            self.diff_g_paths = [self.group_file_sorted]
        self.run_mul_stat()
        if self.option("id_trans"):
            self.change_id()
        self.set_output()

    def sort_group_file_fun(self):
        group_file = self.option("group_file").prop["path"]
        group_map_sample = dict()
        with open(group_file) as f:
            f.readline()
            group_list = list()
            for line in f:
                line = line.strip()
                if line=='':
                    continue
                spline = line.split("\t")
                if spline[1] not in group_list:
                    group_list.append(spline[1])
                    group_map_sample[spline[1]] = list()
                group_map_sample[spline[1]].append(spline[0])
        group_file_sorted = 'sorted_group.xls'
        with open(group_file_sorted, 'w') as fw:
            fw.write("#sample\tgroup_name\n")
            for g in sorted(group_list):
                fw.write('\n'.join([i+'\t'+ g for i in sorted(group_map_sample[g])]) + '\n')
        return group_file_sorted

    def check_exp_file(self):
        exp_file = self.option("exp_file").prop["path"]
        data = pd.read_table(exp_file,sep='\t',header=0)
        data.drop_duplicates(['metab_id'],keep='first',inplace=True)
        self.new_exp = self.work_dir+'/new_exp.xls'
        data.to_csv(self.new_exp,sep='\t',index=False)

    def make_group(self, group_file, diff_group):
        """
        根据差异分组名，从总group表中生成各差异组group
        :param group_file: 总group文件
        :param diff_group: 差异组名， eg："A|B";"A|C"
        :return: 返回生成差异group文件路径diff_g_paths列表
        """
        diff_groups = diff_group.split(";")
        group_dict = {}
        diff_g_paths = []
        with open(group_file, "r") as f:
            for line in f:
                line = line.strip()
                line1 = line.split("\t")
                if not "#" in line:
                    group = line1[1]
                    name = line1[0]
                    if not group_dict.has_key(group):
                        group_dict[group] = []
                    group_dict[group].append(name)
        for each in diff_groups:
            two_groups = each.split("_vs_")
            #diff_name = each.replace("|", "_vs_")
            diff_name = each
            diff_g_path = self.work_dir + "/" + diff_name + "_group"
            diff_g_paths.append(diff_g_path)
            with open(diff_g_path, "w") as gf:
                gf.write("#sample\tgroup\n")
                for i in two_groups:
                    if not group_dict.has_key(i):
                        raise OptionError('分组样本-%s-不在丰度表所拥有的样本内，请检查分组方案', variables=(i), code="34700103")
                    names = group_dict[i]
                    for j in names:
                        gf.write(j + "\t" + i + "\n")
            self.logger.info("生成{}_group表成功".format(each))
        return diff_g_paths

    def run_mul_stat(self):
        """
        pca、plsda、oplsda分析
        """
        count = 0
        for myfile in self.diff_g_paths:
            group = myfile.split("/")[-1].rpartition("_group")[0]
            exp = self.new_exp
            if not self.option("id_trans"):
                out = self.output_dir
            else:
                out = self.work_dir
            if self.option("group_name"):
                out = out + "/" + group
            if not os.path.exists(out):
                os.mkdir(out)
            ci = self.option("confidence")
            data_trans = self.option("data_trans")
            tmp_trans = data_trans.split(";")
            tmp_list = []
            for x in tmp_trans:
                if x == "":
                    tmp_list.append("none")
                else:
                    tmp_list.append(x)
            data_trans = ";".join(tmp_list)
            mul_type = self.option("mul_type")
            perm = self.option("perm")
            diffstat = MulDiffStat()
            diffstat.mul_pca(exp, out, myfile, ci=ci, data_trans=data_trans, mul_type=mul_type, perm=perm)
            mul_type_new = mul_type.replace(";","_")
            rname = " run_{}_{}.r".format(mul_type_new, group)
            cmd = self.r_path + rname
            self.cmds.append(cmd)
            count += 1
            command = self.add_command("diff" + str(count), cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run {} succeed".format(cmd))
            else:
                self.set_error("run %s failed", variables=(cmd), code="34700101")
                raise Exception("run {} failed".format(cmd))
        if count == len(self.cmds):
            self.logger.info("全部命令运行完成！！！")
        else:
            self.set_error("cmd未全部完成", code="34700103")
            raise Exception("cmd未全部完成")

    def set_output(self):
        """
        id 转化
        """
        self.logger.info("进行id转化...")
        pls_dir = self.output_dir
        table_path = self.option("metab_desc").prop["path"]
        self.metab_trans = Relation()
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path)
        groups = os.listdir(self.output_dir)
        if self.option("group_name") != "":
            for eachgroup in groups:
                load_files = glob.glob(pls_dir + "/" + eachgroup + '/*loading.xls')
                vip_files = glob.glob(pls_dir + "/" + eachgroup + '/*vip.xls')
                for eachfile in load_files:
                    newfile = eachfile.replace("loading","loadings")
                    self.trans_data(eachfile, newfile, map_table, id_name_dict)
                for eachfile in vip_files:
                    newfile = eachfile.replace("vip.xls","vips.xls")
                    self.trans_data(eachfile, newfile, map_table, id_name_dict)
                #转splot图的数据 ，目前只有oplsda 画splot图
                splot_files =  glob.glob(pls_dir + "/" + eachgroup + '/OPLS-DA.splot.xls')
                for eachfile in splot_files:
                    newfile = eachfile.replace("splot.xls","splots.xls")
                    tmp_data = pd.read_table(eachfile,sep='\t')
                    tmp_data = tmp_data.reset_index()
                    tmp_data.columns = ['metab id','p','cor_value']
                    tmp_data.to_csv(eachfile,sep='\t', index=False)
                    self.trans_data(eachfile, newfile, map_table, id_name_dict)
        else:
            load_files = glob.glob(pls_dir + '/*loading.xls')
            for eachfile in load_files:
                newfile = eachfile.replace("loading","loadings")
                self.trans_data(eachfile, newfile, map_table, id_name_dict)
        self.end()

    def trans_data(self, oldfile, newfile, map_table, id_name_dict):
        self.metab_trans.get_trans_file(map_table, id_name_dict, oldfile, newfile)
