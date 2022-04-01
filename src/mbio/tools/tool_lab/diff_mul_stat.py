# -*- coding: utf-8 -*-

import os, re, subprocess, glob
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.scripts.mul_diff_stat import MulDiffStat


class DiffMulStatAgent(Agent):

    def __init__(self, parent):
        super(DiffMulStatAgent, self).__init__(parent)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'tool_lab.table'},  # 表达矩阵文件
            {"name": "groups", "type": "string"},
            {"name": "group_file", "type": "infile", "format": "tool_lab.simple"},  # 分组文件
            {'name': 'mul_type', 'type': 'string', 'default': 'pca'},  # 多元统计类型，pca，plsda, oplsda ，可以是“pca;plsda”形式
            {'name': 'confidence', 'type': 'string', 'default': '0.95'},  # 置信度，与mul_type对应的“0.95;0.95”形式
            {'name': 'perm', 'type': 'string', 'default': ''},  # 置换次数，与mul_type对应的“200;100”形式
            {'name': 'data_trans', 'type': 'string', 'default': ''},  # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数

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
        self.diff_g_paths = [self.group_file_sorted]
        self.run_mul_stat()
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
            for g in group_list:
                if g in self.option("groups").split(','):
                    fw.write('\n'.join([i+'\t'+ g for i in sorted(group_map_sample[g])]) + '\n')
        return group_file_sorted

    def check_exp_file(self):
        exp_file = self.option("exp_file").prop["path"]
        data = pd.read_table(exp_file,sep='\t',header=0)
        col_1 = data.columns[0]
        data.drop_duplicates(col_1,keep='first',inplace=True)
        self.new_exp = self.work_dir+'/new_exp.xls'
        data.to_csv(self.new_exp,sep='\t',index=False)



    def run_mul_stat(self):
        """
        pca、plsda、oplsda分析
        """
        count = 0
        for myfile in self.diff_g_paths:
            group = myfile.split("/")[-1].split("_group")[0]
            exp = self.new_exp
            out = self.output_dir
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
        self.end()

