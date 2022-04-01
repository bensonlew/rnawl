# coding=utf-8
#__author__ = 'shaohua.yuan'

import os,re
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from mbio.packages.metabolome.scripts.merge_table import merge_table, merge_mul_metabset
from mbio.packages.metabolome.scripts.profile_select import profile_select
import pandas as pd
import copy
import numpy as np


class CreatTableAgent(Agent):
    """
    代谢项目时workflow中使用阴阳离子合并list
    """
    def __init__(self, parent):
        super(CreatTableAgent, self).__init__(parent)
        options = [
            {'name': 'diff_dir1', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},
            {'name': 'diff_dir2', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},
            {'name': 'exp', 'type': 'infile', 'format': 'metabolome.metab_abun'},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {'name': 'out_venn', 'type': 'outfile', 'format': 'metabolome.mul_metabset,sequence.profile_table'},
            {'name': 'out_abu', 'type': 'outfile', 'format': 'metabolome.express,sequence.profile_table'},
            {'name': 'set_list', 'type': 'outfile', 'format': 'metabolome.metabset,sequence.profile_table'},
            {'name': 'diff_mul_set', 'type': 'outfile', 'format': 'metabolome.mul_metabset,sequence.profile_table'},
            {'name': 'top', 'type': 'int'},
            {'name': 'filter', 'type': 'bool', 'default': True},  # 是否去除pos，neg名字的代谢物
            {"name": "scale", "type": "bool", "default": True},  ## 是否标准化聚类数据
            {'name': 'scale_abu', 'type': 'outfile', 'format': 'metabolome.express,sequence.profile_table'},  # 标准化时使用
            {'name': 'two_sample_diff_dir1', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},
            {'name': 'two_sample_diff_dir2', 'type': 'infile', 'format': 'annotation.mg_anno_dir'}
        ]
        self.add_option(options)
        self.used_metab_select = ""

    def check_options(self):
        if self.option("top"):
            if not int(self.option("top")) > 0:
                raise OptionError("top必须大于0的整数", code="34701701")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('5')

    def end(self):
        super(CreatTableAgent, self).end()


class CreatTableTool(Tool):
    """
    LC项目时workflow中使用阴阳离子合并list
    """
    def __init__(self, config):
        super(CreatTableTool, self).__init__(config)

    ## 合并module的 mix  （一列）
    def run_merge(self,diff_intersection1,diff_intersection2, outfile):
        ##diff_intersection1 = os.path.join(self.option("diff_dir1").prop["path"], "Diff_intersection.metabset.xls")
        ##diff_intersection2 = os.path.join(self.option("diff_dir2").prop["path"], "Diff_intersection.metabset.xls")

        #outfile = os.path.join(self.work_dir, "tmp_merge_diff_intersection.metabset.xls")
        try:
            merge_file = merge_table(str(diff_intersection1), str(diff_intersection2), "metab_id", outfile=outfile,
                                     how="outer")
        except Exception as e:
            self.set_error("创建合并非冗余代谢集失败——%s", variables=(e), code="34701701")

    #2个diff_dir  合并各差异组代谢集  如果filter,根据代谢物名称过滤
    def run_merge_mulset(self, mul_metabset1, mul_metabset2, outfile):
        #mul_metabset1 = os.path.join(self.option("diff_dir1").prop["path"], "mul.metabset.list.xls")
        #mul_metabset2 = os.path.join(self.option("diff_dir2").prop["path"], "mul.metabset.list.xls")
        #outfile = os.path.join(self.output_dir, "merge_mul.metabset.xls")
        if self.option("filter"):
            exp_file = self.option("exp").prop["path"]
            exp_des_file = exp_file.replace("metab_abund.txt", "metab_desc.txt")
            if not os.path.exists(exp_des_file):
                self.set_error("文件不存在——%s", variables=(exp_des_file), code="34701703")
                raise Exception("文件不存在——{}".format(exp_des_file))
        else:
            exp_des_file = None
        try:
            merge_mul_metabset(mul_metabset1, outfile, mul_file2=mul_metabset2, exp_des_file=exp_des_file)
        except Exception as e:
            self.set_error("创建合并venn用非冗余代谢集失败——%s", variables=(e), code="34701704")

    ##只有1个diff_dir1
    def run_mulset(self, mul_metabset1, outfile):
        #mul_metabset1 = os.path.join(self.option("diff_dir1").prop["path"], "mul.metabset.list.xls")
        #outfile = os.path.join(self.output_dir, "merge_mul.metabset.xls")
        if self.option("filter"):
            exp_file = self.option("exp").prop["path"]
            exp_des_file = exp_file.replace("metab_abund.txt", "metab_desc.txt")
            if not os.path.exists(exp_des_file):
                self.set_error("文件不存在——%s", variables=(exp_des_file), code="34701706")
                raise Exception("文件不存在——{}".format(exp_des_file))
        else:
            exp_des_file = None
        try:
            merge_mul_metabset(mul_metabset1, outfile, exp_des_file=exp_des_file)
        except Exception as e:
            self.set_error("创建合并venn用非冗余代谢集失败——%s", variables=(e), code="34701707")

    def run_remove_head(self):
        out_list = os.path.join(self.output_dir, "merge_diff_intersection.metabset.xls")
        set_table = pd.read_table(self.used_metab_select, sep="\t", header=0, index_col=0)
        set_table.to_csv(out_list, sep='\t', header=False, index=True)
        set_table_list = set_table.index.tolist()
        set_table_str = ",".join(set_table_list)
        out_diff_mul = os.path.join(self.output_dir, "diff_set_mul.metabset.xls")  #self.used_metab_select
        with open(out_diff_mul, "w") as outf:
            outf.write("DiffSet_mix\t" + set_table_str + "\n")

    #根据 代谢物，top，样本，筛选exp表
    def select_table(self):
        """
        """
        select = profile_select()
        profile = self.option("exp").prop["path"]
        outfile = os.path.join(self.output_dir, "merge_diff_intersection.metabset.abu.xls")
        select_column = "metab_id"
        select_file = self.used_metab_select
        samples = self.option("group_file").prop["sample"]
        sams = ",".join(samples)
        if self.option("top"):
            top = int(self.option("top"))
        else:
            top = 0
        self.logger.info(select_file)
        if "scale" in self.get_option_object().keys() and self.option("scale"):
            exp_profile = self.scale_data(profile)
            outfile_scale = os.path.join(self.output_dir, "merge_diff_intersection.metabset.abu_scale.xls")
            select.run_select(exp_profile, outfile_scale, 0, select_column=select_column, sam=sams,
                              gene_list=select_file, total="F", top=top)
            table = pd.read_table(outfile_scale, sep='\t', index_col=0)
            origin_table = pd.read_table(profile, sep='\t', index_col=0)
            select_names = table.index
            self.logger.info(profile)
            self.logger.info(outfile_scale)
            selecl_origin = origin_table.loc[select_names,].reindex()
            selecl_origin.to_csv(outfile, index=True, header=True, sep="\t")
        else:
            select.run_select(profile, outfile, 0, select_column=select_column, sam=sams, gene_list=select_file,
                              total="F", top=top)

    def scale_data(self, ori_table):
        table = pd.read_table(ori_table, sep="\t", index_col=0)
        scaled_data = table.apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)
        exp_profile = os.path.join(self.output_dir, "scale_data.xls")
        scaled_data.to_csv(exp_profile, index=True, header=True, sep="\t")
        return exp_profile

    def set_output(self):
        out_abu = os.path.join(self.output_dir, "merge_diff_intersection.metabset.abu.xls")
        out_list = os.path.join(self.output_dir, "merge_diff_intersection.metabset.xls")  #单列的mix 无表头
        out_diff_mul = os.path.join(self.output_dir, "diff_set_mul.metabset.xls")  #2列的mix
        out_venn = os.path.join(self.output_dir, "merge_mul.metabset.xls")  #各比较组的代谢集
        outfile_scale = os.path.join(self.output_dir, "merge_diff_intersection.metabset.abu_scale.xls")
        try:
            self.option('out_abu', out_abu)
            self.option('out_venn', out_venn)
            self.option("set_list", out_list)
            self.option('diff_mul_set', out_diff_mul)
            if "scale" in self.get_option_object().keys() and self.option("scale"):
                self.option('scale_abu', outfile_scale)
            self.logger.info("设置输出结果文件成功")
        except Exception as e:
            self.set_error("输出结果文件异常——%s", variables=(e), code="34701709")
            raise Exception("输出结果文件异常——{}".format(e))

    def run(self):
        super(CreatTableTool, self).run()
        in_dirs = []
        ##Diff_intersection.metabset.xls    #tmp_merge_diff_intersection.metabset.xl
        ##mul.metabset.list.xls

        for  dir in ['diff_dir1', 'diff_dir2', 'two_sample_diff_dir1','two_sample_diff_dir2']:
            if self.option(dir).is_set:
                in_dirs.append(dir)

        if len(in_dirs)==1:
            self.used_metab_select = os.path.join(self.option(in_dirs[0]).prop["path"],"Diff_intersection.metabset.xls")
            mul_metabset1 = os.path.join(self.option(in_dirs[0]).prop["path"], "mul.metabset.list.xls")
            outfile = os.path.join(self.output_dir, "merge_mul.metabset.xls")
            self.run_mulset(mul_metabset1, outfile)
        else:
            dir_1 = in_dirs[0]
            dir_1_mix = os.path.join(self.option(dir_1).prop["path"],"Diff_intersection.metabset.xls")
            dir_1_mul =  os.path.join(self.option(dir_1).prop["path"],"mul.metabset.list.xls")
            for id, o_dir in enumerate(in_dirs[1:],1):
                dir_2_mix = os.path.join(self.option(o_dir).prop["path"],"Diff_intersection.metabset.xls")
                dir_2_mul = os.path.join(self.option(o_dir).prop["path"],"mul.metabset.list.xls")
                out_mix = self.work_dir+'/tmp_merge_mix.'+str(id)
                out_mul = self.work_dir+'/tmp_merge_mul.'+str(id)
                self.run_merge(dir_1_mix, dir_2_mix, out_mix)
                self.run_merge_mulset(dir_1_mul,dir_2_mul, out_mul)
                dir_1_mix = out_mix
                dir_1_mul = out_mul
            os.rename(out_mix, self.work_dir+'/tmp_merge_diff_intersection.metabset.xls')
            self.used_metab_select = self.work_dir+'/tmp_merge_diff_intersection.metabset.xls'
            if os.path.exists(self.output_dir+'/merge_mul.metabset.xls'):
                os.remove(self.output_dir+'/merge_mul.metabset.xls')
            os.link(out_mul,self.output_dir+'/merge_mul.metabset.xls')

        result = self.read_metab_select()
        if result:
            if self.option("filter"):
                self.filter_matabs()
            self.select_table()
            self.run_remove_head()
        else:
            self.set_empty_file()
        self.set_output()
        self.end()

    def read_metab_select(self):
        table = pd.read_table(self.used_metab_select, sep="\t", header=0)
        if len(table) < 1:
            result = False
        else:
            result = True
        return result

    #设置空文件，丰度文件按top筛选
    def set_empty_file(self):
        out_abu = os.path.join(self.output_dir, "merge_diff_intersection.metabset.abu.xls")
        out_list = os.path.join(self.output_dir, "merge_diff_intersection.metabset.xls")
        out_diff_mul = os.path.join(self.output_dir, "diff_set_mul.metabset.xls")
        df_empty = pd.DataFrame(columns=["metab"])
        #df_empty.to_csv(out_abu,sep="\t",index=False,header=True)
        df_empty.to_csv(out_list, sep="\t", index=False, header=True)
        df_empty.to_csv(out_diff_mul, sep="\t", index=False, header=True)
        profile = self.option("exp").prop["path"]
        table = pd.read_table(profile, sep="\t", header=0, index_col=0)
        if self.option("top"):
            top = int(self.option("top"))
        else:
            top = len(table)
        columns_name = table.columns
        self.logger.info(profile)
        table['Total'] = table.loc[:, columns_name].apply(lambda x: x.sum(), axis=1)
        table = table.sort_values(["Total"], ascending=False)[0:top]
        table = table.drop('Total', 1)
        table.to_csv(out_abu, index=True, header=True, sep="\t")


    def filter_matabs(self):
        """
        筛除neg，pos无名称的代谢物
        """
        exp_file = self.option("exp").prop["path"]
        exp_des_file = exp_file.replace("metab_abund.txt", "metab_desc.txt")
        if not os.path.exists(exp_des_file):
            self.set_error("文件不存在——%s", variables=(exp_des_file), code="34701711")
            raise Exception("文件不存在——{}".format(exp_des_file))
        # if self.option("diff_dir2").is_set:
        #     tmp_metabset = os.path.join(self.work_dir, "tmp_merge_diff_intersection.metabset.xls")  #
        # else:
        #     tmp_metabset = os.path.join(self.option("diff_dir1").prop["path"], "Diff_intersection.metabset.xls")
        tmp_metabset = self.used_metab_select
        exp_table = pd.read_table(tmp_metabset, sep="\t", header=0)
        exp_table = exp_table[["metab_id"]]
        exp_des = pd.read_table(exp_des_file, sep="\t", header=0)
        exp_des = exp_des[["metab_id", "Metabolite"]]
        merge_table = pd.merge(exp_table, exp_des, how="inner", on="metab_id")
        merge_table.index = merge_table["Metabolite"]
        exp_des_list = exp_des["Metabolite"].values.tolist()
        metab_list = copy.copy(exp_des_list)
        for each in exp_des_list:
            lower_each = each.lower()
            '''
            if lower_each.startswith('pos') or lower_each.startswith('neg'):
            '''
            m1 = re.match(r'^pos\d+$', lower_each) # modify 20181115 ,代谢物可能带pos前缀，eg. postin
            m2 = re.match(r'^neg\d+$', lower_each)
            if m1 or m2:
                metab_list.remove(each)
        filter_table = merge_table.loc[metab_list].reindex()
        filter_table = filter_table.dropna()
        filter_table = filter_table[["metab_id"]]
        outfile = os.path.join(self.work_dir, "remove_noname.metabset.xls")
        self.used_metab_select = outfile
        filter_table.to_csv(outfile, sep="\t", index=False, quoting=3)
        return outfile


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
            "id": "CreatTable",
            "type": "tool",
            "name": "metabolome.creat_table",
            "instant": True,
            "options": dict(
                diff_dir1="/mnt/ilustre/users/sanger-dev/workspace/20181115/Metabolome_tsg_32888/DiffPls/MergeDiff/output/Metabset/",
                #diff_dir2="/mnt/ilustre/users/sanger-dev/workspace/20180803/Metabolome_tsg_31353/DiffPls1/MergeDiff/output/Metabset/",
                exp="/mnt/ilustre/users/sanger-dev/workspace/20181115/Metabolome_tsg_32888/output/Preprocess/pos/metab_abund.txt",
                group_file="/mnt/ilustre/users/sanger-dev/workspace/20181115/Metabolome_tsg_32888/noQC_group.xls",
                filter=True,
                #top=50
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
