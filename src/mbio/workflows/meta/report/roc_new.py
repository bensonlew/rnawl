# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

"""个性化ROC分析"""

import os
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import pandas as pd


class RocNewWorkflow(Workflow):
    """
    涉及lefse分析、两组比较分析和随机森林分析以及roc分析
    version v1.0
    author: zhouxuan
    last_modify: 2017.05.12
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RocNewWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "lefse_analysis", "type": "bool", 'default': False},  # 判断是否做这三个分析
            {"name": "two_group_analysis", "type": "bool", 'default': False},
            {"name": "ran_for_analysis", "type": "bool", 'default': False},
            {"name": "otu_id", "type": "string", "default": ""},  # 导表使用

            {"name": "lefse_otu", "type": "infile", 'format': "meta.otu.otu_table"},  # lefse 分析的参数设置
            {"name": "lefse_group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "second_group_detail", "type": "string"},
            {"name": "lda_filter", "type": "float", "default": 2.0},
            {"name": "strict", "type": "int", "default": 0},
            {"name": "group_name", "type": "string"},
            {"name": "start_level", "type": "int", "default": 9},
            {"name": "end_level", "type": "int", "default": 9},

            {"name": "two_ran_otu", "type": "infile", 'format': "meta.otu.otu_table"},  # 两组和随机森林的参数设置
            {"name": "two_ran_group", "type": "infile", 'format': "meta.otu.group_table"},
            {"name": "level", "type": "int"},
            {"name": "method_cal", "type": "string", "default": "student"},
            {"name": "ci", "type": "float", "default": 0.99},
            {"name": "q_test", "type": "string", "default": "fdr"},
            {"name": "tree_number", "type": "int", "default": 500},

            {"name": "lefse_cho", "type": "string", "default": "p-value"},  # 筛选条件设置
            {"name": "lefse_num", "type": "int", "default": 50},
            {"name": "two_group_cho", "type": "string", "default": "p-value"},
            {"name": "two_group_num", "type": "int", "default": 50},
            {"name": "Ran_for_num", "type": "int", "default": 50},
            {"name": "intersection", "type": "bool", "default": False},

            {"name": "roc_calc_method", "type": "string", "default": "sum"},  # ROC计算方法以及env表的计算方法
            {"name": "roc_method_1", "type": "string", "default": ""},  # a:0.1
            # 后续需要从网页端传入,必须携带相应的分组名称
            {"name": "roc_method_2", "type": "string", "default": ""},  # b:0.2
            {"name": "env_table", "type": "infile", 'format': "meta.otu.group_table"},  # env表格存在有或没有的状态
            {"name": "env_labs", "type": "string", "default": ""},
            {"name": "env_id", "type": "string", "default": ""},

            {"name": "update_info", "type": "string"},  # workflow更新
            {"name": "main_id", "type": "string", "default": "56a847100e6da953f3ce31e6"},  # 主表id 测试用default
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.two_ran_lefse = self.add_tool('meta.beta_diversity.roc_new')
        self.roc_lefse = self.add_tool('meta.beta_diversity.roc')
        self.roc_ran = self.add_tool('meta.beta_diversity.roc')
        self.roc_two = self.add_tool('meta.beta_diversity.roc')
        self.roc_intersection = self.add_tool('meta.beta_diversity.roc')
        self.roc_env = self.add_tool('meta.beta_diversity.roc')
        self.dict = {}

    def check_options(self):
        if self.option('method_cal') not in ['student', 'welch', 'wilcox']:
            raise OptionError('错误的两组差异检验方法：%s', variables=(self.option('method_cal')), code="12704001")
        if self.option('q_test') not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            raise OptionError('错误的多重检验矫正方法：%s', variables=(self.option('q_test')), code="12704002")
        if self.option('lefse_cho') not in ['p-value', 'LDA']:
            raise OptionError('错误的lefse筛选方法：%s', variables=(self.option('lefse_cho')), code="12704003")
        if self.option('two_group_cho') not in ['p-value', 'Corrected-p-value']:
            raise OptionError('错误的两组差异比较结果筛选方法：%s', variables=(self.option('two_group_cho')), code="12704004")
        if self.option('roc_calc_method') not in ['sum', 'average', 'median', 'MI']:
            raise OptionError('错误的roc计算方法：%s', variables=(self.option('roc_calc_method')), code="12704005")
        if self.option('roc_calc_method') == "MI" and (self.option('roc_method_1') == "" or self.option('roc_method_2') == ""):
            raise OptionError('错误的roc计算方法参数：%s %s', variables=(self.option('roc_method_1'),
                                                               self.option('roc_method_2')), code="12704006")

    def run(self):
        self.run_two_ran_lefse()
        if self.option('roc_calc_method') == "MI":
            value_1 = self.option("roc_method_1").split(":")
            self.dict[value_1[0]] = value_1[1]
            value_2 = self.option("roc_method_2").split(":")
            self.dict[value_2[0]] = value_2[1]
        super(RocNewWorkflow, self).run()

    def run_two_ran_lefse(self):
        if self.option('intersection'):
            start_level = self.option("level")
            end_level = self.option("level")
        else:
            start_level = self.option("start_level")
            end_level = self.option("end_level")
        new_table = self.change_otuname(self.option("two_ran_otu").prop["path"])
        options = {
            "lefse_input": self.option("lefse_otu"),
            "lefse_group": self.option("lefse_group"),
            "otutable": new_table,
            # "otutable": self.option("two_ran_otu"),
            "grouptable": self.option("two_ran_group"),
            "lda_filter": self.option("lda_filter"),
            "strict": self.option("strict"),
            "lefse_gname": self.option("group_name"),
            "start_level": start_level,
            "end_level": end_level,
            "method_cal": self.option("method_cal"),
            "ci": self.option("ci"),
            "q_test": self.option("q_test"),
            "tree_number": self.option("tree_number"),
        }
        self.two_ran_lefse.set_options(options)
        if self.option('intersection'):
            self.two_ran_lefse.on("end", self.run_roc_intersection)
        elif self.option('lefse_analysis'):
            self.two_ran_lefse.on("end", self.run_roc_lefse)
        elif self.option('two_group_analysis'):
            self.two_ran_lefse.on("end", self.run_roc_two)
        elif self.option('ran_for_analysis'):
            self.two_ran_lefse.on("end", self.run_roc_ran)
        elif self.option('env_table').is_set:
            self.two_ran_lefse.on("end", self.run_roc_env)
        else:
            self.set_error('参数不正确，无法正常进行roc分析', code="12704001")
        self.two_ran_lefse.run()

    def run_roc_intersection(self):
        fin_otu = self.get_intersection_otu()
        if self.option('roc_calc_method') == "MI":
            mi_group = os.path.join(self.work_dir, "intersection_mi_otu.xls")
            self.get_mi_otu(fin_otu, self.option('two_ran_group').prop['path'], mi_group, self.dict)
        else:
            mi_group = self.option('two_ran_group')
        # new_table = self.change_otuname(fin_otu)
        options = {
            'otu_table': fin_otu,
            'level': self.option('level'),
            'group_table': mi_group,
            'method': self.option('roc_calc_method'),
            'top_n': 0
        }
        self.roc_intersection.set_options(options)
        if self.option('env_table').is_set:
            self.roc_intersection.on("end", self.run_roc_env)
        else:
            self.roc_intersection.on("end", self.set_db)
        self.roc_intersection.run()

    def run_roc_env(self):
        env_otu = self.get_env_otu(self.option('env_table').prop['path'])
        if self.option('roc_calc_method') == "MI":
            mi_group = os.path.join(self.work_dir, "roc_mi_otu.xls")
            self.get_mi_otu(env_otu, self.option('two_ran_group').prop['path'], mi_group, self.dict)
        else:
            mi_group = self.option('two_ran_group')
        options = {
            'otu_table': env_otu,  # 把传入的env表格进行修改，otu表一个格式的进行
            'level': self.option('level'),
            'group_table': mi_group,  # 分组文件表
            'method': self.option('roc_calc_method'),
            'top_n': 0
        }
        self.roc_env.set_options(options)
        self.roc_env.on("end", self.set_db)
        self.roc_env.run()

    def run_roc_lefse(self):
        lefse_otu = self.get_otu_table(name="lefse")
        if self.option('roc_calc_method') == "MI":
            mi_group = os.path.join(self.work_dir, "lefse_mi_otu.xls")
            self.get_mi_otu(lefse_otu, self.option('two_ran_group').prop['path'], mi_group, self.dict)
        else:
            mi_group = self.option('two_ran_group')
        # new_table = self.change_otuname(lefse_otu)
        options = {
            'otu_table': lefse_otu,
            'level': self.option('level'),
            'group_table': mi_group,
            'method': self.option('roc_calc_method'),
            'top_n': 0
        }
        self.roc_lefse.set_options(options)
        if self.option('two_group_analysis'):
            self.roc_lefse.on("end", self.run_roc_two)
        elif self.option('ran_for_analysis'):
            self.roc_lefse.on("end", self.run_roc_ran)
        elif self.option('env_table').is_set:
            self.roc_lefse.on("end", self.run_roc_env)
        else:
            self.roc_lefse.on("end", self.set_db)
        self.roc_lefse.run()

    def run_roc_two(self):
        two_otu = self.get_otu_table(name="two_group")
        if self.option('roc_calc_method') == "MI":
            mi_group = os.path.join(self.work_dir, "two_mi_otu.xls")
            self.get_mi_otu(two_otu, self.option('two_ran_group').prop['path'], mi_group, self.dict)
        else:
            mi_group = self.option('two_ran_group')
        # new_table = self.change_otuname(two_otu)
        options = {
            'otu_table': two_otu,
            'level': self.option('level'),
            'group_table': mi_group,
            'method': self.option('roc_calc_method'),
            'top_n': 0
        }
        self.roc_two.set_options(options)
        if self.option('ran_for_analysis'):
            self.roc_two.on("end", self.run_roc_ran)
        elif self.option('env_table').is_set:
            self.roc_two.on("end", self.run_roc_env)
        else:
            self.roc_two.on("end", self.set_db)
        self.roc_two.run()

    def run_roc_ran(self):
        ran_otu = self.get_otu_table(name="random_forest")
        if self.option('roc_calc_method') == "MI":
            mi_group = os.path.join(self.work_dir, "ran_mi_otu.xls")
            self.get_mi_otu(ran_otu, self.option('two_ran_group').prop['path'], mi_group, self.dict)
        else:
            mi_group = self.option('two_ran_group')
        # new_table = self.change_otuname(ran_otu)
        options = {
            'otu_table': ran_otu,
            'level': self.option('level'),
            'group_table': mi_group,
            'method': self.option('roc_calc_method'),
            'top_n': 0
        }
        self.roc_ran.set_options(options)
        if self.option('env_table').is_set:
            self.roc_ran.on("end", self.run_roc_env)
        else:
            self.roc_ran.on("end", self.set_db)
        self.roc_ran.run()

    def get_env_otu(self, env_table):
        con = pd.read_table(env_table, header=0, sep="\t")
        cont = con.T
        sample_list = [i for i in cont.ix[0]]
        cont = cont.drop("#SampleID")
        cont.columns = sample_list
        cont.index.name = "OTU ID"
        table = os.path.join(self.work_dir, "env_table")
        cont.to_csv(table, sep="\t", index=True)
        return table

    def get_otu_table(self, name=None):
        if name == "lefse":
            file_path = os.path.join(self.two_ran_lefse.output_dir, "lefse_LDA.xls")
            con = pd.read_table(file_path, header=0, sep="\t")
            con = con[con["pvalue"] != "-"]
            if len(con) == 0:
                self.set_error("lefse表格筛选不出OTU表", code="12704002")
            if self.option('lefse_cho') == "p-value":
                con = con.sort_values(by="pvalue", ascending=False)
            else:
                con = con.sort_values(by="lda", ascending=False)
            con = con.iloc[:self.option('lefse_num')]
            lefse_spe_list_old = [i for i in con['taxon']]  # lda表中的物种
            lefse_spe_list = []
            for i in lefse_spe_list_old:
                c = i.split('.')
                lefse_spe_list.append((';').join(c))
            ori_otu = pd.read_table(os.path.join(self.two_ran_lefse.output_dir, "lefse_otu_table"), header=0, sep="\t")
            spe_list = [i for i in ori_otu['OTU ID']]  # 原始otu表中的物种
            ori_otu.index = spe_list
            index = set(lefse_spe_list) & set(spe_list)
            index_list = [i for i in index]
            if len(index_list) == 0:
                self.set_error("lefse表格筛选不出OTU表", code="12704003")
            # 重置index的物种应该为交集
            new_otu = ori_otu.reindex(index=index_list)  # 重置筛选物种
            table = os.path.join(self.work_dir, 'lefse_roc_input_otu.xls')
            new_otu.to_csv(table, sep="\t", index=False)
        if name == "two_group":
            file_path = os.path.join(self.two_ran_lefse.output_dir, "two_group_table.xls")
            con = pd.read_table(file_path, header=0, sep="\t")
            if self.option('two_group_cho') == "p-value":
                con = con.sort_values(by="pvalue", ascending=False)
            else:
                con = con.sort_values(by="corrected_pvalue", ascending=False)
            con = con.iloc[:self.option('two_group_num')]
            two_group_spe_list = [i for i in con[' ']]
            ori_otu = pd.read_table(os.path.join(self.work_dir, 'change_otu_name'), header=0, sep="\t")
            spe_list = [i for i in ori_otu['OTU ID']]
            ori_otu.index = spe_list
            new_otu = ori_otu.reindex(index=two_group_spe_list)
            table = os.path.join(self.work_dir, 'two_group_otu.xls')
            new_otu.to_csv(table, sep="\t", index=False)
        if name == "random_forest":
            file_path = os.path.join(self.two_ran_lefse.output_dir, "Random_table.xls")
            con = pd.read_table(file_path, header=0, sep="\t")
            con = con.sort_values(by="MeanDecreaseAccuracy", ascending=False)
            con = con.iloc[:self.option('Ran_for_num')]
            index = con.index
            ori_otu = pd.read_table(os.path.join(self.work_dir, 'change_otu_name'), header=0, sep="\t")
            spe_list = [i.split(";")[-1] for i in ori_otu['OTU ID']]
            ori_otu.index = spe_list
            new_otu = ori_otu.reindex(index=index)
            table = os.path.join(self.work_dir, 'random_forest_otu.xls')
            new_otu.to_csv(table, sep="\t", index=False)
        return table

    def get_intersection_otu(self):
        path =[]
        if self.option('lefse_analysis'):
            table_lefse = self.get_otu_table(name="lefse")
            path.append(table_lefse)
        if self.option('two_group_analysis'):
            table_two = self.get_otu_table(name="two_group")
            path.append(table_two)
        if self.option('ran_for_analysis'):
            table_ran = self.get_otu_table(name="random_forest")
            path.append(table_ran)
        if len(path) == 1:
            table = path[0]
        else:
            T = 0
            new = set([])
            for i in path:
                T += 1
                con = pd.read_table(i, header=0, sep="\t")
                spe_list = [m.split(";")[-1] for m in con['OTU ID']]  # 取最后一个分类水平的名称
                old = new
                new = set(spe_list)
                self.logger.info("{}的spe集合{}".format(i, new))
                if T >= 2:
                    new = old & new
            all_spe = [i for i in new]
            self.logger.info("all_spe:{}".format(all_spe))
            if len(all_spe) == 0:
                self.set_error("共有物种数目为0，无法进行roc分析", code="12704004")
            if self.option('two_group_analysis'):
                any_one_con = pd.read_table(table_two, header=0, sep="\t")
            else:
                any_one_con = pd.read_table(table_ran, header=0, sep="\t")
            spe_list = [i.split(";")[-1] for i in any_one_con['OTU ID']]
            any_one_con.index = spe_list
            new_otu = any_one_con.reindex(index=all_spe)
            table = os.path.join(self.work_dir, 'intersection_otu.xls')
            new_otu.to_csv(table, sep="\t", index=False)
        return table

    def get_mi_otu(self, otu_table, group_table, result_path, dict_):
        group = pd.read_table(group_table, header=0, sep="\t")
        group_dict = {}
        for i in group.index:
            group_dict[group["#sample"][i]] = group[group.columns[1]][i]
        table = pd.read_table(otu_table, header=0, sep="\t")
        samp_list = []
        for i in table.columns:
            samp_list.append(i)
        g_list = []
        for i in samp_list[1:]:
            g_list.append(group_dict[i])
        tablet = table.T
        tablet = tablet.drop("OTU ID")
        tablet["g_name"] = g_list
        sum_table = tablet.groupby('g_name').sum().T
        spe_group = []
        for i in sum_table.index:
            if sum_table.ix[i][0] < sum_table.ix[i][1]:
                spe_group.append(sum_table.columns[0])
            else:
                spe_group.append(sum_table.columns[1])
        table["sp_group"] = spe_group
        fin_otu = table.T
        fin_otu = fin_otu.drop("OTU ID")
        a_number = 0
        b_number = 0
        for i in spe_group:
            if i == sum_table.columns[0]:
                a_number = a_number + 1
            else:
                b_number = b_number + 1
        if a_number == 0 or b_number == 0:
            self.logger.error("筛选后物种无法进行MI指数的roc分析{}".format(result_path))
            self.set_error("MI指数的roc分析失败", code="12704005")
        data = {}
        con = pd.DataFrame(data)
        for i in fin_otu.index:
            a = 0
            b = 0
            if i == "sp_group":
                break
            for m in range(len(spe_group)):
                if fin_otu.ix['sp_group'][m] == sum_table.columns[0]:
                    a = a + int(fin_otu.ix[i][m])
                else:
                    b = b + int(fin_otu.ix[i][m])
            a = a / a_number
            b = b / b_number
            if data:
                con[i] = [a, b]
            else:
                data = {i: [a, b]}
                con = pd.DataFrame(data)
        con.index = [sum_table.columns[0], sum_table.columns[1]]
        fin_con = con.T
        fin_value = []
        for i in fin_con.index:
            value = fin_con.ix[i][0] * float(dict_[fin_con.columns[0]]) - fin_con.ix[i][1] * float(dict_[fin_con.columns[1]])
            fin_value.append(value)
        fin_con["value"] = fin_value
        fin_con["group_name"] = g_list
        fin_con = fin_con.drop([sum_table.columns[0], sum_table.columns[1]], axis=1)
        fin_con.index.name = "#Sample"
        fin_con.to_csv(result_path, sep="\t", index=True)

    def end(self):
        """
        上传结果目录
        """
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "LEfSe差异分析结果目录"],
        #     ["./lefse_LDA.cladogram.png", "png", "LEfSe分析cladogram结果图片"],
        #     ["./lefse_LDA.png", "png", "LEfSe分析LDA图片"],
        #     ["./lefse_LDA.xls", "xls", "LEfSe分析lda数据表"]
        # ])
        super(RocNewWorkflow, self).end()

    def change_otuname(self, tablepath):
        newtable = os.path.join(self.work_dir, 'change_otu_name')
        f2 = open(newtable, 'w+')
        with open(tablepath, 'r') as f:
            i = 0
            for line in f:
                if i == 0:
                    i = 1
                    f2.write(line)
                else:
                    line = line.strip().split('\t')
                    line_data = line[0].strip().split(' ')
                    line_he = "".join(line_data)
                    line[0] = line_he
                    # line[0] = line_data[-1]
                    for i in range(0, len(line)):
                        if i == len(line) - 1:
                            f2.write("%s\n" % (line[i]))
                        else:
                            f2.write("%s\t" % (line[i]))
        f2.close()
        return newtable

    def set_db(self):
        """
        设置结果文件，将数据导入mongo库中
        """
        if self.option("intersection"):
            self.linkdir(self.roc_intersection.output_dir, "Intersection_roc")
        else:
            if self.option("lefse_analysis"):
                self.linkdir(self.roc_lefse.output_dir, "Lefse_roc")
            if self.option("two_group_analysis"):
                self.linkdir(self.roc_two.output_dir, "Two_roc")
            if self.option("ran_for_analysis"):
                self.linkdir(self.roc_ran.output_dir, "Ran_roc")
        if self.option('env_table').is_set:
            self.linkdir(self.roc_env.output_dir, "Env_roc")
        api_roc = self.api.roc
        dir_name = os.listdir(self.output_dir)
        for name in dir_name:
            dir_path = os.path.join(self.output_dir, name)
            file_name = os.listdir(dir_path)
            for f_name in file_name:
                if f_name == "roc_auc.xls":
                    api_roc.add_roc_auc(file_path=os.path.join(dir_path, f_name), table_id=self.option("main_id"),
                                        dir_name=name)
                if f_name == "roc_curve.xls":
                    api_roc.add_roc_curve(file_path=os.path.join(dir_path, f_name), table_id=self.option("main_id"),
                                          dir_name=name)
                if f_name == "roc_plot_rocarea.xls":
                    api_roc.add_roc_area(file_path=os.path.join(dir_path, f_name), table_id=self.option("main_id"),
                                         dir_name=name)
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))