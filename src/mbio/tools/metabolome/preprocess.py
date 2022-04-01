# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import pandas as pd
from mbio.packages.metabolome.preprocess import Preprocess
from mbio.packages.metabolome.common import combine_table
from mbio.packages.metabolome.common import check_head
from mbio.packages.metabolome.common import check_index
from mbio.packages.metabolome.contend_name import Contend_Name
import unittest

class PreprocessAgent(Agent):
    """
    代谢组数据预处理
    version 1.0
    author: guhaidong
    last_modify: 2018.6.13
    """

    def __init__(self, parent):
        super(PreprocessAgent, self).__init__(parent)
        options = [
            {'name': 'ana_method', 'type': 'string', 'default': 'GC'},  # 是LC、GC分析
            {'name': 'interactive', 'type': 'string', 'default': 'False'},  # 是否是交互分析
            {'name': 'pos_table', 'type': 'infile', 'format': 'metabolome.metab_table'},  # 工作流的分析输入
            {'name': 'neg_table', 'type': 'infile', 'format': 'metabolome.metab_table'},  # LC的工作流分析输入
            {'name': 'pos_rout', 'type': 'infile', 'format': 'metabolome.metab_table_dir'},  # 交互分析的输入
            {'name': 'neg_rout', 'type': 'infile', 'format': 'metabolome.metab_table_dir'},  # LC的交互分析的输入
            {'name': 'group_table', 'type': 'infile', 'format': 'meta.otu.group_table'},  # 分组文件
            {'name': 'fillna', 'type': 'string', 'default': 'min'},  # 缺失值填充方法: min/median/mean/rf/none
            {'name': 'rsd', 'type': 'string', 'default': '30%'},
            # QC验证RSD，需数字+% 数字范围 25-30 ,允许填"none", 用于QC样本量足够但是不进行rsd质控 add by shaohua.yuan 20180920
            {'name': 'norm', 'type': 'string', 'default': 'sum'},  # 数据归一化方法：median/mean/sum/sample/inner/none
            {'name': 'sample_name', 'type': 'string'},  # 标准化样品， norm = sample时使用
            {'name': 'inner_ref', 'type': 'string'},  # 内参代谢物，norm = inner时使用
            {'name': 'scale', 'type': 'string', 'default': 'none'},  # 取log值，形式log+数字/"none"/"defined"
            {'name': 'log', 'type': 'int'},  # 取log值，当scale = defined时使用
            {'name': 'org_pos_out', 'type': 'outfile', 'format': 'metabolome.metab_table_dir'},  # 处理的阳离子原始表
            {'name': 'org_neg_out', 'type': 'outfile', 'format': 'metabolome.metab_table_dir'},  # 处理的阴离子原始表
            {'name': 'org_set', 'type': 'outfile', 'format': 'metabolome.metabset'},  # 初始代谢集
            {'name': 'has_name_org_set', 'type': 'outfile', 'format': 'metabolome.metabset'}, #有代谢物名称的代谢集
            {'name': 'org_mul_set', 'type': 'outfile', 'format': 'metabolome.mul_metabset'},  # 可能要增加的代谢集文件
            {'name': 'has_name_org_mul_set', 'type': 'outfile', 'format': 'metabolome.mul_metabset'},  #有代谢物名称的原始代谢集
            {'name': 'has_name_origin_mul_set', 'type': 'outfile', 'format': 'metabolome.mul_metabset'},  #有代谢物名称的origin 代谢集
            {'name': 'pos_out', 'type': 'outfile', 'format': 'metabolome.metab_table_dir'},  # 阳离子结果表
            {'name': 'neg_out', 'type': 'outfile', 'format': 'metabolome.metab_table_dir'},  # 阴离子结果表
            {'name': 'mix_out', 'type': 'outfile', 'format': 'metabolome.metab_table_dir'},  # 合并后的结果表
            {'name': 'rm_nan','type':'float','default':50}, #20190603 v2.0  mean: 50%
            {"name": "fill_type", "type":"string","default":"all"},  #v3 202003
            {"name": "raw_cv", "type": "string", "default": "T" },  # T: 计算raw数据的cv 值， 否则不计算 #v3
            {"name": "task_version", "type":"string","default":"3.0"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if self.option('ana_method') not in ["GC", "LC"]:
            raise OptionError("分析流程类型必须为GC或LC，debug: %s" , variables=(self.option('ana_method')), code="34701801")
        if self.option('interactive') not in ["True", "False"]:
            raise OptionError("交互分析判断参数必须为True或False, debug: %s" , variables=(self.option('interactive')), code="34701802")
        if self.option('ana_method') == 'GC' and self.option('interactive') == 'False':
            if not self.option('pos_table').is_set:
                raise OptionError("请传入阳离子原始表", code="34701803")
        elif self.option('ana_method') == 'GC' and self.option('interactive') == 'True':
            if not self.option('pos_rout').is_set:
                raise OptionError("请传入阳离子路径", code="34701804")
        elif self.option('ana_method') == 'LC' and self.option('interactive') == 'False':
            if not self.option('pos_table').is_set:
                raise OptionError("请传入阳离子原始表", code="34701805")
            if not self.option('neg_table').is_set:
                raise OptionError("请传入阴离子原始表", code="34701806")
        elif self.option('ana_method') == 'LC' and self.option('interactive') == 'True':
            if not self.option('pos_rout').is_set:
                raise OptionError("请传入阳离子路径", code="34701807")
            if not self.option('neg_rout').is_set:
                raise OptionError("请传入阴离子路径", code="34701808")
        if not self.option('group_table').is_set:
            raise OptionError("请传入分组文件", code="34701809")
        if self.option('fillna') not in ['min', 'median', 'mean', 'rf', 'none']:
            raise OptionError("fillna 参数不在范围内", code="34701810")
        self.logger.info(self.option('rsd'))
        if self.option('rsd') != "none" and not self.option('rsd').endswith("%"):
            raise OptionError("rsd参数必须以%结尾", code="34701811")
        else:
            if self.option('rsd') != "none":
                try:
                    check_int = self.option('rsd').rstrip("%")
                    check_int = int(check_int)
                except:
                    raise OptionError("rsd参数必须是整数+%", code="34701812")
                if check_int < 20 or check_int > 30:
                    raise OptionError("rsd参数必须在20%-30%范围内", code="34701813")
        if self.option('norm') not in ["median", 'mean', 'sum', 'sample', 'inner', 'none']:
            raise OptionError("norm参数不在范围内", code="34701814")
        elif self.option('norm') == 'sample':
            if not self.option('sample_name'):
                raise OptionError("选择内参样品时，需指定样品", code="34701815")
        elif self.option('norm') == 'inner':
            if not self.option('inner_ref'):
                raise OptionError("须填写内参代谢物", code="34701816")
        if self.option('scale') not in ["none", "defined"]:
            if not self.option('scale').startswith("log"):
                raise OptionError("scale 参数必须为log+数字，或者为none，defined", code="34701817")
            else:
                check_int = self.option('scale').lstrip("log")
                try:
                    check_int = int(check_int)
                except:
                    raise OptionError("scale 参数必须为log + 整数，或者为none, defined", code="34701818")
                if check_int < 2:
                    raise OptionError("scale 参数的整数部分必须大于等于2", code="34701819")
        elif self.option('scale') == 'defined':
            if self.option('log') < 2:
                raise OptionError("log参数值必须大于等于2的整数", code="34701820")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '15G'

    def end(self):
        super(PreprocessAgent, self).end()


class PreprocessTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(PreprocessTool, self).__init__(config)
        self.pos_table = ""
        self.neg_table = ""
        self.pos_desc = ""
        self.neg_desc = ""
        self.metab_list = []
        self.has_name_metab_list = []
        #self.piplist = ["fillna", "control", "norm", "log", "to_table"]
        #self.piplist = ["rm_nan","fillna", "control", "norm", "log", "to_table"]  #add rm_nan 20190604
        self.piplist = ["rm_nan","fillna", "norm","control", "log", "to_table"]
        self.fillna_method = self.option('fillna')  # 需要具体填写
        if self.option("rsd") != "none":              # add by shaohua.yuan 20180920
            self.control_mehod = float(self.option('rsd').rstrip("%")) / 100
        else:
            self.control_mehod = "none"
        self.logger.info(self.control_mehod)
        self.norm_method = self.option('norm')
        if self.norm_method == "sample":
            self.norm_method += ":" + self.option("sample_name")
        elif self.norm_method == "inner":
            self.norm_method += ":" + self.option("inner_ref")
        self.log_method = self.option("scale")
        if self.log_method == "defined":
            self.log_method = self.option("log")
        elif "log" in self.log_method:
            self.log_method = int(self.log_method.lstrip("log"))

        self.get_samples_len()

    def get_samples_len(self):
        with open(self.option('group_table').path, 'r') as fr:
            lines = fr.readlines()
        self.samples_num = len(lines) - 1



    def run_gc_data(self):
        # self.logger.info(">>>>>>in func run_gc_data")
        pos_table = self.check_dir(os.path.join(self.output_dir, "org_pos/metab_abund.txt"))
        pos_desc = self.check_dir(os.path.join(self.output_dir, "org_pos/metab_desc.txt"))
        # data = pd.read_table(self.option('pos_table').path, index_col="Metabolite").fillna("-")  # 表头必须固定，用-防止表格缺失内容
        data = pd.read_table(self.option("pos_table").path).fillna("-")
        data['Metabolite'] = data['Metabolite'].apply(lambda x: re.sub(r"\bpos", "pos", x, flags=re.I))
        data['Mode'] = data["Mode"].apply(lambda x: re.sub(r"\bpos", "pos", x, flags=re.I))
        data.set_index("Metabolite", inplace=True)
        #data = data[data[data.columns[8:]].apply(lambda x: x.max() > 0, axis=1)]  # 去丰度全为0的数据  ##列固定，要根据样本数。改成不为0的比例!!!
        data = data[~data.index.duplicated(keep='first')]  # 去名称重复的数据
        # filter_index = data[data.index.str.startswith('pos') == False].index  # 筛选有名称的数据
        #pos_set = set(data.index)
        pos_set = self.mix_add(data.index)  #生成无重复的代谢物list
        # pos_set = set(filter_index)
        metab_map = {}
        for index, one in enumerate(pos_set):
            metab_map[one] = 'metab_%s' % index
            #if one.startswith('pos'):
            self.metab_list.append('metab_%s' % index)
            if re.match(r"^pos\d+", one, flags=re.I):  # modify 20181115 for metab eg. postion
                continue
            self.has_name_metab_list.append('metab_%s' % index)
        self.write_table(data, metab_map, pos_table, pos_desc)
        # self.logger.info(">>>>>>func run_gc_data end")
        return pos_table, pos_desc

    def run_lc_data(self):
        # self.logger.info(">>>>>>in func run_lc_data")
        pos_table = self.check_dir(os.path.join(self.output_dir, "org_pos/metab_abund.txt"))
        neg_table = self.check_dir(os.path.join(self.output_dir, "org_neg/metab_abund.txt"))
        pos_desc = self.check_dir(os.path.join(self.output_dir, "org_pos/metab_desc.txt"))
        neg_desc = self.check_dir(os.path.join(self.output_dir, "org_neg/metab_desc.txt"))
        self.org_pos_abund_table = pos_table
        self.org_neg_abund_table = neg_table
        self.org_pos_desc_table = pos_desc
        self.org_neg_desc_table = neg_desc
        # pos_data = pd.read_table(self.option('pos_table').path, index_col="Metabolite").fillna("-")
        # if self.option('norm') == "inner":
        #     self.logger.info('Contend_Name %s'%self.option("inner_ref"))
        #     CN = Contend_Name(self.option("pos_table").path, self.option("neg_table").path, not_dispose=self.option("inner_ref")) #zouguanqing 20190604
        # else:
        self.logger.info('Contend_Name')
        CN = Contend_Name(self.option("pos_table").path, self.option("neg_table").path)
        CN.contend_pip(sample_info=self.option("group_table").path) #produce pos.rm_rep.csv and neg.rm_rep.csv
        pos_table_change_name = self.work_dir + '/pos.rm_rep.csv'
        neg_table_change_name = self.work_dir + '/neg.rm_rep.csv'
        #pos_data = pd.read_table(self.option("pos_table").path).fillna("-")
        pos_data = pd.read_table(pos_table_change_name).fillna("-")
        self.id_num = 0
        def change_2_id(x):  #zouguanqing 20190726
            if x['Metabolite']=='-':
                x['Metabolite']= 'id_'+ str(self.id_num)
                self.id_num +=1
            return x
        pos_data = pos_data.apply(change_2_id,axis=1)
        pos_data['Metabolite'] = pos_data['Metabolite'].apply(lambda x: re.sub(r"\bpos", "pos", x, flags=re.I))
        pos_data['Mode'] = pos_data["Mode"].apply(lambda x: re.sub(r"\bpos", "pos", x, flags=re.I))
        pos_data.set_index("Metabolite", inplace=True)
        # neg_data = pd.read_table(self.option("neg_table").path, index_col="Metabolite").fillna("-")
        #neg_data = pd.read_table(self.option("neg_table").path).fillna("-")
        neg_data = pd.read_table(neg_table_change_name).fillna("-")
        neg_data = neg_data.apply(change_2_id,axis=1)
        neg_data['Metabolite'] = neg_data['Metabolite'].apply(lambda x: re.sub(r"\bneg", "neg", x, flags=re.I))

        neg_data['Mode'] = neg_data["Mode"].apply(lambda x: re.sub(r"\bneg", "neg", x, flags=re.I))
        neg_data.set_index("Metabolite", inplace=True)
        ##pos_data = pos_data[pos_data[pos_data.columns[8:]].apply(lambda x: x.max() > 0, axis=1)]
        ##neg_data = neg_data[neg_data[neg_data.columns[8:]].apply(lambda x: x.max() > 0, axis=1)]

        ###pos_data = pos_data[~pos_data.index.duplicated(keep='first')]
        ###neg_data = neg_data[~neg_data.index.duplicated(keep='first')]


        # pos_filter_index = pos_data[pos_data.index.str.startswith('pos') == False].index  # 筛选有名称的数据
        # neg_filter_index = neg_data[neg_data.index.str.startswith('neg') == False].index  # 筛选有名称的数据
        ############# 不区分大小写
        # mix_set = set(pos_data.index) | set(neg_data.index)
        mix_set = self.mix_add(pos_data.index.tolist(), neg_data.index.tolist())
        # filter_mix_set = set(pos_filter_index) | set(neg_filter_index)
        metab_map = {}
        for index, one in enumerate(mix_set):
            metab_map[one] = 'metab_%s' % index
            #if one.startswith("pos") or one.startswith("neg"):
            self.metab_list.append('metab_%s' % index)
            if re.match(r"^pos\d+|^neg\d+|^id_\d+$", one, flags=re.I):  # modify 20181115 for metab eg. postion
                continue
            self.has_name_metab_list.append('metab_%s' % index)

        self.write_table(pos_data, metab_map, pos_table, pos_desc)
        self.write_table(neg_data, metab_map, neg_table, neg_desc)
        # self.logger.info(">>>>>>func run_gc_data end")
        return pos_table, neg_table, pos_desc, neg_desc

    def rm_special(self, metab_name):
        '''
        目前只去除<i>,</i>
        '''
        metab_name = metab_name.replace("<i>","").replace("</i>","")
        return metab_name

    def mix_add(self, pos_list, neg_list=None):
        new_mix_set = []
        for each in pos_list:
            each = each.lower()
            if not each in new_mix_set:
                new_mix_set.append(each)
        if neg_list:
            for each in neg_list:
                each = each.lower()
                if not each in new_mix_set:
                    new_mix_set.append(each)
        return new_mix_set

    def write_table(self, df, metab_map, abund, desc):
        # self.logger.info(">>>>>>in func write_table")
        f1 = open(abund, 'w')
        f2 = open(desc, 'w')
        # f1_head = df.columns[8:]  # 不取Metabolite为index时使用
        col_nums = df.shape[1]
        desc_num = col_nums - self.samples_num  ## 因为有不固定的列，根据样本数来确定不固定列的最大列数
        f1_head = df.columns[desc_num:]
        abund_data = df[f1_head]
        f1.write("metab_id\t" + "\t".join(f1_head.tolist()) + '\n')
        # f2_head = df.columns[:8]  # 不取Metabolite为index时使用
        f2_head = df.columns[:desc_num]
        desc_data = df[f2_head].replace(['_', '--', '—'], '-')

        f2.write("metab_id\tMetabolite\t" + "\t".join(f2_head.tolist()) + '\n')  # 取Metabolite为index时手动写入Metabolite
        for i in df.index:
            i_lower = i.lower()
            f1.write("%s\t%s\n" % (metab_map[i_lower], "\t".join(abund_data.loc[i].astype("string").tolist())))
            if re.match('^id_\d+$',i_lower) :
                meta = metab_map[i_lower]  ##'-',没有 Metabolite用 meta id代替
            else:
                meta = self.rm_special(i)
            f2.write("%s\t%s\t%s\n" % (metab_map[i_lower], meta ,"\t".join(desc_data.loc[i].astype("string").tolist())))
        f1.close()
        f2.close()
        # self.logger.info(">>>>>>func write_table end")

    def run_set(self):
        # self.logger.info(">>>>>>in func run_set")
        set_file = open(self.output_dir + '/metabset.list', 'w')
        mul_set_file = open(self.output_dir + '/mul_metabset.list', 'w')
        mul_set_file.write("origin\t%s" % ','.join(self.metab_list))
        set_file.write("\n".join(self.metab_list))
        set_file.close()
        mul_set_file.close()

        has_name_set_file = open(self.output_dir + '/has_name_metabset.list', 'w')
        has_name_mul_set_file = open(self.output_dir + '/has_name_mul_metabset.list', 'w')
        has_name_mul_set_file.write("origin\t%s" % ','.join(self.has_name_metab_list))
        has_name_set_file.write('\n'.join(self.has_name_metab_list))
        has_name_set_file.close()
        has_name_mul_set_file.close()


    def run_set_origin(self):
        if self.option('ana_method') == 'LC':
            ori_desc_mix = self.output_dir + '/mix/metab_desc.txt'
        else:
            ori_desc_mix = self.output_dir + '/pos/metab_desc.txt'
        ori_meta_list = []
        data = pd.read_table(ori_desc_mix,sep='\t',header=0)
        for i in range(len(data)):
            if not re.match('^metab_\d+$',data['Metabolite'][i]):
                ori_meta_list.append(data['metab_id'][i])
        if len(ori_meta_list) >0:
            with open(self.output_dir + '/has_name_origin_mul_metabset.list', 'w') as has_name_origin_set_file:
                has_name_origin_set_file.write("origin\t%s"% ','.join(ori_meta_list))

        # self.logger.info(">>>>>>func run_set end")

    def get_control_group(self):
        # self.logger.info(">>>>>>in func get_control_group")
        group_data = pd.read_table(self.option('group_table').path)
        group_data.index = group_data[group_data.columns[0]]
        group_data = group_data[group_data[group_data.columns[1]] == "QC"]
        if self.option("rsd") == "none":
            self.piplist.remove("control")
        elif len(group_data) < 2:
            self.piplist.remove("control")
        # self.logger.info(">>>>>>func get_control_group end")
        self.qc_samples = group_data.index.tolist()
        return group_data.index.tolist()

    def run_gc_preprocess(self):
        # self.logger.info(">>>>>>in func run_gc_preprocess")
        group_list = self.get_control_group()
        if self.option("norm") == 'inner':
            prep_obj = Preprocess(self.pos_table, "metab_id", group_list, inner_file=self.pos_inner_file,group_file=self.option('group_table').path)
        else:
            prep_obj = Preprocess(self.pos_table, "metab_id", group_list,group_file=self.option('group_table').path)
        tmp_pos_abund = self.check_dir(os.path.join(self.work_dir, "pos/metab_abund.txt"))
        pos_abund = self.check_dir(os.path.join(self.output_dir, "pos/metab_abund.txt"))
        pipdic = {
            "fillna": self.fillna_method,
            "control": self.control_mehod,
            "norm": self.norm_method,
            "log": self.log_method,
            "to_table": tmp_pos_abund,
            "rm_nan" : float(self.option('rm_nan'))/100,  #zgq 20190604
            "rm_type" : self.option("fill_type") #2020 v3
        }
        # self.logger.info(">>>>>>>>>>>piplist: %s\npipdic: %s" % (self.piplist, pipdic))
        metab_list = prep_obj.pipline(self.piplist, pipdic)
        index_table = os.path.join(self.work_dir, "pos_remain.index")
        with open(index_table, 'w') as f:
            self.logger.info(">>>DEBUG")
            self.logger.info(type(metab_list))
            self.logger.info(metab_list)
            f.write("\n".join(metab_list.tolist()))
        combine_table(self.pos_desc, os.path.join(self.output_dir, "pos/metab_desc.txt"), index_table)
        data = pd.read_table(tmp_pos_abund)
        if "control" in self.piplist:
            data = data.drop(["rsd", "sum"], axis=1)
        data.to_csv(pos_abund, sep="\t", index=False)
        # self.logger.info(">>>>>>func run_gc_preprocess end")

    def run_lc_preprocess(self,pip_type=None):
        # self.logger.info(">>>>>>in func run_lc_preprocess")
        group_list = self.get_control_group()
        if self.option('norm') == 'inner':
            pre_pos_obj = Preprocess(self.pos_table, "metab_id", group_list, inner_file=self.pos_inner_file,group_file=self.option('group_table').path)
            pre_neg_obj = Preprocess(self.neg_table, "metab_id", group_list, inner_file=self.neg_inner_file,group_file=self.option('group_table').path)
        else:
            pre_pos_obj = Preprocess(self.pos_table, "metab_id", group_list,group_file=self.option('group_table').path)
            pre_neg_obj = Preprocess(self.neg_table, "metab_id", group_list,group_file=self.option('group_table').path)
        tmp_pos_abund = self.check_dir(os.path.join(self.work_dir, "pos/metab_abund.txt"))
        tmp_neg_abund = self.check_dir(os.path.join(self.work_dir, "neg/metab_abund.txt"))
        pos_abund = self.check_dir(os.path.join(self.output_dir, "pos/metab_abund.txt"))
        neg_abund = self.check_dir(os.path.join(self.output_dir, "neg/metab_abund.txt"))
        new_pos_desc = self.check_dir(os.path.join(self.output_dir, "pos/metab_desc.txt"))
        new_neg_desc = self.check_dir(os.path.join(self.output_dir, "neg/metab_desc.txt"))
        pipdic = {
            "fillna": self.fillna_method,
            "control": self.control_mehod,
            "norm": self.norm_method,
            "log": self.log_method,
            "rm_nan" : float(self.option('rm_nan'))/100,  #zgq 20190604
            "rm_type" : self.option("fill_type") #2020 v3
        }
        pipdic.update({"to_table": tmp_pos_abund})
        self.logger.info('self.piplist: %s'%(',').join(self.piplist))
        self.logger.info(str(pipdic))
        metab_pos_list = pre_pos_obj.pipline(self.piplist, pipdic)
        pipdic.update({"to_table": tmp_neg_abund})
        metab_neg_list = pre_neg_obj.pipline(self.piplist, pipdic)
        index_pos = os.path.join(self.work_dir, "pos_remain.index")
        index_neg = os.path.join(self.work_dir, "neg_remain.index")
        with open(index_pos, 'w') as f1:
            f1.write("\n".join(metab_pos_list.tolist()))
        with open(index_neg, 'w') as f2:
            f2.write("\n".join(metab_neg_list.tolist()))
        combine_table(self.pos_desc, new_pos_desc, index_pos)
        combine_table(self.neg_desc, new_neg_desc, index_neg)

        repeat_id = set(metab_pos_list.tolist()) & set(metab_neg_list.tolist())
        if self.option('interactive')=='True' and len(repeat_id) != 0:  ##20200303 兼容v1版数据预处理的交互分析
            self.logger.info("V1 interactive and has repeat_id ............ ")
            self.logger.info("repeat : %s"%(','.join(repeat_id)))
            mix_set = self.merge_mix(tmp_pos_abund, tmp_neg_abund, pos_abund, neg_abund,new_pos_desc, new_neg_desc)
            self.merge_mix_desc(new_pos_desc, new_neg_desc, mix_set)

        else:
            self.treat_data_v2(tmp_pos_abund, pos_abund)
            self.treat_data_v2(tmp_neg_abund, neg_abund)
            self.merge_mix_v2(new_pos_desc, new_neg_desc,pos_abund, neg_abund,type='qc')
            if pip_type:  #工作流 增加原始的mix表
                self.merge_mix_v2(self.org_pos_desc_table, self.org_neg_desc_table,self.org_pos_abund_table, self.org_neg_abund_table,type='org')

    def merge_mix_v2(self,pos_desc,neg_desc,pos_abund,neg_abund,type='qc'):
        if type == 'org':
            mix_path = self.check_dir(os.path.join(self.output_dir, "org_mix/metab_abund.txt"))
            mix_desc = self.check_dir(os.path.join(self.output_dir, "org_mix/metab_desc.txt"))
            self.logger.info('org: ' + mix_path)
        else:
            mix_path = self.check_dir(os.path.join(self.output_dir, "mix/metab_abund.txt"))
            mix_desc = self.check_dir(os.path.join(self.output_dir, "mix/metab_desc.txt"))
            self.logger.info('qc: ' + mix_path)
        pos_desc_data = pd.read_table(pos_desc, index_col="metab_id")
        neg_desc_data = pd.read_table(neg_desc, index_col="metab_id")

        mix_desc_data = pd.concat([pos_desc_data,neg_desc_data],axis=0,join='outer')
        mix_desc_data.to_csv(mix_desc,sep='\t', quoting=3)

        pos_abund_data = pd.read_table(pos_abund, index_col="metab_id")
        neg_abund_data = pd.read_table(neg_abund, index_col="metab_id")
        mix_abund_data = pd.concat([pos_abund_data,neg_abund_data],axis=0,join='outer')
        mix_abund_data.to_csv(mix_path,sep='\t')

    """ v1 版"""
    def merge_mix(self, pos_table, neg_table, out_pos, out_neg, new_pos_desc, new_neg_desc):
        # self.logger.info(">>>>>>in func merge_mix")
        pos_head, pos_hash, pos_stat = self.treat_data(pos_table, out_pos)
        neg_head, neg_hash, neg_stat = self.treat_data(neg_table, out_neg)
        mix_set = set(pos_hash) | set(neg_hash)
        mix_path = self.check_dir(os.path.join(self.output_dir, "mix/metab_abund.txt"))
        output = open(mix_path, 'w')
        output.write("metab_id\t%s\n" % '\t'.join(pos_head))
        pos_data = pd.read_table(new_pos_desc, index_col="metab_id").fillna("-").astype("string")
        neg_data = pd.read_table(new_neg_desc, index_col="metab_id").fillna("-").astype("string")
        mix_desc = self.check_dir(os.path.join(self.output_dir, "mix/metab_desc.txt"))
        desc_output = open(mix_desc, 'w')
        desc_output.write("metab_id\t%s\n" % '\t'.join(pos_data.columns.tolist()))
        for one in mix_set:
            if pos_hash.has_key(one) and neg_hash.has_key(one):
                if pos_stat[one]['rsd'] == neg_stat[one]['rsd']:
                    if pos_stat[one]['sum'] == neg_stat[one]['sum']:
                        output.write("%s\t%s\n" % (one, '\t'.join(pos_hash[one])))
                        desc_output.write("%s\t%s\n" % (one, '\t'.join(pos_data.loc[one].tolist())))
                    elif pos_stat[one]['sum'] > neg_stat[one]['sum']:
                        output.write("%s\t%s\n" % (one, '\t'.join(pos_hash[one])))
                        desc_output.write("%s\t%s\n" % (one, '\t'.join(pos_data.loc[one].tolist())))
                    elif neg_stat[one]['sum'] > pos_stat[one]['sum']:
                        output.write("%s\t%s\n" % (one, '\t'.join(neg_hash[one])))
                        desc_output.write("%s\t%s\n" % (one, '\t'.join(neg_data.loc[one].tolist())))
                elif pos_stat[one]['rsd'] > neg_stat[one]['rsd']:
                    # output.write("%s\t%s\n" % (one, '\t'.join(pos_hash[one]))) modify by shaohua,yuan, rsd大的删除
                    output.write("%s\t%s\n" % (one, '\t'.join(neg_hash[one])))
                    desc_output.write("%s\t%s\n" % (one, '\t'.join(neg_data.loc[one].tolist())))
                elif neg_stat[one]['rsd'] > pos_stat[one]['rsd']:
                    # output.write("%s\t%s\n" % (one, '\t'.join(neg_hash[one])))
                    output.write("%s\t%s\n" % (one, '\t'.join(pos_hash[one])))
                    desc_output.write("%s\t%s\n" % (one, '\t'.join(pos_data.loc[one].tolist())))
            elif pos_hash.has_key(one):
                output.write("%s\t%s\n" % (one, '\t'.join(pos_hash[one])))
                desc_output.write("%s\t%s\n" % (one, '\t'.join(pos_data.loc[one].tolist())))
            elif neg_hash.has_key(one):
                output.write("%s\t%s\n" % (one, '\t'.join(neg_hash[one])))
                desc_output.write("%s\t%s\n" % (one, '\t'.join(neg_data.loc[one].tolist())))
        output.close()
        desc_output.close()
        # self.logger.info(">>>>>>func merge_mix end")
        return mix_set


    def merge_mix_desc(self, new_pos_desc, new_neg_desc, mix_set):
        # self.logger.info(">>>>>>in func merge_mix_desc")
        pos_data = pd.read_table(new_pos_desc, index_col="metab_id").fillna("-").astype("string")
        neg_data = pd.read_table(new_neg_desc, index_col="metab_id").fillna("-").astype("string")
        mix_desc = self.check_dir(os.path.join(self.output_dir, "mix/metab_desc.txt"))
        output = open(mix_desc, 'w')
        output.write("metab_id\t%s\n" % '\t'.join(pos_data.columns.tolist()))
        for one in mix_set:
            if one in pos_data.index:
                output.write("%s\t%s\n" % (one, '\t'.join(pos_data.loc[one].tolist())))
            elif one in neg_data.index:
                output.write("%s\t%s\n" % (one, '\t'.join(neg_data.loc[one].tolist())))
        output.close()
        # self.logger.info(">>>>>>func merge_mix_desc end")

    def treat_data_v2(self,table, out_table):
        data = pd.read_table(table, index_col="metab_id")
        if "control" not in self.piplist:
            data.to_csv(out_table, sep="\t")
        else:
            data = data.drop(["rsd", "sum"], axis=1).astype("string")
            data.to_csv(out_table, sep='\t')



    def treat_data(self, table, out_table):
        # self.logger.info(">>>>>>in func treat_data")
        data = pd.read_table(table, index_col="metab_id")
        info_dic = {}
        stat_dic = {}
        if "control" not in self.piplist:
            data.to_csv(out_table, sep="\t")
            head = data.columns.tolist()
            for i in data.index:
                info_dic[i] = data.loc[i].tolist()
                info_dic[i]= map(str,info_dic[i])   # add by shaohua.yuan 20180920
                stat_dic[i] = {
                    "rsd": 0,
                    "sum": 0
                }
            return head, info_dic, stat_dic
        stat_data = data[["rsd", "sum"]]
        data = data.drop(["rsd", "sum"], axis=1).astype("string")
        head = data.columns.tolist()
        data.to_csv(out_table, sep='\t', quoting=3)
        for i in data.index:
            info_dic[i] = data.loc[i].tolist()
            # for i in stat_data.index:
            stat_dic[i] = {
                "rsd": stat_data["rsd"][i],
                "sum": stat_data["sum"][i]
            }
            # self.logger.info("stat_data log: %s" % i)
        #self.logger.info(">>>>>>func treat_data end")
        return head, info_dic, stat_dic

    def check_dir(self, path):
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        return path

    def set_output(self):
        if os.path.isdir(self.output_dir + '/org_pos'):
            # self.option("org_pos_out").set_path(self.output_dir + '/org_pos')
            self.option("org_pos_out", self.output_dir + '/org_pos')
        if os.path.isdir(self.output_dir + '/org_neg'):
            # self.option("org_neg_out").set_path(self.output_dir + '/org_neg')
            self.option("org_neg_out", self.output_dir + '/org_neg')
        if os.path.isfile(self.output_dir + '/metabset.list'):
            # self.option("org_set").set_path(self.output_dir + '/metabset.list')
            self.option("org_set", self.output_dir + '/metabset.list')
        if os.path.isfile(self.output_dir + '/mul_metabset.list'):
            # self.option("org_mul_set").set_path(self.output_dir + '/mul_metabset.list')
            self.option("org_mul_set", self.output_dir + '/mul_metabset.list')

        if os.path.isfile(self.output_dir + '/has_name_metabset.list'):
            self.option("has_name_org_set", self.output_dir + '/has_name_metabset.list')
        if os.path.isfile(self.output_dir + '/has_name_mul_metabset.list'):
            self.option("has_name_org_mul_set", self.output_dir + '/has_name_mul_metabset.list')
        if os.path.isfile(self.output_dir + '/has_name_origin_mul_metabset.list'):
            self.option("has_name_origin_mul_set", self.output_dir + '/has_name_origin_mul_metabset.list')

        if os.path.isdir(self.output_dir + '/pos'):
            # self.option("pos_out").set_path(self.output_dir + '/pos')
            self.option("pos_out", self.output_dir + '/pos')
        if os.path.isdir(self.output_dir + '/neg'):
            # self.option("neg_out").set_path(self.output_dir + '/neg')
            self.option("neg_out", self.output_dir + '/neg')
        if os.path.isdir(self.output_dir + '/mix'):
            # self.option("mix_out").set_path(self.output_dir + '/mix')
            self.option("mix_out", self.output_dir + '/mix')
        self.end()

    def check_ref_sample(self, table):
        value = True
        if self.option('norm') == 'sample':
            value = check_head(table, self.option("sample_name"))
        return value

    def check_ref_metab(self, table):
        value = True
        if self.option('norm') == 'inner':
            value = check_index(table, self.option("inner_ref"), "Metabolite", getid=True)
            if value:
                self.norm_method = 'inner:' + value
        return value


    def get_rid_of_inner(self, inner, infile, inner_file, new_file,search_col='Metabolite'):
        data = pd.read_table(infile, sep='\t',header=0)
        data1 = data[data[search_col] == inner]
        if len(data1) == 0:
            self.set_error("所选择的内参代谢物%s不存在，请重新输入" , variables=(self.option("inner_ref")), code="34701802")
        data2 = data[data[search_col] != inner]
        data2.to_csv(new_file,sep='\t',index=False, quoting=3)

        #inner_low = inner.lower()
        #metab_map = {inner_low: inner}
        #desc = inner_file + '_desc'
        #data1 = data.set_index('Metabolite')
        #self.write_table(data1,metab_map,inner_file,desc)
        if self.option('interactive') == 'False':
            col_nums = data1.shape[1]
            sample_col = col_nums - self.samples_num
            f1_head = data1.columns[sample_col:].tolist()
            data3 = data1[f1_head]
            data3['metab_id'] = data1[search_col]
            data3.to_csv(inner_file,sep='\t',index=False,columns=['metab_id']+f1_head, quoting=3)
        else:
            data1.to_csv(inner_file, sep='\t', index=False, quoting=3)

    def search_metabolite_id_by_desc(self, desc_file, metab_name):
        data = pd.read_table(desc_file, sep='\t',index_col=0)
        try:
            metab_id = data[data['Metabolite']==metab_name].index[0]
            return metab_id
        except Exception as e:
            return False

    def run(self):
        """
        运行
        """
        super(PreprocessTool, self).run()
        if self.option('ana_method') == 'GC' and self.option('interactive') == 'False':
            self.logger.info(">>>in GC process, interactive is false")
            if self.option('norm') == 'inner':
                infile = self.option('pos_table').path
                new_pos_table = self.work_dir + '/new_pos_table.xls'
                inner_file = self.work_dir+'/pos_inner.xls'
                self.pos_inner_file = inner_file
                inner = self.option("inner_ref")
                # if not self.check_ref_metab(infile):
                #     self.set_error("所选择的内参代谢物%s不存在，请重新输入" , variables=(self.option("inner_ref")), code="34701802")
                self.get_rid_of_inner(inner, infile, inner_file, new_pos_table)
                self.norm_method = 'inner:' + self.option("inner_ref")
                self.option('pos_table',new_pos_table)

            self.pos_table, self.pos_desc = self.run_gc_data()
            if not self.check_ref_sample(self.pos_table):
                self.set_error("所选择的内参样品%s不存在，请重新输入" , variables=(self.option('sample_name')), code="34701801")

            self.run_set()
            self.run_gc_preprocess()
            self.run_set_origin()
        elif self.option('ana_method') == 'LC' and self.option('interactive') == 'False':
            self.logger.info(">>> in LC process, interactive is false")
            if self.option('norm') == 'inner':
                pos_infile = self.option('pos_table').path
                neg_infile = self.option('neg_table').path
                self.pos_inner_file = self.work_dir+'/pos_inner.xls'
                self.neg_inner_file = self.work_dir+'/neg_inner.xls'
                new_pos_file = self.work_dir+'/new_pos_table.xls'
                new_neg_file = self.work_dir+'/new_neg_table.xls'
                inner = self.option("inner_ref")
                # if not self.check_ref_metab(pos_infile) or not self.check_ref_metab(neg_infile):
                #     self.set_error("所选择的内参代谢物%s不存在，请重新输入" , variables=(self.option("inner_ref")), code="34701802")
                self.get_rid_of_inner(inner, pos_infile, self.pos_inner_file, new_pos_file)
                self.get_rid_of_inner(inner, neg_infile, self.neg_inner_file, new_neg_file)
                self.norm_method = 'inner:' + self.option("inner_ref")
                self.option('pos_table', new_pos_file)
                self.option('neg_table', new_neg_file)

            self.pos_table, self.neg_table, self.pos_desc, self.neg_desc = self.run_lc_data()
            if not self.check_ref_sample(self.pos_table) or not self.check_ref_sample(self.neg_table):
                self.set_error("所选择的内参样品%s不存在，请重新输入" , variables=(self.option('sample_name')), code="34701803")

            self.run_set()
            self.run_lc_preprocess(pip_type=True)
            self.run_set_origin()
        elif self.option('ana_method') == 'GC' and self.option('interactive') == 'True':
            self.logger.info(">>> in GC process, interactive is true")
            self.pos_table = self.option('pos_rout').path + '/metab_abund.txt'
            self.pos_desc = self.option('pos_rout').path + '/metab_desc.txt'


            if not self.check_ref_sample(self.pos_table):
                self.set_error("所选择的内参样品%s不存在，请重新输入" , variables=(self.option('sample_name')), code="34701805")
            if not self.check_ref_metab(self.pos_desc):
                self.set_error("所选择的内参代谢物%s不存在，请重新输入" , variables=(self.option("inner_ref")), code="34701806")
            ##20200320
            if self.option('norm') == 'inner':
                infile = self.pos_table
                new_pos_table = self.work_dir + '/new_pos_table.xls'
                inner_file = self.work_dir+'/pos_inner.xls'
                self.pos_inner_file = inner_file
                inner = self.option("inner_ref")
                # if not self.check_ref_metab(infile):
                #     self.set_error("所选择的内参代谢物%s不存在，请重新输入" , variables=(self.option("inner_ref")), code="34701802")
                inner_id = self.search_metabolite_id_by_desc(self.pos_desc, inner)
                self.get_rid_of_inner(inner_id, infile, inner_file, new_pos_table, search_col= 'metab_id')
                #self.norm_method = 'inner:' + self.option("inner_ref")
                self.norm_method = 'inner:' + inner_id
                ####self.option('pos_table',new_pos_table)
                self.pos_table = new_pos_table

            self.run_gc_preprocess()
        elif self.option('ana_method') == 'LC' and self.option('interactive') == 'True':
            self.logger.info(">>> in LC process, interactive is true")
            self.pos_table = self.option('pos_rout').path + '/metab_abund.txt'
            self.neg_table = self.option('neg_rout').path + '/metab_abund.txt'
            self.pos_desc = self.option('pos_rout').path + '/metab_desc.txt'
            self.neg_desc = self.option("neg_rout").path + '/metab_desc.txt'
            self.mix_desc =  self.concat_file(self.pos_desc,self.neg_desc,'mix_desc.txt')
            self.mix_table = self.concat_file(self.pos_table, self.neg_table, 'mix_abund.txt')

            if not self.check_ref_sample(self.pos_table) or not self.check_ref_sample(self.neg_table):
                self.set_error("所选择的内参样品%s不存在，请重新输入" , variables=(self.option('sample_name')), code="34701807")
            if not self.check_ref_metab(self.pos_desc) or not self.check_ref_metab(self.neg_desc):
                self.set_error("所选择的内参代谢物%s不存在，请重新输入" , variables=(self.option("inner_ref")), code="34701808")

            ## add inner 20200320
            if self.option('norm') == 'inner':
                pos_infile = self.pos_table
                neg_infile = self.neg_table
                self.inner_file = self.work_dir+'/inner.xls'

                new_pos_file = self.work_dir+'/new_pos_table.xls'
                new_neg_file = self.work_dir+'/new_neg_table.xls'
                inner = self.option("inner_ref")
                inner_id_pos = self.search_metabolite_id_by_desc(self.pos_desc, inner)
                inner_id_neg = self.search_metabolite_id_by_desc(self.neg_desc, inner)
                if inner_id_pos:
                    self.get_rid_of_inner(inner_id_pos, pos_infile, self.inner_file, new_pos_file)
                    inner_id = inner_id_pos
                    self.pos_table = new_pos_file
                elif inner_id_neg:
                    self.get_rid_of_inner(inner_id_neg, neg_infile, self.inner_file, new_neg_file)
                    inner_id = inner_id_neg
                    self.neg_table = new_neg_file
                else:
                    self.set_error("所选择的内参代谢物%s不存在，请重新输入" , variables=(self.option("inner_ref")), code="34701808")
                self.norm_method = 'inner:' + inner_id

            self.run_lc_preprocess()
        ##v3 202003 增加 hmdb等的统计和 cv画图的数据
        if self.option('ana_method') == 'GC':
            if self.option("task_version") !="1.0":
                self.count_fun(pos_desc=self.output_dir+'/pos/metab_desc.txt')
                self.qc_group_cv(pos_abund=self.output_dir+'/pos/metab_abund.txt', has_rsd_txt=True)
                if self.option("raw_cv") == 'T':
                    if self.option('interactive') == 'True':
                        self.count_fun(pos_desc=self.pos_desc,out_name='raw_')
                        self.qc_group_cv(pos_abund=self.pos_table,out_name='raw_')
                    else:
                        self.count_fun(pos_desc=self.output_dir+'/org_pos/metab_desc.txt',out_name='raw_')
                        self.qc_group_cv(pos_abund=self.output_dir+'/org_pos/metab_abund.txt',out_name='raw_')
        else:
            if self.option("task_version") !="1.0":
                self.count_fun(pos_desc=self.output_dir+'/pos/metab_desc.txt',
                               neg_desc=self.output_dir+'/neg/metab_desc.txt',
                               mix_desc=self.output_dir+'/mix/metab_desc.txt')
                self.qc_group_cv(pos_abund=self.output_dir+'/pos/metab_abund.txt',
                                 neg_abund=self.output_dir+'/neg/metab_abund.txt',
                                 mix_abund=self.output_dir+'/mix/metab_abund.txt',
                                 has_rsd_txt=True)
                if self.option("raw_cv") == 'T':
                    if self.option('interactive') == 'True':

                        self.count_fun(pos_desc=self.pos_desc,
                                       neg_desc=self.neg_desc,
                                       mix_desc = self.mix_desc ,
                                       out_name='raw_')
                        self.qc_group_cv(pos_abund=self.pos_table,
                                         neg_abund=self.neg_table,
                                         mix_abund= self.mix_table,
                                         out_name='raw_')
                    else:
                        self.count_fun(pos_desc=self.output_dir+'/org_pos/metab_desc.txt',
                                       neg_desc=self.output_dir+'/org_neg/metab_desc.txt',
                                       mix_desc=self.output_dir+'/org_mix/metab_desc.txt',
                                       out_name='raw_')
                        self.qc_group_cv(pos_abund=self.output_dir+'/org_pos/metab_abund.txt',
                                         neg_abund=self.output_dir+'/org_neg/metab_abund.txt',
                                         mix_abund=self.output_dir+'/org_mix/metab_abund.txt',
                                         out_name='raw_')

        self.set_output()

    def concat_file(self,f1,f2,new_file):
        f1_data = pd.read_table(f1,sep='\t')
        f2_data = pd.read_table(f2,sep='\t')
        cat_data = pd.concat([f1_data,f2_data],axis=0,join='outer')
        cat_data.to_csv(new_file,sep='\t', quoting=3)
        return new_file

    #统计hmdb kegg all等离子数 v3 增加 202003.输入文件是质控后的结果  zouguanqing
    def count_fun(self,pos_desc=None,neg_desc=None,mix_desc=None,out_name=''):
        tmp_dic = {}
        if pos_desc:
            tmp_dic['pos'] = pos_desc
        if neg_desc:
            tmp_dic['neg'] = neg_desc
        if mix_desc:
            tmp_dic['mix'] = mix_desc
        type_l = []
        hmdb_l = []
        kegg_l = []
        all_l = []
        name_l = []


        for type in ['pos','neg','mix']:
            if type not in tmp_dic:
                continue
            desc_file = tmp_dic[type]
            data = pd.read_table(desc_file, sep='\t')


            type_l.append(type)
            all_l.append(data.shape[0])
            if 'Library ID'  in data.columns:
                hmdb_col = 'Library ID'
                hmdb_l.append(data[data[hmdb_col] != "-"].shape[0])
            elif 'HMDB_ID' in data.columns:
                hmdb_col = 'HMDB_ID'
                hmdb_l.append(data[data[hmdb_col] != "-"].shape[0])
            else:
                hmdb_l.append(0)
            kegg_l.append(data[data['KEGG Compound ID']!="-"].shape[0])
            name_l.append(data[data['Metabolite'].map(lambda x: False if  re.match('metab_\d+$',x) else True)].shape[0])
        out_df = pd.DataFrame({'type':type_l, "all": all_l,"has_name": name_l,'hmdb':hmdb_l,'kegg':kegg_l})
        out_df.to_csv(self.output_dir+'/{}statistic.xls'.format(out_name),sep='\t',index=False)

    #qc组计算 cv    zouguanqing
    def qc_group_cv(self,pos_abund=None, neg_abund=None, mix_abund=None,out_name='',has_rsd_txt=False):
        if len(self.qc_samples) < 2:
            return False

        tmp_dic = {}
        if pos_abund:
            tmp_dic['pos'] = pos_abund
        if neg_abund:
            tmp_dic['neg'] = neg_abund
        if mix_abund:
            tmp_dic['mix'] = mix_abund


        if has_rsd_txt:
            for type in ['pos', 'neg']:
                if type in tmp_dic:
                    if os.path.exists(self.output_dir+ '/'+type+'/metab_abund.txt_rsd.txt'):
                        os.remove(self.output_dir+ '/'+type+'/metab_abund.txt_rsd.txt')
                    os.link(self.work_dir+'/'+type+'/metab_abund.txt_rsd.txt', self.output_dir+ '/'+type+'/metab_abund.txt_rsd.txt')

            if mix_abund:
                mix_rsd = mix_abund+ '_rsd.txt'
                f1_data = pd.read_table(pos_abund+'_rsd.txt',sep='\t',index_col=0,header=None)
                f2_data = pd.read_table(neg_abund+'_rsd.txt',sep='\t',index_col=0,header=None)
                cat_data = pd.concat([f1_data,f2_data],axis=0,join='outer')
                cat_data.to_csv(mix_rsd,sep='\t')

        per_l = ['per'+str(i) for i in range(1,11)]

        res = {}

        if has_rsd_txt == False:
            for type in tmp_dic:
                res[type] = []
                abund_file = tmp_dic[type]
                data = pd.read_table(abund_file, sep='\t',index_col='metab_id')
                metab_num = float(data.shape[0])
                if len(self.qc_samples) !=0:
                    qc_data = data[self.qc_samples]
                    qc_data['rsd'] = qc_data.apply(lambda x: x.std()/x.mean(), axis=1)
                    for per in range(1,11):
                        per_c = per*0.1
                        res[type].append(qc_data[qc_data['rsd']<per_c].shape[0]/metab_num)
                    qc_data['rsd'].to_csv(abund_file+'_rsd.txt', sep='\t')
        else:
            for type in tmp_dic:
                res[type] = []
                rsd_file = tmp_dic[type] + '_rsd.txt'
                data = pd.read_table(rsd_file, sep='\t',index_col=0,header=None)
                tc = data.columns[0]
                metab_num = float(data.shape[0])
                for per in range(1,11):
                        per_c = per*0.1
                        res[type].append(data[data[tc]<per_c].shape[0]/metab_num)

        out_df = pd.DataFrame({'percent':per_l})
        for type in res:
            out_df[type] = res[type]

        out_df.to_csv(self.output_dir+'/{}cv_percent.xls'.format(out_name), sep='\t', index=False)
