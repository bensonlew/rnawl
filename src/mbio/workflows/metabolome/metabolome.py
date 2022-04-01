# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# __last_modified__ = '20180620'

"""代谢分析工作流"""

from mainapp.libs.param_pack import group_detail_sort
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
from biocluster.config import Config
from biocluster.api.database.base import ApiManager
import pandas as pd
import numpy as np
import os, glob
import json
import shutil
import time
import datetime
import gevent
import functools
from copy import deepcopy
from  mainapp.models.mongo.metabolome import  Metabolome
import re
from mbio.packages.meta.delete_mongo import DeleteDemoMongo

def time_count(func):  # 统计导表时间
    @functools.wraps(func)
    def wrapper(self, *args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run %s at %s' % (func_name, start_time))
        if func_name != "export_group_and_samples" and func_name != "export_org_metabset" and func_name != "export_diff_metabset":
            self.main_collection_delete(func_name)
        func(self, *args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，该阶段导表已进行{}s".format(func_name, end - start))
    return wrapper

def tryforgood(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

class MetabolomeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        代谢workflow option参数设置
        """
        self._sheet = wsheet_object
        super(MetabolomeWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'test', 'type': 'bool', 'default': False},  # 是否为测试workflow
            {'name': 'project_type', 'type': 'string', 'default': '', 'required': True},  # GC、LC类型
            {'name': 'pos_table', 'type': 'infile', 'format': 'metabolome.metab_table,sequence.profile_table',
             'required': True},  # 阳离子输入表
            {'name': 'neg_table', 'type': 'infile', 'format': 'metabolome.metab_table,sequence.profile_table'},
            # 阴离子输入表，LC使用
            {'name': 'group_table', 'type': 'infile', 'format': 'meta.otu.group_table', 'required': True},  # 分组文件
            {'name': 'diff_group', 'type': 'infile', 'format': 'sample.control_table', 'required': True},
            {'name': 'other_organism', 'type': 'string'},
            {'name': 'organism_type', 'type': 'string', 'default': '', 'required': True},  # KEGG物种类型
            {'name': 'organism', 'type': 'string', 'default': ''},  # KEGG物种名称
            {'name': 'fillna', 'type': 'string', 'default': 'none'},  # 缺失值填充方法: min/median/mean/rf/none
            {'name': 'rsd', 'type': 'string', 'default': '30'},  # QC验证RSD，需数字+% 数字范围 25-30
            {'name': 'norm', 'type': 'string', 'default': 'none'},  # 数据归一化方法：median/mean/sum/sample/inner/none
            {'name': 'sample_name', 'type': 'string'},  # 标准化样品， norm = sample时使用
            {'name': 'inner_ref', 'type': 'string'},  # 内参代谢物，norm = inner时使用
            {'name': 'scale', 'type': 'string', 'default': 'none'},  # 取log值，形式log+数字/"none"/"defined"
            {'name': 'log', 'type': 'int'},  # 取log值，当scale = undefined时使用,
            {'name': 'hmdb', 'type': 'string', 'default': "T"},
            #{'name': 'rm_nan','type':'float','default':50}, #20190603 v2.0
            {'name': 'fill_percent', 'type': 'int', 'default': 50},
            {'name': 'fill_type', 'type': 'string', 'default': 'all'},  #add v3 20200326
            {'name': 'diff_check', 'type': 'string', 'default': 'Student'}, #Student,Welch,Wilcox  # "t-test", 'welch', 'wilcox'
            {'name': "p_fdr",'type': 'string', 'default': ''},  #P_value,FDR
            {'name': "p_fdr_condition", 'type': 'string','default' : '<'},  #＜,≤,=,＞,≥
            {'name': 'p_fdr_value','type': 'string', 'default': '0.05'},
            {'name': 'vip', 'type': 'string', 'default':''}, #VIP_PLS-DA,VIP_OPLS-DA
            {'name': 'vip_condition','type':'string', 'default':'<'}, #＜,≤,=,＞,≥
            {'name': 'vip_value', 'type':'float', 'default' : 0},
            {'name': 'up_condition', 'type': 'string', 'default': ""}, # ＞,≥
            {'name': 'up_value', 'type': 'float','default':1}, #
            {'name': 'down_condition', 'type': 'string', 'default':''},
            {'name': 'down_value', 'type': 'float', 'default':1},
            {'name': 'pca_data','type': 'string', 'default': 'UV'}, #UV,Ctr,Par,none
            {'name': 'plsda_data', 'type': 'string', 'default': 'UV'}, #UV,Ctr,Par,none
            {'name': 'oplsda_data', 'type': 'string', 'default': 'UV'}, #UV,Ctr,Par,none
            {'name': 'mix_table', 'type': 'string', 'default':'F'}, # LC时是否合并阴阳离子表 T or F
            {'name': 's2_diff_check', 'type': 'string', 'default':'fisher'},  #add v3
            {'name': 's2_p_value_fdr', 'type': 'string', 'default': 'P_value'},
            {'name': 's2_p_fdr_condition', 'type': 'string','default': '$gt'},
            {'name': 's2_p_fdr_value', 'type': 'float'},
            {'name': 's2_up_condition', 'type': 'string','default': '$gt'},
            {'name': 's2_up_value', 'type': 'float', 'default': 1.2},
            {'name': 's2_down_condition', 'type': 'string', 'default': '$lt'},
            {'name': 's2_down_value', 'type': 'float', 'default': 0.8},
            {'name': 'kegg_version', 'type': 'string', 'default': "KEGG"},## kegg数据库注释版本
            {'name': 'save_pdf', 'type': 'int', 'default': 0},
            {'name': 'tsanger', 'type': 'string', 'default': ''},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.is_mix_table = self.option('mix_table')
        self.check_group_conditon()  #v3  检查是否有 2样本的差异组 ,2样本以上的 组是否大于2组， 是否全是单样本组

        # 添加module和tool
        if self.option("project_type") == "GC":
            self.run_type1 = ["pos"]  # pca, 差异分析用
            self.run_type2 = ["pos"]  # 其他分析用
            self.run_type3 = ["pos"]
        else:
            if self.is_mix_table == 'T':
                self.run_type3 = ['mix']  #diff_pls,metabset_vip
            else:
                self.run_type3 = ["pos","neg"]
            self.run_type1 = ["pos","neg"]
            self.run_type2 = ["mix"]


        if self.option('diff_check') == 'Student':
            self.option('diff_check','t-test')
        elif self.option('diff_check') == 'Welch':
            self.option('diff_check', 'welch')
        elif self.option('diff_check') == 'Wilcox':
            self.option('diff_check', 'wilcox')

        self.task_id = self._sheet.id
        ## 不涉及数据表类型
        self.preprocess = self.add_tool('metabolome.preprocess')
        self.metabset_venn = self.add_tool('metabolome.metabset.venn')
        self.metabset_keggc = self.add_tool('metabolome.metabset.keggc')
        self.metabset_keggp = self.add_tool('metabolome.metabset.keggp')
        #self.metabset_enrich = self.add_tool('metabolome.metabset.enrich')
        self.metabset_enrich = self.add_module('metabolome.enrich_topo')
        #self.metabset_ipath = self.add_tool('metabolome.metabset.ipath')
        self.metabset_ipath = self.add_tool('metabolome.metabset.ipath3')
        self.creat_analysis_table = self.add_tool('metabolome.creat_table')
        self.diff_metabset_pip = self.add_module("metabolome.metabsets")  #差异代谢集的多个分析的集成modoule v2.0添加

        self.pre_metabset_tools = []  # 代谢集前分析，主要包括注释、pca、相关性和差异分析

        self.metabset_tools = [self.metabset_venn, self.metabset_keggc, self.metabset_keggp, self.metabset_enrich,
                               self.metabset_ipath]  # 代谢集分析
        ## 涉及数据表，LC和GC tool个数使用可能不一致
        self.add_lc_tools(self.run_type2, "annotation", "metabolome.anno", tool_type="module")
        #self.add_lc_tools(self.run_type1, "diff_pls", "metabolome.diff_pls", tool_type="module")
        if len(self.run_type3)== 2:   ## 即pos和neg 2种情况，生成mix的结果。用于代谢集的vip分析
            #self.diff_pls_types  = ['mix']
            #self.diff_pls_types.extend(self.run_type3)
            self.diff_pls_types = self.run_type3
        else:
            self.diff_pls_types = self.run_type3

        if self.more2group_group:
            self.add_lc_tools(self.run_type3, "groups_diff","metabolome.groups_diff", tool_type='module')  #add v3
        if self.vs_more_diff:
            self.add_lc_tools(self.diff_pls_types, "diff_pls", "metabolome.diff_pls", tool_type="module")  #20190722 run_type1 change to run_type3
        if self.s1vs1_diff:
            self.add_lc_tools(self.run_type3, "two_sample", "metabolome.two_sample_diff", tool_type="module") #add v3

        self.add_lc_tools(self.run_type3, "pca", "metabolome.diff.diff_mul_stat")
        self.add_lc_tools(self.run_type3, "sample_corr", "metabolome.compare.corr_tree")
        self.add_lc_tools(self.run_type3, "sample_venn", "metabolome.sample_venn", tool_type="module")  #add v3
        self.add_lc_tools(self.run_type3, "sample_plsda", "metabolome.diff.diff_mul_stat") #add v3
        self.add_lc_tools(self.run_type2, "metabset_cluster", "metabolome.metabset.metab_cluster")
        #self.add_lc_tools(self.run_type1, "metabset_vip", "metabolome.metabset.metab_vip")
        if self.vs_more_diff:
            self.add_lc_tools(self.run_type2, "metabset_vip", "metabolome.metabset.metab_vip")  ##run_type1 change to run_type3
        self.add_lc_tools(self.run_type2, "metabset_corr", "metabolome.compare.corr_tree")
        if self.option("hmdb") == "T":
            self.add_lc_tools(self.run_type2, "anno_hmdb", "metabolome.annotation.anno_hmdb")
            self.metabset_hmdb = self.add_tool('metabolome.annotation.anno_hmdb')
            self.metabset_tools.append(self.metabset_hmdb)
        # 定义公共参数
        self.project_type = self.option("project_type")
        self.group_table = self.option("group_table")
        self.raw_pos_table = self.option("pos_table").prop["path"]

        '''add_steps'''
        self.step.add_steps('preprocess', 'annotation', 'exp_pca', 'sample_corr', "sample_venn","sample_plsda",
                            'creat_table', 'metabset')
        if self.more2group_group:
            self.step.add_steps('groups_diff')

        if self.vs_more_diff:
            self.step.add_steps('diff_pls')

        if self.s1vs1_diff:
            self.step.add_steps('two_sample')

        '''初始化自定义变量'''

        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.api_dic = {}  # 存放導表函數
        self.group_detail = {}  # 分组对应样品id{group1: [id1,id2,id3], group2: [id4,id5,id6]}
        self.metabset_analysis_file = ""
        self.start = {}  ## set_step 使用
        self.venn_sets = {}
        self.default_top = 50
        self.hmdb_overview = ''
        self.metabolome_api = Metabolome()  # 20190717
        '''
        sanger_type, sanger_path = self._sheet.output.split(':')
        self.logger.info(self._sheet.output)
        self.logger.info(sanger_path)
        sanger_prefix = Config().get_netdata_config(sanger_type)
        self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        self.my_bucket = Config().get_project_region_bucket("metabolome")
        '''
        self.objid2dir = {}  # 每一个主表ID对应的结果文件夹相对路径
        # output 结果 转到 upload文件夹的对应文件夹名称
        self.link_up_map = {
            'Preprocess': '1.Preprocess',
            'ExpCorr' : '2.SampleComp/01.ExpCorr',
            'ExpPCA': '2.SampleComp/02.ExpPCA',
            'ExpVenn' : "2.SampleComp/03.ExpVenn",   #add v3
            'ExpPLSDA' :"2.SampleComp/04.ExpPLSDA",  #add v3
            'AnnoKeggc': '3.Anno/01.AnnoKeggc',
            'AnnoKeggp' : '3.Anno/02.AnnoKeggp',
            'AnnoHmdb' : '3.Anno/03.AnnoHmdb',
            'AnnoOverview': '3.Anno/04.AnnoOverview',
            #'ExpDiff':'4.ExpDiff',
            'ExpDiff' : '4.ExpDiff/01.TwoGroupExpDiff',
            'TwoSample' : '4.ExpDiff/02.TwoSamExpDiff',  #add v3
            'GroupsDiff' : '4.ExpDiff/03.MultiGroupExpDiff',  #add v3
            'MetabsetVenn':'5.Metabset/01.MetabsetVenn',
            'MetabsetCluster':'5.Metabset/02.MetabsetCluster',
            'MetabsetVip':'5.Metabset/03.MetabsetVip',
            'MetabsetCorr': '5.Metabset/04.MetabsetCorr' ,
            'MetabsetKeggc':'5.Metabset/05.MetabsetKeggc',
            'MetabsetKeggp':'5.Metabset/06.MetabsetKeggp',
            'MetabsetEnrich':'5.Metabset/07.MetabsetEnrich',
            'MetabsetHmdb':'5.Metabset/08.MetabsetHmdb',
            'MetabsetIpath':'5.Metabset/09.MetabsetIpath',
            'MetabsetRoc':'5.Metabset/10.MetabsetRoc',
            'MetabsetIntergratedCorr':'5.Metabset/19.IntergratedCorr',
            'MetabsetIntergratedProc':'5.Metabset/20.IntergratedProc',
        }
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'metabolome')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    ## 动态添加变量名
    def add_lc_tools(self, run_type, analysis, tool, tool_type="tool"):
        for each in run_type:
            toolname = analysis + "_" + each
            self.logger.info(toolname)
            if tool_type == "tool":
                globals()[toolname] = self.add_tool(tool)
            elif tool_type == "module":
                globals()[toolname] = self.add_module(tool)
            if analysis in ["pca", "annotation", "sample_corr", "diff_pls","sample_venn", "sample_plsda","groups_diff","two_sample"]:  #v3 add sample_venn sample_plsda, two_sample
                self.pre_metabset_tools.append(globals()[toolname])
                self.logger.info(self.pre_metabset_tools)
            elif analysis == "anno_hmdb":
                self.logger.info("add anno hmdb tool")
            else:
                self.metabset_tools.append(globals()[toolname])

            if toolname in ['diff_pls_pos', 'diff_pls_mix'] :
                self.diff_pls_tool = globals()[toolname]
            elif toolname in ['diff_pls_neg']:
                self.diff_pls_tool_neg =  globals()[toolname]


    def check_options(self):
        """
        检查参数设置
        """
        self.logger.info("check_sheet_data...")
        if not self._sheet.id:
            raise OptionError('需要输入task_id', code="14700101")
        if not self._sheet.member_type:
            raise OptionError("需要输入member_type", code="14700102")
        if not self._sheet.cmd_id:
            raise OptionError("需要输入cmd_id", code="14700103")
        self.logger.info("check options...")
        if self.option("project_type") == "LC":
            if not self.option("neg_table").is_set:
                raise OptionError("当项目类型为LC时必须输入负离子表格！", code="14700104")
        if self.option("fillna") not in ["min", "median", "mean", "rf", "none"]:
            raise OptionError("缺失值填充方法错误，只能为min, median,mean，rf，none！", code="14700105")
        if self.option('rsd') != "none":
            try:
                rsd_value = int(self.option('rsd'))
            except:
                raise OptionError("rsd type must be int or string 'none'!")
            if int(self.option('rsd')) < 25 or int(self.option('rsd')) > 30:
                raise OptionError("rsd参数必须在25%-30%范围内", code="14700106")
        if self.option("norm") not in ["median", "mean", "sum", "sample", "inner", "none"]:
            raise OptionError("数据归一化方法，只能为median,mean,sum,sample,inner,none！", code="14700107")
        elif self.option('norm') == 'sample':
            if not self.option('sample_name'):
                raise OptionError("选择内参样品时，需指定样品", code="14700108")
            self.check_norm(sample=self.option('sample_name'))
        elif self.option('norm') == 'inner':
            if not self.option('inner_ref'):
                raise OptionError("须填写内参代谢物", code="14700109")
            self.check_norm(inner=self.option('inner_ref'))
        self.logger.info(self.option('scale'))
        if self.option('scale') not in ["none", "defined"]:
            if not self.option('scale').startswith("log"):
                raise OptionError("scale 参数必须为log+数字，或者为none，defined", code="14700110")
            else:
                check_int = self.option('scale').lstrip("log")
                try:
                    check_int = int(check_int)
                except:
                    raise OptionError("scale 参数必须为log + 整数，或者为none, defined", code="14700111")
                if check_int < 2:
                    raise OptionError("scale 参数的整数部分必须大于等于2", code="14700112")
        elif self.option('scale') == 'defined':
            if self.option('log') < 2:
                self.logger.info(self.option('log'))
                raise OptionError("log参数值必须大于等于2的整数", code="14700113")
        return True

    def check_specimen_names(self):
        """
        检查丰度表中是否存在group_file中对应样本名
        """
        self.logger.info("check group file...")
        other_header = ["Metabolite", "Formula", "Mass", "RT (min)", "KEGG Compound ID", "HMDB_ID",
                        "CAS number", "Mode", 'm/z']
        raw_pos_file = self.option("pos_table").prop["path"]
        with open(raw_pos_file, "r") as f1:
            pos_header = f1.next().strip().split("\t")
            for item in other_header:
                if item in pos_header:
                    pos_header.remove(item)
        self.logger.info(pos_header)
        if self.option("project_type") == "LC":
            raw_neg_file = self.option("neg_table").prop["path"]
            with open(raw_neg_file, "r") as f2:
                neg_header = f2.next().strip().split("\t")
                for item in other_header:
                    if item in neg_header:
                        neg_header.remove(item)
            if sorted(pos_header) != sorted(neg_header):
                self.logger.error("正负离子表的样本不匹配！")
                raise OptionError("正负离子表的样本不匹配！", code="14700114")
        #if sorted(pos_header) != sorted(self.get_group()):
        for g_sample in sorted(self.get_group()):
            if g_sample not in pos_header:
                self.logger.info(sorted(pos_header))
                self.logger.info(sorted(self.get_group()))
                self.logger.error("离子表和group表的样本不匹配！")
                raise OptionError("离子表和group表的样本不匹配！", code="14700115")

    # def get_diff_group(self):
    #     """
    #     获取对照组组名
    #     """
    #     self.logger.info("get_diff_group...")
    #     group_diff_name = []
    #     diff_file = self.option("diff_group").prop["path"]
    #     with open(diff_file, "r") as dg:
    #         for line in dg:
    #             line = line.strip().split("\t")
    #             if "#" not in line[0]:
    #                 g1 = line[0]
    #                 g2 = line[1]
    #                 diffname = g2 + "_vs_" + g1
    #                 group_diff_name.append(diffname)
    #     group_diff_str = ";".join(group_diff_name)
    #     return group_diff_str

    def get_group(self):
        samples = self.option("group_table").prop["sample"]
        for eachsample in samples:
            if "-" in eachsample:
                self.logger.error("请将样品名中中划线改为下划线！")
                raise OptionError("请将样品名中中划线改为下划线！", code="14700116")
        return samples

    def get_noQC_group(self):
        sub_group_file = self.work_dir + "/" + "noQC_group.xls"
        group_file = self.option("group_table").prop["path"]
        with open(group_file, "r") as f, open(sub_group_file, "w") as outf:
            for line in f:
                line = line.strip().split("\t")
                sample = line[0]
                group = line[1]
                if not group == "QC":
                    outf.write(sample + "\t" + group + "\n")
        self.no_qc_group_table = sub_group_file


    def check_group(self):
        cmp_list = self.option("diff_group").prop["cmp_list"]
        group_dict = self.option("group_table").get_group_spname()
        self.logger.info(cmp_list)
        group_list = self.option('group_table').get_group_spname()
        for each in cmp_list:
            for eachname in each:
                if eachname not in group_list:
                    self.logger.error("Control table 中group--{} 不在Group 表中，请检查！".format(eachname))
                    raise OptionError("Control table 中group--%s 不在Group 表中，请检查！", variables=(eachname), code="14700117")
                else:
                    if "_vs_" in each:
                        self.logger.error("分组名中不允许带'_vs_',会与差异分组冲突，请修改！")
                        raise OptionError("分组名中不允许带'_vs_',会与差异分组冲突，请修改！", code="14700118")
        for eachgroup in group_dict:
            # if eachgroup != "QC":
            #     len_sam = len(group_dict[eachgroup])
            #     if len_sam < 3:
            #         self.logger.error("除QC外每组需含有三个以上样本，组--{}不满足！".format(eachgroup))
            #         ###raise OptionError("除QC外每组需含有三个以上样本，组--%s不满足！", variables=(eachgroup), code="14700119")
            if eachgroup == "QC":
                len_sam = len(group_dict[eachgroup])
                if len_sam < 2:
                    raise OptionError("if has QC, QC samples must >1！")

        #检查表达文件中有QC组, 但分组文件没有QC组
        if 'QC' not in group_dict:
            raw_pos_file = self.option("pos_table").prop["path"]
            with open(raw_pos_file) as f:
                heads = f.readline().strip().split('\t')
            for h in heads:
                if re.match('QC\d+$',h):
                    raise OptionError("丰度表有QC组，但分组文件没有QC组！")


    ##检查差异组是否存在2样本的情况, 2样本比较情况，多组比较的情况
    def check_group_conditon(self):
        vs_more_diff = []
        s1vs1_diff = []
        self.s1vs1_group_detail = {}  # 做2样本比较的group detail
        self.vs_more_group_detail = {}   #做2组比较的group detail
        self.more2_group_detail = {}    #做多组方差分析的group detail
        diff_file = self.option("diff_group").prop["path"]
        group_file = self.option("group_table").prop["path"]
        group = pd.read_table(group_file,sep='\t')
        group.columns = ['sample','group']
        group_detail = {}
        tmp = dict(list(group.groupby(by=['group'])))
        for k in tmp.keys():
            group_detail[k] = sorted(tmp[k]['sample'].tolist())

        diff = pd.read_table(diff_file, sep='\t')
        col1 = diff.columns[0]
        col2 = diff.columns[1]
        for i in diff.index:
            g1 = diff.loc[i][col1]
            g2 = diff.loc[i][col2]
            if g1 not in group_detail:
                self.set_error("%s 不在分组文件中"%g1)
            elif g2 not in group_detail:
                self.set_error("%s 不在分组文件中"%g2)

            if len(group_detail[g1])==1 and len(group_detail[g2])==1:
                diff_name = g2 + '_vs_' + g1
                s1vs1_diff.append(diff_name)
                if g1 not in self.s1vs1_group_detail:
                    self.s1vs1_group_detail[g1] = group_detail[g1]
                if g2 not in self.s1vs1_group_detail:
                    self.s1vs1_group_detail[g2] = group_detail[g2]
            else:
                diff_name = g2 + '_vs_' + g1
                vs_more_diff.append(diff_name)
                if g1 not in self.vs_more_group_detail:
                    self.vs_more_group_detail[g1] = group_detail[g1]
                if g2 not in self.vs_more_group_detail:
                    self.vs_more_group_detail[g2] = group_detail[g2]


        ###记录 满足多vs 1或1vs多， 或多比多的 差异组. 用于2组差异分析
        if vs_more_diff:
            self.vs_more_diff = ';'.join(vs_more_diff)
        else:
            self.vs_more_diff = ''

        ##记录比较组的2组都是1个样本
        if s1vs1_diff:
            self.s1vs1_diff = ';'.join(s1vs1_diff)
        else:
            self.s1vs1_diff = ''

        ####
        self.raw_group_detail = group_detail

        ##检查组内样本数大于2的组，生成group 文件用于多组比较分析
        more_2_group = []
        for g in group_detail:
            if len(group_detail[g]) > 2:
                if g not in ['QC']:
                    more_2_group.append(g)
                    self.more2_group_detail[g] = group_detail[g]

        if len(more_2_group) > 2:
            self.more2group_group = self.work_dir+'/groups_diff_group.txt'  #组样本数大于2的，这种组的数量大于2组，的分组文件
            with open(self.more2group_group, 'w') as fw:
                fw.write('#sample\tgroup\n')
                for g in more_2_group:
                    for s in group_detail[g]:
                        fw.write(s+'\t'+g+'\n')

        else:
            self.more2group_group = ''



    def check_norm(self, sample=None, inner=None):
        if sample:
            samples = self.get_group()
            if not sample in samples:
                raise OptionError("内参样本-%s-不 在group样本表范围内！", variables=(sample), code="14700120")
        if inner:
            raw_pos_file = self.option("pos_table").prop["path"]
            inner_table = pd.read_table(raw_pos_file, sep="\t", header=0, index_col='Metabolite')
            metab = inner_table.index.tolist()
            if not inner in metab:
                raise OptionError("内参代谢物-%s-不在丰度表内！", variables=(inner), code="14700121")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.logger.info(event['data'])
        self.logger.info(event['data'].values())
        if event['data'].values() == "creat_table":
            self.step.annotation.finish()
            self.step.exp_pca.finish()
            self.step.sample_corr.finish()
            self.step.diff_pls.finish()
        self.step.update()

    def set_step_stat(self, opts, module, event, step=None, start=False, run=False):
        #if step != self.step.metabset :
        #    module.on('start', self.set_step, {'start': step})
        if step:
            module.on('start', self.set_step, {'start': step})
            module.on('end', self.set_step, {'end': step})

        module.on('end', self.set_output, event)
        if start:
            module.set_options(opts)
        if run:
            module.run()



    def run_preprocess(self):
        self.logger.info("-----start run_preprocess----")
        if self.option("rsd") == "none":
            rsd = "none"
        else:
            rsd = self.option("rsd") + "%"
        opts = {
            'ana_method': self.project_type,
            'pos_table': self.raw_pos_table,
            'group_table': self.group_table,
            'fillna': self.option("fillna"),
            'rsd': rsd,
            'norm': self.option("norm"),
            'scale': self.option('scale'),
            'rm_nan': self.option('fill_percent'),  #20190603 v2.0
            'fill_type' :self.option('fill_type')
        }
        if self.project_type == "LC":
            opts['neg_table'] = self.option("neg_table").prop["path"]
        if self.option('norm') == 'sample':
            opts['sample_name'] = self.option("sample_name")
        if self.option('norm') == 'inner':
            opts['inner_ref'] = self.option("inner_ref")
        if self.option('scale') == 'defined':
            opts['log'] = self.option("log")
        self.set_step_stat(opts, self.preprocess, 'preprocess', self.step.preprocess, start=True, run=True)

    def run_annotation(self):
        self.logger.info("-----start annotation----")
        if self.option("organism"):
            kegg_species = self.option("organism_type") + ";" + self.option("organism")
        elif self.option("organism_type") != "All":
            kegg_species = self.option("organism_type")
        else:
            kegg_species = "False"
        org_set = self.preprocess.option("org_set")
        org_metab_pos = self.preprocess.option("org_pos_out").prop["metab_desc"]
        opts = {
            #'metab_table': metab_table,
            'organism': kegg_species,
            "org_metab_pos": org_metab_pos,
            #"org_metab_neg": org_metab_neg,
            'anno_list': 'keggc,keggp,overview',
            'metabset': org_set,
            'work_type': self.project_type
        }
        if self.option("project_type") == "LC":
            org_metab_neg = self.preprocess.option("org_neg_out").prop["metab_desc"]
            opts["org_metab_neg"] = org_metab_neg
        self.add_pos_neg_mix_tool(opts, self.run_type2, 'annotation', self.step.annotation, parm_file="metab_table",
                                  file_name="metab_desc.txt")
    def run_anno_hmdb(self):
        self.logger.info("-----start annotation HMDB----")
        if self.option("project_type") == "LC":
            annotation_overview = annotation_mix.option("overview_out").prop["path"]
        else:
            annotation_overview = annotation_pos.option("overview_out").prop["path"]
        self.logger.info(annotation_overview)
        opts = {
            "anno_overview": annotation_overview,
            "type": "anno"
        }
        self.add_pos_neg_mix_tool(opts, self.run_type2, 'anno_hmdb', self.step.annotation, parm_file="anno_hmdb",
                                  file_name="")

    def run_pca(self):
        self.logger.info("-----start pca----")
        opts = {
            #'exp_file': exp_file,
            "group_file": self.group_table,
            'mul_type': 'pca',
            'confidence': '0.95',
            'data_trans': 'UV',
        }
        self.add_pos_neg_mix_tool(opts, self.run_type3, 'pca', self.step.exp_pca, parm_file="pca")

    #阴阳离子情况 添加表达量表 或描述表
    def add_pos_neg_mix_tool(self, opts, run_type, analysis, step=None, parm_file="exp_file", file_name="metab_abund.txt"):
        """
        动态添加tool的set_option和set_step
        """
        for each in run_type:
            tmp_option = each + "_out"
            toolname = analysis + "_" + each
            if analysis == "annotation":
                db_v = self.database_version()
                self.logger.info("start old kegg database：{}".format(db_v))
                opts['database_version'] = db_v ## add by qingchen.zhang用于区分新老工作流
            if parm_file == "diff_dir":   # 只有vip 用
                ##################
                if len(self.run_type3) == 2:
                    diff_pos_tool = "diff_pls_pos"
                    diff_neg_tool = "diff_pls_neg"
                    diff_dir = globals()[diff_pos_tool].option("id_diffStat_dir").prop['path']
                    diff_dir2 = globals()[diff_neg_tool].option("id_diffStat_dir").prop['path']
                    all_dir = self.work_dir + '/tmp_all_diff_dir/'
                    if not os.path.exists(all_dir):
                        os.mkdir(all_dir)
                    files = os.listdir(diff_dir)
                    files2 = os.listdir(diff_dir2)
                    share_file = set(files) & set(files2)
                    for f in share_file:
                        d1 = pd.read_table(diff_dir+'/'+f,header=0)
                        d2 = pd.read_table(diff_dir2+'/'+f,header=0)
                        data_cat = pd.concat([d1,d2],axis=0)  #LC 阴阳离子表合并
                        data_cat.to_csv(all_dir+'/'+f,sep='\t' ,index=False)
                    self.pos_neg_cat_diffstat_dir = all_dir
                    opts[parm_file] = all_dir
                else:
                    diff_tool = "diff_pls_" + each
                    diff_dir = globals()[diff_tool].option("id_diffStat_dir")
                    opts[parm_file] = diff_dir
                opts['metab_trans'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_desc.txt")
                opts['metab_table'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_abund.txt")
            elif parm_file in ["pca" ,"diff_pls", "sample_plsda"]:
                opts['metab_desc'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_desc.txt")
                opts['exp_file'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_abund.txt")
            elif parm_file == "anno_hmdb":
                opts['metabset'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_abund.txt")
            elif parm_file == 'sample_venn':  #add v3
                opts['metab_table'] =  os.path.join(self.preprocess.output_dir , "org_"+each+"/metab_abund.txt")  #使用原始表
            elif parm_file in  ['groups_diff', "two_sample"]:  #add v3
                opts['metab_table'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_abund.txt")
                opts['metab_desc'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_desc.txt")
            elif parm_file != "":
                new_exp_file = os.path.join(self.preprocess.option(tmp_option).prop["path"], file_name)
                opts[parm_file] = new_exp_file
            elif parm_file == "":
                opts['metab_trans'] = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_desc.txt")
            globals()[toolname].set_options(opts)
            self.logger.info(("set {} ok").format(toolname))

            if not self.start.has_key(analysis):
                if step:
                    globals()[toolname].on('start', self.set_step, {'start': step})
                self.start[analysis] = 1
            globals()[toolname].on('end', self.set_output, analysis + "_" + each)

            if analysis == "anno_hmdb":
                globals()[toolname].run()


    def run_sample_corr(self):
        self.logger.info("----start samole corr------")
        opts = {
            #'exp': exp_file,
            'scm': "complete",
            'scd': "euclidean",
            'corr_method': "pearson",
            'sct': "hierarchy",
            'transform' : 'none'

        }
        self.add_pos_neg_mix_tool(opts, self.run_type3, 'sample_corr', self.step.sample_corr, parm_file="exp")

    #add v3 20200326
    def run_sample_venn(self):
        self.logger.info('-----start sample venn ----')
        opts = {
            "threshold" : "50",
            "group_table" : self.no_qc_group_table
        }
        self.add_pos_neg_mix_tool(opts, self.run_type3, 'sample_venn', self.step.sample_venn, parm_file="sample_venn")


    #add v3  20200326
    def run_sample_plsda(self):
        self.logger.info('-----start >2 groups plasda ----')
        opts = {
            "group_file":  self.option("group_table").prop['path'] , #self.no_qc_group_table,
            "mul_type" : 'plsda',
            "confidence" : "0.95",
            "perm" : "200",
            "data_trans" : self.option('plsda_data')
        }
        self.add_pos_neg_mix_tool(opts, self.run_type3, 'sample_plsda', self.step.sample_plsda, parm_file="diff_pls")

    #add v3 20200326
    def run_groups_diff(self):
        self.logger.info('-----start  groups diff ----')
        opts = {
            "group": self.more2group_group,
            "test_method": 'ow',
            "group_name": "group",   #和 group table 第二列的列名一致
            "post_hoc": "scheffe",
            "coverage": 0.95
        }
        self.add_pos_neg_mix_tool(opts, self.run_type3, 'groups_diff', self.step.groups_diff, parm_file="groups_diff")

    def run_diff_pls(self):

        self.logger.info("-----start diff_pls----")
        group_table = self.no_qc_group_table
        #exp_file = os.path.join(self.option("pos_out"), "metab_abund.txt")
        diff_group_name =  self.vs_more_diff  #self.get_diff_group()
        opts = {
            #'exp_file': exp_file,
            "group_file": group_table,
            'group_name': diff_group_name,
            'mul_type': 'pca;plsda;oplsda',
            'confidence': '0.95;0.95;0.95',
            'perm': '0;200;200',
            #'data_trans': 'UV;Par;Par',
            #'test_method': 't-test',
            'side_type': 'two-tailed',
            'creat_metabset': True
        }
        ## 20190702 version v2
        opts['test_method'] = self.option('diff_check')
        opts['data_trans'] = ';'.join([self.option('pca_data'), self.option('plsda_data'), self.option('oplsda_data')])
        if self.filter_k != []:
            opts['filter_k'] = self.filter_k
            opts['filter_t'] = self.filter_t
            opts['filter_v'] = self.filter_v

        #self.add_pos_neg_mix_tool(opts, self.run_type1, 'diff_pls', self.step.diff_pls, parm_file="diff_pls")
        self.add_pos_neg_mix_tool(opts, self.diff_pls_types, 'diff_pls', self.step.diff_pls, parm_file="diff_pls")

    def run_two_sample(self):
        self.logger.info("-----start two_sample-----")
        opts = {
            "group_detail" : str(self.s1vs1_group_detail),
            "diff_group_name" : self.s1vs1_diff,
            "test_method" : self.option("s2_diff_check"),
            "side_type" : "two.side",
            "correct" : "bonferroni",
            "s2_p_value_fdr":self.option("s2_p_value_fdr"),
            "s2_p_fdr_condition" : self.option("s2_p_fdr_condition"),
            "s2_p_fdr_value" : self.option("s2_p_fdr_value"),
            "s2_up_condition" : self.option("s2_up_condition"),
            "s2_up_value" : self.option("s2_up_value"),
            "s2_down_condition" : self.option("s2_down_condition"),
            "s2_down_value" : self.option("s2_down_value")
        }
        self.add_pos_neg_mix_tool(opts, self.run_type3, 'two_sample', self.step.two_sample, parm_file='two_sample')

    def run_metab_venn(self):
        self.logger.info("start metab_venn")
        self.logger.info(self.creat_analysis_table.option("out_venn").prop["path"])
        if self.do_venn:
            opts = {
                "list_file": self.venn_input
            }
            self.set_step_stat(opts, self.metabset_venn, 'metabset_venn', self.step.metabset, start=True, run=False)
        else:
            pass

    def run_metab_cluster(self):
        self.logger.info("start metab_cluster")
        exp_file = self.used_out_abu
        self.logger.info(exp_file)
        opts = {
            'exp': exp_file,  ## 所有差异代谢物丰度矩阵
            'sct': 'hierarchy',
            'scd': 'euclidean',
            'scm': 'complete',
            'mct': 'hierarchy',
            'mcd': 'euclidean',
            'mcm': 'complete',
            'n_cluster':10
        }
        if self.scale_abu:
            opts["exp"] = self.scale_abu
            opts["before_scale"] = self.used_out_abu
        self.add_pos_neg_mix_tool(opts, self.run_type2, 'metabset_cluster', parm_file="")

    def run_metab_vip(self):
        self.logger.info("start metab_vip")
        metab_set_table = self.used_set_list
        opts = {
            #'diff_dir': diff_dir,
            "metab_set_table": metab_set_table,
            'vip_type': 'oplsda',
            'vip_cut': 1.0,
            'vip_top': 30,
            'mct': 'hierarchy',
            'mcd': 'euclidean',
            'mcm': 'complete',
            'scale': True
        }
        #self.add_pos_neg_mix_tool(opts, self.run_type1, 'metabset_vip', self.step.metabset, parm_file="diff_dir")
        self.add_pos_neg_mix_tool(opts, self.run_type2, 'metabset_vip',  parm_file="diff_dir")
    def run_metab_corr(self):
        self.logger.info("start metab_corr")
        exp_file = self.used_out_abu
        opts = {
            'exp': exp_file,
            'scm': "complete",
            'scd': "euclidean",
            'sct': "hierarchy",
            'corr_method': "pearson",
            'file_tran': True
        }
        self.add_pos_neg_mix_tool(opts, self.run_type2, 'metabset_corr',  parm_file="")

    def run_metab_hmdb(self):
        self.logger.info("start metab_hmdb")
        if self.option("project_type") == "LC":
            hmdb_overview = anno_hmdb_mix.option("level_out_origin").prop["path"]
        else:
            hmdb_overview = anno_hmdb_pos.option("level_out_origin").prop["path"]
        self.logger.info(hmdb_overview)
        self.hmdb_overview = hmdb_overview
        use_set = self.used_set_list
        opts = {
            'anno_overview': hmdb_overview,
            'metabset': use_set,
            "type": "metabsetanno"
        }
        self.set_step_stat(opts, self.metabset_hmdb, 'metabset_hmdb', start=True, run=False)

    def run_metab_keggc(self):
        self.logger.info("start metab_keggc")
        anno_overview = self.anno_overview
        use_set = self.used_set_list
        opts = {
            'anno_overview': anno_overview,
            'metabset': use_set
        }
        self.set_step_stat(opts, self.metabset_keggc, 'metabset_keggc', start=True, run=False)

    def run_metab_keggp(self):
        self.logger.info("start metab_keggp")
        ko_overview = self.ko_overview
        org_mul_set = self.used_diff_mul_set
        db_v = self.database_version()
        opts = {
            'anno_overview': ko_overview,
            'metabset': org_mul_set,
            "database_version": db_v
        }
        self.set_step_stat(opts, self.metabset_keggp, 'metabset_keggp', start=True, run=False)

    def run_metab_enrich(self):
        self.logger.info("start metab_enrich")
        anno_overview = self.anno_overview

        org_set = self.used_set_list
        opts = {
            'anno_overview': anno_overview,
            "ko_overview": self.ko_overview,  # add by ghd @20191016
            'metabset': org_set,
            'correct': 'BH',
            'bg': 'species',
            'method' : 'rbc'  #20190730
        }
        if self.option("organism"):
            kegg_species = self.option("organism_type") + ";" + self.option("organism")
        elif self.option("organism_type") != "All":
            kegg_species = self.option("organism_type")
        else:
            kegg_species = "All"
        opts["species"] = kegg_species
        db_v = self.database_version()
        opts['database_version'] = db_v
        self.logger.info("###w database_version: {}".format(opts))
        self.set_step_stat(opts, self.metabset_enrich, 'metabset_enrich', start=True, run=False)

    def run_metab_ipath(self):
        self.logger.info("start metab_ipath")
        anno_overview = self.anno_overview

        org_mul_set = self.used_diff_mul_set
        opts = {
            'anno_overview': anno_overview,
            'metabset': org_mul_set
        }
        self.set_step_stat(opts, self.metabset_ipath, 'metabset_ipath',  start=True, run=False)

    def run_pre_metabset(self):
        """
        代谢集前分析，主要包括注释、pca、相关性和差异分析
        """
        self.logger.info("start run annotation,sample_corr,sample_venn,sample_plada,pca,diff_pls")
        self.run_annotation()
        self.run_sample_corr()
        self.run_pca()
        self.run_sample_venn()  #add v3
        self.run_sample_plsda()    #add v3
        if self.more2group_group:
            self.run_groups_diff()  #add v3

        if self.s1vs1_diff:
            self.run_two_sample()  #add v3
        if self.vs_more_diff:
            self.run_diff_pls()
        for eachtool in self.pre_metabset_tools:
            eachtool.run()

    def run_metabset_analysis(self):
        #self.step.creat_table.finish()
        self.logger.info("start run metabset analysis...")
        self.determine_venn()
        self.determine_file()
        self.run_metab_venn()
        self.run_metab_cluster()
        if self.vs_more_diff:  #add v3
            self.run_metab_vip()
        self.run_metab_corr()
        self.run_metab_keggc()
        self.run_metab_keggp()
        self.run_metab_enrich()
        self.run_metab_ipath()
        if self.option("hmdb") == "T":
            self.run_metab_hmdb()
        for eachtool in self.metabset_tools:
            if eachtool._name == "Venn":
                if self.do_venn:
                    eachtool.run()
                else:
                    eachtool._end = True
            elif eachtool._name == "Keggp" and not self.ko:
                eachtool._end = True
            else:
                eachtool.run()
                #self.step.metabset.finish()

    def run_diff_metabset_pip(self):   # add zouguanqing v2.0
        if self.project_type == 'GC':
            diff_tool = "diff_pls_pos"
            tmp_option = 'pos_out'
            if self.vs_more_diff:
                diff_dir = globals()[diff_tool].option("id_diffStat_dir")
        else:
            if self.vs_more_diff:
                if len(self.run_type3) == 2:
                    diff_dir = self.pos_neg_cat_diffstat_dir
                else:
                    diff_tool =  "diff_pls_mix"
                    diff_dir = globals()[diff_tool].option("id_diffStat_dir")
            tmp_option = 'mix_out'
        metab_trans = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_desc.txt")
        metab_table = os.path.join(self.preprocess.option(tmp_option).prop["path"], "metab_abund.txt")

        if self.option("organism"):
            kegg_species = self.option("organism_type") + ";" + self.option("organism")
        elif self.option("organism_type") != "All":
            kegg_species = self.option("organism_type")
        else:
            kegg_species = "All"


        opts = {
            "anno_overview": self.anno_overview,
            "mul_metabset": self.creat_analysis_table.option("out_venn").prop["path"],
            "ko_overview":self.ko_overview,
            # "hmdb_overview": self.hmdb_overview,
            ###"diff_dir": diff_dir,
            "metab_trans":metab_trans,
            "metab_table": metab_table,
            "species": kegg_species
        }

        if self.hmdb_overview != '':
            opts['hmdb_overview'] = self.hmdb_overview

        if self.vs_more_diff:
            opts['diff_dir'] = diff_dir
        db_v = self.database_version()
        opts['database_version'] = db_v
        self.diff_metabset_pip.set_options(opts)
        self.diff_metabset_pip.run()


    def determine_venn(self):
        profile = self.creat_analysis_table.option("out_venn").prop["path"]
        self.logger.info(profile)
        if os.path.getsize(profile) > 0:
            table = pd.read_table(profile, sep="\t", header=None)
            if len(table) < 2:
                self.do_venn = False
            else:
                self.do_venn = True
            if len(table) > 5:
                table = pd.read_table(profile, sep="\t", header=None, nrows=5)
                tmp_venn = os.path.join(self.work_dir, "venn_input.xls")
                table.to_csv(tmp_venn, index=False, header=False, sep="\t")
                self.venn_input = tmp_venn
            else:
                self.venn_input = self.creat_analysis_table.option("out_venn")
        else:
            self.do_venn = False

    def creat_metabset_analysis_file(self):
        group_file = self.work_dir + "/" + "noQC_group.xls"
        opts = {
            #"diff_dir1": diff_pls_pos.option("metabset_dir"),
            "group_file": group_file,
            "top": self.default_top
        }
        ###20190722
        if self.option("project_type") == "LC":
            if self.vs_more_diff:
                if self.option('mix_table') == 'T':
                    opts['diff_dir1'] = diff_pls_mix.option('metabset_dir')
                else:
                    opts['diff_dir1'] = diff_pls_pos.option('metabset_dir')
                    opts['diff_dir2'] = diff_pls_neg.option('metabset_dir')
            if self.s1vs1_diff:  # add v3
                if self.option('mix_table') == 'T':
                    if globals()['two_sample_mix'].option("mix_metabset").is_set:
                        opts['two_sample_diff_dir1'] = globals()['two_sample_mix'].work_dir + '/Metabset'
                else:
                    if globals()['two_sample_pos'].option("mix_metabset").is_set:
                        opts['two_sample_diff_dir1'] = globals()['two_sample_pos'].work_dir + '/Metabset'
                    if globals()['two_sample_neg'].option("mix_metabset").is_set:
                        opts['two_sample_diff_dir2'] = globals()['two_sample_neg'].work_dir + '/Metabset'

            opts['exp'] = self.preprocess.option("mix_out").prop["metab_abun"]
        else:
            if self.vs_more_diff:
                opts['diff_dir1'] = diff_pls_pos.option('metabset_dir')
            if self.s1vs1_diff:  #add v3
                if globals()['two_sample_pos'].option("mix_metabset").is_set:
                    opts['two_sample_diff_dir1'] = globals()['two_sample_pos'].work_dir + '/Metabset'
            opts['exp'] = self.preprocess.option("pos_out").prop["metab_abun"]
        ###
        if self.option("project_type") == "LC":
            self.anno_overview = annotation_mix.option("overview_out")
            self.ko_overview = annotation_mix.option("overview_ko")
        else:
            self.anno_overview = annotation_pos.option("overview_out")
            self.ko_overview = annotation_pos.option("overview_ko")
        self.set_step_stat(opts, self.creat_analysis_table, 'creat_table', self.step.creat_table, start=True, run=True)

    def determine_file(self):
        self.logger.info("start determine_file")
        self.used_set_list = self.creat_analysis_table.option("set_list").prop["path"]
        self.used_out_abu = self.creat_analysis_table.option("out_abu").prop["path"]
        self.used_diff_mul_set = self.creat_analysis_table.option("diff_mul_set").prop["path"]
        table = pd.read_table(self.used_set_list, sep="\t", header=None)
        if len(table) < 3:
            self.used_set_list = self.preprocess.option("org_set").prop["path"]
            self.used_diff_mul_set = self.preprocess.option("org_mul_set").prop["path"]
            if self.project_type == "LC":
                self.used_out_abu = os.path.join(self.preprocess.option("mix_out").prop["path"], "metab_abund.txt")
            else:
                self.used_out_abu = os.path.join(self.preprocess.option("pos_out").prop["path"], "metab_abund.txt")
            self.scale_abu = self.scale_data(self.used_out_abu)
        else:
            self.scale_abu = self.creat_analysis_table.option("scale_abu").prop["path"]
        ko_overview = self.ko_overview.prop["path"]
        ko_table = pd.read_table(ko_overview, sep="\t", header=0)
        if len(ko_table) < 1:
            self.ko = False
        else:
            self.ko = True
        if self.option("project_type") == "LC":
            keggc_level = annotation_mix.option('keggc_level').path
        else:
            keggc_level = annotation_pos.option('keggc_level').path
        keggc_level_table = pd.read_table(keggc_level, sep="\t", header=0)
        if len(keggc_level_table) < 1:
            self.keggc_result = False
        else:
            self.keggc_result = True

    def determine_hmdb(self):
        if self.option("project_type") == "LC":
            hmdb_level = anno_hmdb_mix.option("level_out").prop["path"]
        else:
            hmdb_level = anno_hmdb_pos.option("level_out").prop["path"]
        metab_level = self.metabset_hmdb.option('level_out').path
        hmdb_level_table = pd.read_table(hmdb_level, sep="\t", header=0)
        metab_level_table = pd.read_table(metab_level, sep="\t", header=0)
        self.hmdb_level = True
        self.metab_level = True
        if len(hmdb_level_table) < 1:
            self.hmdb_level = False
        if len(metab_level_table) < 1:
            self.metab_level = False

    def scale_data(self, file_path):
        table = pd.read_table(file_path, sep="\t", index_col=0)
        scaled_data = table.apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)
        exp_profile = os.path.join(self.work_dir, "scale_data.xls")
        scaled_data.to_csv(exp_profile, index=True, header=True, sep="\t")
        return exp_profile

    def t_end(self):
        self.upload_results()
        super(MetabolomeWorkflow, self).end()

    def t_run(self):
        if os.path.exists(self.work_dir + "/FigSave"):
            os.system("rm -rf " + self.work_dir + "/FigSave")
        self.figsave = self.add_tool("metabolome.fig_save")
        self.figsave.on('end', self.t_end)
        self.figsave.set_options({
            "task_id": self.sheet.id,
            "project": "metabolome",
            "interaction": 0,
        })
        self.figsave.run()

    def get_rerun_task(self):
        task_ids = []
        list_path = os.path.abspath(os.path.dirname(__file__)) + "/rerun.txt"
        self.logger.info("list_path:" + list_path)
        if os.path.exists(list_path):
            with open(list_path, "r") as r:
                for l in r:
                    task_ids.extend(l.strip('\n').split())
        return task_ids

    def run(self):
        """
        运行 metabolome workflow
        :return:
        """
        if self.sheet.id in self.get_rerun_task() and self.option("save_pdf"):
            self.logger.info("rerun figsave")
            self.t_run()
            super(MetabolomeWorkflow, self).run()
            return
        task_info = self.api.api('task_info.metabolome_task_info')
        sg_task_diff_params = {
            'p_fdr': self.option('p_fdr'),
            'p_fdr_condition': self.option('p_fdr_condition'),
            'p_fdr_value': self.option('p_fdr_value'),
            'vip': self.option('vip'),
            'vip_condition': self.option('vip_condition'),
            'vip_value': self.option('vip_value'),
            'up_condition': self.option('up_condition'),
            'up_value': self.option('up_value'),
            'down_condition': self.option('down_condition'),
            'down_value': self.option('down_value')
        }

        sg_task_two_sample_diff_params = {
            's2_p_value_fdr' : self.option('s2_p_value_fdr'),
            's2_p_fdr_condition' : self.option('s2_p_fdr_condition'),
            's2_p_fdr_value' :  self.option("s2_p_fdr_value"),
            's2_up_condition' : self.option('s2_up_condition'),
            's2_up_value' : self.option('s2_up_value'),
            's2_down_condition' : self.option('s2_down_condition'),
            's2_down_value' : self.option('s2_down_value')
        }
        task_info.add_task_info(diff_params=sg_task_diff_params, two_sample_diff_params=sg_task_two_sample_diff_params)
        ## 下面对数据库的版本进行更新，便于进行新老版本兼容，老版本的数据库不会存在database字段
        if self.option("kegg_version").lower() in ['KEGG_V94.2', 'KEGG_V94', "kegg_v94", "kegg_v94.2"]:
            database_type = json.dumps({
                "kegg": "v94.2",
            },sort_keys=True, separators=(',', ':'))
        elif self.option("kegg_version").lower() in ["kegg_v2021.09.18"]:
            database_type = json.dumps({
                "kegg": "v2021.09.18",
            },sort_keys=True, separators=(',', ':'))
        else:
            self.logger.info("本次分析用老的KEGG数据库进行分析，数据库版本未知，请知悉！")
        task_info.update_mongo('sg_task',{"task_id":self.task_id}, {"database":database_type})


        k_list = []
        t_list = []
        v_list = []

        if self.option('p_fdr'):
            if self.option('p_fdr') in ['P_value']:
                k_list.append('P_value')
            else:
                k_list.append('fdr')
            v_list.append(str(self.option('p_fdr_value')))

            if self.option('p_fdr_condition') in ['>','》','$gt']:
                t_list.append('gt')
            elif self.option('p_fdr_condition') in ['<','《','$lt']:
                t_list.append('lt')
            elif self.option('p_fdr_condition') in ['=']:
                t_list.append('eq')
            elif self.option('p_fdr_condition') in ['≤','$lte']:
                t_list.append('lte')
            elif self.option('p_fdr_condition') in ['≥','$gte']:
                t_list.append('gte')
            else:
                raise Exception('WRONG: %s' %self.option('p_fdr_condition'))

        if self.option('vip'):
            if self.option('vip') in ['VIP_PLS-DA']:
                k_list.append('Vip_plsda')
            else:
                k_list.append('Vip_oplsda')
            v_list.append(str(self.option('vip_value')))

            if self.option('vip_condition') in ['>','》','$gt']:
                t_list.append('gt')
            elif self.option('vip_condition') in ['<','《','$lt']:
                t_list.append('lt')
            elif self.option('vip_condition') in ['=']:
                t_list.append('eq')
            elif self.option('vip_condition') in ['≤','$lte']:
                t_list.append('lte')
            elif self.option('vip_condition') in ['≥','$gte']:
                t_list.append('gte')
            else:
                raise Exception('WRONG: %s' %self.option('vip_condition'))


        if self.option('up_condition'):
            k_list.append('FC_up')
            if self.option('up_condition') in ['>','》','$gt']:
                t_list.append('gt')
            else:
                t_list.append('gte')
            v_list.append(str(self.option('up_value')))

        if self.option('down_condition'):
            k_list.append('FC_down')
            if self.option('down_condition') in ['<','《','＜','$lt']:
                t_list.append('lt')
            else:
                t_list.append('lte')
            v_list.append(str(self.option('down_value')))

        self.filter_k = ','.join(k_list)
        self.filter_t = ','.join(t_list)
        self.filter_v = ','.join(v_list)
        self.logger.info('print metabset filter:')
        self.logger.info(self.filter_k)
        self.logger.info(self.filter_t)
        self.logger.info(self.filter_v)


        hmdb = self.option("hmdb")
        if self.project_type == 'GC':
            task_info.add_project_type(self.task_id, self._sheet.project_sn, self.project_type, hmdb=hmdb)
        else:
            task_info.add_project_type(self.task_id, self._sheet.project_sn, self.project_type, hmdb=hmdb,mix_table=self.is_mix_table) #20190717

        if self.option('test'):
            '''
            self.creat_analysis_table.on("end", self.run_metabset_analysis)
            self.on_rely(self.metabset_tools, self.run_api)
            self.creat_metabset_analysis_file()
            '''
            self.run_api()
        else:
            self.preprocess.on("end", self.run_pre_metabset)
            if self.option("hmdb") == "T":
                self.on_rely(self.pre_metabset_tools, self.run_anno_hmdb)
                if self.option("project_type") == "LC":
                    anno_hmdb_mix.on("end", self.creat_metabset_analysis_file)
                else:
                    anno_hmdb_pos.on("end", self.creat_metabset_analysis_file)
            else:
                self.on_rely(self.pre_metabset_tools, self.creat_metabset_analysis_file)
            self.creat_analysis_table.on("end", self.run_metabset_analysis)

            self.on_rely(self.metabset_tools, self.run_diff_metabset_pip)
            self.diff_metabset_pip.on('end',self.run_api)  ## v2.0 添加

            self.check_group()
            self.get_noQC_group()   #获取 self.no_qc_group_table
            self.check_specimen_names()
            self.run_preprocess()

        super(MetabolomeWorkflow, self).run()

    def run_api(self, test=False):  # 原run_api_and_set_outpu
        self.logger.info("开始导表")

        if self.option("test"):
            pass
        if self.option("hmdb") == "T":
            self.determine_hmdb()
        self.logger.info("开始样本，分组和代谢集导表")
        self.export_group_and_samples()
        self.set_params(set_analysis="preprocess")
        self.export_preprocess()
        self.group_detail_noQC = deepcopy(self.group_detail)
        if self.group_detail_noQC.has_key('QC'):
            self.group_detail_noQC.pop('QC')
        self.export_org_metabset()
        self.export_diff_metabset()
        self.logger.info("开始代谢集前导表")
        self.set_params(set_analysis="pre_metabset")
        self.export_pca()
        self.export_sample_corr()
        self.export_sample_venn() #add v3
        self.export_sample_plsda()  #add v3
        if self.more2group_group:
            self.export_groups_diff()  #add v3
        if self.vs_more_diff:  #add v3
            self.export_diff_pls()
        if self.s1vs1_diff:  #add v3
            self.export_two_sample()
        if self.keggc_result:
            self.export_anno_keggc()
        if self.ko:
            self.export_anno_keggp()
        if self.option("hmdb") == "T" and self.hmdb_level:
            self.export_anno_hmdb()
        self.export_anno_overview()
        self.logger.info("开始代谢集分析导表")
        self.set_params(set_analysis="metabset")
        if self.do_venn:
            self.export_metab_venn()
        self.export_metab_cluster()
        if self.vs_more_diff:
            self.export_metab_vip()
        self.export_metab_corr()
        self.export_metab_keggc()
        if self.ko:
            self.export_metab_keggp()
        self.export_metab_enrich()
        self.export_metab_ipath()
        #self.export_metab_roc()
        if self.option("hmdb") == "T" and self.metab_level:
            self.export_metab_hmdb()

        self.export_diff_metabset_pip_result() #zouguanqing  导多个代谢集分析的pip分析结果

        #工作流有代谢集参与的分析，记录到sg_status中，在代谢集管理页面中会用到  #add v3
        task_info = self.api.api('task_info.metabolome_task_info')
        updata_collect_names = ['metabset_vip','metabset_venn','metabset_cluster','metabset_kegg_enrich',
                                'metabset_corr','metabset_keggp','metabset_keggc','metabset_hmdb','metabset_ipath']

        task_info.workflow_metabset_add_sg_status(updata_collect_names, self.task_id,is_workflow_use=True)


        project_stat = self.api.api("metabolome.project_statistics")
        if self.option("project_type") == "LC":
            origin_exp_file = os.path.join(self.preprocess.option("mix_out").prop["path"],'metab_abund.txt')
        else:
            origin_exp_file = os.path.join(self.preprocess.option("pos_out").prop["path"],'metab_abund.txt')

        if len(self.run_type3) == 2:
            if self.s1vs1_diff:
                project_stat.project_check_pip(self.output_dir, self.option('group_table').path, self.creat_analysis_table.output_dir, None,None,abund=origin_exp_file)
            else:
                pls_dir = self.diff_pls_tool.option("pls_dir").prop["path"]
                pls_dir_neg = self.diff_pls_tool_neg.option("pls_dir").prop["path"]
                project_stat.project_check_pip(self.output_dir, self.option('group_table').path, self.creat_analysis_table.output_dir, pls_dir,pls_dir_neg,abund=origin_exp_file)
        else:
            if self.s1vs1_diff:
                project_stat.project_check_pip(self.output_dir, self.option('group_table').path, self.creat_analysis_table.output_dir, None,abund=origin_exp_file)
            else:
                if self.run_type3[0] == 'pos':
                    pls_dir = self.diff_pls_tool.option("pls_dir").prop["path"]
                else:
                    pls_dir = self.diff_pls_tool.option("pls_dir").prop["path"]

                project_stat.project_check_pip(self.output_dir, self.option('group_table').path, self.creat_analysis_table.output_dir, pls_dir,abund=origin_exp_file)
        task_info.update_mongo('sg_task',{"task_id":self.task_id}, {'version':'3.5'})

        if self.option("save_pdf"):
            task_info.update_mongo('sg_task', {"task_id": self.task_id}, {"save_pdf": 1})
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": self.sheet.id,
                "project": "metabolome",
                "interaction": 0,
            })
            self.figsave.run()
        else:
            self.end()

    def set_output(self, event):
        """
        将各个模块的结果输出至output
        """
        self.logger.info("set_output...")
        obj = event['bind_object']
        self.logger.info(event)
        o_dirs = {
            "preprocess": 'Preprocess',
            "sample_corr_pos": "ExpCorr/pos",
            "sample_corr_neg": "ExpCorr/neg",
            "sample_corr_mix": "ExpCorr/mix", ###20190717
            "pca_pos": "ExpPCA/pos",
            "pca_neg": "ExpPCA/neg",
            "pca_mix": "ExpPCA/mix",   ###20190717
            "diff_pls_pos": "ExpDiff/pos",
            "diff_pls_neg": "ExpDiff/neg",
            "diff_pls_mix": "ExpDiff/mix",  ###20190717
            "metabset_venn": "MetabsetVenn/",
            "metabset_roc": "MetabsetRoc/",
            "metabset_keggc": "MetabsetKeggc/",
            "metabset_keggp": "MetabsetKeggp/",
            "metabset_enrich": "MetabsetEnrich/",
            "metabset_ipath": "MetabsetIpath/",
            "metabset_cluster_mix": "MetabsetCluster/",
            "metabset_vip_pos": "MetabsetVip/pos",
            "metabset_vip_neg": "MetabsetVip/neg",
            "metabset_vip_mix": "MetabsetVip/mix",    ###20190717
            "metabset_corr_mix": "MetabsetCorr/",
            "anno_hmdb_mix": "AnnoHmdb/",
            "metabset_hmdb": "MetabsetHmdb/",
            #### GC用
            "metabset_corr_pos": "MetabsetCorr/",
            "metabset_cluster_pos": "MetabsetCluster/",
            "anno_hmdb_pos": "AnnoHmdb/",
            ##add v3
            "sample_venn_pos": "ExpVenn/pos",
            "sample_venn_neg": "ExpVenn/neg",
            "sample_venn_mix": "ExpVenn/mix",
            "sample_plsda_pos" : "ExpPLSDA/pos",
            "sample_plsda_neg" : "ExpPLSDA/neg",
            "sample_plsda_mix" : "ExpPLSDA/mix",
            "two_sample_pos" : "TwoSample/pos",
            "two_sample_neg" : "TwoSample/neg",
            "two_sample_mix" : "TwoSample/mix",
            "groups_diff_pos" : "GroupsDiff/pos",
            "groups_diff_neg" : "GroupsDiff/neg",
            "groups_diff_mix" : "GroupsDiff/mix",
        }
        #list1 = ['preprocess', "sample_corr_pos", "sample_corr_neg", "pca_pos", "pca_neg", "metabset_venn",
        #         "metabset_roc", "metabset_keggc", "metabset_keggp", "metabset_enrich", "metabset_ipath",
        #         "metabset_cluster", "metabset_vip"]
        list2 = ["annotation_mix", "annotation_pos", 'diff_pls_pos', 'diff_pls_neg','diff_pls_mix', 'creat_table']  #add diff_pls_mix
        if event['data'] not in list2:
            self.move_dir(obj.output_dir, o_dirs[event['data']], event['data'])
        if event['data'] in ['diff_pls_pos', 'diff_pls_neg','diff_pls_mix']:  #20190717 add diff_pls_mix
            self.move_pls(obj.work_dir, o_dirs[event['data']])
        if event['data'] == 'annotation_mix' or event['data'] == 'annotation_pos':
            self.move_anno(obj.work_dir)

    def move_dir(self, olddir, newname, event_data):  # 原函数名move2outputdir
        """
        移动一个目录下所有文件/文件夹到workflow输出路径下，供set_output调用
        """
        start = time.time()
        if not os.path.isdir(olddir):
            #raise Exception('需要移动到output目录的文件夹不存在。', code="14700101")
            self.set_error('需要移动到output目录的文件夹不存在。', code="14700101")
        newdir = os.path.join(self.output_dir, newname)
        self.logger.info("newdir is : " + newdir)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = []
        newfiles = []
        if "anno_hmdb" in event_data:
            allfiles.remove("anno.xls")
        for i in allfiles:
            oldfile = os.path.join(olddir, i)
            oldfiles.append(oldfile)
            newfilename = self.change_name(i, event_data)
            newfile = os.path.join(newdir, newfilename)
            newfiles.append(newfile)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(oldfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}移动到{},耗时{}s".format(olddir, newdir, duration))


    def change_name(self, filename, event_data):
        if "sample_corr" in event_data:
            prename = "sample"
        if "metabset_corr" in event_data:
            prename = "metab"
        if filename == "corr.cluster_tree.xls":
            newname = prename + "_corr_tree.xls"
        elif filename == "corr.xls":
            newname = prename + "_corr.xls"
        elif filename == "pvalue.xls":
            newname = "corr_pvalue.xls"
        else:
            newname = filename
        return newname

    def move_anno(self, olddir):
        start = time.time()
        if not os.path.isdir(olddir):
            #raise Exception('需要移动到output目录的文件夹不存在。', code="14700102")
            self.set_error('需要移动到output目录的文件夹不存在。', code="14700102")
        all_dirs = os.listdir(olddir)
        oldfiles = []
        newfiles = []
        for eachdir in all_dirs:
            if eachdir in ["AnnoKeggc", "AnnoKeggp", "AnnoOverview"]:
                newdir = os.path.join(self.output_dir, eachdir)
                self.logger.info("newdir is : " + newdir)
                if not os.path.exists(newdir):
                    os.makedirs(newdir)
                eachdirpath = olddir + "/" + eachdir + "/output"
                allfiles = os.listdir(eachdirpath)
                for i in allfiles:
                    oldfile = os.path.join(eachdirpath, i)
                    oldfiles.append(oldfile)
                    newfile = os.path.join(newdir, i)
                    newfiles.append(newfile)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(oldfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("注释文件夹{}已链接,耗时{}s".format(olddir, duration))

    def move_pls(self, olddir, newdir):
        start = time.time()
        if not os.path.isdir(olddir):
            #raise Exception('需要移动到output目录的文件夹不存在。', code="14700103")
            self.set_error('需要移动到output目录的文件夹不存在。', code="14700103")
        all_dirs = os.listdir(olddir)
        oldfiles = []
        newfiles = []
        newdir = os.path.join(self.output_dir, newdir)
        # for eachdir in all_dirs:
        #     if eachdir == "DiffMulStat":
        #         eachdirpath = olddir + "/DiffMulStat/output"
        #         allfiles = os.listdir(eachdirpath)
        #         newdir_pls = os.path.join(newdir, "DiffMulStat")
        #         for i in allfiles:
        #             oldfile = os.path.join(eachdirpath, i)
        #             oldfiles.append(oldfile)
        #             newfile = os.path.join(newdir_pls, i)
        #             newfiles.append(newfile)
        #     if eachdir == "MergeDiff":
        #         DiffTestdir = olddir + "/MergeDiff/output/DiffStat"
        if "DiffMulStat" in all_dirs:
            eachdirpath = olddir + "/DiffMulStat/output"
            allfiles = os.listdir(eachdirpath)
            newdir_pls = os.path.join(newdir, "DiffMulStat")
            for i in allfiles:
                oldfile = os.path.join(eachdirpath, i)
                oldfiles.append(oldfile)
                newfile = os.path.join(newdir_pls, i)
                newfiles.append(newfile)

        for eachdir in ['MergeDiff__2','MergeDiff__1','MergeDiff']:
            if eachdir in all_dirs:
                DiffTestdir = olddir + "/"+eachdir+"/output/DiffStat"
                if os.path.exists(DiffTestdir):
                    allfiles = os.listdir(DiffTestdir)
                    newdir_test = os.path.join(newdir, "DiffTest")
                    for i in allfiles:
                        oldfile = os.path.join(DiffTestdir, i)
                        oldfiles.append(oldfile)
                        newfile = os.path.join(newdir_test, i)
                        newfiles.append(newfile)
                break


        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(oldfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}已链接,耗时{}s".format(olddir, duration))

    def move_file(self, old_file, new_file):
        """
        递归移动文件夹的内容，供move_dir调用
        """
        if os.path.isfile(old_file):
            if not os.path.isdir(os.path.dirname(new_file)):
                os.makedirs(os.path.dirname(new_file))
            old_file_name = old_file.split("/")[-1]
            if not old_file_name in ["PCA.ellipse.xls", "OPLS-DA.ellipse.xls", "OPLS-DA.intercept.xls", "ko.xls",
                                     "OPLS-DA.loading.xls", "OPLS-DA.vip.xls", "PCA.loading.xls", "PLS-DA.ellipse.xls",
                                     "PLS-DA.intercept.xls", "PLS-DA.loading.xls", "PLS-DA.vip.xls", "metabset.list",
                                     "mul_metabset.list", "gene_ipath_input.xls"]:
                os.link(old_file, new_file)
        elif os.path.isdir(old_file):
            os.makedirs(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)
        else:
            self.logger.info("导出失败：请检查{}".format(old_file))

    def end(self):
        self.up_files()
        super(MetabolomeWorkflow, self).end()

    def table_len(self, file):
        with open(file, "r") as f1:
            lines = f1.readlines()
        table_length = len(lines)
        return table_length

    def set_params(self, set_analysis=None):
        ########### params 导表有的内部已加 json.dumps
        self.dic_params = {}
        #self.metabset_table_type = "mix"
        self.metabset_table_type = (self.run_type2)[0]
        if set_analysis == "preprocess":
            if self.option("rsd") == "none":
                rsd = ""
            else:
                rsd = self.option("rsd")
            self.logger.info("set preprocess,pca,anno,diff_pls,sample_corr params")
            para_preprocess = {
                "fillna": self.option("fillna"),
                "norm": self.option("norm"),
                "rsd": rsd,
                "scale": self.option("scale"),
                'task_id': self.task_id,
                'submit_location': "exp",
                'task_type': 2,
                'fill_type': self.option('fill_type')  #v3 202003
            }
            if self.option("norm") == "sample":
                para_preprocess["sample_name"] = self.option("sample_name")
            if self.option("norm") == "inner":
                para_preprocess["inner_ref"] = self.option("inner_ref")
            if self.option("scale") == 'defined':
                para_preprocess['log'] = self.option("log")
            ## add 20190722  rm_nan
            para_preprocess['fill_percent'] = str(self.option("fill_percent"))

            self.dic_params["preprocess"] = para_preprocess
            para_preprocess_raw = {
                "fillna": "none",
                "norm": "none",
                "rsd": "",
                "fill_percent": "100",
                "scale": "none",
                'task_id': self.task_id,
                'submit_location': "exp",
                'task_type': 2,
            }

            self.dic_params["preprocess_raw"] = para_preprocess_raw
        elif set_analysis == "pre_metabset":
            para_anno_keggc = {
                'metab_table': str(self.metab_table_id),
                ###'table_type': self.metabset_table_type,  v2.0
                'task_id': self.task_id,
                'submit_location': "annokeggc",
                'task_type': 2,
                'database':'CBR'  #2019 v2
            }

            self.dic_params["anno_keggc"] = para_anno_keggc
            if self.option("organism"):
                kegg_species = self.option("organism_type") + ";" + self.option("organism")
            elif self.option("organism_type") != "All":
                kegg_species = self.option("organism_type")
            else:
                kegg_species = "All"
            para_anno_keggp = {
                'metab_table': str(self.metab_table_id),
                ###'table_type': self.metabset_table_type,   v2.0
                'task_id': self.task_id,
                'organism': kegg_species,
                'submit_location': "annokeggp",
                'task_type': 2,

            }
            self.dic_params["anno_keggp"] = para_anno_keggp
            para_anno_overview = {
                'set_id': str(self.org_set_id),
                'task_id': self.task_id,
                'organism': kegg_species,
                'submit_location': "anno_overview",
                #'task_type': 1
            }
            self.dic_params["anno_overview"] = para_anno_overview
            para_sample_corr = {
                "metab_table": str(self.metab_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.group_detail,
                "dist": "euclidean",
                "coefficient": "pearson",
                "cluster": "complete",
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "expcorr",
                "transform" : 'None'    #None
            }
            self.dic_params["sample_corr"] = para_sample_corr
            para_pca = {
                "metab_table": str(self.metab_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.group_detail,
                "transform": "UV",
                "confidence": 0.95,
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "exppca",
            }
            self.dic_params["pca"] = para_pca

            #add v3 params
            para_sample_venn = {
                "metab_table": str(self.raw_main_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.group_detail_noQC,
                'threshold' : "50",
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "expvenn",
            }
            self.dic_params["sample_venn"] = para_sample_venn

            para_sample_plsda = {
                "metab_table": str(self.metab_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.group_detail,  #self.group_detail_noQC,
                "confidence": 0.95,
                "replace" : 200,
                "trans" : self.option('plsda_data'),
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "expplsda",
            }
            self.dic_params["sample_plsda"] = para_sample_plsda

            #add v3 params
            para_groups_diff = {
                "metab_table": str(self.metab_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.more2_group_detail,
                "post_hoc": "scheffe",
                "coverage" : "0.95",
                "test_method" : "ow",
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "diff_multigroup",
            }
            self.dic_params["groups_diff"] = para_groups_diff

            #add v3 2020
            two_sample = {
                "metab_table" : str(self.metab_table_id),
                "group_id" :str(self.group_id),
                "group_detail": self.s1vs1_group_detail ,
                "diff_group_id": str(self.diff_group_id),
                "diff_group_detail" : sorted(self.s1vs1_diff.split(';')),
                "test_method" : self.option("s2_diff_check"),
                "correct" : "bonferroni",
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "diff_twosample"  ##
            }
            if self.option("s2_diff_check") == 'fisher':
                two_sample['tail'] = 'two-tailed'
            self.dic_params["two_sample"] = two_sample

            para_diff_pls = {
                "metab_table": str(self.metab_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.vs_more_group_detail,  #self.group_detail_noQC
                "diff_group_id": str(self.diff_group_id),
                "diff_group_detail": sorted(self.vs_more_diff.split(";")),   # ori self.diff_detail
                "test_method": self.option('diff_check'),
                "tail": "two-tailed",
                "pca_trans": self.option("pca_data"),
                "pca_confidence": "0.95",
                "plsda_trans": self.option('plsda_data'),
                "plsda_confidence": "0.95",
                "plsda_replace": "200",
                "oplsda_trans": self.option('oplsda_data'),
                "oplsda_confidence": "0.95",
                "oplsda_replace": "200",
                "task_id": self.task_id,
                "task_type": 2,
                "submit_location": "diff_twogroup"  #v2 v3
            }
            self.dic_params["diff_pls"] = para_diff_pls
            para_anno_hmdb = {
                'metab_table': str(self.metab_table_id),
                # 'table_type': self.metabset_table_type,  # version v2 modified by ghd @191024
                'task_id': self.task_id,
                'submit_location': "annohmdb",
                'task_type': 2
            }
            self.dic_params["anno_hmdb"] = para_anno_hmdb
            ##### 代谢集params
        else:
            para_metab_venn = {
                'metabset': ",".join(self.diff_group_set_id[0:5]),
                "task_id": self.task_id,
                "task_type": 2,
                "submit_location": "metabsetvenn"
            }
            self.dic_params["metab_venn"] = para_metab_venn
            para_metab_cluster = {
                "metab_table": str(self.metab_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.group_detail_noQC,
                "metab_cluster_method": "hierarchy",
                "metab_dist": "euclidean",
                "metab_cluster": "complete",
                "sam_cluster_method": "hierarchy",
                "sam_dist": "euclidean",
                "sam_cluster": "complete",
                "scale": "scale",
                "top_meta": self.top_diff,
                "metab_set": str(self.diff_set_id),
                "group_method": "no",
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "metabsetcluster",
                "n_cluster" : 10  #20190729
            }
            self.dic_params["metab_cluster"] = para_metab_cluster
            para_metab_vip = {
                "metab_table": str(self.metab_table_id),  #preprocess id
                "group_id": str(self.group_id),
                "group_detail": self.group_detail_noQC,
                "group_method": "none",
                "top_vip": "30",
                "metab_cluster": "complete",
                "metab_cluster_method": "hierarchy",
                "metab_dist": "euclidean",
                "diff_group_id": str(self.diff_group_id),
                "diff_group_detail": self.diff_detail,
                "metab_set": str(self.diff_set_id),
                "vip": "1",
                "scale": "scale",
                "vip_from": "oplsda",
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "metabsetvip",
            }
            self.dic_params["metab_vip"] = para_metab_vip
            para_metab_corr = {
                "metab_table": str(self.metab_table_id),
                "group_id": str(self.group_id),
                "group_detail": self.group_detail_noQC,
                "coefficient": "pearson",
                "dist": "euclidean",
                "cluster_type": "hierarchy",
                "cluster": "complete",
                "top": self.top_diff,
                "metab_set": str(self.diff_set_id),
                'task_id': self.task_id,
                "task_type": 2,
                "submit_location": "metabsetcorr",
            }
            self.dic_params["metab_corr"] = para_metab_corr
            para_metab_keggc = {
                'metabset': str(self.diff_set_id),
                'task_type': 2,
                'task_id': self.task_id,
                'submit_location': "metabsetkeggc",
                'compound':'CBR'  #2019 v2
            }
            self.dic_params["metab_keggc"] = para_metab_keggc
            para_metab_hmdb = {
                'metabset': str(self.diff_set_id),
                'task_type': 2,
                'task_id': self.task_id,
                'submit_location': "metabsethmdb"
            }
            self.dic_params["metab_hmdb"] = para_metab_hmdb
            para_metab_keggp = {
                'metabset': str(self.diff_set_id),
                'task_type': 2,
                'task_id': self.task_id,
                'submit_location': "metabsetkeggp"
            }
            self.dic_params["metab_keggp"] = para_metab_keggp
            para_metab_enrich = {
                'metabset': str(self.diff_set_id),
                'correct': "BH",
                'bg': "species",
                'task_type': 2,
                'task_id': self.task_id,
                'submit_location': "metabsetkeggenrich",
                'topo' : 'rc'  #20190730
            }
            self.dic_params["metab_enrich"] = para_metab_enrich
            para_metab_ipath = {
                'metabset': str(self.diff_set_id),
                'task_type': 2,
                'task_id': self.task_id,
                'submit_location': "metabsetipath"
            }
            self.dic_params["metab_ipath"] = para_metab_ipath

    @time_count
    def export_group_and_samples(self):
        self.api_sam_group = self.api.api('metabolome.specimen_group')
        file_path = self.group_table.prop["path"]
        rm_repeat = self.api.api('metabolome.repeat')
        rm_repeat.rm_each("specimen_group", name="group_name", name_detail="GroupOrigin")
        group_info = self.api_sam_group.add_group_table(file_path, task_id=self.task_id, is_used="0",
                                                        origin_name="GroupOrigin")
        self.api_sam_group.add_sample_table(file_path)
        diff_file_path = self.option("diff_group").prop["path"]
        group_info = group_info[0]
        self.group_id = group_info["group_id"]
        self.group_detail = group_info["group_detail"]
        # modify 解决交互分析重运行问题
        items = self.group_detail.items()
        items = sorted(items)
        self.group_detail = {}
        for key, value  in items:
            value = sorted(value)
            self.group_detail[key] = value
        rm_repeat.rm_each("specimen_group", name="compare_group_name", name_detail="DiffOrigin")
        self.diff_group_id, self.diff_detail = self.api_sam_group.add_compare_group(diff_file_path, self.group_id,
                                                                                    diff_group_name="DiffOrigin")
        self.diff_detail = sorted(self.diff_detail)

    @time_count
    def export_preprocess(self):
        # rere_path0 = self._sheet.output.split(':')[1]
        self.api_dic["preprocess"] = self.api.api('metabolome.preprocess')
        # rere_path = rere_path0 + "/Preprocess/org_pos"
        rere_path = os.path.join(self._sheet.output, self.link_up_map["Preprocess"]+"/org_pos/")
        table_path = self.preprocess.option("org_pos_out").prop["path"]
        project_type = self.option('project_type')
        if project_type == "LC":
            rere_path += "," + os.path.join(self._sheet.output, self.link_up_map["Preprocess"]+"/org_neg/")
            rere_path += "," + os.path.join(self._sheet.output, self.link_up_map["Preprocess"]+"/org_mix/")
            table_path += "," + self.preprocess.option("org_neg_out").prop["path"]
        name = "raw"
        params_raw = self.dic_params["preprocess_raw"]
        raw_main_table_id = self.api_dic["preprocess"].add_metab_table(name, project_type, rere_path, is_raw=1,
                                                                   params=params_raw)
        self.raw_main_table_id = raw_main_table_id
        self.api_dic["preprocess"].add_metab_table_detail(raw_main_table_id, table_path, raw_main_table_id)
        name2 = "MetabTable_Origin"
        # rere_path2 = rere_path0 + "/Preprocess/pos"
        rere_path2 = os.path.join(self._sheet.output, self.link_up_map["Preprocess"]+"/pos/")
        table_path2 = self.preprocess.option("pos_out").prop["path"]
        if project_type == "LC":    ###LC preprocess 始终分3种情况
            rere_path2 += "," + os.path.join(self._sheet.output, self.link_up_map["Preprocess"]+"/neg/")
            rere_path2 += ',' + os.path.join(self._sheet.output, self.link_up_map["Preprocess"]+"/mix/")
            table_path2 += "," + self.preprocess.option("neg_out").prop["path"]
            table_path2 += ',' + self.preprocess.option("mix_out").prop["path"]
        params = self.dic_params["preprocess"]
        self.logger.info(table_path2)
        main_table_id2 = self.api_dic["preprocess"].add_metab_table(name2, project_type, rere_path2, is_raw=0,
                                                                    params=params)
        self.api_dic["preprocess"].add_metab_table_detail(main_table_id2, table_path2, raw_main_table_id)
        self.logger.info(main_table_id2)
        self.metab_table_id = main_table_id2
        self.objid2dir[str(raw_main_table_id)] = self.link_up_map["Preprocess"]
        self.objid2dir[str(main_table_id2)] = self.link_up_map["Preprocess"]
        return main_table_id2

    @time_count
    def export_anno_keggc(self):
        self.api_dic["anno_keggc"] = self.api.api('metabolome.anno_keggc')
        name = "AnnoKeggc_Origin"
        params = self.dic_params["anno_keggc"]
        main_table_id = self.api_dic["anno_keggc"].add_anno_keggc_main(name, params=params)
        if self.option("project_type") == "LC":
            level_out = annotation_mix.option('keggc_level').path
            stat_out = annotation_mix.option('keggc_stat').path
        else:
            level_out = annotation_pos.option('keggc_level').path
            stat_out = annotation_pos.option('keggc_stat').path
        self.logger.info("stat_out_file:{}".format(stat_out))
        self.api_dic["anno_keggc"].add_level_detail(main_table_id, level_out)
        self.api_dic["anno_keggc"].add_stat_detail(main_table_id, stat_out)
        self.objid2dir[str(main_table_id)] = self.link_up_map["AnnoKeggc"]

    @time_count
    def export_anno_hmdb(self):
        self.api_dic["anno_hmdb"] = self.api.api('metabolome.anno_hmdb')
        name = "AnnoHmdb_Origin"
        params = self.dic_params["anno_hmdb"]
        if self.option("project_type") == "LC":
            #hmdb_overview_origin = anno_hmdb_mix.option("level_out_origin").prop["path"]
            hmdb_overview_origin = os.path.join(self._sheet.output, self.link_up_map["AnnoHmdb"]+"/HmdbLevel_Origin.xls")
            level_out = anno_hmdb_mix.option('level_out').path
            out_dir = anno_hmdb_mix.output_dir
        else:
            #hmdb_overview_origin = anno_hmdb_pos.option("level_out_origin").prop["path"]
            hmdb_overview_origin = os.path.join(self._sheet.output, self.link_up_map["AnnoHmdb"]+"/HmdbLevel_Origin.xls")
            level_out = anno_hmdb_pos.option('level_out').path
            out_dir = anno_hmdb_pos.output_dir
        main_table_id = self.api_dic["anno_hmdb"].add_anno_hmdb_main(name, params=params,
                                                                     anno_hmdb=hmdb_overview_origin)
        self.api_dic["anno_hmdb"].add_level_detail(main_table_id, level_out, database="annohmdb")
        self.api_dic["anno_hmdb"].add_stat_detail(main_table_id, out_dir, database="annohmdb")
        self.objid2dir[str(main_table_id)] = self.link_up_map["AnnoHmdb"]

    @time_count
    def export_anno_keggp(self):
        self.api_dic["anno_keggp"] = self.api.api('metabolome.anno_keggp')
        name = "AnnoKeggp_Origin"
        params = self.dic_params["anno_keggp"]
        if self.option("project_type") == "LC":
            level_out = annotation_mix.option('keggp_level').path
            stat_out = annotation_mix.option('keggp_stat').path
        else:
            level_out = annotation_pos.option('keggp_level').path
            stat_out = annotation_pos.option('keggp_stat').path
        # result_dir = self.remote_dir + "/AnnoKeggp/pathway_img/"
        result_dir = os.path.join(self._sheet.output, self.link_up_map["AnnoKeggp"]+"/pathway_img/")
        # final_result_dir = "rerewrweset" + result_dir.split("rerewrweset")[1]
        main_table_id = self.api_dic["anno_keggp"].add_anno_keggp_main(name, result_dir, params=params)
        self.api_dic["anno_keggp"].add_level_detail(main_table_id, level_out)
        self.api_dic["anno_keggp"].add_stat_detail(main_table_id, stat_out)
        self.objid2dir[str(main_table_id)] = self.link_up_map["AnnoKeggp"]

    @time_count
    def export_anno_overview(self):
        if self.option("project_type") == "LC":
            final_anno_overview = os.path.join(self.output_dir, "AnnoOverview/anno.xls")
            if self.option("hmdb") == "T":
                new_anno_overview = anno_hmdb_mix.option("new_overview").prop["path"]
        else:
            final_anno_overview = os.path.join(self.output_dir, "AnnoOverview/anno.xls")
            if self.option("hmdb") == "T":
                new_anno_overview = anno_hmdb_pos.option("new_overview").prop["path"]
        if self.option("hmdb") == "T":
            if os.path.exists(final_anno_overview):
                os.remove(final_anno_overview)
            os.link(new_anno_overview, final_anno_overview)
        self.api_dic["anno_overview"] = self.api.api('metabolome.anno_overview')
        name = "Overview_Origin"
        params = self.dic_params["anno_overview"]
        main_table_id = self.api_dic["anno_overview"].add_overview(params, name=name)
        if self.option("project_type") == "LC":
            #overview_table = annotation_mix.option("overview_out").prop["path"]
            ko_table = annotation_mix.option("overview_ko").prop["path"]
            origin_exp_file = os.path.join(self.preprocess.option("mix_out").prop["path"],'metab_abund.txt')
        else:
            #overview_table = annotation_pos.option("overview_out").prop["path"]
            ko_table = annotation_pos.option("overview_ko").prop["path"]
            origin_exp_file = os.path.join(self.preprocess.option("pos_out").prop["path"],'metab_abund.txt')
        self.api_dic["anno_overview"].add_overview_detail(main_table_id, final_anno_overview,origin_exp_file=origin_exp_file)
        if self.ko:
            self.api_dic["anno_overview"].add_overview_ko(main_table_id, ko_table)
        self.objid2dir[str(main_table_id)] = self.link_up_map["AnnoOverview"]

    @time_count
    def export_pca(self):
        self.api_dic["pca"] = self.api.api('metabolome.exp_pca')
        for each in self.run_type3:
            toolname = "pca_" + each
            out_dir = globals()[toolname].output_dir
            params = self.dic_params["pca"]
            self.logger.info(params)
            self.logger.info(type(params))
            if self.option("project_type") == "LC" and self.is_mix_table == 'F':   # LC 且不合并才有
                params["table_type"] = each
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            name = "ExpPca_Origin_" + each
            group = self.group_table.prop["path"]
            pca_file = os.path.join(out_dir, "PCA.sites.xls")
            model_file = os.path.join(out_dir, "PCA.model.xls")
            ellipse = os.path.join(out_dir, "PCA.ellipse.xls")
            main_table_id = self.api_dic["pca"].add_exp_pca(name=name, params=params,table_type=each)
            self.api_dic["pca"].add_exp_pca_detail(main_table_id, pca_file, ellipse, group_file=group)
            self.api_dic["pca"].add_exp_pca_model(main_table_id, model_file)
            self.objid2dir[str(main_table_id)] = self.link_up_map["ExpPCA"] + '/' + each

    @time_count
    def export_sample_corr(self):
        self.api_dic["sample_corr"] = self.api.api("metabolome.exp_corr")
        for each in self.run_type3:
            toolname = "sample_corr_" + each
            out_dir = globals()[toolname].output_dir
            name = "ExpCorr_Origin_" + each
            params = self.dic_params["sample_corr"]
            if self.option("project_type") == "LC" and self.is_mix_table == 'F':   # LC 且不合并才有
                params["table_type"] = each
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            tree_file = os.path.join(out_dir, "corr.cluster_tree.xls")
            p_file = os.path.join(out_dir, "pvalue.xls")
            corr_file = os.path.join(out_dir, "corr.xls")
            main_table_id = self.api_dic["sample_corr"].add_exp_corr(name=name, tree_file=tree_file, params=params,table_type=each)
            self.api_dic["sample_corr"].add_exp_corr_detail(main_table_id, corr_file, p_file)
            self.objid2dir[str(main_table_id)] = self.link_up_map["ExpCorr"] + '/' + each

    #add v3
    @time_count
    def export_sample_venn(self):
        self.api_dic["sample_venn"] = self.api.api("metabolome.sample_venn")
        for each in self.run_type3:
            toolname = "sample_venn_" + each
            out_dir = globals()[toolname].output_dir
            name = "SampleVenn_Origin_" + each
            params = self.dic_params["sample_venn"]
            if self.option("project_type") == "LC" and self.is_mix_table == 'F':   # LC 且不合并才有
                params["table_type"] = each
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_table_id = self.api_dic["sample_venn"].add_venn(name, params)
            venn_table_result = out_dir + "/venn_table.xls"
            desc = os.path.join(self.preprocess.output_dir , "org_"+each+"/metab_desc.txt")
            self.api_dic["sample_venn"].add_venn_detail(venn_table_result, main_table_id, len(self.group_detail_noQC.keys()), desc)
            self.objid2dir[str(main_table_id)] = self.link_up_map["ExpVenn"] + '/' + each

    @time_count
    def export_sample_plsda(self):
        self.api_dic["sample_plsda"] = self.api.api("metabolome.sample_plsda")
        for each in self.run_type3:
            toolname = "sample_plsda_" + each
            pls_dir = globals()[toolname].output_dir
            name = "SamplePlsda_Origin_" + each
            params = self.dic_params["sample_plsda"]
            if self.option("project_type") == "LC" and self.is_mix_table == 'F':   # LC 且不合并才有
                params["table_type"] = each
            metab_table_id = params['metab_table']
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = self.api_dic["sample_plsda"].add_sample_plsda(name=name, params=params,metab_table_id=metab_table_id)
            self.api_dic["sample_plsda"].add_exp_diff_bar(main_id, pls_dir)
            self.api_dic["sample_plsda"].add_exp_diff_comp(main_id, pls_dir, self.option("group_table").prop['path'])
            self.api_dic["sample_plsda"].add_exp_diff_model(main_id, pls_dir)
            self.api_dic["sample_plsda"].add_exp_diff_scatter(main_id, pls_dir)
            self.objid2dir[str(main_id)] = self.link_up_map["ExpPLSDA"] + '/' + each

    #add v3
    @time_count
    def export_groups_diff(self):
        self.api_dic["groups_diff"] = self.api.api("metabolome.groups_diff")
        name = "GroupsDiff_Origin"
        params = self.dic_params["groups_diff"]
        metab_table_id = params["metab_table"]
        params = json.dumps(params, sort_keys=True, separators=(",", ":"))
        main_id = self.api_dic["groups_diff"].add_groups_diff(name=name, params=params, metab_table_id=metab_table_id)
        for each in self.run_type3:
            toolname = "groups_diff_" + each
            target_dir = globals()[toolname].output_dir
            self.api_dic["groups_diff"].add_groups_diff_detail(target_dir,main_id,table_type=each)
        self.objid2dir[str(main_id)] = self.link_up_map["GroupsDiff"]

    @time_count
    def export_diff_pls(self):
        self.api_dic["diff_pls"] = self.api.api("metabolome.exp_diff")
        name = "ExpDiff_Origin"
        params = self.dic_params["diff_pls"]
        params = json.dumps(params, sort_keys=True, separators=(",", ":"))
        if len(self.run_type3) ==1:
            tmp = self.link_up_map["ExpDiff"]+"/"+self.run_type3[0]+"/DiffTest/"
            diff_result_dir = os.path.join(self._sheet.output, tmp)
        else:
            tmp = self.link_up_map["ExpDiff"]+"/"+self.run_type3[0]+"/DiffTest/"
            tmp2 = self.link_up_map["ExpDiff"]+"/"+self.run_type3[1]+"/DiffTest/"
            diff_result_dir = os.path.join(self._sheet.output, tmp) + ','+ os.path.join(self._sheet.output, tmp2)
        main_table_id = self.api_dic["diff_pls"].add_exp_diff(self.metab_table_id, '', name=name, params=params,
                                                                  diff_dir=diff_result_dir)
        self.diff_pls_main_id = main_table_id
        self.objid2dir[str(main_table_id)] = self.link_up_map["ExpDiff"]

        for each in self.run_type3:    # 20190722 run_type1 change to run_type3
            toolname = "diff_pls_" + each
            out_dir = globals()[toolname].output_dir

            pls_dir = globals()[toolname].option("pls_dir").prop["path"]
            diff_file_dir = globals()[toolname].option("diffStat_dir").prop["path"]
            id_diff_file_dir = globals()[toolname].option("id_diffStat_dir").prop["path"]
            id_diff_plot_dir = globals()[toolname].option("diff_plot_dir").prop["path"]
            # diff_result_dir = rere_path + "/ExpDiff/" + each + "/DiffTest/"

            diff_main_id = "diff_pls_" + each
            globals()[diff_main_id] = main_table_id
            self.api_dic["diff_pls"].add_exp_diff_detail(main_table_id, id_diff_file_dir, each)
            self.api_dic["diff_pls"].add_exp_diff_detail_plot(main_table_id, id_diff_plot_dir, each)
            self.api_dic["diff_pls"].add_exp_diff_model(main_table_id, pls_dir,type=each)
            self.api_dic["diff_pls"].add_exp_diff_bar(main_table_id, pls_dir,type=each)
            self.api_dic["diff_pls"].add_exp_diff_comp(main_table_id, pls_dir, self.group_table.prop["path"],type=each)
            self.api_dic["diff_pls"].add_exp_diff_scatter(main_table_id, pls_dir,type=each)
            self.api_dic["diff_pls"].add_exp_diff_load(main_table_id, pls_dir, each) #20200316
            self.api_dic["diff_pls"].add_exp_diff_splot(main_table_id, pls_dir, each) #20200316

    @time_count
    def export_two_sample(self):
        self.api_dic["two_sample"] = self.api.api("metabolome.two_sample")
        params = self.dic_params['two_sample']
        params = json.dumps(params, sort_keys=True, separators=(",", ":"))
        name = 'TwoSample_Origin'
        main_id = self.api_dic['two_sample'].add_two_sample(name, params, self.metab_table_id)
        for each in self.run_type3:
            tool_name = 'two_sample_'+each
            output_dir = globals()[tool_name].output_dir
            self.api_dic["two_sample"].add_two_sample_detail(output_dir,main_id,table_type=each)
        self.objid2dir[str(main_id)] = self.link_up_map["TwoSample"]

    @time_count
    def export_metab_cluster(self):
        self.api_dic["metab_cluster"] = self.api.api("metabolome.metabset_cluster")
        for each in self.run_type2:
            name = "MetabsetCluster_Origin"
            params = self.dic_params["metab_cluster"]
            params_ori = params
            ###params["table_type"] = each  v2.0
            params_ori = params
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            toolname = "metabset_cluster_" + each
            out_dir = globals()[toolname].output_dir
            sam_tree = os.path.join(out_dir, "sample.cluster_tree.xls")
            matab_tree = os.path.join(out_dir, "metab.cluster_tree.xls")
            matab_tree_id = os.path.join(out_dir, "metab_id.cluster_tree.xls")
            expression_matrix = os.path.join(out_dir, "cluster_scale_exp.xls")
            main_table_id = self.api_dic["metab_cluster"].add_metabset_cluster(name=name, sam_tree=sam_tree,
                                                        matab_tree=matab_tree_id, params=params,table_type=each,metab_set=ObjectId(params_ori['metab_set']))
            self.api_dic["metab_cluster"].add_metabset_cluster_detail(main_table_id, expression_matrix)
            self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetCluster"] + "/DiffSet_mix"

    @time_count
    def export_metab_vip(self):
        self.api_dic["metab_vip"] = self.api.api("metabolome.metabset_vip")
        for each in self.run_type2:   #20190722 run_type1 change to run_type3
            name = "MetabsetVip_Origin_" + each
            params = self.dic_params["metab_vip"]
            params_ori = params
            diff_main_id = "diff_pls_" + each
            params["diff_table"] =  str(self.diff_pls_main_id)  #str(globals()[diff_main_id])
            ###params["table_type"] = each  v2.0
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            toolname = "metabset_vip_" + each
            out_dir = globals()[toolname].output_dir
            diff_detail = json.dumps(self.group_detail)
            diff_main_id = "diff_pls_" + each
            diff_id = self.diff_pls_main_id     #globals()[diff_main_id]
            main_table_id = self.api_dic["metab_vip"].add_metabset_vip(self.metab_table_id, params=params, name=name,
                                                            diff_detail=diff_detail, diff_id=diff_id,table_type=each,metab_set=ObjectId(params_ori['metab_set']))
            if each == 'pos':
                metab_desc = self.preprocess.option('pos_out').prop['path']+'/metab_desc.txt'
            elif each == 'neg':
                metab_desc = self.preprocess.option('neg_out').prop['path']+'/metab_desc.txt'
            else:
                metab_desc = self.preprocess.option('mix_out').prop['path']+'/metab_desc.txt'

            self.api_dic["metab_vip"].add_metabset_vip_detail(main_table_id, out_dir, "oplsda", metab_desc, scale="scale")
            self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetVip"] + "/DiffSet_mix/" + each

    @time_count
    def export_metab_corr(self):
        self.api_dic["metabset_corr"] = self.api.api("metabolome.metabset_corr")
        for each in self.run_type2:
            name = "MetabsetCorr_Origin"
            toolname = "metabset_corr_" + each
            params = self.dic_params["metab_corr"]
            params_ori = params
            ###params["table_type"] = each v2.0
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            out_dir = globals()[toolname].output_dir
            tree_file = os.path.join(out_dir, "corr.cluster_tree.xls")
            tree_id_file = os.path.join(out_dir, "metab_id.corr_tree.xls")
            p_file = os.path.join(out_dir, "pvalue.xls")
            corr_file = os.path.join(out_dir, "corr.xls")
            main_table_id = self.api_dic["metabset_corr"].add_metabset_corr(params=params, name=name,
                                                            tree_file=tree_id_file, list_file=corr_file,table_type=each,metab_set=ObjectId(params_ori['metab_set']))
            self.api_dic["metabset_corr"].add_metabset_corr_detail(main_table_id, corr_file, p_file)
            self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetCorr"] + "/DiffSet_mix"


    @time_count
    def export_metab_keggc(self):
        self.api_dic["metab_keggc"] = self.api.api('metabolome.metabset_keggc')
        name = "MetabsetKeggc_Origin"
        params = self.dic_params["metab_keggc"]
        main_table_id = self.api_dic["metab_keggc"].add_metabsetc(params, name=name)
        self.logger.info(self.metabset_keggc.option('stat_out'))
        stat_out = self.metabset_keggc.option('stat_out').prop["path"]
        self.api_dic["metab_keggc"].add_metabsetc_detail(main_table_id, stat_out)
        self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetKeggc"] + "/DiffSet_mix"

    @time_count
    def export_metab_hmdb(self):
        self.api_dic["metab_hmdb"] = self.api.api('metabolome.anno_hmdb')
        name = "MetabsetHmdb_Origin"
        params = self.dic_params["metab_hmdb"]
        params_ori = params
        main_table_id = self.api_dic["metab_hmdb"].add_anno_hmdb_main(name, params=params, database="metabsethmdb", table_type=self.run_type2[0],metab_set=ObjectId(params_ori['metabset']))
        level_out = self.metabset_hmdb.option('level_out').path
        out_dir = self.metabset_hmdb.output_dir
        self.api_dic["metab_hmdb"].add_level_detail(main_table_id, level_out, database="metabsethmdb")
        self.api_dic["metab_hmdb"].add_stat_detail(main_table_id, out_dir, database="metabsethmdb")
        self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetHmdb"] + "/DiffSet_mix"

    @time_count
    def export_metab_keggp(self):
        self.api_dic["metab_keggp"] = self.api.api('metabolome.metabset_keggp')
        params = self.dic_params["metab_keggp"]
        out_dir = self.metabset_keggp.output_dir
        # img_dir = self.remote_dir + "/MetabsetKeggp/pathway_img/"
        # final_img_dir = "rerewrweset" + img_dir.split("rerewrweset")[1]
        final_img_dir = os.path.join(self._sheet.output, self.link_up_map["MetabsetKeggp"]+"/pathway_img/")
        set_name = "DiffSet_mix"
        name = "MetabsetKeggp_Origin"
        main_table_id = self.api_dic["metab_keggp"].add_metabsetp(params, final_img_dir, set_name, name=name)
        self.logger.info(self.metabset_keggp.option('level_out'))
        self.logger.info(self.metabset_keggp.option('level_out').prop["path"])
        self.api_dic["metab_keggp"].add_metabsetp_level(main_table_id, self.metabset_keggp.option('level_out').path)
        stat_path = os.path.join(out_dir, 'stat.xls')
        self.api_dic["metab_keggp"].add_metabsetp_stat(main_table_id, set_name, stat_path)
        self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetKeggp"] + "/DiffSet_mix"

    @time_count
    def export_metab_enrich(self):
        self.api_dic["metab_enrich"] = self.api.api('metabolome.enrich')
        name = "MetabsetEnrich_Origin"
        params = self.dic_params["metab_enrich"]
        out_dir = self.metabset_enrich.output_dir
        if self.option("project_type") == "LC":
            overview_table = annotation_mix.option("overview_out").prop["path"]
        else:
            overview_table = annotation_pos.option("overview_out").prop["path"]
        enrich_table = os.path.join(out_dir, 'DE.list.check.kegg_enrichment.xls')
        if int(self.table_len(enrich_table)) > 1:
            main_table_id = self.api_dic["metab_enrich"].add_enrich(name, params,metab_set=self.diff_set_id)
            self.api_dic["metab_enrich"].add_enrich_detail(main_table_id, overview_table, enrich_table)
            topo_table = os.path.join(out_dir, 'topology_png/kegg_topology.xls')
            if os.path.exists(topo_table):
                remote_path = os.path.join(self._sheet.output,self.link_up_map["MetabsetEnrich"] +'/DiffSet_mix/topology_png')
                self.api_dic["metab_enrich"].add_topology(main_table_id,topo_table,remote_path)
            self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetEnrich"] + "/DiffSet_mix"

    @time_count
    def export_metab_ipath(self):
        self.api_dic["metab_ipath"] = self.api.api('metabolome.ipath')
        name = "MetabsetIpath_Origin"
        params = self.dic_params["metab_ipath"]
        # rere_path = self.remote_dir + "/MetabsetIpath/"
        # final_rere_path = "rerewrweset" + rere_path.split("rerewrweset")[1]
        final_rere_path = os.path.join(self._sheet.output, self.link_up_map["MetabsetIpath"] + '/DiffSet_mix')
        ipath_dir = self.metabset_ipath.output_dir
        insert_file = os.path.join(ipath_dir, 'gene_ipath_input.xls')
        if int(self.table_len(insert_file)) > 1:
            main_table_id = self.api_dic["metab_ipath"].add_ipath(name, final_rere_path, params=params)
            self.api_dic["metab_ipath"].add_ipath_detail(main_table_id, ipath_dir,
                                                         self.preprocess.option("org_set").prop["path"])
            self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetIpath"] + "/DiffSet_mix"

    @time_count
    def export_metab_venn(self):
        self.api_dic["metab_venn"] = self.api.api('metabolome.metabset_venn')
        name = "MetabsetVenn_Origin"
        params = self.dic_params["metab_venn"]
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        venn_path = os.path.join(self.metabset_venn.output_dir, "venn_table.xls")
        venn_graph_result = os.path.join(self.metabset_venn.output_dir, "metabset_detail.xls")
        main_table_id = self.api_dic["metab_venn"].add_venn(name=name, params=params)
        self.api_dic["metab_venn"].add_venn_detail(venn_path, main_table_id)
        self.api_dic["metab_venn"].add_venn_graph(venn_graph_result, main_table_id)
        self.objid2dir[str(main_table_id)] = self.link_up_map["MetabsetEnrich"] + "/DiffSet_mix"

    @time_count
    def export_org_metabset(self):
        self.api_dic["org_set"] = self.api.api('metabolome.metabset')
        org_metab_file = self.preprocess.option("has_name_org_mul_set").prop["path"]
        rm_repeat = self.api.api('metabolome.repeat')
        rm_repeat.rm_each("metab_set")
        with open(org_metab_file, "r") as f:
            line = f.readline().strip().split("\t")
            metabs = line[1]
            metablist = metabs.split(",")
        me_num = len(metablist)
        self.org_set_id = self.api_dic["org_set"].add_metab_set("Set_Raw", me_num,not_delete=True)
        self.api_dic["org_set"].add_metab_set_detail(self.org_set_id, metablist)
        self.logger.info("raw metabset finish!")

        if self.preprocess.option("has_name_origin_mul_set").is_set:
            origin_metab_file = self.preprocess.option("has_name_origin_mul_set").prop["path"]
            with open(origin_metab_file, "r") as f:
                line = f.readline().strip().split("\t")
                metabs = line[1]
                metablist = metabs.split(",")
                me_num = len(metablist)
                origin_set_id = self.api_dic["org_set"].add_metab_set("Set_Origin", me_num,not_delete=True)
                self.api_dic["org_set"].add_metab_set_detail(origin_set_id, metablist)
                self.logger.info("origin metabset finish!")


    @time_count
    def export_diff_metabset(self):
        self.diff_group_set_id = []
        self.api_dic["diff_set"] = self.api.api('metabolome.metabset')
        diff_metab_file = os.path.join(self.creat_analysis_table.output_dir, "merge_mul.metabset.xls")
        diff_set_list = self.creat_analysis_table.option("set_list").path
        self.diff_set_ids = {}  ##important

        if os.path.getsize(diff_metab_file) > 0:
            with open(diff_metab_file, "r") as f:
                for line in f:
                    line = line.strip().split("\t")
                    if len(line) < 2:
                        continue
                    metabs = line[1]
                    metablist = metabs.split(",")
                    if len(metablist) == 0:
                        continue
                    me_num = len(metablist)
                    set_name = line[0]
                    self.diff_each_set_id = self.api_dic["diff_set"].add_metab_set(set_name, me_num)
                    self.diff_group_set_id.append(str(self.diff_each_set_id))
                    self.diff_set_ids[set_name] = str(self.diff_each_set_id)
                    self.api_dic["diff_set"].add_metab_set_detail(self.diff_each_set_id, metablist)

        table = pd.read_table(diff_set_list, sep="\t",index_col=0,header=None)
        if len(table) < 1: #20201126  无表头，所以小于1是没有代谢物
            self.diff_set_id = self.org_set_id
            self.top_diff = self.default_top
        else:
            all_diff_metab = table.index.tolist()
            self.top_diff = len(all_diff_metab)
            if self.top_diff > self.default_top:
                self.top_diff = int(self.default_top)
            set_name = "DiffSet_mix"
            self.diff_set_id = self.api_dic["diff_set"].add_metab_set(set_name, len(all_diff_metab),not_delete=True)
            self.api_dic["diff_set"].add_metab_set_detail(self.diff_set_id, all_diff_metab)

        self.logger.info("Diff metabset finish!")

    def export_diff_metabset_pip_result(self):
        if self.project_type == 'GC':
            table_type = 'pos'
        else:
            table_type = 'mix'
        for set_name in self.diff_set_ids:
            self._export_one_diff_metabset_pip_result(set_name,table_type, self.diff_metabset_pip)


    def _export_one_diff_metabset_pip_result(self,set_name, table_type, metabset_tool):  #zouguanqing
        #export emtab_cluster

        out_dir = metabset_tool.output_dir + '/MetabsetCluster/' + set_name
        if os.listdir(out_dir):
            self.api_dic["metab_cluster"] = self.api.api("metabolome.metabset_cluster")
            name = "MetabsetCluster_Origin_" + set_name
            params = self.dic_params["metab_cluster"]
            params['metab_set'] = self.diff_set_ids[set_name]
            params_ori = params
            sam_tree = os.path.join(out_dir, "sample.cluster_tree.xls")
            matab_tree = os.path.join(out_dir, "metab.cluster_tree.xls")
            matab_tree_id = os.path.join(out_dir, "metab_id.cluster_tree.xls")
            expression_matrix = os.path.join(out_dir, "cluster_scale_exp.xls")
            exp_data = pd.read_table(expression_matrix,sep='\t',header=0)
            top_num = len(exp_data)
            params['top_meta'] = top_num
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_table_id = self.api_dic["metab_cluster"].add_metabset_cluster(name=name, sam_tree=sam_tree,
                                      matab_tree=matab_tree_id, params=params,table_type=table_type,metab_set=ObjectId(params_ori['metab_set']))
            self.api_dic["metab_cluster"].add_metabset_cluster_detail(main_table_id, expression_matrix)

        ##vip
        out_dir = metabset_tool.output_dir + '/MetabsetVip/' + set_name
        if os.listdir(out_dir):
            self.api_dic["metab_vip"] = self.api.api("metabolome.metabset_vip")
            name = "MetabsetVip_Origin_" + set_name
            params = self.dic_params["metab_vip"]
            params['metab_set'] = self.diff_set_ids[set_name]
            params_ori = params
            params["diff_table"] =  str(self.diff_pls_main_id)
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            diff_detail = json.dumps(self.group_detail)
            diff_id = self.diff_pls_main_id
            main_table_id = self.api_dic["metab_vip"].add_metabset_vip(self.metab_table_id, params=params, name=name,
                                    diff_detail=diff_detail, diff_id=diff_id,table_type=table_type,metab_set=ObjectId(params_ori['metab_set']))
            if table_type == 'pos':
                metab_desc = self.preprocess.option('pos_out').prop['path']+'/metab_desc.txt'
            elif table_type == 'neg':
                metab_desc = self.preprocess.option('neg_out').prop['path']+'/metab_desc.txt'
            else:
                metab_desc = self.preprocess.option('mix_out').prop['path']+'/metab_desc.txt'
            self.api_dic["metab_vip"].add_metabset_vip_detail(main_table_id, out_dir, "oplsda", metab_desc, scale="scale",out_dir=out_dir)

        ##corr
        out_dir = metabset_tool.output_dir + '/MetabsetCorr/' + set_name
        if os.listdir(out_dir):
            self.api_dic["metabset_corr"] = self.api.api("metabolome.metabset_corr")
            name = "MetabsetCorr_Origin_" + set_name
            params = self.dic_params["metab_corr"]
            params['metab_set'] = self.diff_set_ids[set_name]
            corr_file = os.path.join(out_dir, "corr.xls")
            corr_data = pd.read_table(corr_file, sep='\t', header=0)
            params['top'] = len(corr_data)
            params_ori = params
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))

            tree_file = os.path.join(out_dir, "corr.cluster_tree.xls")
            tree_id_file = os.path.join(out_dir, "metab_id.corr_tree.xls")
            p_file = os.path.join(out_dir, "pvalue.xls")

            main_table_id = self.api_dic["metabset_corr"].add_metabset_corr(params=params, name=name,
                                                            tree_file=tree_id_file, list_file=corr_file,table_type=table_type,metab_set=ObjectId(params_ori['metab_set']))
            self.api_dic["metabset_corr"].add_metabset_corr_detail(main_table_id, corr_file, p_file)

        ##keggc
        out_dir = metabset_tool.output_dir + '/MetabsetKeggc/' + set_name
        if os.listdir(out_dir):
            self.api_dic["metab_keggc"] = self.api.api('metabolome.metabset_keggc')
            name = "MetabsetKeggc_Origin_" + set_name
            params = self.dic_params["metab_keggc"]
            params['metabset'] = self.diff_set_ids[set_name]
            main_table_id = self.api_dic["metab_keggc"].add_metabsetc(params, name=name)
            out_dir = metabset_tool.output_dir + '/MetabsetKeggc/' + set_name
            stat_out = out_dir+ '/stat.xls'
            self.api_dic["metab_keggc"].add_metabsetc_detail(main_table_id, stat_out)

        ##hmdb
        out_dir = metabset_tool.output_dir + '/MetabsetHmdb/' + set_name
        if os.listdir(out_dir):
            self.api_dic["metab_hmdb"] = self.api.api('metabolome.anno_hmdb')
            name = "MetabsetHmdb_Origin_" + set_name
            params = self.dic_params["metab_hmdb"]
            params['metabset'] = self.diff_set_ids[set_name]
            params_ori = params
            main_table_id = self.api_dic["metab_hmdb"].add_anno_hmdb_main(name, params=params, database="metabsethmdb", table_type=table_type,metab_set=ObjectId(params_ori['metabset']))

            level_out = out_dir +'/HmdbLevel.xls'
            self.api_dic["metab_hmdb"].add_level_detail(main_table_id, level_out, database="metabsethmdb")
            self.api_dic["metab_hmdb"].add_stat_detail(main_table_id, out_dir, database="metabsethmdb")

        ##keggp
        out_dir = metabset_tool.output_dir + '/MetabsetKeggp/' + set_name
        if os.listdir(out_dir):
            self.api_dic["metab_keggp"] = self.api.api('metabolome.metabset_keggp')
            params = self.dic_params["metab_keggp"]
            params['metabset'] = self.diff_set_ids[set_name]

            final_img_dir = os.path.join(self._sheet.output, self.link_up_map["MetabsetKeggp"]+"/"+set_name+"/pathway_img/")
            name = "MetabsetKeggp_Origin_" + set_name
            main_table_id = self.api_dic["metab_keggp"].add_metabsetp(params, final_img_dir, set_name, name=name)
            level_out = out_dir + '/level.xls'
            self.api_dic["metab_keggp"].add_metabsetp_level(main_table_id, level_out)
            stat_path = os.path.join(out_dir, 'stat.xls')
            self.api_dic["metab_keggp"].add_metabsetp_stat(main_table_id, set_name, stat_path)


        ##enrich
        out_dir = metabset_tool.output_dir + '/MetabsetEnrich/' + set_name
        if os.listdir(out_dir):
            self.api_dic["metab_enrich"] = self.api.api('metabolome.enrich')
            name = "MetabsetEnrich_Origin_" + set_name
            params = self.dic_params["metab_enrich"]
            params['metabset'] = self.diff_set_ids[set_name]
            out_dir = metabset_tool.output_dir + '/MetabsetEnrich/' + set_name

            overview_table = self.anno_overview.path
            enrich_table = os.path.join(out_dir, 'DE.list.check.kegg_enrichment.xls')
            if int(self.table_len(enrich_table)) > 1:
                main_table_id = self.api_dic["metab_enrich"].add_enrich(name, params,metab_set=self.diff_set_ids[set_name])
                self.api_dic["metab_enrich"].add_enrich_detail(main_table_id, overview_table, enrich_table)
                topo_table = os.path.join(out_dir, 'topology_png/kegg_topology.xls')
                if os.path.exists(topo_table):
                    remote_path = os.path.join(self._sheet.output,self.link_up_map["MetabsetEnrich"]+'/'+set_name +'/topology_png')
                    self.api_dic["metab_enrich"].add_topology(main_table_id,topo_table,remote_path)

        ##ipath
        ipath_dir = metabset_tool.output_dir + '/MetabsetIpath/' + set_name
        out_dir = ipath_dir
        if os.listdir(out_dir):
            self.api_dic["metab_ipath"] = self.api.api('metabolome.ipath')
            name = "MetabsetIpath_Origin_" + set_name
            params = self.dic_params["metab_ipath"]
            params['metabset'] = self.diff_set_ids[set_name]
            final_rere_path = os.path.join(self._sheet.output, self.link_up_map["MetabsetIpath"],set_name)

            insert_file = os.path.join(ipath_dir, 'gene_ipath_input.xls')
            if int(self.table_len(insert_file)) > 1:
                main_table_id = self.api_dic["metab_ipath"].add_ipath(name, final_rere_path, params=params)
                self.api_dic["metab_ipath"].add_ipath_detail(main_table_id, ipath_dir,self.preprocess.option("org_set").prop["path"])
                self.objid2dir[str(main_table_id)] = self.link_up_map[out_dir.split('/')[-2]] + set_name

    def main_collection_delete(self, func_name):
        """
        重导表检查，去除重复主表
        """
        #myApiManager =  ApiManager(self)
        rm_repeat = self.api.api('metabolome.repeat')
        collect_relation = {
            "export_preprocess": "exp",
            "export_anno_keggc": "anno_keggc",
            "export_anno_keggp": "anno_keggp",
            "export_anno_overview": "anno_overview",
            "export_pca": "exp_pca",
            "export_sample_corr": "exp_corr",
            "export_diff_pls": "exp_diff",
            "export_metab_cluster": "metabset_cluster",
            "export_metab_vip": "metabset_vip",
            "export_metab_corr": "metabset_corr",
            "export_metab_keggc": "metabset_keggc",
            "export_metab_keggp": "metabset_keggp",
            "export_metab_enrich": "metabset_kegg_enrich",
            "export_metab_ipath": "metabset_ipath",
            "export_metab_venn": "metabset_venn",
            "export_anno_hmdb": "anno_hmdb",
            "export_metab_hmdb": "metabset_hmdb",
            "export_sample_venn": "sample_venn",  #add v3
            "export_sample_plsda": "sample_plsda",
            "export_two_sample" : "diff_sample",
            "export_groups_diff" : "groups_diff",

        }
        collect_name = collect_relation[func_name]
        if func_name == "export_preprocess":
            rm_repeat.rm_main(collect_name, name="MetabTable_Origin")
            rm_repeat.rm_main(collect_name, name="raw")
        else:
            #params_name = func_name.replace("export_","")
            #params = self.dic_params[params_name]
            rm_repeat.rm_main(collect_name)
            #print "repeat main_collection {} delete".format(collect_name)
            #self.logger.info("repeat main_collection {} delete".format(collect_name))

    def database_version(self):
        """
        功能：根据注释版本返回任务id
        工作流默认跑最新的，返回任务id为空
        如果跑老数据库，则需要返回任务id
        :return:
        """
        if self.option("kegg_version") in ['KEGG_V94.2', 'KEGG_V94', 'kegg_v94.2', "kegg_v94"]:
            db_v = "v94.2"
        elif self.option("kegg_version") in ['KEGG', 'kegg']:
            db_v = 'kegg'
        elif self.option("kegg_version").lower() in ['kegg_v2021.09.18', 'v2021.09.18']:
            db_v = 'v2021.09.18'
        return db_v


    def move_dir_upload(self, olddir,  newdir):  #zouguanqing 20190802
        all = os.listdir(olddir)
        for each in all:
            file = os.path.join(olddir,each)
            new_file = os.path.join(newdir,each)
            if os.path.isdir(file):
                if not os.path.exists(new_file):
                    os.mkdir(new_file)
                    self.move_dir_upload(file,new_file)
            else:
                if os.path.exists(new_file):
                    os.remove(new_file)
                os.link(file,new_file)

    def up_files(self):
        #dir_up = self.output_dir
        link_up_map = self.link_up_map
        metab_pip_dir = [
            'MetabsetCluster',
            'MetabsetVip',
            'MetabsetCorr',
            'MetabsetKeggc',
            'MetabsetKeggp',
            'MetabsetEnrich',
            'MetabsetHmdb' ,
            'MetabsetIpath'
        ]


        dir_up = self.work_dir +'/upload_dir'
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.makedirs(dir_up)
        output_dirs = os.listdir(self.output_dir)
        for dir in output_dirs:
            if dir in link_up_map.keys():
                new_dir = os.path.join(dir_up,link_up_map[dir])
                if dir in metab_pip_dir + ['MetabsetVenn']:
                    new_dir = new_dir + '/DiffSet_mix'
                shutil.copytree(self.output_dir+'/'+dir, new_dir)

        ##将metaset 的pip 结果 run_diff_metabset_pip 结果链接到 upload文件夹
        for dir in metab_pip_dir:
            for set_name  in self.diff_set_ids:
                ori_dir = os.path.join(self.diff_metabset_pip.output_dir, dir,set_name)
                tar_dir = os.path.join(dir_up,self.link_up_map[dir],set_name)
                if not os.path.exists(ori_dir):
                    continue
                shutil.copytree(ori_dir,tar_dir)
                if dir == 'MetabsetKeggc':
                    os.remove(tar_dir+'/level.xls')
                    os.rename(tar_dir+'/stat.xls',tar_dir+'/level.xls')  #保持文件夹文件和原来的统一

        #删除不上传的文件
        metab_lists = glob.glob(dir_up +'/1.Preprocess/*.list')
        qc_xls = glob.glob(dir_up+'/1.Preprocess/*.xls')
        for f in metab_lists + qc_xls:
            os.remove(f)

        # 打包压缩pathway_img
        anno_img = dir_up + '/3.Anno/02.AnnoKeggp/pathway_img'
        anno_set_imgs =glob.glob(dir_up + '/5.Metabset/06.MetabsetKeggp/*/pathway_img')
        anno_set_imgs.append(anno_img)
        wf_dir = os.path.abspath(os.curdir)
        try:
            for img in anno_set_imgs:
                os.chdir(img + '/..')
                os.system('tar -zcf pathway_img.tar.gz pathway_img')
                os.system('rm -r pathway_img')
                os.chdir(wf_dir)
        except:
            os.chdir(wf_dir)
            #self.set_error("pathway_img 压缩失败")
            self.logger.info("pathway_img 压缩失败, 请检查原因")
        self.upload_results()

    def upload_results(self):
        dir_up = self.work_dir +'/upload_dir'
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, dir_up))

        repaths = [
            [".", "", "流程分析结果目录"],
            ["1.Preprocess", "", "预处理结果文件夹", 0, "150001"],
            ["1.Preprocess/QCcv/", "", "QC CV评估图", 0, ""],
            ["1.Preprocess/pos/", "", "预处理阳离子结果", 0, "150002"],
            ["1.Preprocess/pos/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["1.Preprocess/pos/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["1.Preprocess/neg/", "", "预处理阴离子结果", 0, "150005"],
            ["1.Preprocess/neg/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["1.Preprocess/neg/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["1.Preprocess/mix/", "", "预处理混合离子结果", 0, "150006"],
            ["1.Preprocess/mix/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["1.Preprocess/mix/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["1.Preprocess/org_neg/", "", "原始阴离子表", 0, "150007"],
            ["1.Preprocess/org_neg/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["1.Preprocess/org_neg/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["1.Preprocess/org_pos/", "", "原始阳离子表", 0, "150008"],
            ["1.Preprocess/org_pos/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["1.Preprocess/org_pos/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["1.Preprocess/org_mix/", "", "原始合并表", 0],
            ["1.Preprocess/org_mix/metab_abund.txt", "txt", "离子强度表", 0, "150003"],
            ["1.Preprocess/org_mix/metab_desc.txt", "txt", "代谢物信息表", 0, "150004"],
            ["2.SampleComp", "", "样本相关分析", 0],
            ["2.SampleComp/01.ExpCorr", "", "样本相关性结果文件夹", 0, "150009"],
            ["2.SampleComp/01.ExpCorr/pos/sample_corr_tree.xls", "xls", "样本相关性树文件", 0, "150010"],
            ["2.SampleComp/01.ExpCorr/pos/sample_corr.xls", "xls", "样本相关性表", 0, "150011"],
            ["2.SampleComp/01.ExpCorr/pos/corr_pvalue.xls", "xls", "样本相关性P值表", 0, "150012"],
            ["2.SampleComp/01.ExpCorr/neg/sample_corr_tree.xls", "xls", "样本相关性树文件", 0, "150010"],
            ["2.SampleComp/01.ExpCorr/neg/sample_corr.xls", "xls", "样本相关性表", 0, "150011"],
            ["2.SampleComp/01.ExpCorr/neg/corr_pvalue.xls", "xls", "样本相关性P值表", 0, "150012"],
            ['2.SampleComp/02.ExpPCA', '', 'PCA结果文件夹', 0, "150013"],
            ['2.SampleComp/02.ExpPCA/pos/PCA.loadings.xls', 'xls', 'PCA代谢物主成分贡献度表', 0, "150014"],
            ['2.SampleComp/02.ExpPCA/pos/PCA.model.xls', 'xls', 'PCA模型参数表', 0, "150015"],
            ["2.SampleComp/02.ExpPCA/pos/PCA.sites.xls", "xls", "PCA样本各维度坐标", 0, "150016"],
            ["2.SampleComp/02.ExpPCA/neg/PCA.loadings.xls", "xls", "PCA代谢物主成分贡献度表", 0, "150014"],
            ["2.SampleComp/02.ExpPCA/neg/PCA.sites.xls", "xls", "样本各维度坐标", 0, "150016"],
            ["2.SampleComp/02.ExpPCA/neg/PCA.model.xls", "xls", "PCA模型参数表", 0, "150015"],
            ["2.SampleComp/03.ExpVenn","", '样本Venn结果文件夹',0],  #add v3
            ["2.SampleComp/03.ExpVenn/pos", "",'阳离子Venn结果',0],
            ["2.SampleComp/03.ExpVenn/neg", "",'阴离子Venn结果',0],
            ["2.SampleComp/04.ExpPLSDA","", '样本Plsda结果文件夹',0],
            ["2.SampleComp/04.ExpPLSDA/pos","", '阳离子Plsda结果文件夹',0],
            ["2.SampleComp/04.ExpPLSDA/neg","", '阴离子Plsda结果文件夹',0],
            ["3.Anno", "", "注释结果", 0],
            ["3.Anno/01.AnnoKeggc", "", "KEGG化合物分类注释结果", 0, "150017"],
            ["3.Anno/01.AnnoKeggc/KEGG_Compound_Classification.pdf", "pdf", "化合物分类统计柱状图", 0, ""],
            ["3.Anno/01.AnnoKeggc/level.xls", "xls", "化合物分类层级表", 0, "150018"],
            ["3.Anno/01.AnnoKeggc/stat.xls", "xls", "化合物分类统计表", 0, "150019"],
            ["3.Anno/02.AnnoKeggp/", "", "KEGG功能通路结果", 0, "150020"],
            ["3.Anno/02.AnnoKeggp/KEGG_Pathway.pdf", "pdf", "KEGG通路统计图", 0, ""],
            ["3.Anno/02.AnnoKeggp/Top20_KEGG_pathway.pdf", "pdf", "重要通路统计图", 0, ""],
            ["3.Anno/02.AnnoKeggp/level.xls", "xls", "Pathway层级表", 0, "150062"],
            ["3.Anno/02.AnnoKeggp/stat.xls", "xls", "KEGG通路统计表", 0, "150063"],
            ["3.Anno/02.AnnoKeggp/pathway_img/", "", "KEGG通路图", 0, "150021"],
            ["3.Anno/04.AnnoOverview/", "", "代谢物总览", 0, "150022"],
            ["3.Anno/04.AnnoOverview/anno.xls", "xls", "代谢物信息总览表", 0, "150023"],
            ["3.Anno/03.AnnoHmdb", "", "HMDB化合物分类结果", 0, "150082"],
            ["3.Anno/03.AnnoHmdb/HmdbLevel.xls", "xls", "HMDB化合物分类层级表", 0, "150083"],
            ["3.Anno/03.AnnoHmdb/HmdbSuperclass.xls", "xls", "HMDB Superclas统计表", 0, "150084"],
            ["3.Anno/03.AnnoHmdb/HmdbClass.xls", "xls", "HMDB CLass统计表", 0, "150085"],
            ["3.Anno/03.AnnoHmdb/HmdbSubclass.xls", "xls", "HMDb Subclass统计表", 0, "150086"],
            ["4.ExpDiff", "", "差异代谢物计算与统计结果文件夹", 0, "150048"],
            ["4.ExpDiff/01.TwoGroupExpDiff","","两组差异代谢物计算结果文件夹",0],
            ["4.ExpDiff/02.TwoSamExpDiff","","两样本差异代谢物计算结果文件夹",0],
            ["4.ExpDiff/03.MultiGroupExpDiff","","多组差异代谢物计算结果文件夹",0],
            ["4.ExpDiff/03.MultiGroupExpDiff/pos","","阳离子结果文件夹",0],
            ["4.ExpDiff/03.MultiGroupExpDiff/neg","","阴离子结果文件夹",0],
            ["4.ExpDiff/03.MultiGroupExpDiff/mix","","Mix结果文件夹",0],
            ["5.Metabset", "", "代谢集结果", 0],
            ["5.Metabset/01.MetabsetVenn/", "", "Venn分析结果给", 0, "150024"],
            ["5.Metabset/02.MetabsetCluster", "", "代谢物聚类结果文件夹", 0, "150027"],
            ["5.Metabset/03.MetabsetVip", "", "VIP分析结果文件夹", 0, "150040"],
            ["5.Metabset/04.MetabsetCorr", "", "代谢物相关性结果文件", 0, "150031"],
            ["5.Metabset/05.MetabsetKeggc", "", "代谢集KEGG化合物分类结果", 0, "150035"],
            ["5.Metabset/06.MetabsetKeggp", "xls", "代谢集KEGG功能通路结果", 0, "150037"],
            ["5.Metabset/07.MetabsetEnrich", "", "代谢集KEGG富集分析结果", 0, "150042"],
            ["5.Metabset/08.MetabsetHmdb", "", "代谢集HMDB化合物分类结果", 0, "150087"],
            ["5.Metabset/09.MetabsetIpath", "", "iPath代谢通路分析结果", 0, "150044"],
            ["5.Metabset/10.MetabsetRoc", "", "ROC结果", 0 ],
            ["5.Metabset/19.IntergratedCorr", "", "Corr结果", 0],
            ["5.Metabset/20.IntergratedProc", "", "Proc结果", 0],

        ]
        regexps = [
            #[r"AnnoKeggp/pathway_img/.*\.pdf", "pdf", "KEGG通路图"],
            #[r"MetabsetKeggp/pathway_img/.*\.pdf", "pdf", "代谢集KEGG通路图"],
            [r"2.SampleComp/01.ExpCorr/.*/SampleCorr.*.pdf", "pdf", "样本相关性热图", 0, ""],
            [r"2.SampleComp/02.ExpPCA/.*/SamplePCA.*.pdf", "pdf", "样本PCA分析图", 0, ""],
            [r"2.SampleComp/03.ExpVenn/.*/SampleVenn.*.pdf","pdf", '样本比较Venn图',0], #add v3
            [r"2.SampleComp/03.ExpVenn/.*/venn_table.xls","xls", 'venn结果',0], #add v3
            [r"2.SampleComp/03.ExpVenn/.*/metabset_detail.xls","xls", '样本代谢集',0],
            [r"2.SampleComp/04.ExpPLSDA/.*/SamplesScorePLSDA.pdf", "pdf",'代谢样本比较PLSDA图',0],
            [r"2.SampleComp/04.ExpPLSDA/.*/SamplesPLSDA_ModelOverview.pdf", "pdf",'样本比较PLSDA模型总览图',0],
            [r"2.SampleComp/04.ExpPLSDA/.*/SamplesOPLSDA_ModelOverview.pdf", "pdf",'样本比较PLSDA模型总览图',0],
            [r"2.SampleComp/04.ExpPLSDA/.*/.*loadings.xls", "xls",'代谢物载荷图数据',0],
            [r"2.SampleComp/04.ExpPLSDA/.*/.*model.xls","xls", '模型参数',0],
            [r"2.SampleComp/04.ExpPLSDA/.*/.*permMN.xls","xls",'模型参数',0],
            [r"2.SampleComp/04.ExpPLSDA/.*/.*sites.xls","xls", '散点图数据',0],
            [r"3.Anno/03.AnnoHmdb/HMDB_Compound_.*.pdf", "pdf", "HMDB化合物分类统计图", 0, "150083"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/", "", "多元统计模型结果", 0, "150050"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*ScorePCA.pdf ", "pdf", "PCA分析散点图", 0, "150051"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*ScorePCAloading.pdf ", "pdf", "PCA分析载荷图", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*ScorePLSDA.pdf ", "pdf", "Plsda分析散点图", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*ScorePLSDA_loading.pdf ", "pdf", "Plsda分析载荷图", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*PLSDA_ModelOverview.pdf ", "pdf", "Plsda模型概览", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*PLSDA_Permutation.pdf ", "pdf", "Plsda置换检验", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*ScoreOPLSDA.pdf ", "pdf", "OPlsda分析散点图", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*ScoreOPLSDA_loading.pdf ", "pdf", "OPlsda分析载荷图", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*OPLSDA_ModelOverview.pdf ", "pdf", "OPlsda模型概览", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*OPLSDA_Permutation.pdf ", "pdf", "OPlsda置换检验", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/.*ScoreOPLSDA_Splot.pdf ", "pdf", "OPlsda分析splot图", 0, ""],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/OPLS-DA\.model\.xls", "xls", "OPLS-DA模型参数表", 0, "150051"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/OPLS-DA\.permMN\.xls", "xls", "OPLS-DA响应排序检验结果表", 0, "150052"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/OPLS-DA\.loadings\.xls", "xls", "OPLS-DA代谢物主成分贡献度表", 0, "150053"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/OPLS-DA\.sites\.xls", "xls", "OPLS-DA样本各维度坐标", 0, "150054"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/OPLS-DA\.vips\.xls", "xls", "OPLS-DA的VIP值表", 0, "150055"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/PCA\.loadings\.xls", "xls", "PCA代谢物主成分贡献度表", 0, "150014"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/PCA\.model\.xls", "xls", "PCA模型参数表", 0, "150015"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/PCA\.sites\.xls", "xls", "PCA样本各维度坐标", 0, "150016"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/PLS-DA\.loadings\.xls", "xls", "PLS-DA代谢物主成分贡献度表", 0, "150056"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/PLS-DA\.model\.xls", "xls", "PLS-DA模型参数表", 0, "150057"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/PLS-DA\.permMN\.xls", "xls", "PLS-DA响应排序检验结果表", 0, "150058"],
            [r"4.ExpDiff/01.TwoGroupExpDifff/.+/DiffMulStat/.*_vs_.*/PLS-DA\.sites\.xls", "xls", "PLS-DA样本各维度坐标", 0, "150059"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffMulStat/.*_vs_.*/PLS-DA\.vips\.xls", "xls", "PLS-DA的VIP值表", 0, "150060"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffTest/.*_vs_.*\.diff\.exp\.xls", "xls", "差异检验结果表", 0, "150061"],
            [r"4.ExpDiff/01.TwoGroupExpDiff/.+/DiffTest/.*Volcano.pdf", "pdf", "Volcano 火山图", 0, "150061"],
            [r"4.ExpDiff/02.TwoSamExpDiff/pos", "","阳离子表结果",0],
            [r"4.ExpDiff/02.TwoSamExpDiff/neg", "","阴离子表结果",0],
            [r"4.ExpDiff/02.TwoSamExpDiff/.*/.*/", "","比较组结果文件夹",0],
            [r"4.ExpDiff/02.TwoSamExpDiff/.*/.*/.*.xls", "","结果文件",0],
            [r"4.ExpDiff/03.MultiGroupExpDiff/.*/.*.xls","xls","结果文件",0],
            ["5.Metabset/01.MetabsetVenn/.*/venn_table.xls", "xls", "Venn分析结果表", 0, "150025"],
            ["5.Metabset/01.MetabsetVenn/.*/metabset_detail.xls", "xls", "各代谢集包含的代谢物列表", 0, "150026"],
            ["5.Metabset/02.MetabsetCluster/.*", "", "", 0],
            ["5.Metabset/02.MetabsetCluster/.*/metab.cluster_tree.xls", "xls", "代谢物聚类树文件", 0, "150028"],
            ["5.Metabset/02.MetabsetCluster/.*/sample.cluster_tree.xls", "xls", "样本聚类树文件", 0, "150029"],
            ["5.Metabset/02.MetabsetCluster/.*/cluster_exp.xls", "xls", "代谢物表达量", 0, "150030"],
            ["5.Metabset/02.MetabsetCluster/.*/cluster_scale_exp.xls", "xls", "标准化后代谢物表达量", 0, "150077"],
            ["5.Metabset/02.MetabsetCluster/.*/metab_id.cluster_tree.xls", "xls", "代谢物id聚类树文件", 0, "150075"],
            [r"5.Metabset/02.MetabsetCluster/.*/.*Heatmap.pdf", "pdf", "聚类热图", 0, ""],
            [r"5.Metabset/02.MetabsetCluster/.*/.*Heatmap.+.pdf", "pdf", "子聚类趋势图", 0, ""],
            ["5.Metabset/03.MetabsetVip/.*", "", "", 0],
            [r"5.Metabset/03.MetabsetVip/.+/.*_vs_.*/.*VIPbar.pdf ", "pdf", "VIP分析bar图", 0, "150049"],
            [r"5.Metabset/03.MetabsetVip/.+/.*_vs_.*/.*VIPplot.pdf ", "pdf", "VIP分析plot图", 0, "150049"],
            [r"5.Metabset/03.MetabsetVip/.+/.*_vs_.*/Vip_exp\.xls", "xls", "代谢物VIP值表", 0, "150049"],
            [r"5.Metabset/03.MetabsetVip/.+/.*_vs_.*/metab\.cluster_tree\.xls", "xls", "代谢物聚类树文件", 0, "150028"],
            [r"5.Metabset/03.MetabsetVip/.+/.*_vs_.*/metab_id\.cluster_tree\.xls", "xls", "代谢物id聚类树文件", 0, "150075"],
            [r"5.Metabset/03.MetabsetVip/.+/.*_vs_.*/Vip_scale_exp\.xls", "xls", "代谢物VIP值与标准化后表达量", 0, "150078"],
            ["5.Metabset/04.MetabsetCorr/.*", "", "", 0,],
            ["5.Metabset/04.MetabsetCorr/.*/.*Correlation.pdf", "pdf", "代谢物相关性热图", 0, ""],
            ["5.Metabset/04.MetabsetCorr/.*/metab_corr_tree.xls", "xls", "代谢物相关性树文件", 0, "150032"],
            ["5.Metabset/04.MetabsetCorr/.*/metab_corr.xls", "xls", "代谢物相关性表", 0, "150033"],
            ["5.Metabset/04.MetabsetCorr/.*/corr_pvalue.xls", "xls", "代谢物相关性P值表", 0, "150034"],
            ["5.Metabset/04.MetabsetCorr/.*/metab_id.corr_tree.xls", "xls", "代谢物id聚类树文件", 0, "150076"],
            ["5.Metabset/05.MetabsetKeggc/.*", "", 0],
            ["5.Metabset/05.MetabsetKeggc/.*/.*KEGG_Compound_Classification.pdf", "pdf", "化合物分类统计图", 0, ""],
            ["5.Metabset/05.MetabsetKeggc/.*/level.xls", "xls", "代谢集KEGG化合物分类层级表", 0, "150036"],
            ["5.Metabset/06.MetabsetKeggp/.*", "xls", "代谢集KEGG功能通路结果", 0, "150037"],
            ["5.Metabset/06.MetabsetKeggp/.*/.*KEGG_Pathway.pdf", "pdf", "KEGG统计图", 0, "150038"],
            ["5.Metabset/06.MetabsetKeggp/.*/level.xls", "xls", "代谢集KEGG通路分类层级表", 0, "150038"],
            ["5.Metabset/06.MetabsetKeggp/.*/stat.xls", "xls", "代谢集KEGG通路统计表", 0, "150039"],
            ["5.Metabset/06.MetabsetKeggp/.*/pathway_img", "", "代谢集KEGG通路图", 0, "150041"],
            ["5.Metabset/07.MetabsetEnrich.*", "", "", 0],
            ["5.Metabset/07.MetabsetEnrich/.*/DE.list.check.kegg_enrichment.xls", "xls", "代谢集KEGG富集分析结果表", 0, "150043"],
            ["5.Metabset/07.MetabsetEnrich/.*/.*KEGG_topology.pdf", "pdf", "KEGG拓扑学分析气泡图", 0, ""],
            ["5.Metabset/07.MetabsetEnrich/.*/.*KEGG_enrich_bar_p-ER.pdf", "pdf", "KEGG富集分析图", 0, ""],
            ["5.Metabset/07.MetabsetEnrich/.*/.*KEGG_enrich_bar_ER.pdf", "pdf", "KEGG富集分析图", 0, ""],
            ["5.Metabset/07.MetabsetEnrich/.*/.*KEGG_enrich_bar_p.pdf", "pdf", "KEGG富集分析图", 0, ""],
            ["5.Metabset/07.MetabsetEnrich/.*/.*KEGG_enrich_bubble-p.pdf", "pdf", "KEGG富集分析图", 0, ""],
            ["5.Metabset/07.MetabsetEnrich/.*/.*KEGG_enrich_bubble-ER.pdf", "pdf", "KEGG富集分析图", 0, ""],
            ["5.Metabset/09.MetabsetIpath/.*", "", 0],
            ["5.Metabset/09.MetabsetIpath/.*/Biosynthesis_of_secondary_metabolities.svg", "svg",
             "Biosynthesis of secondary metabolities", 0, "150045"],
            ["5.Metabset/09.MetabsetIpath/.*/Metabolic_pathways.svg", "svg", "Metabolic pathways", 0, "150046"],
            ["5.Metabset/09.MetabsetIpath/.*/Regulatory_pathways.svg", "svg", "Regulatory pathways", 0, "150047"],
            ["5.Metabset/08.MetabsetHmdb/.*", "", 0],
            ["5.Metabset/08.MetabsetHmdb/.*/.*HMDB_Compound_Classification.*.pdf", "pdf", "HMDB化合物分类统计图", 0, "150088"],
            ["5.Metabset/08.MetabsetHmdb/.*/HmdbLevel.xls", "xls", "代谢集HMDB化合物分类层级表", 0, "150088"],
            ["5.Metabset/08.MetabsetHmdb/.*/HmdbSuperclass.xls", "xls", "代谢集HMDB Superclas统计表", 0, "150089"],
            ["5.Metabset/08.MetabsetHmdb/.*/HmdbClass.xls", "xls", "代谢集HMDB CLass统计表", 0, "150090"],
            ["5.Metabset/08.MetabsetHmdb/.*/HmdbSubclass.xls", "xls", "代谢集HMDb Subclass统计表", 0, "150091"],
            ["5.Metabset/10.MetabsetRoc/.*/.*.pdf", "pdf", "代谢物ROC曲线图", 0, ""],
            [r"5.Metabset/11.MetabsetEnrichHeatmap/.*/.*Heatmap.pdf", "pdf", "聚类热图", 0, ""],
            [r"5.Metabset/11.MetabsetEnrichHeatmap/.*/.*Heatmap.+.pdf", "pdf", "子聚类趋势图", 0, ""],
            [r"5.Metabset/19.IntergratedCorr/.*.pdf", "pdf", "代谢与表型相关性热图", 0],
            [r"5.Metabset/20.IntergratedProc/.*.pdf", "pdf", "代谢与表型普式分析图", 0],
        ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
