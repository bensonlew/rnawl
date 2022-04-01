# -*- coding:utf-8 -*-
# __author__ = 'fengyitong'
"""dia工作流"""

from biocluster.workflow import Workflow
from mbio.workflows.labelfree.labelfree import LabelfreeWorkflow
# from mbio.workflows.dia.labelfree import LabelfreeWorkflow
import os
import shutil
import re
import pandas as pd
from biocluster.wpm.client import *
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.ref_rna_v2.functions import tryforgood
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo


class DiaWorkflow(LabelfreeWorkflow):
    def __init__(self, wsheet_object):
        """
        workflow option参数设置
        """
        self._sheet = wsheet_object
        super(LabelfreeWorkflow, self).__init__(wsheet_object)
        options = [
            ##选择输入文件
            {"name": "protein", "type": "infile", "format": "labelfree.common"},  # 输入蛋白鉴定表
            {"name": "ratio_exp", "type": "infile", 'format': "labelfree.ratio_exp"},  # 输入表达量ratio表
            {"name": "protein_fasta", "type": "infile", "format": "labelfree.common"},  # 输入蛋白fasta序列文件
            {"name": "protein_group", "type": "infile", "format": "labelfree.group_table"},  # 输入分组文件
            {"name": "protein_control", "type": "infile", "format": "labelfree.compare_table"},  # 输入对照组文件

            ##功能注释参数
            {"name": "data_source", "type": "string", "default": "Uniprot"}, # FASTA序列来源
            {"name": "go_evalue", "type": "string", "default": "1e-5"}, # go注释evalue
            {"name": "go_identity", "type": "float", "default": 0.98}, # go注释identity
            {"name": "cog_evalue", "type": "string", "default": "1e-5"}, # cog注释evalue
            {"name": "cog_identity", "type": "float", "default": 0}, # cog注释identity
            {"name": "kegg_class", "type": "string", "default": ""}, # kegg注释物种分类
            {"name": "kegg_org", "type": "string", "default": ""}, # kegg注释具体物种
            {"name": "kegg_evalue", "type": "string", "default": "1e-5"}, # kegg注释evalue
            {"name": "kegg_identity", "type": "float", "default": 0.98}, # kegg注释identity
            {"name": "sub_loc", "type": "string", "default": "Plants"}, # 亚细胞定位 物种分类
            {"name": "pfam_evalue", "type": "string", "default": "1e-5"}, # pfam注释evalue
            {"name": "database", "type": "string", "default": 'go,kegg,pfam,cog'},
            {"name": "nr_database", "type": "string", "default": "All"},  # nr库类型
            {"name": "nr_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "pfam_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-5},
            {"name": "gram", "type": "string", "default": 'neg'}, # 如果是原核的话，格兰阴氏和格兰阳氏会选择不同的训练集

            ##差异表达分析参数
            {"name": "fc_up", "type": "float", "default": 1.2},
            {"name": "fc_down", "type": "float", "default": 0.83},
            {"name": "pvalue", "type": "float", "default": 0.05},
            {"name": "correct_method", "type": "string", "default":"two.sided"}, #卡方检验的页面没有单双尾检验
            {"name": "mul_test", "type": "string", "default": "none"},# 统一用p.adjust来矫正p值
            # param mul_test: 多重检验方法选择，默认为none，包括: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
            {"name": "method_type", "type": "string", "default": "student"},
            # 升级后新增了填充缺失值和选择每组有效样本数的功能
            {"name": "fillna", "type": "string", "default": "none"},
            {"name": "cutoffs", "type": "string", "default": "none"},

            ##样本比较分析
            {"name": "sam_analysis", "type": "bool", "default": True},# 设置页面端是否显示样本比较分析
            {"name": "network_analyse", "type": "bool", "default": True},  # 设置页面端是否显示ppi分析
            {"name": "string_analyse", "type": "bool", "default": True},  # 设置页面端是否显示string爬虫分析
            {"name": "WGCNA_analyse", "type": "bool", "default": True},  # 设置页面端是否显示WGCNA分析

            ##蛋白互作网络
            {"name": "ppi_category", "type": "string", "default": "All"},# 设置页面端是否显示样本比较分析
            # {"name": "ppi_species", "type": "string", "default": "Homo ""sapiens"}, # 设置PPI物种
            {"name": "ppi_species", "type": "int", "default": 0}, # 设置PPI物种

            ##是否是DIA数据
            {"name": "DIA", "type": "bool", "default": True},#如果选择DIA会传入这个参数，默认为False
            {"name": "report_dia", "type": "infile", "format": "labelfree.common"},#选择dia会需要上传的report文件
            {"name": "exp_ana_dia", "type": "infile", "format": "labelfree.common"},#选择dia会需要上传的实验配置文件，用来生成蛋白信息
            {"name": "psm", "type": "string", "default": ""},#用来欺骗labelfree的文件检查，心累
            {"name": "protein_information", "type": "string", "default": ""},#用来欺骗labelfree的文件检查，心累
            {"name": "peptide", "type": "string", "default": ""},#用来欺骗labelfree的文件检查，心累
            # {"name": "protein_dia", "type": "infile", "format": "labelfree.common"},#选择dia会需要上传的搜库文件

            {"name": "change_des", "type": "bool", "default": False},#是否选择用nr注释结果替换description信息
            {"name": "useblast", "type": "bool", "default": True},#是否选择用string注释结果爬取官网
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_202007"},
        ]

        #获取输出目录
        self.workflow_output_tmp = self._sheet.output
        try:
            if re.match(r'tsanger:',self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
            elif re.match(r'sanger:',self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
            elif re.match(r'^\w+://\S+/.+$',self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp
            else:
                self.set_error("json output wrong")
        except:
            self.workflow_output = 'just_test'
        self.project_sn = self._sheet.project_sn #获取project_sn
        self.task_id = self._sheet.id #获取task_id
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False

        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()
        # 骗过文件检查
        self.option('psm', self.option('report_dia').prop['path'])
        self.option('protein_information', self.option('report_dia').prop['path'])
        self.option('peptide', self.option('report_dia').prop['path'])

        # 事先处理那3个表格，如果是excel，则转化成txt
        try:
            protein_df = pd.read_csv(self.option('protein').prop['path'])
        except:
            protein_df = pd.read_excel(self.option('protein').prop['path'], sep='\t')
            protein_path = os.path.join(os.path.dirname(self.option('protein').prop['path']), 'protein_.xls')
            protein_df.to_csv(protein_path,sep='\t',header=True,index=False)
            self.option('protein', protein_path)
        del protein_df

        # labelfree更换搜库软件后要兼容两种格式
        self.new = False
        with open(self.option('protein').prop['path'], 'r') as pro:
            header = pro.readline()
            if header.startswith('Checked'):
                self.new = True

        #添加tool/module
        self.filecheck = self.add_tool("labelfree.filecheck_labelfree")
        if self.new:
            self.searchdb = self.add_tool("itraq_and_tmt.searchdb")
        else:
            self.searchdb = self.add_tool("labelfree.searchdb")
        self.sam_corr = self.add_tool("labelfree.exp_corr")
        # self.sam_pca = self.add_tool("labelfree.exp_pca")
        self.sam_pca = self.add_tool("labelfree.exp_pca_meta")
        # self.sam_pca = self.add_tool("labelfree.pca")
        # 置信圈
        self.ellipse = self.add_tool("graph.ellipse")
        self.diff_pep = self.add_tool("labelfree.diff")
        if self.new:
            self.annotation = self.add_module("itraq_and_tmt.protein_annotation")
        else:
            self.annotation = self.add_module("labelfree.protein_annotation")

        # 工作流中跑蛋白集相关内容
        self.exp_venn = self.add_tool("labelfree.exp_venn")
        self.diff_cluster = self.add_module("labelfree.diff_cluster")
        self.diff_cog_class = self.add_module("labelfree.diff_cog_class")
        self.diff_go_class = self.add_module("labelfree.diff_go_class")
        self.diff_kegg_class = self.add_module("labelfree.diff_kegg_class")
        self.diff_pfam_stat = self.add_module("labelfree.diff_pfam_stat")
        self.diff_subloc_stat = self.add_module("labelfree.diff_subloc_stat")
        self.diff_ipath = self.add_module("labelfree.diff_ipath")
        self.diff_ppi = self.add_module("labelfree.diff_ppi")
        self.diff_string_picture = self.add_module("labelfree.diff_string_pictures")
        self.diff_go_enrich = self.add_module("labelfree.diff_go_enrich")
        self.diff_kegg_enrich = self.add_module("labelfree.diff_kegg_enrich")
        self.diff_circle = self.add_module("labelfree.diff_circle")

        #判断流程结束tool/module list
        self.final_tools = [self.searchdb, self.sam_corr, self.sam_pca, self.exp_venn, self.diff_cluster, self.diff_cog_class,
                            self.diff_go_class, self.diff_kegg_class, self.diff_pfam_stat, self.diff_subloc_stat, self.diff_ipath, self.diff_ppi, self.diff_circle]
        # self.final_tools = [self.searchdb, self.sam_corr, self.sam_pca, self.exp_venn, self.diff_cluster,
        #                     self.diff_go_class, self.diff_kegg_class, self.diff_ipath, self.diff_ppi, self.diff_circle]

        # 添加step，显示在页面进度条
        all_steps = ["filecheck", "searchdb", "exp_venn", "sam_corr", "sam_pca", "diff_pep",
                     "annotation", "diff_cluster", "diff_cog_class", "diff_go_class",
                     "diff_kegg_class", "diff_pfam_stat", "diff_subloc_stat", "diff_ipath", "diff_ppi", "diff_string_picture",
                     "diff_go_enrich", "diff_kegg_enrich", "diff_circle"]
        # all_steps = ["filecheck", "searchdb", "exp_venn", "sam_corr", "sam_pca", "diff_pep",
        #              "annotation", "diff_cluster", "diff_go_class",
        #              "diff_kegg_class", "diff_ipath", "diff_ppi",
        #              "diff_go_enrich", "diff_kegg_enrich", "diff_circle"]
        # self.step.add_steps("filecheck", "searchdb", "sam_corr", "sam_pca", "diff_pep","annotation")

        group_spname = self.option("protein_group").prop['group_dict']
        group_snum = [len(group_spname[g]) for g in group_spname]
        self.min_group_num = min(group_snum)
        if self.min_group_num >= 3:
            self.final_tools[self.final_tools.index(self.sam_pca)] = self.ellipse
            all_steps.append("ellipse")

        if self.option("protein_group").prop["sample_number"] < 3:
            all_steps.remove('sam_pca')
            all_steps.remove('sam_corr')
            self.final_tools.remove(self.sam_pca)
            self.final_tools.remove(self.sam_corr)
            try:
                all_steps.remove("ellipse")
                self.final_tools.remove(self.ellipse)
            except:
                pass
        for step in all_steps:
            self.step.add_steps(step)

        # 事先处理表达量表
        self.ref_file =self.treat()

        qc = os.path.join(self.work_dir, 'qc')
        if os.path.exists(qc):
            shutil.rmtree(qc)
        os.makedirs(qc)

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self.task_id, 'labelfree')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))

import unittest
import datetime
class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        data = {
            "id": "dia_workflow_fff_" + id,
            "type": "workflow",
            "name": "dia.dia",
            "options": dict(
                protein="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/protein.xlsx",
                ratio_exp="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/exp.txt",
                protein_fasta="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/exp.fasta",
                protein_group="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/group.txt",
                protein_control="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/control.txt",
                report_dia="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/report.xls",
                exp_ana_dia="/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/20190421_185436_SLD_ALL_SUM_ExperimentAnalysis.txt",
                )
        }

        info = worker.add_task(data)
        print(info)


if __name__ == '__main__':
    unittest.main()
