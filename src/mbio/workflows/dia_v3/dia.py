# -*- coding:utf-8 -*-
# __author__ = 'xuxi'
""" dia v2工作流 ; last modified 20201117 """

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import glob
import json
import shutil
import re
import time
import gevent
import functools
from bson.son import SON
from biocluster.config import Config
from collections import OrderedDict
import pandas as pd
import numpy as np
from biocluster.wpm.client import *
from mbio.packages.rna.annot_config import AnnotConfig
import tarfile
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo


# 定义用于统计导表时间的装饰器
def time_count(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run ' + func_name + ' at ' + start_time)
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End ' + func_name + ' at ' + end_time)
        print("{}函数执行时间约为{}s".format(func.__name__, end - start))
    return wrapper

class DiaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        workflow option参数设置
        """
        self._sheet = wsheet_object
        super(DiaWorkflow, self).__init__(wsheet_object)
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
            {"name": "method_type", "type": "string", "default": "t.test"},
            # 升级后新增了填充缺失值和选择每组有效样本数的功能
            {"name": "fillna", "type": "string", "default": "none"},
            {"name": "cutoffs", "type": "string", "default": "none"},            
            # {"name": "remove_eigenvalues", "type": "int", "default": 90},
            # {"name": "ingroup_eigenvalues", "type": "int", "default": 50},
            # {"name": "fillna_samples", "type": "string", "default": "ingroup"},
            {"name": "log", "type": "string", "default": "none"},   #  新增参数 20201125
            ## 数据预处理参数
            {"name": "all_eliminate", "type": "string", "default": "all"},  #
            {"name": "all_percent", "type": "string", "default": "90"},
            {"name": "if_group", "type": "string", "default": "yes"},   # if perform specific
            {"name": "group_specific", "type": "string", "default": "any"},
            {"name": "group_percent", "type": "string", "default": "50"},
            {"name": "fillna", "type": "string", "default": "seqknn"},  # method for fill na
            {"name": "fill_type", "type": "string", "default": "group"},    # based on all/group  ；预处理和run_diff_pep共用参数


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
            {"name": "report_img", "type": "bool", "default": True},
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

        # 骗过文件检查
        self.option('psm', self.option('report_dia').prop['path'])
        self.option('protein_information', self.option('report_dia').prop['path'])
        self.option('peptide', self.option('report_dia').prop['path'])

        # 项目重运行时，清除mongo库中已有数据
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False

        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

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

        # self.filecheck = self.add_tool("labelfree.filecheck_labelfree")
        self.filecheck = self.add_tool("dia_v3.filecheck_dia")
        self.preprocess = self.add_tool("dia_v3.preprocess")
        # if self.new:
        #     self.searchdb = self.add_tool("itraq_and_tmt.searchdb")
        # else:
        self.searchdb = self.add_tool("dia_v3.searchdb")
        self.sam_corr = self.add_tool("labelfree.exp_corr")
        # self.sam_pca = self.add_tool("labelfree.exp_pca")
        self.sam_pca = self.add_tool("labelfree.exp_pca_meta")
        # self.sam_pca = self.add_tool("labelfree.pca")
        # 置信圈
        self.ellipse = self.add_tool("graph.ellipse")
        self.diff_pep = self.add_tool("dia_v3.diff")

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
        self.final_tools = [self.searchdb, self.preprocess, self.sam_corr, self.sam_pca, self.exp_venn, self.diff_cluster, self.diff_cog_class,
                            self.diff_go_class, self.diff_kegg_class, self.diff_pfam_stat, self.diff_subloc_stat, self.diff_ipath, self.diff_ppi, self.diff_circle]
        # self.final_tools = [self.searchdb, self.sam_corr, self.sam_pca, self.exp_venn, self.diff_cluster,
        #                     self.diff_go_class, self.diff_kegg_class, self.diff_ipath, self.diff_ppi, self.diff_circle]

        # 添加step，显示在页面进度条
        all_steps = ["filecheck", "preprocess", "searchdb", "exp_venn", "sam_corr", "sam_pca", "diff_pep",
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
        self.raw_ref_file =self.treat()

        qc = os.path.join(self.work_dir, 'qc')
        if os.path.exists(qc):
            shutil.rmtree(qc)
        os.makedirs(qc)

    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self.task_id, 'dia')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        """
        检查选项
        """
        # 注释相关参数
        try:
            go_evalue = float(self.option("go_evalue"))
            kegg_evalue = float(self.option("kegg_evalue"))
            pfam_evalue = float(self.option("pfam_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范")
        else:
            self.option("nr_blast_evalue", go_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("pfam_blast_evalue", pfam_evalue)
        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("GO比对的E值超出范围")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围")
        if not self.option("pfam_blast_evalue") > 0 and not self.option("pfam_blast_evalue") < 1:
            raise OptionError("Pfam比对的E值超出范围")
        if not self.option("go_identity") >= 0 and not self.option("go_identity") <= 1:
            raise OptionError("GO identity值超出范围")
        if not self.option("kegg_identity") >= 0 and not self.option("kegg_identity") <= 1:
            raise OptionError("KEGG identity值超出范围")
        if self.option("nr_database") == "All":
            self.option("nr_database", "nr")
        else:
            nr = self.option("nr_database").lower()
            self.option("nr_database", nr)

        # 使用分库代替nr的总库
        if self.option('kegg_class') == 'Animals':
            self.option('nr_database', 'animal')
        elif self.option('kegg_class') == 'Plants':
            self.option('nr_database', 'plant')
        elif self.option('kegg_class') == 'Fungi':
            self.option('nr_database', 'fungi')
        elif self.option('kegg_class') == 'Protists':
            self.option('nr_database', 'protist')


        if self.option("method_type") == "t.test":
            self.option("method_type", "student")
        elif self.option("method_type") == "chisq.test":
            self.option("method_type", "chi")
        elif self.option("method_type") == "fisher.test":
            self.option("method_type", "fisher")
        # if self.option('DIA'):
        #     if not all(
        #             os.path.exists(self.option('report_dia').prop['path']),
        #             os.path.exists(self.option('exp_ana_dia').prop['path']),
        #             os.path.exists(self.option('protein_dia').prop['path']),
        #                ):
        #         raise OptionError("DIA分析需要的文件不全")

        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        """
        labelfree workflow run方法
        """
        self.filecheck.on('end', self.run_preprocess)
        self.preprocess.on('end', self.run_searchdb)
        if not self.option("protein_group").prop["sample_number"] < 3:
            self.preprocess.on('end', self.run_sam_pca)
            self.preprocess.on('end', self.run_sam_corr)
        self.preprocess.on('end', self.run_diff_pep)
        # self.preprocess.on('end', self.run_diamond)
        self.diff_pep.on('end', self.run_diamond)
        self.preprocess.on('end', self.run_exp_venn)
        self.on_rely([self.diff_kegg_enrich, self.diff_go_enrich], self.run_diff_circle)
        self.on_rely(self.final_tools, self.run_chart)
        self.run_filecheck()
        super(DiaWorkflow, self).run()

    def run_plot_static_peperr(self):
        # 将上行生成的dmass文件绘制成静态图(dMass.pdf,self.plot_static_peperr.output_dir+dMass.png)，因为页面的图片加载太大太慢了
        self.plot_static_peperr = self.add_tool("itraq_and_tmt.draw_staticplot")
        self.logger.info("开始运行绘制peperr静态图")
        opts = {
            "inputfile" : os.path.join(self.work_dir, 'tmp_qc','dMass.xls')
        }
        self.plot_static_peperr.set_options(opts)
        self.plot_static_peperr.on('end', self.end)
        self.plot_static_peperr.run()

    def run_chart(self):
        '''
        绘图步骤插入在导表前
        '''

        qc_dir = os.path.join(self.work_dir, 'tmp_qc')
        if not os.path.exists(qc_dir):
            os.makedirs(qc_dir)
        self.tmp_add_proteinmw(proteinmw_exp=self.option("protein").prop['path'])
        self.tmp_add_coverage(coverage_exp=self.option("protein").prop['path'])
        self.tmp_add_proteininfo(proteininfo_exp=self.option('exp_ana_dia').prop['path'])
        self.tmp_add_pepnum(pepnum_exp=self.option("protein").prop['path'])
        self.tmp_add_peplen(peplen_exp=self.option('report_dia').prop['path'])
        # self.tmp_add_peperror(peperror_exp=self.option("psm").prop['path'])
        self.chart = self.add_tool('dia_v3.chart')
        chart_dict = {
            # "type": "workflow",
            # "diff_group":[i[1]+"_vs_"+i[0] for i in self.option("protein_control").prop["cmp_list"]],
            "protein_mw_distribution": self.work_dir + "/tmp_qc/Protein_mw_distribution.xls",
            "protein_seq_cover_distribution": self.work_dir + "/tmp_qc/Protein_seq_cover_distribution.xls",
            "protein_infomation": self.work_dir + "/tmp_qc/Protein_infomation.xls",
            "peptite_number_distribution": self.work_dir + "/tmp_qc/Peptite_number_distribution.xls",
            "peptite_length_distribution": self.work_dir + "/tmp_qc/Peptite_length_distribution.xls",
            # "peptide_dmass": self.work_dir + "/tmp_qc/dMass.xls",
            "all_annotation_stat": self.annotation.output_dir + "/anno_stat/all_annotation_statistics.xls",
            "go12level_statistics": self.annotation.output_dir + "/go/go12level_statistics.xls",
            "go123level_statistics": self.annotation.output_dir + "/go/go123level_statistics.xls",
            "go1234level_statistics": self.annotation.output_dir + "/go/go1234level_statistics.xls",
            "kegg_layer": self.annotation.output_dir + "/kegg/kegg_layer.xls",
            "kegg_pathway": self.annotation.output_dir + "/kegg/pathway_table.xls",
            "cog_summary": self.annotation.output_dir + "/cog/cog_summary.xls",
            "pfam_domain": self.annotation.output_dir + "/blast_xml/pfam_domain",
            "subloc_stat": self.annotation.output_dir + "/subloc/multiloc_stat.xls",
            "num_summary": self.diff_pep.work_dir + "/num_summary.xls",
            "allsummary": self.diff_pep.work_dir + "/allsummary.xls",
            "diff": self.diff_pep.work_dir + "/output/{group_name}_diff.xls",

            "expression_matrix": self.diff_cluster.output_dir + "/{group_name}/expression_matrix.xls",
            "seq_tree": self.diff_cluster.output_dir + "/{group_name}/seq.cluster_tree.txt",
            "sample_tree": self.diff_cluster.output_dir + "/{group_name}/sample.cluster_tree.txt",

            "proteinsetgo_class_all": self.diff_go_class.output_dir + "/{group_name}/{group_name}_all_go_class_table.xls",
            "proteinsetgo_class_updown": self.diff_go_class.output_dir + "/{group_name}/{group_name}_up_down_go_class_table.xls",
            #
            "go_enrich_all": self.diff_go_enrich.output_dir + "/{group_name}_all/go_enrich_{group_name}_all_protein.xls",
            "go_enrich_down": self.diff_go_enrich.output_dir + "/{group_name}_down/go_enrich_{group_name}_down_protein.xls",
            "go_enrich_up": self.diff_go_enrich.output_dir + "/{group_name}_up/go_enrich_{group_name}_up_protein.xls",
            #
            "kegg_class_level": self.diff_kegg_class.work_dir + "/protein_kegg_level_table.xls",
            "kegg_class_stat_all": self.diff_kegg_class.work_dir + "/output/{group_name}_all/kegg_stat.xls",
            "kegg_class_stat_updown": self.diff_kegg_class.work_dir + "/output/{group_name}_up_down/kegg_stat.xls",
            #
            "kegg_enrich_all": self.diff_kegg_enrich.output_dir + "/{group_name}_all/{group_name}_all_protein.list.DE.list.check.kegg_enrichment.xls",
            "kegg_enrich_down": self.diff_kegg_enrich.output_dir + "/{group_name}_down/{group_name}_down_protein.list.DE.list.check.kegg_enrichment.xls",
            "kegg_enrich_up": self.diff_kegg_enrich.output_dir + "/{group_name}_up/{group_name}_up_protein.list.DE.list.check.kegg_enrichment.xls",
            #
            "geneset_circ_choose": self.diff_circle.output_dir + "/{{anno_type}}_circle_{group_name}_{{all_up_down}}/{{anno_type}}_enrich_choose.table",
            "geneset_circ_zscore": self.diff_circle.output_dir + "/{{anno_type}}_circle_{group_name}_{{all_up_down}}/enrich_zscore",
            #
            "ppi_centrality_all": self.diff_ppi.output_dir + "/{group_name}_all/ppinetwork_topology/protein_interaction_network_centrality.txt",
            "ppi_centrality_up": self.diff_ppi.output_dir + "/{group_name}_up/ppinetwork_topology/protein_interaction_network_centrality.txt",
            "ppi_centrality_down": self.diff_ppi.output_dir + "/{group_name}_down/ppinetwork_topology/protein_interaction_network_centrality.txt",
            "ppi_degree_all": self.diff_ppi.output_dir + "/{group_name}_all/ppinetwork_topology/protein_interaction_network_degree_distribution.txt",
            "ppi_degree_up": self.diff_ppi.output_dir + "/{group_name}_up/ppinetwork_topology/protein_interaction_network_degree_distribution.txt",
            "ppi_degree_down": self.diff_ppi.output_dir + "/{group_name}_down/ppinetwork_topology/protein_interaction_network_degree_distribution.txt",
            # #
            "cog_all": self.diff_cog_class.output_dir + "/{group_name}/{group_name}_all_cog_class_table.xls",
            "cog_up_down": self.diff_cog_class.output_dir + "/{group_name}/{group_name}_up_down_cog_class_table.xls",
            #
            "pfam_all": self.diff_pfam_stat.output_dir + "/{group_name}_all/pfam_stat.xls",
            "pfam_up_down": self.diff_pfam_stat.output_dir + "/{group_name}_up_down/pfam_stat.xls",
            #
            "subloc_all": self.diff_subloc_stat.output_dir + "/{group_name}_all/subloc_stat.xls",
            "subloc_up_down": self.diff_subloc_stat.output_dir + "/{group_name}_up_down/subloc_stat.xls",
            # # 下面画的两个图所在的页面 位于 样本比较分析（样本相关性热图，pca分析），这个分析模块的展示是可选的
            "sample_corr_tree": self.sam_corr.work_dir + "/sample.cluster_tree.txt",
            "sample_corr_matrix": self.sam_corr.work_dir + "/sample_correlation.xls",
            #
            "sample_pca_sites": self.sam_pca.output_dir + "/pca_sites.xls",
            "sample_pca_propor": self.sam_pca.output_dir + "/pca_importance.xls",

            "exp_venn": self.exp_venn.output_dir + '/venn_graph.xls',

            # 新增 数据预处理
            'preprocess_assess': os.path.join(self.preprocess.work_dir, '{}_nrmse.txt'.format(self.option('fillna'))),
            'preprocess_box_raw': os.path.join(self.preprocess.work_dir, 'boxdata_raw.xls'),
            'preprocess_box_origin': os.path.join(self.preprocess.work_dir, 'boxdata.xls'),
            'preprocess_summary_raw': os.path.join(self.preprocess.work_dir, '{}_cv_summary_raw.xls'.format(self.option('fillna'))),
            'preprocess_summary_origin': os.path.join(self.preprocess.work_dir, '{}_cv_summary.xls'.format(self.option('fillna'))),

            'group_file': self.option('protein_group').prop['path'],


            # # 下面这个venn图是新增的，以前的页面上都是没有这个图的
            # "diff_path_for_proteinset_for_venn": self.diff_pep.work_dir,

            # # #
            # above is pdf generated by main workflow
        }
        # if hasattr(self, 'ellipse'):
        #     chart_dict.update({
        #         "exp_pca_ellipse": "{table_dir}/ellipse_out.xls".format(table_dir=self.ellipse.output_dir)
        #     })
        chart_dict2 = {k: v for k, v in chart_dict.items() if os.path.exists(v) or "{" in v}
        chart_dict2.update({"diff_path_for_proteinset_for_venn": self.diff_pep.work_dir, "type": "workflow",
                            "diff_group": [i[1] + "_vs_" + i[0] for i in
                                           self.option("protein_control").prop["cmp_list"]]})
        if self.option('protein_group').is_set:
            chart_dict2.update({'protein_group': json.dumps(self.option('protein_group').prop['group_dict'], sort_keys=True, separators=(',', ':'))})
        with open(self.work_dir + "/chart_workflow.json", 'w') as json_f:
            json.dump(chart_dict2, json_f, sort_keys=True, indent=4)
        self.chart.set_options({
            "file_json": self.work_dir + "/chart_workflow.json"
        })
        self.chart.on('end', self.end)
        # self.chart.on('end', self.run_plot_static_peperr)
        self.chart.run()

    def run_filecheck(self):
        self.logger.info("开始运行文件检查")
        opts = {
            'protein_group': self.option('protein_group'),
            'protein_control': self.option('protein_control'),
            'ratio_exp': self.option('ratio_exp'),
            #'scaled_exp': self.option('scaled_exp'),
            'protein': self.option('protein'),
            'psm': self.option('psm'),
            'protein_information': self.option('protein_information'),
            'peptide': self.option('peptide'),
            'protein_fasta': self.option('protein_fasta'),
        }
        self.filecheck.set_options(opts)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_preprocess(self):
        self.logger.info("开始运行数据预处理")
        self.preprocess.set_options({
            'raw_path': self.work_dir + "/raw_treat_ref",
            'group_table': self.option("protein_group").prop["path"],
            'all_eliminate':self.option('all_eliminate'),
            'all_percent':float(self.option('all_percent')),
            'if_group':self.option('if_group'),
            'group_specific':self.option('group_specific'),
            'group_percent':float(self.option('group_percent')),
            'fillna':self.option('fillna'),
            'fill_type':self.option('fill_type'),
        })
        self.preprocess.on('end', self.set_preprocess_output)
        self.preprocess.on('start', self.set_step, {'start': self.step.preprocess})
        self.preprocess.on('end', self.set_step, {'end': self.step.preprocess})
        self.preprocess.run()

    def set_preprocess_output(self):
        if not os.path.exists(self.work_dir + "/treat_ref"):
            os.link(self.preprocess.option('preprocess_exp').prop['path'], self.work_dir + "/treat_ref")
        self.ref_file = self.work_dir + "/treat_ref"
        if not os.path.exists(self.work_dir + "/treat_ref"):
            raise OptionError("预处理失败，没有找到预处理结果文件")


    def run_searchdb(self):
        self.logger.info("开始运行搜库")
        self.searchdb.set_options({
            'protein': self.option('protein').prop['path'],
            # 'ratio_exp': self.preprocess.option('preprocess_exp').prop['path'],
            'ratio_exp': self.work_dir + "/raw_treat_ref",
            'protein_fasta':self.option('protein_fasta').prop['path'],
            'report_dia':self.option('report_dia').prop['path'],
        })
        self.searchdb.on('end', self.set_output, 'searchdb')
        self.searchdb.on('start', self.set_step, {'start': self.step.searchdb})
        self.searchdb.on('end', self.set_step, {'end': self.step.searchdb})
        self.searchdb.run()

    def treat(self):
        # 先暂定不取log
        df = pd.read_table(self.option("ratio_exp").prop["path"], header=0, index_col=0)
        df = df.fillna(0)
        df_scaled = df.replace([r'-', r'_', 'Filtered'], 0)
        df_scaled = df_scaled.ix[~(df_scaled == 0).all(axis=1), :]
        ref_file = self.work_dir + "/raw_treat_ref"
        df_scaled.to_csv(ref_file, sep = '\t')
        ### 注释下面，用preprocess代替其之功能
        # 如果选择填充缺失值的话，这边会进行，为了少些点代码，这里直接在工作流里跑了，幸亏数据肯定不会很大
        # if self.option('fillna').lower() not in ['no', 'none']:
        #     self.logger.info('选择了对表达量数据进行预处理，有效的蛋白的0值将被填充%s' % self.option('fillna').lower())
        #     from mbio.packages.labelfree.preprocess import Preprocess
        #     prep_obj = Preprocess(ref_file, self.option('protein_group').prop['path'],
        #                           method=self.option('fillna'), cutoffs=self.option('cutoffs'), out=ref_file)
        #     prep_obj.fillzero()
        #     prep_obj.to_out()
        return ref_file

    def run_sam_corr(self):
        self.logger.info("开始运行样本相关性分析")
        self.sam_corr.set_options({'exp': self.preprocess.option('preprocess_exp').prop['path']})
        self.sam_corr.on('end', self.set_output, 'sam_corr')
        self.sam_corr.on('start', self.set_step, {'start': self.step.sam_corr})
        self.sam_corr.on('end', self.set_step, {'end': self.step.sam_corr})
        self.sam_corr.run()

    def run_sam_pca(self):
        self.logger.info("开始运行样本PCA分析")
        # self.sam_pca.set_options({'otutable': self.ref_file})
        # self.sam_pca.set_options({'exp': self.ref_file})
        self.sam_pca.set_options(dict(
            exp_file=self.preprocess.option('preprocess_exp').prop['path'],
            group_file=self.option('protein_group').prop['path'],
        ))
        self.sam_pca.on('end', self.set_output, 'sam_pca')
        if self.min_group_num >= 3:
            self.sam_pca.on('end', self.run_ellipse)
        self.sam_pca.on('start', self.set_step, {'start': self.step.sam_pca})
        self.sam_pca.on('end', self.set_step, {'end': self.step.sam_pca})
        self.sam_pca.run()

    def run_ellipse(self):
        self.logger.info("开始运行pca ellipse")
        opts = {
            'analysis': "pca",
            # 'pc_table': self.sam_pca.output_dir + "/PCA.xls",
            'pc_table': self.sam_pca.output_dir + "/pca_sites.xls",
            'group_table':self.option('protein_group').prop['path']
        }

        self.ellipse.set_options(opts)
        self.ellipse.on("end", self.set_output, "ellipse")
        self.ellipse.on('start', self.set_step, {'start': self.step.ellipse})
        self.ellipse.on('end', self.set_step, {'end': self.step.ellipse})
        self.ellipse.run()

    def run_exp_venn(self):
        import itertools
        self.logger.info("开始运行组间venn分析")
        group_file = self.option('protein_group').prop['path']
        group2samples = self.option('protein_group').prop['group_dict']
        if len(group2samples) > 5:
            group_file += '_for_exp_venn'
            groups = group2samples.keys()
            with open(group_file, 'w') as gw:
                gw.write('#sample\tgroup\n')
                for n in itertools.count(1):
                    if n == 6:
                        break
                    group = groups[n]
                    for sample in group2samples[group]:
                        gw.write(sample + '\t' + group + '\n')

        opt = {
            "express_matrix" : self.preprocess.option('preprocess_exp').prop['path'],
            "group_table" : group_file
        }
        self.exp_venn.set_options(opt)
        self.exp_venn.on('end', self.set_output, 'exp_venn')
        self.exp_venn.on('start', self.set_step, {'start': self.step.exp_venn})
        self.exp_venn.on('end', self.set_step, {'end': self.step.exp_venn})
        self.exp_venn.run()

    def run_diff_pep(self):
        # 自己在diff函数中replace
        self.logger.info("开始运行差异蛋白分析")
        opts = {
            'group': self.option('protein_group'),
            'cmp': self.option('protein_control'),
            'ratio_exp': self.preprocess.option('preprocess_exp').prop['path'],
            'method_type': self.option('method_type'),
            'correct_method': self.option('correct_method'),
            'pvalue': self.option('pvalue'),
            'fc_up': self.option('fc_up'),
            'fc_down': self.option('fc_down'),
            'cutoffs': self.option('cutoffs'),
            'log':self.option('log'),
            'fill_type':self.option('fill_type'),
        }
        self.diff_pep.set_options(opts)
        self.diff_pep.on('end', self.run_diff_cluster)
        self.diff_pep.on('end', self.set_output, 'diff_pep')
        self.diff_pep.on('start', self.set_step, {'start': self.step.diff_pep})
        self.diff_pep.on('end', self.set_step, {'end': self.step.diff_pep})
        self.diff_pep.run()

    def run_diamond(self):
        self.logger.info("开始运行diamond比对")
        self.blast_modules = []
        blast_opts = {
            'query': self.option("protein_fasta"),
            'query_type': 'prot',
            'database': None,
            'blast': 'blastp',
            'evalue': None,
            'outfmt': 5,
            "diamond_version": "v2.0.4",
        }
        if self.option("annot_group") in ["REFRNA_GROUP_202110"]:
            blast_opts.update({
                "diamond_version": "v2.0.13"
            })
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            self.diamond_nr = self.add_module("labelfree.diamond")
            blast_opts.update(
                {
                    'database': self.option("nr_database"),
                    'evalue': self.option('go_evalue'),
                    'identity': self.option('go_identity'),
                    'nr_version' :self.annot_config_dict['nr']['version'],
                    'version' :self.annot_config_dict['diamond_all']['version']

                }
            )
            self.diamond_nr.set_options(blast_opts)
            self.blast_modules.append(self.diamond_nr)
            self.diamond_nr.on('end', self.set_output, 'diamond_nr')
            self.diamond_nr.run()
        if 'cog' in self.option('database'):
            self.diamond_eggnog = self.add_module("labelfree.diamond")
            blast_opts.update(
                {
                    'database': 'eggnog',
                    'evalue': self.option('cog_evalue'),
                    'identity': self.option('cog_identity'),
                    'eggnog_version' :self.annot_config_dict['eggnog']['version'],
                    'version' :self.annot_config_dict['diamond_all']['version']
                }
            )
            self.diamond_eggnog.set_options(blast_opts)
            self.blast_modules.append(self.diamond_eggnog)
            self.diamond_eggnog.on('end', self.set_output, 'diamond_eggnog')
            self.diamond_eggnog.run()

            self.diamond_string = self.add_module("labelfree.diamond")
            blast_opts.update(
                {
                    'database': 'string',
                    'evalue': self.option('cog_evalue'),
                    'identity': self.option('cog_identity'),
                    'string_version' :self.annot_config_dict['string']['version'],
                    'version' :self.annot_config_dict['diamond_all']['version']

                }
            )
            self.diamond_string.set_options(blast_opts)
            self.blast_modules.append(self.diamond_string)
            self.diamond_string.on('end', self.set_output, 'diamond_string')
            self.diamond_string.run()
        if 'kegg' in self.option('database'):
            self.diamond_kegg = self.add_module("labelfree.diamond")
            blast_opts.update(
                {
                    'database': 'kegg',
                    'evalue': self.option('kegg_evalue'),
                    'kegg_species' : self.option('kegg_org'),
                    'identity': self.option('kegg_identity'),
                    'kegg_version' :self.annot_config_dict['kegg']['version'],
                    'version' :self.annot_config_dict['diamond_all']['version']
                 }
            )
            self.diamond_kegg.set_options(blast_opts)
            self.blast_modules.append(self.diamond_kegg)
            self.diamond_kegg.on('end', self.set_output, 'diamond_kegg')
            self.diamond_kegg.run()
        if 'pfam' in self.option("database"):
            self.pfam = self.add_tool("labelfree.annotation.pfam")
            opts = {
                "fasta": self.option('protein_fasta'),
                "e_value": self.option('pfam_evalue'),
                "pfam_version": self.annot_config_dict['pfam']['version'],
            }
            self.pfam.set_options(opts)
            self.pfam.on("end", self.set_output, "pfam")
            self.blast_modules.append(self.pfam)
            self.pfam.run()
        self.on_rely([self.diamond_nr, self.diamond_kegg, self.pfam, self.diamond_eggnog, self.diamond_string], self.run_annotation)

    def run_annotation(self):
        self.logger.info("开始运行注释统计")
        anno_opts = {
            "des" : self.option("protein").prop['path'],
            "fa" : self.option('protein_fasta').prop['path'],
            "go_annot" : True,
            "nr_annot" : True,
            'kegg_species' : self.option('kegg_org'),
            "taxonomy" : self.option("kegg_class"),
            "blast_nr_xml" : self.diamond_nr.option('outxml'),
            "blast_kegg_xml" : self.diamond_kegg.option('outxml'),
            "blast_string_xml" : self.diamond_string.option('outxml'),
            "blast_eggnog_xml" : self.diamond_eggnog.option('outxml'),
            "pfam_domain" : self.pfam.output_dir + "/pfam_domain",
            "kegg_version" : self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "go_version" : self.annot_config_dict['go']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],
            'gram': self.option('gram'),
            'sub_loc': self.option('sub_loc'),
        }
        if self.option("annot_group") == "GROUP_2017":
            anno_opts.update({"cog_type": "string"})

        self.annotation.set_options(anno_opts)
        self.annotation.on('end', self.set_output, 'annotation')
        self.annotation.on('start', self.set_step, {'start': self.step.annotation})
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        class_prefix = [self.annotation, self.diff_pep]
        for func in [self.run_diff_cog_class, self.run_diff_go_class, self.run_diff_kegg_class, self.run_diff_pfam_stat, self.run_diff_subloc_stat,
                     self.run_diff_go_enrich, self.run_diff_kegg_enrich, self.run_diff_ipath, self.run_diff_ppi, self.run_diff_string_picture]:
            # self.on_rely(class_prefix, func)
            self.annotation.on('end', func)
            gevent.sleep(1)
        self.annotation.run()

    def run_diff_cluster(self):
        self.logger.info("开始运行差异蛋白cog分类注释")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "exp" : self.preprocess.option('preprocess_exp').prop['path'],
            "group" : self.option('protein_group'),
            "scd" : "euclidean"
        }
        self.diff_cluster.set_options(opts)
        self.diff_cluster.on('end', self.set_output, 'diff_cluster')
        self.diff_cluster.on('start', self.set_step, {'start': self.step.diff_cluster})
        self.diff_cluster.on('end', self.set_step, {'end': self.step.diff_cluster})
        self.diff_cluster.run()

    def run_diff_cog_class(self):
        self.logger.info("开始运行差异蛋白cog分类注释")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "cog_stat" : os.path.join(self.annotation.output_dir, 'cog', 'cog_summary.xls'),
        }
        self.diff_cog_class.set_options(opts)
        self.diff_cog_class.on('end', self.set_output, 'diff_cog_class')
        self.diff_cog_class.on('start', self.set_step, {'start': self.step.diff_cog_class})
        self.diff_cog_class.on('end', self.set_step, {'end': self.step.diff_cog_class})
        try:
            self.diff_cog_class.run()
        except Exception:
            try:
                self.diff_cog_class.end()
            except Exception:
                pass

    def run_diff_go_class(self):
        self.logger.info("开始运行差异蛋白go分类注释")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "go_version" : self.annot_config_dict['go']['version'],
            "go_stat" : os.path.join(self.annotation.output_dir, 'go', 'go12level_statistics.xls'),
        }
        self.diff_go_class.set_options(opts)
        self.diff_go_class.on('end', self.set_output, 'diff_go_class')
        self.diff_go_class.on('start', self.set_step, {'start': self.step.diff_go_class})
        self.diff_go_class.on('end', self.set_step, {'end': self.step.diff_go_class})
        try:
            self.diff_go_class.run()
        except Exception:
            # self.diff_go_class.end()
            try:
                self.diff_go_class.end()
            except Exception:
                pass

    def run_diff_kegg_class(self):
        self.logger.info("开始运行差异蛋白kegg分类注释")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "png_dir" : os.path.join(self.annotation.output_dir, 'kegg','pathways'),
            "kegg_table" : os.path.join(self.annotation.output_dir, 'kegg','kegg_table.xls'),
            "pathway_table" : os.path.join(self.annotation.output_dir, 'kegg','pathway_table.xls'),
            'kegg_version' :self.annot_config_dict['kegg']['version'],
        }
        self.diff_kegg_class.set_options(opts)
        self.diff_kegg_class.on('end', self.set_output, 'diff_kegg_class')
        self.diff_kegg_class.on('start', self.set_step, {'start': self.step.diff_kegg_class})
        self.diff_kegg_class.on('end', self.set_step, {'end': self.step.diff_kegg_class})
        try:
            self.diff_kegg_class.run()
        except Exception:
            # self.diff_kegg_class.end()
            try:
                self.diff_kegg_class.end()
            except Exception:
                pass

    def run_diff_pfam_stat(self):
        self.logger.info("开始运行差异蛋白pfam分类注释")
        opts = {
            "diff_path": self.diff_pep.output_dir,
            "pfam_file": os.path.join(self.annotation.output_dir, 'blast_xml', 'pfam_domain'),
        }
        self.diff_pfam_stat.set_options(opts)
        self.diff_pfam_stat.on('end', self.set_output, 'diff_pfam_stat')
        self.diff_pfam_stat.on('start', self.set_step, {'start': self.step.diff_pfam_stat})
        self.diff_pfam_stat.on('end', self.set_step, {'end': self.step.diff_pfam_stat})
        try:
            self.diff_pfam_stat.run()
        except Exception:
            # self.diff_pfam_stat.end()
            try:
                self.diff_pfam_stat.end()
            except Exception:
                pass

    def run_diff_subloc_stat(self):
        self.logger.info("开始运行差异蛋白subloc分类注释")
        opts = {
            "diff_path": self.diff_pep.output_dir,
            "subloc_file": os.path.join(self.annotation.output_dir, 'subloc', 'multiloc.xls'),
        }
        self.diff_subloc_stat.set_options(opts)
        self.diff_subloc_stat.on('end', self.set_output, 'diff_subloc_stat')
        self.diff_subloc_stat.on('start', self.set_step, {'start': self.step.diff_subloc_stat})
        self.diff_subloc_stat.on('end', self.set_step, {'end': self.step.diff_subloc_stat})
        try:
            self.diff_subloc_stat.run()
        except Exception:
            # self.diff_subloc_stat.end()
            try:
                self.diff_subloc_stat.end()
            except Exception:
                pass

    def run_diff_ipath(self):
        self.logger.info("开始运行差异蛋白ipath分析")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "kegg_table" : os.path.join(self.annotation.output_dir, 'kegg', 'kegg_table.xls'),
        }
        self.diff_ipath.set_options(opts)
        self.diff_ipath.on('end', self.set_output, 'diff_ipath')
        self.diff_ipath.on('start', self.set_step, {'start': self.step.diff_ipath})
        self.diff_ipath.on('end', self.set_step, {'end': self.step.diff_ipath})
        try:
            self.diff_ipath.run()
        except Exception:
            # self.diff_ipath.end()
            try:
                self.diff_ipath.end()
            except Exception:
                pass

    def run_diff_ppi(self):
        self.logger.info("开始运行差异蛋白ppi分析")
        self.logger.info(self.option('ppi_species'))
        opts = {
            "diff_path": self.diff_pep.output_dir,
            "seq": self.option('protein_fasta'),
            "species": self.option('ppi_species'),
            "string_xml": os.path.join(self.annotation.output_dir, 'blast_xml', 'string.xml'),
        }
        self.diff_ppi.set_options(opts)
        self.diff_ppi.on('end', self.set_output, 'diff_ppi')
        self.diff_ppi.on('start', self.set_step, {'start': self.step.diff_ppi})
        self.diff_ppi.on('end', self.set_step, {'end': self.step.diff_ppi})
        try:
            self.diff_ppi.run()
        except Exception:
            # self.diff_ppi.end()
            try:
                self.diff_ppi.end()
            except Exception:
                pass

    def run_diff_string_picture(self):
        self.logger.info("开始运行差异蛋白string数据库图片等的爬取")
        useblast = 'no'
        if self.option('useblast'):
            useblast = 'yes'
        opts = {
            "diff_path": self.diff_pep.output_dir,
            "useblast": useblast,
            "species": self.option('ppi_species'),
            "string_xml": os.path.join(self.annotation.output_dir, 'blast_xml', 'string.xml'),
        }
        self.diff_string_picture.set_options(opts)
        self.diff_string_picture.on('end', self.set_output, 'diff_string_picture')
        self.diff_string_picture.on('start', self.set_step, {'start': self.step.diff_string_picture})
        self.diff_string_picture.on('end', self.set_step, {'end': self.step.diff_string_picture})
        try:
            self.diff_string_picture.run()
        except Exception:
            # self.diff_string_picture.end()
            try:
                self.diff_string_picture.end()
            except Exception:
                pass

    def run_diff_kegg_enrich(self):
        self.logger.info("开始运行差异蛋白kegg富集分析")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "kegg_table" : os.path.join(self.annotation.output_dir, 'kegg', 'kegg_table.xls'),
            "png_dir" : os.path.join(self.annotation.output_dir, 'kegg','pathways'),
            "pathway_table" : os.path.join(self.annotation.output_dir, 'kegg', 'pathway_table.xls'),
            'kegg_version' :self.annot_config_dict['kegg']['version'],
        }
        self.diff_kegg_enrich.set_options(opts)
        self.diff_kegg_enrich.on('end', self.set_output, 'diff_kegg_enrich')
        self.diff_kegg_enrich.on('start', self.set_step, {'start': self.step.diff_kegg_enrich})
        self.diff_kegg_enrich.on('end', self.set_step, {'end': self.step.diff_kegg_enrich})
        try:
            self.diff_kegg_enrich.run()
        except Exception:
            # self.diff_kegg_enrich.end()
            try:
                self.diff_kegg_enrich.end()
            except Exception:
                pass

    def run_diff_go_enrich(self):
        self.logger.info("开始运行差异蛋白go富集分析")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "go_version" : self.annot_config_dict['go']['version'],
            "go_list" : os.path.join(self.annotation.output_dir, 'go', 'query_gos.list'),
        }
        self.diff_go_enrich.set_options(opts)
        self.diff_go_enrich.on('end', self.set_output, 'diff_go_enrich')
        self.diff_go_enrich.on('start', self.set_step, {'start': self.step.diff_go_enrich})
        self.diff_go_enrich.on('end', self.set_step, {'end': self.step.diff_go_enrich})
        try:
            self.diff_go_enrich.run()
        except Exception:
            # self.diff_go_enrich.end()
            try:
                self.diff_go_enrich.end()
            except Exception:
                pass

    def run_diff_circle(self):
        self.logger.info("开始运行差异蛋白go富集分析")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "kegg_enrich_path" : self.diff_kegg_enrich.output_dir,
            "go_enrich_path" : self.diff_go_enrich.output_dir,
        }
        self.diff_circle.set_options(opts)
        self.diff_circle.on('end', self.set_output, 'diff_circle')
        self.diff_circle.on('start', self.set_step, {'start': self.step.diff_circle})
        self.diff_circle.on('end', self.set_step, {'end': self.step.diff_circle})
        try:
            self.diff_circle.run()
        except Exception:
            # self.diff_circle.end()
            try:
                self.diff_circle.end()
            except Exception:
                pass

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
        newdir = os.path.join(self.output_dir, newname)
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}到{}移动耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        if os.path.isfile(old_file):
            os.link(old_file, new_file)
        else:
            os.mkdir(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == 'searchdb':
            self.move2outputdir(obj.output_dir, 'SearchDB')
        if event['data'] == 'sam_pca':
            self.move2outputdir(obj.output_dir, 'SamPca')
        if event['data'] == 'sam_corr':
            self.move2outputdir(obj.output_dir, 'SamCorr')
        if event['data'] == 'annotation':
            self.move2outputdir(obj.output_dir, 'Annotation')
        if event['data'] == 'diff_pep':
            self.move2outputdir(obj.output_dir, 'DiffPep')
        if event['data'] == 'ellipse':
            self.move2outputdir(obj.output_dir, 'Ellipse')
        if event['data'] == 'exp_venn':
            self.move2outputdir(obj.output_dir, 'ExpVenn')
        if event['data'] == 'diff_cluster':
            self.move2outputdir(obj.output_dir, '5_Proteinset/01_PsetCluster')
        if event['data'] == 'diff_go_class':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_PsetAnno/01_PsetGO')
        if event['data'] == 'diff_kegg_class':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_PsetAnno/02_PsetKEGG')
            kegg_dir = os.path.join(self.output_dir, '5_Proteinset/03_PsetAnno/02_PsetKEGG')
            for sub_dir in os.listdir(kegg_dir):
                dir_to_tar = os.path.join(kegg_dir, sub_dir, 'pathways')
                target_tar = os.path.join(kegg_dir, sub_dir, 'pathways.tar')
                if os.path.exists(dir_to_tar):
                    self.make_tar(target_tar,dir_to_tar)
                if os.path.exists(target_tar) and os.path.exists(dir_to_tar):
                    shutil.rmtree(dir_to_tar)
                dir_to_tar2 = os.path.join(kegg_dir, sub_dir, 'ko')
                target_tar2 = os.path.join(kegg_dir, sub_dir, 'ko.tar')
                if os.path.exists(dir_to_tar2):
                    self.make_tar(target_tar2,dir_to_tar2)
                if os.path.exists(target_tar2) and os.path.exists(dir_to_tar2):
                    shutil.rmtree(dir_to_tar2)
        if event['data'] == 'diff_cog_class':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_PsetAnno/03_PsetCOG')
        if event['data'] == 'diff_pfam_stat':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_PsetAnno/04_PsetPfam')
        if event['data'] == 'diff_subloc_stat':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_PsetAnno/05_PsetSubLoc')
        if event['data'] == 'diff_go_enrich':
            self.move2outputdir(obj.output_dir, '5_Proteinset/04_PsetEnrich/01_EnrichGO')
        if event['data'] == 'diff_kegg_enrich':
            self.move2outputdir(obj.output_dir, '5_Proteinset/04_PsetEnrich/02_EnrichKEGG')
            kegg_dir = os.path.join(self.output_dir, '5_Proteinset/04_PsetEnrich/02_EnrichKEGG')
            for sub_dir in os.listdir(kegg_dir):
                dir_to_tar = os.path.join(kegg_dir, sub_dir, 'pathways')
                target_tar = os.path.join(kegg_dir, sub_dir, 'pathways.tar')
                if os.path.exists(dir_to_tar):
                    self.make_tar(target_tar, dir_to_tar)
                if os.path.exists(target_tar) and os.path.exists(dir_to_tar):
                    shutil.rmtree(dir_to_tar)
        if event['data'] == 'diff_circle':
            self.move2outputdir(obj.output_dir, '5_Proteinset/04_PsetEnrich/03_EnrichChord')
            circ_dir = os.path.join(self.output_dir, '5_Proteinset/04_PsetEnrich/03_EnrichChord')
            for each in os.listdir(circ_dir):
                old_dir = os.path.join(circ_dir, each)
                new_dir = os.path.join(circ_dir, each.split('circle_')[1])
                if not os.path.exists(new_dir):
                    os.mkdir(new_dir)
                for i in os.listdir(old_dir):
                    if i == 'enrich_zscore':
                        j = each.split('circle_')[0] + i
                    else:
                        j = i
                    os.link(os.path.join(old_dir, i), os.path.join(new_dir, j))
                shutil.rmtree(old_dir)
        if event['data'] == 'diff_ppi':
            self.move2outputdir(obj.output_dir, '5_Proteinset/05_PsetPPI')
        if event['data'] == 'diff_string_picture':
            self.move2outputdir(obj.output_dir, '5_Proteinset/05_PsetStringPic')
        if event['data'] == 'diff_ipath':
            self.move2outputdir(obj.output_dir, '5_Proteinset/06_PsetIpath')

    def end(self):
        self.build_seq_database() # 创建序列数据库
        self.run_api() # 运行导表函数
        self.add_anno_to_diff() # 给差异表达文件增加注释
        # self.run_chart() # 运行chart生成静态图，放在这里是因为qc文件必须在run_api()之后

        ## 导表后，修改文件的绝对路径
        db = Config().get_mongo_client(mtype="dia")[Config().get_mongo_dbname("dia")]
        col1 = db["sg_annotation_stat"]
        col1.update({"task_id" : self.task_id}, {"$set": {"result_dir": os.path.join(self.workflow_output, "Annotation")}}, upsert=True)
        col2 = db["sg_task"]
        col2.update({"task_id" : self.task_id}, {"$set": {"seq_db": os.path.join(self.workflow_output, "SequenceDatabase/seq_db.sqlite3")}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"protein_sliced": os.path.join(self.workflow_output, "SearchDB/protein_sliced.xls")}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"protein_fa": os.path.join(self.workflow_output, "ProteinFasta/protein.fa")}}, upsert=True)
        # 声明是V3
        col2.update({"task_id" : self.task_id}, {"$set": {"version": "V3.1"}}, upsert=True)
        annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
        if self.annot_config_dict['kegg']['version'] > "2020":
            if self.option("kegg_org") not in [None, ""]:
                annot_version_dict['kegg'] += "_spe"
        else:
            del annot_version_dict['kegg']
        col2.update({'task_id': self.task_id},
                    {'$set': {'database_version': annot_version_dict,
                              'annot_group': self.option("annot_group")}}, upsert=True)
        self.modify_output() # 修改文件目录结构
        super(DiaWorkflow, self).end()

    def make_tar(self, output_filename, source_dir):
        with tarfile.open(output_filename, "w") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))

    def modify_output(self):
        target_dir = self.output_dir
        if os.path.exists(target_dir + "/SequenceDatabase"):
            shutil.rmtree(target_dir + "/SequenceDatabase")
        os.mkdir(target_dir + "/SequenceDatabase")
        seq_db = self.work_dir  + "/seq_db.sqlite3"
        os.link(seq_db, target_dir + "/SequenceDatabase/seq_db.sqlite3")
        if os.path.exists(target_dir + "/ProteinFasta"):
            shutil.rmtree(target_dir + "/ProteinFasta")
        os.mkdir(target_dir + "/ProteinFasta")
        fa_file = self.option('protein_fasta').prop["path"]
        os.link(fa_file, target_dir + "/ProteinFasta/protein.fa")
        if os.path.exists(target_dir + "/SearchDB/protein_sliced.xls"):
            os.remove(target_dir + "/SearchDB/protein_sliced.xls")
        protein_splice = self.searchdb.work_dir + "/protein_sliced.xls"
        os.link(protein_splice, target_dir + "/SearchDB/protein_sliced.xls")
        if os.path.exists(target_dir + "/Express"):
            shutil.rmtree(target_dir + "/Express")
        os.mkdir(target_dir + "/Express")
        ratio_exp = os.path.join(self.work_dir, "treat_ref")
        os.link(ratio_exp, target_dir + "/Express/exp.xls")
        sam_pca_dir = target_dir + "/SamPca"
        if os.path.exists(sam_pca_dir + "/pca_rotation_all.xls"):
            os.remove(sam_pca_dir + "/pca_rotation.xls")
            os.rename(sam_pca_dir + "/pca_rotation_all.xls", sam_pca_dir + "/pca_rotation.xls")
        if os.path.exists(target_dir + "/Annotation/anno_stat/blast"):
            shutil.rmtree(target_dir + "/Annotation/anno_stat/blast")

        # V2版文件结构调整，多生成了一个交付给客户的结果文件目录
        for dir in ['1_DataInfo', '2_SampleComp', '3_Annotation', '4_ExpDiff']:
            dir = os.path.join(target_dir, dir)
            if os.path.exists(dir):
                shutil.rmtree(dir)
            os.makedirs(dir)

        tmp = os.path.join(target_dir, '1_DataInfo', '01_SearchDB')
        os.makedirs(tmp)
        for file in glob.glob(os.path.join(target_dir, 'SearchDB/*xls')):
            os.link(file, os.path.join(tmp, os.path.basename(file)))
        tmp = os.path.join(target_dir, '1_DataInfo', '02_QcInfo')
        os.makedirs(tmp)
        for file in glob.glob(os.path.join(self.work_dir, 'qc/*xls')):
            os.link(file, os.path.join(tmp, os.path.basename(file)))
        tmp = os.path.join(target_dir, '1_DataInfo', '03_ExpInfo')
        os.makedirs(tmp)
        for file in glob.glob(os.path.join(target_dir, 'Express/*xls')):
            os.link(file, os.path.join(tmp, os.path.basename(file)))

        tmp = os.path.join(target_dir, '2_SampleComp', '01_SamPca')
        os.makedirs(tmp)
        for file in glob.glob(os.path.join(target_dir, 'SamPca/*xls')):
            os.link(file, os.path.join(tmp, os.path.basename(file)))
        tmp = os.path.join(target_dir, '2_SampleComp', '02_SamCorr')
        os.makedirs(tmp)
        for file in glob.glob(os.path.join(target_dir, 'SamCorr/*xls')):
            os.link(file, os.path.join(tmp, os.path.basename(file)))
        tmp = os.path.join(target_dir, '2_SampleComp', '03_ExpVenn')
        os.makedirs(tmp)
        for file in glob.glob(os.path.join(target_dir, 'ExpVenn/*xls')):
            os.link(file, os.path.join(tmp, os.path.basename(file)))

        tmp = os.path.join(target_dir, '3_Annotation')
        dir_to_tar = os.path.join(target_dir, 'Annotation', 'kegg','pathways')
        target_tar = os.path.join(target_dir, 'Annotation', 'kegg','pathways.tar')
        if os.path.exists(dir_to_tar):
            self.make_tar(target_tar,dir_to_tar)
        if os.path.exists(target_tar) and os.path.exists(dir_to_tar):
            shutil.rmtree(dir_to_tar)
        # os.makedirs(tmp)
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'anno_stat'), os.path.join(tmp, '01_Stat'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'go'), os.path.join(tmp, '02_GO'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'kegg'), os.path.join(tmp, '03_KEGG'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'cog'), os.path.join(tmp, '04_COG'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'subloc'), os.path.join(tmp, '06_SubLoc'))
        tmp = os.path.join(target_dir, '3_Annotation', '05_Pfam')
        os.makedirs(tmp)
        os.link(os.path.join(target_dir, 'Annotation', 'blast_xml', 'pfam_domain'), os.path.join(tmp, 'pfam_domain'))

        tmp = os.path.join(target_dir, '4_ExpDiff')
        # os.makedirs(tmp)
        diff_path = os.path.join(target_dir, 'DiffPep')
        for file in os.listdir(diff_path):
            if file.endswith('_diff.xls'):
                cmp = file.split('_diff.xls')[0]
                os.makedirs(os.path.join(tmp, cmp))
                for file_ in glob.glob(os.path.join(diff_path, cmp+'*')):
                    os.link(file_, os.path.join(tmp, cmp, os.path.basename(file_)))
            try:
                os.link(os.path.join(diff_path, 'diff_protein_summary.xls'), os.path.join(tmp, 'diff_protein_summary.xls'))
            except:
                pass
        for dir in ['Express', 'SamCorr', 'SamPca', 'DiffPep', 'ExpVenn', 'Ellipse']:
            try:
                shutil.rmtree(os.path.join(target_dir, dir))
            except:
                pass

        # #下面的fat_dir中的目录的文件个数太多，会导致接口传输失败不能上传结果目录，需要打包减少文件个数。
        # fat_dir=os.path.join(self.output_dir,'5_Proteinset','03_PsetAnno','02_PsetKEGG')
        # for dir_name in os.listdir(fat_dir):
        #     dir_to_tar = os.path.join(fat_dir,dir_name,'pathways')
        #     target_tar = os.path.join(fat_dir,dir_name,'pathways.tar')
        #     if os.path.exists(dir_to_tar):
        #         self.make_tar(target_tar,dir_to_tar)
        #     if os.path.exists(target_tar):
        #         shutil.rmtree(dir_to_tar)
        #
        # for dir_name3 in os.listdir(fat_dir):
        #     dir_to_tar3 = os.path.join(fat_dir,dir_name3,'ko')
        #     target_tar3 = os.path.join(fat_dir,dir_name3,'ko.tar')
        #     if os.path.exists(dir_to_tar3):
        #         self.make_tar(target_tar3,dir_to_tar3)
        #     if os.path.exists(target_tar3):
        #         shutil.rmtree(dir_to_tar3)
        #
        # fat_dir4=os.path.join(self.output_dir,'5_Proteinset','04_PsetEnrich','02_EnrichKEGG')
        # for dir_name4 in os.listdir(fat_dir4):
        #     dir_to_tar4 = os.path.join(fat_dir4,dir_name4,'pathways')
        #     target_tar4 = os.path.join(fat_dir4,dir_name4,'pathways.tar')
        #     if os.path.exists(dir_to_tar4):
        #         self.make_tar(target_tar4,dir_to_tar4)
        #     if os.path.exists(target_tar4):
        #         shutil.rmtree(dir_to_tar4)
        #
        # fat_dir2=os.path.join(self.output_dir,'3_Annotation','03_KEGG')
        # dir_to_tar2 = os.path.join(fat_dir2,'pathways')
        # target_tar2 = os.path.join(fat_dir2,'pathways.tar')
        # if os.path.exists(dir_to_tar2):
        #     self.make_tar(target_tar2,dir_to_tar2)
        # if os.path.exists(target_tar2):
        #     shutil.rmtree(dir_to_tar2)

        ### 增加静态图文件
        all_cmp_names_list = [i[1]+"_vs_"+i[0] for i in self.option("protein_control").prop["cmp_list"]]

        def os_link(a,b):
            if os.path.exists(b):
                os.remove(b)
            os.link(a,b)
        tmp_work_dir_name = time.strftime('%Y%m%d_%H%M%S', time.localtime(time.time())) +"_tmp_pdf_dir_by_mainworkflow"
        os.makedirs(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name))
        for pdf_file in glob.glob(os.path.join(self.work_dir, 'Chart/output/*.pdf')):
            os.system('cp {} {}'.format(pdf_file, os.path.join(os.path.dirname(pdf_file), tmp_work_dir_name, ".".join(os.path.basename(pdf_file).split('.')[:-2])+'.pdf')))

        pdf_files = [os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, i) for i in ["pep_len.pdf", "info.pdf","mw.pdf"]]
        pdf_files += glob.glob(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, '*inter-cv.pdf'))
        pdf_files += glob.glob(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, '*intra-cv.pdf'))
        pdf_files += glob.glob(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, '*assessment.pdf'))
        for each in pdf_files:
            if os.path.exists(each):
                os_link(each, os.path.join(self.output_dir, '1_DataInfo', '02_QcInfo', os.path.basename(each)))
        # #
        # if os.path.exists(os.path.join(self.plot_static_peperr.output_dir, "dMass.pdf")):
        #     if os.path.exists(os.path.join(self.output_dir, '1_DataInfo', '02_QcInfo', "pep_error.pdf")):
        #         os.remove(os.path.join(self.output_dir, '1_DataInfo', '02_QcInfo', "pep_error.pdf"))
        #     os_link(os.path.join(self.plot_static_peperr.output_dir, "dMass.pdf"), os.path.join(self.output_dir, '1_DataInfo', '02_QcInfo', "pep_error.pdf"))
        # if os.path.exists(os.path.join(self.plot_static_peperr.output_dir, "dMass.png")):
        #     os_link(os.path.join(self.plot_static_peperr.output_dir, "dMass.png"), os.path.join(self.output_dir, '1_DataInfo', '02_QcInfo', "pep_error.png"))
        # #
        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "sam_pca.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '2_SampleComp', '01_SamPca', "sam_pca.pdf"))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "exp.heatmap.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '2_SampleComp', '02_SamCorr', "exp.heatmap.pdf"))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "exp.tree.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '2_SampleComp', '02_SamCorr', "exp.tree.pdf"))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "exp_venn.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '2_SampleComp', '03_ExpVenn', "exp_venn.pdf"))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "anno_stat.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '3_Annotation', '01_Stat', "anno_stat.pdf"))

        pdf_files = glob.glob(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, 'go_lev*.pdf'))
        for each in pdf_files:
            if os.path.exists(each):
                os_link(each, os.path.join(self.output_dir, '3_Annotation', '02_GO', os.path.basename(each)))

        for i in ["path_class.pdf","key_path.pdf"]:
            pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, i)
            if os.path.exists(pdf_file):
                os_link(pdf_file, os.path.join(self.output_dir, '3_Annotation', '03_KEGG', i))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "cog_bar.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '3_Annotation', '04_COG', "cog_bar.pdf"))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "pfam_bar.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '3_Annotation', '05_Pfam', "pfam_bar.pdf"))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "subloc_bar.pdf")
        if os.path.exists(pdf_file):
            os_link(pdf_file, os.path.join(self.output_dir, '3_Annotation', '06_SubLoc', "subloc_bar.pdf"))

        pdf_files = glob.glob(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, '*_volcano.pdf'))
        pdf_files += glob.glob(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, '*_scatter.pdf'))
        for each in pdf_files:
            cmp = os.path.basename(each).split('_volcano')[0].split('_scatter')[0]
            if not os.path.exists(os.path.join(self.output_dir, '4_ExpDiff', cmp)):
                os.makedirs(os.path.join(self.output_dir, '4_ExpDiff', cmp))
            os_link(each, os.path.join(self.output_dir, '4_ExpDiff', cmp, os.path.basename(each)))

        for i in ["cluster_group_","cluster_sample_"]:
            for file in glob.glob(os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, i+'*.pdf')):
                one_pdf = os.path.basename(file)
                if "____" in one_pdf:
                    if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '01_PsetCluster', one_pdf.split('____')[0])):
                        os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '01_PsetCluster', one_pdf.split('____')[0]))
                    os_link(file, os.path.join(self.output_dir, '5_Proteinset', '01_PsetCluster', one_pdf.split('____')[0], one_pdf.split('____')[1]))

        pdf_file = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, "venn.pdf")
        if os.path.exists(pdf_file):
            if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '02_PsetVenn')):
                os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '02_PsetVenn'))
            os_link(pdf_file, os.path.join(self.output_dir, '5_Proteinset', '02_PsetVenn', "venn.pdf"))

        for i in ['_all_go.pdf','_up_down_go.pdf','_up_down_go_neg.pdf', '_up_go.pdf', '_down_go.pdf']:
            for comp in all_cmp_names_list:
                one_pdf = comp+'____'+comp+i
                one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                if os.path.exists(one_pdf_full):
                    if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','01_PsetGO', comp)):
                        os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','01_PsetGO', comp))
                    os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno', '01_PsetGO', comp, comp+i))

        for comp in all_cmp_names_list:
            one_pdf = comp+'_all'+'____'+comp+"_path_all.pdf"
            one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
            if os.path.exists(one_pdf_full):
                if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','02_PsetKEGG', comp+'_all')):
                    os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','02_PsetKEGG', comp+'_all'))
                os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno', '02_PsetKEGG', comp+'_all', comp+'_path_all.pdf'))
            two_pdf = comp+'_up_down'+'____'+comp+"_path_down.pdf"
            two_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, two_pdf)
            if os.path.exists(two_pdf_full):
                if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','02_PsetKEGG', comp+'_up_down')):
                    os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','02_PsetKEGG', comp+'_up_down'))
                os_link(two_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno', '02_PsetKEGG', comp+'_up_down', comp+'_path_down.pdf'))
            three_pdf = comp+'_up_down'+'____'+comp+"_path_up.pdf"
            three_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, three_pdf)
            if os.path.exists(three_pdf_full):
                if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','02_PsetKEGG', comp+'_up_down')):
                    os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','02_PsetKEGG', comp+'_up_down'))
                os_link(three_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno', '02_PsetKEGG', comp+'_up_down', comp+'_path_up.pdf'))

        for i in ['_all_cog.pdf','_up_down_cog.pdf']:
            for comp in all_cmp_names_list:
                one_pdf = comp+'____'+comp+i
                one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                if os.path.exists(one_pdf_full):
                    if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','03_PsetCOG', comp)):
                        os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','03_PsetCOG', comp))
                    os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno', '03_PsetCOG', comp, comp+i))

        for i in ['_all', '_up_down']:
            for comp in all_cmp_names_list:
                one_pdf = comp+'____'+comp+i + '_pfam.pdf'
                one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                if os.path.exists(one_pdf_full):
                    if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','04_PsetPfam', comp+i)):
                        os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','04_PsetPfam', comp+i))
                    os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno', '04_PsetPfam', comp+i, comp+i+'_pfam.pdf'))

        for i in ['_all','_up_down']:
            for comp in all_cmp_names_list:
                one_pdf = comp+'____'+comp+i+'_subloc.pdf'
                one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                if os.path.exists(one_pdf_full):
                    if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','05_PsetSubLoc', comp+i)):
                        os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno','05_PsetSubLoc', comp+i))
                    os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '03_PsetAnno', '05_PsetSubLoc', comp+i, comp+i+'_subloc.pdf'))

        for i in ['_go_enrich_bar.pdf','_go_enrich_bubble.pdf','_go_enrich_bubble2.pdf']:
            for comp in all_cmp_names_list:
                for dir_suf in ['_all','_up','_down']:
                    one_pdf = comp+dir_suf+'____'+comp+i
                    one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                    if os.path.exists(one_pdf_full):
                        if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich','01_EnrichGO', comp+dir_suf)):
                            os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich','01_EnrichGO', comp+dir_suf))
                        os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich', '01_EnrichGO', comp+dir_suf, comp+i))

        for i in ['_kegg_enrich_bar.pdf','_kegg_enrich_bubble.pdf','_kegg_enrich_bubble2.pdf']:
            for comp in all_cmp_names_list:
                for dir_suf in ['_all','_up','_down']:
                    one_pdf = comp+dir_suf+'____'+comp+i
                    one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                    if os.path.exists(one_pdf_full):
                        if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich','02_EnrichKEGG', comp+dir_suf)):
                            os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich','02_EnrichKEGG', comp+dir_suf))
                        os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich', '02_EnrichKEGG', comp+dir_suf, comp+i))

        for i in ['_kegg_chord.pdf','_go_chord.pdf']:
            for comp in all_cmp_names_list:
                for dir_suf in ['_all','_up','_down']:
                    one_pdf = comp+dir_suf+'____'+comp+i
                    one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                    if os.path.exists(one_pdf_full):
                        if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich','03_EnrichChord', comp+dir_suf)):
                            os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich','03_EnrichChord', comp+dir_suf))
                        os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '04_PsetEnrich', '03_EnrichChord', comp+dir_suf, comp+i))

        for i in ['ppi.centrality.line.pdf','ppi.degree.line.pdf']:
            for comp in all_cmp_names_list:
                for dir_suf in ['_all','_up','_down']:
                    one_pdf = comp+dir_suf+'____ppinetwork_topology____'+i
                    one_pdf_full = os.path.join(self.work_dir, 'Chart/output/', tmp_work_dir_name, one_pdf)
                    if os.path.exists(one_pdf_full):
                        if not os.path.exists(os.path.join(self.output_dir, '5_Proteinset', '05_PsetPPI', comp+dir_suf, 'ppinetwork_topology')):
                            os.makedirs(os.path.join(self.output_dir, '5_Proteinset', '05_PsetPPI', comp+dir_suf, 'ppinetwork_topology'))
                        os_link(one_pdf_full, os.path.join(self.output_dir, '5_Proteinset', '05_PsetPPI', comp+dir_suf, 'ppinetwork_topology', i))

        ####

        repaths = [
            [".", "", "流程分析结果目录"],
            ["SequenceDatabase", "", "蛋白序列数据库文件"],
            ["SequenceDatabase/seq_db.sqlite3", "", "", 1],
            ["ProteinFasta", "", "蛋白序列FASTA文件"],
            ["SearchDB", "", "搜库结果目录"],
            ["SearchDB/protein_sliced.xls", "", "搜库得到的蛋白详细信息"],
            ["SearchDB/search_db.xls", "", "搜库得到的蛋白详细信息以及在各个样本中的丰度"],
            ["SamPca", "", "样本PCA分析结果目录"],
            ["SamPca/pca_rotation.xls", "", "蛋白相关主成分贡献度"],
            ["SamPca/pca_sites.xls", "", "PCA分析样本坐标表"],
            ["SamPca/pca_importance.xls", "", "PCA主成分解释表"],
            ["SamCorr", "", "样本相关性分析结果目录"],
            ["SamCorr/sample_correlation.xls", "", "样本间相关性系数表"],
            ["SamCorr/expression_matrix.xls", "", "蛋白表达丰度矩阵表"],
            ["Annotation", "", "蛋白注释结果中间文件"],
            ["Annotation/anno_stat", "", "蛋白注释统计结果目录"],
            ["Annotation/anno_stat/proteins_anno_detail.xls", "", "鉴定蛋白附加6数据库详细信息的表"],
            ["Annotation/anno_stat/all_annotation_statistics.xls", "", "6个数据库中相关蛋白数以及总蛋白占比"],
            ["Annotation/anno_stat/venn", "", "鉴定结果中蛋白在各个数据库注释结果统计VENN图目录"],
            ["Annotation/anno_stat/venn/pfam_venn.txt", "", "有与pfam数据库相关信息的蛋白"],
            ["Annotation/anno_stat/venn/kegg_venn.txt", "", "有与kegg数据库相关信息的蛋白"],
            ["Annotation/anno_stat/venn/go_venn.txt", "", "有与go数据库相关信息的蛋白"],
            ["Annotation/anno_stat/venn/cog_venn.txt", "", "有与cog数据库相关信息的蛋白"],
            ["Annotation/anno_stat/venn/nr_venn.txt", "", "有与nr数据库相关信息的蛋白"],
            ["Annotation/anno_stat/venn/subloc_venn.txt", "", "有与subloc数据库相关信息的蛋白"],
            ["Annotation/go", "", "鉴定蛋白与GO数据库比对结果详情目录"],
            ["Annotation/go/blast2go.annot", "", "鉴定蛋白与GO的对应关系以及GOterm名称"],
            ["Annotation/go/query_gos.list", "", "鉴定蛋白与GO的对应关系"],
            ["Annotation/go/go12level_statistics.xls", "", "鉴定蛋白的GO1、2层级的详细信息"],
            ["Annotation/go/go123level_statistics.xls", "", "鉴定蛋白的GO1、2、3层级的详细信息"],
            ["Annotation/go/go1234level_statistics.xls", "", "鉴定蛋白的GO1、2、3、4层级的详细信息"],
            ["Annotation/cog", "", "鉴定蛋白与COG数据库比对结果详情目录"],
            ["Annotation/cog/cog_table.xls", "", "鉴定蛋白序列文件与COG比对结果表"],
            ["Annotation/cog/cog_summary.xls", "", "COG数据库注释结果汇总"],
            ["Annotation/cog/cog_list.xls", "", "鉴定蛋白与COG的对应关系"],
            ["Annotation/kegg", "", "鉴定蛋白与KEGG数据库比对结果详情目录"],
            ["Annotation/kegg/pathways", "", "KEGG通路注释结果图片目录"],
            ["Annotation/kegg/pathways/map*.png", "", "通路注释结果图片-位图"],
            ["Annotation/kegg/pathways/map*.pdf", "", "通路注释结果图片-矢量图"],
            ["Annotation/kegg/pathway_table.xls", "", "每行以通路为单位显示通路注释详情表"],
            ["Annotation/kegg/kegg_table.xls", "", "每行以蛋白为单位显示通路注释详情表"],
            ["Annotation/kegg/pid.txt", "", "每个注释通路和映射上去的蛋白联系表"],
            ["Annotation/subloc", "", "蛋白亚细胞定位注释结果目录"],
            ["Annotation/subloc/multiloc_stat.xls", "", "亚细胞位置对应蛋白的数量表"],
            ["Annotation/subloc/multiloc.xls", "", "亚细胞定位结果表"],
            ["Annotation/blast_xml", "", "蛋白比对注释结果目录", 1],
            ["Annotation/blast_xml/kegg.xml", "", "鉴定蛋白与数据库比对结果文件目录"],
            ["Annotation/blast_xml/string.xml", "", "鉴定蛋白与STRING数据库比对结果xml形式"],
            ["Annotation/blast_xml/pfam_domain", "", "鉴定蛋白与pfam数据库比对结果xml形式"],
            ["Annotation/blast_xml/nr.xml", "", "鉴定蛋白与NR数据库比对结果xml形式"],
            ["Annotation/blast_xml/protein.xls", "", "鉴定蛋白详情表"],
            ["Annotation/blast_xml/protein.fa", "", "蛋白序列"],
            ["DiffPep", "", "全部样本差异表达详情目录"],
            ["DiffPep/diff_protein_summary.xls", "", "所有组差异信息汇总表"],
            ["DiffPep/*.diff.xls", "", "样本间差异比较信息表"],
            ["DiffPep/*.diff.detail.xls", "", "样本间差异比较信息与蛋白详细注释信息表"],
            ["Express", "", "蛋白定量分析结果目录"],
            ["Express/exp.xls", "", "蛋白在各个样本间质谱峰面积表"],
            ["1_DataInfo", "", "搜库数据信息"],
            ["1_DataInfo/01_SearchDB", "", "搜库结果目录"],
            ["1_DataInfo/01_SearchDB/protein_sliced.xls", "", "搜库得到的蛋白详细信息"],
            ["1_DataInfo/01_SearchDB/search_db.xls", "", "搜库得到的蛋白详细信息以及在各个样本中的丰度"],
            ["1_DataInfo/02_QcInfo", "", "质控结果目录"],
            ["1_DataInfo/02_QcInfo/pep_len.pdf", "", "肽段长度分布图"],
            ["1_DataInfo/02_QcInfo/pep_num.pdf", "", "肽段数量分布图"],
            ["1_DataInfo/02_QcInfo/pep_error.pdf", "", "匹配误差分布图"],
            ["1_DataInfo/02_QcInfo/pep_error.png", "", "匹配误差分布图"],
            ["1_DataInfo/02_QcInfo/info.pdf", "", "蛋白质信息统计图"],
            ["1_DataInfo/02_QcInfo/mw.pdf", "", "蛋白质分子量分布图"],
            ["1_DataInfo/02_QcInfo/seq_cov.pdf", "", "蛋白质覆盖度分布图"],
            ["2_SampleComp/01_SamPca", "", "样本PCA分析结果目录"],
            ["2_SampleComp", "", "样本间比较分析"],
            ["2_SampleComp/01_SamPca/pca_rotation.xls", "", "蛋白相关主成分贡献度"],
            ["2_SampleComp/01_SamPca/pca_sites.xls", "", "PCA分析样本坐标表"],
            ["2_SampleComp/01_SamPca/pca_importance.xls", "", "PCA主成分解释表"],
            ["2_SampleComp/01_SamPca/sam_pca.pdf", "", "PCA分析图"],
            ["2_SampleComp/02_SamCorr", "", "样本相关性分析结果目录"],
            ["2_SampleComp/02_SamCorr/sample_correlation.xls", "", "样本间相关性系数表"],
            ["2_SampleComp/02_SamCorr/expression_matrix.xls", "", "蛋白表达丰度矩阵表"],
            ["2_SampleComp/02_SamCorr/all.exp.pdf", "", "相关性热图"],
            ["2_SampleComp/03_ExpVenn", "", "样本间venn分析结果目录"],
            ["2_SampleComp/03_ExpVenn/venn_graph.xls", "", "各样本有效蛋白列表"],
            ["2_SampleComp/03_ExpVenn/exp_venn.pdf", "", "样本间venn图"],
            ["3_Annotation", "", "蛋白注释结果目录"],
            ["3_Annotation/01_Stat", "", "蛋白注释统计结果目录"],
            ["3_Annotation/01_Stat/proteins_anno_detail.xls", "", "鉴定蛋白附加6数据库详细信息的表"],
            ["3_Annotation/01_Stat/all_annotation_statistics.xls", "", "6个数据库中相关蛋白数以及总蛋白占比"],
            ["3_Annotation/01_Stat/anno_stat.pdf", "", "注释结果柱状图"],
            ["3_Annotation/01_Stat/venn", "", "鉴定结果中蛋白在各个数据库注释结果统计VENN图目录"],
            ["3_Annotation/01_Stat/venn/pfam_venn.txt", "", "有与pfam数据库相关信息的蛋白"],
            ["3_Annotation/01_Stat/venn/kegg_venn.txt", "", "有与kegg数据库相关信息的蛋白"],
            ["3_Annotation/01_Stat/venn/go_venn.txt", "", "有与go数据库相关信息的蛋白"],
            ["3_Annotation/01_Stat/venn/cog_venn.txt", "", "有与cog数据库相关信息的蛋白"],
            ["3_Annotation/01_Stat/venn/nr_venn.txt", "", "有与nr数据库相关信息的蛋白"],
            ["3_Annotation/01_Stat/venn/subloc_venn.txt", "", "有与subloc数据库相关信息的蛋白"],
            ["3_Annotation/02_GO", "", "鉴定蛋白与GO数据库比对结果详情目录"],
            ["3_Annotation/02_GO/blast2go.annot", "", "鉴定蛋白与GO的对应关系以及GOterm名称"],
            ["3_Annotation/02_GO/query_gos.list", "", "鉴定蛋白与GO的对应关系"],
            ["3_Annotation/02_GO/go12level_statistics.xls", "", "鉴定蛋白的GO1、2层级的详细信息"],
            ["3_Annotation/02_GO/go123level_statistics.xls", "", "鉴定蛋白的GO1、2、3层级的详细信息"],
            ["3_Annotation/02_GO/go1234level_statistics.xls", "", "鉴定蛋白的GO1、2、3、4层级的详细信息"],
            ["3_Annotation/04_COG", "", "鉴定蛋白与COG数据库比对结果详情目录"],
            ["3_Annotation/04_COG/cog_table.xls", "", "鉴定蛋白序列文件与COG比对结果表"],
            ["3_Annotation/04_COG/cog_summary.xls", "", "COG数据库注释结果汇总"],
            ["3_Annotation/04_COG/cog_list.xls", "", "鉴定蛋白与COG的对应关系"],
            ["3_Annotation/04_COG/cog_bar.pdf", "", "COG分类统计柱状图"],
            ["3_Annotation/03_KEGG", "", "鉴定蛋白与KEGG数据库比对结果详情目录"],
            ["3_Annotation/03_KEGG/pathways", "", "KEGG通路注释结果图片目录"],
            ["3_Annotation/03_KEGG/pathways/.*png", "", "通路注释结果图片-位图"],
            ["3_Annotation/03_KEGG/pathways/.*pdf", "", "通路注释结果图片-矢量图"],
            ["3_Annotation/03_KEGG/pathway_table.xls", "", "每行以通路为单位显示通路注释详情表"],
            ["3_Annotation/03_KEGG/kegg_table.xls", "", "每行以蛋白为单位显示通路注释详情表"],
            ["3_Annotation/03_KEGG/pid.txt", "", "每个注释通路和映射上去的蛋白联系表"],
            ["3_Annotation/03_KEGG/path_class.pdf", "", "Pathway分类统计柱状图"],
            ["3_Annotation/03_KEGG/key_path.pdf", "", "重要通路统计图"],
            ["3_Annotation/06_SubLoc", "", "蛋白亚细胞定位注释结果目录"],
            ["3_Annotation/06_SubLoc/multiloc_stat.xls", "", "亚细胞位置对应蛋白的数量表"],
            ["3_Annotation/06_SubLoc/multiloc.xls", "", "亚细胞定位结果表"],
            ["3_Annotation/06_SubLoc/subloc_bar.pdf", "", "Subloc注释柱状图"],
            ["3_Annotation/05_Pfam/", "", "鉴定蛋白与Pfam数据库比对结果文件目录"],
            ["3_Annotation/05_Pfam/pfam_domain", "", "鉴定蛋白与pfam数据库比对结果xml形式"],
            ["3_Annotation/05_Pfam/pfam_bar.pdf", "", "Pfam注释柱状图"],
            ["4_ExpDiff", "", "全部样本差异表达详情目录"],
            ["4_ExpDiff/diff_protein_summary.xls", "", "所有组差异信息汇总表"],
            ["4_ExpDiff/*.diff.xls", "", "样本间差异比较信息表"],
            ["4_ExpDiff/*.diff.detail.xls", "", "样本间差异比较信息与蛋白详细注释信息表"],
            ["1_DataInfo/03_ExpInfo", "", "蛋白定量分析结果目录"],
            ["1_DataInfo/03_ExpInfo/exp.xls", "", "蛋白在各个样本间质谱峰面积表"],
            ["5_Proteinset", "", "差异蛋白相关分析结果目录"],
            ["5_Proteinset/01_PsetCluster", "", "差异蛋白聚类分析结果目录"],
            ["5_Proteinset/02_PsetVenn", '', '蛋白集Venn分析结果目录'],
            ["5_Proteinset/02_PsetVenn/venn.pdf", '', 'Venn图'],
            ["5_Proteinset/03_PsetAnno", "", "差异蛋白分类注释相关分析结果目录"],
            ["5_Proteinset/03_PsetAnno/01_PsetGO", "", "差异蛋白GO分类注释分析结果目录"],
            ["5_Proteinset/03_PsetAnno/02_PsetKEGG", "", "差异蛋白kegg分类注释分析结果目录"],
            ["5_Proteinset/03_PsetAnno/03_PsetCOG", "", "差异蛋白cog分类注释分析结果目录"],
            ["5_Proteinset/03_PsetAnno/04_PsetPfam", "", "差异蛋白Pfam分类注释分析结果目录"],
            ["5_Proteinset/03_PsetAnno/05_PsetSubLoc", "", "差异蛋白亚细胞定位分类注释分析结果目录"],
            ["5_Proteinset/04_PsetEnrich", "", "差异蛋白富集相关分析结果目录"],
            ["5_Proteinset/04_PsetEnrich/01_EnrichGO", "", "差异蛋白GO富集分析结果目录"],
            ["5_Proteinset/04_PsetEnrich/02_EnrichKEGG", "", "差异蛋白KEGG富集分析结果目录"],
            ["5_Proteinset/04_PsetEnrich/03_EnrichChord", "", "差异蛋白弦图结果目录"],
            ["5_Proteinset/05_PsetPPI", "", "差异蛋白蛋白互作网络结果目录"],
            ["5_Proteinset/05_PsetStringPic", "", "差异蛋白String数据库爬虫结果目录"],
            ["5_Proteinset/06_PsetIpath", "", "差异蛋白Ipath分析结果目录"],
        ]

        regpaths = [
            [r"1_DataInfo/02_QcInfo/.*_assessment.pdf", '', '预处理效果图'],
            [r"1_DataInfo/02_QcInfo/.*inter-cv.pdf", '', '组间CV分布图'],
            [r"1_DataInfo/02_QcInfo/.*_intra-cv.pdf", '', '组内CV分布图'],
            [r"3_Annotation/02_GO/go_lev.*_pie.pdf", "", "饼图"],
            [r"3_Annotation/02_GO/go_lev.*_bar.pdf", "", "直方图"],
            [r"4_ExpDiff/.*_vs_.*/.*_vs_.*scatter.pdf", '', '差异散点图'],
            [r"4_ExpDiff/.*_vs_.*/.*_vs_.*volcano.pdf", '', '差异火山图'],
            [r"5_Proteinset/01_PsetCluster/cluster_.*/cluster_group.*pdf", '', '聚类热图'],
            [r"5_Proteinset/01_PsetCluster/cluster_.*/sub.*pdf", '', '折线图'],
            [r"5_Proteinset/03_PsetAnno/01_PsetGO/.*_vs_.*/.*_vs_.*pdf", '', 'GO注释柱状图'],
            [r"5_Proteinset/03_PsetAnno/02_PsetKEGG/.*_vs_.*/.*_vs_.*pdf", '', 'Pathway分类统计柱状图'],
            [r"5_Proteinset/03_PsetAnno/03_PsetCOG/.*_vs_.*/.*_vs_.*pdf", '', 'COG注释统计图'],
            [r"5_Proteinset/03_PsetAnno/04_PsetPfam/.*_vs_.*/.*_vs_.*pdf", '', 'Pfam注释统计柱状图'],
            [r"5_Proteinset/03_PsetAnno/05_PsetSubLoc/.*_vs_.*/.*_vs_.*pdf", '', 'SubLoc注释统计柱状图'],
            [r"5_Proteinset/04_PsetEnrich/01_EnrichGO/.*_vs_.*/.*_vs_.*bar\.pdf", '', '富集直方图'],
            [r"5_Proteinset/04_PsetEnrich/01_EnrichGO/.*_vs_.*/.*_vs_.*bubble.*pdf", '', '富集气泡图'],
            [r"5_Proteinset/04_PsetEnrich/02_EnrichKEGG/.*_vs_.*/.*_vs_.*bar\.pdf", '', '富集直方图'],
            [r"5_Proteinset/04_PsetEnrich/02_EnrichKEGG/.*_vs_.*/.*_vs_.*bubble.*pdf", '', '富集气泡图'],
            [r"5_Proteinset/04_PsetEnrich/03_EnrichChord/.*_vs_.*/.*_vs_.*pdf", '', '富集弦图'],
            [r"5_Proteinset/05_PsetPPI/.*_vs_.*/ppinetwork_topology/ppi\.centrality\.line\.pdf", '', '网络中心系数分布图'],
            [r"5_Proteinset/05_PsetPPI/.*_vs_.*/ppinetwork_topology/ppi\.degree\.line\.pdf", '', '网络节点度分布图'],
        ]

        sdir = self.add_upload_dir(target_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regpaths)
        if self.option("report_img"):
            self.logger.info("开始上传报告中的png图片")
            s3 = self._sheet.output.split(":")[0]
            report_img_dir = self.chart.work_dir + '/png/'
            report_img_s3 = s3 + "://commonbucket/files/report_img/diav3/" + self.task_id + "/"
            self.upload_to_s3(report_img_dir, report_img_s3)
            self.logger.info("结束上传报告中的png图片")

    def run_api(self):
        import traceback
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.logger.info("导表开始")
        self.stop_timeout_check()
        task_info = self.api.api('task_info.dia_task_info')
        task_info.add_task_info()
        self.export_specimen()
        self.export_searchdb()
        self.export_qc()
        if not self.option("protein_group").prop["sample_number"] < 3:
            self.export_sam_corr()
            self.export_sam_pca()
        self.export_exp_venn()
        self.export_annotation()
        self.export_pep_diff()
        self.export_proteinset()
        self.export_proteinset_venn()
        self.logger.info("开始导差异蛋白分析相关表，如果错误会尽量忽略")
        try:
            self.export_diff_cluster()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_cog_class()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_go_class()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_kegg_class()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_pfam_stat()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_subloc_stat()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_go_enrich()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_kegg_enrich()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_circle()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_ppi()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_ipath()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        try:
            self.export_diff_string_picture()
        except:
            estr = traceback.format_exc()
            self.logger.info(estr)
        # self.export_diff_cluster()
        # self.export_diff_cog_class()
        # self.export_diff_go_class()
        # self.export_diff_kegg_class()
        # self.export_diff_pfam_stat()
        # self.export_diff_subloc_stat()
        # self.export_diff_go_enrich()
        # self.export_diff_kegg_enrich()
        # self.export_diff_ppi()
        # self.export_diff_string_picture()
        # self.export_diff_circle()
        # self.export_diff_ipath()
        # super(DiaWorkflow, self).end()
        self.export_project_overview()
        if self.option("report_img"):
            self.export_report_img()
        self.logger.info("导表完成")

    def export_report_img(self):
        report_config = os.path.join(self.chart.work_dir, 'report_config.json')
        api =  self.api.api('dia.report_model')
        s3 = self._sheet.output.split(":")[0]
        report_img_s3 = s3 + ":commonbucket/files/report_img/diav3/" + self.task_id
        api.add_report_image(self.task_id, report_config, report_img_s3)

    @time_count
    def export_project_overview(self):
        self.api_project_overview = self.api.api("dia.project_overview")
        params = dict(
            task_id=self.task_id,
            specimen_count=str(len([i for item in self.group_detail for i in item]))+'/'+str(len(self.group_detail)),
            database_type=self.option('data_source'),
            kegg_class=self.option("kegg_class"),
            kegg_name=self.option('kegg_org'),
            p_value_fdr=self.option('pvalue'),
            fc_1=self.option('fc_up'),
            fc_2=self.option('fc_down'),
            spectrum=self.proteininfo_data_dic[0]['identified_spectrum'],
            spectrum_total=self.proteininfo_data_dic[0]['total_spectrum'],
            spectrum_percent='{:.2%}'.format(float(self.proteininfo_data_dic[0]['identified_spectrum'])/float(\
                self.proteininfo_data_dic[0]['total_spectrum'])),
            protein=self.proteininfo_data_dic[0]['protein_num']
        )
        self.api_project_overview.add_project_overview(params=params)
        self.api_project_overview.add_project_protein(task_id=self.task_id)
        self.api_project_overview.add_project_error(task_id=self.task_id)
        self.api_project_overview.add_project_expcorr(task_id=self.task_id, \
            express_file=self.preprocess.option('preprocess_exp').prop['path'],\
            group_file=self.option('protein_group').prop['path'])


    @time_count
    def export_specimen(self):
        gevent.sleep()
        self.api_specimen = self.api.api("dia.specimen")
        self.api_specimen.add_specimen(group_txt=self.option("protein_group").prop["path"], params=None,
                project_sn=self.project_sn, task_id=self.task_id, type="dia")
        self.group_id, self.group_detail, self.group_category = self.api_specimen.add_specimen_group(self.option("protein_group").prop["path"])
        self.control_id, self.compare_detail = self.api_specimen.add_control_group(self.option("protein_control").prop["path"], self.group_id)

    @time_count
    def export_searchdb(self):
        protein_selected = self.searchdb.work_dir + "/protein_sliced.xls"
        if self.option('change_des'):
            nr_file = os.path.join(self.annotation.output_dir, 'anno_stat', 'blast', 'nr.xls')
            nr_df = pd.read_csv(nr_file, sep='\t')
            acc2des = dict()
            for acc, des in zip(nr_df['Query-Name'], nr_df['Hit-Description']):
                if not acc in acc2des:
                    acc2des[acc] = des
            sel_df = pd.read_csv(protein_selected, sep='\t', index_col=0)
            for ind in sel_df.index:
                if ind in acc2des:
                    sel_df.loc[ind, 'Description'] = acc2des[ind]
            sel_df.to_csv(protein_selected, index=True, header=True, sep='\t')
        gevent.sleep()
        # if self.new:
        #     self.api_searchdb = self.api.api("dia.searchdb_new")
        # else:
        self.api_searchdb = self.api.api("dia.searchdb")
        self.api_searchdb.add_searchdb(project_sn=self.project_sn,
                task_id=self.task_id, name=None, protein_selected = protein_selected)

    @time_count
    def export_sam_corr(self):
        gevent.sleep()
        self.api_sam_corr = self.api.api("dia.all_exp")
        params = dict(
            task_id=self.task_id,
            submit_location='expresscorr',
            group_id=str(self.group_id),
            group_dict=self.option('protein_group').prop['group_dict'],
            scm="complete",
            scd="euclidean",
            corr_method="pearson",
            task_type=2,
            express_id=str(self.express_id)
        )
        corr_output = self.sam_corr.work_dir
        self.api_sam_corr.add_exp_corr3(corr_output, params=params, main_id=None,
                      project_sn=self.project_sn, task_id=self.task_id)

    @time_count
    def export_qc(self):
        gevent.sleep()
        if self.option('DIA'):
            self.api_dia = self.api.api("dia.dia_quality_control")
            # self.api_peperror = self.api.api("dia.peperror_new")
            params = {"software": "DIA"}
            exp_params = {
            # 'raw_path': self.work_dir + "/raw_treat_ref",
            # 'group_table': self.option("protein_group").prop["path"],
            'all_eliminate':self.option('all_eliminate'),
            'all_percent':self.option('all_percent'),
            'if_group':self.option('if_group'),
            'fillna':self.option('fillna'),
            'fill_type':self.option('fill_type'),
            'submit_location':'preprocess',
            'task_id':self.task_id,
            'task_type':2
            }
            if self.option('if_group') == "yes":
                exp_params.update({
                    "group_specific": self.option('group_specific'),
                    "group_percent": self.option('group_percent')
                    })
            # self.api_dia.add_peplen(self.option('report_dia').prop['path'], params=params)
            # self.api_dia.add_pepnum(self.option('report_dia').prop['path'], params=params)
            # self.api_dia.add_proteininfo(self.option('exp_ana_dia').prop['path'], params=params)
            # self.api_dia.add_proteinmw(self.option('protein_dia').prop['path'], params=params)
            # self.api_dia.add_coverage(self.option('protein_dia').prop['path'], params=params)
            self.api_dia.add_peplen(self.option('report_dia').prop['path'], params=params)
            self.api_dia.add_pepnum(self.option('report_dia').prop['path'], params=params)
            self.proteininfo_data_dic = self.api_dia.add_proteininfo(self.option('exp_ana_dia').prop['path'], params=params)
            # self.api_dia.add_proteinmw(self.option('protein_dia'), params=params)
            self.api_dia.add_proteinmw(self.option('protein').prop['path'], params=params)
            # self.api_dia.add_coverage(self.option('protein_dia'), params=params)
            self.api_dia.add_coverage(self.option('protein').prop['path'], params=params)
            # self.api_peperror.add_peperror(params=params, main_id=None,
            #                                project_sn=self.project_sn, task_id=self.task_id,
            #                                peperror_exp=self.option("psm").prop['path'],
            #                                s3_png_path=os.path.join(self._sheet.output,"1_DataInfo/02_QcInfo/pep_error.png"))
            self.api_express = self.api.api("dia.express")
            self.api_preprocess = self.api.api("dia.preprocess")
            # self.express_id_raw = self.api_express.add_express(params=params, main_id=None,
            #                              project_sn=self.project_sn, task_id=self.task_id,
            #                              express_exp=self.work_dir + "/raw_treat_ref",name="Raw")
            self.express_id_raw = self.api_preprocess.add_preprocess_exp(
                exp_path=self.work_dir + "/raw_treat_ref",
                cv_path=os.path.join(self.preprocess.work_dir, "{}_cv_stats_raw.xls".format(self.option("fillna"))),
                cv_summary=os.path.join(self.preprocess.work_dir, "{}_cv_summary_raw.xls".format(self.option("fillna"))),
                searchdb=self.option("protein").prop['path'], params=params, main_id=None, project_sn=self.project_sn,
                task_id=self.task_id, exp_type="raw")
            time.sleep(5) #确保上面add_express存的name为raw的表 在下面add_preprocess_exp存的名为Protein_Table_Origin的表 的时间点之前.
            self.express_id = self.api_preprocess.add_preprocess_exp(\
                exp_path=self.preprocess.option('preprocess_exp').prop['path'], \
                cv_path=os.path.join(self.preprocess.work_dir, "{}_cv_stats.xls".format(self.option("fillna"))),\
                nrmse_path=os.path.join(self.preprocess.work_dir, "{}_nrmse.txt".format(self.option("fillna"))),\
                na_path=os.path.join(self.preprocess.work_dir, "{}_na_stats.xls".format(self.option("fillna"))),\
                cv_summary=os.path.join(self.preprocess.work_dir, "{}_cv_summary.xls".format(self.option("fillna"))),\
                searchdb=self.option("protein").prop['path'],
                main_id=None,project_sn=self.project_sn,task_id=self.task_id,params=exp_params,raw_exp_id=str(self.express_id_raw))
            # self.api_preprocess.add_cv_box(\
            #     cv_matrix=os.path.join(self.preprocess.work_dir, "{}_cv_stats.xls".format(self.option("fillna"))),\
            #     express_id=self.express_id)
            return
        if self.new:
            self.api_coverage = self.api.api("dia.coverage_new")
            self.api_express = self.api.api("dia.express")
            self.api_preprocess = self.api.api("dia.preprocess")
            self.api_peperror = self.api.api("dia.peperror_new")
            self.api_peplen = self.api.api("dia.peplen_new")
            self.api_pepnum = self.api.api("dia.pepnum_new")
            self.api_proteininfo = self.api.api("dia.proteininfo_new")
            self.api_proteinmw = self.api.api("dia.proteinmw_new")
        else:
            self.api_coverage = self.api.api("dia.coverage")
            self.api_express = self.api.api("dia.express")
            self.api_preprocess = self.api.api("dia.preprocess")
            self.api_peperror = self.api.api("dia.peperror")
            self.api_peplen = self.api.api("dia.peplen")
            self.api_pepnum = self.api.api("dia.pepnum")
            self.api_proteininfo = self.api.api("dia.proteininfo")
            self.api_proteinmw = self.api.api("dia.proteinmw")
        params = {"software": "peaks 8.5"}
        exp_params = {
        # 'raw_path': self.work_dir + "/raw_treat_ref",
        # 'group_table': self.option("protein_group").prop["path"],
        'all_eliminate':self.option('all_eliminate'),
        'all_percent':self.option('all_percent'),
        'if_group':self.option('if_group'),
        'fillna':self.option('fillna'),
        'fill_type':self.option('fill_type'),
        'submit_location':'preprocess',
        'task_id':self.task_id,
        'task_type':2
        }
        if self.option('if_group') == "yes":
            exp_params.update({
                "group_specific": self.option('group_specific'),
                "group_percent": self.option('group_percent')
                })
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        corr_output = self.sam_corr.work_dir
        self.api_coverage.add_coverage(params=params, main_id=None,
                                       project_sn=self.project_sn, task_id=self.task_id,
                                       coverage_exp=self.option("protein").prop['path'])
        # self.express_id_raw = self.api_express.add_express(params=params, main_id=None,
        #                              project_sn=self.project_sn, task_id=self.task_id,
        #                              express_exp=self.work_dir + "/raw_treat_ref",name="Raw")
        self.express_id_raw = self.api_preprocess.add_preprocess_exp(
            exp_path=self.work_dir + "/raw_treat_ref",
            cv_path=os.path.join(self.preprocess.work_dir, "{}_cv_stats_raw.xls".format(self.option("fillna"))),
            cv_summary=os.path.join(self.preprocess.work_dir, "{}_cv_summary_raw.xls".format(self.option("fillna"))),
            searchdb=self.option("protein").prop['path'], params=params, main_id=None, project_sn=self.project_sn,
            task_id=self.task_id, exp_type="raw")
        self.express_id = self.api_preprocess.add_preprocess_exp(\
            exp_path=self.preprocess.option('preprocess_exp').prop['path'], \
            cv_path=os.path.join(self.preprocess.work_dir, "{}_cv_stats.xls".format(self.option("fillna"))),\
            nrmse_path=os.path.join(self.preprocess.work_dir, "{}_nrmse.txt".format(self.option("fillna"))),\
            na_path=os.path.join(self.preprocess.work_dir, "{}_na_stats.xls".format(self.option("fillna"))),\
            cv_summary=os.path.join(self.preprocess.work_dir, "{}_cv_summary.xls".format(self.option("fillna"))),\
            searchdb=self.option("protein").prop['path'],
            main_id=None,project_sn=self.project_sn,task_id=self.task_id,params=exp_params,raw_exp_id=str(self.express_id_raw))
        # self.api_preprocess.add_cv_box(\
        #     cv_matrix=os.path.join(self.preprocess.work_dir, "{}_cv_stats.xls".format(self.option("fillna"))),\
        #     express_id=self.express_id)
        self.api_peplen.add_peplen(params=params, main_id=None,
                                   project_sn=self.project_sn, task_id=self.task_id,
                                   peplen_exp=self.option("peptide").prop['path'])
        self.api_pepnum.add_pepnum(params=params, main_id=None,
                                   project_sn=self.project_sn, task_id=self.task_id,
                                   pepnum_exp=self.option("protein").prop['path'])
        self.proteininfo_data_dic = self.api_proteininfo.add_proteininfo(params=params, main_id=None,
                                             project_sn=self.project_sn, task_id=self.task_id,
                                             proteininfo_exp=self.option("protein_information").prop['path'])
        self.api_proteinmw.add_proteinmw(params=params, main_id=None,
                                             project_sn=self.project_sn, task_id=self.task_id,
                                             proteinmw_exp=self.option("protein").prop['path'])
        # self.api_peperror.add_peperror(params=params, main_id=None,
        #                                project_sn=self.project_sn, task_id=self.task_id,
        #                                peperror_exp=self.option("psm").prop['path'],
        #                                s3_png_path=os.path.join(self._sheet.output,"1_DataInfo/02_QcInfo/pep_error.png"))

    @time_count
    def export_sam_pca(self):
        gevent.sleep()
        # self.api_sam_pca = self.api.api("labelfree.pca")
        all_exp = self.api.api("dia.all_exp")
        params = dict(
            task_id=self.task_id,
            submit_location="expresspca",
            group_id=str(self.group_id),
            group_dict=self.option('protein_group').prop['group_dict'],
            task_type=2,
            express_id=str(self.express_id)
        )
        #os.remove(self.sam_pca.output_dir + "/pca_rotation.xls")
        #os.rename(self.sam_pca.output_dir + "/pca_rotation_all.xls", self.sam_pca.output_dir + "/pca_rotation.xls")
        # pca_output = self.sam_pca.output_dir
        # self.api_sam_pca.add_pca(pca_output, exp_id=None, params=params,
        #             project_sn=self.project_sn, task_id=self.task_id, main_id=None)
        # pca_output = self.sam_pca.work_dir
        pca_output = self.sam_pca.output_dir
        pca_main_id = all_exp.add_exp_pca2(pca_output, params=params,
                                  project_sn=self.project_sn, task_id=self.task_id, main_id=None)
        if self.min_group_num >= 3:
            all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(pca_main_id))

    @time_count
    def export_exp_venn(self):
        gevent.sleep()
        venn_group = self.option('protein_group').prop['path']+'_for_exp_venn'
        if os.path.exists(venn_group):
            group_id, group_detail, group_category = self.api_specimen.add_specimen_group(venn_group)
            group_dict = OrderedDict()
            with open(venn_group, "r") as f:
                f.readline()
                for line in f:
                    tmp = line.strip().split("\t")
                    group_dict.setdefault(tmp[1], list())
                    if tmp[0] not in group_dict[tmp[1]]:
                        group_dict[tmp[1]].append(tmp[0])
        else:
            group_dict = self.option('protein_group').prop['group_dict']
            group_id = self.group_id
        api_exp_venn = self.api.api("dia.all_exp")
        params = dict(
            task_id=self.task_id,
            submit_location="expvenn",
            group_id=str(group_id),
            group_dict=group_dict,
            task_type=2,
            express_id=str(self.express_id)
        )
        graph_table = os.path.join(self.exp_venn.output_dir, 'venn_graph.xls')
        api_exp_venn.add_exp_venn(graph_table, params=params, main_id=None)

    @time_count
    def export_pep_diff(self):
        gevent.sleep()
        self.api_diff = self.api.api("dia.diff")

        params = dict(
            task_id=self.task_id,
            submit_location="diff",
            group_id=str(self.group_id),
            group_dict=self.option('protein_group').prop['group_dict'],
            type= "origin", 
            correct_method=self.option('correct_method'),
            padjust_method='fdr',
            fc_up=str(self.option('fc_up')), 
            fc_down=str(self.option('fc_down')),
            pvalue=str(self.option('pvalue')), 
            diff_method=self.option('method_type'),
            control_id=str(self.control_id), 
            task_type=2,
            express_id=str(self.express_id),
            log=self.option('log'),
            sig_type='pvalue'
            # fill_type=self.option('fill_type')
        )

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        work_dir = self.diff_pep.work_dir
        compare_dict_xls = work_dir + "/compare_dict.xls"
        num_summary_xls = work_dir + "/num_summary.xls"
        allsummary_xls = work_dir + "/all_summary.xls"
        self.api_diff.add_diff(work_dir, 
            compare_dict_xls=compare_dict_xls,
            num_summary_xls=num_summary_xls,
            allsummary_xls=allsummary_xls,
            params=params,
            project_sn=self.project_sn,
            task_id=self.task_id,
            method_type=self.option("method_type"),
            fill_type=self.option('fill_type')
            )

    @time_count
    def export_annotation(self):
        '''
        导入注释结果表格 liubinxu
        '''
        gevent.sleep()
        if self.new:
            self.annotation_api = self.api.api("dia.protein_annotation")
        else:
            self.annotation_api = self.api.api("dia.protein_annotation")
        annotation_result_dir = self.annotation.output_dir
        params = {
            "go_evalue": str(self.option("go_evalue")),
            "cog_evalue": str(self.option("cog_evalue")),
            "kegg_evalue": str(self.option("kegg_evalue")),
            "pfam_evalue": str(self.option("pfam_evalue")),
            "go_identity": str(self.option("go_identity")),
            "cog_identity": str(self.option("cog_identity")),
            "kegg_identity": str(self.option("kegg_identity")),
            "taxon": str(self.option("kegg_class")),
            "kegg_species":str(self.option("kegg_org")),
            "submit_location": "annotationstat",
            "task_id": self.task_id,
            "task_type": 2
        }
        self.type = self.annotation_api.run(annotation_result_dir, params, change_des=self.option('change_des'))

    @time_count
    def export_proteinset(self):
        gevent.sleep()
        self.logger.info("导入蛋白集")
        self.export_proteinset = self.api.api("dia.proteinset")
        self.export_proteinset.add_proteinset(diff_work_dir=self.diff_pep.work_dir, group_id=str(self.group_id),
                                              task_id=self.task_id, project_sn=self.project_sn)
    @time_count
    def export_proteinset_venn(self):
        gevent.sleep()
        self.logger.info("导入蛋白集以绘制页面初始venn图(只在蛋白集数目>=2情况下，超过6个则取前6个蛋白集绘制venn)")
        self.export_proteinset = self.api.api("dia.proteinset")
        self.export_proteinset.add_proteinset_venn(task_id=self.task_id, project_sn=self.project_sn)

    @time_count
    def export_diff_cluster(self):
        self.logger.info("导入差异蛋白聚类结果")
        all_exp = self.api.api("dia.all_exp")
        cluster_out = self.diff_cluster.output_dir
        for dir in os.listdir(cluster_out):
            if os.path.isdir(os.path.join(cluster_out, dir)):
                all_exp.add_geneset_cluster(os.path.join(cluster_out, dir), main_id='', express_id=str(self.express_id))

    @time_count
    def export_diff_cog_class(self):
        self.logger.info("导入差异蛋白cog分类注释")
        api_proteinset = self.api.api('dia.proteinset')
        for cog_result in glob.glob(os.path.join(self.diff_cog_class.output_dir, '*', '*cog_class_table.xls')):
            api_proteinset.add_proteinset_cog_detail(cog_result, '')

    @time_count
    def export_diff_go_class(self):
        self.logger.info("导入差异蛋白go分类注释")
        api_proteinset = self.api.api('dia.proteinset')
        for go_result in glob.glob(os.path.join(self.diff_go_class.output_dir, '*', '*go_class_table.xls')):
            api_proteinset.add_go_regulate_detail(go_result, '')
            if '_up_down' not in go_result:
                continue
            api_proteinset.add_go_regulate_detail2(go_result, '')

    @time_count
    def export_diff_pfam_stat(self):
        self.logger.info("导入差异蛋白pfam分类注释")
        api_proteinset = self.api.api('dia.proteinset')
        for pfam_stat in glob.glob(os.path.join(self.diff_pfam_stat.output_dir, '*', '*stat.xls')):
            api_proteinset.add_pfam_stat('', pfam_stat)

    @time_count
    def export_diff_subloc_stat(self):
        self.logger.info("导入差异蛋白subloc分类注释")
        api_proteinset = self.api.api('dia.proteinset')
        for subloc_stat in glob.glob(os.path.join(self.diff_subloc_stat.output_dir, '*', '*stat.xls')):
            api_proteinset.add_subloc_stat('', subloc_stat)

    @time_count
    def export_diff_kegg_class(self):
        self.logger.info("导入差异蛋白kegg分类注释")
        api_proteinset = self.api.api('dia.proteinset')
        kegg_out = self.diff_kegg_class.output_dir
        kegg_work = self.diff_kegg_class.work_dir
        for dir in os.listdir(kegg_out):
            if os.path.isdir(os.path.join(kegg_out, dir)):
                kegg_stat = os.path.join(kegg_out, dir, 'kegg_stat.xls')
                pathway_file = os.path.join(kegg_out, dir, 'pathways')
                protein_kegg = os.path.join(kegg_work, dir+'_protein.list')
                kegg_table_2 = os.path.join(kegg_work, 'protein_kegg_level_table.xls')
                record_id = api_proteinset.add_kegg_regulate_new2("", protein_kegg, kegg_stat, kegg_table_2)
                api_proteinset.add_kegg_regulate_pic(record_id, kegg_table_2, pathway_file)

    @time_count
    def export_diff_go_enrich(self):
        self.logger.info("导入差异蛋白go富集结果")
        api_proteinset = self.api.api('dia.proteinset')
        go_out = self.diff_go_enrich.output_dir
        for dir in os.listdir(go_out):
            if os.path.isdir(os.path.join(go_out, dir)):
                go_adjust_png = os.path.join(go_out, dir, 'adjust_lineage.png')
                go_adjust_pdf = os.path.join(go_out, dir, 'adjust_lineage.pdf')
                enrich_file = glob.glob(os.path.join(go_out, dir, '*.xls'))[0]
                record_id = api_proteinset.add_go_enrich_detail('', enrich_file)
                api_proteinset.update_directed_graph(record_id, go_adjust_png, go_adjust_pdf)

    @time_count
    def export_diff_kegg_enrich(self):
        self.logger.info("导入差异蛋白kegg富集结果")
        api_proteinset = self.api.api('dia.proteinset')
        kegg_out = self.diff_kegg_enrich.output_dir
        for dir in os.listdir(kegg_out):
            if os.path.isdir(os.path.join(kegg_out, dir)):
                enrich_file = glob.glob(os.path.join(kegg_out, dir, '*.xls'))[0]
                record_id = api_proteinset.add_kegg_enrich_detail('', enrich_file)
                pathway_file = os.path.join(kegg_out, dir, 'pathways')
                api_proteinset.add_kegg_enrich_pic(record_id, enrich_file, pathway_file)


    @time_count
    def export_diff_circle(self):
        self.logger.info("导入差异蛋白kegg富集炫图和GO富集炫图结果")
        api_proteinset = self.api.api('dia.proteinset')
        circle_out = self.diff_circle.output_dir
        self.logger.info(circle_out)
        for dir in os.listdir(circle_out):
            if os.path.isdir(os.path.join(circle_out, dir)):
                enrich_zscore_file = os.path.join(circle_out, dir, 'enrich_zscore')
                if dir.startswith('go_circle'):
                    enrich_circ_file = os.path.join(circle_out, dir, 'go_enrich_choose.table')
                    enrich_detail_file = os.path.join(circle_out, dir, 'go_enrich_detail.table')
                    record_id = api_proteinset.add_circ_graph('', enrich_circ_file, 'go')
                    api_proteinset.add_circ_detail(record_id, enrich_detail_file, 'go')
                    api_proteinset.update_circ_main(record_id, enrich_zscore_file, 'go')
                else:
                    enrich_circ_file = os.path.join(circle_out, dir, 'kegg_enrich_choose.table')
                    enrich_detail_file = os.path.join(circle_out, dir, 'kegg_enrich_detail.table')
                    record_id = api_proteinset.add_circ_graph('', enrich_circ_file, 'kegg')
                    api_proteinset.add_circ_detail(record_id, enrich_detail_file, 'kegg')
                    api_proteinset.update_circ_main(record_id, enrich_zscore_file, 'kegg')

    @time_count
    def export_diff_ipath(self):
        self.logger.info("导入差异蛋白KeggIpath结果")
        api_proteinset = self.api.api('dia.proteinset')
        ipath_out = self.diff_ipath.output_dir
        ipath_work = self.diff_ipath.work_dir
        for dir in os.listdir(ipath_out):
            if os.path.isdir(os.path.join(ipath_out, dir)):
                ipath_file = os.path.join(ipath_out, dir, 'gene_ipath_input.xls')
                protein_file = os.path.join(ipath_work, dir+'_protein.list')
                api_proteinset.add_ipath_detail('', ipath_file, protein_file)

    @time_count
    def export_diff_ppi(self):
        self.logger.info("导入差异蛋白互作结果")
        api_ppinetwork = self.api.api('dia.proteinset_ppi')
        ppi_out = self.diff_ppi.output_dir
        for dir in os.listdir(ppi_out):
            if os.path.isdir(os.path.join(ppi_out, dir)):
                all_nodes_path = os.path.join(ppi_out, dir, 'ppinetwork_predict', 'all_nodes.txt')  # 画图节点属性文件
                interaction_path = os.path.join(ppi_out, dir, 'ppinetwork_predict', 'interaction_detail.txt')  # 画图的边文件
                network_stats_path = os.path.join(ppi_out, dir, 'ppinetwork_predict', 'network_stats.txt')  # 网络全局属性统计
                network_centrality_path = os.path.join(ppi_out, dir, 'ppinetwork_topology', 'protein_interaction_network_centrality.txt')
                degree_distribution_path = os.path.join(ppi_out, dir, 'ppinetwork_topology', 'protein_interaction_network_degree_distribution.txt')
                table_id = api_ppinetwork.add_node_table(file_path=all_nodes_path, group_id=self.group_id,
                                              table_id='')  # 节点的属性文件（画网络图用）
                api_ppinetwork.add_edge_table(file_path=interaction_path, table_id=table_id)  # 边信息
                api_ppinetwork.add_network_attributes(file2_path=network_stats_path,
                                                      table_id=table_id)  # 网络全局属性
                api_ppinetwork.add_network_centrality(file_path=network_centrality_path, file1_path=all_nodes_path,
                                                      table_id=table_id)  # 中心信息
                api_ppinetwork.add_degree_distribution(file_path=degree_distribution_path,
                                                       table_id=table_id)

    @time_count
    def export_diff_string_picture(self):
        self.logger.info("导入差异蛋白STRING数据库爬虫结果")
        api_proteinset = self.api.api('dia.proteinset')
        string_out = self.diff_string_picture.output_dir
        for dir in os.listdir(string_out):
            if os.path.isdir(os.path.join(string_out, dir)):
                api_proteinset.add_string_picture('', os.path.join(string_out, dir))

    def build_seq_database(self):
        self.logger.info("创建蛋白数据库")
        self.export_seq = self.api.api("dia.protein_seq")
        pep = self.option("protein_fasta").prop["path"]
        seq_db = os.path.join(self.work_dir, 'seq_db.sqlite3')
        self.export_seq.build_seq_database(seq_db, pep, task_id=self.task_id)

    def add_anno_to_diff(self):
        self.logger.info("给差异蛋白文件添加注释")
        protein_anno = self.annotation.output_dir + '/anno_stat/proteins_anno_detail.xls'
        protein_anno_pd = pd.read_table(protein_anno, header=0, index_col=0)
        diff_output = self.diff_pep.output_dir
        target_files = glob.glob(diff_output + '/' + '*_vs_*_diff.xls')
        treat_ref_df = pd.read_table(self.work_dir + "/raw_treat_ref", sep="\t", header=0, index_col=0)
        for each in target_files:
            protein_pd = pd.read_table(each, header=0, index_col=0)
            protein_anno_result = pd.concat([protein_pd, protein_anno_pd, treat_ref_df], axis=1)
            # protein_anno_out = each.split('.xls')[0] + '.detail.xls'
            protein_anno_out = self.output_dir + '/DiffPep/' + os.path.basename(each).split('.xls')[0] + '.detail.xls'
            protein_anno_result.to_csv(protein_anno_out, header=True, index=True, sep='\t')

    ############################################################################
    # 下面的函数为了生成qc中的文件，提前生成这样的文件为了画静态pdf,代码参考的是导表中的函数
    def tmp_add_proteinmw(self, proteinmw_exp=None):
        df = pd.read_table(proteinmw_exp, sep ="\t", header = 0)
        mylist = df.loc[:, 'MW [kDa]'].tolist()
        num_1_20 =  num_21_40 = num_41_60 = num_61_80 = num_81_100 = \
            num_101_120 = num_121_140\
        = num_141_160 = num_161_180 = num_181_200 = num_201_220 = \
            num_221_240 = num_241_260 = num_261_280 = num_281_300 = num_300_inf = 0

        for i in mylist:
            if i <= 21:
                num_1_20 += 1
            elif i <= 41 and i > 21:
                num_21_40 += 1
            elif i <= 61 and i > 41:
                num_41_60 += 1
            elif i <= 81 and i > 61:
                num_61_80 += 1
            elif i <= 101 and i > 81:
                num_81_100 += 1
            elif i <= 121 and i > 101:
                num_101_120 += 1
            elif i <= 141 and i > 121:
                num_121_140 += 1
            elif i <= 161 and i > 141:
                num_141_160 += 1
            elif i <= 181 and i > 161:
                num_161_180 += 1
            elif i <= 201 and i > 181:
                num_181_200 += 1
            elif i <= 221 and i > 201:
                num_201_220 += 1
            elif i <= 241 and i > 221:
                num_221_240 += 1
            elif i <= 261 and i > 241:
                num_241_260 += 1
            elif i <= 281 and i > 261:
                num_261_280 += 1
            elif i <= 301 and i > 281:
                num_281_300 += 1
            elif  i > 301:
                num_300_inf += 1


        mw_dict =OrderedDict()
        mw_dict['1-21'] = num_1_20
        mw_dict['21-41'] = num_21_40
        mw_dict['41-61'] = num_41_60
        mw_dict['61-81'] = num_61_80
        mw_dict['81-101'] = num_81_100
        mw_dict['101-121'] = num_101_120
        mw_dict['121-141'] = num_121_140
        mw_dict['141-161'] = num_141_160
        mw_dict['161-181'] = num_161_180
        mw_dict['181-201'] = num_181_200
        mw_dict['201-221'] = num_201_220
        mw_dict['221-241'] = num_221_240
        mw_dict['241-261'] = num_241_260
        mw_dict['261-281'] = num_261_280
        mw_dict['281-301'] = num_281_300
        mw_dict['>301'] = num_300_inf
        qc_dir = os.path.join(self.work_dir, 'tmp_qc')
        if os.path.exists(qc_dir):
            mw_file = os.path.join(qc_dir, 'Protein_mw_distribution.xls')
            with open(mw_file, 'w') as ew:
                ew.write('moleweight' + '\t' + 'number' + '\n')
                for mw, n in mw_dict.items():
                    ew.write(mw + '\t' + str(n) + '\n')

    def tmp_add_peplen(self, peplen_exp=None):
        try:
            report_df = pd.read_csv(peplen_exp, sep='\t')
        except:
            report_df = pd.read_excel(peplen_exp, sep='\t')
        peps = list()
        len2num = dict()
        if 'PEP.GroupingKey' in report_df.columns:
            PEP_GroupingKey = 'PEP.GroupingKey'
        else:
            PEP_GroupingKey = 'PEP.StrippedSequence'
        for pep in report_df[PEP_GroupingKey]:
            pep = ''.join(list(filter(str.isalpha, pep)))
            len_ = len(pep)
            if pep in peps or not len_:
                continue
            peps.append(pep)
            if len_ not in len2num:
                len2num[len_] = 0
            len2num[len_] += 1
        del report_df
        len2num_s = OrderedDict()
        for l in sorted(len2num.keys()):
            len2num_s[str(l)] = len2num[l]
        qc_dir = os.path.join(self.work_dir, 'tmp_qc')
        if os.path.exists(qc_dir):
            len_file = os.path.join(qc_dir, 'Peptite_length_distribution.xls')
            with open(len_file, 'w') as ew:
                ew.write('length' + '\t' + 'number' + '\n')
                for l, n in len2num_s.items():
                    ew.write(l + '\t' + str(n) + '\n')

    def tmp_add_proteininfo(self, proteininfo_exp=None):
        info2num = OrderedDict()
        with open(proteininfo_exp, 'r') as er:
            for line in er:
                if not line.strip():
                    continue
                if u'Sparse Profiles - ' in line:
                    info, num = line.strip().split('\t')
                    info = info.split('Sparse Profiles - ')[1]
                    num = num.split(' of ')[0]
                    info2num[info] = int(num.replace(',', ''))
        new_columns = ['total_spectrum', 'identified_spectrum', 'peptide_sequences', 'protein_num', 'protein_groups']
        info_df = pd.DataFrame([info2num])
        tmp = info_df.columns.tolist()
        tmp[-1], tmp[-2] = tmp[-2], tmp[-1]
        info_df = info_df[tmp]
        info_df.rename(columns=dict(zip(info_df.columns.tolist(), new_columns)), inplace=True)
        info_df.drop(columns=['identified_spectrum'], inplace=True)
        qc_dir = os.path.join(self.work_dir, 'tmp_qc')
        if os.path.exists(qc_dir):
            try:
                info_df.to_csv(os.path.join(qc_dir, 'Protein_infomation.xls'), index=False, header=True)
            except:
                pass

    def tmp_add_coverage(self, coverage_exp=None):
        df = pd.read_table(coverage_exp, sep ="\t", header = 0)
        try:
            mylist = df.loc[:, 'Coverage'].tolist()
        except:
            mylist = df.loc[:, 'Coverage [%]'].tolist()
        num_1 =  num_1_5 = num_5_10 = num_10_20 = num_20_40 = \
            num_40_60 = num_60_80  = num_80_inf = 0

        for i in mylist:
            if i <= 1:
                num_1 += 1
            elif i <= 5 and i > 1:
                num_1_5 += 1
            elif i <= 10 and i > 5:
                num_5_10 += 1
            elif i <= 20 and i > 10:
                num_10_20 += 1
            elif i <= 40 and i > 20:
                num_20_40 += 1
            elif i <= 60 and i > 40:
                num_40_60 += 1
            elif i <= 80 and i > 60:
                num_60_80 += 1
            elif i > 80:
                num_80_inf += 1

        mw_dict =OrderedDict()
        mw_dict['<1'] = num_1
        mw_dict['1-5'] = num_1_5
        mw_dict['5-10'] = num_5_10
        mw_dict['10-20'] = num_10_20
        mw_dict['20-40'] = num_20_40
        mw_dict['40-60'] = num_40_60
        mw_dict['60-80'] = num_60_80
        mw_dict['>80'] = num_80_inf
        qc_dir = os.path.join(self.work_dir, 'tmp_qc')
        if os.path.exists(qc_dir):
            cover_file = os.path.join(qc_dir, 'Protein_seq_cover_distribution.xls')
            with open(cover_file, 'w') as ew:
                ew.write('coverage' + '\t' + 'number' + '\n')
                for c, n in mw_dict.items():
                    ew.write(c + '\t' + str(n) + '\n')

    def tmp_add_pepnum(self, pepnum_exp=None):
        df = pd.read_table(pepnum_exp, sep ="\t", header = 0)
        #num = dict(df.loc[:, '# Peptides'].value_counts())
        num = OrderedDict(df.loc[:, '# Peptides'].value_counts(ascending=True).sort_index())
        keys_1 = [str(int(x)) for x in num.keys()]
        values_1 = num.values()
        num = OrderedDict(zip(keys_1, num.values()))
        qc_dir = os.path.join(self.work_dir, 'tmp_qc')
        if os.path.exists(qc_dir):
            num_file = os.path.join(qc_dir, 'Peptite_number_distribution.xls')
            with open(num_file, 'w') as ew:
                ew.write('pep_num' + '\t' + 'number' + '\n')
                for l, n in num.items():
                    ew.write(l + '\t' + str(n) + '\n')

    def tmp_add_peperror(self,peperror_exp=None):
        df = pd.read_table(peperror_exp, sep ="\t", header = 0)
        try:
            mz = df.loc[:, 'm/z [Da]'].tolist()
            mz = [round(x, 6) for x in mz]
            delta = df.loc[:, 'DeltaM [ppm]'].tolist()
            delta = [round(x, 6) for x in delta]
            mz_delta = zip(mz, delta)
        except:
            mz_delta = list()
        if len(mz_delta) > 300000:
            mz_delta = mz_delta[0: 100000]
        # if len(mz_delta) > 15000:
        #     new_all_dot = [[int('%.f' % (i[0]/10))*10,float('%.1f' % i[1])] for i in mz_delta]
        #     new_new_all_dot = []
        #     for i in new_all_dot:
        #         if i not in new_new_all_dot:
        #             new_new_all_dot.append(i)
        #     mz_delta = [[i[0]+random.randint(-3, 3),i[1]] for i in new_new_all_dot]
        qc_dir = os.path.join(self.work_dir, 'tmp_qc')
        if os.path.exists(qc_dir):
            error_file = os.path.join(qc_dir, 'dMass.xls')
            with open(error_file, 'w') as ew:
                ew.write('m/z [Da]' + '\t' + 'DeltaM [ppm]' + '\n')
                for m, d in mz_delta:
                    ew.write(str(m) + '\t' + str(d) + '\n')

