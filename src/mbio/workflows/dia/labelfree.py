# -*- coding:utf-8 -*-
# __author__ = 'shicaiping'
"""labelfree工作流"""

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
import pandas as pd
from collections import OrderedDict
import pandas as pd
import numpy as np

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


class LabelfreeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        workflow option参数设置
        """
        self._sheet = wsheet_object
        super(LabelfreeWorkflow, self).__init__(wsheet_object)
        options = [
            ##选择输入文件
            {"name": "protein", "type": "infile", "format": "labelfree.common"},  # 输入蛋白鉴定表
            {"name": "psm", "type": "infile", "format": "labelfree.common"},  # 输入PSM表
            {"name": "peptide", "type": "infile", "format": "labelfree.common"},  # 输入肽段表
            {"name": "protein_information", "type": "infile", "format": "labelfree.common"},  # 输入蛋白信息表
            {"name": "ratio_exp", "type": "infile", 'format': "labelfree.ratio_exp"},  # 输入表达量ratio表
            #{"name": "scaled_exp", "type": "infile", "format": "labelfree.scaled_exp"},  # 输入表达量Abundance表
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
            {"name": "sub_loc", "type": "string", "default": ""}, # 亚细胞定位 物种分类
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

            ##蛋白互作网络
            {"name": "ppi_category", "type": "string", "default": "All"},# 设置页面端是否显示样本比较分析
            # {"name": "ppi_species", "type": "string", "default": "Homo ""sapiens"}, # 设置PPI物种
            {"name": "ppi_species", "type": "int", "default": 3702}, # 设置PPI物种

            ##是否是DIA数据
            {"name": "DIA", "type": "bool", "default": False},#如果选择DIA会传入这个参数，默认为False
            # {"name": "report_dia", "type": "infile", "format": "labelfree.common", "default": '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/report.xls'},#选择dia会需要上传的report文件
            # {"name": "exp_ana_dia", "type": "infile", "format": "labelfree.common", "default": '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/20190421_185436_SLD_ALL_SUM_ExperimentAnalysis.txt'},#选择dia会需要上传的实验配置文件，用来生成蛋白信息
            # {"name": "protein_dia", "type": "infile", "format": "labelfree.common", "default": '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/protein.xlsx'},#选择dia会需要上传的搜库文件
            {"name": "report_dia", "type": "string",
             "default": '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/report.xls'},
            # 选择dia会需要上传的report文件
            {"name": "exp_ana_dia", "type": "string",
             "default": '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/20190421_185436_SLD_ALL_SUM_ExperimentAnalysis.txt'},
            # 选择dia会需要上传的实验配置文件，用来生成蛋白信息
            {"name": "protein_dia", "type": "string",
             "default": '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/labelfree_dia/dia_qc/protein.xlsx'},
            # 选择dia会需要上传的搜库文件

            {"name": "change_des", "type": "bool", "default": True},#是否选择用nr注释结果替换description信息
            {"name": "useblast", "type": "bool", "default": True},#是否选择用string注释结果爬取官网
        ]

        #获取输出目录
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")
        self.project_sn = self._sheet.project_sn #获取project_sn
        self.task_id = self._sheet.id #获取task_id
        self.add_option(options)
        self.set_options(self._sheet.options())

        # 事先处理那3个表格，如果是excel，则转化成txt
        try:
            protein_df = pd.read_csv(self.option('protein').prop['path'])
        except:
            protein_df = pd.read_excel(self.option('protein').prop['path'], sep='\t')
            protein_path = os.path.join(os.path.dirname(self.option('protein').prop['path']), 'protein_.xls')
            protein_df.to_csv(protein_path,sep='\t',header=True,index=False)
            self.option('protein', protein_path)
        del protein_df

        if not self.option('DIA'):
            try:
                peptide_df = pd.read_csv(self.option('peptide').prop['path'])
            except:
                peptide_df = pd.read_excel(self.option('peptide').prop['path'], sep='\t')
                peptide_path = os.path.join(os.path.dirname(self.option('peptide').prop['path']), 'peptide_.xls')
                peptide_df.to_csv(peptide_path,sep='\t',header=True,index=False)
                self.option('peptide', peptide_path)
            del peptide_df

            try:
                psm_df = pd.read_csv(self.option('psm').prop['path'])
            except:
                psm_df = pd.read_excel(self.option('psm').prop['path'], sep='\t')
                psm_path = os.path.join(os.path.dirname(self.option('psm').prop['path']), 'psm_.xls')
                psm_df.to_csv(psm_path,sep='\t',header=True,index=False)
                self.option('psm', psm_path)
            del psm_df

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
        self.sam_pca = self.add_tool("labelfree.pca")
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
        self.final_tools = [self.searchdb, self.sam_corr, self.sam_pca, self.exp_venn, self.diff_cluster,
                            self.diff_go_class, self.diff_kegg_class, self.diff_ipath, self.diff_ppi, self.diff_circle]

        # 添加step，显示在页面进度条
        all_steps = ["filecheck", "searchdb", "exp_venn", "sam_corr", "sam_pca", "diff_pep",
                     "annotation", "diff_cluster", "diff_go_class",
                     "diff_kegg_class", "diff_ipath", "diff_ppi",
                     "diff_go_enrich", "diff_kegg_enrich", "diff_circle"]
        # self.step.add_steps("filecheck", "searchdb", "sam_corr", "sam_pca", "diff_pep","annotation")

        group_spname = self.option("protein_group").prop['group_dict']
        group_snum = [len(group_spname[g]) for g in group_spname]
        self.min_group_num = min(group_snum)
        # if self.min_group_num >= 3:
        #     self.final_tools[self.final_tools.index(self.sam_pca)] = self.ellipse
        #     all_steps.append("ellipse")

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
        self.filecheck.on('end', self.run_searchdb)
        if not self.option("protein_group").prop["sample_number"] < 3:
            self.filecheck.on('end', self.run_sam_pca)
            self.filecheck.on('end', self.run_sam_corr)
        self.filecheck.on('end', self.run_diff_pep)
        # self.filecheck.on('end', self.run_diamond)
        self.diff_pep.on('end', self.run_diamond)
        self.filecheck.on('end', self.run_exp_venn)
        self.on_rely([self.diff_kegg_enrich, self.diff_go_enrich], self.run_diff_circle)
        self.on_rely(self.final_tools, self.end)
        self.run_filecheck()
        super(LabelfreeWorkflow, self).run()


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

    def run_searchdb(self):
        self.logger.info("开始运行搜库")
        self.searchdb.set_options({
            'protein': self.option('protein').prop['path'],
            'ratio_exp': self.work_dir + "/treat_ref"
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
        ref_file = self.work_dir + "/treat_ref"
        df_scaled.to_csv(ref_file, sep = '\t')
        # 如果选择填充缺失值的话，这边会进行，为了少些点代码，这里直接在工作流里跑了，幸亏数据肯定不会很大
        if self.option('fillna').lower() not in ['no', 'none']:
            self.logger.info('选择了对表达量数据进行预处理，有效的蛋白的0值将被填充%s' % self.option('fillna').lower())
            from mbio.packages.labelfree.preprocess import Preprocess
            prep_obj = Preprocess(ref_file, self.option('protein_group').prop['path'],
                                  method=self.option('fillna'), cutoffs=self.option('cutoffs'), out=ref_file)
            prep_obj.fillzero()
            prep_obj.to_out()
        return ref_file

    def run_sam_corr(self):
        self.logger.info("开始运行样本相关性分析")
        self.sam_corr.set_options({'exp': self.ref_file})
        self.sam_corr.on('end', self.set_output, 'sam_corr')
        self.sam_corr.on('start', self.set_step, {'start': self.step.sam_corr})
        self.sam_corr.on('end', self.set_step, {'end': self.step.sam_corr})
        self.sam_corr.run()

    # def run_sam_pca(self):
    #     self.logger.info("开始运行样本PCA分析")
    #     # self.sam_pca.set_options({'otutable': self.ref_file})
    #     self.sam_pca.set_options({'exp': self.ref_file})
    #     self.sam_pca.on('end', self.set_output, 'sam_pca')
    #     if self.min_group_num >= 3:
    #         self.sam_pca.on('end', self.run_ellipse)
    #     self.sam_pca.on('start', self.set_step, {'start': self.step.sam_pca})
    #     self.sam_pca.on('end', self.set_step, {'end': self.step.sam_pca})
    #     self.sam_pca.run()

    def run_sam_pca(self):
        self.logger.info("开始运行样本PCA分析")
        self.sam_pca.set_options({'otutable': self.ref_file})
        self.sam_pca.on('end', self.set_output, 'sam_pca')
        self.sam_pca.on('start', self.set_step, {'start': self.step.sam_pca})
        self.sam_pca.on('end', self.set_step, {'end': self.step.sam_pca})
        self.sam_pca.run()

    def run_ellipse(self):
        self.logger.info("开始运行pca ellipse")
        opts = {
            'analysis': "pca",
            'pc_table': self.sam_pca.output_dir + "/PCA.xls",
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
            "express_matrix" : self.ref_file,
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
            'ratio_exp': self.ref_file,
            'method_type': self.option('method_type'),
            'correct_method': self.option('correct_method'),
            'pvalue': self.option('pvalue'),
            'fc_up': self.option('fc_up'),
            'fc_down': self.option('fc_down'),
            # 'cutoffs': self.option('cutoffs'),
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
        }
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            self.diamond_nr = self.add_module("labelfree.diamond")
            blast_opts.update(
                {
                    'database': self.option("nr_database"),
                    'evalue': self.option('go_evalue'),
                    'identity': self.option('go_identity')
                }
            )
            self.diamond_nr.set_options(blast_opts)
            self.blast_modules.append(self.diamond_nr)
            self.diamond_nr.on('end', self.set_output, 'diamond_nr')
            self.diamond_nr.run()
        if 'cog' in self.option('database'):
            self.diamond_string = self.add_module("labelfree.diamond")
            blast_opts.update(
                {
                    'database': 'string',
                    'evalue': self.option('cog_evalue'),
                    'identity': self.option('cog_identity')
                }
            )
            self.diamond_string.set_options(blast_opts)
            self.blast_modules.append(self.diamond_string)
            self.diamond_string.on('end', self.set_output, 'diamond_string')
            self.diamond_string.run()
        if 'kegg' in self.option('database'):
            self.diamond_kegg = self.add_module("labelfree.diamond")
            blast_opts.update(
                {'database': 'kegg',
                 'evalue': self.option('kegg_evalue'),
                 'kegg_species' : self.option('kegg_org'),
                 'identity': self.option('kegg_identity')
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
            }
            self.pfam.set_options(opts)
            self.pfam.on("end", self.set_output, "pfam")
            self.blast_modules.append(self.pfam)
            self.pfam.run()
        self.on_rely([self.diamond_nr, self.diamond_kegg, self.pfam, self.diamond_string], self.run_annotation)

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
            "pfam_domain" : self.pfam.output_dir + "/pfam_domain",
            'gram' : self.option('gram'),
            'sub_loc' : self.option('sub_loc'),
        }
        self.annotation.set_options(anno_opts)
        self.annotation.on('end', self.set_output, 'annotation')
        self.annotation.on('start', self.set_step, {'start': self.step.annotation})
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        class_prefix = [self.annotation, self.diff_pep]
        for func in [self.run_diff_go_class, self.run_diff_kegg_class,
                     self.run_diff_go_enrich, self.run_diff_kegg_enrich, self.run_diff_ipath, self.run_diff_ppi]:
            # self.on_rely(class_prefix, func)
            self.annotation.on('end', func)
            gevent.sleep(1)
        self.annotation.run()

    def run_diff_cluster(self):
        self.logger.info("开始运行差异蛋白cog分类注释")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "exp" : self.ref_file,
            "group" : self.option('protein_group')
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
        self.diff_cog_class.run()

    def run_diff_go_class(self):
        self.logger.info("开始运行差异蛋白go分类注释")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "go_stat" : os.path.join(self.annotation.output_dir, 'go', 'go12level_statistics.xls'),
        }
        self.diff_go_class.set_options(opts)
        self.diff_go_class.on('end', self.set_output, 'diff_go_class')
        self.diff_go_class.on('start', self.set_step, {'start': self.step.diff_go_class})
        self.diff_go_class.on('end', self.set_step, {'end': self.step.diff_go_class})
        self.diff_go_class.run()

    def run_diff_kegg_class(self):
        self.logger.info("开始运行差异蛋白kegg分类注释")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "png_dir" : os.path.join(self.annotation.output_dir, 'kegg','pathways'),
            "kegg_table" : os.path.join(self.annotation.output_dir, 'kegg','kegg_table.xls'),
            "pathway_table" : os.path.join(self.annotation.output_dir, 'kegg','pathway_table.xls'),
        }
        self.diff_kegg_class.set_options(opts)
        self.diff_kegg_class.on('end', self.set_output, 'diff_kegg_class')
        self.diff_kegg_class.on('start', self.set_step, {'start': self.step.diff_kegg_class})
        self.diff_kegg_class.on('end', self.set_step, {'end': self.step.diff_kegg_class})
        self.diff_kegg_class.run()

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
        self.diff_pfam_stat.run()

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
        self.diff_subloc_stat.run()

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
        self.diff_ipath.run()

    def run_diff_ppi(self):
        self.logger.info("开始运行差异蛋白ppi分析")
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
        self.diff_ppi.run()

    def run_diff_string_picture(self):
        self.logger.info("开始运行差异蛋白STRING数据库图片等的爬取")
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
        self.diff_string_picture.run()

    def run_diff_kegg_enrich(self):
        self.logger.info("开始运行差异蛋白kegg富集分析")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "kegg_table" : os.path.join(self.annotation.output_dir, 'kegg', 'kegg_table.xls'),
            "pathway_table" : os.path.join(self.annotation.output_dir, 'kegg', 'pathway_table.xls'),
        }
        self.diff_kegg_enrich.set_options(opts)
        self.diff_kegg_enrich.on('end', self.set_output, 'diff_kegg_enrich')
        self.diff_kegg_enrich.on('start', self.set_step, {'start': self.step.diff_kegg_enrich})
        self.diff_kegg_enrich.on('end', self.set_step, {'end': self.step.diff_kegg_enrich})
        self.diff_kegg_enrich.run()

    def run_diff_go_enrich(self):
        self.logger.info("开始运行差异蛋白go富集分析")
        opts = {
            "diff_path" : self.diff_pep.output_dir,
            "go_list" : os.path.join(self.annotation.output_dir, 'go', 'query_gos.list'),
        }
        self.diff_go_enrich.set_options(opts)
        self.diff_go_enrich.on('end', self.set_output, 'diff_go_enrich')
        self.diff_go_enrich.on('start', self.set_step, {'start': self.step.diff_go_enrich})
        self.diff_go_enrich.on('end', self.set_step, {'end': self.step.diff_go_enrich})
        self.diff_go_enrich.run()

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
        self.diff_circle.run()


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
            self.move2outputdir(obj.output_dir, '5_Proteinset/01_ProteinsetCluster')
        if event['data'] == 'diff_go_class':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_ProteinsetAnno/01_ProteinsetGO')
        if event['data'] == 'diff_kegg_class':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_ProteinsetAnno/02_ProteinsetKEGG')
        if event['data'] == 'diff_cog_class':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_ProteinsetAnno/03_ProteinsetCOG')
        if event['data'] == 'diff_pfam_stat':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_ProteinsetAnno/04_ProteinsetPfam')
        if event['data'] == 'diff_subloc_stat':
            self.move2outputdir(obj.output_dir, '5_Proteinset/03_ProteinsetAnno/05_ProteinsetSubLoc')
        if event['data'] == 'diff_go_enrich':
            self.move2outputdir(obj.output_dir, '5_Proteinset/04_ProteinsetEnrich/01_ProteinsetEnrichGO')
        if event['data'] == 'diff_kegg_enrich':
            self.move2outputdir(obj.output_dir, '5_Proteinset/04_ProteinsetEnrich/02_ProteinsetEnrichKEGG')
        if event['data'] == 'diff_circle':
            self.move2outputdir(obj.output_dir, '5_Proteinset/04_ProteinsetEnrich/03_ProteinsetEnrichChord')
        if event['data'] == 'diff_ppi':
            self.move2outputdir(obj.output_dir, '5_Proteinset/05_ProteinsetPPI')
        if event['data'] == 'diff_string_picture':
            self.move2outputdir(obj.output_dir, '5_Proteinset/05_ProteinsetStingPictures')
        if event['data'] == 'diff_ipath':
            self.move2outputdir(obj.output_dir, '5_Proteinset/06_ProteinsetIpath')

    def end(self):
        self.build_seq_database() # 创建序列数据库
        self.add_anno_to_diff() # 给差异表达文件增加注释
        self.run_api() # 运行导表函数

        ## 导表后，修改文件的绝对路径
        db = Config().get_mongo_client(mtype="labelfree")[Config().get_mongo_dbname("labelfree")]
        col1 = db["sg_annotation_stat"]
        col1.update({"task_id" : self.task_id}, {"$set": {"result_dir": self.workflow_output + "/Annotation"}}, upsert=True)
        col2 = db["sg_task"]
        col2.update({"task_id" : self.task_id}, {"$set": {"seq_db": self.workflow_output + "/SequenceDatabase/seq_db.sqlite3"}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"protein_sliced": self.workflow_output + "/SearchDB/protein_sliced.xls"}}, upsert=True)
        col2.update({"task_id" : self.task_id}, {"$set": {"protein_fa": self.workflow_output + "/ProteinFasta/protein.fa"}}, upsert=True)
        # 声明是V2
        col2.update({"task_id" : self.task_id}, {"$set": {"version": "V2"}}, upsert=True)
        self.modify_output() # 修改文件目录结构
        super(LabelfreeWorkflow, self).end()

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
        # dia项目暂时不要质控信息
        if not self.option('DIA'):
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
        # os.makedirs(tmp)
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'anno_stat'), os.path.join(tmp, '01_AnnoStat'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'go'), os.path.join(tmp, '02_AnnoGO'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'kegg'), os.path.join(tmp, '03_AnnoKEGG'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'cog'), os.path.join(tmp, '04_AnnoCOG'))
        shutil.copytree(os.path.join(target_dir, 'Annotation', 'subloc'), os.path.join(tmp, '06_SubLoc'))
        tmp = os.path.join(target_dir, '3_Annotation', '05_AnnoPfam')
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
        for dir in ['Express', 'SamCorr', 'SamPca', 'DiffPep', 'ExpVenn']:
            shutil.rmtree(os.path.join(target_dir, dir))

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
            ["Annotation", "", "蛋白注释结果目录"],
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
            ["DiffPep/diff_protein_summary_labelfree.xls", "", "所有组差异信息汇总表"],
            ["DiffPep/*.diff.xls", "", "样本间差异比较信息表"],
            ["DiffPep/*.diff.detail.xls", "", "样本间差异比较信息与蛋白详细注释信息表"],
            ["Express", "", "蛋白定量分析结果目录"],
            ["Express/exp.xls", "", "蛋白在各个样本间质谱峰面积表"],
            ["1_DataInfo", "", "搜库数据信息"],
            ["1_DataInfo/01_SearchDB", "", "搜库结果目录"],
            ["1_DataInfo/01_SearchDB/protein_sliced.xls", "", "搜库得到的蛋白详细信息"],
            ["1_DataInfo/01_SearchDB/search_db.xls", "", "搜库得到的蛋白详细信息以及在各个样本中的丰度"],
            ["2_SampleComp", "", "样本间比较分析"],
            ["2_SampleComp/01_SamPca", "", "样本PCA分析结果目录"],
            ["2_SampleComp/01_SamPca/pca_rotation.xls", "", "蛋白相关主成分贡献度"],
            ["2_SampleComp/01_SamPca/pca_sites.xls", "", "PCA分析样本坐标表"],
            ["2_SampleComp/01_SamPca/pca_importance.xls", "", "PCA主成分解释表"],
            ["2_SampleComp/02_SamCorr", "", "样本相关性分析结果目录"],
            ["2_SampleComp/02_SamCorr/sample_correlation.xls", "", "样本间相关性系数表"],
            ["2_SampleComp/02_SamCorr/expression_matrix.xls", "", "蛋白表达丰度矩阵表"],
            ["2_SampleComp/03_ExpVenn", "", "样本间venn分析结果目录"],
            ["2_SampleComp/03_ExpVenn/venn_graph.xls", "", "各样本有效蛋白列表"],
            ["3_Annotation", "", "蛋白注释结果目录"],
            ["3_Annotation/01_AnnoStat", "", "蛋白注释统计结果目录"],
            ["3_Annotation/01_AnnoStat/proteins_anno_detail.xls", "", "鉴定蛋白附加6数据库详细信息的表"],
            ["3_Annotation/01_AnnoStat/all_annotation_statistics.xls", "", "6个数据库中相关蛋白数以及总蛋白占比"],
            ["3_Annotation/01_AnnoStat/venn", "", "鉴定结果中蛋白在各个数据库注释结果统计VENN图目录"],
            ["3_Annotation/01_AnnoStat/venn/pfam_venn.txt", "", "有与pfam数据库相关信息的蛋白"],
            ["3_Annotation/01_AnnoStat/venn/kegg_venn.txt", "", "有与kegg数据库相关信息的蛋白"],
            ["3_Annotation/01_AnnoStat/venn/go_venn.txt", "", "有与go数据库相关信息的蛋白"],
            ["3_Annotation/01_AnnoStat/venn/cog_venn.txt", "", "有与cog数据库相关信息的蛋白"],
            ["3_Annotation/01_AnnoStat/venn/nr_venn.txt", "", "有与nr数据库相关信息的蛋白"],
            ["3_Annotation/01_AnnoStat/venn/subloc_venn.txt", "", "有与subloc数据库相关信息的蛋白"],
            ["3_Annotation/02_AnnoGO", "", "鉴定蛋白与GO数据库比对结果详情目录"],
            ["3_Annotation/02_AnnoGO/blast2go.annot", "", "鉴定蛋白与GO的对应关系以及GOterm名称"],
            ["3_Annotation/02_AnnoGO/query_gos.list", "", "鉴定蛋白与GO的对应关系"],
            ["3_Annotation/02_AnnoGO/go12level_statistics.xls", "", "鉴定蛋白的GO1、2层级的详细信息"],
            ["3_Annotation/02_AnnoGO/go123level_statistics.xls", "", "鉴定蛋白的GO1、2、3层级的详细信息"],
            ["3_Annotation/02_AnnoGO/go1234level_statistics.xls", "", "鉴定蛋白的GO1、2、3、4层级的详细信息"],
            ["3_Annotation/04_AnnoCOG", "", "鉴定蛋白与COG数据库比对结果详情目录"],
            ["3_Annotation/04_AnnoCOG/cog_table.xls", "", "鉴定蛋白序列文件与COG比对结果表"],
            ["3_Annotation/04_AnnoCOG/cog_summary.xls", "", "COG数据库注释结果汇总"],
            ["3_Annotation/04_AnnoCOG/cog_list.xls", "", "鉴定蛋白与COG的对应关系"],
            ["3_Annotation/03_AnnoKEGG", "", "鉴定蛋白与KEGG数据库比对结果详情目录"],
            ["3_Annotation/03_AnnoKEGG/pathways", "", "KEGG通路注释结果图片目录"],
            ["3_Annotation/03_AnnoKEGG/pathways/map*.png", "", "通路注释结果图片-位图"],
            ["3_Annotation/03_AnnoKEGG/pathways/map*.pdf", "", "通路注释结果图片-矢量图"],
            ["3_Annotation/03_AnnoKEGG/pathway_table.xls", "", "每行以通路为单位显示通路注释详情表"],
            ["3_Annotation/03_AnnoKEGG/kegg_table.xls", "", "每行以蛋白为单位显示通路注释详情表"],
            ["3_Annotation/03_AnnoKEGG/pid.txt", "", "每个注释通路和映射上去的蛋白联系表"],
            ["3_Annotation/06_SubLoc", "", "蛋白亚细胞定位注释结果目录"],
            ["3_Annotation/06_SubLoc/multiloc_stat.xls", "", "亚细胞位置对应蛋白的数量表"],
            ["3_Annotation/06_SubLoc/multiloc.xls", "", "亚细胞定位结果表"],
            ["3_Annotation/05_AnnoPfam/", "", "鉴定蛋白与Pfam数据库比对结果文件目录"],
            ["3_Annotation/05_AnnoPfam/pfam_domain", "", "鉴定蛋白与pfam数据库比对结果xml形式"],
            ["4_ExpDiff", "", "全部样本差异表达详情目录"],
            ["4_ExpDiff/diff_protein_summary_labelfree.xls", "", "所有组差异信息汇总表"],
            ["4_ExpDiff/*.diff.xls", "", "样本间差异比较信息表"],
            ["4_ExpDiff/*.diff.detail.xls", "", "样本间差异比较信息与蛋白详细注释信息表"],
            ["1_DataInfo/03_ExpInfo", "", "蛋白定量分析结果目录"],
            ["1_DataInfo/03_ExpInfo/exp.xls", "", "蛋白在各个样本间质谱峰面积表"],
            ["5_Proteinset", "", "差异蛋白相关分析结果目录"],
            ["5_Proteinset/01_ProteinsetCluster", "", "差异蛋白聚类分析结果目录"],
            ["5_Proteinset/03_ProteinsetAnno", "", "差异蛋白分类注释相关分析结果目录"],
            ["5_Proteinset/03_ProteinsetAnno/01_ProteinsetGO", "", "差异蛋白GO分类注释分析结果目录"],
            ["5_Proteinset/03_ProteinsetAnno/02_ProteinsetKEGG", "", "差异蛋白kegg分类注释分析结果目录"],
            ["5_Proteinset/03_ProteinsetAnno/03_ProteinsetCOG", "", "差异蛋白cog分类注释分析结果目录"],
            ["5_Proteinset/03_ProteinsetAnno/04_ProteinsetPfam", "", "差异蛋白Pfam分类注释分析结果目录"],
            ["5_Proteinset/03_ProteinsetAnno/05_ProteinsetSubLoc", "", "差异蛋白亚细胞定位分类注释分析结果目录"],
            ["5_Proteinset/04_ProteinsetEnrich", "", "差异蛋白富集相关分析结果目录"],
            ["5_Proteinset/04_ProteinsetEnrich/01_ProteinsetEnrichGO", "", "差异蛋白GO富集分析结果目录"],
            ["5_Proteinset/04_ProteinsetEnrich/02_ProteinsetEnrichKEGG", "", "差异蛋白KEGG富集分析结果目录"],
            ["5_Proteinset/04_ProteinsetEnrich/03_ProteinsetEnrichChord", "", "差异蛋白弦图结果目录"],
            ["5_Proteinset/05_ProteinsetPPI", "", "差异蛋白蛋白互作网络结果目录"],
            ["5_Proteinset/05_ProteinsetStingPictures", "", "差异蛋白String数据库爬虫结果目录"],
            ["5_Proteinset/06_ProteinsetIpath", "", "差异蛋白Ipath分析结果目录"],
        ]
        sdir = self.add_upload_dir(target_dir)
        sdir.add_relpath_rules(repaths)

    def run_api(self):
        import traceback
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.logger.info("导表开始")
        self.stop_timeout_check()
        task_info = self.api.api('task_info.labelfree_task_info')
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
        self.logger.info("开始导差异蛋白分析相关表，如果错误会尽量忽略")
        # try:
        #     self.export_diff_cluster()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_cog_class()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_go_class()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_kegg_class()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_pfam_stat()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_subloc_stat()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_go_enrich()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_kegg_enrich()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_circle()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_ppi()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # try:
        #     self.export_diff_ipath()
        # except:
        #     estr = traceback.format_exc()
        #     self.logger.info(estr)
        # self.export_diff_cluster()
        # self.export_diff_cog_class()
        self.export_diff_go_class()
        self.export_diff_kegg_class()
        # self.export_diff_pfam_stat()
        # self.export_diff_subloc_stat()
        self.export_diff_go_enrich()
        self.export_diff_kegg_enrich()
        self.export_diff_ppi()
        # self.export_diff_string_picture()
        self.export_diff_circle()
        self.export_diff_ipath()
        # super(LabelfreeWorkflow, self).end()
        self.logger.info("导表完成")

    @time_count
    def export_specimen(self):
        gevent.sleep()
        self.api_specimen = self.api.api("labelfree.specimen")
        self.api_specimen.add_specimen(group_txt=self.option("protein_group").prop["path"], params=None,
                project_sn=self.project_sn, task_id=self.task_id, type="labelfree")
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
        if self.new:
            self.api_searchdb = self.api.api("itraq_and_tmt.searchdb")
            self.api_searchdb._project_type = 'labelfree'
        else:
            self.api_searchdb = self.api.api("labelfree.searchdb")
        self.api_searchdb.add_searchdb(project_sn=self.project_sn,
                task_id=self.task_id, name=None, protein_selected = protein_selected)

    @time_count
    def export_sam_corr(self):
        gevent.sleep()
        self.api_sam_corr = self.api.api("labelfree.all_exp")
        params = dict(
            task_id=self.task_id,
            submit_location='expresscorr',
            group_id=str(self.group_id),
            group_dict=self.option('protein_group').prop['group_dict'],
            scm="complete",
            scd="euclidean",
            corr_method="pearson",
            task_type=2
        )
        corr_output = self.sam_corr.work_dir
        self.api_sam_corr.add_exp_corr3(corr_output, params=params, main_id=None,
                      project_sn=self.project_sn, task_id=self.task_id)

    @time_count
    def export_qc(self):
        gevent.sleep()
        if self.option('DIA'):
            self.api_dia = self.api.api("labelfree.dia_quality_control")
            params = {"software": "DIA"}
            # self.api_dia.add_peplen(self.option('report_dia').prop['path'], params=params)
            # self.api_dia.add_pepnum(self.option('report_dia').prop['path'], params=params)
            # self.api_dia.add_proteininfo(self.option('exp_ana_dia').prop['path'], params=params)
            # self.api_dia.add_proteinmw(self.option('protein_dia').prop['path'], params=params)
            # self.api_dia.add_coverage(self.option('protein_dia').prop['path'], params=params)
            self.api_dia.add_peplen(self.option('report_dia').prop['path'], params=params)
            self.api_dia.add_pepnum(self.option('report_dia').prop['path'], params=params)
            self.api_dia.add_proteininfo(self.option('exp_ana_dia').prop['path'], params=params)
            # self.api_dia.add_proteinmw(self.option('protein_dia'), params=params)
            self.api_dia.add_proteinmw(self.option('protein').prop['path'], params=params)
            # self.api_dia.add_coverage(self.option('protein_dia'), params=params)
            self.api_dia.add_coverage(self.option('protein').prop['path'], params=params)
            self.api_express = self.api.api("labelfree.express")
            self.api_express.add_express(params=params, main_id=None,
                                         project_sn=self.project_sn, task_id=self.task_id, type='ratio',
                                         express_exp=self.work_dir + "/treat_ref")
            self.api_express.add_express(params=params, main_id=None,
                                         project_sn=self.project_sn, task_id=self.task_id, type='scaled',
                                         express_exp=self.work_dir + "/treat_ref")
            return
        if self.new:
            self.api_coverage = self.api.api("itraq_and_tmt.coverage")
            self.api_coverage._project_type = 'labelfree'
            self.api_express = self.api.api("labelfree.express")
            self.api_peperror = self.api.api("itraq_and_tmt.peperror")
            self.api_peperror._project_type = 'labelfree'
            self.api_peplen = self.api.api("itraq_and_tmt.peplen")
            self.api_peplen._project_type = 'labelfree'
            self.api_pepnum = self.api.api("itraq_and_tmt.pepnum")
            self.api_pepnum._project_type = 'labelfree'
            self.api_proteininfo = self.api.api("labelfree.proteininfo")
            self.api_proteinmw = self.api.api("itraq_and_tmt.proteinmw")
            self.api_proteinmw._project_type = 'labelfree'
        else:
            self.api_coverage = self.api.api("labelfree.coverage")
            self.api_express = self.api.api("labelfree.express")
            self.api_peperror = self.api.api("labelfree.peperror")
            self.api_peplen = self.api.api("labelfree.peplen")
            self.api_pepnum = self.api.api("labelfree.pepnum")
            self.api_proteininfo = self.api.api("labelfree.proteininfo")
            self.api_proteinmw = self.api.api("labelfree.proteinmw")
        params = {"software": "peaks 8.5"}
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        corr_output = self.sam_corr.work_dir
        self.api_coverage.add_coverage(params=params, main_id=None,
                                       project_sn=self.project_sn, task_id=self.task_id,
                                       coverage_exp=self.option("protein").prop['path'])
        self.api_express.add_express(params=params, main_id=None,
                                     project_sn=self.project_sn, task_id=self.task_id, type='ratio',
                                     express_exp=self.work_dir + "/treat_ref")
        self.api_express.add_express(params=params, main_id=None,
                                     project_sn=self.project_sn, task_id=self.task_id, type='scaled',
                                     express_exp=self.work_dir + "/treat_ref")
        self.api_peplen.add_peplen(params=params, main_id=None,
                                   project_sn=self.project_sn, task_id=self.task_id,
                                   peplen_exp=self.option("peptide").prop['path'])
        self.api_pepnum.add_pepnum(params=params, main_id=None,
                                   project_sn=self.project_sn, task_id=self.task_id,
                                   pepnum_exp=self.option("protein").prop['path'])
        self.api_proteininfo.add_proteininfo(params=params, main_id=None,
                                             project_sn=self.project_sn, task_id=self.task_id,
                                             proteininfo_exp=self.option("protein_information").prop['path'])
        self.api_proteinmw.add_proteinmw(params=params, main_id=None,
                                             project_sn=self.project_sn, task_id=self.task_id,
                                             proteinmw_exp=self.option("protein").prop['path'])

    # @time_count
    # def export_sam_pca(self):
    #     gevent.sleep()
    #     # self.api_sam_pca = self.api.api("labelfree.pca")
    #     all_exp = self.api.api("labelfree.all_exp")
    #     params = dict(
    #         task_id=self.task_id,
    #         submit_location="expresspca",
    #         group_id=str(self.group_id),
    #         group_dict=self.option('protein_group').prop['group_dict'],
    #         task_type=2
    #     )
    #     #os.remove(self.sam_pca.output_dir + "/pca_rotation.xls")
    #     #os.rename(self.sam_pca.output_dir + "/pca_rotation_all.xls", self.sam_pca.output_dir + "/pca_rotation.xls")
    #     # pca_output = self.sam_pca.output_dir
    #     # self.api_sam_pca.add_pca(pca_output, exp_id=None, params=params,
    #     #             project_sn=self.project_sn, task_id=self.task_id, main_id=None)
    #     pca_output = self.sam_pca.work_dir
    #     pca_main_id = all_exp.add_exp_pca2(pca_output, params=params,
    #                               project_sn=self.project_sn, task_id=self.task_id, main_id=None)
    #     if self.min_group_num >= 3:
    #         all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(pca_main_id))

    @time_count
    def export_sam_pca(self):
        gevent.sleep()
        self.api_sam_pca = self.api.api("labelfree.pca")
        params = dict(
            task_id=self.task_id,
            submit_location="expresspca",
            group_id=str(self.group_id),
            group_dict=self.option('protein_group').prop['group_dict'],
            task_type=2
        )
        # os.remove(self.sam_pca.output_dir + "/pca_rotation.xls")
        # os.rename(self.sam_pca.output_dir + "/pca_rotation_all.xls", self.sam_pca.output_dir + "/pca_rotation.xls")
        pca_output = self.sam_pca.output_dir
        self.api_sam_pca.add_pca(pca_output, exp_id=None, params=params,
                                 project_sn=self.project_sn, task_id=self.task_id, main_id=None)

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
        api_exp_venn = self.api.api("labelfree.all_exp")
        params = dict(
            task_id=self.task_id,
            submit_location="expvenn",
            group_id=str(group_id),
            group_dict=group_dict,
            task_type=2
        )
        graph_table = os.path.join(self.exp_venn.output_dir, 'venn_graph.xls')
        api_exp_venn.add_exp_venn(graph_table, params=params, main_id=None)

    @time_count
    def export_pep_diff(self):
        gevent.sleep()
        self.api_diff = self.api.api("labelfree.diff")

        params = dict(
            task_id=self.task_id,
            submit_location="diff",
            group_id=str(self.group_id),
            group_dict=self.option('protein_group').prop['group_dict'],
            type= "origin", correct_method=self.option('correct_method'), fc_up=self.option('fc_up'), fc_down=self.option('fc_down'),
            pvalue=self.option('pvalue'), diff_method=self.option('method_type'),control_id=str(self.control_id), task_type=2)

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        work_dir = self.diff_pep.work_dir
        compare_dict_xls = work_dir + "/compare_dict.xls"
        num_summary_xls = work_dir + "/num_summary.xls"
        allsummary_xls = work_dir + "/all_labelfree_summary.xls"
        self.api_diff.add_diff(work_dir, compare_dict_xls=compare_dict_xls, num_summary_xls=num_summary_xls,
            allsummary_xls=allsummary_xls,params=params, project_sn=self.project_sn,
            task_id=self.task_id, method_type=self.option("method_type"))

    @time_count
    def export_annotation(self):
        '''
        导入注释结果表格 liubinxu
        '''
        gevent.sleep()
        if self.new:
            self.annotation_api = self.api.api("labelfree.protein_annotation")
        else:
            self.annotation_api = self.api.api("labelfree.protein_annotation")
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
        self.export_proteinset = self.api.api("labelfree.proteinset")
        self.export_proteinset.add_proteinset(diff_work_dir=self.diff_pep.work_dir, group_id=str(self.group_id),
                                              task_id=self.task_id, project_sn=self.project_sn)

    @time_count
    def export_diff_cluster(self):
        self.logger.info("导入差异蛋白聚类结果")
        all_exp = self.api.api("labelfree.all_exp")
        cluster_out = self.diff_cluster.output_dir
        for dir in os.listdir(cluster_out):
            if os.path.isdir(os.path.join(cluster_out, dir)):
                all_exp.add_geneset_cluster(os.path.join(cluster_out, dir), main_id='')

    @time_count
    def export_diff_cog_class(self):
        self.logger.info("导入差异蛋白cog分类注释")
        api_proteinset = self.api.api('labelfree.proteinset')
        for cog_result in glob.glob(os.path.join(self.diff_cog_class.output_dir, '*', '*cog_class_table.xls')):
            api_proteinset.add_proteinset_cog_detail(cog_result, '')

    @time_count
    def export_diff_go_class(self):
        self.logger.info("导入差异蛋白go分类注释")
        api_proteinset = self.api.api('labelfree.proteinset')
        for go_result in glob.glob(os.path.join(self.diff_go_class.output_dir, '*', '*go_class_table.xls')):
            api_proteinset.add_go_regulate_detail(go_result, '')
            if '_up_down' not in go_result:
                continue
            api_proteinset.add_go_regulate_detail2(go_result, '')

    @time_count
    def export_diff_pfam_stat(self):
        self.logger.info("导入差异蛋白pfam分类注释")
        api_proteinset = self.api.api('labelfree.proteinset')
        for pfam_stat in glob.glob(os.path.join(self.diff_pfam_stat.output_dir, '*', '*stat.xls')):
            api_proteinset.add_pfam_stat('', pfam_stat)

    @time_count
    def export_diff_subloc_stat(self):
        self.logger.info("导入差异蛋白subloc分类注释")
        api_proteinset = self.api.api('labelfree.proteinset')
        for subloc_stat in glob.glob(os.path.join(self.diff_subloc_stat.output_dir, '*', '*stat.xls')):
            api_proteinset.add_subloc_stat('', subloc_stat)

    @time_count
    def export_diff_kegg_class(self):
        self.logger.info("导入差异蛋白kegg分类注释")
        api_proteinset = self.api.api('labelfree.proteinset')
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
        api_proteinset = self.api.api('labelfree.proteinset')
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
        api_proteinset = self.api.api('labelfree.proteinset')
        kegg_out = self.diff_kegg_enrich.output_dir
        for dir in os.listdir(kegg_out):
            if os.path.isdir(os.path.join(kegg_out, dir)):
                enrich_file = glob.glob(os.path.join(kegg_out, dir, '*.xls'))[0]
                api_proteinset.add_kegg_enrich_detail('', enrich_file)

    @time_count
    def export_diff_circle(self):
        self.logger.info("导入差异蛋白kegg富集炫图和GO富集炫图结果")
        api_proteinset = self.api.api('labelfree.proteinset')
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
        api_proteinset = self.api.api('labelfree.proteinset')
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
        api_ppinetwork = self.api.api('labelfree.proteinset_ppi')
        ppi_out = self.diff_ppi.output_dir
        for dir in os.listdir(ppi_out):
            if os.path.isdir(os.path.join(ppi_out, dir)):
                all_nodes_path = os.path.join(ppi_out, dir, 'ppinetwork_predict', 'all_nodes.txt')  # 画图节点属性文件
                interaction_path = os.path.join(ppi_out, dir, 'ppinetwork_predict', 'interaction_detail.txt')  # 画图的边文件
                network_stats_path = os.path.join(ppi_out, dir, 'ppinetwork_predict', 'network_stats.txt')  # 网络全局属性统计
                network_centrality_path = os.path.join(ppi_out, dir, 'ppinetwork_topology', 'protein_interaction_network_centrality.txt')
                degree_distribution_path = os.path.join(ppi_out, dir, 'ppinetwork_topology', 'protein_interaction_network_degree_distribution.txt')
                table_id = api_ppinetwork.add_node_table(file_path=all_nodes_path,
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
        api_proteinset = self.api.api('labelfree.proteinset')
        string_out = self.diff_string_picture.output_dir
        for dir in os.listdir(string_out):
            if os.path.isdir(os.path.join(string_out, dir)):
                api_proteinset.add_string_picture('', os.path.join(string_out, dir))

    def build_seq_database(self):
        self.logger.info("创建蛋白数据库")
        self.export_seq = self.api.api("labelfree.protein_seq")
        pep = self.option("protein_fasta").prop["path"]
        seq_db = os.path.join(self.work_dir, 'seq_db.sqlite3')
        self.export_seq.build_seq_database(seq_db, pep, task_id=self.task_id)

    def add_anno_to_diff(self):
        self.logger.info("给差异蛋白文件添加注释")
        protein_anno = self.annotation.output_dir + '/anno_stat/proteins_anno_detail.xls'
        protein_anno_pd = pd.read_table(protein_anno, header=0, index_col=0)
        diff_output = self.diff_pep.output_dir
        target_files = glob.glob(diff_output + '/' + '*_vs_*_diff.xls')
        for each in target_files:
            protein_pd = pd.read_table(each, header=0, index_col=0)
            protein_anno_result = pd.concat([protein_pd, protein_anno_pd], axis=1)
            # protein_anno_out = each.split('.xls')[0] + '.detail.xls'
            protein_anno_out = self.output_dir + '/DiffPep/' + os.path.basename(each).split('.xls')[0] + '.detail.xls'
            protein_anno_result.to_csv(protein_anno_out, header=True, index=True, sep='\t')
