# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# __last_modified__ = '20180301'

"""宏基因组分析工作流"""
# 不质控暂无法完成导表(缺少base_info碱基统计文件，无法导入raw_data)

from mainapp.libs.param_pack import group_detail_sort
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.id_convert import name2id, id2name
from mbio.packages.metagenomic.common import get_old_mongo, get_mongo
from mainapp.models.mongo.metagenomic import Metagenomic
from mbio.packages.metagenomic.copy_demo import CopyDemo
from mbio.packages.meta.delete_mongo import DeleteDemoMongo
from bson import ObjectId
import os, re
import glob
import json
import shutil
import time
import datetime
import gevent
import functools


def time_count(func):  # 统计导表时间
    @functools.wraps(func)
    def wrapper(self, *args, **kw):
        start = time.time()
        func_name = func.__name__
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

def get_insertsize(file_path, output):
    """
    生成insert_size文件至workspace下
    :param file_path: 原始数据文件路径
    :return: insert_size文件路径
    """
    if os.path.isdir(output):
        output = os.path.join(output, 'insertsize')
    fw = open(output, 'w')
    with open(file_path, 'rb') as f:
        lines = f.readlines()
        for line in lines[1:]:
            line = line.strip().split('\t')
            fw.write(line[0] + "\t" + line[2] + "\n")
    fw.close()
    return output


class MetaGenomicWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        宏基因组workflow option参数设置
        """
        self._sheet = wsheet_object
        super(MetaGenomicWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'test', 'type': 'bool', 'default': False},  # 是否为测试workflow
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'speciman_info', 'type': 'infile', 'format': 'meta_genomic.specimen_info'},  # 样本集信息表
            {'name': 'raw_info', 'type': 'infile', 'format': 'sequence.profile_table'},  # 原始序列的信息表
            {'name': 'qc_info', 'type': 'infile', 'format': 'sequence.profile_table'},  # 质控后的信息表
            {'name': 'qc', 'type': 'bool', 'default': False},  # 是否需要质控
            {'name': 'insertsize', 'type': 'infile', 'format': 'sample.insertsize_table'},  # 插入片段长度表
            {'name': 'qc_quality', 'type': 'int', 'default': 20},  # 质控质量值标准
            {'name': 'qc_length', 'type': 'int', 'default': 30},  # 质控最短序列长度
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            # 未定义的宿主序列所在文件夹，多个宿主cat到一个文件，并作为tool:align.bwa的输入文件，可不提供
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            # {'name': 'assemble_tool', 'type': 'string', 'default': 'idba'},  # 选择拼接工具，soapdenovo OR idba
            {'name': 'assemble_type', 'type': 'string', 'default': 'Megahit'},
            # 选择拼接策略，soapdenovo OR idba OR megahit OR multiple
            # 选择拼接策略，SOAPdenovo2 OR IDBA_UD OR Megahit OR Multiple OR Multiple_Megahit OR Multiple_IDBA_UD
            {'name': 'use_newbler', 'type': 'string', 'default': 'False'}, # 增加use_newbler参数，控制拼接模块是否调用newbler拼接 by ghd @ 20181023
            {'name': 'min_contig', 'type': 'int', 'default': 300},  # 拼接序列最短长度
            {'name': 'gene_predictor', 'type': 'string', 'default': "prodigal"},  # 基因预测软件 prodigal metagene metagenemaker
            {'name': 'min_gene', 'type': 'int', 'default': 100},  # 预测基因最短长度
            {'name': 'cdhit_identity', 'type': 'float', 'default': 0.95},  # 基因序列聚类相似度
            {'name': 'cdhit_coverage', 'type': 'float', 'default': 0.9},  # 基因序列聚类覆盖度
            {'name': 'soap_identity', 'type': 'float', 'default': 0.95},  # 基因丰度计算相似度
            {'name': 'anno_nr', 'type': 'string', 'default': ''},  # nr 注释使用e值，不注释为空
            #注释使用e值，不注释为空， 逗号分割，依次对应 'kegg,cog,cazy,ardb,card,vfdb'
            {'name': 'anno_list', 'type': 'string', 'default': ''},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},  # 物种/功能分析输入环境因子表
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 物种/功能分析输入group表
            {'name': 'personal_anno', 'type': 'string', 'default': 'none'}, ##个性化注释模块调用参数，如果为None不做
            {"name": "unique_fa", "type": "infile", "format": "sequence.fasta"},  # 工作流2使用的非冗余基因集fa
            {'name': 'clean_list', 'type': 'infile', 'format': 'sequence.profile_table'}, # 工作流2使用的clean_data的list.txt
            {'name': 'fastp', 'type': 'string', 'default': 'fastp'}, #宏基因组质控 add by qingchen.zhang 20190507
            {'name': 'old_task_id', 'type': 'string', 'default': ''},  # 工作流2通过老项目task_id拉取数据
            {'name': 'old_task_db', 'type': 'int', 'default': 1},  # 老项目tasks所在的mongo服务器版本
            {'name': 'pipeline', 'type': 'string', 'default': '1'},  # 指定流程类型 1，2，2.1，3
            # metaphlan3.0 参数
            {'name': 'metaphlan', 'type': 'bool', 'default': False},
            {'name': 'mph_min_cu_len', 'type': 'float', 'default': 0.1},
            {'name': 'mph_stat', 'type': 'string', 'default': 'tavg_g'},
            {'name': 'mph_stat_q', 'type': 'float', 'default': 0.2},
            # kraken2.0 参数
            {'name': 'kraken', 'type': 'bool', 'default': False},
            {'name': 'kk_confidence', 'type': 'float', 'default': 0.1},
            # 基因组物种与功能注释 evalue 设置
            {'name': 'nr_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'cog_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'kegg_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'cazy_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'ardb_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'card_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'vfdb_evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'save_pdf', 'type': 'int', 'default': 0}
        ]
        self.wf_version = "3.1"
        self.add_option(options)
        self.set_options(self._sheet.options())
        '''初始化module/tool'''
        #self.sequence = self.add_module('sequence.meta_genomic')
        if self.option("gene_predictor") in ["metagene", "MetaGene"]:
            self.gene_predict = self.add_module('gene_structure.gene_predict')
        else:
            self.gene_predict = self.add_module("gene_structure.metag_predict")
        self.qc = self.add_module('meta.qc.qc_and_stat')
        self.rm_host = self.add_module('meta.qc.bwa_remove_host')
        self.assem_soapdenovo = self.add_module('assemble.mg_ass_soapdenovo')
        self.assem_idba = self.add_module('assemble.mg_ass_idba')
        self.gene_set = self.add_module('cluster.uni_gene')
        self.nr = self.add_module('align.meta_diamond')
        self.cog = self.add_module('align.meta_diamond')
        self.kegg = self.add_module('align.meta_diamond')
        self.anno = self.add_module('annotation.mg_common_anno_stat')
        self.cazy = self.add_module('annotation.cazy_annotation')
        self.ardb = self.add_module('annotation.ardb_annotation')
        self.card = self.add_module('annotation.card_annotation')
        self.vfdb = self.add_module('annotation.vfdb_annotation')
        self.overview = self.add_tool('annotation.database_merge')
        self.table = dict()
        self.table['check'] = self.add_tool('meta.create_abund_table')
        self.composition = self.add_module('meta.composition.composition_analysis')
        self.compare = self.add_module('meta.beta_diversity.beta_diversity')
        self.correlation = self.add_tool('statistical.pearsons_correlation')
        #### version新增模块
        self.anno_nr_lca = self.add_module('annotation.mg_common_anno_stat')
        self.anno_nr_deunclass = self.add_module('annotation.mg_common_anno_stat')
        self.wf_soft = []  # 工作流中使用导的软件列表
        self.wf_database = []  # 工作起来流程使用到的数据库列表
        self.taxon_outset = {}

        ####
        '''add_steps'''
        if self.pipeline == 'pipeline3':
            self.sequence = self.add_module('metagenome.input_process')
            self.step.add_steps('sequence')
        elif self.pipeline == 'pipeline2.1':
            self.step.add_steps('copy_data', 'nr_', 'cog', 'kegg','anno', 'cazy', 'vfdb', 'ardb', 'card',
                                'overview', 'table', 'composition', 'compare','correlation')
        elif self.pipeline == "pipeline2":
            self.sequence_deal = self.add_module('sequence.metagenomic_check')
            self.step.add_steps('sequence', 'gene_set', 'nr_', 'cog', 'kegg','anno', 'cazy', 'vfdb', 'ardb', 'card',
                                'overview', 'table', 'composition', 'compare','correlation')
        elif self.pipeline == "pipeline1":
            # self.sequence = self.add_module('sequence.meta_genomic')
            self.sequence = self.add_module("metagenome.input_process")
            self.step.add_steps('sequence', 'qc_', 'rm_host', 'assem', 'gene_predict', 'gene_set', 'nr_', 'cog', 'kegg',
                            'anno', 'cazy', 'vfdb', 'ardb', 'card', 'overview', 'table', 'composition', 'compare',
                            'correlation')
        if self.option("personal_anno") != "none":
            self.step.add_steps('sequence', 'qc_', 'rm_host', 'assem', 'gene_predict', 'gene_set', 'nr_', 'cog', 'kegg',
                            'anno', 'cazy', 'vfdb', 'ardb', 'card', 'overview', 'table', 'composition', 'compare',
                            'correlation','anno_personal')
        self.logger.info("test pipeline>>>>>>>>>>>>>>>")
        '''初始化自定义变量'''
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.api_dic = {}  # 存放導表函數
        self.anno_tool = []  # nr/kegg/cog注释记录
        self.all_anno = []  # 全部的注释记录(用于依赖关系)
        self.new_table = []  # 构建新丰度表模块(module)
        self.table_pointer = 0  # 查看对应的丰度表运行任务输出路径
        self.analysis = []  # 分析模块具体分析内容(module/tool)
        self.nr_dir = ''  # nr注释结果文件路径，导表时用
        self.cog_dir = ''
        self.kegg_dir = ''
        self.anno_table = dict()  # 注释结果表(含所有注释水平，含丰度结果表)
        self.profile_table1 = dict()  # 注释丰度表(用于组成分析，相关性heatmap图)
        self.profile_table2 = dict()  # 注释丰度表(用于样品比较分析、rda、cca、db_rda分析)
        self.soft_db_info = {}
        self.default_level1 = {
            'nr': 'Genus',
            'cog': 'Function',
            'kegg': 'Level1',
            'cazy': 'Class',
            'vfdb': 'Level1',
            'ardb': 'Type',
            'card': 'ARO',
        }
        self.default_level2 = {
            'nr': 'Genus',
            'cog': 'NOG',
            'kegg': 'Level3',
            'cazy': 'Family',
            'vfdb': 'VFs',
            'ardb': 'ARG',
            'card': 'ARO',
        }
        # self.composition_dir2anno = {}  # 输出结果和导表时根据此值判断数据库类型
        # self.compare_dir2anno = {}
        # self.correlation_dir2anno = {}
        # self.anno2correlation_tree = {}
        self.sample_in_group = []  # 根据group文件获取样本名
        self.specimen_group = ""  # group_id
        self.group_detail = {}  # 分组对应样品id{group1: [id1,id2,id3], group2: [id4,id5,id6]}
        self.env_id = ""  # 环境因子的id
        self.env_labs = ""  # 环境因子，","分割
        self.geneset_id = ""  # geneset主表_id
        self.spname_spid = {}  # 样本对应原始数据表_id
        self.anno_id = {}  # 注释表对应mongo数据库主表_id
        # sanger_type, sanger_path = self._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.remote_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # sanger_prefix = Config().get_project_region_bucket(project_type="metagenomic")
        # self.remote_dir = os.path.join(sanger_prefix, self._sheet.output.split(":")[-1].lstrip("/"))  # region://bucket/files/...
        self.remote_dir = self._sheet.output
        self.logger.info(">>>>>>>>>>>>>>>>>>>>>>>>><<<test workflow new version<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        if self.option('test'):
            '''
            self.gene_set.option('uni_fasta').set_path(
                "/mnt/ilustre/users/sanger-dev/workspace/20171212/MetaGenomic_metagenome_debug6/output/geneset/uniGeneset/gene.uniGeneset.fa")
            self.gene_set.option('uni_fastaa').set_path(
                "/mnt/ilustre/users/sanger-dev/workspace/20171212/MetaGenomic_metagenome_debug6/output/geneset/uniGeneset/gene.uniGeneset.faa")
            self.gene_set.option('reads_abundance').set_path(
                "/mnt/ilustre/users/sanger-dev/workspace/20171212/MetaGenomic_metagenome_debug6/output/geneset/gene_profile/reads_number.xls")

            self.anno_table = {
                'geneset': '/mnt/ilustre/users/sanger-dev/workspace/20171212/MetaGenomic_metagenome_debug6/output/geneset/gene_profile/RPKM.xls',
                'ardb': '/mnt/ilustre/users/sanger-dev/workspace/20171212/MetaGenomic_metagenome_debug6/output/ardb/gene_ardb_anno.xls',
                'card': '/mnt/ilustre/users/sanger-dev/workspace/20171212/MetaGenomic_metagenome_debug6/output/card/gene_card_anno.xls',
                # 'vfdb': '/mnt/ilustre/users/sanger-dev/workspace/20170928/MetaGenomic_metagenome/output/vfdb/gene_vfdb_predict_anno.xls',
            }
            '''
            pass

        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'metagenomic')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        """
        检查参数设置
        """
        self.logger.info("check_sheet_data...")
        if not self._sheet.id:
            raise OptionError('需要输入task_id', code="12800101")
        if not self._sheet.member_type:
            raise OptionError("需要输入member_type", code="12800102")
        if not self._sheet.cmd_id:
            raise OptionError("需要输入cmd_id", code="12800103")
        self.logger.info("check options...")
        self.pipeline = 'pipeline' + self.option('pipeline')
        if self.option('old_task_id'):
            if self.option('old_task_db') == 0:
                db = get_old_mongo()[1]
            else:
                db = get_mongo()[1]
            old_task_info = db['sg_task'].find_one({'task_id': self.option('old_task_id')})
            if not old_task_info:
                self.set_error('输入的task_id为{}, {}中没有相关记录，请确保task_id正确且运行成功'.format(self.option('old_task_id'), db))
            if 'wf_version' in old_task_info:
                self.wf_version = old_task_info["wf_version"]
            else:
                self.wf_version = ''
        if self.pipeline == "pipeline3":
            self.option("qc", True)
        if self.pipeline == "pipeline1":
            if not self.option('in_fastq'):
                raise OptionError('需要输入原始fastq序列', code="12800104")
            if self.option('qc') and not self.option('speciman_info').is_set:
                raise OptionError('质控需提供样本集信息表', code="12800105")
            # if not self.option('qc') and not self.option('raw_info').is_set:
            #    raise OptionError('需进行质控，或者输入原始数据统计表', code="12800106")
            #if not self.option('qc') and not self.option('qc_info').is_set:
            #    raise OptionError('需进行质控，或者输入质控后数据统计表', code="12800107")
            if not self.option('qc_quality') > 0 and not self.option('qc_quality') < 42:
                raise OptionError('qc最小质量值超出范围，应在0~42之间', code="12800108")
            if not self.option('qc_length') > 0:
                raise OptionError('qc最小长度值超出范围，应大于0', code="12800109")
            if self.option('rm_host'):
                if self.option('ref_database') == '' and not self.option('ref_undefined').is_set:
                    raise OptionError('已选择去宿主，需输入参考数据库或参考序列', code="12800110")
                # if self.option('ref_database') not in ['', 'Custom'] and self.option('ref_undefined').is_set:
                #     raise OptionError('去宿主不可同时提供参考数据库及参考序列', code="12800111")
            if not self.option('assemble_type') in ['SOAPdenovo2', 'IDBA_UD', 'Megahit', 'Multiple', 'Multiple_Megahit',
                                                    'Multiple_IDBA_UD']:
                raise OptionError('拼接策略参数错误，应为SOAPdenovo2/IDBA_UD/Megahit/Multiple/Multiple_Megahit/Multiple_IDBA_UD', code="12800112")
            if self.option('min_contig') < 200 or self.option('min_contig') > 1000:
                raise OptionError('最小Contig长度参数超出范围200~1000', code="12800113")
            if self.option('min_gene') < 0:
                raise OptionError('最小基因长度参数不能为负', code="12800114")
            if not 0.75 <= self.option("cdhit_identity") <= 1:
                raise OptionError("cdhit identity必须在0.75，1之间", code="12800115")
            if not 0 <= self.option("cdhit_coverage") <= 1:
                raise OptionError("cdhit coverage必须在0,1之间", code="12800116")
        if not 0 < self.option("soap_identity") < 1:
            raise OptionError("soap identity必须在0，1之间", code="12800117")
        anno_list = ["kegg", "cog", "cazy", "ardb", "card", "vfdb"]
        new_list = []
        anno_list_e = self.option("anno_list").split(',')
        for i in range(len(anno_list_e)):
            if anno_list_e[i]:
                anno = anno_list[i]
                new_list.append(anno)
                self.option(anno + '_evalue', float(anno_list_e[i]))
        if not new_list:
            new_list = ["none"]
        self.option("anno_list", ','.join(new_list))
        if self.option("anno_list") == '':
            raise OptionError('anno_list 为空字符，不选比对数据库时应为none', code="12800118")
        elif self.option("anno_list") == 'none':
            self.choose_anno = []
            self.choose_nr = []
        else:
            self.choose_anno = self.option("anno_list").split(",")
            self.choose_nr = []
        if not self.option("anno_nr"):
            nr = "none"
        else:
            nr = "nr"
            nr_e = float(self.option("anno_nr"))
            self.option("nr_evalue", nr_e)
        self.option("anno_nr", nr)
        if self.option("anno_nr") == "nr":
            self.choose_anno.append("nr")
            self.choose_nr.append("nr_lca")
            self.choose_nr.append("nr_deunclass")
        elif self.option("anno_nr") != "none":
            raise OptionError('nr注释参数必须为nr或none', code="12800119")
        if not set(self.choose_anno).issubset(set(['kegg', 'cog', 'cazy', 'ardb', 'card', 'vfdb', 'nr'])):
            raise OptionError('数据库比对参数错误: %s, 请检查参数', variables=(self.choose_anno), code="12800120")
        if not self.option('group').is_set:
            if self.option('envtable').is_set:
                raise OptionError('设置环境因子请先设置group参数', code="12800121")
            if self.option('speciman_info').is_set:
                self.inner_group = self.get_group(self.option('speciman_info').prop['path'])
                self.logger.info("get inner_group !!!!!!!%s" % self.inner_group)
            else:
                if self.pipeline == "pipeline1":
                    self.inner_group = self.get_group(self.option('qc_info').prop['path'])
                elif self.pipeline == 'pipeline2.1':
                    self.inner_group = self.get_group_pipe2_1()
                else:
                    self.inner_group = self.get_group_pipe2(self.option('insertsize').prop['path'])
        if self.option('envtable').is_set:
            return_code = self.option('envtable').env_check()  # 对环境因子的数据类型进行检测 @20180301
        if self.pipeline == "pipeline1":
            check_speimen_names_messages = self.check_specimen_names()
            if check_speimen_names_messages != '':
                # check_speimen_names_messages = "以下文件样本名称不一致：" + check_speimen_names_messages
                raise OptionError("以下文件样本名称不一致：%s", variables=(check_speimen_names_messages), code="12800122")
        else:
            pass
        return True

    def check_specimen_names(self):
        list_set = self.read_file_to_set(self.option('in_fastq').prop['path'] + '/list.txt', header=False, list_index=1)
        error_message = ""
        if self.option('speciman_info').is_set:
            speciman_info_set = self.read_file_to_set(self.option('speciman_info').prop['path'], header=True)
            if list_set == speciman_info_set:
                self.logger.info("speciman_info sample check pass!!!")
            else:
                self.logger.error("\n speciman_info sample check has error:\n\t in_fastq samples: \t%s, \n\tspeciman_info samples: \t%s" % (list_set, speciman_info_set))
                error_message += "speciman_info, "
        if self.option('group').is_set:
            group_set = self.read_file_to_set(self.option('group').prop['path'], header=True)
            if list_set == group_set:
                self.logger.info("group sample check pass!!!")
            else:
                self.logger.error("\n group sample check has error:\n\t in_fastq samples: \t%s, \n\tsamples in group: \t%s" % (list_set, group_set))
                error_message += "group, "
        if self.option('envtable').is_set:
            env_set = self.read_file_to_set(self.option('envtable').prop['path'], header=True)
            if list_set == env_set:
                self.logger.info("envtable sample check pass!!!")
            else:
                self.logger.error("\n envtable sample check has error:\n\t in_fastq samples: \t%s, \n\tsamples in envtable: \t%s" % (list_set, env_set))
                error_message += "envtable, "
        if self.option('raw_info').is_set:
            raw_set = self.read_file_to_set(self.option('raw_info').prop['path'], header=True)
            if list_set == raw_set:
                self.logger.info("rawinfo sample check pass!!!")
            else:
                self.logger.error("\n raw_info sample check has error:\n\t in_fastq samples: \t%s, \n\tsamples in raw_info: \t%s" % (list_set, raw_set))
                error_message += "raw_info, "
        if self.option('qc_info').is_set:
            qc_set = self.read_file_to_set(self.option('qc_info').prop['path'], header=True)
            if list_set == qc_set:
                self.logger.info("qc_info sample check pass!!!")
            else:
                self.logger.error("\n qc_info sample check has error:\n\t in_fastq samples: \t%s, \n\tsamples in qc_info: \t%s" % (list_set, qc_set))
                error_message += "qc_info, "
        return error_message

    def read_file_to_set(self, file, header=True, list_index=0):
        specimen_set = set()
        try:
            list_index = int(list_index)
        except:
            raise OptionError("list_index is not integer!", code="12800123")
        with open(file, 'r') as rf:
            lines = rf.readlines()
            if header:
                lines.pop(0)
            for line in lines:
                t = line.split('\t')
                if len(t) <= list_index:
                    raise OptionError("%s 只有 %s 列, 不存在第 %s 列", variables=(file, len(t), list_index + 1), code="12800124")
                specimen_set.add(t[list_index])
        return specimen_set

    def get_group(self, path):
        inner_group_path = self.work_dir + '/inner_group'
        f = open(inner_group_path, 'w')
        f.write("#sample\tgroup_name\n")
        with open(path, 'r') as rf:
            lines = rf.readlines()
            for line in lines[1:]:
                t = line.split('\t')
                f.write(t[0] + '\t' + t[0] + '\n')
        f.close()
        return inner_group_path

    def get_group_pipe2(self, path):
        inner_group_path = self.work_dir + '/inner_group'
        f = open(inner_group_path, 'w')
        f.write("#sample\tgroup_name\n")
        with open(path, 'r') as rf:
            lines = rf.readlines()
            for line in lines[0:]:
                t = line.split('\t')
                f.write(t[0] + '\t' + t[0] + '\n')
        f.close()
        return inner_group_path

    def get_group_pipe2_1(self):
        # 导出老项目的分组文件
        inner_group_path = self.work_dir + '/inner_group'
        metag_api = Metagenomic()
        metag_api._config = self.config
        if self.option('old_task_db') == 0:
            db = get_old_mongo()[1]
        else:
            db = get_mongo()[1]
        i2n = id2name(self.option('old_task_id'), 'task', self.option('old_task_db'))
        sp_group = db['specimen_group'].find_one({'task_id': self.option('old_task_id')})
        self.group_name = None
        self.logger.info('sp_group: {}\n old_task_id: {}\n'.format(sp_group, self.option('old_task_id')))
        with open(inner_group_path, 'w') as f:
            f.write("#sample\tgroup_name\n")
            if sp_group:
                self.group_name = sp_group['group_name']
                category, specimen = sp_group['category_names'], sp_group['specimen_names']
                for i in range(len(category)):
                    for sp in specimen[i]:
                        f.write(i2n[sp] + '\t' + category[i] + '\n')
        return inner_group_path

    def get_sample(self):
        samp_list = []
        if self.option('group').is_set:
            path = self.option('group').prop['path']
        else:
            path = self.inner_group
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                t = line.split('\t')
                samp_list.append(t[0])
        return samp_list

    def get_env(self):
        path = self.option('envtable').prop['path']
        env_list = open(path, 'r').readline().rstrip('\n').split('\t')[1:]
        return env_list

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_run(self, opts, module, event, step, start=True):
        self.run_info(event)
        module.set_options(opts)
        module.on('start', self.set_step, {'start': step})
        module.on('end', self.set_step, {'end': step})
        module.on('end', self.set_output, event)
        if start:
            module.run()

    def run_info(self, run_type):
        if run_type == "sequence":
            self.wf_soft.append("fastp")
            self.soft_db_info["qc"] = {"soft": 18}
        elif run_type == "rm_host":
            self.wf_soft.append("bwa")
            self.wf_database.append("ensembl")
            self.soft_db_info["rm_host"] = {"soft": 19}
        elif run_type == "assem":
            if self.option('assemble_type') == 'SOAPdenovo2':
                self.wf_soft.append("SOAPdenovo2")
                soft = 26
            elif 'IDBA_UD' in self.option('assemble_type'):
                self.wf_soft.append("idba")
                soft = 34
            elif 'Megahit' in self.option('assemble_type'):
                self.wf_soft.append("megahit")
                soft = 20
            if 'Multiple' in self.option('assemble_type'):
                # self.wf_soft.append("newbler")
                self.wf_soft.append("bowtie")
                soft = str(soft) + ",30"
            self.soft_db_info["assem"] = {"soft": soft, "params": [str(self.option("min_contig"))]}
        elif run_type == "gene_predict":
            # prodigal metagene metagenemak
            p = {"MetaProdigal": "prodigal", "Prodigal": "prodigal", "MetaGene": "metagene", "MetaGeneMark": "metagenemark"}
            s = {"MetaProdigal": 22, "Prodigal": 22, "MetaGene": 23, "MetaGeneMark": 24}
            self.wf_soft.append(p[self.option("gene_predictor")])
            self.soft_db_info["gene_predict"] = {"soft": s[self.option("gene_predictor")], "params": [str(self.option("min_gene"))]}
        elif run_type == "gene_set":
            self.wf_soft.append("SOAPaligner")
            self.wf_soft.append("cdhit")
            self.soft_db_info["gene_set_n"] = {"soft": 25,
                                               "params": [str(self.option("cdhit_identity")),
                                                          str(self.option("cdhit_coverage"))]}
            self.soft_db_info["gene_set_a"] = {"soft": 26,
                                               "params": [str(self.option("soap_identity"))]}
        elif run_type == "nr":
            "diamond" in self.wf_soft or self.wf_soft.append("diamond")
            "nr" in self.wf_database or self.wf_database.append("nr")
            self.soft_db_info["nr"] = {"soft": 27, "db": 1,
                                         "params": ["blastp", "E-value ≤ " + str(self.option("nr_evalue"))]}
        elif run_type == "kegg":
            "diamond" in self.wf_soft or self.wf_soft.append("diamond")
            self.wf_database.append("kegg")
            self.soft_db_info["kegg"] = {"soft": 27, "db": 3,
                                           "params": ["blastp", "E-value ≤ " + str(self.option("kegg_evalue"))]}
        elif run_type == "cog":
            "diamond" in self.wf_soft or self.wf_soft.append("diamond")
            self.wf_database.append("cog")
            self.soft_db_info["cog"] = {"soft": 27, "db": 2,
                                          "params": ["blastp", "E-value ≤ " + str(self.option("cog_evalue"))]}
        elif run_type == "cazy":
            "hmmer" in self.wf_soft or self.wf_soft.append("hmmer")
            self.wf_database.append("cazy")
            self.soft_db_info["cazy"] = {"soft": 28, 'db': 36,
                                           "params": ["hmmscan", "E-value ≤ " + str(self.option("cazy_evalue"))]}
        elif run_type == "vfdb":
            "diamond" in self.wf_soft or self.wf_soft.append("diamond")
            self.wf_database.append("vfdb")
            self.soft_db_info["vfdb"] = {"soft": 27, 'db': 5,
                                           "params": ["blastp", "E-value ≤ " + str(self.option("vfdb_evalue"))]}
        elif run_type == "ardb":
            "diamond" in self.wf_soft or self.wf_soft.append("diamond")
            self.wf_database.append("ardb")
            self.soft_db_info["ardb"] = {"soft": 27, 'db': 4,
                                           "params": ["blastp", "E-value ≤ " + str(self.option("ardb_evalue"))]}
        elif run_type == "card":
            "diamond" in self.wf_soft or self.wf_soft.append("diamond")
            self.wf_database.append("card")
            self.soft_db_info["card"] = {"soft": 27, 'db': 37,
                                           "params": ["blastp", "E-value ≤ " + str(self.option("card_evalue"))]}
        elif run_type == "kraken":
            self.wf_soft.append("kraken2")
            self.wf_database.append("kraken2_db")
            self.soft_db_info["kraken2"] = {"soft": 32, 'db': 16,
                                             "params": ["confidence: " + str(self.option("kk_confidence"))]}
        elif run_type == "metaphlan":
            self.wf_soft.append("metaphlan3")
            self.wf_database.append("metaphlan3_db")
            self.soft_db_info["metaphlan3"] = {"soft": 31, 'db': 15,
                                             "params": ["min_cu_len : " + str(self.option("mph_min_cu_len ")),
                                                        "stat: " + self.option("mph_stat")]}
            if self.option("mph_stat") in ["tavg_g", "tavg_l", "wavg_g", "wavg_l"]:
                self.soft_db_info["metaphlan3"]["params"].append("stat_q: "+ str(self.option("mph_stat_q")))

    def run_sequence(self):
        opts = {
            'input_dir': self.option("in_fastq").path,
            'qc_quality': self.option('qc_quality'),
            'qc_length': self.option('qc_length'),
            'qc': self.option('qc'),
            'sample_info': self.option('speciman_info')
        }
        self.set_run(opts, self.sequence, 'sequence', self.step.sequence)

    def run_qc(self):
        opts = {
            'fastq_dir': self.sequence.output_dir + '/data',
            'insert_size': self.option('speciman_info').prop['path'],
            'stat_dir': self.sequence.output_dir + '/base_info',
            'qc_quality': self.option('qc_quality'),
            'qc_length': self.option('qc_length'),
            'fastp': self.option("fastp")
        }
        self.set_run(opts, self.qc, 'qc', self.step.qc_)

    def get_ref_path(self, ensembl_path):
        file_path = os.path.join(ensembl_path, 'ref_path.txt')
        ref_path = {}
        with open(file_path, 'r') as r:
            for l in r:
                line = l.strip().split('\t')
                key = "{}:{}".format(line[0], line[1])
                ref_path[key] = line[2]
        return ref_path

    def run_rm_host(self):
        ref_database = []
        ensembl_path = self.config.SOFTWARE_DIR + '/database/ensembl/'
        ref_fna_path = self.get_ref_path(ensembl_path)
        for one in self.option("ref_database").split(','):
            if one != "Custom":
                ref_database.append(ref_fna_path[one])
        opts = {"fq_type": "PE", "fastq_dir": self.sequence.fq_path}
        ref_path = None
        if not os.path.exists(self.work_dir + "/host_fasta"):
            os.mkdir(self.work_dir + "/host_fasta")
        if len(ref_database) == 1:
            ref_path = os.path.join(ensembl_path, ref_database[0])
        elif len(ref_database) > 1:
            for one in ref_database:
                self.link(os.path.join(ensembl_path, one), "host_fasta/")
            ref_path = self.work_dir + "/host_fasta"
        if self.option("ref_undefined").is_set:
            if ref_path:
                for f in os.listdir(self.option("ref_undefined").path):
                    self.link(os.path.join(self.option("ref_undefined").path, f), "host_fasta/")
                if "host_fasta" not in ref_path:
                    self.link(os.path.join(ensembl_path, ref_database[0]), "host_fasta/")
                ref_path = self.work_dir + "/host_fasta"
            else:
                ref_path = self.option("ref_undefined")
        if len(ref_database) == 1 and "Custom" not in self.option("ref_database"):
            opts["ref_database"] = ref_path
        else:
            opts["ref_database"] = ''
            opts["ref_undefined"] = ref_path
        self.set_run(opts, self.rm_host, 'rm_host', self.step.rm_host)

    def run_assem(self):
        if self.option('qc'):
            opts = {
                'qc_stat': self.sequence.qc_stat,
                'raw_stat': self.sequence.raw_stat,
                'QC_dir': self.sequence.fq_path
            }
        else:
            opts = {
                'qc_stat': self.option('qc_info'),
                'raw_stat': self.option('raw_info'),
                'QC_dir': self.option('in_fastq')
            }
            if not self.option('raw_info').is_set:
                opts["raw_stat"] = self.sequence.raw_stat
            if not self.option('qc_info').is_set:
                opts["qc_stat"] = self.sequence.qc_stat
        if self.option('rm_host'):
            opts['QC_dir'] = self.rm_host.option('result_fq_dir')

        opts['min_contig'] = self.option('min_contig')
        if self.option('assemble_type') == 'SOAPdenovo2':
            self.set_run(opts, self.assem_soapdenovo, 'assem', self.step.assem)
        else:
            if 'IDBA_UD' in self.option('assemble_type'):
                opts['assemble_tool'] = 'idba'
            elif 'Megahit' in self.option('assemble_type'):
                opts['assemble_tool'] = 'megahit'
            if 'Multiple' in self.option('assemble_type'):
                opts['method'] = 'multiple'
                opts['use_newbler'] = self.option('use_newbler')  # 加入参数use_newbler @ 20181023
            else:
                opts['method'] = 'simple'
                opts['use_newbler'] = 'False'
            self.set_run(opts, self.assem_idba, 'assem', self.step.assem)

    def run_gene_predict(self):
        opts = {
            'min_gene': str(self.option('min_gene')),
        }
        if self.option('assemble_type') == 'SOAPdenovo2':
            opts['input_fasta'] = self.assem_soapdenovo.option('contig').prop['path']
        else:
            opts['input_fasta'] = self.assem_idba.option('contig').prop['path']
        if self.option("gene_predictor") == "metagene":
            #self.gene_predict = self.add_module('gene_structure.gene_predict')
            pass
        else:
            opts["predict_meth"] = self.option("gene_predictor")
        self.set_run(opts, self.gene_predict, 'gene_predict', self.step.gene_predict)

    def run_gene_set(self):
        if self.pipeline == "pipeline1":
            opts = {
                'insertsize': self.option('insertsize'),
                'cdhit_identity': self.option('cdhit_identity'),
                'cdhit_coverage': self.option('cdhit_coverage'),
                'soap_identity': self.option('soap_identity'),
                'ana_type': 'prot',
            }
            if self.option('insertsize').is_set:
                opts['insertsize'] = self.option('insertsize')
            elif self.option('raw_info').is_set:
                opts['insertsize'] = get_insertsize(self.option('raw_info').prop['path'], self.gene_set.work_dir)
            else:
                opts['insertsize'] = get_insertsize(self.option("speciman_info").path, self.gene_set.work_dir)
            if self.option('rm_host'):
                opts['QC_dir'] = self.rm_host.option('result_fq_dir')
            elif self.option('qc'):
                opts['QC_dir'] = self.sequence.fq_path
            else:
                opts['QC_dir'] = self.option('in_fastq')
            # 混拼，单拼序列合起来做基因集
            # if "Multiple" in self.option('assemble_type'):
            #     opts['gene_tmp_fa'] = self.gene_predict.option('out_fa')
            #     opts['gene_tmp_faa'] = self.gene_predict.option('out_faa')
            #     opts['gene_tmp_fa_mix'] = self.gene_predict.option('out_fa_mix')
            #     opts['gene_tmp_faa_mix'] = self.gene_predict.option('out_faa_mix')
            # else:
            #     opts['gene_tmp_faa'] = self.gene_predict.option('out')
            #     opts['gene_tmp_fa'] = self.gene_predict.option('out_fa')
            opts['gene_tmp_faa'] = self.gene_predict.option('out')
            opts['gene_tmp_fa'] = self.gene_predict.option('out').path[:-1]
        elif self.pipeline == "pipeline2":
            opts = {
                'insertsize': self.option('insertsize'),
                'QC_dir': os.path.join(self.sequence_deal.output_dir,"clean_dir"),
                'uni_fasta': self.option('unique_fa'),
                'soap_identity': self.option('soap_identity'),
                "map_type": 2
            }
        self.set_run(opts, self.gene_set, 'gene_set', self.step.gene_set)
        self.anno_table['geneset'] = os.path.join(self.gene_set.output_dir,
                                                  'gene_profile/RPKM.xls')

    def run_nr(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'query_type': "prot",
            'database': 'nr_v20200604',
            'lines': '50000',
            "target_num": 5,
            "evalue": self.option("nr_evalue"),
        }
        self.set_run(opts, self.nr, 'nr', self.step.nr_)

    def run_kegg(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'query_type': "prot",
            'database': 'kegg_v94.2',
            'lines': '50000',
            "evalue": self.option("kegg_evalue"),
        }
        self.set_run(opts, self.kegg, 'kegg', self.step.kegg)

    def run_cog(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'query_type': "prot",
            'database': 'eggnog',
            'lines': '50000',
            "evalue": self.option("cog_evalue"),
        }
        self.set_run(opts, self.cog, 'cog', self.step.cog)

    def run_anno(self):
        opts = {
            'reads_profile_table': self.gene_set.option('reads_abundance'),
        }
        if 'nr' in self.choose_anno:
            opts['nr_xml_dir'] = self.nr.option('outxml_dir')
            opts['out_type'] = 1  ##zouguanqing 20190314
        if 'kegg' in self.choose_anno:
            opts['kegg_xml_dir'] = self.kegg.option('outxml_dir')
        if 'cog' in self.choose_anno:
            opts['cog_xml_dir'] = self.cog.option('outxml_dir')
        self.set_run(opts, self.anno, 'anno', self.step.anno, False)
        if 'nr' in self.choose_anno:
            self.nr_dir = os.path.join(self.anno.output_dir, 'nr_tax_level')
            self.anno_table['nr'] = os.path.join(self.nr_dir, 'gene_nr_anno.xls')
            self.anno_table['nr_blastout'] = os.path.join(self.nr_dir, 'nr_align_table.xls')
        if 'cog' in self.choose_anno:
            self.cog_dir = os.path.join(self.anno.output_dir, 'cog_result_dir')
            self.anno_table['cog'] = os.path.join(self.cog_dir, 'gene_cog_anno.xls')
        if 'kegg' in self.choose_anno:
            self.kegg_dir = os.path.join(self.anno.output_dir, 'kegg_result_dir')
            self.anno_table['kegg'] = os.path.join(self.kegg_dir, 'gene_kegg_anno.xls')
        self.anno.run()

    def run_cazy(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            "evalue": self.option("cazy_evalue"),
        }
        self.set_run(opts, self.cazy, 'cazy', self.step.cazy, False)
        self.anno_table['cazy'] = os.path.join(self.cazy.output_dir, 'anno_result', 'gene_cazy_anno.xls')
        self.cazy.run()

    def run_vfdb(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            'lines': '400000',
            "evalue": self.option("vfdb_evalue"),
        }
        self.set_run(opts, self.vfdb, 'vfdb', self.step.vfdb, False)
        self.anno_table['vfdb'] = os.path.join(self.vfdb.output_dir, 'gene_vfdb_total_anno.xls')
        self.vfdb.run()

    def run_ardb(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            'lines': '400000',
            "evalue": self.option("ardb_evalue"),
        }
        self.set_run(opts, self.ardb, 'ardb', self.step.ardb, False)
        self.anno_table['ardb'] = os.path.join(self.ardb.output_dir, 'gene_ardb_anno.xls')
        self.ardb.run()

    def run_card(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            'lines': '400000',
            "evalue": self.option("card_evalue"),
        }
        self.set_run(opts, self.card, 'card', self.step.card, False)
        self.anno_table['card'] = os.path.join(self.card.output_dir, 'gene_card_anno.xls')
        self.card.run()

    def run_overview(self):
        opts = {
            'gene_length_table': self.gene_set.option('gene_length_list')
        }
        for anno in self.choose_anno:
            opt_k = 'gene_' + anno + '_anno'
            opts[opt_k] = self.anno_table[anno]
        self.set_run(opts, self.overview, 'overview', self.step.overview)

    def run_analysis(self, event):
        if type(event) is not str:
            event = event['data']

        for db in self.choose_anno:
            ana_str = '1'
            if event == 'all' and self.default_level2[db] == self.default_level1[db]:
                ana_str = '1,2'
            elif event == 'all' and self.default_level2[db] != self.default_level1[db]:
                # ana_str = '1'
                count = 0
                for index, line in enumerate(open(self.anno_table[db], 'r')):
                    count += 1
                self.logger.info("card 注释的结果含行数：{}".format(count))
                if count > 1:
                    self.run_new_table(self.anno_table[db], self.anno_table['geneset'], self.default_level2[db], db=db,ana_str='2')
            count = 0
            for index, line in enumerate(open(self.anno_table[db], 'r')):
                count += 1
            self.logger.info("{}注释的结果含行数:{}".format(db, count))
            if count > 1:
                self.run_new_table(self.anno_table[db], self.anno_table['geneset'], self.default_level1[db], db=db,ana_str=ana_str)
        if len(self.new_table) != 1:
            self.on_rely(self.new_table, self.run_analysis3)
            for module in self.new_table:
                module.run()
        else:
            self.new_table[0].on('end', self.run_analysis3)
            self.new_table[0].run()

    def run_analysis3(self):
        composition_path = os.path.join(self.work_dir, 'new_abund_file')
        other_path = os.path.join(self.work_dir, 'new_abund_file2')
        check_path = self.check_abund_file(composition_path)
        for db in check_path:
            profile_table1 = os.path.join(self.work_dir, 'new_abund_file', db, 'new_abund_table.xls')
            if self.option('group').is_set:
                group_table = self.option('group').prop['path']  # add_batch调用的参数必须可实例化 20180323
            else:
                group_table = self.inner_group
            self.logger.info(">>>report group:")
            self.logger.info(group_table)
            if db in ['nr', 'gene', 'ardb', 'card']:
                self.func_composition(profile_table1, db, group_table)
            elif db in ['cog', 'kegg', 'cazy', 'vfdb']:
                self.func_composition(profile_table1, db, group_table, analysis='heatmap,circos')
                # self.composition_dir2anno[self.composition.output_dir] = db  # batch 对象没有output_dir
                # self.compare_dir2anno[self.composition.work_dir] = db
                self.func_composition(profile_table1, db, group_table, analysis='bar', others=0)
            # self.composition_dir2anno[self.composition.output_dir] = db  # batch 对象没有output_dir
            ##self.compare_dir2anno[self.composition.work_dir] = db
        if os.path.exists(other_path):
            check_path = self.check_abund_file(other_path)
            for db in check_path:
                profile_table2 = os.path.join(self.work_dir, 'new_abund_file2', db, 'new_abund_table.xls')
                if self.option('group').is_set:
                    self.logger.info("report compare group:")
                    self.logger.info(self.option('group').prop['path'])
                    self.func_compare(profile_table2, db, self.option('group').prop['path'])  # add_batch调用的参数必须可实例化 20180323
                else:
                    self.func_compare(profile_table2, db)
                # self.compare_dir2anno[self.compare.output_dir] = db  # batch 对象没有output_dir
                ##self.compare_dir2anno[self.compare.work_dir] = db
                if self.option('envtable').is_set:
                    env_list = self.get_env()
                    if len(env_list) >= 2:
                        self.logger.info("report cor_heatmap env:")
                        self.logger.info(self.option('envtable').prop['path'])
                        self.func_correlation(profile_table2, db, self.option('envtable').prop['path'])  # add_batch调用的参数必须可实例化 20180323
                        # self.correlation_dir2anno[self.correlation.output_dir] = db
                        ##self.correlation_dir2anno[self.correlation.work_dir] = db  # batch 对象没有output_dir
                        # self.logger.info("report %s database correlation workdir: %s" % (db,self.correlation.work_dir))
                        # self.anno2correlation_tree[db] = self.correlation.work_dir
        if len(self.analysis) != 1:
            if self.option("personal_anno") != "none":
                self.on_rely(self.analysis, self.version2_add)
            else:
                self.on_rely(self.analysis, self.end)
            gevent_list = []
            for module in self.analysis:
                gevent_list.append(gevent.spawn(module.run))
            gevent.joinall(gevent_list)
        else:
            if self.option("personal_anno") != "none":
                self.composition.on("end", self.version2_add)
            else:
                self.composition.on("end", self.end)
            self.composition.run()

    ###### version 新增函数
    def run_personal_anno(self):
        gene_anno = self.anno_table['nr']
        gene_anno_lca = self.anno_table['nr_lca']
        gene_anno_deunclass = self.anno_table['nr_deunclass']
        if self.option('group').is_set:
            group_table = self.option('group').prop['path']
        else:
            group_table = self.inner_group
        if self.option("personal_anno") != "none":
            database = self.option("personal_anno")
            self.personal_anno = self.add_module('annotation.personal_anno')
            self.logger.info("test>>>>>>>>>>>>>>>>>>>>>>>>>")
            self.logger.info(database)
            #database = "qs;probio;p450;pfam;mvirdb;phi;tcdb;go"
            opts = {
                'database': database,
                'query': self.gene_set.option('uni_fastaa'),
                'reads_profile_table': self.gene_set.option('reads_abundance'),
                'blastout': self.anno_table['nr_blastout'],
                'nr_gene_anno': gene_anno,
                'nr_gene_anno_de': gene_anno_deunclass,
                'nr_gene_anno_lca': gene_anno_lca,
                'group_table': group_table,
            }
            self.set_run(opts, self.personal_anno, 'anno_personal', self.step.anno_personal)
        else:
            self.end()

    def version2_add(self):
        self.personal_anno.on("end", self.end)
        self.run_personal_anno()

    def run_anno_new(self):
        self.run_anno()
        self.run_anno_nr_lca()
        self.run_anno_nr_deunclass()

    def run_anno_nr_lca(self):
        opts = {
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            "nr_method": "lca",
            "nr_xml_dir": self.nr.option('outxml_dir'),
            'out_type': 1   ##zouguanqing
        }
        self.set_run(opts, self.anno_nr_lca, 'anno_lca', self.step.anno, False)
        #self.anno_nr_lca.set_options(opts)
        self.nr_lca_dir = os.path.join(self.anno_nr_lca.output_dir, 'nr_tax_level')
        self.anno_table['nr_lca'] = os.path.join(self.nr_lca_dir, 'gene_nr_anno.xls')
        self.anno_nr_lca.run()

    def run_anno_nr_deunclass(self):
        opts = {
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            "nr_method": "deunclassied",
            "nr_xml_dir": self.nr.option('outxml_dir'),
            'out_type' : 1   ##zouguanqing
        }
        self.set_run(opts, self.anno_nr_deunclass, 'anno_deunclass', self.step.anno, False)
        #self.anno_nr_deunclass.set_options(opts)
        self.nr_deunclass_dir = os.path.join(self.anno_nr_deunclass.output_dir, 'nr_tax_level')
        self.anno_table['nr_deunclass'] = os.path.join(self.nr_deunclass_dir, 'gene_nr_anno.xls')
        self.anno_nr_deunclass.run()

    def run_sequence_deal(self):
        opts = {
            'fastq_dir': self.option('in_fastq'),
            'clean_list': self.option('clean_list'),
            'insertsize': self.option('insertsize')
        }
        self.set_run(opts, self.sequence_deal, 'sequence_deal', self.step.sequence)

    def run_kraken(self):
        self.kraken = self.add_module("metagenome.kraken2")
        self.kraken.on("end", self.run_kraken_outset)
        self.step.add_steps("kraken")
        opts = {
            "confidence": self.option("kk_confidence")
        }
        if self.option("qc"):
            opts['fq_dir'] = self.sequence.fq_path
        else:
            opts['fq_dir'] = self.option('in_fastq')
        self.set_run(opts, self.kraken, 'kraken', self.step.kraken)

    def run_kraken_outset(self):
        kraken_outset = self.add_tool("metagenomic.taxon_outset")
        self.taxon_outset["kraken"] = kraken_outset
        opts = {"result_dir": self.kraken.output_dir,
                "name2id": json.dumps(self.spname_spid)}
        kraken_outset.set_options(opts)
        kraken_outset.run()

    def run_metaphlan(self):
        self.metaphlan = self.add_module("metagenome.metaphlan3")
        self.metaphlan.on("end", self.run_metaphlan_outset)
        self.step.add_steps("metaphlan")
        opts = {
            "min_cu_len": self.option("mph_min_cu_len"),
            "stat": self.option("mph_stat"),
            "stat_q": self.option("mph_stat_q"),
        }
        if self.option("qc"):
            opts['fq_dir'] = self.sequence.fq_path
        else:
            opts['fq_dir'] = self.option('in_fastq')
        self.set_run(opts, self.metaphlan, 'metaphlan', self.step.metaphlan)

    def run_metaphlan_outset(self):
        metaphlan_outset = self.add_tool("metagenomic.taxon_outset")
        self.taxon_outset["metaphlan"] = metaphlan_outset
        opts = {"result_dir": self.metaphlan.output_dir,
                "name2id": json.dumps(self.spname_spid)}
        metaphlan_outset.set_options(opts)
        metaphlan_outset.run()

    def check_abund_file(self, path):
        """
        此函数是为了检查计算的丰度表，如果该数据库的丰度表没有获得到，一直等待直到获得结果；
        如果注释结果为空，那么在进行组成分析时应该跳过；
        :param path:
        :return:
        """
        check_path = os.listdir(path)
        interrupt = 'wait'
        self.logger.info("choose databse:%s"%(self.choose_anno))
        while interrupt in ["wait"]:
            interrupt = 'run'
            for db in self.choose_anno:  # 注： 不含gene
                if db not in check_path:
                    self.logger.info("db: %s" %db)
                    count = 0
                    for index, line in enumerate(open(self.anno_table[db], 'r')):
                        count += 1
                    if count > 1:
                        interrupt = 'wait'
                        time.sleep(10)
                    else: ##说明结果为空
                        interrupt = 'run'
                        time.sleep(10)
                    break
                else:
                    file = os.path.join(path, db, 'new_abund_table.xls')
                    self.logger.info("file_path:%s"%file)
                    if not os.path.isfile(file):
                        interrupt = 'wait'
                        time.sleep(10)
                        break
        return check_path

    def wait_file(self, path, wait_times=1):
        '''
        增加等待文件结果方法
        :param path: 结果文件路径
        :param wait_times: 等待次数
        :return: 文件路径
        :time: 20180425
        '''
        while wait_times < 50:
            self.logger.info("等待第%s次" % wait_times)
            if not os.path.isfile(path):
                time.sleep(10)
                wait_times += 1
                self.wait_file(path, wait_times=wait_times)
            return path
        self.logger.info("超过文件等待次数，需检查文件%s" % path)
        return

    def run_new_table(self, anno, gene, level, db, ana_str='1,2'):
        if level != "":
            opts = {
                'anno_table': anno,
                'geneset_table': gene,
                'level_type': level,
            }
        else:
            opts = {
                'gene_list': anno,
                'geneset_table': gene,
            }
        self.table[self.table_pointer] = self.add_tool('meta.create_abund_table')
        self.set_run(opts, self.table[self.table_pointer], 'table', self.step.table, False)
        ana_list = ana_str.split(',')
        for ana in ana_list:
            db_ana_pointer = db + ',' + ana + ',' + str(self.table_pointer)
            self.table[self.table_pointer].on('end', self.set_abund_file, db_ana_pointer)
        self.new_table.append(self.table[self.table_pointer])
        self.table_pointer += 1

    def set_abund_file(self, event):
        db_ana_pointer = event['data']
        db = db_ana_pointer.split(',')[0]
        ana = db_ana_pointer.split(',')[1]
        pointer = int(db_ana_pointer.split(',')[2])
        abund_dir = ''
        if ana == '1':
            abund_dir = os.path.join(self.work_dir, 'new_abund_file')
        elif ana == '2':
            abund_dir = os.path.join(self.work_dir, 'new_abund_file2')
        if not os.path.isdir(abund_dir):
            os.mkdir(abund_dir)
        abund_file_path = os.path.join(abund_dir, db)
        if not os.path.isdir(abund_file_path):
            os.mkdir(abund_file_path)
        old_table_file = os.path.join(self.table[pointer].output_dir, 'new_abund_table.xls')
        target_path = os.path.join(abund_file_path, 'new_abund_table.xls')
        if os.path.exists(target_path):
            os.remove(target_path)
        os.link(old_table_file, target_path)

    def func_composition(self, abund, anno, group, analysis='bar,heatmap,circos', others=0.01):
        opts = {
            'analysis': analysis,
            'group': group,
            'abundtable': abund,
            'species_number': '50',
            'others': others,
        }
        # self.composition = self.add_module('meta.composition.composition_analysis')
        self.composition = self.add_batch('meta.composition.composition_analysis', ignore_error=True, batch_type="module")  # 出错也没关系 @ 20180321
        self.set_run(opts, self.composition, 'composition__' + anno + '__' + analysis, self.step.composition, False)
        self.analysis.append(self.composition)

    def func_compare(self, abund, anno, group=''):
        opts = {
            'dis_method': 'bray_curtis',
            'otutable': abund,
        }
        if group != '':
            self.logger.info("group exists!!!! ->>|%s|<<-" % group)
            opts['group'] = group
        if self.option('envtable').is_set:
            self.logger.info("envtalbe exists!!! ->>|%s|<<-" % self.option('envtable').prop['path'])
            env_list = self.get_env()
            if len(env_list) <= len(self.sample_in_group):
                opts['envtable'] = self.option('envtable').prop['path']
                opts['analysis'] = 'distance,pca,pcoa,nmds,rda_cca,dbrda,hcluster'
            else:
                opts['analysis'] = 'distance,pca,pcoa,nmds,hcluster'
        else:
            opts['analysis'] = 'distance,pca,pcoa,nmds,hcluster'
        # self.compare = self.add_module('meta.beta_diversity.beta_diversity')
        self.compare = self.add_batch('meta.beta_diversity.beta_diversity', ignore_error=True, batch_type="module")  # 出错也没关系 @ 20180321
        self.set_run(opts, self.compare, 'compare__' + anno + '__' + opts['analysis'], self.step.compare, False)
        self.analysis.append(self.compare)

    def func_correlation(self, abund, anno, envtable):
        opts = {
            'method': 'spearmanr',
            'otutable': abund,
            'envtable': envtable,
            "top_species": 50,
        }
        # self.correlation = self.add_tool('statistical.pearsons_correlation')
        self.correlation = self.add_batch('statistical.pearsons_correlation', ignore_error=True, batch_type="tool")  # 出错也没关系 @ 20180321
        self.set_run(opts, self.correlation, 'correlation__' + anno + '__cor_heatmap', self.step.correlation, False)
        self.analysis.append(self.correlation)

    '''处理输出文件'''



    def set_output(self, event):
        """
        将各个模块的结果输出至output
        """
        obj = event['bind_object']
        if event['data'] == 'sequence':
            self.move_dir(obj.rawdata, 'rawdata')
            if self.option('qc'):
                self.move_dir(obj.fq_path, 'qc/after_qc_dir')
                if os.path.exists(self.output_dir + "/qc/reads.cleanData.stat.xls"):
                    os.remove(self.output_dir + "/qc/reads.cleanData.stat.xls")
                self.move_file(obj.qc_stat, os.path.join(self.output_dir, "qc/reads.cleanData.stat.xls"))
        elif event['data'] == "sequence_deal":
            self.move_dir(obj.output_dir, 'optimize_data')
        elif event['data'] == 'qc':
            self.move_dir(os.path.join(obj.output_dir, 'after_qc_dir'), 'qc/after_qc_dir')
            self.move_dir(os.path.join(obj.output_dir, 'qc_stat'), 'qc/qc_stat')
        elif event['data'] == 'rm_host':
            self.move_dir(obj.output_dir, 'rm_host')
        elif event['data'] == 'assem':
            self.move_dir(obj.output_dir, 'assemble')
        elif event['data'] == 'gene_predict':
            self.move_dir(obj.output_dir, 'predict')
        elif event['data'] == 'gene_set':
            self.move_dir(obj.output_dir, 'geneset')
        elif event['data'] == 'anno':
            # if self.option('nr'):
            if 'nr' in self.choose_anno:
                self.move_dir(self.nr_dir, 'nr')
            # if self.option('cog'):
            if 'cog' in self.choose_anno:
                self.move_dir(self.cog_dir, 'cog')
            if 'kegg' in self.choose_anno:
                self.move_dir(self.kegg_dir, 'kegg')
        elif event['data'] == 'anno_lca':
            self.move_dir(self.nr_lca_dir, 'nr_lca')
        elif event['data'] == 'anno_deunclass':
            self.move_dir(self.nr_deunclass_dir, 'nr_deunclassied')
        elif event['data'] == 'anno_personal':
            self.move_dir(obj.output_dir, 'anno_personal')
        elif event['data'] == 'cazy':
            self.move_dir(os.path.join(obj.output_dir, 'anno_result'), 'cazy')
        elif event['data'] == 'vfdb':
            self.move_dir(obj.output_dir, 'vfdb')
        elif event['data'] == 'ardb':
            self.move_dir(obj.output_dir, 'ardb')
        elif event['data'] == 'card':
            self.move_dir(obj.output_dir, 'card')
        elif event["data"] == "kraken":
            if not os.path.exists(os.path.join(self.output_dir, "kraken")):
                os.mkdir(os.path.join(self.output_dir, "kraken"))
            for f in os.listdir(obj.output_dir):
                if f.endswith(".taxon.xls") or f.endswith(".taxon_bracken.xls"):
                    file_path = os.path.join(obj.output_dir, f)
                    self.link(file_path, "output/kraken/")
        elif event["data"] == "metaphlan":
            if not os.path.exists(os.path.join(self.output_dir, "metaphlan")):
                os.mkdir(os.path.join(self.output_dir, "metaphlan"))
            for f in os.listdir(obj.output_dir):
                if f.endswith(".taxon.xls"):
                    file_path = os.path.join(obj.output_dir, f)
                    self.link(file_path, "output/metaphlan/")
        elif event['data'] == 'overview':
            self.logger.info("DEBUG: set overview output")
            old_file = self.overview.option('gene_overview').prop['path']
            geneset_path = os.path.join(self.output_dir, 'geneset', 'anno_overview.xls')
            if os.path.isfile(geneset_path):
                os.remove(geneset_path)
            os.link(old_file, geneset_path)
        elif len(event['data'].split('__')) == 3:
            analysis, anno, allfiles = event['data'].split('__')  # 需检查具体分析的大小写,双下划线
            dir_dic = {  # 上传文件名:module结果文件夹名
                'bar': 'bar',
                'heatmap': 'heatmap',
                'circos': 'circos',
                'distance': 'Distance',
                'pca': 'Pca',
                'pcoa': 'Pcoa',
                'nmds': 'Nmds',
                'hcluster': 'Hcluster',
                'rda_cca': 'Rda',
                'dbrda': 'Dbrda',
                'cor_heatmap': '',
            }
            analysis_dic = {  # 新路径名:module名
                'composition': 'CompositionAnalysis',
                'compare': 'BetaDiversity',
                'correlation': 'PearsonsCorrelation'
            }
            for dir in allfiles.split(','):
                dir_path = os.path.join(obj.work_dir, analysis_dic[analysis], "output", dir_dic[dir])
                if os.path.isdir(dir_path):
                    if dir in ['bar', 'heatmap', 'circos']:
                        self.move_dir(dir_path, os.path.join(analysis, dir_dic[dir], anno))  # 组成分析
                    elif dir in ['distance','pca','pcoa','nmds','hcluster']:  # 比较分析
                        self.move_dir(dir_path, os.path.join(analysis, dir_dic[dir], anno))
                    elif dir in ['rda_cca','dbrda']:
                        self.move_dir(dir_path, os.path.join('correlation', dir_dic[dir], anno))  # 环境因子相关性分析
                    else:
                        new_dir = os.path.join('correlation', dir, anno)  # 相关性heatmap图
                        self.move_dir(dir_path, new_dir)  # 相关性heatmap图
                        os.rename(os.path.join(self.output_dir, new_dir, "pearsons_correlation.xls"),
                                  os.path.join(self.output_dir, new_dir, "correlation.xls"))
                        os.rename(os.path.join(self.output_dir, new_dir, "pearsons_pvalue.xls"),
                                  os.path.join(self.output_dir, new_dir, "pvalue.xls"))
                    self.logger.info("在%s水平下的%s分析，移动输出路径。from：%s" % (anno, dir, dir_path))
                else:
                    self.logger.info("在%s水平下的%s分析未能正常结束，跳过output输出\n对应路径：%s" % (anno, dir, dir_path))
        else:
            self.logger.error("set_output传递的参数不正确：%s" % event['data'])
        """
        if event['data'] == 'composition':
            # anno = self.composition_dir2anno[obj.output_dir]  # batch对象没有output_dir
            # anno = self.composition_dir2anno[obj.work_dir]
            # allfiles = os.listdir(obj.output_dir)  # batch对象没有output_dir
            allfiles = os.listdir(os.path.join(obj.work_dir, "output"))
            if not allfiles:  # 防止此分析未完成而进行结果输出 @ 20180321
                self.logger.info("%s的composition无结果，跳过此处的output输出" % anno)
                return
            for dir in allfiles:
                # self.move_dir(os.path.join(obj.output_dir, dir), os.path.join('composition', dir, anno))  # batch对象没有output_dir
                self.move_dir(os.path.join(obj.work_dir, "output", dir), os.path.join('composition', dir, anno))
        if event['data'] == 'compare':
            # anno = self.compare_dir2anno[obj.output_dir]  # batch对象没有output_dir
            # anno = self.compare_dir2anno[obj.work_dir]
            # allfiles = os.listdir(obj.output_dir)  # batch对象没有output_dir
            allfiles = os.listdir(os.path.join(obj.work_dir, "output"))
            if not allfiles:  # 防止此分析未完成而进行结果输出 @ 20180321
                self.logger.info("%s的compare无结果，跳过此处的output输出" % anno)
                return
            for dir in allfiles:
                if dir in ['Pca', 'Pcoa', 'Hcluster', 'Nmds', 'Distance']:
                    # self.move_dir(os.path.join(obj.output_dir, dir), os.path.join('compare', dir, anno))  # batch对象没有output_dir
                    self.move_dir(os.path.join(obj.work_dir, "output", dir), os.path.join('compare', dir, anno))
                else:
                    # self.move_dir(os.path.join(obj.output_dir, dir), os.path.join('correlation', dir, anno))  # batch对象没有output_dir
                    self.move_dir(os.path.join(obj.work_dir, "output", dir), os.path.join('correlation', dir, anno))
        if event['data'] == 'correlation':
            # anno = self.correlation_dir2anno[obj.output_dir]  # batch对象没有output_dir
            # anno = self.correlation_dir2anno[obj.work_dir]
            # if not os.listdir(obj.output_dir):  # 防止此分析未完成而进行结果输出 @ 20180321
            if not os.listdir(os.path.join(obj.work_dir, "output")):  # batch对象没有output_dir
                self.logger.info("%s的correlation无结果，跳过此处的output输出" % anno)
                return
            new_dir = os.path.join('correlation', 'cor_heatmap', anno)
            # self.move_dir(obj.output_dir, new_dir)  # batch对象没有output_dir
            self.move_dir(os.path.join(obj.work_dir, "output"), new_dir)
            os.rename(os.path.join(self.output_dir, new_dir, "pearsons_correlation.xls"),
                      os.path.join(self.output_dir, new_dir, "correlation.xls"))
            os.rename(os.path.join(self.output_dir, new_dir, "pearsons_pvalue.xls"),
                      os.path.join(self.output_dir, new_dir, "pvalue.xls"))
        """

    def rm_tmp_file(self):
        """
        去除中间过程文件
        :return:
        """
        self.rm_dir(os.path.join(self.work_dir, 'MetaGenomic'), [])
        self.rm_dir(os.path.join(self.work_dir, 'QcAndStat'), [])
        self.rm_dir(os.path.join(self.work_dir, 'BwaRemoveHost'), [])
        self.rm_dir(os.path.join(self.work_dir, 'MgAssIdba'), [])  # 检查混拼
        self.rm_dir(os.path.join(self.work_dir, 'MgAssSoapdenovo'), [])  # 需检查
        self.rm_dir(os.path.join(self.work_dir, 'GenePredict'), ['output'])  # 中间文件需确定
        self.rm_dir(os.path.join(self.work_dir, 'UniGene'), ['gene.uniGeneset.fa.cd-hit-para-tmp'])
        self.rm_dir(os.path.join(self.work_dir, 'MetaDiamond'), [])
        self.rm_dir(os.path.join(self.work_dir, 'MetaDiamond1'), [])
        self.rm_dir(os.path.join(self.work_dir, 'MetaDiamond2'), [])
        self.rm_dir(os.path.join(self.work_dir, 'output'), ['qc', 'rm_host'])

    def rm_dir(self, dir, expect_list):
        if not os.path.exists(dir):
            self.logger.info('there is no dir named:' + dir)
            return
        elif os.path.isdir(dir):
            dir_list = os.listdir(dir)
            for file in dir_list:
                self.logger.info(dir + 'contains file:' + 'file')
                if file in expect_list:
                    self.logger.info('file: ' + file + 'in expect_list, list is :')
                    self.logger.info(expect_list)
                    continue
                else:
                    '''
                    for expect_file in expect_list:
                        if expect_file in file:
                            self.logger.info('file ' + file + 'is in expect_file' + expect_file)
                            rm_ = 0
                            break
                        else:
                            rm_ = 1
                    '''
                    rm_ = 1
                    if rm_ == 1:
                        next_dir = os.path.join(dir, file)
                        self.logger.info('now go into next dir' + next_dir)
                        self.rm_dir(next_dir, expect_list)
            check_dir = os.listdir(dir)
            if len(check_dir) == 0:
                self.logger.info('dir is empty , remove it :' + dir)
                shutil.rmtree(dir)
        elif os.path.isfile(dir):
            if dir in expect_list:
                self.logger.info("file is in expect_list, file is : " + dir + " expect_list is ")
                self.logger.info(expect_list)
                return
            else:
                file = os.path.basename(dir)
                rm_ = 1
                for expect_file in expect_list:
                    if expect_file in file:
                        self.logger.info('file' + file + 'is in expect_file' + expect_file)
                        rm_ = 0
                        break
                    else:
                        rm_ = 1
                if rm_ == 1:
                    filesize = os.path.getsize(dir) / 1024 / 1024
                    if filesize > 50:
                        os.remove(dir)
                        self.logger.info('now remove file: ' + dir + 'its size is : ')
                        self.logger.info(filesize)

    def set_output_all(self):
        """
        将所有结果一起导出至output，暂不需要
        """
        pass

    def move_dir(self, olddir, newname):  # 原函数名move2outputdir
        """
        移动一个目录下所有文件/文件夹到workflow输出路径下，供set_output调用
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="12800101")
        newdir = os.path.join(self.output_dir, newname)
        self.logger.info("newdir is : " + newdir)
        if not os.path.exists(newdir):
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
        self.logger.info("文件夹{}移动到{},耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        """
        递归移动文件夹的内容，供move_dir调用
        """
        if os.path.isfile(old_file):
            if not os.path.isdir(os.path.dirname(new_file)):
                os.makedirs(os.path.dirname(new_file))
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
        if self.taxon_outset:
            wait_steps = []
            while 1:
                self.logger.info("检查物种注释运行状态: {}".format(self.taxon_outset))
                status = [True if t.is_end else False for t in self.taxon_outset.values()]
                if all(status):
                    break
                else:
                    gevent.sleep(5)
        #self.wf_end()
        if self.option('test') == 'True':
            self.run_api(test=True)
        else:
            self.run_api()
        self.save_pdf()


    def wf_end(self):
        if self.sheet.id in ["i-sanger_183654"]:
            self.send_files()
            super(MetaGenomicWorkflow, self).end()
            return
        self.send_files()
        # self.rm_tmp_file() #暂时停用，防止上传结果出错退出
        super(MetaGenomicWorkflow, self).end()

    def save_pdf(self):
        if self.option("save_pdf") == 1:
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.wf_end)
            self.figsave.set_options({
                "task_id": self.sheet.id,
                "project": "metagenomic",
                "interaction": 0,
            })
            self.figsave.run()
        else:
            self.wf_end()

    def send_files(self):  # 原modify_output
        """
        结果放置到/upload_results
        """
        if self.option("save_pdf") == 1:
        # if self.option("save_pdf"):
            os.system("cp -r {}/* {}/".format(self.figsave.output_dir, self.output_dir))
        dir_o = self.output_dir
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        if os.path.exists(os.path.join(dir_o, "anno_personal")):
            self.move_dir(os.path.join(dir_o, "anno_personal"), os.path.join(dir_up, "anno_personal"))
        if os.path.exists(os.path.join(dir_o, "optimize_data")):
            self.move_dir(os.path.join(dir_o, "optimize_data"), os.path.join(dir_up, "optimize_data"))
        if os.path.exists(os.path.join(dir_o, "rawdata")) and self.option("qc"):
            self.move_file(os.path.join(dir_o, "rawdata/base_info.txt"), os.path.join(dir_up, "rawdata/base_info.txt"))
            self.move_file(os.path.join(dir_o, "rawdata/reads.rawData.stat"),
                           os.path.join(dir_up, "rawdata/reads.rawData.stat.xls"))
            self.move_file(os.path.join(dir_o, "base_quality_distribution.pdf.tar.gz"),
                           os.path.join(dir_up, "rawdata/base_quality_distribution.pdf.tar.gz"))  # pdf压缩文件夹
            self.move_file(os.path.join(dir_o, "qc/reads.cleanData.stat.xls"),
                           os.path.join(dir_up, "qc/qc_stat/reads.cleanData.stat.xls"))
        if os.path.exists(os.path.join(dir_o, "rm_host")):
            self.move_file(os.path.join(dir_o, "rm_host/stat.list.txt"), os.path.join(dir_up, "rm_host/stat.list.txt"))
        if os.path.exists(os.path.join(dir_o, "assemble")):
            self.move_file(os.path.join(dir_o, "assemble/assembly.stat"),
                           os.path.join(dir_up, "assemble/assembly.stat.xls"))
            if os.path.exists(os.path.join(dir_o, "contig.length.pdf.tar.gz")):
                self.move_file(os.path.join(dir_o, "contig.length.pdf.tar.gz"),
                               os.path.join(dir_up, "assemble/contig.length.pdf.tar.gz"))  # pdf压缩文件夹
            else:
                self.move_file(os.path.join(dir_o, "All.contig.length.pdf"),
                               os.path.join(dir_up, "assemble/All.contig.length.pdf"))  # pdf图片
            files = os.listdir(os.path.join(dir_o, "assemble"))
            for file in files:
                if file.endswith("contig.fa.tar.gz"):
                    self.move_file(os.path.join(dir_o, "assemble", file), os.path.join(dir_up, "assemble", file))
            with open (os.path.join(dir_up, "assemble", "list.txt"), "w") as f:
                f.write("name\tfiles\n")
                for file in files:
                    if file.endswith("contig.fa.tar.gz"):
                        m = re.search('(.*).contig.fa.tar.gz', file)
                        f.write("{}\t{}\n".format(m.group(1), file))
        if os.path.exists(os.path.join(dir_o, "predict")):
            self.move_file(os.path.join(dir_o, "predict/sample.metagene.stat"),
                           os.path.join(dir_up, "predict/genePredict_stat.xls"))
            if os.path.exists(os.path.join(dir_o, "gene.length.pdf.tar.gz")):
                self.move_file(os.path.join(dir_o, "gene.length.pdf.tar.gz"),
                               os.path.join(dir_up, "predict/gene.length.pdf.tar.gz"))  # pdf压缩文件夹
            else:
                self.move_file(os.path.join(dir_o, "Total.gene.length.pdf"),
                               os.path.join(dir_up, "predict/Total.gene.length.pdf"))  # pdf图片
            files = os.listdir(os.path.join(dir_o, "predict"))
            for file in files:
                new_file = file.replace("metagene.more100.fa", "genePredict.fa")
                if file.endswith("tar.gz"):
                    self.move_file(os.path.join(dir_o, "predict", file), os.path.join(dir_up, "predict", new_file))
        if os.path.exists(os.path.join(dir_o, "metaphlan")):
            self.move_dir(os.path.join(dir_o, "metaphlan"), os.path.join(dir_up, "metaphlan"))
        if os.path.exists(os.path.join(dir_o, "kraken")):
            self.move_dir(os.path.join(dir_o, "kraken"), os.path.join(dir_up, "kraken"))
        if os.path.exists(os.path.join(dir_o, "geneset")):
            self.move_file(os.path.join(dir_o, "geneset/uniGeneset"), os.path.join(dir_up, "geneset/uniGeneset"))
            self.move_file(os.path.join(dir_o, "geneset/gene_profile"), os.path.join(dir_up, "geneset/gene_profile"))
            if os.path.exists(os.path.join(dir_o, "geneCatalog.length.pdf")):
                self.move_file(os.path.join(dir_o, "geneCatalog.length.pdf"),
                               os.path.join(dir_up, "geneset/geneCatalog.length.pdf"))  # pdf压缩文件夹
            if self.sheet.id in ["i-sanger_183654"]:
                overviewfile = self.wait_file(os.path.join(dir_o, "geneset/anno_overview.xls"))  # 等待anno_overview转移至结果目录下
            else:
                overviewfile = self.overview.option('gene_overview').prop['path']  # 使用tool的结果至接过目录下 by GHD @20180806
            self.move_file(overviewfile, os.path.join(dir_up, "geneset/anno_overview.xls"))  # modified by GHD @ 20180425
            # self.move_file(os.path.join(dir_o, "geneset/anno_overview.xls"),
            #                os.path.join(dir_up, "geneset/anno_overview.xls"))
        for file in ["nr", "kegg", "cog", "cazy", "vfdb", "ardb", "card", "composition", "compare", "correlation","nr_lca","nr_deunclassied"]:
            if os.path.exists(os.path.join(dir_o, file)):
                self.move_file(os.path.join(dir_o, file), os.path.join(dir_up, file))
            if file == "vfdb":
                if os.path.exists(os.path.join(dir_o, "vfdb_level2.pdf.tar.gz")):
                    self.move_file(os.path.join(dir_o, "vfdb_level2.pdf.tar.gz"),
                                   os.path.join(dir_up, "vfdb/vfdb_level2.pdf.tar.gz"))  # pdf压缩文件夹
            if file == "kegg":
                # for rm_file in ["kegg_merge.xml", "pathway.kgml", "pathway.png"]:  # 保留kegg_merge.xml
                for rm_file in ["pathway.kgml", "pathway.png", "kegg_merge.xml"]:  # 不再保留
                    rm_file_path = os.path.join(dir_up, file, rm_file)
                    if os.path.isfile(rm_file_path):
                        os.remove(rm_file_path)
            if file == 'cazy':
                rm_file_path = os.path.join(dir_up, file, "gene_cazy_parse_anno.xls")
                if os.path.isfile(rm_file_path):
                    os.remove(rm_file_path)
            """
            if file == "composition":
                file_names = ["nr_pdf", "cog_pdf", "kegg_pdf", "cazy_pdf", "ardb_pdf", "card_pdf", "vfdb_pdf", "gene_pdf", "personal_pdf"]
                for tmp_file in file_names:
                    anno = tmp_file[:-4]
                    if os.path.exists(os.path.join(dir_o, tmp_file)):
                        if os.path.exists(os.path.join(dir_o, tmp_file, "one_sample_bar.pdf")):
                            shutil.rmtree(os.path.join(dir_o, tmp_file, "one_sample_bar.pdf"))
                        elif os.path.exists(os.path.join(dir_o, tmp_file, "pie.pdf")):
                            shutil.rmtree(os.path.join(dir_o, tmp_file, "pie.pdf"))
                        self.move_dir(os.path.join(dir_o, tmp_file), os.path.join(dir_up, file, "bar", anno))
            """


        repaths = [
            [".", "", "流程分析结果目录", 0, "120001"],
            ["rawdata", "", "原始序列目录", 0, "120002"],
            ["rawdata/base_info.txt", "", "原始序列质量统计表", 0, "120003"],
            ["rawdata/reads.rawData.stat.xls", "", "各样品原始数据统计表", 0, "120006"],
            ["rawdata/base_quality_distribution.pdf.tar.gz", "", "碱基质量图及碱基分布图"],
            ["qc", "", "质控结果目录", 0, "120007"],
            ["qc/reads.cleanData.stat.xls", "", "各样品质控后数据统计表", 0, "120009"],
            ["rm_host", "", "去宿主后序列目录", 0, "120010"],
            ["rm_host/stat.list.txt", "", "各样品去宿主后序列统计表", 0, "120011"],
            ["assemble", "", "拼接结果目录", 0, "120012"],
            ["assemble/assembly.stat.xls", "", "各步骤组装结果统计表", 0, "120013"],
            ["assemble/contig.length.pdf.tar.gz", "", "Contig长度分布图"],
            ["assemble/All.contig.length.pdf", "pdf", "Contig长度分布图"],
            ["predict", "", "基因预测结果目录", 0, "120015"],
            ["predict/genePredict_stat.xls", "", "基因预测结果统计表", 0, "120016"],
            ["predict/gene.length.pdf.tar.gz", "", "各样本预测基因长度分布图"],
            ["predict/Total.gene.length.pdf", "pdf", "预测基因长度分布图"],
            ["geneset", "", "非冗余基因集结果目录", 0, "120018"],
            ["geneset/anno_overview.xls", "", "物种与功能注释总览表", 0, "120019"],
            ["geneset/gene_profile", "", "非冗余集因丰度目录", 0, "120020"],
            ["geneset/gene_profile/gene.uniGeneset.fa.length.txt", "", "基因长度表", 0, "120034"],
            ["geneset/gene_profile/RPKM.xls", "", "基因在各个样品中的RPKM丰度表", 0, "120025"],
            ["geneset/gene_profile/reads_length_ratio_relative.xls", "", "基因在各个样品中的相对丰度/基因长度表", 0, "120026"],
            ["geneset/gene_profile/reads_length_ratio.xls", "", "基因在各个样品中的丰度/基因长度表", 0, "120027"],
            ["geneset/gene_profile/reads_number_relative.xls", "", "基因在各个样品中的相对丰度表", 0, "120028"],
            ["geneset/gene_profile/reads_number.xls", "", "基因在各个样品中的丰度表", 0, "120029"],
            ["geneset/gene_profile/TPM.xls", "", "基因在各个样品中的TPM丰度表", 0, "120030"],
            ['geneset/gene_profile/top100_reads_number.xls', '', '丰度前100的基因丰度表', 0, "120031"],
            ['geneset/gene_profile/top100_reads_number_relative.xls', '', '丰度前100的基因相对丰度表', 0, "120032"],
            ['geneset/gene_profile/reads_profile.tar.gz', '', '基因reads数相对丰度与基因reads数丰度的压缩文件', 0, "120033"],
            ["geneset/uniGeneset", "", "非冗余基因集序列统计目录", 0, "120021"],
            ["geneset/uniGeneset/geneCatalog_stat.xls", "", "非冗余基因数目和长度统计表", 0, "120022"],
            ["geneset/uniGeneset/gene.uniGeneset.fa", "", "非冗余基因集核酸序列", 0, "120023"],
            ["geneset/uniGeneset/gene.uniGeneset.faa", "", "非冗余基因集蛋白序列", 0, "120024"],
            ["geneset/geneCatalog.length.pdf", "pdf", "非冗余基因集长度分布图"],
            ["nr", "", "NR功能注释结果目录", 0, "120035"],
            ["nr/gene_nr_anno.xls", "", "每条基因的物种注释表", 0, "120036"],
            ["nr/nr_align_table.xls", "", "物种序列比对结果", 0, "120037"],
            ["nr/tax_d.xls", "", "域注释丰度表", 0, "120038"],
            ["nr/tax_k.xls", "", "界注释丰度表", 0, "120039"],
            ["nr/tax_p.xls", "", "门注释丰度表", 0, "120040"],
            ["nr/tax_c.xls", "", "纲注释丰度表", 0, "120041"],
            ["nr/tax_o.xls", "", "目注释丰度表", 0, "120042"],
            ["nr/tax_f.xls", "", "科注释丰度表", 0, "120043"],
            ["nr/tax_g.xls", "", "属注释丰度表", 0, "120044"],
            ["nr/tax_s.xls", "", "种注释丰度表", 0, "120045"],
            ["kegg", "", "KEGG功能注释结果目录", 0, "120046"],
            ["kegg/gene_kegg_anno.xls", "", "每条基因的KEGG功能注释表", 0, "120047"],
            ["kegg/kegg_align_table.xls", "", "KEGG序列比对结果表", 0, "120048"],
            ["kegg/kegg_pathway_eachmap.xls", "", "Pathway在各个样品中的丰度表", 0, "120049"],
            ["kegg/kegg_enzyme_profile.xls", "", "各样品KEGG酶丰度表", 0, "120050"],
            ["kegg/kegg_gene_profile.xls", "", "各样品KEGG基因丰度表", 0, "120051"],
            ["kegg/kegg_KO_profile.xls", "", "各样品KO丰度表", 0, "120052"],
            ["kegg/kegg_level1_profile.xls", "", "各样品Pathway level1丰度表", 0, "120053"],
            ["kegg/kegg_level2_profile.xls", "", "各样品Pathway level2丰度表", 0, "120054"],
            ["kegg/kegg_level3_profile.xls", "", "各样品Pathway level3丰度表", 0, "120055"],
            ["kegg/kegg_module_profile.xls", "", "各样品KEGG Module丰度表", 0, "120056"],
            ["kegg/kegg_pathway_profile.xls", "", "各样品Pathway丰度表", 0, "120057"],
            ["kegg/pathway_img.tar.gz", "", "", 1, ""],  # 图片路径，临时添加
            ["kegg/pathway_img.tar.gz/pathway_img.pdf", "pdf", ""],
            ["cog", "", "COG功能注释结果目录", 0, "120058"],
            ["cog/cog_align_table.xls", "", "COG序列比对结果表", 0, "120059"],
            ["cog/cog_nog_profile.xls", "", "各样品COG NOG丰度表", 0, "120060"],
            ["cog/cog_category_profile.xls", "", "各样品COG Category丰度表", 0, "120061"],
            ["cog/cog_function_profile.xls", "", "各样品COG Function丰度表", 0, "120062"],
            ["cog/gene_cog_anno.xls", "", "每条基因的COG功能注释表", 0, "120063"],
            ["cazy", "", "CAZy碳水化合物活性酶注释结果目录", 0, "120064"],
            ["cazy/cazy_class_profile.xls", "", "各样品CAZY Class丰度表", 0, "120065"],
            ["cazy/gene_cazy_anno.xls", "", "每条基因的CAZY功能注释表", 0, "120066"],
            ["cazy/gene_cazy_class_stat.xls", "", "Class基因信息统计表", 0, "120069"],
            ["cazy/gene_cazy_family_stat.xls", "", "Family基因信息统计表", 0, "120070"],
            ["cazy/cazy_family_profile.xls", "", "各样品CAZY Family丰度表", 0, "120067"],
            ["cazy/gene_dbCAN.hmmscan.out.dm.ds", "", "CAZY序列比对结果表", 0, "120068"],
            ["vfdb", "", "VFDB毒力因子注释结果目录", 0, "120071"],
            ["vfdb/gene_vfdb_core_anno.xls", "", "每条基因的VFDB核心库功能注释表", 0, "120072"],
            ["vfdb/gene_vfdb_predict_anno.xls", "", "每条基因的VFDB预测库功能注释表", 0, "120073"],
            ["vfdb/gene_vfdb_total_anno.xls", "", "每条基因的VFDB功能注释表", 0, "120074"],
            ["vfdb/vfdb_all_Gi_profile.xls", "", "各样品VFDB基因丰度表", 0, "120075"],
            ["vfdb/vfdb_all_VF_profile.xls", "", "各样品VFDB毒力因子丰度表", 0, "120076"],
            ["vfdb/vfdb_core_align_table.xls", "", "VFDB核心库序列比对结果", 0, "120077"],
            ["vfdb/vfdb_core_Gi_profile.xls", "", "各样品VFDB核心库基因丰度表", 0, "120078"],
            ["vfdb/vfdb_core_VF_profile.xls", "", "各样品VFDB核心库毒力因子丰度表", 0, "120079"],
            ["vfdb/vfdb_predict_align_table.xls", "", "VFDB预测库序列比对结果表", 0, "120081"],
            ["vfdb/vfdb_predict_Gi_profile.xls", "", "各样品VFDB预测库基因丰度表", 0, "120082"],
            ["vfdb/vfdb_predict_VF_profile.xls", "", "各样品VFDB预测库毒力因子丰度表", 0, "120083"],
            ["vfdb/vfdb_level_pie.xls", "", "VFDB两级分类的丰度统计表", 0, "120080"],
            ["vfdb/vfdb_level2.pdf.tar.gz", "", "毒力基因预测分类统计饼图"],
            ["ardb", "", "ARDB抗性基因功能注释结果目录", 0, "120084"],
            ["ardb/ardb_align_table.xls", "", "ARDB序列比对结果表", 0, "120085"],
            ["ardb/ardb_class_profile.xls", "", "各样品ARDB Class丰度表", 0, "120086"],
            ["ardb/gene_ardb_anno.xls", "", "每条基因的ARDB功能注释表", 0, "120087"],
            ["ardb/ardb_ARG_profile.xls", "", "各样品ARDB ARG丰度表", 0, "120088"],
            ["ardb/ardb_type_profile.xls", "", "各样品ARDB type丰度表", 0, "120089"],
            ['ardb/gene_ardb_class_stat.xls', '', 'class基因信息统计表', 0, "120090"],
            ['ardb/ardb_Antibiotic_class_profile.xls', '', '各样品Antibiotic_class丰度表', 0, ""],
            ["card", "", "CARD抗性基因功能注释结果目录", 0, "120091"],
            ["card/card_align_table.xls", "", "CARD序列比对结果", 0, "120092"],
            ["card/card_class_profile.xls", "", "各样品CARD Class丰度表", 0, "120093"],
            ["card/card_ARO_gene_number.xls", "", "CARD每个ARO比对基因信息表", 0, "120094"],
            ["card/gene_card_anno.xls", "", "每条基因的CARD功能注释表", 0, "120095"],
            ["card/card_ARO_profile.xls", "", "各样品CARD ARO丰度表", 0, "120096"],
            ["card/card_ARO_profile_all.xls", "", "各样品ARO丰度总表",0,"120350"],
            ["card/card_Resistance_Mechanism_profile.xls", "", "各样品抗性机制丰度表",0,"120351"],
            ["card/card_Antibiotic_class_profile.xls", "", "各样品Antibiotic class丰度表",0,"120352"],
            ["card/card_Drug_class_profile.xls", "", "各样品Drug_class丰度表",0,"120353"],
            ["composition", "", "物种与功能组成分析结果目录", 0 , "120097"],
            ["composition/bar", "", "柱形图结果目录", 0, "120098"],
            ["composition/heatmap", "", "Heatmap图结果目录", 0, "120099"],
            ["composition/circos", "", "Circos样本与物种或功能关系图结果目录", 0, "120100"],
            ["compare", "", "物种与功能比较分析结果目录", 0, "120103"],
            ["compare/Hcluster", "", "样本层次聚类分析结果目录", 0, "120104"],
            ["compare/Distance", "", "距离矩阵计算结果目录", 0, "120106"],
            ["compare/Pca", "", "PCA分析结果目录", 0, "120108"],
            ["compare/Pcoa", "", "PCoA分析结果目录", 0, "120109"],
            ["compare/Nmds", "", "NMDS分析结果目录", 0, "120110"],
            ["correlation", "", "环境因子关联分析结果目录", 0, "120124"],
            ["correlation/Rda", "", "RDA_CCA分析结果目录", 0, "120125"],
            ["correlation/Dbrda", "", "db_RDA分析结果目录", 0, "120126"],
            ["correlation/cor_heatmap", "", "相关性Heatmap分析结果目录", 0, "120127"],
            ["metaphlan", "", "metaphlan3注释结果目录", 0, ""],
            ["metaphlan/taxon_D.xls", "", "域注释丰度表", 0, ""],
            ["metaphlan/taxon_K.xls", "", "界注释丰度表", 0, ""],
            ["metaphlan/taxon_P.xls", "", "门注释丰度表", 0, ""],
            ["metaphlan/taxon_C.xls", "", "纲注释丰度表", 0, ""],
            ["metaphlan/taxon_O.xls", "", "目注释丰度表", 0, ""],
            ["metaphlan/taxon_F.xls", "", "科注释丰度表", 0, ""],
            ["metaphlan/taxon_G.xls", "", "属注释丰度表", 0, ""],
            ["metaphlan/taxon_S.xls", "", "种注释丰度表", 0, ""],
            ["kraken", "", "kraken2注释结果目录", 0, ""],
            ["kraken/taxon_D.xls", "", "域注释丰度表", 0, ""],
            ["kraken/taxon_K.xls", "", "界注释丰度表", 0, ""],
            ["kraken/taxon_P.xls", "", "门注释丰度表", 0, ""],
            ["kraken/taxon_C.xls", "", "纲注释丰度表", 0, ""],
            ["kraken/taxon_O.xls", "", "目注释丰度表", 0, ""],
            ["kraken/taxon_F.xls", "", "科注释丰度表", 0, ""],
            ["kraken/taxon_G.xls", "", "属注释丰度表", 0, ""],
            ["kraken/taxon_S.xls", "", "种注释丰度表", 0, ""],
        ]
        regexps = [
            [r"rawdata/base_info.txt", "", "原始序列碱基质量统计文件", 0, "120004"],
            [r"assemble/[^/]+\.contig\.fa", "", "拼接contig结果", 0, "120014"],
            [r"predict/.+\.genePredict.fa.tar.gz", "", "长度大于等于100bp的基因的核酸序列", 0, "120017"],
            [r"geneset/gene_profile/top100\..+xls", "", "丰度前100的基因丰度表", 0, "120031"],
            [r"composition/.+/.+/taxa\.percents\.table\.xls", "", "物种/基因/功能相对丰度结果表", 0, "120101"],
            [r"composition/.+/.+/taxa\.table\.xls", "", "物种/基因/功能丰度结果表", 0, "120102"],
            [r"composition/.+/.+/bar\.pdf", "", "群落柱形图", 0, ""],
            [r"composition/.+/.+/one_sample_bar\.pdf\.tar\.gz", "", "单样本柱形图", 0, ""],
            [r"composition/.+/.+/pie\.pdf\.tar\.gz", "", "单样本饼图", 0, ""],
            [r"compare/Hcluster/.+/hcluster\.tre", "graph.newick_tree", "样本层次聚类树结果表", 0, "120105"],
            [r"compare/Distance/.+/bray_curtis.+\.xls$", "meta.beta_diversity.distance_matrix", "样本距离矩阵文件", 0, "120107"],
            [r"compare/Nmds/.+/nmds_sites\.xls$", "xls", "NMDS样本各维度坐标", 0, "120111"],
            [r"compare/Nmds/.+/nmds_stress\.xls$", "xls", "NMDS样本特征拟合度值", 0, "120112"],
            [r"compare/Pca/.+/pca_importance\.xls$", "xls", "PCA主成分解释度表", 0, "120113"],
            [r"compare/Pca/.+/pca_sites\.xls$", "xls", "PCA样本各成分轴坐标", 0, "120114"],
            [r"compare/Pca/.+/pca_rotation\.xls$", "xls", "PCA主成分贡献度表", 0, "120115"],
            [r"compare/Pca/.+/pca_rotation_all\.xls$", "xls", "PCA全部主成分贡献度表", 0, "120116"],
            [r"compare/Pca/.+/pca_envfit_vector\.xls$", "xls", "PCA数量型环境因子的显著性检验值", 0, "120120"],
            [r"compare/Pca/.+/pca_envfit_vector_scores\.xls$", "xls", "PCA数量型环境因子坐标表", 0, "120119"],
            [r"compare/Pca/.+/pca_envfit_factor\.xls$", "xls", "PCA哑变量环境因子的显著性检验值",0,"120231"],
            [r"compare/Pca/.+/pca_envfit_factor_scores\.xls$", "xls", "PCA哑变量环境因子坐标表", 0, "120117"],
            [r"compare/Pcoa/.+/pcoa_sites\.xls", "xls", "PCoA样本坐标表", 0, "120123"],
            [r"compare/Pcoa/.+/pcoa_eigenvalues\.xls", "xls", "PCoA矩阵特征值", 0, "120121"],
            [r"compare/Pcoa/.+/pcoa_eigenvaluespre\.xls", "xls", "PCoA特征解释度百分比", 0, "120122"],
            [r"correlation/Rda/.+/dca\.xls", "xls", "判断用RDA还是CCA的DCA文件", 0, "120128"],
            [r"correlation/Rda/.+/.*rda_biplot\.xls", "xls", "RDA_CCA数量型环境因子坐标表", 0, "120129"],
            [r"correlation/Dbrda/.+/.*rda_biplot\.xls", "xls", "dbRDA数量型环境因子坐标表", 0, "120135"],
            [r"correlation/Rda/.+/rda_envfit\.xls", "xls", "RDA_CCA各环境因子的显著性检验值", 0, "120130"],
            [r"correlation/Dbrda/.+/.*rda_envfit\.xls", "xls", "dbRDA各环境因子的显著性检验值", 0, "120136"],
            [r"correlation/Rda/.+/.*rda_importance\.xls", "xls", "RDA_CCA主成分解释度表", 0, "120131"],
            [r"correlation/Dbrda/.+/.*rda_importance\.xls", "xls", "dbRDA主成分解释度表", 0, "120137"],
            [r"correlation/Rda/.+/.*rda_plot_species_data\.xls", "xls", "RDA_CCA绘图物种_功能坐标表", 0, "120132"],
            [r"correlation/Dbrda/.+/.*rda_plot_species_data\.xls", "xls", "dbRDA绘图物种_功能坐标表", 0, "120138"],
            [r"correlation/Rda/.+/.*rda_sites\.xls", "xls", "RDA_CCA样本坐标表", 0, "120133"],
            [r"correlation/Dbrda/.+/.*rda_sites\.xls", "xls", "dbRDA样本坐标表", 0, "120139"],
            [r"correlation/Rda/.+/.*rda_species\.xls", "xls", "RDA_CCA物种_功能坐标表", 0, "120134"],
            [r"correlation/Dbrda/.+/.*rda_species\.xls", "xls", "dbRDA物种_功能坐标表", 0, "120140"],
            [r"correlation/cor_heatmap/.+/correlation\.xls", "xls", "相关系数表", 0, "120141"],
            [r"correlation/cor_heatmap/.+/pvalue\.xls", "xls", "相关系数对应p值表", 0, "120142"],
            [r"metaphlan/.*\.taxon\.xls", "xls", "样本metaphlan3注释文件", 0, ""],
            [r"kraken/.*\.taxon\.xls", "xls", "样本kraken2注释文件", 0, ""],
            [r"kraken/.*\.taxon_bracken\.xls", "xls", "样本bracken结果文件", 0, ""],
        ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

    '''导表'''

    def run_api(self, test=False):  # 原run_api_and_set_output
        self.logger.info("开始导表")
        if self.pipeline == "pipeline3":
            gevent.sleep(10)
            self.export_qc()
        elif self.pipeline != 'pipeline2.1':
            self.export_qc()
            if self.pipeline == "pipeline1":
                self.export_assem()
                self.export_predict()
            self.export_geneset()
        self.update_sg_task()  # 新增数据库版本的更新时间用以判断数据库更新，每个对应的版本在主表中有存
        self.insert_soft_db()
        if len(self.all_anno) == 0:
            self.logger.info("导表完成")
            return
        if "nr" in self.choose_anno:
            self.export_nr_deunclass()
            self.export_nr_lca()
            self.export_nr()
        if "kegg" in self.choose_anno:
            self.export_kegg()
        if "cog" in self.choose_anno:
            self.export_cog()
        if "cazy" in self.choose_anno:
            self.export_cazy()
        if "vfdb" in self.choose_anno:
            self.export_vfdb()
        if "ardb" in self.choose_anno:
            self.export_ardb()
        if "card" in self.choose_anno:
            self.export_card()
        if self.option("kraken"):
            self.export_kraken()
        if self.option("metaphlan"):
            self.export_metaphlan()
        self.export_overview()
        self.geneset_id = str(self.geneset_id)
        self.env_id = str(self.env_id)
        ana_list = ["composition", "compare", "correlation"]
        for ana in ana_list:
            ana_path = os.path.join(self.output_dir, ana)
            if os.path.isdir(ana_path):
                analysis_list = os.listdir(ana_path)
                for analysis in analysis_list:
                    anno_list = os.listdir(os.path.join(self.output_dir, ana, analysis))
                    if ana == "composition":
                        self.export_composition(analysis, anno_list)
                    if ana == "compare":
                        if analysis == "Hcluster":
                            self.export_hclust(anno_list)
                        elif analysis in ["Nmds", "Pca", "Pcoa"]:
                            self.export_beta(analysis, anno_list)
                    if ana == "correlation":
                        if analysis == "cor_heatmap":
                            self.export_cor_heatmap(anno_list)
                        elif analysis in "Dbrda, Rda":
                            self.logger.info("analysis is :" + analysis)
                            self.export_beta(analysis, anno_list)
        if self.option("personal_anno") != "none":
            # 导入个性化注释主表
            self.export_personal_anno()
        self.export_software()
        self.logger.info("导表完成")

    def insert_soft_db(self):
        """
        插入本次工作流运行使用的软件和数据库信息
        """
        soft_db = self.api.api("metagenomic.common_api")
        soft_db.insert_soft_db_task(self.sheet.id, self.wf_soft, self.wf_database)

    def update_sg_task(self):
        """
        兼容数据库版本，先更新sg_task版本，以此作为判断
        :return:
        """
        self.logger.info("开始更新sg_task表!")
        task_id = "_".join(self._sheet.id.split("_")[0:2])
        software = self.api.api('metagenomic.software')
        database_type = json.dumps({"update_time": "202011"},sort_keys=True, separators=(',', ':'))
        software.update_mongo('sg_task',{"task_id":task_id}, {"database":database_type,
                                                              "wf_version": self.wf_version,
                                                              "sample_num": len(self.spname_spid), "save_pdf": self.option('save_pdf')})
        self.logger.info("更新sg_task表完成!")

    @time_count
    def export_specimen_group(self):
        '''
        工作流2单独导specimen_group
        '''
        if self.option('group').is_set:
            self.api_dic["group"] = self.api.api("metagenomic.specimen_group")
            group_info = self.api_dic["group"].add_ini_group_table(self.option('group').prop['path'],"")
            if len(group_info) == 1:
                self.specimen_group = group_info[0]["specimen_group"]
                self.group_detail = group_info[0]["group_detail"]
            else:
                self.logger.info("group_info length is %s" % len(group_info))
                return
        else:
            self.specimen_group = 'all'
            self.group_detail = {'all': sorted(self.spname_spid.values())}
        if self.option('envtable').is_set:
            self.api_dic["envtable"] = self.api.api("metagenomic.env_metagenomic")
            self.env_id = self.api_dic["envtable"].add_env_table(self.option('envtable').prop['path'], self.spname_spid)
            self.env_labs = self.api_dic["envtable"].get_env_lab()

    @time_count
    def export_qc(self):
        self.api_dic["data_stat"] = self.api.api("metagenomic.data_stat")
        if self.pipeline in ["pipeline1", "pipeline3"]:
            soft_db = {"soft_db_info": self.soft_db_info["qc"]}
            if self.option('qc'):
                data_stat_id = self.api_dic["data_stat"].add_data_stat("raw",
                                                                       self.sequence.raw_stat,
                                                                       self.output_dir + "/rawdata/base_info.txt", "null",
                                                                       software_ver=soft_db)
                self.api_dic["data_stat"].add_data_stat("clean", self.output_dir + "/qc/reads.cleanData.stat.xls",
                                                        "null", data_stat_id, software_ver=soft_db)
            else:
                #data_stat_id = self.api_dic["data_stat"].add_data_stat('raw', self.option('raw_info').prop['path'],
                                                                       #self.base_info, "null")  # 暂不存在碱基统计路径
                if self.option('raw_info').is_set:
                    data_stat_id = self.api_dic["data_stat"].add_data_stat('raw', self.option('raw_info').prop['path'],
                                                                           "null", "null", software_ver=soft_db)  # 暂不存在碱基统计路径
                else:
                    data_stat_id = self.api_dic["data_stat"].add_data_stat('raw', self.sequence.raw_stat,
                                                                           "null", "null", software_ver=soft_db)
                if self.option("qc_info").is_set:
                    self.api_dic["data_stat"].add_data_stat("clean", self.option('qc_info').prop['path'],
                                                            "null", data_stat_id, software_ver=soft_db)
        elif self.pipeline == "pipeline2":
            stat_path = self.output_dir + "/optimize_data/optimize_reads.stat.xls"
            base_info_path = self.output_dir + "/optimize_data/base_info"
            #data_stat_id = self.api_dic["data_stat"].add_data_stat("unigene", self.option("insertsize").prop["path"],"null","null")
            data_stat_id = self.api_dic["data_stat"].add_data_stat("unigene", stat_path,base_info_path,"null",
                                                                   software_ver=soft_db)
        self.spname_spid = name2id(data_stat_id, type="raw")
        if self.pipeline == "pipeline3":
            return
        if self.option('rm_host'):
            if self.option('ref_database') in ["", "Custom"]:
                params = "Custom," + self.option('ref_undefined_name')
            else:
                params = self.option('ref_database')
            software_ver = {'soft_db_info': {"soft":19}}  # 增加软件版本信息
            self.api_dic["data_stat"].add_data_stat("optimised", self.output_dir + "/rm_host/stat.list.txt", "null",
                                                    data_stat_id, software_ver=software_ver, params=params)
        if self.option('group').is_set:
            self.api_dic["group"] = self.api.api("metagenomic.specimen_group")
            group_info = self.api_dic["group"].add_ini_group_table(self.option('group').prop['path'], self.spname_spid)
            if len(group_info) == 1:
                self.specimen_group = group_info[0]["specimen_group"]
                self.group_detail = group_info[0]["group_detail"]
            else:
                self.logger.info("group_info length is %s" % len(group_info))
                return
        else:
            self.specimen_group = 'all'
            self.group_detail = {'all': sorted(self.spname_spid.values())}
        self.logger.info("specimen_group is :")
        self.logger.info(self.specimen_group)
        self.logger.info("group_detail is :")
        self.logger.info(self.group_detail)
        if self.option('envtable').is_set:
            self.api_dic["envtable"] = self.api.api("metagenomic.env_metagenomic")
            self.env_id = self.api_dic["envtable"].add_env_table(self.option('envtable').prop['path'], self.spname_spid)
            self.env_labs = self.api_dic["envtable"].get_env_lab()

    @time_count
    def export_assem(self):
        self.api_dic["assem_gene"] = self.api.api("metagenomic.assemble_gene")
        assemble_type = self.option('assemble_type')
        soft_db = {"soft_db_info": self.soft_db_info["assem"]}
        assem_id = self.api_dic["assem_gene"].add_assemble_stat(assemble_type, self.option("min_contig"),
                                                                software_ver=soft_db)
        self.api_dic["assem_gene"].add_assemble_stat_detail(assem_id, self.output_dir + "/assemble/assembly.stat")
        self.api_dic["assem_gene"].add_assemble_stat_bar(assem_id, self.output_dir + "/assemble/len_distribute")

    @time_count
    def export_predict(self):
        if "assem_gene" not in self.api_dic.keys():
            self.api_dic["assem_gene"] = self.api.api("metagenomic.assemble_gene")
        assemble_type = self.option('assemble_type')
        soft_db = {"soft_db_info": self.soft_db_info["gene_predict"]}
        gene_id = self.api_dic["assem_gene"].add_predict_gene(assemble_type, self.option("min_gene"),
                                                              software_ver=soft_db)
        self.api_dic["assem_gene"].add_predict_gene_detail(gene_id, self.output_dir + "/predict/sample.metagene.stat")
        self.api_dic["assem_gene"].add_predict_gene_bar(gene_id, self.output_dir + "/predict/len_distribute")
        self.api_dic["assem_gene"].add_predict_gene_total(self.output_dir + "/predict/sample.metagene.stat")

    @time_count
    def export_geneset(self):
        self.api_dic["geneset"] = self.api.api("metagenomic.geneset")
        geneset_path = os.path.join(self.remote_dir, 'geneset')
        soft_db = {"soft_db_info": self.soft_db_info["gene_set_n"]}
        soft_db2 = {"soft_db_info2": self.soft_db_info["gene_set_a"]}
        self.geneset_id = self.api_dic["geneset"].add_geneset(self.output_dir + "/geneset", geneset_path, 1,
                                                              identity=self.option('cdhit_identity'),
                                                              software_ver=soft_db,
                                                              software_ver2=soft_db2)  # mark add by qingchen.zhang@20190621
        self.api_dic["geneset"].add_geneset_bar(self.geneset_id, self.output_dir + "/geneset/length_distribute")
        self.api_dic["geneset"].add_geneset_readsn(self.geneset_id,
                                                   self.output_dir + "/geneset/gene_profile/top100_reads_number.xls")
        self.api_dic["geneset"].add_geneset_readsr(self.geneset_id,
                                                   self.output_dir + "/geneset/gene_profile/top100_reads_number_relative.xls")
        self.anno_id["gene"] = ""

    @time_count
    def export_nr(self):
        self.api_dic["nr"] = self.api.api("metagenomic.mg_anno_nr")
        nr_path = os.path.join(self.remote_dir, 'nr/gene_nr_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["nr"]}
        return_id = self.api_dic["nr"].add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group,
                                                   group_detail=self.group_detail,
                                                   software_ver=soft_db)
        self.api_dic["nr"].add_anno_nr_detail(return_id, self.output_dir + "/nr")
        self.anno_id["nr"] = str(return_id)

    @time_count
    def export_nr_lca(self):
        self.api_dic["nr"] = self.api.api("metagenomic.mg_anno_nr")
        nr_path = os.path.join(self.remote_dir, 'nr_lca/gene_nr_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["nr"]}
        return_id = self.api_dic["nr"].add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group,
                                                   group_detail=self.group_detail, name = "NR_Origin_LCA",
                                                   software_ver=soft_db)
        self.api_dic["nr"].add_anno_nr_detail(return_id, self.output_dir + "/nr_lca")
        self.anno_id["nr"] = str(return_id)

    @time_count
    def export_nr_deunclass(self):
        self.api_dic["nr"] = self.api.api("metagenomic.mg_anno_nr")
        nr_path = os.path.join(self.remote_dir, 'nr_deunclassied/gene_nr_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["nr"]}
        return_id = self.api_dic["nr"].add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group,
                                                   group_detail=self.group_detail, name = "NR_Origin_Deunclassified",
                                                   software_ver=soft_db)
        self.api_dic["nr"].add_anno_nr_detail(return_id, self.output_dir + "/nr_deunclassied")
        self.anno_id["nr"] = str(return_id)

    @time_count
    def export_cog(self):
        self.api_dic["cog"] = self.api.api("metagenomic.mg_anno_cog")
        cog_path = os.path.join(self.remote_dir, 'cog/gene_cog_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["cog"]}
        return_id = self.api_dic["cog"].add_anno_cog(self.geneset_id, "All", cog_path, group_id=self.specimen_group,
                                                     group_detail=self.group_detail, software_ver=soft_db)
        self.api_dic["cog"].add_anno_cog_nog(return_id, self.output_dir + "/cog")
        self.api_dic["cog"].add_anno_cog_function(return_id, self.output_dir + "/cog")
        self.api_dic["cog"].add_anno_cog_category(return_id, self.output_dir + "/cog")
        self.anno_id["cog"] = str(return_id)

    @time_count
    def export_kegg(self):
        self.api_dic["kegg"] = self.api.api("metagenomic.mg_anno_kegg")
        # xml_path = os.path.join(self.output_dir, "kegg/kegg_merge.xml")
        # xml_path = os.path.join(self.remote_dir, "kegg/kegg_merge.xml")  # 此处已将文件路径上传，导表函数只用来保存路径，不再上传kegg_merge.xml
        kegg_web_path = os.path.join(self.remote_dir, 'kegg/gene_kegg_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["kegg"]}
        # return_id = self.api_dic["kegg"].add_anno_kegg(self.geneset_id, 'All',
        #                                                kegg_web_path, xml_file=xml_path, group_id=self.specimen_group,
        #                                                group_detail=self.group_detail, software_ver=soft_db
        #                                                )
        return_id = self.api_dic["kegg"].add_anno_kegg(self.geneset_id, 'All',
                                                       kegg_web_path, group_id=self.specimen_group,
                                                       group_detail=self.group_detail, software_ver=soft_db
                                                       )
        self.api_dic["kegg"].add_anno_kegg_gene(return_id, self.output_dir + "/kegg")
        self.api_dic["kegg"].add_anno_kegg_orthology(return_id, self.output_dir + "/kegg")
        self.api_dic["kegg"].add_anno_kegg_module(return_id, self.output_dir + "/kegg")
        self.api_dic["kegg"].add_anno_kegg_enzyme(return_id, self.output_dir + "/kegg")
        # graph_path = os.path.join('/'.join(self.remote_dir.split('/')[:5]), 'metag/KEGG_Pathway', str(self._sheet.id),
        #                           str(return_id) + '/')
        # graph_path = os.path.join('metag/KEGG_Pathway', str(self._sheet.id), str(return_id) + '/')
        graph_path = os.path.join(self.remote_dir, "kegg/pathway_img.tar.gz")
        self.api_dic["kegg"].add_anno_kegg_pathway(return_id, self.output_dir + "/kegg", graph_path)
        self.anno_id["kegg"] = str(return_id)

    @time_count
    def export_cazy(self):
        self.api_dic["cazy"] = self.api.api("metagenomic.mg_anno_cazy")
        cazy_path = os.path.join(self.remote_dir, 'cazy/gene_cazy_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["cazy"]}
        return_id = self.api_dic["cazy"].add_anno_cazy(self.geneset_id, "All",
                                                       cazy_path, group_id=self.specimen_group,
                                                       group_detail=self.group_detail, software_ver=soft_db)
        self.api_dic["cazy"].add_anno_cazy_family(return_id, self.output_dir + "/cazy")
        self.api_dic["cazy"].add_anno_cazy_class(return_id, self.output_dir + "/cazy")
        self.anno_id["cazy"] = str(return_id)

    @time_count
    def export_vfdb(self):
        self.api_dic["vfdb"] = self.api.api("metagenomic.mg_anno_vfdb")
        vfdb_path = os.path.join(self.remote_dir, 'vfdb/gene_vfdb_total_anno.xls')
        vfdb_group_detail = {'all': sorted(self.spname_spid.values())}
        soft_db = {"soft_db_info": self.soft_db_info["vfdb"]}
        return_id = self.api_dic["vfdb"].add_anno_vfdb(self.geneset_id, "All",
                                                       vfdb_path, group_id="all",
                                                       group_detail=vfdb_group_detail, software_ver=soft_db)
        self.api_dic["vfdb"].add_anno_vfdb_vfs(return_id, self.output_dir + "/vfdb", "core")
        self.api_dic["vfdb"].add_anno_vfdb_vfs(return_id, self.output_dir + "/vfdb", "predict")
        self.api_dic["vfdb"].add_anno_vfdb_pie(return_id, self.output_dir + "/vfdb")
        self.anno_id["vfdb"] = str(return_id)

    @time_count
    def export_ardb(self):
        self.api_dic["ardb"] = self.api.api("metagenomic.mg_anno_ardb")
        ardb_path = os.path.join(self.remote_dir, 'ardb/gene_ardb_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["ardb"]}
        return_id = self.api_dic["ardb"].add_anno_ardb(self.geneset_id, "All",
                                                       ardb_path, group_id=self.specimen_group,
                                                       group_detail=self.group_detail, software_ver=soft_db)
        count = 0
        new_ardb_path = os.path.join(self.output_dir, 'ardb/gene_ardb_anno.xls')
        if not os.path.exists(new_ardb_path):
            self.logger.info("ardb注释结果为空")
            return
        for index, line in enumerate(open(new_ardb_path, 'r')):
            count += 1
        if count > 1:
            self.api_dic["ardb"].add_anno_ardb_arg(return_id, self.output_dir + "/ardb")
            self.api_dic["ardb"].add_anno_ardb_type(return_id, self.output_dir + "/ardb")
            self.api_dic["ardb"].add_anno_ardb_class(return_id, self.output_dir + "/ardb")
            self.anno_id["ardb"] = str(return_id)

    @time_count
    def export_card(self):
        self.api_dic["card"] = self.api.api("metagenomic.mg_anno_card")
        card_path = os.path.join(self.remote_dir, 'card/gene_card_anno.xls')
        soft_db = {"soft_db_info": self.soft_db_info["card"]}
        return_id = self.api_dic["card"].add_anno_card(self.geneset_id, "All",
                                                       card_path, group_id=self.specimen_group,
                                                       group_detail=self.group_detail, software_ver=soft_db)
        count = 0
        new_card_path = os.path.join(self.output_dir, 'card/gene_card_anno.xls')
        if not os.path.exists(new_card_path):
            self.logger.info("card注释结果为空")
            return
        for index, line in enumerate(open(new_card_path, 'r')):
            count += 1
        self.logger.info("card 注释的结果含行数：{}".format(count))
        if count > 1:
            self.api_dic["card"].add_anno_card_aro(return_id, self.output_dir + "/card")
            # self.api_dic["card"].add_anno_card_class(return_id, self.output_dir + "/card") # 新版本的数据库不再提供这张表 20201027 qingchen.zhang
            self.anno_id["card"] = str(return_id)
        else:
            self.logger.info("card注释结果为空")

    @time_count
    def export_overview(self):
        """
        注釋表信息總覽
        :return:
        """
        self.api_dic["overview"] = self.api.api("metagenomic.mg_anno_overview")
        return_id = self.api_dic["overview"].add_anno_overview(self.geneset_id)
        self.api_dic["overview"].add_anno_overview_detail(return_id, self.overview.option('gene_overview').prop['path'])
        if self.pipeline == 'pipeline2.1':
            self.api_dic["overview"].db['anno_overview'].update({'_id': return_id},
                    {'$set': {'new_sum': os.path.join(self.remote_dir, 'geneset', 'anno_overview.xls')}})

    @time_count
    def export_kraken(self):
        self.api_dic["kraken"] = self.api.api("metagenomic.common_api")
        name = "Kraken_Origin_" + datetime.datetime.now().strftime("%Y%m%d")
        self.spname_spid
        other = {"specimen": ','.join(sorted(self.spname_spid.values())),
                 "anno_file": os.path.join(self.sheet.output, "kraken/")}
        params_json = {
            "confidence": float(self.option("kk_confidence")),
            "task_id": self.sheet.id,
            "qc": self.option("qc"),
        }
        if int(self.option("qc")):
            params_json.update({
                'qc_quality': int(self.option("qc_quality")),
                "qc_length": int(self.option("qc_length"))
            })
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        soft_db = {"soft_db_info": self.soft_db_info["kraken2"]}
        return_id = self.api_dic["kraken"].add_main("anno_kraken", name=name, others=other,
                                                    params=params, software_ver=soft_db)
        key_map = {
            "Domain": 'd__' , "Kingdom": 'k__', 'Phylum': "p__", 'Class': "c__",
            'Order': "o__", 'Family': "f__", 'Genus': "g__", 'Species': "s__"
        }
        key_map.update(self.spname_spid)
        for f in os.listdir(self.taxon_outset["kraken"].output_dir):
            file_path = os.path.join(self.taxon_outset["kraken"].output_dir, f)
            self.link(file_path, "output/kraken/")
            self.api_dic["kraken"].add_detail(file_path, 'anno_kraken_detail',
                                              return_id, 'anno_id',
                                              key_map=key_map)

    @time_count
    def export_metaphlan(self):
        self.api_dic["metaphlan"] = self.api.api("metagenomic.common_api")
        name = "Metaphlan_Origin_" + datetime.datetime.now().strftime("%Y%m%d")
        self.spname_spid
        other = {"specimen": ','.join(sorted(self.spname_spid.values())),
                 "anno_file": os.path.join(self.sheet.output, "metaphlan/")}
        params_json = {
            "min_cu_len": int(self.option("mph_min_cu_len")),
            "stat": self.option("mph_stat"),
            "stat_q": float(self.option("mph_stat_q")),
            "qc": int(self.option("qc")),
        }
        if int(self.option("qc")):
            params_json.update({
                'qc_quality': int(self.option("qc_quality")),
                "qc_length": int(self.option("qc_length"))
            })
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        soft_db = {"soft_db_info": self.soft_db_info["metaphlan3"]}
        return_id = self.api_dic["metaphlan"].add_main("anno_metaphlan", name=name, others=other,
                                                    params=params, software_ver=soft_db)
        key_map = {
            "Domain": 'd__' , "Kingdom": 'k__', 'Phylum': "p__", 'Class': "c__",
            'Order': "o__", 'Family': "f__", 'Genus': "g__", 'Species': "s__"
        }
        key_map.update(self.spname_spid)
        for f in os.listdir(self.taxon_outset["metaphlan"].output_dir):
            file_path = os.path.join(self.taxon_outset["metaphlan"].output_dir, f)
            self.link(file_path, "output/metaphlan/")
            self.api_dic["metaphlan"].add_detail(file_path, 'anno_metaphlan_detail',
                                              return_id, 'anno_id',
                                              key_map=key_map)

    def get_param(self, anno):
        params = {
            "geneset_id": str(self.geneset_id),
            "group_detail": group_detail_sort(self.group_detail),
            "group_id": str(self.specimen_group),
            "method": "rpkm",
            "task_type": 2,
            "anno_id": str(self.anno_id[anno]),
            "anno_type": anno,
        }
        return params

    @time_count
    def export_composition(self, analysis, anno_list):
        self.api_dic["composition"] = self.api.api("metagenomic.composition")
        level_id = {
            "gene": "",
            "nr": 7,
            "cog": 10,
            "kegg": 12,
            "cazy": 18,
            "vfdb": 26,
            "ardb": 21,
            "card": 24,
        }
        params = {
            "graphic_type": analysis,
            "group_method": "",
        }
        analysis_dic = {
            "heatmap": "Heatmap",
            "circos": "Circos",
            "bar": "Bar",
        }
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if analysis == "heatmap":
            params["species_cluster"] = ""
            params["specimen_cluster"] = ""
            params["top"] = 50
        # for select_anno in self.anno_table.keys():
        for select_anno in anno_list:
            if select_anno == 'geneset':
                continue
            if select_anno != 'nr':
                params['second_level'] = ""
            elif 'second_level' in params.keys():
                del params['second_level']
            submit_location = "composition_" + analysis + "_" + select_anno
            #name = analysis_dic[analysis] + "_" + select_anno.upper() + "_" + self.default_level1[
            #    select_anno] + "_" + time
            name = analysis_dic[analysis] + "_" + select_anno.upper() + "_" + self.default_level1[
                select_anno] + "_Origin"  # 工作流运行的分析，导表用Origin后缀
            self.main_collection_delete('export_composition', name)  # 导表前删除重复的主表
            if analysis in ["bar", "circos"]:
                params["combine_value"] = "0.01" if select_anno in ['nr', 'gene', 'ardb', 'card'] else "0"
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            return_id = self.api_dic["composition"].add_composition(analysis, select_anno, name=name,
                                                                    params=params)
            if analysis == 'heatmap':
                detail_file = os.path.join(self.output_dir, "composition", analysis, select_anno, "taxa.table.xls")
            else:
                detail_file = os.path.join(self.output_dir, "composition", analysis, select_anno,
                                           "taxa.percents.table.xls")
            self.api_dic["composition"].add_composition_detail(detail_file, return_id, species_tree="",
                                                               specimen_tree="", group_method="")  # 导表函数增加了group_method参数

    @time_count
    def export_beta(self, analysis, anno_list):
        self.api_dic["beta"] = self.api.api("metagenomic.beta_diversity")
        level_id = {
            "gene": "",
            "nr": 7,
            "cog": 11,
            "kegg": 14,
            "cazy": 19,
            "vfdb": 28,
            "ardb": 23,
            "card": 25,
        }
        params = {}
        analysis_type = analysis.lower()
        analysis_dic = {
            "pca": "PCA",
            "pcoa": "PCoA",
            "nmds": "NMDS",
            "rda": "RDACCA",
            "dbrda": "dbRDA",
        }
        params["analysis_type"] = "rda_cca" if analysis_type == "rda" else analysis_type
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if analysis_type in ["nmds", "pcoa", "dbrda"]:
            params["distance_method"] = "bray_curtis"
        if analysis_type in ["pca", "rda", "dbrda"] and self.option('envtable').is_set:
            params["env_id"] = str(self.env_id)  #
            params["env_labs"] = self.env_labs  #
        # for select_anno in self.anno_table.keys():
        for select_anno in anno_list:
            self.logger.info("anno is " + select_anno)
            if select_anno == 'geneset':
                continue
            if analysis_type == 'rda':
                anno_type = 'rdacca'  # var anno_type added by zouxuan
            else:
                anno_type = analysis_type
            submit_location = 'betadiversity_' + anno_type + '_' + select_anno
            # name = analysis_dic[analysis_type] + "_" + select_anno.upper() + "_" + self.default_level2[
            #     select_anno] + "_" + time
            name = analysis_dic[analysis_type] + "_" + select_anno.upper() + "_" + self.default_level2[
                select_anno] + "_Origin"  # 工作流运行的分析，导表用Origin后缀
            self.main_collection_delete('export_beta', name)  # 导表前删除重复的主表
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            self.logger.info(params)
            if analysis_type in ["nmds", "pcoa"]:
                path = os.path.join(self.output_dir, "compare", analysis, select_anno)
                self.api_dic["beta"].add_beta_diversity(path, analysis_type, main="True", anno_type=select_anno,
                                                        name=name, params=params)
            elif analysis_type == 'pca':
                path = os.path.join(self.output_dir, "compare", analysis, select_anno)
                pca_path = os.path.join(self.remote_dir, 'compare/Pca/', select_anno)
                if self.env_id != '':
                    self.api_dic["beta"].add_beta_diversity(path, analysis_type, main="True", web_path=pca_path,
                                                            anno_type=select_anno,
                                                            env_id=self.env_id, name=name, params=params)
                else:
                    self.api_dic["beta"].add_beta_diversity(path, analysis_type, main="True", web_path=pca_path,
                                                            anno_type=select_anno,
                                                            name=name, params=params)
            elif analysis_type == "rda":
                path = os.path.join(self.output_dir, "correlation", analysis, select_anno)
                if not os.path.exists(path):
                    self.logger.info("不存在路径%s，跳过%s的rda分析" % (path, select_anno))
                    return
                self.api_dic["beta"].add_beta_diversity(path, "rda_cca", main="True", anno_type=select_anno,
                                                        env_id=self.env_id, name=name, params=params)
            elif analysis_type == "dbrda":
                path = os.path.join(self.output_dir, "correlation", analysis, select_anno)
                if not os.path.exists(path):
                    self.logger.info("不存在路径%s，跳过%s的dbrda分析" % (path, select_anno))
                    return
                self.api_dic["beta"].add_beta_diversity(path, "dbrda", main="True", anno_type=select_anno,
                                                        env_id=self.env_id, name=name, params=params)

    @time_count
    def export_hclust(self, anno_list):
        self.api_dic["hcluster"] = self.api.api("metagenomic.hcluster_tree")
        self.api_dic["dist"] = self.api.api("metagenomic.distance_metagenomic")
        level_id = {
            "gene": "",
            "nr": 7,
            "cog": 11,
            "kegg": 14,
            "cazy": 19,
            "vfdb": 28,
            "ardb": 23,
            "card": 25,
        }
        params = {
            "hcluster_method": "average",
            "distance_method": "bray_curtis",
        }
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        # for select_anno in self.anno_table.keys():
        for select_anno in anno_list:
            if select_anno == 'geneset':
                continue
            else:
                submit_location = 'hclustertree_' + select_anno
                # name = "Hcluster_" + select_anno.upper() + "_" + self.default_level2[
                #     select_anno] + "_" + time
                name = "Hcluster_" + select_anno.upper() + "_" + self.default_level2[
                    select_anno] + "_Origin"  # 工作运行的分析，导表用Origin后缀
            self.main_collection_delete('export_hclust', name)  # 导表前删除重复的主表
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            matrix_path = os.path.join(self.output_dir, "compare/Distance", select_anno,
                                       "bray_curtis_new_abund_table.xls")
            newick_path = os.path.join(self.output_dir, "compare/Hcluster", select_anno, "hcluster.tre")
            dist_id = self.api_dic["dist"].add_dist_table(matrix_path, main=True, name=name, params=params)
            self.api_dic["hcluster"].add_hcluster_tree(newick_path, main=True, anno_type=select_anno,
                                                       name=name, params=params, update_dist_id=dist_id)

    @time_count
    def export_cor_heatmap(self, anno_list):
        self.api_dic["cor_heatmap"] = self.api.api("metagenomic.heatmap_cor")
        level_id = {
            "gene": "",
            "nr": 7,
            "cog": 11,
            "kegg": 14,
            "cazy": 19,
            "vfdb": 28,
            "ardb": 23,
            "card": 25,
        }
        params = {
            "top": 50,
            "env_id": str(self.env_id),
            "env_labs": self.env_labs,
            "coefficient": "spearmanr",
            "species_cluster": "",
            "env_cluster": "",
        }
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        # for select_anno in self.anno_table.keys():
        for select_anno in anno_list:
            if select_anno == 'geneset':
                continue
            submit_location = "heatmapcor_" + select_anno
            # name = "CorrHeatmapSpearman_" + select_anno.upper() + "_" + self.default_level2[
            #     select_anno] + "_" + time
            name = "CorrHeatmapSpearman_" + select_anno.upper() + "_" + self.default_level2[
                select_anno] + "_Origin"  # 工作流运行的分析，导表用Origin后缀
            self.main_collection_delete('export_cor_heatmap', name)  # 导表前删除重复的主表
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            return_id = self.api_dic["cor_heatmap"].add_heatmap_cor(select_anno, name=name, params=params)
            corr_path = os.path.join(self.output_dir, "correlation/cor_heatmap", select_anno,
                                     "correlation.xls")
            pvalue_path = os.path.join(self.output_dir, "correlation/cor_heatmap", select_anno, "pvalue.xls")
            self.api_dic["cor_heatmap"].add_heatmap_cor_detail(corr_path, "correlation", return_id, species_tree="",
                                                               env_tree="")
            self.api_dic["cor_heatmap"].add_heatmap_cor_detail(pvalue_path, "pvalue", return_id,
                                                               species_tree="", env_tree="")

    @time_count
    def export_personal_anno(self):
        ###################### 测试用
        #self.geneset_id = ObjectId("5bf4f7f8a4e1af4253f48c45")
        common_param = {}
        self.personal_api = {}
        common_param["geneset_id"] = self.geneset_id
        common_param["group_detail"] = self.group_detail
        common_param["anno_dir"] = os.path.join(self.remote_dir, "personal_anno")
        self.personal_api_anno = self.api.api("metagenomic.personal_anno")
        self.personal_api["go"] = self.api.api("metagenomic.anno_go")
        self.personal_api["qs"] = self.api.api("metagenomic.qs_anno")
        self.personal_api["probio"] = self.api.api("metagenomic.probiotics")
        self.personal_api["pfam"] = self.api.api("metagenomic.mg_anno_pfam")
        self.personal_api["p450"] = self.api.api("metagenomic.mg_anno_cyps")
        self.personal_api["tcdb"] = self.api.api("metagenomic.tcdb")
        self.personal_api["mvirdb"] = self.api.api("metagenomic.mvirdb")
        self.personal_api["sec"] = self.api.api("metagenomic.mg_anno_sec")
        self.personal_api["t3ss"]  = self.api.api("metagenomic.mg_anno_ttss")
        self.personal_api["phi"] = self.api.api("metagenomic.mg_anno_phi")
        self.personal_params = {
            "geneset_id": str(self.geneset_id),
            "group_detail": group_detail_sort(self.group_detail),
            "group_id": str(self.specimen_group),
            "database": "",
            "submit_location": "",
            "task_type": 2
        }
        database = self.option("personal_anno")
        #############################测试用
        #database = "qs;probio;p450;pfam;mvirdb;tcdb"
        nr_method = "best_hit,de_unclassied,lca"
        self.personal_api_anno.add_main(database, self.personal_api, self.personal_params, common_param,nr_method=nr_method)
        self.logger.info("personal main table finish")

    def export_software(self):
        """
        软件列表导表
        :return:
        """
        software = self.api.api('metagenomic.software')
        task_id = "_".join(self._sheet.id.split("_"))
        new_main_id = software.add_software(task_id=task_id, params="")
        software.add_software_detail(new_main_id)
        ## 更新组装软件名称
        if 'IDBA_UD' in self.option("assemble_type"):
            software_config = {"name" : "IDBA_UD",
                                "source" : "http://www.cs.hku.hk/~alse/idba_ud",
                                "en_item" : "Assembly",
                                "item" : "拼接",
                               "soft_id": ObjectId(new_main_id),
                               "sort": 3,
                                "version" : "Version 1.1.1"}
        elif 'Megahit' in self.option('assemble_type'):
            software_config = {"name" : "Megahit",
                                "source" : "http://www.l3-bioinfo.com/products/megahit.htm",
                                "en_item" : "Assembly",
                                "item" : "拼接",
                               "soft_id": ObjectId(new_main_id),
                               "sort": 3,
                                "version" : "Version 1.1.2"}
        elif 'SOAPdenovo2' in self.option('assemble_type'):
            software_config = {"name" : "SOAPdenovo2",
                                "source" : "https://github.com/aquaskyline/SOAPdenovo2",
                                "en_item" : "Assembly",
                                "item" : "拼接",
                               "soft_id": ObjectId(new_main_id),
                               "sort": 3,
                                "version" : "Version 2.04"}
        software.update_mongo('soft_detail',{"soft_id":ObjectId(new_main_id), "name": "Megahit"}, software_config)
        self.logger.info("software finish")

    def main_collection_delete(self, func_name, name=None):
        """
        重导表检查，去除重复主表
        新增一个name参数，用于兼容删除重复'export_composition', 'export_beta', 'export_hclust', 'export_cor_heatmap'四个主表时误删的的问题
        """
        rm_repeat = self.api.api('metagenomic.repeat')
        collect_relation = {
            "export_specimen_group": "specimen_group",
            "export_qc": "data_stat",
            "export_assem": "assemble_stat",
            "export_predict": "predict_gene",
            "export_geneset": "geneset",
            "export_nr": "anno_nr",
            "export_nr_lca": "anno_nr",
            "export_nr_deunclass": "anno_nr",
            "export_cog": "anno_cog",
            "export_kegg": "anno_kegg",
            "export_cazy": "anno_cazy",
            "export_vfdb": "anno_vfdb",
            "export_ardb": "anno_ardb",
            "export_card": "anno_card",
            "export_overview": "anno_overview",
            "export_composition": "composition",
            "export_beta": "beta_diversity",
            "export_hclust": "hcluster_tree",
            "export_cor_heatmap": "heatmap_cor",
            "export_personal_anno": "anno_cyps,anno_go,anno_mvir,anno_pfam,anno_phi,anno_probio,anno_qs,anno_sec,anno_tcdb,anno_ttss",
            "export_kraken": "anno_kraken",
            "export_metaphlan": "anno_metaphlan"
        }
        # rm_each_list = ["export_qc","export_assem", "export_predict","export_composition","export_beta","export_hclust","export_specimen_group"]
        rm_each_list = ["export_qc","export_assem", "export_predict","export_composition","export_beta","export_hclust","export_specimen_group"]
        ana_list = ['export_composition', 'export_beta', 'export_hclust', 'export_cor_heatmap']
        collect_name = collect_relation[func_name]
        if func_name in ana_list:
            if not name:
                return
            else:
                rm_repeat.rm_main(collect_name, name=name)
        elif func_name == "export_nr":
            rm_repeat.rm_main(collect_name, name="NR_Origin")
        elif func_name == "export_nr_lca":
            rm_repeat.rm_main(collect_name, name="NR_Origin_LCA")
        elif func_name == "export_nr_deunclass":
            rm_repeat.rm_main(collect_name, name="NR_Origin_deunclassied")
        elif  func_name in rm_each_list:
            rm_repeat.rm_each(collect_name)
        elif func_name == "export_personal_anno":
            person_main_list = collect_name.split(",")
            for each in person_main_list:
                rm_repeat.rm_each(each)
        else:
            #params_name = func_name.replace("export_","")
            #params = self.dic_params[params_name]
            rm_repeat.rm_main(collect_name)
            #print "repeat main_collection {} delete".format(collect_name)
            #self.logger.info("repeat main_collection {} delete".format(collect_name))

    def run_from_task(self):
        '''
        pipeline2.1 从老的task_id获取数据用于后续流程
        因为后续的模块都依赖于 self.gene_set 模块，所以需要手动设置此模块的监听，start和end事件
        同时配置后续模块依赖的变量
        '''
        # 手动开启监听并启动start事件
        self.gene_set.start_listener()
        self.gene_set.fire('start')
        self.step.copy_data.start()
        self.step.update()
        # 从老项目复制mongo表数据
        copy_data = CopyDemo(self.option('old_task_id'), self._sheet.id, self._sheet.project_sn,
                             self._sheet.member_id, self._sheet.member_type,
                             old_db_version=self.option('old_task_db'))
        if not copy_data.from_task():
            self.set_error('从老的task_id: {}获取mongo数据失败'.format(self.option('old_task_id')))

        # 下载注释需要用的基因组蛋白序列文件和reads_num文件
        try:
            # 获取复制后的分组信息
            metag_api = Metagenomic()
            metag_api._config = self.config
            group_info = metag_api.common_find_one('specimen_group',
                                                   {'task_id': self._sheet.id,
                                                    'group_name': self.group_name})
            if group_info:
                self.specimen_group = str(group_info['_id'])
                self.spname_spid = name2id(self._sheet.id, 'task')
                category, specimen = group_info['category_names'], group_info['specimen_names']
                for i in range(len(category)):
                    for sp in specimen[i]:
                        if category[i] not in self.group_detail:
                            self.group_detail[category[i]] = []
                        self.group_detail[category[i]].append(sp)
            else:
                data_stat = metag_api.common_find_one('data_stat',
                                                      {'task_id': self._sheet.id,
                                                       'type': 'raw'})
                data_stat_id = str(data_stat['_id'])
                self.spname_spid = name2id(data_stat_id, type="raw")
                self.specimen_group = 'all'
                self.group_detail = {'all': sorted(self.spname_spid.values())}

            # 获取老项目基因集信息，获取相应的文件地址并下载
            if self.option('old_task_db') == 0:
                db = get_old_mongo()[1]
            else:
                db = get_mongo()[1]
            geneset_info = db['geneset'].find_one({'task_id': self.option('old_task_id'),
                                                   'type': 1})
            reads_num = geneset_info['reads_num']
            uni_faa = '/'.join(reads_num.split('/')[:-2]) + '/uniGeneset/gene.uniGeneset.faa'
            from biocluster.api.file.remote import RemoteFileManager
            reads_num = RemoteFileManager(reads_num)
            uni_faa = RemoteFileManager(uni_faa)
            rpkm = RemoteFileManager(geneset_info['rpkm'])
            gene_list_length = RemoteFileManager(geneset_info['gene_list_length'])
            reads_num.download(os.path.join(self.work_dir, 'remote_input/'))
            uni_faa.download(os.path.join(self.work_dir, 'remote_input/'))
            rpkm.download(os.path.join(self.work_dir, 'remote_input/'))
            gene_list_length.download(os.path.join(self.work_dir, 'remote_input/'))

            self.gene_set.option('reads_abundance', reads_num.local_path)
            self.gene_set.option('uni_fastaa', uni_faa.local_path)
            self.gene_set.option('gene_length_list', gene_list_length.local_path)
            self.anno_table['geneset'] = rpkm.local_path  # 后续依赖的变量
            # 上传老项目中用于此次分析中的文件
            if not os.path.exists(os.path.join(self.output_dir, 'geneset')):
                os.mkdir(os.path.join(self.output_dir, 'geneset'))
                os.mkdir(os.path.join(self.output_dir, 'geneset/uniGeneset'))
                os.mkdir(os.path.join(self.output_dir, 'geneset/gene_profile'))
            self.link(uni_faa.local_path, 'output/geneset/uniGeneset/')
            self.link(reads_num.local_path, 'output/geneset/gene_profile/')
            self.link(gene_list_length.local_path, 'output/geneset/gene_profile/')
            # 更新geneset信息
            gene_set_update = {
                'reads_num': os.path.join(self.remote_dir,
                                          'geneset/gene_profile/' + os.path.basename(reads_num.local_path)),
                'gene_list_length': os.path.join(self.remote_dir,
                                                 'geneset/gene_profile/' + os.path.basename(gene_list_length.local_path)),
            }
            gene_set_col = metag_api.db['geneset']
            gene_set_col.update_one({'task_id': self._sheet.id, 'name': 'GENESET_Origin'},
                                    {'$set': gene_set_update})
            new_geneset_info = metag_api.find_origin_geneset_info(self._sheet.id)
            self.geneset_id = new_geneset_info['_id']  # 后续依赖的变量
        except Exception as e:
            self.set_error('下载老项目task_id {}文件时出错: {}'.format(self.option('old_task_id'), e))
        # 手动出发结束事件
        self.gene_set.set_end()
        self.gene_set.fire('end')
        self.step.copy_data.finish()
        self.step.update()

    def run(self):
        """
        运行 meta_genomic workflow
        :return:
        """
        task_info = self.api.api('task_info.mg_task_info')
        task_info.add_task_info()
        if self.pipeline == "pipeline3":
            self.sequence.on('end', self.end)
            self.run_sequence()
        else:
            if self.pipeline == "pipeline1":
                if self.option('rm_host'):
                    self.sequence.on('end', self.run_rm_host)
                    self.rm_host.on('end', self.run_assem)
                else:
                    self.sequence.on('end', self.run_assem)
                self.add_event("read_taxon_anno")
                if self.option("kraken"):
                    self.sequence.on("end", self.run_kraken)
                if self.option("metaphlan"):
                    self.sequence.on("end", self.run_metaphlan)
                self.assem_soapdenovo.on('end', self.run_gene_predict)
                self.assem_idba.on('end', self.run_gene_predict)
                self.gene_predict.on('end', self.run_gene_set)
            elif self.pipeline == "pipeline2":
                self.sequence_deal.on('end', self.run_gene_set)
            if 'nr' in self.choose_anno:
                self.gene_set.on('end', self.run_nr)
                self.anno_tool.append(self.nr)
            if 'kegg' in self.choose_anno:
                self.gene_set.on('end', self.run_kegg)
                self.anno_tool.append(self.kegg)
            if 'cog' in self.choose_anno:
                self.gene_set.on('end', self.run_cog)
                self.anno_tool.append(self.cog)
            if 'cazy' in self.choose_anno:
                self.gene_set.on('end', self.run_cazy)
                self.all_anno.append(self.cazy)
            if 'ardb' in self.choose_anno:
                self.gene_set.on('end', self.run_ardb)
                self.all_anno.append(self.ardb)
            if 'card' in self.choose_anno:
                self.gene_set.on('end', self.run_card)
                self.all_anno.append(self.card)
            if 'vfdb' in self.choose_anno:
                self.gene_set.on('end', self.run_vfdb)
                self.all_anno.append(self.vfdb)
            if len(self.anno_tool) != 0:
                #self.on_rely(self.anno_tool, self.run_anno)
                self.on_rely(self.anno_tool, self.run_anno_new)
                self.all_anno.append(self.anno)
                self.all_anno.append(self.anno_nr_lca)
                self.all_anno.append(self.anno_nr_deunclass)
            self.sample_in_group = self.get_sample()
            if len(self.all_anno) == 0:
                # self.gene_set.on('end', self.end)
                if self.option("personal_anno") != "none":
                    self.gene_set.on('end', self.version2_add)
                else:
                    self.gene_set.on('end', self.end)
            else:
                self.logger.info("test>>>>>>>>>>>>>>>>>>>>>>>>>")
                self.logger.info(self.all_anno)
                self.on_rely(self.all_anno, self.run_overview)
                if len(self.sample_in_group) < 2:
                    if self.option("personal_anno") != "none":
                        self.overview.on('end', self.version2_add)
                    else:
                        self.overview.on('end', self.end)
                elif len(self.sample_in_group) == 2:
                    self.overview.on('end', self.run_analysis, 'composition')
                elif len(self.sample_in_group) > 2:
                    self.overview.on('end', self.run_analysis, 'all')
            if self.option('test'):
                pass
            elif self.pipeline == "pipeline1":
                self.run_sequence()
            elif self.pipeline == "pipeline2":
                self.logger.info("start pipeline2")
                self.run_sequence_deal()
            elif self.pipeline == 'pipeline2.1':
                self.run_from_task()
        super(MetaGenomicWorkflow, self).run()

