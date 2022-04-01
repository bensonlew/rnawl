# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

"""宏基因组分析工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from mbio.packages.metagenomic.id_convert import name2id
from bson import ObjectId
import os
import json
import shutil
import time
import datetime
import functools


def time_count(func):  # 统计导表时间
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run %s at %s' % (func_name, start_time))
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，该阶段导表已进行{}s".format(func_name, end - start))

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
            # {'name': 'main_id', 'type': 'string'},  # 原始序列主表_id
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            # {'name': 'fq_type', 'type': 'string', 'default': 'PE'},  # PE OR SE
            {'name': 'speciman_info', 'type': 'infile', 'format': 'sequence.profile_table'},  # 样本集信息表
            {'name': 'raw_info', 'type': 'infile', 'format': 'sequence.profile_table'},  # 原始序列的信息表
            {'name': 'qc_info', 'type': 'infile', 'format': 'sequence.profile_table'},  # 质控后的信息表
            {'name': 'insertsize', 'type': 'infile', 'format': 'sample.insertsize_table'},  # 插入片段长度表
            {'name': 'qc', 'type': 'bool', 'default': True},  # 是否需要质控
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
            # 选择拼接策略，SOAPdenovo2 OR IDBA_UD OR Megahit OR Multiple
            {'name': 'min_contig', 'type': 'int', 'default': 300},  # 拼接序列最短长度
            {'name': 'min_gene', 'type': 'int', 'default': 100},  # 预测基因最短长度
            {'name': 'cdhit_identity', 'type': 'float', 'default': 0.95},  # 基因序列聚类相似度
            {'name': 'cdhit_coverage', 'type': 'float', 'default': 0.9},  # 基因序列聚类覆盖度
            {'name': 'soap_identity', 'type': 'float', 'default': 0.95},  # 基因丰度计算相似度
            {'name': 'nr', 'type': 'bool', 'default': True},  # 是否进行nr注释
            {'name': 'cog', 'type': 'bool', 'default': True},  # 是否进行cog注释
            {'name': 'kegg', 'type': 'bool', 'default': True},  # 是否进行kegg注释
            {'name': 'cazy', 'type': 'bool', 'default': True},  # 是否进行cazy注释
            {'name': 'ardb', 'type': 'bool', 'default': True},  # 是否进行ardb注释
            {'name': 'card', 'type': 'bool', 'default': True},  # 是否进行card注释
            {'name': 'vfdb', 'type': 'bool', 'default': True},  # 是否进行vfdb注释
            {'name': 'anno_nr', 'type': 'string', 'default': 'nrnr'},  # 适应页面参数
            {'name': 'anno_list' , 'type': 'string', 'default': 'function'},  # 功能注释列表，逗号分割，适应页面参数
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},  # 物种/功能分析输入环境因子表
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 物种/功能分析输入group表
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        '''获取数据库信息'''
        # self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.json"  # 尚未指定具体的数据库
        # self.json_dict = self.get_json()
        '''初始化module/tool'''
        self.sequence = self.add_module('sequence.meta_genomic')
        self.qc = self.add_module('meta.qc.qc_and_stat')
        self.rm_host = self.add_module('meta.qc.bwa_remove_host')
        self.assem_soapdenovo = self.add_module('assemble.mg_ass_soapdenovo')
        self.assem_idba = self.add_module('assemble.mg_ass_idba')
        self.gene_predict = self.add_module('gene_structure.gene_predict')
        self.gene_set = self.add_module('cluster.uni_gene')
        self.nr = self.add_module('align.meta_diamond')
        self.cog = self.add_module('align.meta_diamond')
        self.kegg = self.add_module('align.meta_diamond')
        self.anno = self.add_module('annotation.mg_common_anno_stat')
        self.cazy = self.add_module('annotation.cazy_annotation')
        self.ardb = self.add_module('annotation.ardb_annotation')
        self.card = self.add_module('annotation.card_annotation')
        self.vfdb = self.add_module('annotation.vfdb_annotation')
        self.table = dict()
        self.table['check'] = self.add_tool('meta.create_abund_table')
        self.composition = self.add_module('meta.composition.composition_analysis')
        self.compare = self.add_module('meta.beta_diversity.beta_diversity')
        self.correlation = self.add_tool('statistical.pearsons_correlation')
        # self.XXX = self.add_module("XXX")
        # self.XXX = self.add_tool("XXX")
        '''add_steps'''
        self.step.add_steps('sequence', 'qc_', 'rm_host', 'assem', 'gene_predict', 'gene_set', 'nr_', 'cog',
                            'kegg', 'anno', 'cazy', 'vfdb', 'ardb', 'card', 'table', 'composition', 'compare',
                            'correlation')
        '''初始化自定义变量'''
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.api_dic = {}  # 存放導表函數
        self.anno_tool = []  # nr/kegg/cog注释记录
        self.all_anno = []  # 全部的注释记录(用于依赖关系)
        self.choose_anno = []  # 全部注释记录(字符型，用于物种与功能分析, 不含geneset)
        self.new_table = []  # 构建新丰度表模块(module)
        self.table_pointer = 0 # 查看对应的丰度表运行任务输出路径
        self.analysis = []  # 分析模块具体分析内容(module/tool)
        self.nr_dir = ''  # nr注释结果文件路径，导表时用
        self.cog_dir = ''
        self.kegg_dir = ''
        self.anno_table = dict()  # 注释结果表(含所有注释水平，含丰度结果表)
        self.profile_table1 = dict()  # 注释丰度表(用于组成分析，相关性heatmap图)
        self.profile_table2 = dict()  # 注释丰度表(用于样品比较分析、rda、cca、db_rda分析)
        self.default_level1 = {
            'nr': 'Genus',
            'cog': 'Function',
            'kegg': 'Level1',
            'cazy': 'Class',
            'vfdb': 'Level1',
            'ardb': 'Type',
            'card': 'Class',
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
        self.composition_dir2anno = {}  # 输出结果和导表时根据此值判断数据库类型
        self.compare_dir2anno = {}
        self.correlation_dir2anno = {}
        self.anno2correlation_tree = {}
        self.sample_in_group = []  # 根据group文件获取样本名
        self.specimen_group = "" #  group_id
        self.group_detail = {}  # 分组对应样品id{group1: [id1,id2,id3], group2: [id4,id5,id6]}
        self.env_id = ""  # 环境因子的id
        self.env_labs = ""  # 环境因子，","分割
        self.geneset_id = ""  # geneset主表_id
        self.spname_spid = {}  # 样本对应原始数据表_id
        self.anno_id = {}  # 注释表对应mongo数据库主表_id
        self.remote_dir = '/mnt/ilustre/' + self._sheet.output.replace(':','-data/') # /mnt/../m_188/188_5a1e1bc312416/tsg_26035/workflow_results
        # self.sanger_path = self.get_sanger_path()
        if self.option('test'):
            self.gene_set.option('uni_fasta').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/geneset/uniGeneset/gene.uniGeneset.fa")
            self.gene_set.option('uni_fastaa').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/geneset/uniGeneset/gene.uniGeneset.faa")
            self.gene_set.option('reads_abundance').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/geneset/gene_profile/reads_number.xls")

            self.anno_table = {
                'geneset': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/geneset/gene_profile/RPKM.xls',
                'ardb': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/ardb/gene_ardb_anno.xls',
                'card': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/card/gene_card_anno.xls',
                'vfdb': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/vfdb/gene_vfdb_predict_anno.xls',
                'cazy': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/cazy/gene_cazy_anno.xls',
                'nr': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/nr/gene_nr_anno.xls',
                'kegg': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/kegg/gene_kegg_anno.xls',
                'cog': '/mnt/ilustre/users/sanger-dev/workspace/20171213/MetaGenomic_tsg_26265/output/cog/gene_cog_anno.xls',
            }
            # self.gene_predict.option('out').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171205/MetaGenomic_tsg_26030/GenePredict/MetageneStat/output/Total.metagene.fa")
            # self.qc.option('before_qc_stat').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171205/MetaGenomic_tsg_26030/output/qc/qc_stat/new_reads.rawData.stat")
            # self.qc.option('after_qc_dir').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171205/MetaGenomic_tsg_26030/output/qc/after_qc_dir")
            '''
            self.gene_set.option('uni_fasta').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/geneset/uniGeneset/gene.uniGeneset.fa")
            self.gene_set.option('uni_fastaa').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/geneset/uniGeneset/gene.uniGeneset.faa")
            self.gene_set.option('reads_abundance').set_path("/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/geneset/gene_profile/reads_number.xls")

            self.anno_table = {
                'geneset': '/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/geneset/gene_profile/RPKM.xls',
                'ardb': '/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/ardb/gene_ardb_anno.xls',
                'card': '/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/card/gene_card_anno.xls',
                'vfdb': '/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/vfdb/gene_vfdb_predict_anno.xls',
                'cazy': '/mnt/ilustre/users/sanger-dev/workspace/20171023/MetaGenomic_metagenome/output/cazy/anno_result/gene_cazy_anno.xls',
                'nr': '/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/output/nr/gene_nr_anno.xls',
                'kegg': '/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/output/kegg/gene_kegg_anno.xls',
                'cog': '/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/output/cog/gene_cog_anno.xls',
            }
            '''

            # self.qc_fastq = self.qc.option('in_fastq')  # 暂未加入质控步骤，输入质控序列

    def check_options(self):
        """
        检查参数设置
        """
        # if not self.option('main_id') and self.option('test') == False:
        #      raise OptionError('缺少主表id')
        # self.logger.info("anno_nr is %s" % self.option('anno_nr'))
        # self.logger.info("anno_list is %s" % self.option('anno_list'))
        # self.logger.info("ref_undefined_name is %s" % self.option('ref_undefined_name'))
        self.logger.info("check_sheet_data...")
        if not self._sheet.id:
            raise OptionError('需要输入task_id')
        if not self._sheet.member_type:
            raise OptionError("需要输入member_type")
        if not self._sheet.cmd_id:
            raise OptionError("需要输入cmd_id")
        self.logger.info("check_option...")
        if not self.option('in_fastq'):
            raise OptionError('需要输入原始fastq序列')
        # if not self.option('fq_type') in ['PE', 'SE']:
        #    raise OptionError('fq序列应为PE或SE')
        if self.option('qc') and not self.option('speciman_info').is_set:
            raise OptionError('质控需提供样本集信息表')
        if not self.option('qc') and not self.option('raw_info'):
            raise OptionError('需进行质控，或者输入原始数据统计表')
        if not self.option('qc') and not self.option('qc_info'):
            raise OptionError('需进行质控，或者输入质控后数据统计表')
        # if not self.option('insertsize'):
        #     raise OptionError('需要输入insertsize表')
        if not self.option('qc_quality') > 0 and not self.option('qc_quality') < 42:
            raise OptionError('qc最小质量值超出范围，应在0~42之间')
        if not self.option('qc_length') > 0:
            raise OptionError('qc最小长度值超出范围，应大于0')
        if self.option('rm_host'):
            if self.option('ref_database') == '' and not self.option('ref_undefined').is_set:
                raise OptionError('已选择去宿主，需输入参考数据库或参考序列')
            if self.option('ref_database') != '' and self.option('ref_undefined').is_set:
                raise OptionError('去宿主不可同时提供参考数据库及参考序列')
        # if not self.option('assemble_tool') in ['soapdenovo', 'idba']:
        #     raise OptionError('请检查拼接工具是否输入正确')
        # if not self.option('assemble_type') in ['soapdenovo', 'idba', 'megahit', 'multiple']:
        if not self.option('assemble_type') in ['SOAPdenovo2', 'IDBA_UD', 'Megahit', 'Multiple']:
            raise OptionError('拼接策略参数错误，应为SOAPdenovo2/IDBA_UD/Megahit/Multiple')
        # if self.option('assemble_tool') == 'soapdenovo' and self.option('assemble_type') == 'multiple':
        #    raise OptionError('不支持SOAPdenovo混拼流程')
        if self.option('min_contig') < 200 or self.option('min_contig') > 1000:
            raise OptionError('最小Contig长度参数超出范围200~1000')
        self.logger.info("check_min_contig...")
        if self.option('min_gene') < 0:
            raise OptionError('最小基因长度参数不能为负')
        if not 0.75 <= self.option("cdhit_identity") <= 1:
            raise OptionError("cdhit identity必须在0.75，1之间")
        if not 0 <= self.option("cdhit_coverage") <= 1:
            raise OptionError("cdhit coverage必须在0,1之间")
        if not 0 < self.option("soap_identity") < 1:
            raise OptionError("soap identity必须在0，1之间")
        if not self.option('group').is_set:
            raise OptionError('必须输入分组文件')
        self.logger.info("check_option over...")
        return True

    # def get_json(self):
    #     f = open(self.json_path, 'r')
    #     json_dict = json.loads(f.read())
    #     return json_dict

    def get_sample(self):
        samp_list = []
        with open(self.option('group').prop['path'], 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                t = line.split('\t')
                samp_list.append(t)
        return samp_list

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_run(self, opts, module, event, step, start=True):
        module.set_options(opts)
        module.on('start', self.set_step, {'start': step})
        module.on('end', self.set_step, {'end': step})
        module.on('end', self.set_output, event)
        if start:
            module.run()

    def run_sequence(self):
        opts = {
            'fastq_dir': self.option('in_fastq'),
        }
        self.set_run(opts, self.sequence, 'sequence', self.step.sequence)

    def run_qc(self):
        opts = {
            'fastq_dir': self.sequence.output_dir + '/data',
            'stat_dir': self.sequence.output_dir + '/base_info',
            'insert_size': self.option('speciman_info'),
        }
        self.set_run(opts, self.qc, 'qc', self.step.qc_)

    def run_rm_host(self):
        opts = {
            'fq_type': 'PSE',
            'ref_database': self.option('ref_database'),
            'ref_undefined': self.option('ref_undefined'),
        }
        if self.option('qc'):
            opts['fastq_dir'] = self.qc.option('after_qc_dir')
        else:
            opts['fastq_dir'] = self.option('in_fastq')
        self.set_run(opts, self.rm_host, 'rm_host', self.step.rm_host)

    def run_assem(self):
        if self.option('qc'):
            opts = {
                'qc_stat': self.qc.option('after_qc_stat'),
                'raw_stat': self.qc.option('before_qc_stat'),
                'QC_dir': self.qc.option('after_qc_dir')
            }
            self.logger.info('have qc ,raw_stat is :')
            self.logger.info(self.qc.option('before_qc_stat'))
            self.logger.info(self.qc.option('before_qc_stat').prop.keys())
            self.logger.info(self.qc.option('before_qc_stat').prop['path'])
        else:
            opts = {
                'qc_stat': self.option('qc_info'),
                'raw_stat': self.option('raw_info'),
            }
            self.logger.info('havenot qc ')
            if self.option('rm_host'):
                opts['QC_dir'] = self.rm_host.option('result_fq_dir')
            else:
                opts['QC_dir'] = self.option('in_fastq')
        opts['min_contig'] = self.option('min_contig')
        self.logger.info('begin set_run ...')
        if self.option('assemble_type') == 'SOAPdenovo2':
            self.set_run(opts, self.assem_soapdenovo, 'assem', self.step.assem)
        else:
            if self.option('assemble_type') == 'IDBA_UD':
                opts['assemble_tool'] = 'idba'
                opts['method'] = 'simple'
            if self.option('assemble_type') == 'Megahit':
                opts['assemble_tool'] = 'megahit'
                opts['method'] = 'simple'
            if self.option('assemble_type') == 'Multiple':
                opts['method'] = 'multiple'
            self.logger.info('set_run idba ...')
            self.set_run(opts, self.assem_idba, 'assem', self.step.assem)
        self.logger.info('set_run over ...')

    def run_gene_predict(self):
        opts = {
            'min_gene': str(self.option('min_gene')),
        }
        if self.option('assemble_type') == 'SOAPdenovo2':
            opts['input_fasta'] = self.assem_soapdenovo.option('contig').prop['path']
        else:
            opts['input_fasta'] = self.assem_idba.option('contig').prop['path']
        self.logger.info("metagene input_fasta")
        self.logger.info(opts['input_fasta'])
        self.set_run(opts, self.gene_predict, 'gene_predict', self.step.gene_predict)

    def run_gene_set(self):
        opts = {
            'gene_tmp_fa': self.gene_predict.option('out'),
            'insertsize': self.option('insertsize'),
            'cdhit_identity': self.option('cdhit_identity'),
            'cdhit_coverage': self.option('cdhit_coverage'),
            'soap_identity': self.option('soap_identity'),
        }
        if self.option('insertsize').is_set:
            opts['insertsize'] = self.option('insertsize')
        elif self.option('qc'):
            opts['insertsize'] = get_insertsize(self.qc.option('before_qc_stat').prop['path'], self.gene_set.work_dir)
        else:
            opts['insertsize'] = get_insertsize(self.option('raw_info').prop['path'], self.gene_set.work_dir)
        if self.option('rm_host'):
            opts['QC_dir'] = self.rm_host.option('result_fq_dir')
        elif self.option('qc'):
            opts['QC_dir'] = self.qc.option('after_qc_dir')
        else:
            opts['QC_dir'] = self.option('in_fastq')
        self.set_run(opts, self.gene_set, 'gene_set', self.step.gene_set)
        self.anno_table['geneset'] = os.path.join(self.gene_set.output_dir,
                                                  'gene_profile/RPKM.xls')  # self.gene_set.option('rpkm_abundance').prop['path']

    def run_nr(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            # '/mnt/ilustre/users/sanger-dev/workspace/20170921/MetaGenomic_metagenome/UniGene/output/uniGeneset/gene.uniGeneset.faa',
            'query_type': "prot",
            'database': 'nr',
            'lines': '50000',
        }
        self.set_run(opts, self.nr, 'nr', self.step.nr_)

    def run_kegg(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            # '/mnt/ilustre/users/sanger-dev/workspace/20170921/MetaGenomic_metagenome/UniGene/output/uniGeneset/gene.uniGeneset.faa',
            'query_type': "prot",
            'database': 'kegg',
            'lines': '100000',
        }
        self.set_run(opts, self.kegg, 'kegg', self.step.kegg)

    def run_cog(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            # '/mnt/ilustre/users/sanger-dev/workspace/20170921/MetaGenomic_metagenome/UniGene/output/uniGeneset/gene.uniGeneset.faa',
            'query_type': "prot",
            'database': 'eggnog',
            'lines': '100000',
        }
        self.set_run(opts, self.cog, 'cog', self.step.cog)

    def run_anno(self):
        opts = {
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            # self.gene_set.option('rpkm_abundance'),  # '/mnt/ilustre/users/sanger-dev/workspace/20170921/MetaGenomic_metagenome/UniGene/output/gene_profile/RPKM.xls'
        }
        if self.option('nr'):
            opts['nr_xml_dir'] = self.nr.option('outxml_dir')
        if self.option('kegg'):
            opts['kegg_xml_dir'] = self.kegg.option('outxml_dir')
        if self.option('cog'):
            opts['cog_xml_dir'] = self.cog.option('outxml_dir')
        self.set_run(opts, self.anno, 'anno', self.step.anno, False)
        if self.option('nr'):
            self.nr_dir = os.path.join(self.anno.output_dir, 'nr_tax_level')
            self.anno_table['nr'] = os.path.join(self.nr_dir, 'gene_nr_anno.xls')
        if self.option('cog'):
            self.cog_dir = os.path.join(self.anno.output_dir, 'cog_result_dir')
            self.anno_table['cog'] = os.path.join(self.cog_dir, 'gene_cog_anno.xls')
        if self.option('kegg'):
            self.kegg_dir = os.path.join(self.anno.output_dir, 'kegg_result_dir')
            self.anno_table['kegg'] = os.path.join(self.kegg_dir, 'gene_kegg_anno.xls')
        self.anno.run()

    def run_cazy(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
        }
        self.set_run(opts, self.cazy, 'cazy', self.step.cazy, False)
        self.anno_table['cazy'] = os.path.join(self.cazy.output_dir, 'anno_result', 'gene_cazy_anno.xls')
        self.cazy.run()

    def run_vfdb(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            'lines': '400000',
        }
        self.set_run(opts, self.vfdb, 'vfdb', self.step.vfdb, False)
        self.anno_table['vfdb'] = os.path.join(self.vfdb.output_dir, 'gene_vfdb_total_anno.xls')
        self.vfdb.run()

    def run_ardb(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            'lines': '400000',
        }
        self.set_run(opts, self.ardb, 'ardb', self.step.ardb, False)
        self.anno_table['ardb'] = os.path.join(self.ardb.output_dir, 'gene_ardb_anno.xls')
        self.ardb.run()

    def run_card(self):
        opts = {
            'query': self.gene_set.option('uni_fastaa'),
            'reads_profile_table': self.gene_set.option('reads_abundance'),
            'lines': '400000',
        }
        self.set_run(opts, self.card, 'card', self.step.card, False)
        self.anno_table['card'] = os.path.join(self.card.output_dir, 'gene_card_anno.xls')
        self.card.run()

    """
    def run_analysis(self, event):
        for db in self.choose_anno:
            if type(event) is not str:
                event = event['data']
            self.profile_table1[db] = self.run_new_table(self.anno_table[db], self.anno_table['geneset'],
                                                         self.default_level1[db])
            if self.default_level2[db] == self.default_level1[db] and event == 'all':
                self.profile_table2[db] = self.profile_table1[db]
            elif self.default_level2[db] != self.default_level1[db] and event == 'all':
                self.profile_table2[db] = self.run_new_table(self.anno_table[db], self.anno_table['geneset'],
                                                         self.default_level2[db])
        self.profile_table1['gene'] = self.run_new_table(self.output_dir + '/geneset/gene_profile/gene.uniGeneset.fa.length.txt', self.anno_table['geneset'], "")
        if event == 'all':
            self.profile_table2['gene'] = self.profile_table1['gene']
        if len(self.new_table) != 1:
            self.on_rely(self.new_table, self.run_analysis2)
            for module in self.new_table:
                module.run()
        else:
            self.table.on('end', self.run_analysis2)
            self.table.run()
    """

    def run_analysis(self, event):
        if type(event) is not str:
                event = event['data']
        for db in self.choose_anno:
            ana_str = '1'
            # if event != 'all':
            #     ana_str = '1'
            if event == 'all' and self.default_level2[db] == self.default_level1[db]:
                ana_str = '1,2'
            elif event == 'all' and self.default_level2[db] != self.default_level1[db]:
                # ana_str = '1'
                self.run_new_table(self.anno_table[db], self.anno_table['geneset'], self.default_level2[db], db=db, ana_str='2')
            self.run_new_table(self.anno_table[db], self.anno_table['geneset'], self.default_level1[db], db=db, ana_str=ana_str)
        # gene_str = '1,2' if event == 'all' else '1'
        # self.run_new_table(self.gene_set.output_dir + '/gene_profile/gene.uniGeneset.fa.length.txt', self.anno_table['geneset'], '', db='gene', ana_str=gene_str)
        if len(self.new_table) != 1:
            self.on_rely(self.new_table, self.run_analysis3)
            for module in self.new_table:
                module.run()
        else:
            self.table.on('end', self.run_analysis3)
            self.table.run()

    def run_analysis2(self):
        # self.profile_table1['gene'] = self.anno_table['geneset']
        # self.profile_table2['gene'] = self.anno_table['geneset']
        for db in self.profile_table1.keys():
            self.func_composition(self.profile_table1[db], self.option('group'))
            self.composition_dir2anno[self.composition.output_dir] = db
        for db in self.profile_table2.keys():
            self.func_compare(self.profile_table2[db], self.option('group'))
            self.compare_dir2anno[self.compare.output_dir] = db
            if self.option('envtable').is_set:
                self.func_correlation(self.profile_table2[db], self.option('envtable'))
                self.correlation_dir2anno[self.correlation.output_dir] = db
                self.anno2correlation_tree[db] = self.correlation.work_dir
        if len(self.analysis) != 1:
            self.on_rely(self.analysis, self.end)
            for module in self.analysis:
                module.run()
        else:
            self.composition.on('end', self.end)
            self.composition.run()

    def run_analysis3(self):
        composition_path = os.path.join(self.work_dir, 'new_abund_file')
        other_path = os.path.join(self.work_dir, 'new_abund_file2')
        check_path = self.check_abund_file(composition_path)
        self.logger.info('--------check_path--------')
        self.logger.info(check_path)
        for db in check_path:
            profile_table1 = os.path.join(self.work_dir, 'new_abund_file', db, 'new_abund_table.xls')
            if db in ['nr', 'gene', 'ardb', 'card']:
                self.func_composition(profile_table1, self.option('group'))
            elif db in ['cog', 'kegg', 'cazy', 'vfdb']:
                self.func_composition(profile_table1, self.option('group'), analysis='heatmap,circos')
                self.composition_dir2anno[self.composition.output_dir] = db
                self.func_composition(profile_table1, self.option('group'), analysis='bar', others=0)
            self.composition_dir2anno[self.composition.output_dir] = db
        if os.path.exists(other_path):
            check_path = self.check_abund_file(other_path)
            for db in check_path:
                profile_table2 = os.path.join(self.work_dir, 'new_abund_file2', db, 'new_abund_table.xls')
                self.func_compare(profile_table2, self.option('group'), db)
                self.compare_dir2anno[self.compare.output_dir] = db
                if self.option('envtable').is_set:
                    self.func_correlation(profile_table2, self.option('envtable'))
                    self.correlation_dir2anno[self.correlation.output_dir] = db
                    self.anno2correlation_tree[db] = self.correlation.work_dir
        if len(self.analysis) != 1:
            self.on_rely(self.analysis, self.end)
            for module in self.analysis:
                module.run()
        else:
            self.composition.on("end", self.end)
            self.composition.run()

    def check_abund_file(self, path):
        check_path = os.listdir(path)
        self.logger.info('>>>in check_abund_file, check_path is ')
        self.logger.info(check_path)
        interrupt = 'wait'
        while interrupt == "wait":
            interrupt = 'run'
            for db in self.choose_anno:  # 注： 不含gene
                if db not in check_path:
                    interrupt = 'wait'
                    time.sleep(10)
                    break
                else:
                    file = os.path.join(path, db, 'new_abund_table.xls')
                    if not os.path.isfile(file):
                        interrupt = 'wait'
                        time.sleep(10)
                        break
        return check_path

    def run_new_table(self, anno, gene, level, db, ana_str='1,2'):
        self.logger.info(">>>in run_new_table db is %s ; ana_str is %s" % (db, ana_str))
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
            # self.logger.info("db is %s, ana is %s" % (db,ana))
            db_ana_pointer = db + ','+ ana + ',' + str(self.table_pointer)
            # self.logger.info("db_ana is "  + db_ana_pointer)
            self.table[self.table_pointer].on('end', self.set_abund_file, db_ana_pointer)
        self.new_table.append(self.table[self.table_pointer])
        # self.logger.info("pointer %s has route %s" % (self.table_pointer, self.table[self.table_pointer].output_dir))
        self.table_pointer += 1
        # new_table_file = self.table.output_dir + '/new_abund_table.xls'
        # return new_table_file

    def set_abund_file(self, event):
        # self.logger.info("set_abund_file logger")
        # self.logger.info(event)
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
        # if ana == 1:
        #     target_path = os.path.join(abund_file_path, 'new_abund_table.xls')
        #     path_2 = os.path.join(abund_file_path, 'new_abund_table2.xls')
        #     if not os.path.exists(path_2):
        #         os.link(old_table_file, path_2)
        # elif ana == 2:
        #     target_path = os.path.join(abund_file_path, 'new_abund_table2.xls')
        if os.path.exists(target_path):
            os.remove(target_path)
        os.link(old_table_file, target_path)

    def func_composition(self, abund, group, analysis='bar,heatmap,circos', others=0.01):
        opts = {
            # 'analysis': 'bar,heatmap,circos',
            'analysis': analysis,
            'abundtable': abund,
            'group': group,
            'species_number': '50',
            'others': others,
        }
        self.logger.info('abundtable is :' + abund)
        self.logger.info('group is : ' + group.prop['path'])
        self.composition = self.add_module('meta.composition.composition_analysis')
        self.set_run(opts, self.composition, 'composition', self.step.composition, False)
        self.analysis.append(self.composition)

    def func_compare(self, abund, group, db):
        opts = {
            'dis_method': 'bray_curtis',
            'otutable': abund,
            'group': group,
        }
        if self.option('envtable').is_set and db != 'gene':
            opts['envtable'] = self.option('envtable')
            opts['analysis'] = 'distance,pca,pcoa,nmds,rda_cca,dbrda,hcluster'
        else:
            opts['analysis'] = 'distance,pca,pcoa,nmds,hcluster'
        self.compare = self.add_module('meta.beta_diversity.beta_diversity')
        self.set_run(opts, self.compare, 'compare', self.step.compare, False)
        self.analysis.append(self.compare)

    def func_correlation(self, abund, envtable):
        opts = {
            'method': 'spearmanr',
            'otutable': abund,
            'envtable': envtable,
            "top_species": 50,
        }
        self.correlation = self.add_tool('statistical.pearsons_correlation')
        self.set_run(opts, self.correlation, 'correlation', self.step.correlation, False)
        self.analysis.append(self.correlation)

    '''处理输出文件'''

    def set_output(self, event):
        """
        将各个模块的结果输出至output
        """
        obj = event['bind_object']
        if event['data'] == 'sequence':
            self.move_dir(obj.output_dir, 'rawdata')
        if event['data'] == 'qc':
            self.move_dir(os.path.join(obj.output_dir, 'after_qc_dir'), 'qc/after_qc_dir')
            self.move_dir(os.path.join(obj.output_dir, 'qc_stat'), 'qc/qc_stat')
        if event['data'] == 'rm_host':
            self.move_dir(obj.output_dir, 'rm_host')
        if event['data'] == 'assem':
            self.move_dir(obj.output_dir, 'assemble')
        if event['data'] == 'gene_predict':
            self.move_dir(obj.output_dir, 'predict')
        if event['data'] == 'gene_set':
            self.move_dir(obj.output_dir, 'geneset')
        if event['data'] == 'anno':
            if self.option('nr'):
                # self.nr_dir = os.path.join(obj.output_dir, 'nr_tax_level')
                # self.anno_table['nr'] = os.path.join(self.nr_dir, 'gene_nr_anno.xls')
                self.move_dir(self.nr_dir, 'nr')
            if self.option('cog'):
                # self.cog_dir = os.path.join(obj.output_dir, 'cog_result_dir')
                # self.anno_table['cog'] = os.path.join(self.cog_dir, 'gene_cog_anno.xls')
                self.move_dir(self.cog_dir, 'cog')
            if self.option('kegg'):
                # self.kegg_dir = os.path.join(obj.output_dir, 'kegg_result_dir')
                # self.anno_table['kegg'] = os.path.join(self.kegg_dir, 'gene_kegg_anno.xls')
                self.move_dir(self.kegg_dir, 'kegg')
        if event['data'] == 'cazy':
            # self.anno_table['cazy'] = os.path.join(obj.output_dir, 'anno_result', 'gene_cazy_anno.xls')
            self.move_dir(os.path.join(obj.output_dir, 'anno_result'), 'cazy')
        if event['data'] == 'vfdb':
            # self.anno_table['vfdb'] = os.path.join(obj.output_dir, 'gene_vfdb_total_anno.xls')
            self.move_dir(obj.output_dir, 'vfdb')
        if event['data'] == 'ardb':
            # self.anno_table['ardb'] = os.path.join(obj.output_dir, 'gene_ardb_anno.xls')
            self.move_dir(obj.output_dir, 'ardb')
        if event['data'] == 'card':
            # self.anno_table['card'] = os.path.join(obj.output_dir, 'gene_card_anno.xls')
            self.move_dir(obj.output_dir, 'card')
        if event['data'] == 'composition':
            anno = self.composition_dir2anno[obj.output_dir]
            allfiles = os.listdir(obj.output_dir)
            for dir in allfiles:
                self.move_dir(os.path.join(obj.output_dir, dir), os.path.join('composition', dir, anno))
        if event['data'] == 'compare':
            anno = self.compare_dir2anno[obj.output_dir]
            allfiles = os.listdir(obj.output_dir)
            for dir in allfiles:
                if dir in ['Pca', 'Pcoa', 'Hcluster', 'Nmds', 'Distance']:
                    self.move_dir(os.path.join(obj.output_dir, dir), os.path.join('compare', dir, anno))
                else:
                    self.move_dir(os.path.join(obj.output_dir, dir), os.path.join('correlation', dir, anno))
        if event['data'] == 'correlation':  # ouput里面是一个路径？还是一组文件？
            anno = self.correlation_dir2anno[obj.output_dir]
            new_dir = os.path.join('correlation', 'cor_heatmap', anno)
            self.move_dir(obj.output_dir, new_dir)
            os.rename(os.path.join(self.output_dir, new_dir, "pearsons_correlation.xls"), os.path.join(self.output_dir, new_dir, "correlation.xls"))
            os.rename(os.path.join(self.output_dir, new_dir, "pearsons_pvalue.xls"), os.path.join(self.output_dir, new_dir, "pvalue.xls"))

    def set_output_all(self):
        """
        将所有结果一起导出至output
        """
        pass

    def move_dir(self, olddir, newname):  # 原函数名move2outputdir
        """
        移动一个目录下所有文件/文件夹到workflow输出路径下，供set_output调用
        """
        start = time.time()
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
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
        self.logger.info('start export api')
        self.logger.info('IMPORT_DATA_AFTER_END is :')
        self.logger.info(self.IMPORT_REPORT_AFTER_END)
        self.logger.info("test is : " )
        self.logger.info(self.option('test'))
        self.run_api(test=self.option('test'))
        self.logger.info("start send_files")
        self.send_files()
        # self.rm_tmp_file()
        super(MetaGenomicWorkflow, self).end()

    def send_files(self):  # 原modify_output
        """
        结果放置到/upload_results
        """
        dir_o = self.output_dir
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        if os.path.exists(os.path.join(dir_o, "rawdata")):
            self.move_file(os.path.join(dir_o, "rawdata/base_info"), os.path.join(dir_up, "rawdata/base_info"))
            self.move_file(os.path.join(dir_o, "qc/qc_stat/reads.rawData.stat"),
                           os.path.join(dir_up, "rawdata/data/reads.rawData.stat.xls"))
            self.move_file(os.path.join(dir_o, "qc/qc_stat/reads.cleanData.stat"),
                           os.path.join(dir_up, "qc/qc_stat/reads.cleanData.stat.xls"))
        if os.path.exists(os.path.join(dir_o, "qc")):
            self.move_file(os.path.join(dir_o, "assemble/assembly.stat"),
                           os.path.join(dir_up, "assemble/assembly.stat.xls"))
            files = os.listdir(os.path.join(dir_o, "qc"))
            for file in files:
                if file.endswith("contig.fa"):
                    self.move_file(os.path.join(dir_o, "qc", file), os.path.join(dir_up, "qc", file))
        if os.path.exists(os.path.join(dir_o, "predict")):
            self.move_file(os.path.join(dir_o, "predict/sample.metagene.stat"),
                           os.path.join(dir_up, "predict/genePredict_stat.xls"))
            files = os.listdir(os.path.join(dir_o, "predict"))
            for file in files:
                new_file = file.replace("metagene.fa", "genePredict.fa")
                if file.endswith("metagene.fa"):
                    self.move_file(os.path.join(dir_o, "predict", file), os.path.join(dir_up, "predict", new_file))
        if os.path.exists(os.path.join(dir_o, "geneset")):
            self.move_file(os.path.join(dir_o, "geneset/uniGeneset"), os.path.join(dir_up, "geneset/uniGeneset"))
            self.move_file(os.path.join(dir_o, "geneset/gene_profile"), os.path.join(dir_up, "geneset/gene_profile"))
            # os.remove(os.path.join(dir_up, "geneset/gene_profile/gene.uniGeneset.fa.length.txt"))
        for file in ["nr", "kegg", "cog", "cazy", "vfdb", "ardb", "card", "composition", "compare", "correlation"]:
            if os.path.exists(os.path.join(dir_o, file)):
                self.move_file(os.path.join(dir_o, file), os.path.join(dir_up, file))
            if file == "kegg":
                for rm_file in ["kegg_merge.xml", "pathway.kgml", "pathway.png"]:
                    rm_file_path = os.path.join(dir_up, file, rm_file)
                    if os.path.isfile(rm_file_path):
                        os.remove(rm_file_path)
        self.logger.info("repaths...")
        repaths = [
            [".", "", "流程分析结果目录"],
            ["rawdata", "", "原始序列目录"],
            ["rawdata/base_info", "", "原始序列质量统计目录"],
            ["rawdata/data", "", "原始序列文件路径"],
            ["rawdata/data/reads.rawData.stat.xls", "", "各样品原始数据统计表"],
            ["qc", "", "质控结果目录"],
            ["qc/qc_stat", "", "质控后质量统计结果"],
            ["qc/qc_stat/reads.cleanData.stat.xls", "", "各样品过滤后数据统计表"],
            ["rm_host", "", "去宿主后序列目录"],
            ["rm_host/stat.list.txt", "", "各样品去宿主后序列统计表"],
            ["assemble", "", "拼接结果目录"],
            ["assemble/assembly.stat.xls", "", "各步骤组装结果统计表"],
            ["predict", "", "基因预测结果目录"],
            ["predict/genePredict_stat.xls", "", "基因预测结果统计表"],
            ["geneset", "", "非冗余基因集结果目录"],
            ["geneset/gene_profile", "", "非冗余集因丰度目录"],
            ['geneset/gene_profile/gene.uniGeneset.fa.length.txt', '', '基因长度表'],
            ["geneset/gene_profile/RPKM.xls", "", "基因在各个样品中的RPKM丰度表"],
            ["geneset/gene_profile/reads_length_ratio_relative.xls", "", "基因在各个样品中的相对丰度/基因长度表"],
            ["geneset/gene_profile/reads_length_ratio.xls", "", "基因在各个样品中的丰度/基因长度表"],
            ["geneset/gene_profile/reads_number_relative.xls", "", "基因在各个样品中的相对丰度表"],
            ["geneset/gene_profile/reads_number.xls", "", "基因在各个样品中的丰度表"],
            ["geneset/gene_profile/TPM.xls", "", "基因在各个样品中的TPM丰度表"],
            ['geneset/gene_profile/top100_reads_number.xls', '', '丰度前100的基因在各个样品中的丰度表'],
            ['geneset/gene_profile/top100_reads_number_relative.xls', '', '丰度前100的基因在各个样品中的相对丰度表'],
            ['geneset/gene_profile/reads_profile.tar.gz', '', '基因reads数相对丰度与基因reads数丰度的压缩文件'],
            ["geneset/uniGeneset", "", "非冗余基因集序列统计目录"],
            ["geneset/uniGeneset/geneCatalog_stat.xls", "", "去冗余前后基因数目和长度统计表"],
            ["geneset/uniGeneset/gene.uniGeneset.fa", "", "非冗余基因集核酸序列"],
            ["geneset/uniGeneset/gene.uniGeneset.faa", "", "非冗余基因集蛋白序列"],
            ["nr", "", "NR功能注释结果目录"],
            ["nr/gene_nr_anno.xls", "", "每条基因的物种注释表"],
            ["nr/nr_align_table.xls", "", "物种序列比对结果"],
            ["nr/tax_d.xls", "", "域注释丰度表"],
            ["nr/tax_k.xls", "", "界注释丰度表"],
            ["nr/tax_p.xls", "", "门注释丰度表"],
            ["nr/tax_c.xls", "", "纲注释丰度表"],
            ["nr/tax_o.xls", "", "目注释丰度表"],
            ["nr/tax_f.xls", "", "科注释丰度表"],
            ["nr/tax_g.xls", "", "属注释丰度表"],
            ["nr/tax_s.xls", "", "种注释丰度表"],
            ["kegg", "", "KEGG功能注释结果目录"],
            ["kegg/gene_kegg_anno.xls", "", "每条基因的KEGG功能注释表"],
            ["kegg/kegg_align_table.xls", "", "KEGG序列比对结果表"],
            ["kegg/kegg_pathway_eachmap.xls", "", "各样本Pathway map表"],
            ["kegg/kegg_enzyme_profile.xls", "", "各样品KEGG酶丰度表"],
            ["kegg/kegg_gene_profile.xls", "", "各样品KEGG基因丰度表"],
            ["kegg/kegg_KO_profile.xls", "", "各样品KO丰度表"],
            ["kegg/kegg_level1_profile.xls", "", "各样品Pathway level1丰度表"],
            ["kegg/kegg_level2_profile.xls", "", "各样品Pathway level2丰度表"],
            ["kegg/kegg_level3_profile.xls", "", "各样品Pathway level3丰度表"],
            ["kegg/kegg_module_profile.xls", "", "各样品KEGG Module丰度表"],
            ["kegg/kegg_pathway_profile.xls", "", "各样品Pathway丰度表"],
            ["cog", "", "COG功能注释结果目录"],
            ["cog/cog_align_table.xls", "", "COG序列比对结果表"],
            ["cog/cog_nog_profile.xls", "", "各样品COG NOG丰度表"],
            ["cog/cog_category_profile.xls", "", "各样品COG Category丰度表"],
            ["cog/cog_function_profile.xls", "", "各样品COG Function丰度表"],
            ["cog/gene_cog_anno.xls", "", "每条基因的COG功能注释表"],
            ["cazy", "", "CAZy碳水化合物活性酶注释结果目录"],
            ["cazy/cazy_class_profile.xls", "", "各样品CAZY Class丰度表"],
            ["cazy/gene_cazy_anno.xls", "", "每条基因的CAZY功能注释表"],
            ["cazy/cazy_family_profile.xls", "", "各样品CAZY Family丰度表"],
            ["cazy/gene_dbCAN.hmmscan.out.dm.ds", "", "CAZY序列比对结果表"],
            ["vfdb", "", "VFDB毒力因子注释结果目录"],
            ["vfdb/gene_vfdb_core_anno.xls", "", "每条基因的VFDB核心库功能注释表"],
            ["vfdb/gene_vfdb_predict_anno.xls", "", "每条基因的VFDB预测库功能注释表"],
            ["vfdb/gene_vfdb_total_anno.xls", "", "每条基因的VFDB功能注释表"],
            ["vfdb/vfdb_all_Gi_profile.xls", "", "各样品VFDB基因丰度表"],
            ["vfdb/vfdb_all_VF_profile.xls", "", "各样品VFDB毒力因子丰度表"],
            ["vfdb/vfdb_core_align_table.xls", "", "VFDB核心库序列比对结果"],
            ["vfdb/vfdb_core_Gi_profile.xls", "", "各样品VFDB核心库基因丰度表"],
            ["vfdb/vfdb_core_VF_profile.xls", "", "各样品VFDB核心库毒力因子丰度表"],
            ["vfdb/vfdb_predict_align_table.xls", "", "VFDB预测库序列比对结果表"],
            ["vfdb/vfdb_predict_Gi_profile.xls", "", "各样品VFDB预测库基因丰度表"],
            ["vfdb/vfdb_predict_VF_profile.xls", "", "各样品VFDB预测库毒力因子丰度表"],
            ["vfdb/vfdb_level_pie.xls", "", "VFDB两级分类的丰度统计表"],
            ["ardb", "", "ARDB抗性基因功能注释结果目录"],
            ["ardb/ardb_align_table.xls", "", "ARDB序列比对结果表"],
            ["ardb/ardb_class_profile.xls", "", "各样品ARDB Class丰度表"],
            ["ardb/gene_ardb_anno.xls", "", "每条基因的ARDB功能注释表"],
            ["ardb/ardb_ARG_profile.xls", "", "各样品ARDB ARG丰度表"],
            ["ardb/ardb_type_profile.xls", "", "各样品ARDB type丰度表"],
            ['ardb/gene_ardb_class_stat.xls', '', 'class基因信息统计表'],
            ["card", "", "CARD抗性基因功能注释结果目录"],
            ["card/card_align_table.xls", "", "CARD序列比对结果"],
            ["card/card_class_profile.xls", "", "各样品CARD Class丰度表"],
            ["card/card_ARO_gene_number.xls", "", "CARD每个ARO比对基因信息表"],
            ["card/gene_card_anno.xls", "", "每条基因的CARD功能注释表"],
            ["card/card_ARO_profile.xls", "", "各样品CARD ARO丰度表"],
            ["composition", "", "物种与功能组成分析结果目录"],
            ["composition/bar", "", "柱形图结果目录"],
            ["composition/heatmap", "", "Heatmap图结果目录"],
            ["composition/circos", "", "Circos样本与物种或功能关系图结果目录"],
            ["compare", "", "物种与功能比较分析结果目录"],
            ["compare/Hcluster", "", "样本层次聚类分析结果目录"],
            ["compare/Distance", "", "距离矩阵计算结果目录"],
            ["compare/Pca", "", "PCA分析结果目录"],
            ["compare/Pcoa", "", "PCoA分析结果目录"],
            ["compare/Nmds", "", "NMDS分析结果目录"],
            ["correlation", "", "环境因子关联分析结果目录"],
            ["correlation/Rda", "", "RDA_CCA分析结果目录"],
            ["correlation/Dbrda", "", "db_RDA分析结果目录"],
            ["correlation/cor_heatmap", "", "相关性Heatmap分析结果目录"],
        ]
        self.logger.info("regexps...")
        regexps = [
            [r"rawdata/base_info/.+\.fastq\.fastxstat\.txt", "", "原始序列碱基质量统计文件"],
            [r"assemble/[^/]+\.contig\.fa", "", "拼接contig结果"],
            [r"predict/.+\.genePredict.fa", "", "长度大于等于100bp的基因的核酸序列"],
            [r"geneset/gene_profile/top100\..+xls", "", "丰度前100的基因丰度表"],
            [r"composition/.+/.+/taxa\.percents\.table\.xls", "", "物种/基因/功能相对丰度结果表"],
            [r"composition/.+/.+/taxa\.table\.xls", "", "物种/基因/功能丰度结果表"],
            [r"compare/Hcluster/.+/hcluster\.tre", "graph.newick_tree", "样本层次聚类树结果表"],
            [r"compare/Distance/.+/bray_curtis.+\.xls$", "meta.beta_diversity.distance_matrix", "样本距离矩阵文件"],
            [r"compare/Nmds/.+/nmds_sites\.xls$", "xls", "NMDS样本各维度坐标"],
            [r"compare/Nmds/.+/nmds_stress\.xls$", "xls", "NMDS样本特征拟合度值"],
            [r"compare/Pca/.+/pca_importance\.xls$", "xls", "PCA主成分解释度表"],
            [r"compare/Pca/.+/pca_sites\.xls$", "xls", "PCA样本各成分轴坐标"],
            [r"compare/Pca/.+/pca_rotation\.xls$", "xls", "PCA主成分贡献度表"],
            [r"compare/Pca/.+/pca_rotation_all\.xls$", "xls", "PCA全部主成分贡献度表"],
            [r"compare/Pca/.+/pca_envfit_vector\.xls$", "xls", "PCA数量型环境因子的显著性检验值"],
            [r"compare/Pca/.+/pca_envfit_vector_scores\.xls$", "xls", "PCA数量型环境因子坐标表"],
            [r"compare/Pca/.+/pca_envfit_factor\.xls$", "xls", "PCA哑变量环境因子的显著性检验值"],
            [r"compare/Pca/.+/pca_envfit_factor_scores\.xls$", "xls", "PCA哑变量环境因子坐标表"],
            [r"compare/Pcoa/.+/pcoa_sites\.xls", "xls", "PCoA样本坐标表"],
            [r"compare/Pcoa/.+/pcoa_eigenvalues\.xls", "xls", "PCoA矩阵特征值"],
            [r"compare/Pcoa/.+/pcoa_eigenvaluespre\.xls", "xls", "PCoA特征解释度百分比"],
            [r"correlation/Rda/.+/dca\.xls", "xls", "判断用RDA还是CCA的DCA文件"],
            [r"correlation/Rda/.+/.*rda_biplot\.xls", "xls", "RDA_CCA数量型环境因子坐标表"],
            [r"correlation/Dbrda/.+/.*rda_biplot\.xls", "xls", "dbRDA数量型环境因子坐标表"],
            [r"correlation/Rda/.+/rda_envfit\.xls", "xls", "RDA_CCA各环境因子的显著性检验值"],
            [r"correlation/Rda/.+/.*rda_importance\.xls", "xls", "RDA_CCA主成分解释度表"],
            [r"correlation/Dbrda/.+/.*rda_importance\.xls", "xls", "dbRDA主成分解释度表"],
            [r"correlation/Rda/.+/.*rda_plot_species_data\.xls", "xls", "RDA_CCA绘图物种_功能坐标表"],
            [r"correlation/Dbrda/.+/.*rda_plot_species_data\.xls", "xls", "dbRDA绘图物种_功能坐标表"],
            [r"correlation/Rda/.+/.*rda_sites\.xls", "xls", "RDA_CCA样本坐标表"],
            [r"correlation/Dbrda/.+/.*rda_sites\.xls", "xls", "dbRDA样本坐标表"],
            [r"correlation/Rda/.+/.*rda_species\.xls", "xls", "RDA_CCA物种_功能坐标表"],
            [r"correlation/Dbrda/.+/.*rda_species\.xls", "xls", "dbRDA物种_功能坐标表"],
            [r"correlation/cor_heatmap/.+/correlation\.xls", "xls", "相关性Heatmap图的Pearson相关系数表"],
            [r"correlation/cor_heatmap/.+/pvalue\.xls", "xls", "相关性Heatmap图的Pearson相关系数对应p值表"],
        ]
        self.logger.info("add_upload_dir")
        sdir = self.add_upload_dir(dir_up)
        self.logger.info("add_relpath_rules...")
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

    '''导表'''

    def run_api(self, test=False):  # 原run_api_and_set_output
        self.logger.info("run_api")
        self.logger.info("--------")
        # greenlets_list_first = []  # 一阶段导表
        # greenlets_list_sec = []  # 二阶段导表
        # greenlets_list_third = []  # 三阶段导表
        if test:
            self.logger.info("test is %s" % test)
            self.export_qc()  # self.specimen_group  self.group_detail self.spname_spid
            self.export_assem()
            self.export_predict()
            self.export_geneset()  # self.geneset_id
            self.export_ardb()
            self.export_nr()
            self.export_kegg()
            self.export_cog()
            self.export_vfdb()
            self.export_cazy()
            self.export_card()
            self.export_overview()
            self.anno_table = {
                "geneset": "",
                "card": "",
                "nr": "",
                "kegg": "",
                "cog": "",
                "vfdb": "",
                "cazy": "",
                "ardb": "",
            }
            self.geneset_id = str(self.geneset_id)
            self.env_id = str(self.env_id)
            self.anno2correlation_tree = {
                "gene": "/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/PearsonsCorrelation2",
                "nr": "/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/PearsonsCorrelation6",
                "kegg": "/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/PearsonsCorrelation1",
                "cog": "/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/PearsonsCorrelation4",
                "vfdb": "/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/PearsonsCorrelation3",
                "cazy": "/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/PearsonsCorrelation7",
                "card": "/mnt/ilustre/users/sanger-dev/workspace/20171024/MetaGenomic_metagenome/PearsonsCorrelation8",
                "ardb": "/mnt/ilustre/users/sanger-dev/workspace/20171106/MetaGenomic_metagenome_table/PearsonsCorrelation1",
            }
            ana_list = ["composition", "compare", "correlation"]
            for ana in ana_list:
                ana_path = os.path.join(self.output_dir, ana)
                if os.path.isdir(ana_path):
                    analysis_list = os.listdir(ana_path)
                    for analysis in analysis_list:
                        if ana == "composition":
                            self.export_composition(analysis)
                        if ana == "compare":
                            if analysis == "Hcluster":
                                self.export_hclust()
                            elif analysis in ["Nmds", "Pca", "Pcoa"]:
                                self.export_beta(analysis)
                        if ana == "correlation":
                            if analysis == "cor_heatmap":
                                self.export_cor_heatmap()
                            elif analysis in "Dbrda, Rda":
                                self.logger.info("analysis is :" + analysis)
                                self.export_beta(analysis)
            self.logger.info("导表测试完成")
            return
        else:
            self.logger.info("what's wrong?")
        self.logger.info("test is not, run export_qc")
        self.export_qc()
        self.logger.info("export_qc end")
        self.export_assem()
        self.export_predict()
        self.export_geneset()
        self.export_nr()
        self.export_kegg()
        self.export_cog()
        self.export_cazy()
        self.export_vfdb()
        self.export_ardb()
        self.export_card()
        self.export_overview()
        self.geneset_id = str(self.geneset_id)
        self.env_id = str(self.env_id)
        ana_list = ["composition", "compare", "correlation"]
        for ana in ana_list:
            ana_path = os.path.join(self.output_dir, ana)
            if os.path.isdir(ana_path):
                analysis_list = os.listdir(ana_path)
                for analysis in analysis_list:
                    if ana == "composition":
                        self.export_composition(analysis)
                    if ana == "compare":
                        if analysis == "Hcluster":
                            self.export_hclust()
                        elif analysis in ["Nmds", "Pca", "Pcoa"]:
                            self.export_beta(analysis)
                    if ana == "correlation":
                        if analysis == "cor_heatmap":
                            self.export_cor_heatmap()
                        elif analysis in "Dbrda, Rda":
                            self.logger.info("analysis is :" + analysis)
                            self.export_beta(analysis)
        self.logger.info("导表完成")

    @time_count
    def export_qc(self):
        self.logger.info("in export_qc_func")
        self.api_dic["data_stat"] = self.api.api("metagenomic.data_stat")
        if self.option('qc'):
            self.logger.info("add_data_stat")
            data_stat_id = self.api_dic["data_stat"].add_data_stat("raw",
                                                                   self.output_dir + "/qc/qc_stat/new_reads.rawData.stat",
                                                                   self.output_dir + "/rawdata/base_info", "null")
            self.logger.info("add_data_stat clean")
            self.api_dic["data_stat"].add_data_stat("clean", self.output_dir + "/qc/qc_stat/reads.cleanData.stat",
                                                    "null", data_stat_id)
        else:
            data_stat_id = self.api_dic["data_stat"].add_data_stat('raw', self.option('raw_info').prop['path'],
                                                                   self.base_info, "null")
            self.api_dic["data_stat"].add_data_stat("clean", self.option('qc_info').prop['path'], "null", data_stat_id)
        # self.spname_spid = self.api_dic["data_stat"].get_spname_spid()
        self.logger.info('data_stat_id is : %s :' % data_stat_id)
        self.logger.info('get spanme_spid')
        self.spname_spid = name2id(data_stat_id, type="raw")
        if self.option('rm_host'):
            self.logger.info('rm_host')
            self.api_dic["data_stat"].add_data_stat("optimised", self.output_dir + "/rm_host/stat.list.txt", "null",
                                                    data_stat_id)
        self.logger.info('get group')
        if self.option('group').is_set:
            self.api_dic["group"] = self.api.api("metagenomic.specimen_group")
            group_info = self.api_dic["group"].add_ini_group_table(self.option('group').prop['path'], self.spname_spid)
            if len(group_info) == 1:
                self.specimen_group = group_info[0]["specimen_group"]
                self.group_detail = group_info[0]["group_detail"]
                self.logger.info("specimen_group is :")
                self.logger.info(self.specimen_group)
                self.logger.info("group_detail is :")
                self.logger.info(self.group_detail)
            else:
                self.logger.info("group_info length is %s" % len(group_info))
                return
        if self.option('envtable').is_set:
            self.api_dic["envtable"] = self.api.api("metagenomic.env_metagenomic")
            self.env_id = self.api_dic["envtable"].add_env_table(self.option('envtable').prop['path'], self.spname_spid)
            self.env_labs = self.api_dic["envtable"].get_env_lab()

    @time_count
    def export_assem(self):
        self.api_dic["assem_gene"] = self.api.api("metagenomic.assemble_gene")
        if self.option('assemble_type') == 'SOAPdenovo2':
            assemble_type = 'soapdenovo'
        elif self.option('assemble_type') == 'IDBA_UD':
            assemble_type = 'idba'
        else:
            assemble_type = self.option('assemble_type').lower()
        assem_id = self.api_dic["assem_gene"].add_assemble_stat(assemble_type, self.option("min_contig"))
        self.api_dic["assem_gene"].add_assemble_stat_detail(assem_id, self.output_dir + "/assemble/assembly.stat")
        self.api_dic["assem_gene"].add_assemble_stat_bar(assem_id, self.output_dir + "/assemble/len_distribute")

    @time_count
    def export_predict(self):
        if "assem_gene" not in self.api_dic.keys():
            self.api_dic["assem_gene"] = self.api.api("metagenomic.assemble_gene")
        if self.option('assemble_type') == 'SOAPdenovo2':
            assemble_type = 'soapdenovo'
        elif self.option('assemble_type') == 'IDBA_UD':
            assemble_type = 'idba'
        else:
            assemble_type = self.option('assemble_type').lower()
        gene_id = self.api_dic["assem_gene"].add_predict_gene(assemble_type, self.option("min_gene"))
        self.api_dic["assem_gene"].add_predict_gene_detail(gene_id, self.output_dir + "/predict/sample.metagene.stat")
        self.api_dic["assem_gene"].add_predict_gene_bar(gene_id, self.output_dir + "/predict/len_distribute")
        self.api_dic["assem_gene"].add_predict_gene_total(self.output_dir + "/predict/sample.metagene.stat")

    @time_count
    def export_geneset(self):
        self.api_dic["geneset"] = self.api.api("metagenomic.geneset")
        geneset_path = os.path.join(self.remote_dir, 'geneset')
        self.geneset_id = self.api_dic["geneset"].add_geneset(self.output_dir + "/geneset", geneset_path, 1)  # mark
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
        # return_id = self.api_dic["nr"].add_anno_nr(self.geneset_id, "All", self.output_dir + "/nr/gene_nr_anno.xls", group_id=self.specimen_group, group_detail=self.group_detail)
        return_id = self.api_dic["nr"].add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group, group_detail=self.group_detail)
        self.api_dic["nr"].add_anno_nr_detail(return_id, self.output_dir + "/nr")
        self.anno_id["nr"] = str(return_id)

    @time_count
    def export_cog(self):
        self.api_dic["cog"] = self.api.api("metagenomic.mg_anno_cog")
        cog_path = os.path.join(self.remote_dir, 'cog/gene_cog_anno.xls')
        # return_id = self.api_dic["cog"].add_anno_cog(self.geneset_id, "All", self.output_dir + "/cog/gene_cog_anno.xls", group_id=self.specimen_group, group_detail=self.group_detail)
        return_id = self.api_dic["cog"].add_anno_cog(self.geneset_id, "All", cog_path, group_id=self.specimen_group, group_detail=self.group_detail)
        self.api_dic["cog"].add_anno_cog_nog(return_id, self.output_dir + "/cog")
        self.api_dic["cog"].add_anno_cog_function(return_id, self.output_dir + "/cog")
        self.api_dic["cog"].add_anno_cog_category(return_id, self.output_dir + "/cog")
        self.anno_id["cog"] = str(return_id)

    @time_count
    def export_kegg(self):
        self.api_dic["kegg"] = self.api.api("metagenomic.mg_anno_kegg")
        xml_path = os.path.join(self.output_dir, "kegg/kegg_merge.xml")
        kegg_web_path = os.path.join(self.remote_dir, 'kegg/gene_kegg_anno.xls')
        # return_id = self.api_dic["kegg"].add_anno_kegg(self.geneset_id, "All",
        #                                                self.output_dir + "/kegg/gene_kegg_anno.xls", group_id=self.specimen_group, group_detail=self.group_detail)
        return_id = self.api_dic["kegg"].add_anno_kegg(self.geneset_id, 'All',
                                                       kegg_web_path, xml_file=xml_path , group_id=self.specimen_group, group_detail=self.group_detail)
        self.api_dic["kegg"].add_anno_kegg_gene(return_id, self.output_dir + "/kegg")
        self.api_dic["kegg"].add_anno_kegg_orthology(return_id, self.output_dir + "/kegg")
        self.api_dic["kegg"].add_anno_kegg_module(return_id, self.output_dir + "/kegg")
        self.api_dic["kegg"].add_anno_kegg_enzyme(return_id, self.output_dir + "/kegg")
        graph_path = os.path.join('/'.join(self.remote_dir.split('/')[:5]), 'metag/KEGG_Pathway', str(self._sheet.id), str(return_id) + '/')
        # graph_path += '/rerewrweset/metag/KEGG_Pathway' + str(self._sheet.id) + '/' + str(return_id) + '/'
        self.api_dic["kegg"].add_anno_kegg_pathway(return_id, self.output_dir + "/kegg", graph_path)
        self.anno_id["kegg"] = str(return_id)

    @time_count
    def export_cazy(self):
        self.api_dic["cazy"] = self.api.api("metagenomic.mg_anno_cazy")
        cazy_path = os.path.join(self.remote_dir, 'cazy/gene_cazy_anno.xls')
        # return_id = self.api_dic["cazy"].add_anno_cazy(self.geneset_id, "All",
        #                                                self.output_dir + "/cazy/gene_cazy_anno.xls", group_id=self.specimen_group, group_detail=self.group_detail)
        return_id = self.api_dic["cazy"].add_anno_cazy(self.geneset_id, "All",
                                                       cazy_path, group_id=self.specimen_group, group_detail=self.group_detail)
        self.api_dic["cazy"].add_anno_cazy_family(return_id, self.output_dir + "/cazy")
        self.api_dic["cazy"].add_anno_cazy_class(return_id, self.output_dir + "/cazy")
        self.anno_id["cazy"] = str(return_id)

    @time_count
    def export_vfdb(self):
        self.api_dic["vfdb"] = self.api.api("metagenomic.mg_anno_vfdb")
        vfdb_path = os.path.join(self.remote_dir, 'vfdb/gene_vfdb_total_anno.xls')
        # return_id = self.api_dic["vfdb"].add_anno_vfdb(self.geneset_id, "All",
        #                                                self.output_dir + "/vfdb/gene_vfdb_total_anno.xls", group_id=self.specimen_group, group_detail=self.group_detail)
        return_id = self.api_dic["vfdb"].add_anno_vfdb(self.geneset_id, "All",
                                                       vfdb_path, group_id=self.specimen_group, group_detail=self.group_detail)
        self.api_dic["vfdb"].add_anno_vfdb_vfs(return_id, self.output_dir + "/vfdb", "core")
        self.api_dic["vfdb"].add_anno_vfdb_vfs(return_id, self.output_dir + "/vfdb", "predict")
        self.api_dic["vfdb"].add_anno_vfdb_pie(return_id, self.output_dir + "/vfdb")
        self.anno_id["vfdb"] = str(return_id)

    @time_count
    def export_ardb(self):
        self.api_dic["ardb"] = self.api.api("metagenomic.mg_anno_ardb")
        ardb_path = os.path.join(self.remote_dir, 'ardb/gene_ardb_anno.xls')
        # return_id = self.api_dic["ardb"].add_anno_ardb(self.geneset_id, "All",
        #                                                self.output_dir + "/ardb/gene_ardb_anno.xls", group_id=self.specimen_group, group_detail=self.group_detail)
        return_id = self.api_dic["ardb"].add_anno_ardb(self.geneset_id, "All",
                                                       ardb_path, group_id=self.specimen_group, group_detail=self.group_detail)
        self.api_dic["ardb"].add_anno_ardb_arg(return_id, self.output_dir + "/ardb")
        self.api_dic["ardb"].add_anno_ardb_type(return_id, self.output_dir + "/ardb")
        self.api_dic["ardb"].add_anno_ardb_class(return_id, self.output_dir + "/ardb")
        self.anno_id["ardb"] = str(return_id)

    @time_count
    def export_card(self):
        self.api_dic["card"] = self.api.api("metagenomic.mg_anno_card")
        card_path = os.path.join(self.remote_dir, 'card/gene_card_anno.xls')
        # return_id = self.api_dic["card"].add_anno_card(self.geneset_id, "All",
        #                                                self.output_dir + "/card/gene_card_anno.xls", group_id=self.specimen_group, group_detail=self.group_detail)
        return_id = self.api_dic["card"].add_anno_card(self.geneset_id, "All",
                                                       card_path, group_id=self.specimen_group, group_detail=self.group_detail)
        self.api_dic["card"].add_anno_card_aro(return_id, self.output_dir + "/card")
        self.api_dic["card"].add_anno_card_class(return_id, self.output_dir + "/card")
        self.anno_id["card"] = str(return_id)

    @time_count
    def export_overview(self):
        """
        注釋表信息總覽
        :return:
        """
        self.api_dic["overview"] = self.api.api("metagenomic.mg_anno_overview")
        return_id = self.api_dic["overview"].add_anno_overview(self.geneset_id)
        select_table = {}
        for select_anno in ["nr", "kegg", "cog", "cazy", "vfdb", "ardb", "card"]:
            if select_anno in self.anno_table.keys():
                select_table[select_anno] = self.anno_table[select_anno]
            else:
                select_table[select_anno] = None
        # self.api_dic["overview"].add_anno_overview_detail(return_id, self.anno_table["geneset"], nr=select_table["nr"],
        #                                                   kegg=select_table["kegg"], cog=select_table["cog"],
        #                                                  cazy=select_table["cazy"], vfdb=select_table["vfdb"],
        #                                                  ardb=select_table["ardb"], card=select_table["card"])
        self.api_dic["overview"].add_anno_overview_detail(return_id,
                                                          self.output_dir + "/geneset/gene_profile/gene.uniGeneset.fa.length.txt",
                                                          nr=select_table["nr"],
                                                          kegg=select_table["kegg"], cog=select_table["cog"],
                                                          cazy=select_table["cazy"], vfdb=select_table["vfdb"],
                                                          ardb=select_table["ardb"], card=select_table["card"])

    def get_param(self, anno):
        params = {
            "geneset_id": str(self.geneset_id),
            "group_detail": self.group_detail,
            "group_id": str(self.specimen_group),
            "method": "rpkm",
            "task_type": 2,
            "anno_id": str(self.anno_id[anno]),
            "anno_type": anno,
        }
        return params

    @time_count
    def export_composition(self, analysis):
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
        for select_anno in self.anno_table.keys():
            if select_anno == 'geneset':
                continue
                # select_anno = 'gene'
                # submit_location = analysis_dic[analysis] + "_" + select_anno.upper() + "_" + time
            else:
                submit_location = analysis_dic[analysis] + "_" + select_anno.upper() + "_" + self.default_level1[select_anno] + "_" + time
            if analysis in ["circos", "heatmap"]:
                params["combine_value"] = "0.01" if select_anno in ['nr', 'gene', 'ardb', 'card'] else "0"
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            # self.logger.info("export composition-" + analysis + "-" + select_anno)
            # return_id = self.api_dic["composition"].add_composition(analysis, select_anno, self.geneset_id,
            #                                                         self.anno_id[select_anno], level_id[select_anno],
            #                                                         self.specimen_group,
            #                                                         name=analysis_dic[analysis] + "_" + select_anno.upper() + "_" + self.default_level1[select_anno] + "_" + time,
            #                                                         params=params)
            return_id = self.api_dic["composition"].add_composition(analysis, select_anno, name=submit_location,
                                                                    params=params)
            # self.logger.info("next")
            if analysis == 'heatmap':
                detail_file = os.path.join(self.output_dir, "composition", analysis, select_anno, "taxa.table.xls")
            else:
                detail_file = os.path.join(self.output_dir, "composition", analysis, select_anno, "taxa.percents.table.xls")
            self.api_dic["composition"].add_composition_detail(detail_file, return_id, species_tree="",
                                                               specimen_tree="")

    @time_count
    def export_beta(self, analysis):
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
        params = {
            "analysis_type": analysis,
            "second_level": "",
        }
        analysis_type = analysis.lower()
        analysis_dic = {
            "pca": "PCA",
            "pcoa": "PCoA",
            "nmds": "NMDS",
            "rda": "RDACCA",
            "dbrda": "dbRDA",
        }
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if analysis_type in ["nmds", "pcoa", "dbrda"]:
            params["distance_method"] = "bray_curtis"
        elif analysis_type in ["pca", "rda", "dbrda"] and self.option('envtable').is_set:
            params["env_id"] = str(self.env_id)  #
            params["env_labs"] = self.env_labs  #
        for select_anno in self.anno_table.keys():
            self.logger.info("anno is " + select_anno)
            if select_anno == 'geneset':
                continue
                # select_anno = 'gene'
                # submit_location = analysis_dic[analysis_type] + "_" + select_anno.upper() + "_" + time
            else:
                submit_location = analysis_dic[analysis_type] + "_" + select_anno.upper() + "_" + \
                                        self.default_level2[select_anno] + "_" + time
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            self.logger.info(params)
            if analysis_type in ["nmds", "pcoa"]:
                # self.logger.info("analysis:" + analysis + ",\nanno:" + select_anno)
                path = os.path.join(self.output_dir, "compare", analysis, select_anno)
                # self.api_dic["beta"].add_beta_diversity(path, analysis_type, main="True", anno_type=select_anno,
                #                                         params=params, group_id=self.specimen_group,
                #                                         geneset_id=self.geneset_id, anno_id=self.anno_id[select_anno],
                #                                         level_id=level_id[select_anno])
                self.api_dic["beta"].add_beta_diversity(path, analysis_type, main="True", anno_type=select_anno,
                                                            name=submit_location, params=params)
            elif analysis_type == 'pca':
                path = os.path.join(self.output_dir, "compare", analysis, select_anno)
                pca_path = os.path.join(self.remote_dir, 'compare/Pca/', select_anno)
                if self.env_id != '':
                    self.api_dic["beta"].add_beta_diversity(path, analysis_type, main="True", web_path=pca_path, anno_type=select_anno,
                                                            env_id=self.env_id, name=submit_location, params=params)
                else:
                    self.api_dic["beta"].add_beta_diversity(path,  analysis_type, main="True", web_path=pca_path, anno_type=select_anno,
                                                            name=submit_location, params=params)
            elif analysis_type == "rda":
                path = os.path.join(self.output_dir, "correlation", analysis, select_anno)
                if not os.path.exists(path):
                    self.logger.info("不存在路径%s，跳过%s的rda分析" % (path, select_anno))
                    return
                # self.api_dic["beta"].add_beta_diversity(path, "rda_cca", main="True", anno_type=select_anno,
                #                                         params=params, group_id=self.specimen_group,
                #                                         geneset_id=self.geneset_id, anno_id=self.anno_id[select_anno],
                #                                         level_id=level_id[select_anno])
                self.api_dic["beta"].add_beta_diversity(path, "rda_cca", main="True", anno_type=select_anno,
                                                        env_id=self.env_id, name=submit_location,params=params)
            elif analysis_type == "dbrda":
                path = os.path.join(self.output_dir, "correlation", analysis, select_anno)
                if not os.path.exists(path):
                    self.logger.info("不存在路径%s，跳过%s的dbrda分析" % (path, select_anno))
                    return
                # self.api_dic["beta"].add_beta_diversity(path, "dbrda", main="True", anno_type=select_anno,
                #                                         params=params, group_id=self.specimen_group,
                #                                         geneset_id=self.geneset_id, anno_id=self.anno_id[select_anno],
                #                                         level_id=level_id[select_anno])
                self.api_dic["beta"].add_beta_diversity(path, "dbrda", main="True", anno_type=select_anno,
                                                        env_id=self.env_id, name=submit_location, params=params)

    @time_count
    def export_hclust(self):
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
            "second_level": "",
            "hcluster_method": "average",
            "distance_method": "bray_curties",
        }
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        for select_anno in self.anno_table.keys():
            if select_anno == 'geneset':
                continue
                # select_anno = 'gene'
                # submit_location = "Hcluster_" + select_anno.upper() + "_" + time
            else:
                submit_location = "Hcluster_" + select_anno.upper() + "_" + self.default_level2[select_anno] + "_" + time
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            matrix_path = os.path.join(self.output_dir, "compare/Distance", select_anno,
                                       "bray_curtis_new_abund_table.xls")
            newick_path = os.path.join(self.output_dir, "compare/Hcluster", select_anno, "hcluster.tre")
            # dist_id = self.api_dic["dist"].add_dist_table(matrix_path, main=True, anno_id=self.anno_id[select_anno],
            #                                               name="bray_curtis_" + select_anno, params=params,
            #                                               geneset_id=self.geneset_id)
            dist_id = self.api_dic["dist"].add_dist_table(matrix_path, main=True, name=submit_location, params=params)
            # self.api_dic["hcluster"].add_hcluster_tree(newick_path, main=True, anno_type=select_anno,
            #                                            level_id=level_id[select_anno],
            #                                            params=params, specimen_group=self.specimen_group,
            #                                            geneset_id=self.geneset_id, anno_id=self.anno_id[select_anno],
            #                                            update_dist_id=dist_id)
            self.api_dic["hcluster"].add_hcluster_tree(newick_path, main=True, anno_type=select_anno,
                                                       name=submit_location, params=params, update_dist_id=dist_id)

    @time_count
    def export_cor_heatmap(self):
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
            "second_level": "",
            "top": "50",
            "env_id": str(self.env_id),
            "env_labs": self.env_labs,
            "coefficent": "spearmanr",
            "species_cluster": "",
            "env_cluster": "",
        }
        time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        for select_anno in self.anno_table.keys():
            if select_anno == 'geneset':
                continue
                # select_anno = 'gene'
                # submit_location = "CorrHeatmapSpearman_" + select_anno.upper() + "_" + time
            else:
                submit_location = "CorrHeatmapSpearman_" + select_anno.upper() + "_" + self.default_level2[select_anno] + "_" + time
            params = dict(params, **self.get_param(select_anno))
            params["level_id"] = level_id[select_anno]
            params["submit_location"] = submit_location
            # return_id = self.api_dic["cor_heatmap"].add_heatmap_cor(select_anno, self.geneset_id,
            #                                                         self.anno_id[select_anno],
            #                                                         level_id[select_anno], self.env_id,
            #                                                         self.specimen_group,
            #                                                         name="cor_heatmap_" + select_anno, params=params)
            return_id = self.api_dic["cor_heatmap"].add_heatmap_cor(select_anno, name=submit_location, params=params)
            corr_path = os.path.join(self.output_dir, "correlation/cor_heatmap", select_anno,
                                     "correlation.xls")
            pvalue_path = os.path.join(self.output_dir, "correlation/cor_heatmap", select_anno, "pvalue.xls")
            species_tree = os.path.join(self.anno2correlation_tree[select_anno], "final_species_tree.tre")
            env_tree = os.path.join(self.anno2correlation_tree[select_anno], "final_env_tree.tre")
            self.api_dic["cor_heatmap"].add_heatmap_cor_detail(corr_path, "correlation", return_id, species_tree="",
                                                               env_tree="")
            self.api_dic["cor_heatmap"].add_heatmap_cor_detail(pvalue_path, "pvalue", return_id,
                                                               species_tree=species_tree, env_tree=env_tree)

    def run(self):
        """
        运行 meta_genomic workflow
        :return:
        """
        self.logger.info("在run函数中。。。")
        if self._sheet.id:
            self.logger.info("sheet_id is %s" % self._sheet.id)
        if self._sheet.member_type:
            self.logger.info("member_type is %s" % self._sheet.member_type)
        if self._sheet.cmd_id:
            self.logger.info("cmd_id is %s" % self._sheet.cmd_id)
        task_info = self.api.api('task_info.mg_task_info')
        task_info.add_task_info()
        self.sequence.on('end', self.run_qc)
        self.logger.info("运行了第一部分")
        if self.option('rm_host'):
            self.qc.on('end', self.run_rm_host)
            self.rm_host.on('end', self.run_assem)
        else:
            self.qc.on('end', self.run_assem)
        self.assem_soapdenovo.on('end', self.run_gene_predict)
        self.assem_idba.on('end', self.run_gene_predict)
        self.gene_predict.on('end', self.run_gene_set)
        if self.option('nr'):
            self.gene_set.on('end', self.run_nr)
            self.anno_tool.append(self.nr)
            self.choose_anno.append('nr')
        if self.option('kegg'):
            self.gene_set.on('end', self.run_kegg)
            self.anno_tool.append(self.kegg)
            self.choose_anno.append('kegg')
        if self.option('cog'):
            self.gene_set.on('end', self.run_cog)
            self.anno_tool.append(self.cog)
            self.choose_anno.append('cog')
        if self.option('cazy'):
            self.gene_set.on('end', self.run_cazy)
            self.all_anno.append(self.cazy)
            self.choose_anno.append('cazy')
        if self.option('ardb'):
            self.gene_set.on('end', self.run_ardb)
            self.all_anno.append(self.ardb)
            self.choose_anno.append('ardb')
        if self.option('card'):
            self.gene_set.on('end', self.run_card)
            self.all_anno.append(self.card)
            self.choose_anno.append('card')
        if self.option('vfdb'):
            self.gene_set.on('end', self.run_vfdb)
            self.all_anno.append(self.vfdb)
            self.choose_anno.append('vfdb')
        if len(self.anno_tool) != 0:
            self.on_rely(self.anno_tool, self.run_anno)
            self.all_anno.append(self.anno)
        self.sample_in_group = self.get_sample()
        if len(self.all_anno) == 0:
            self.gene_set.on('end', self.run_analysis, 'composition')
        else:
            if len(self.sample_in_group) < 2:
                self.on_rely(self.all_anno, self.end)
            elif len(self.sample_in_group) == 2:
                self.on_rely(self.all_anno, self.run_analysis, 'composition')
                # self.on_rely(self.analysis, self.end)
            elif len(self.sample_in_group) > 2:
                self.on_rely(self.all_anno, self.run_analysis, 'all')
                # self.on_rely(self.analysis, self.end)
        if self.option('test'):
            # self.run_gene_set()
            # self.run_nr()
            # self.run_kegg()
            # self.run_cog()
            # self.run_ardb()
            # self.run_card()
            # self.run_vfdb()
            # self.run_cazy()
            self.run_analysis('all')
            # self.end()
            super(MetaGenomicWorkflow, self).run()
            return True
            # pass
        if self.option('qc'):
            self.run_sequence()
        elif self.option('rm_host'):
            self.run_rm_host()
        else:
            self.run_assem()
        super(MetaGenomicWorkflow, self).run()
