# -*- coding:utf-8 -*-
# __author__ = 'fwy'

import glob
import json
import os
import re
import shutil
from collections import OrderedDict
import datetime
import pandas as pd
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.wpm.client import worker_client
from biocluster.workflow import Workflow
import unittest
from mbio.packages.ref_rna_v3.gene_fusion.extract_pos_bygtf import extract_all_gene_pos,extract_all_chr_length,gtf_check
from mbio.packages.medical_transcriptome.copy_file import CopyFile
from mbio.packages.medical_transcriptome.extract_annotation_result import RefAnnotation
from mbio.packages.medical_transcriptome.functions import workfuncdeco
# from mbio.packages.medical_transcriptome.upload import set_assembly, set_rmats ,set_gene_fusion
from mbio.packages.medical_transcriptome.upload import Upload
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo
from mbio.packages.medical_transcriptome.functions import tryforgood
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple
from mbio.packages.dna_evolution.send_email import SendEmail
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class MedicalTranscriptomeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MedicalTranscriptomeWorkflow, self).__init__(wsheet_object)
        options = [
            ### 基础参数设置
            ## 分析方式
            #分析策略 ['total', 'mapping'，'quantification']
            {'name': 'analysis_strategy', 'type': 'string', 'default': 'total'},
            # 分析水平 ['Transcript', 'Gene']
            {'name': 'level', 'type': 'string', 'default': 'Transcript'},
            # 大于 %的样本
            {'name': 'mapping_sample_percent', 'type': 'float', 'default': 50.0},
            # Mapping Ratio 小于
            {'name': 'mapping_ratio', 'type': 'float', 'default': 65.0},
            # 大于 %的样本
            {'name': 'rrna_sample_percent', 'type': 'float', 'default': 50.0},
            # rRNA Ratio 小于 %
            {'name': 'rrna_ratio', 'type': 'float', 'default': 15.0},

            ## 样本选择
            # 样本个数 ['multiple', 'single']
            {'name': 'sample_num', 'type': 'string', 'default': 'multiple'},
            # 生物学重复 [True, False]
            {'name': 'is_duplicate', 'type': 'bool', 'default': True},
            # 肿瘤样本 ["yes", "no"]
            {'name': 'is_tumor', 'type': 'string', 'default': None},
            # 分组方案
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            # 上机名称
            {'name': 'productive_table', 'type': 'infile', 'format': 'sample.group_table'},
            # 配对信息表
            {'name': 'pair_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'tumor_pair_table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # 对照组文件
            {'name': 'control_file', 'type': 'infile', 'format': 'sample.control_table'},

            ## 物种基因组信息
            # 具体物种 sg_genome_db.name
            {'name': 'ref_genome', 'type': 'string', 'default': None},
            # 基因组版本 sg_genome_db.assembly
            {'name': 'genome_version', 'type': 'string', 'default': None},
            # 基因组注释版本 sg_genome_db.annot_version
            {'name': 'genome_annot_version', 'type': 'string', 'default': None},
            # 基因组编号 sg_genome_db.genome_id
            {'name': 'genome_id', 'type': 'string', 'default': None},

            ## 测序数据信息
            # 数据类型 ['rawdata', 'cleandata']
            {'name': 'datatype', 'type': 'string', 'default': 'rawdata'},
            # 原始序列文件
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            # 质控序列文件
            # {'name': 'qc_dir', 'type': 'string', 'default': None},
            {'name': 'qc_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            # 链特异性 [False, True]
            {'name': 'strand_specific', 'type': 'bool', 'default': False},
            # 链特异性方向 ['PE': {'RF', 'FR'}, 'SE': {'R', 'F'}]
            {'name': 'strand_dir', 'type': 'string', 'default': 'forward'},
            # 测序质量 ['phred 33', 'phred 64']
            {'name': 'quality_score_system', 'type': 'string', 'default': 'phred 33'},
            # 测序类型 ['PE', 'SE']
            {'name': 'fq_type', 'type': 'string', 'default': 'PE'},
            #接口序列(文库类型)
            {"name": "lib_type", "type": "string", "default": 'type I'},
            # "接头序列"    # --adapter_sequence,the adapter for read1
            {'name': 'adapter_a', 'type': 'string', 'default': 'AGATCGGAAGAGCACACGTC'},
            # "接头序列"    # --adapter_sequence,the adapter for read2 -only in PE
            {'name': 'adapter_b', 'type': 'string', 'default': 'AGATCGGAAGAGCGTCGTGT'},
            {'name': 'adapter_s', 'type': 'string', 'default': 'AGATCGGAAGAGCACACGTC'},

            ### 高级参数设置
            # 质控软件 ['fastp', 'seqprep']
            {'name': 'qc_soft', 'type': 'string', 'default': 'fastp'},
            # 长度需求(质控)[未开放]
            {'name': 'length_required', 'type': 'string', 'default': '30'},
            # 比对软件 ['Hisat', 'Tophat', 'STAR']
            {'name': 'align_method', 'type': 'string', 'default': 'Hisat'},
            # 转录组质量评估[未开放]
            {'name': 'map_assess_method', 'type': 'string', 'default': 'saturation,coverage,distribution,chr_stat'},
            # 是否进行组装 [True, False]
            {'name': 'is_assemble', 'type': 'bool', 'default': True},
            # 拼接软件 ['stringtie', 'cufflinks']
            {'name': 'assemble_method', 'type': 'string', 'default': 'stringtie'},

            ## 转录组功能注释
            #基础数据库注释版本
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_202007"},
            # NR库分类 ['Animal, Plant', 'Protist', 'Fungi', 'All']
            {'name': 'nr_database', 'type': 'string', 'default': "metazoa"},
            # KEGG库分类 ['Animal, Plant', 'Protist', 'Fungi', 'All']
            {'name': 'kegg_database', 'type': 'string', 'default': "Animals"},
            # KEGG二级分类
            {"name": "kegg_subtax1", "type": "string", "default": "Vertebrates"},
            # KEGG三级分类
            {"name": "kegg_subtax2", "type": "string", "default": "Mammals"},
            # KEGG物种
            {"name": "kegg_specific", "type": "string", "default": None},

            #其他注释模块未开放参数[注释相关未开放]
            # NR(GO)-Evalue
            {'name': 'nr_evalue', 'type': 'float', 'default': 1e-5},
            # Swissprot-Evalue
            {'name': 'swissprot_evalue', 'type': 'float', 'default': 1e-5},
            # KEGG-Evalue
            {'name': 'kegg_evalue', 'type': 'float', 'default': 1e-5},
            # COG-Evalue
            {'name': 'cog_evalue', 'type': 'float', 'default': 1e-5},
            # Pfam-Evalue
            {'name': 'pfam_evalue', 'type': 'float', 'default': 1e-5},
            # 其它注释过滤参数
            {'name': 'nr_identity', 'type': 'float', 'default': 0},
            {'name': 'nr_similarity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_identity', 'type': 'float', 'default': 0},
            {'name': 'swissprot_similarity', 'type': 'float', 'default': 0},
            {'name': 'kegg_identity', 'type': 'float', 'default': 0},
            {'name': 'kegg_similarity', 'type': 'float', 'default': 0},
            {"name": "kegg_version", "type": "string", "default": "202003"},
            {'name': 'cog_identity', 'type': 'float', 'default': 0},
            {'name': 'cog_similarity', 'type': 'float', 'default': 0},
            {'name': 'database', 'type': 'string', 'default': 'go,nr,cog,kegg,swissprot,pfam'},

            ##表达量分析
            # 表达定量软件 ['RSEM', 'Kallisto', 'Salmon']
            {'name': 'express_method', 'type': 'string', 'default': 'RSEM'},
            # 表达定量指标 ['tpm', 'fpkm']
            {'name': 'exp_way', 'type': 'string', 'default': 'tpm'},

            ##表达差异分析
            # 差异分析软件 ['DESeq2', 'edgeR', 'DEGseq', 'limma', 'NOIseq']
            {'name': 'diff_method', 'type': 'string', 'default': 'DESeq2'},
            # FC[未开放]
            {'name': 'fc', 'type': 'float', 'default': 2},
            # 显著性指标[未开放]
            {'name': 'pvalue_padjust', 'type': 'string', 'default': 'padjust'},
            # 显著性水平[未开放]
            {'name': 'diff_fdr_ci', 'type': 'float', 'default': 0.05},
            # 多重验证校正方法 ['BH', 'Bonferroni', 'Holm', 'BY'][未开放]
            {'name': 'padjust_way', 'type': 'string', 'default': 'BH'},


            ##基因结构分析
            # SNP分析 ['True', 'False', 'Skip']
            {'name': 'is_snp', 'type': 'string', 'default': 'False'},
            # SNP分析方法 ['Samtools', 'GATK', 'Sentieon']
            {'name': 'snp_method', 'type': 'string', 'default': 'samtools'},
            # 可变剪切分析 ['True', 'False', 'Skip']
            {'name': 'is_as', 'type': 'string', 'default': 'True'},
            # 可变剪接分析方法 ['rmats', 'as_pfofile']
            {'name': 'as_method', 'type': 'string', 'default': 'rmats'},
            # 是否进行somatic分析['True', 'False', 'Skip']
            {'name': 'is_somatic', 'type': 'string', 'default': 'False'},
            # SNP分析方法 ['Samtools', 'GATK', 'Sentieon']
            {'name': 'somatic_method', 'type': 'string', 'default': 'sentieon'},
            # 是否将质控结果和比对结果文件传输至线下服务器
            {'name': 'upload_offline', 'type': 'bool', 'default': False},
            # 是否存入报告图片
            {'name': 'report_img', 'type': 'bool', 'default': False},
            {'name': 'get_run_log', 'type': 'bool', 'default': False},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        if self._sheet.rerun:
            self.logger.info("果然有rerun标记")
        try:
            self.rerun = self._sheet.rerun
            self.logger.info("果然有rerun标记")
        except:
            self.logger.info("我就说没rerun标记吧")
        try:
            if self.rerun:
                self.logger.info("就是重运行的")
            else:
                self.logger.info("不是重运行的")
        except:
            self.logger.info("没标记瞎起劲")

        # self.database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname(mtype='ref_rna_v2', dydb_forbid=True)]
        # collection = self.database['sg_genome_db']
        # genome_info = collection.find_one({'genome_id': self.option('genome_id')})
        sg_genome_db = os.path.join(self.config.SOFTWARE_DIR, 'database/mongo/ref_rna_v2/sg_genome_db.json')
        self.sg_genome_db_list = json.load(open(sg_genome_db, 'r'))

        for genome_db in self.sg_genome_db_list:
            # print(self.option('genome_id'))
            if genome_db['genome_id'] == self.option('genome_id'):
                # print(genome_db)
                genome_info = genome_db
                break
        # else:
            # self.set_error('数据库中不存在该物种的注释信息，程序退出', code='13700320')
        genome_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish')
        print("genome_path :%s", genome_path)
        # print("genome_info :%s", genome_info)
        self.ref_annot_dir = os.path.join(genome_path, genome_info['anno_path_v2'])
        print("ref_annot_dir :%s", self.ref_annot_dir)
        if not os.path.isdir(self.ref_annot_dir):
            self.set_error('数据库中不存在该物种的参考注释文件，程序退出', code='13700321')
        self.ref_genome = os.path.join(genome_path, genome_info['dna_fa'])
        if not os.path.isfile(self.ref_genome):
            self.set_error('数据库中不存在该物种的参考基因组序列FASTA文件，程序退出', code='13700323')
        self.ref_gtf = os.path.join(genome_path, genome_info['gtf'])
        if not os.path.isfile(self.ref_gtf):
            self.set_error('数据库中不存在该物种的参考基因组注释GTF文件，程序退出', code='13700324')
        self.genome_id = genome_info['genome_id']
        if not self.genome_id:
            self.set_error('数据库中不存在该物种的参考基因组ID编号，程序退出', code='13700325')
        self.des = os.path.join(genome_path, genome_info['bio_mart_annot'])
        if not os.path.isfile(self.des):
            self.set_error('数据库中不存在该物种的功能描述信息（BIOMART）', code='13700326')
        self.des_type = genome_info['biomart_gene_annotype']
        if not self.des_type:
            self.set_error('数据库中不存在该物种的功能描述信息（ANNOTYPE）', code='13700327')
        self.known_cds = os.path.join(genome_path, genome_info['cds'])
        if not os.path.isfile(self.known_cds):
            self.set_error('数据库中不存在该物种的CDS序列FASTA文件', code='13700328')
        self.known_pep = os.path.join(genome_path, genome_info['pep'])
        if not os.path.isfile(self.known_pep):
            self.set_error('数据库中不存在该物种的PEP序列FASTA文件', code='13700329')
        self.entrez = os.path.join(genome_path, genome_info['ensemble2entrez'])
        if not os.path.isfile(self.entrez):
            self.set_error('数据库中不存在该物种的ENTREZ对应文件', code='13700330')
        self.genome_stat = os.path.join(genome_path, genome_info['gene_stat'])
        if not os.path.isfile(self.genome_stat):
            self.set_error('数据库中不存在该物种的参考基因组统计文件', code='13700331')
        self.g2t2p = os.path.join(genome_path, genome_info['g2t2p'])
        if not os.path.isfile(self.g2t2p):
            self.set_error('数据库中不存在该物种的G2T2P对应文件', code='13700332')
        self.species = genome_info['taxon_id']
        if self.species == None or not self.species:
            self.species = ""
            # self.set_error('数据库中不存在该物种的TAXON编号', code='13700333')
        self.hyperlink = genome_info['ensemble_web']
        if not self.hyperlink:
            self.set_error('数据库中不存在该物种的参考基因组来源超链接', code='13700334')
        self.known_ko = os.path.join(genome_path, genome_info['kegg'])
        if not os.path.isfile(self.known_ko):
            self.set_error('数据库中不存在该种的PATHWAY对应文件', code='13700335')
        self.genome_version = genome_info['assembly']
        self.annot_version = genome_info['annot_version']
        self.organism_name = genome_info['organism_name']
        ## 定义分析方案
        if self.option("analysis_strategy") == "mapping":
            self.analysis_content = ["mapping"]
        elif self.option("analysis_strategy") == "quantification":
            self.analysis_content = ["mapping", "annotation","quantification"]
        elif self.option("analysis_strategy") == "total":
            self.analysis_content = ["mapping", "annotation", "quantification", "other"]
        else:
            self.analysis_content = ["mapping", "annotation", "quantification", "other"]

        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        # data = os.path.join(self.work_dir, 'data.json')
        # if os.path.exists(data):
        #     with open(data, 'r') as load_f:
        #         load_dict = json.load(load_f)
        #         if 'rerun' in load_dict and load_dict['rerun']:
        #             self.logger.info("该项目重运行中，先删除mongo库中已有数据")
        #             self.delete_mongo_data()

        if self._sheet.rerun:
            pass
            # self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            # self.delete_mongo_data()

        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))

    @tryforgood
    def delete_mongo_data(self):
        # self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        # self.program = os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin/python')
        delete = DeleteDemoMongo(self.task_id, 'medical_transcriptome')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format(self.task_id, 'medical_transcriptome')
        # code = os.system(cmd)
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))

    @workfuncdeco
    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        else:
            self.logger.info('succeed in displaying incoming options')
        if self.option('sample_num') == 'single':
            self.option('is_duplicate', False)
            self.option('is_as', 'False')
        dup = False
        # 判断样本个数选择、分组方案以及is_duplicate是否正确
        if self.option('group_table').is_set:
            group_dict = self.option("group_table").prop['group_dict']
            for g in group_dict.keys():
                if len(group_dict[g]) > 1:
                    dup = True
            if self.option('is_duplicate') != dup:
                self.logger.info(self.option('is_duplicate'))
                self.logger.info(dup)
                raise OptionError('生物学重复参数选择和分组方案不匹配，请检查输入是否有误')
            if self.option('group_table').prop['sample_number'] >= 2:
                if self.option('sample_num') == 'single':
                    raise OptionError('分组方案表中有{}个样本，与样本个数为{}相互矛盾，请重新确认参数.'.format(
                        self.option('group_table').prop['sample_number'], self.option('sample_num')))
        if self.option('sample_num') == 'single':
            self.option('is_duplicate', False)
            self.option('is_as', 'False')
        else:
            if self.option("analysis_strategy") in ["quantification", "total"] and not self.option('group_table').is_set:
                raise OptionError('必须上传样本分组表')
            elif not self.option('group_table').is_set:
                samples = list()
                if self.option('datatype') == 'rawdata':
                    fastq_dir = self.option('fastq_dir').path
                else:
                    fastq_dir = self.option('qc_dir').path
                for line in open(os.path.join(fastq_dir, 'list.txt')):
                    sample = line.strip().split('\t')[1]
                    if sample not in samples:
                        samples.append(sample)
                group_table = os.path.join(self.work_dir, 'group.txt')
                with open(group_table, 'w') as w:
                    w.write('#sample\tgroup\n')
                    for sample in samples:
                        w.write('{}\t{}\n'.format(sample, sample))
                self.option("group_table", group_table)
        if not self.option('fq_type') in ['PE', 'SE']:
            raise OptionError('测序类型应为PE或SE', code='13700309')
        if not 0 < self.option('nr_evalue') < 1:
            raise OptionError('NR比对的E值超出范围', code='13700311')
        if not 0 < self.option('cog_evalue') < 1:
            raise OptionError('COG比对的E值超出范围', code='13700312')
        if not 0 < self.option('kegg_evalue') < 1:
            raise OptionError('KEGG比对的E值超出范围', code='13700313')
        if not 0 < self.option('swissprot_evalue') < 1:
            raise OptionError('SWISSPROT比对的E值超出范围', code='13700314')
        if not 0 < self.option('pfam_evalue') < 1:
            raise OptionError('PFAM比对的E值超出范围', code='13700315')
        if not self.option('align_method').lower() in ['tophat', 'hisat', 'star']:
            raise OptionError('比对软件应为Tophat或Hisat或STAR', code='13700316')
        for i in self.option('map_assess_method').split(','):
            if i.lower() not in ['saturation', 'distribution', 'coverage', 'chr_stat']:
                raise OptionError('比对质量评估分析没有%s，请检查', variables=(i), code='13700317')
        if self.option('is_assemble'):
            if self.option('assemble_method').lower() not in ['cufflinks', 'stringtie']:
                raise OptionError('拼接软件应在cufflinks和stringtie中选择', code='13700318')
        if self.option('nr_database') == 'All':
            self.option('nr_database', 'nr')
        elif self.option('nr_database') == 'Animal':
            self.option('nr_database', 'metazoa')
        elif self.option('nr_database') == 'Plant':
            self.option('nr_database', 'viridiplantae')
        elif self.option('nr_database') == 'Protist':
            self.option('nr_database', 'protist')
        elif self.option('nr_database') == 'Fungi':
            self.option('nr_database', 'fungi')
        if self.option('kegg_database') == 'All':
            self.option('kegg_database', 'All')
        elif self.option('kegg_database') == 'Animal':
            self.option('kegg_database', 'Animals')
        elif self.option('kegg_database') == 'Plant':
            self.option('kegg_database', 'Plants')
        elif self.option('kegg_database') == 'Fungi':
            self.option('kegg_database', 'Fungi')
        elif self.option('kegg_database') == 'Protist':
            self.option('kegg_database', 'Protists')
        if not isinstance(self.option('rrna_sample_percent'), float):
                raise OptionError('核糖体判断中止条件的样本比例值应为浮点数', code='13700319')
        if not isinstance(self.option('rrna_ratio'), float):
                raise OptionError('核糖体判断中止条件的阈值比例值应为浮点数', code='13700320')
        if not isinstance(self.option('mapping_sample_percent'), float):
                raise OptionError('比对结果判断中止条件的样本比例值应为浮点数', code='13700321')
        if not isinstance(self.option('mapping_ratio'), float):
                raise OptionError('比对结果判断中止条件的阈值比例值应为浮点数', code='13700322')
        if self.option('pair_table').is_set:
            if self.option('diff_method').lower() not in ['deseq2', 'edger', 'limma']:
                raise OptionError('做配对分析时，需要从DESeq2, edgeR, limma三款软件中选择')
        if self.option('sample_num') == 'multiple':
            group_size = list()
            group_dict = self.option("group_table").prop['group_dict']
            for key in group_dict:
                group_size.append(len(group_dict[key]))
            group_size.sort()
            if group_size[0] == group_size[-1] == 1:
                if self.option("diff_method").lower() == "deseq2":
                    self.option("diff_method", "DEGseq")
                    self.option("diff_fdr_ci", 0.001)
                    self.logger.info("该项目没有生物学重复,不可以使用DESeq2,修改方法为DEGseq,阈值设置为0.001")
            elif group_size[0] == 1 and group_size[-1] >= 2:
                if self.option("diff_method").lower() == "deseq2":
                    self.option("diff_method", "edgeR")
                    self.option("diff_fdr_ci", 0.05)
                    self.logger.info("该项目部分组别有生物学重复,部分组别无生物学重复,不可以使用DESeq2,修改方法为edgeR,阈值设置为0.05")
            elif group_size[0] >= 2 and group_size[-1] >= 2:
                if self.option("diff_method").lower() == "degseq":
                    self.option("diff_method", "DESeq2")
                    self.option("diff_fdr_ci", 0.05)
                    self.logger.info("该项目有生物学重复,不可以使用DEGseq,修改方法为DESeq2,阈值设置为0.05")
            else:
                pass
        # if self.option("diff_method").lower() == "degseq":
        #     self.option("diff_fdr_ci",0.001)

        # if self.option("is_somatic") == "True":
        #     if not self.option('pair_table').is_set:
        #         raise OptionError('肿瘤样本需上传配对分析表')
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        else:
            return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @workfuncdeco
    def run(self):
        # 
        self.add_steps()

        print("start su")
        super(MedicalTranscriptomeWorkflow, self).run()
       


    @workfuncdeco
    def add_steps(self):
        if self.option('datatype') == 'rawdata':
            if self.option('qc_soft') == 'fastp':
                self.step.add_steps('fastp_rna')
            elif self.option('qc_soft') == 'seqprep':
                self.step.add_steps('hiseq_qc')
        if "mapping" in self.analysis_content:
            self.step.add_steps('file_check', 'hiseq_reads_stat_raw', 'hiseq_reads_stat_use', 'rnaseq_mapping',
                                'map_assessment')
        if "quantification" in self.analysis_content:
            self.step.add_steps('quant')
            if self.option('sample_num') == 'multiple':
                group_spname = self.option("group_table").get_group_spname()
                group_snum = [len(group_spname[g]) for g in group_spname]
                min_group_num = min(group_snum)
                if self.option('group_table').prop['sample_number'] > 2:
                    self.step.add_steps('exp_pca')
                    self.step.add_steps('exp_corr')
                if len(self.option("group_table").prop['group_dict']) > 1:
                    self.step.add_steps('exp_venn')
                if min_group_num >= 3:
                    self.step.add_steps('ellipse')
        if "annotation" in self.analysis_content:
            # self.step.add_steps('annot_filter_ref', 'annot_class_beta_ref', 'gene_fa', 'detail')
            self.step.add_steps('annot_filter_ref','gene_fa', 'detail')
            if self.option('is_assemble'):
                self.step.add_steps('refrna_assemble', 'annot_mapdb', 'annot_orfpfam', 'annot_filter_new',
                                    'annot_class_beta_new', 'annot_merge')
            else:
                self.step.add_steps('transcript_abstract')
        if "other" in self.analysis_content:
            if self.option('sample_num') == 'multiple':
                self.step.add_steps('diffexp')
                self.step.add_steps('diff_geneset_analysis')
            if self.option('is_snp') == 'True':
                if self.option('snp_method').lower() == 'samtools':
                    self.step.add_steps('sam_rna')
                elif self.option('snp_method').lower() == 'gatk':
                    self.step.add_steps('snp_rna')
                elif self.option('snp_method').lower() == 'sentieon':
                    self.step.add_steps('call_snp_indel')
            if self.option('is_somatic') == 'True':
                if self.option('somatic_method').lower() == 'samtools':
                    self.step.add_steps('somatic_sam')
                elif self.option('somatic_method').lower() == 'gatk':
                    self.step.add_steps('somatic_gatk')
                elif self.option('somatic_method').lower() == 'sentieon':
                    self.step.add_steps('somatic_sentieon')
            if self.option('is_as') == 'True':
                self.step.add_steps('rmats')

        self.prepare_module_tools()

    @workfuncdeco
    def prepare_module_tools(self):
        self.file_check = self.add_tool('medical_transcriptome.file_check')
        self.hiseq_reads_stat_use = self.add_module('medical_transcriptome.hiseq_reads_stat')
        self.rely = [self.hiseq_reads_stat_use]
        if self.option('datatype') == 'rawdata':
            self.hiseq_reads_stat_raw = self.add_module('medical_transcriptome.hiseq_reads_stat')
            self.on_rely([self.hiseq_reads_stat_raw, self.hiseq_reads_stat_use], self.check_rrna)
            self.rely.append(self.hiseq_reads_stat_raw)
            #根据选择的质控软件决定Module
            if self.option('qc_soft') == 'fastp':
                self.fastp_rna = self.add_module('datasplit.fastp_rna')
            elif self.option('qc_soft') == 'seqprep':
                self.hiseq_qc = self.add_module('medical_transcriptome.hiseq_qc')
        self.rnaseq_mapping = self.add_module('medical_transcriptome.rnaseq_mapping')
        self.rely.append(self.rnaseq_mapping)
        if "mapping" in self.analysis_content:
            self.map_assessment = self.add_module('medical_transcriptome.map_assessment')
            self.rely.append(self.map_assessment)
        if "quantification" in self.analysis_content:
            self.quant = self.add_module('medical_transcriptome.quant.quant')
            self.rely.append(self.quant)
            if self.option('sample_num') == 'multiple':
                group_spname = self.option("group_table").get_group_spname()
                group_snum = [len(group_spname[g]) for g in group_spname]
                min_group_num = min(group_snum)
                group_dict = self.option("group_table").prop['group_dict']
                if len(group_dict) > 1:
                    self.exp_venn = self.add_tool('medical_transcriptome.exp_venn')
                    self.rely.append(self.exp_venn)
                if self.option('group_table').prop['sample_number'] > 2:
                    self.exp_pca = self.add_tool('medical_transcriptome.exp_pca')
                    self.exp_corr = self.add_tool('medical_transcriptome.exp_corr')
                    self.rely.extend([self.exp_corr, self.exp_pca])
                if min_group_num >= 3:
                    self.ellipse = self.add_tool('medical_transcriptome.ellipse')
                    self.rely.append(self.ellipse)
        if "annotation" in self.analysis_content:
            # self.annot_filter_ref = self.add_module('ref_rna_v2.annot_filter')
            # self.annot_class_beta_ref = self.add_module('ref_rna_v2.annot_class_beta')
            self.gene_fa = self.add_tool('medical_transcriptome.gene_fa')
            self.detail = self.add_tool('medical_transcriptome.database.detail')
            self.annot_merge = self.add_tool('medical_transcriptome.annotation.annot_merge')
            self.rely.append(self.annot_merge)
            if self.option('is_assemble'):
                self.refrna_assemble = self.add_module('medical_transcriptome.refrna_assemble')
                self.annot_mapdb = self.add_module('medical_transcriptome.annot_mapdb')
                self.annot_orfpfam = self.add_module('medical_transcriptome.annot_orfpfam')
                self.annot_filter_new = self.add_module('medical_transcriptome.annot_filter')
                self.annot_class_beta_new = self.add_module('medical_transcriptome.annot_class_beta')
                self.rely.append(self.refrna_assemble)
                self.on_rely([self.annot_mapdb, self.annot_orfpfam], self.run_annot_class_beta_new)
                self.on_rely([self.gene_fa, self.annot_orfpfam], self.run_detail)
                self.on_rely([self.annot_class_beta_new], self.run_annot_merge)
                self.rely.append(self.annot_merge)
                # self.run_refrna_assemble()annot_class_medical
            else:
                self.transcript_abstract = self.add_tool('medical_transcriptome.transcript_abstract')
                # self.run_transcript_abstract()
                self.on_rely([ self.detail], self.run_annot_merge)
        if "other" in self.analysis_content:
            if self.option('sample_num') == 'multiple':
                self.diffexp = self.add_tool('medical_transcriptome.batch.diffexp_batch')
                self.diff_geneset_analysis = self.add_module("medical_transcriptome.workflow_diffgt.diff_geneset_all_pipline")
                self.rely.append(self.diffexp)
                self.rely.append(self.diff_geneset_analysis)
                self.on_rely([self.annot_merge,self.quant,self.diffexp], self.run_diff_geneset_analysis)
            if self.option('is_snp') == 'True':
                if self.option('snp_method').lower() == 'samtools':
                    self.sam_rna = self.add_module('medical_transcriptome.snp.call_snp_indel_samtools')
                    self.rely.append(self.sam_rna)
                elif self.option('snp_method').lower() == 'gatk':
                    self.snp_rna = self.add_module('medical_transcriptome.snp.call_snp_indel_gatk')
                    self.rely.append(self.snp_rna)
                elif self.option('snp_method').lower() == 'sentieon':
                    self.call_snp_indel = self.add_module('medical_transcriptome.snp.call_snp_indel_sentieon')
                    self.rely.append(self.call_snp_indel)
            self.gene_fusion = self.add_module("medical_transcriptome.gene_fusion")
            if  self.option('fq_type') == 'PE':
                self.rely.append(self.gene_fusion)
            if self.option('is_as') == 'True':
                if self.option("as_method").lower() == "as_profile":
                    self.as_analysis = self.add_module("medical_transcriptome.asprofile.asprofile")
                else:
                    self.as_analysis = self.add_module('medical_transcriptome.rmats.rmats')
                self.rely.append(self.as_analysis)
            if self.option("is_somatic") == "True":
                self.somatic_analysis = self.add_module('medical_transcriptome.somatic.call_somatic_sentieon')
                self.rely.append(self.somatic_analysis)
        self.logger.info("analysis_content: {}".format(self.analysis_content))
        self.logger.info("self.rely: {}".format(self.rely))
        self.on_rely(self.rely, self.run_chart)
        # self.on_rely(self.rely,self.set_upload)
        self.run_file_check()
        # self.set_upload()
        # self.set_db()

    def run_chart(self):
        '''
        绘图步骤插入在导表前
        '''
        self.chart = self.add_tool('medical_transcriptome.chart')

        # group_dict = self.option("group_table").prop["group_dict"]
        if not self.option('group_table').is_set:
            samples = list()
            if self.option('datatype') == 'rawdata':
                fastq_dir = self.option('fastq_dir').path
            else:
                fastq_dir = self.option('qc_dir').path
            for line in open(os.path.join(fastq_dir, 'list.txt')):
                sample = line.strip().split('\t')[1]
                if sample not in samples:
                    samples.append(sample)
            group_table = os.path.join(self.work_dir, 'group.txt')
            with open(group_table, 'w') as w:
                w.write('#sample\tgroup\n')
                for sample in samples:
                    w.write('{}\t{}\n'.format(sample, sample))
            self.option("group_table", group_table)
        samples = self.option("group_table").prop["sample"]
        # cmp_list = self.option("control_file").prop["cmp_list"]
        chart_dict = {
            "type": "workflow",
            "samples": samples,
            "qc_file_raw": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                table_dir=self.hiseq_reads_stat_raw.output_dir, sample_name='{sample_name}'),
            "qc_file_use": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                table_dir=self.hiseq_reads_stat_use.output_dir, sample_name='{sample_name}'),
            "align_satu_r": "{table_dir}/saturation/satur_{sample_name}.eRPKM.xls.saturation.R".format(
                table_dir=self.map_assessment.output_dir, sample_name='{sample_name}'),
            "align_satu_p": "{table_dir}/saturation/satur_{sample_name}.eRPKM.xls.cluster_percent.xls".format(
                table_dir=self.map_assessment.output_dir, sample_name='{sample_name}'),
            "align_coverage": "{table_dir}/coverage/{sample_name}.geneBodyCoverage.txt".format(
                table_dir=self.map_assessment.output_dir, sample_name='{sample_name}'),
            "align_pos": "{table_dir}/distribution/{sample_name}.reads_distribution.txt".format(
                table_dir=self.map_assessment.output_dir, sample_name='{sample_name}'),
            "align_chr": "{table_dir}/chr_stat/{sample_name}.bam_chr_stat.xls".format(
                table_dir=self.map_assessment.output_dir, sample_name='{sample_name}')

        }
        if self.option('fq_type') == "SE":
            chart_dict.update({
                "qc_file_raw": "{table_dir}/qualityStat/{sample_name}.qual_stat".format(
                    table_dir=self.hiseq_reads_stat_raw.output_dir, sample_name='{sample_name}'),
                "qc_file_use": "{table_dir}/qualityStat/{sample_name}.qual_stat".format(
                    table_dir=self.hiseq_reads_stat_use.output_dir, sample_name='{sample_name}')

            })
        if self.option("group_table").is_set:
            group_dict = self.option("group_table").prop["group_dict"]
            chart_dict.update({
                "group_dict": group_dict
            })
        if self.option("control_file").is_set:
            cmp_list = self.option("control_file").prop["cmp_list"]
            chart_dict.update({
                "cmp_list": cmp_list
            })


        if "other" in self.analysis_content and self.option('sample_num') == 'multiple':
            chart_dict.update({
                "diff_exp_summary": "{table_dir}/diff_summary_{soft}.xls".format(
                    table_dir=self.diffexp.output_dir, soft=self.option("diff_method").lower()),
                "diff_exp": "{table_dir}/{control}_vs_{test}.{soft}.xls".format(
                    table_dir=self.diffexp.output_dir, control='{control}', test='{test}', soft=self.option("diff_method").lower())
            })

        if "other" in self.analysis_content:
            # 基因融合
            if self.option('fq_type') == 'PE':
                sample_names = os.listdir(os.path.join(self.gene_fusion.output_dir, "star_fusion"))
                extract_all_gene_pos(self.ref_gtf, os.path.join(self.gene_fusion.work_dir, "gene_pos"))
                assemble_level_file = os.path.join(os.path.dirname(self.ref_gtf), "assembly_level.txt")
                chr_length_path = os.path.join(os.path.dirname(self.ref_genome), "star_27",
                                               "ctat_genome_lib_build_dir", "ref_genome.fa.star.idx", "chrNameLength.txt")
                extract_all_chr_length(assemble_level_file, chr_length_path,
                                       os.path.join(self.gene_fusion.work_dir, "chr_length"))
                genefusion = {
                    "result_dir": self.gene_fusion.output_dir + "/star_fusion/{sample_name}/star-fusion.fusion_predictions.abridged.tsv",
                    "chr_length_path": os.path.join(self.gene_fusion.work_dir, "chr_length"),
                    "gene_pos": os.path.join(self.gene_fusion.work_dir, "gene_pos"),
                    "samples": sample_names}
                chart_dict.update({"genefusion":genefusion})
        # if "annotation" in self.analysis_content:
        #     chart_dict.update({
        #         "annot_stat": "{table_dir}/allannot_class/all_stat.xls".format(
        #             table_dir=self.annot_merge.output_dir)
        #     })
        #     if self.option('is_assemble'):
        #         chart_dict.update({
        #             "assemble_step": "{table_dir}/Statistics/trans_count_stat_200.txt".format(
        #                 table_dir=self.refrna_assemble.output_dir),
        #             "assemble_stat": "{table_dir}/Statistics/".format(
        #                 table_dir=self.refrna_assemble.output_dir),
        #             "assemble_new": "{table_dir}/Statistics/code_num.txt".format(
        #                 table_dir=self.refrna_assemble.output_dir)
        #         })
        if "quantification" in self.analysis_content:
            chart_dict.update({
                "gene_exp_ref": "{table_dir}/ref.gene.tpm.matrix".format(
                    table_dir=self.quant.work_dir),
                "tran_exp_ref": "{table_dir}/ref.transcript.tpm.matrix".format(
                    table_dir=self.quant.work_dir),
                "gene_exp_all": "{table_dir}/gene.tpm.matrix".format(
                    table_dir=self.quant.work_dir),
                "tran_exp_all": "{table_dir}/transcript.tpm.matrix".format(
                    table_dir=self.quant.work_dir)
            })

            if self.option('sample_num') == 'multiple':
                if self.option('group_table').prop['sample_number'] > 2:
                    chart_dict.update({
                        "exp_pca_file": "{table_dir}/PCA.xls".format(
                            table_dir=self.exp_pca.output_dir),
                        "exp_pca_var_file": "{table_dir}/Explained_variance_ratio.xls".format(
                            table_dir=self.exp_pca.output_dir),
                        "exp_corr_file": "{table_dir}/sample_correlation.xls".format(
                            table_dir=self.exp_corr.work_dir),
                        "exp_corr_tree_file": "{table_dir}/sample.cluster_tree.txt".format(
                            table_dir=self.exp_corr.work_dir)
                })
                if len(self.option("group_table").prop['group_dict']) > 1:
                    chart_dict.update({
                        "exp_venn" : self.exp_venn.output_dir + "/venn_graph.xls"
                    })

        if hasattr(self, 'ellipse'):
            chart_dict.update({
                "exp_pca_ellipse": "{table_dir}/ellipse_out.xls".format(table_dir=self.ellipse.output_dir)
            })
        if self.option('is_somatic') == 'True' and 'other' in self.analysis_content:
            somatic_result = self.somatic_analysis.output_dir
            chart_dict.update({
                "somatic_distribution": "{table_dir}/predeal/snp_position_distribution.xls".format(table_dir=somatic_result),
                "somatic_stat": "{table_dir}/predeal/snp_transition_tranversion_statistics.xls".format(table_dir=somatic_result),
                "somatic_depth": "{table_dir}/predeal/snp_depth_statistics.xls".format(table_dir=somatic_result)
            })
        if self.option('is_snp') == 'True' and 'other' in self.analysis_content:
            if self.option('snp_method').lower() == 'samtools':
                snp_result = self.sam_rna.output_dir
            elif self.option('snp_method').lower() == 'gatk':
                snp_result = self.snp_rna.output_dir
            elif self.option('snp_method').lower() == 'sentieon':
                snp_result = os.path.join(self.call_snp_indel.output_dir,"predeal")
            chart_dict.update({
                "snp_distribution": "{table_dir}/snp_position_distribution.xls".format(table_dir=snp_result),
                "snp_stat" : "{table_dir}/snp_transition_tranversion_statistics.xls".format(table_dir=snp_result),
                "snp_depth" : "{table_dir}/snp_depth_statistics.xls".format(table_dir=snp_result)
            })
        if self.option('is_as') == 'True' and 'other' in self.analysis_content:
            if self.option("as_method").lower() == "rmats":
                chart_dict.update({
                    "splice_stat": "{table_dir}/sample.event.count.{splicestat}.txt".format(table_dir=self.as_analysis.output_dir, splicestat='{splicestat}'),
                    "splice_diff": "{table_dir}/{control}_vs_{test}/event_stats.file.txt".format(table_dir=self.as_analysis.output_dir, control='{control}', test='{test}'),
                    "splice_psi": "{table_dir}/{control}_vs_{test}/psi_stats.file.txt".format(table_dir=self.as_analysis.output_dir, control='{control}', test='{test}')
                })
            else:
                chart_dict.update({
                    "as_profile_stat": "{table_dir}/AS_statistics_merge.txt".format(
                        table_dir=self.as_analysis.output_dir),

                })

        if 'other' in self.analysis_content and self.option('sample_num') == 'multiple':
            genesets = sorted([os.path.basename(i) for i in glob.glob("{}/*vs*".format(self.diff_geneset_analysis.output_dir))])
            chart_dict.update({"genesets":genesets})
            if len(genesets) > 1:
                chart_dict.update({"diff_geneset_venn":self.diff_geneset_analysis.output_dir})
            cluster_geneset_name = os.listdir(os.path.join(self.diff_geneset_analysis.output_dir,"cluster"))[0]
            cluster_dir = os.path.join(self.diff_geneset_analysis.output_dir,"cluster",cluster_geneset_name)
            chart_dict.update({
                "cluster_geneset_name":cluster_geneset_name,
                "cluster_exp" : os.path.join(cluster_dir,"expression_matrix.xls"),
                "cluster_tree" : os.path.join(cluster_dir,"seq.cluster_tree.txt"),
                "sample_tree" : os.path.join(cluster_dir, "sample.cluster_tree.txt"),
                "subcluster_list" : glob.glob(cluster_dir + "/*subcluster_*.xls"),
                "gene_annot_file" :  os.path.join(self.annot_merge.output_dir, "allannot_class", "all_annot_tran.xls")
            })
            file_json_path = os.path.join(self.diff_geneset_analysis.file_prepare.output_dir, "prepare_json")
            with open(file_json_path, "r") as j:
                file_dict = json.load(j)
            kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
            chart_dict.update({
                "go_class":"{table_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(table_dir = self.diff_geneset_analysis.output_dir,geneset_name = "{geneset_name}"),
                "go_enrich": "{table_dir}/{geneset_name}/diff_go_enrich/go_enrich_geneset_list_gene.xls".format(
                table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                "kegg_class": "{table_dir}/{geneset_name}/diff_kegg_class/kegg_stat.xls".format(
                        table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                "kegg_enrich": "{table_dir}/{geneset_name}/diff_kegg_enrich/enrich/{geneset_name}_gene.list.DE.list.check.kegg_enrichment.xls".format(
                            table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                "kegg_level":kegg_level_path,
                "reactome_class": "{table_dir}/{geneset_name}/diff_reactome_class/reactome_class.xls".format(
                                   table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                "reactome_enrich": "{table_dir}/{geneset_name}/diff_reactome_enrich/{geneset_name}_gene.list.reactome_enrichment.xls".format(
                                   table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}")
            })
            if self.organism_name == "Homo_sapiens":
                chart_dict.update({
                    "do_class": "{table_dir}/{geneset_name}/diff_do_class/do_level2_class.xls".format(
                        table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}"),
                    "do_enrich": "{table_dir}/{geneset_name}/diff_do_enrich/do_enrichment.xls".format(
                        table_dir=self.diff_geneset_analysis.output_dir, geneset_name="{geneset_name}")
                })
        with open(self.work_dir + "/chart_workflow.json", 'w') as json_f:
                json.dump(chart_dict, json_f, sort_keys=True, indent=4)
        self.chart.set_options({
            "file_json": self.work_dir + "/chart_workflow.json"
        })
        self.chart.on('end', self.set_db)
        self.chart.run()

    def run_file_check(self):
        if self.option('datatype') == 'rawdata':
            fastq_dir = self.option('fastq_dir')
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir')
        opts = {
            'fq_type': self.option('fq_type'),
            'fastq_dir': fastq_dir,
            'in_gtf': self.ref_gtf,
            'sample_num': self.option('sample_num'),
            'is_duplicate': self.option('is_duplicate')
        }
        if self.option('control_file').is_set:
            opts.update({'control_file': self.option('control_file')})
        if self.option('group_table').is_set:
            opts.update({'group_table': self.option('group_table')})
        self.file_check.set_options(opts)
        self.file_check.on('start', self.set_step, {'start': self.step.file_check})
        self.file_check.on('end', self.set_step, {'end': self.step.file_check})
        if self.option('datatype') == 'rawdata':
            self.file_check.on('end', self.run_hiseq_reads_stat_raw)
            if self.option('qc_soft') == 'fastp':
                self.file_check.on('end', self.run_fastp_rna)
            elif self.option('qc_soft') == 'seqprep':
                self.file_check.on('end', self.run_hiseq_qc)
        elif self.option('datatype') == 'cleandata':
            self.file_check.on('end', self.run_hiseq_reads_stat_use)
        self.file_check.run()

    def run_hiseq_reads_stat_raw(self):
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        self.hiseq_reads_stat_raw.set_options({
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'quality': quality
        })
        self.hiseq_reads_stat_raw.on('start', self.set_step, {'start': self.step.hiseq_reads_stat_raw})
        self.hiseq_reads_stat_raw.on('end', self.set_step, {'end': self.step.hiseq_reads_stat_raw})
        self.hiseq_reads_stat_raw.on('end', self.set_output, 'hiseq_reads_stat_raw')
        self.hiseq_reads_stat_raw.run()

    def run_fastp_rna(self):
        fq_dir = self.option('fastq_dir').path
        sample_path = os.path.join(fq_dir, 'abs.list.txt')
        open(sample_path, 'w').writelines(
            '{}/{}'.format(fq_dir, line) for line in open(os.path.join(fq_dir, 'list.txt'))
        )
        if self.option('fq_type') == "PE":
            self.fastp_rna.set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'length_required': self.option('length_required'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_sequence': self.option('adapter_a'),
                'adapter_sequence_r2': self.option('adapter_b')
            })
        else:
            self.fastp_rna.set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'length_required': self.option('length_required'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_sequence_s': self.option('adapter_a'),
            })
        self.fastp_rna.on('start', self.set_step, {'start': self.step.fastp_rna})
        self.fastp_rna.on('end', self.set_step, {'end': self.step.fastp_rna})
        self.fastp_rna.on('end', self.set_output, 'fastp_rna')
        self.fastp_rna.on('end', self.run_hiseq_reads_stat_use)
        self.fastp_rna.run()

    def run_hiseq_qc(self):
        if self.option('fq_type') == "PE":
            self.hiseq_qc.set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_a': self.option('adapter_a'),
                'adapter_b': self.option('adapter_b')
            })
        else:
            self.hiseq_qc.set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_a': self.option('adapter_a'),
            })
        self.hiseq_qc.on('start', self.set_step, {'start': self.step.hiseq_qc})
        self.hiseq_qc.on('end', self.set_step, {'end': self.step.hiseq_qc})
        self.hiseq_qc.on('end', self.set_output, 'hiseq_qc')
        self.hiseq_qc.on('end', self.run_hiseq_reads_stat_use)
        self.hiseq_qc.run()

    def run_hiseq_reads_stat_use(self):
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        if self.option('datatype') == 'rawdata':
            if self.option('qc_soft') == 'fastp':
                fastq_dir = self.fastp_rna.option('sickle_dir')
            elif self.option('qc_soft') == 'seqprep':
                fastq_dir = self.hiseq_qc.option('sickle_dir')
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir')
        self.hiseq_reads_stat_use.set_options({
            'fastq_dir': fastq_dir,
            'fq_type': self.option('fq_type'),
            'quality': quality,
            'dup': True,
            'rfam': True,
            'rfam_version': self.annot_config_dict['rfam']['version'],
        })
        self.hiseq_reads_stat_use.on('start', self.set_step, {'start': self.step.hiseq_reads_stat_use})
        self.hiseq_reads_stat_use.on('end', self.set_step, {'end': self.step.hiseq_reads_stat_use})
        self.hiseq_reads_stat_use.on('end', self.set_output, 'hiseq_reads_stat_use')
        if self.option('datatype') == 'cleandata':
            self.hiseq_reads_stat_use.on('end', self.check_rrna)
        self.hiseq_reads_stat_use.run()

    def check_rrna(self):
        df = pd.read_table(os.path.join(self.hiseq_reads_stat_use.output_dir, 'stat_results'))
        rrna_sample_percent = len(
            [i for i in df.iloc[:, -1] if i > self.option('rrna_ratio')]
        ) / float(df.shape[0]) * 100
        if rrna_sample_percent > self.option('rrna_sample_percent'):
            self.stop('rrna')
        else:
            self.run_rnaseq_mapping()

    def run_rnaseq_mapping(self):
        if self.option('datatype') == 'rawdata':
            if self.option('qc_soft') == 'fastp':
                fastq_dir = self.fastp_rna.option('sickle_dir')
            elif self.option('qc_soft') == 'seqprep':
                fastq_dir = self.hiseq_qc.option('sickle_dir')
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir')
        self.rnaseq_mapping.set_options({
            'ref_genome': self.option('ref_genome'),
            'genome_version': self.genome_version,
            'genome_annot_version': self.annot_version,
            'mapping_method': self.option('align_method').lower(),
            'seq_method': self.option('fq_type'),
            'fastq_dir': fastq_dir,
            'assemble_method': self.option('assemble_method').lower(),
            'strand_specific': self.option('strand_specific')
        })
        self.rnaseq_mapping.on('start', self.set_step, {'start': self.step.rnaseq_mapping})
        self.rnaseq_mapping.on('end', self.set_step, {'end': self.step.rnaseq_mapping})
        self.rnaseq_mapping.on('end', self.set_output, 'rnaseq_mapping')
        self.rnaseq_mapping.on('end', self.check_mapping)
        self.rnaseq_mapping.on('end', self.run_map_assessment)
        self.rnaseq_mapping.run()

    def check_mapping(self):
        if self.option('align_method').lower() == "hisat":
            n = 0.0
            for i, stat in enumerate(glob.glob(os.path.join(self.rnaseq_mapping.output_dir, 'stat/*.stat'))):
                if float(open(stat).readlines()[-2].split('%')[0]) < float(self.option('mapping_ratio')):
                    n += 1
            else:
                if n / (i + 1) * 100 > float(self.option('mapping_sample_percent')):
                    self.stop('mapping')
                else:
                    if "annotation" in self.analysis_content:
                        if self.option('is_assemble'):
                            self.run_refrna_assemble()
                        else:
                            self.run_transcript_abstract()
                    else:
                        pass
        elif self.option('align_method').lower() == "tophat":
            n = 0.0
            for i, stat in enumerate(glob.glob(os.path.join(self.rnaseq_mapping.output_dir, 'stat/*.stat'))):
                if float(open(stat).readlines()[8].split('%')[0]) < float(self.option('mapping_ratio')):
                    n += 1
            else:
                if n / (i + 1) * 100 > float(self.option('mapping_sample_percent')):
                    self.stop('mapping')
                else:
                    if "annotation" in self.analysis_content:
                        if self.option('is_assemble'):
                            self.run_refrna_assemble()
                        else:
                            self.run_transcript_abstract()
                    else:
                        pass
        elif self.option('align_method').lower() == "star":
            n = 0.0
            i = 0
            stat = os.path.join(self.rnaseq_mapping.output_dir, 'Comparison_results')
            with open(stat, "r") as f1:
                header = f1.readline()
                for line in f1:
                    i += 1
                    if float(line.split("(")[1].split('%')[0]) < float(self.option('mapping_ratio')):
                        n += 1
                else:
                    if n / i * 100 > float(self.option('mapping_sample_percent')):
                        self.stop('mapping')
                    else:
                        if "annotation" in self.analysis_content:
                            if self.option('is_assemble'):
                                self.run_refrna_assemble()
                            else:
                                self.run_transcript_abstract()
                        else:
                            pass


    def run_map_assessment(self):
        self.map_assessment.set_options({
            'bam': self.rnaseq_mapping.option('bam_output'),
            'bed': self.file_check.option('bed'),
            'analysis': self.option('map_assess_method')
        })
        self.map_assessment.on('start', self.set_step, {'start': self.step.map_assessment})
        self.map_assessment.on('end', self.set_step, {'end': self.step.map_assessment})
        self.map_assessment.on('end', self.set_output, 'map_assessment')
        self.map_assessment.run()


    def run_refrna_assemble(self):
        self.step.add_steps('refrna_assemble')
        if self.option('strand_specific'):
            fr_stranded = 'fr-stranded'
            if self.option('strand_dir').startswith('R'):
                strand_direct = 'firststrand'
            elif self.option('strand_dir').startswith('F'):
                strand_direct = 'secondstrand'
        else:
            fr_stranded = 'fr-unstranded'
            strand_direct = 'none'
        self.refrna_assemble.set_options({
            'sample_bam_dir': self.rnaseq_mapping.option('bam_output'),
            'assemble_method': self.option('assemble_method'),
            'ref_gtf': self.file_check.option('gtf'),
            'ref_fa': self.ref_genome,
            'fr_stranded': fr_stranded,
            'strand_direct': strand_direct
        })
        self.refrna_assemble.on('start', self.set_step, {'start': self.step.refrna_assemble})
        self.refrna_assemble.on('end', self.set_step, {'end': self.step.refrna_assemble})
        self.refrna_assemble.on('end', self.set_output, 'refrna_assemble')
        self.refrna_assemble.on('end', self.run_gene_fa)
        self.refrna_assemble.on('end', self.run_annot_mapdb)
        self.refrna_assemble.on('end', self.run_annot_orfpfam)
        if "quantification" in self.analysis_content:
            self.refrna_assemble.on('end', self.run_quant)
        if self.option('is_snp') == 'True' and 'other' in self.analysis_content:
            if self.option('snp_method').lower() == 'samtools':
                self.refrna_assemble.on('end', self.run_sam_rna)
            elif self.option('snp_method').lower() == 'gatk':
                self.refrna_assemble.on('end', self.run_snp_rna)
            elif self.option('snp_method').lower() == 'sentieon':
                self.refrna_assemble.on('end', self.run_call_snp_indel)
        if self.option("is_somatic") == 'True' and 'other' in self.analysis_content:
            self.refrna_assemble.on('end', self.run_somatic_sentieon)
        if self.option('is_as') == 'True':
            if 'other' in self.analysis_content and self.option("as_method") == "rmats":
                self.refrna_assemble.on('end', self.run_rmats)
            elif 'other' in self.analysis_content and self.option("as_method").lower() == "as_profile":
                self.refrna_assemble.on('end', self.run_as_profile)
        if "other" in  self.analysis_content and self.option('fq_type') == 'PE':
            self.refrna_assemble.on('end', self.run_gene_fusion)
        self.refrna_assemble.run()

    def run_transcript_abstract(self):
        self.transcript_abstract.set_options({
            'ref_genome': self.ref_genome,
            'ref_genome_gtf': self.ref_gtf
        })
        self.transcript_abstract.on('start', self.set_step, {'start': self.step.transcript_abstract})
        self.transcript_abstract.on('end', self.set_step, {'end': self.step.transcript_abstract})
        self.transcript_abstract.on('end', self.set_output, 'transcript_abstract')
        self.transcript_abstract.on('end', self.run_gene_fa)
        self.transcript_abstract.on('end', self.run_quant)
        if self.option('sample_num') == 'multiple':
            if self.option('is_as') == 'True' and 'other' in self.analysis_content and self.option("as_method") == "rmats":
                self.transcript_abstract.on('end', self.run_rmats)
        if self.option('is_snp') == 'True' and 'other' in self.analysis_content:
            if self.option('snp_method').lower() == 'samtools':
                self.transcript_abstract.on('end', self.run_sam_rna)
            elif self.option('snp_method').lower() == 'gatk':
                self.transcript_abstract.on('end', self.run_snp_rna)
            elif self.option('snp_method').lower() == 'sentieon':
                self.transcript_abstract.on('end', self.run_call_snp_indel)
        if self.option("is_somatic") == 'True' and 'other' in self.analysis_content:
            self.transcript_abstract.on('end', self.run_somatic_sentieon)
        if "other" in  self.analysis_content:
            self.transcript_abstract.on('end', self.run_gene_fusion)
        self.transcript_abstract.run()

    def run_gene_fusion(self):
        self.step.add_steps('gene_fusion')
        opts = {
            'ref_genome_custom': self.ref_genome,
            'bam_list': os.path.join(self.rnaseq_mapping.output_dir, "bamlist"),
            'ref_gtf': self.ref_gtf,
            'circos' :"yes"
        }
        self.gene_fusion.set_options(opts)
        self.gene_fusion.on('start', self.set_step, {'start': self.step.gene_fusion})
        self.gene_fusion.on('end', self.set_step, {'end': self.step.gene_fusion})
        self.gene_fusion.on('end', self.set_output, 'gene_fusion')
        self.gene_fusion.run()


    def run_gene_fa(self):
        self.step.add_steps('gene_fa')
        if self.option('is_assemble'):
            ref_new_gtf = self.refrna_assemble.option('ref_and_new_gtf').path
        else:
            ref_new_gtf = self.ref_gtf
        self.gene_fa.set_options({
            'ref_new_gtf': ref_new_gtf,
            'ref_genome_custom': self.ref_genome
        })
        self.gene_fa.on('start', self.set_step, {'start': self.step.gene_fa})
        self.gene_fa.on('end', self.set_step, {'end': self.step.gene_fa})
        self.gene_fa.on('end', self.set_output, 'gene_fa')
        if not self.option('is_assemble'):
            self.gene_fa.on('end', self.run_detail)
        self.gene_fa.run()

    def run_detail(self):
        self.step.add_steps('detail')
        opts = {
            'gene_fa': self.gene_fa.option('gene_fa').path,
            'biomart_file': self.des,
            'biomart_type': self.des_type,
            'known_cds': self.known_cds,
            'known_pep': self.known_pep
        }
        if self.option('is_assemble'):
            opts.update({
                'txpt_fa': self.refrna_assemble.option('all_transcripts_fa').path,
                't2g': self.refrna_assemble.option('trans2gene').path,
                'novel_cds': os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.cds'),
                'novel_pep': os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.pep')
            })
        else:
            opts.update({
                'txpt_fa': self.transcript_abstract.option('trans_fa').path,
                't2g': self.transcript_abstract.option('trans2gene').path
            })
        self.detail.set_options(opts)
        self.detail.on('start', self.set_step, {'start': self.step.detail})
        self.detail.on('end', self.set_step, {'end': self.step.detail})
        self.detail.on('end', self.set_output, 'detail')
        self.detail.run()

    def run_annot_mapdb(self):
        self.step.add_steps('annot_mapdb')
        self.annot_mapdb.set_options({
            'query': self.refrna_assemble.option('new_transcripts_fa').path,
            'nr_db': self.option('nr_database'),
            'kegg_version': self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version']
        })
        self.annot_mapdb.on('start', self.set_step, {'start': self.step.annot_mapdb})
        self.annot_mapdb.on('end', self.set_step, {'end': self.step.annot_mapdb})
        self.annot_mapdb.on('end', self.set_output, 'annot_mapdb')
        # self.annot_mapdb.on('end',self.run_annot_class_beta_new)
        self.annot_mapdb.run()

    def run_annot_orfpfam(self):
        self.step.add_steps('annot_orfpfam')
        self.annot_orfpfam.set_options({
            'fasta': self.refrna_assemble.option('new_transcripts_fa').path,
            'gtf': self.refrna_assemble.option('new_transcripts_gtf').path,
            "pfam_version": self.annot_config_dict['pfam']['version']
        })
        self.annot_orfpfam.on('start', self.set_step, {'start': self.step.annot_orfpfam})
        self.annot_orfpfam.on('end', self.set_step, {'end': self.step.annot_orfpfam})
        self.annot_orfpfam.on('end', self.set_output, 'annot_orfpfam')
        self.annot_orfpfam.run()

    # def run_annot_filter_new(self):
    #     self.step.add_steps('annot_filter_new')
    #     self.annot_filter_new.set_options({
    #         'blast_nr_xml': os.path.join(self.annot_mapdb.output_dir, 'nr/blast.xml'),
    #         'blast_eggnog_xml': os.path.join(self.annot_mapdb.output_dir, 'eggnog/blast.xml'),
    #         'blast_kegg_xml': os.path.join(self.annot_mapdb.output_dir, 'kegg/blast.xml'),
    #         'blast_swissprot_xml': os.path.join(self.annot_mapdb.output_dir, 'swissprot/blast.xml'),
    #         'pfam_domain': os.path.join(self.annot_orfpfam.output_dir, 'pfam_domain'),
    #         'blast2go_annot': os.path.join(self.annot_mapdb.output_dir, 'GO/blast2go_merge.xls'),
    #         'nr_evalue': self.option('nr_evalue'),
    #         'nr_identity': self.option('nr_identity'),
    #         'nr_similarity': self.option('nr_similarity'),
    #         'swissprot_evalue': self.option('swissprot_evalue'),
    #         'swissprot_identity': self.option('swissprot_identity'),
    #         'swissprot_similarity': self.option('swissprot_similarity'),
    #         'eggnog_evalue': self.option('cog_evalue'),
    #         'eggnog_identity': self.option('cog_identity'),
    #         'eggnog_similarity': self.option('cog_similarity'),
    #         'kegg_evalue': self.option('kegg_evalue'),
    #         'kegg_identity': self.option('kegg_identity'),
    #         'kegg_similarity': self.option('kegg_similarity'),
    #         'pfam_evalue': self.option('pfam_evalue')
    #     })
    #     self.annot_filter_new.on('start', self.set_step, {'start': self.step.annot_filter_new})
    #     self.annot_filter_new.on('end', self.set_step, {'end': self.step.annot_filter_new})
    #     self.annot_filter_new.on('end', self.set_output, 'annot_filter_new')
    #     self.annot_filter_new.on('end', self.run_annot_class_beta_new)
    #     self.annot_filter_new.run()

    def run_annot_class_beta_new(self):
        self.step.add_steps('annot_class_beta_new')
        self.annot_class_beta_new.set_options({
            'type': 'new',
            'gene2trans': os.path.join(self.annot_orfpfam.output_dir, 'all_tran2gen.txt'),
            'des': self.des,
            'des_type': self.des_type,
            'enterz': self.entrez,
            'blast_nr_xml': os.path.join(self.annot_mapdb.output_dir, 'nr/blast.xml'),
            'blast_kegg_xml': os.path.join(self.annot_mapdb.output_dir, 'kegg/blast.xml'),
            'taxonomy': self.option('kegg_database'),
            'link_bgcolor': 'green',
            'png_bgcolor': '00CD00',
            'blast_eggnog_xml': os.path.join(self.annot_mapdb.output_dir, 'eggnog/blast.xml'),
            'blast_uniprot_xml': os.path.join(self.annot_mapdb.output_dir, 'uniprot/blast.xml'),
            'pfam_domain': os.path.join(self.annot_orfpfam.output_dir, 'pfam_domain'),
            'blast2go_annot': os.path.join(self.annot_mapdb.output_dir, 'GO/blast2go_merge.xls'),
            'kegg_version': self.annot_config_dict['kegg']['version'],
            'kegg_subtax1': self.option('kegg_subtax1'),
            'kegg_subtax2': self.option('kegg_subtax2'),
            'kegg_species': self.option('kegg_specific'),
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],
        })
        self.annot_class_beta_new.on('start', self.set_step, {'start': self.step.annot_class_beta_new})
        self.annot_class_beta_new.on('end', self.set_step, {'end': self.step.annot_class_beta_new})
        self.annot_class_beta_new.on('end', self.set_output, 'annot_class_beta_new')
        self.annot_class_beta_new.run()

    def run_annot_merge(self):
        self.step.add_steps('annot_merge')
        opts = {
            'ref_class_dir': os.path.join(self.ref_annot_dir, 'annot_class_medical'),
            'is_assemble': self.option('is_assemble'),
            'kegg_version': self.annot_config_dict['kegg']['version'],
            'kegg_species': self.option('kegg_specific')
        }
        if self.option('is_assemble'):
            opts.update({
                'new_class_dir': self.annot_class_beta_new.output_dir,
                'new_mapdb_dir': self.annot_mapdb.output_dir
            })
        self.annot_merge.set_options(opts)
        self.annot_merge.on('start', self.set_step, {'start': self.step.annot_merge})
        self.annot_merge.on('end', self.set_step, {'end': self.step.annot_merge})
        self.annot_merge.on('end', self.set_output, 'annot_merge')
        self.annot_merge.run()

    def run_quant(self):
        if self.option('datatype') == 'rawdata':
            if self.option('qc_soft') == 'fastp':
                fastq = self.fastp_rna.option('fq_list').path
            elif self.option('qc_soft') == 'seqprep':
                fastq = self.hiseq_qc.option('fq_list').path
        elif self.option('datatype') == 'cleandata':
            fastq = self.get_fq_list()
        if self.option('strand_specific'):
            if self.option('strand_dir').startswith('R'):
                if self.option('fq_type') == 'PE':
                    libtype = 'rf'
                elif self.option('fq_type') == 'SE':
                    libtype = 'r'
            elif self.option('strand_dir').startswith('F'):
                if self.option('fq_type') == 'PE':
                    libtype = 'fr'
                elif self.option('fq_type') == 'SE':
                    libtype = 'f'
        else:
            libtype = None
        if self.option('is_assemble'):
            transcriptome = self.refrna_assemble.option('all_transcripts_fa').path
            t2g = self.refrna_assemble.option('trans2gene').path
        else:
            transcriptome = self.transcript_abstract.option('trans_fa')
            t2g = self.transcript_abstract.option('trans2gene').path
        self.quant.set_options({
            'fastq': fastq,
            'method': self.option('express_method'),
            'libtype': libtype,
            'transcriptome': transcriptome,
            't2g': t2g,
            'ref_gtf': self.file_check.option('gtf').path
        })
        self.quant.on('start', self.set_step, {'start': self.step.quant})
        self.quant.on('end', self.set_step, {'end': self.step.quant})
        self.quant.on('end', self.set_output, 'quant')
        if self.option('sample_num') == 'multiple':
            if len(self.option('group_table').prop['group_dict']) > 1:
                self.quant.on('end', self.run_exp_venn)
            if self.option('group_table').prop['sample_number'] > 2:
                self.quant.on('end', self.run_exp_pca)
                self.quant.on('end', self.run_exp_corr)
            if "other" in self.analysis_content:
                self.quant.on('end', self.run_diffexp)
        self.quant.run()

    def run_exp_venn(self):
        if self.option('exp_way').lower() == 'fpkm':
            express_matrix = self.quant.option('ref_gene_fpkm').path
        else:
            express_matrix = self.quant.option('ref_gene_tpm').path
        group_dict = self.option('group_table').prop['group_dict']
        if len(group_dict) <= 6:
            group_table = self.option('group_table').path
        else:
            group_table = '{}.6'.format(self.option('group_table').path)
            lines = ['#sample\tgroup\n']
            for i, (g, ss) in enumerate(group_dict.items()):
                if i < 6:
                    lines.extend('{}\t{}\n'.format(s, g) for s in ss)
            else:
                open(group_table, 'w').writelines(lines)
        self.exp_venn.set_options({
            'express_matrix': express_matrix,
            'group_table': group_table,
            'threshold': 1
        })
        self.exp_venn.on('start', self.set_step, {'start': self.step.exp_venn})
        self.exp_venn.on('end', self.set_step, {'end': self.step.exp_venn})
        self.exp_venn.on('end', self.set_output, 'exp_venn')
        self.exp_venn.run()

    def run_exp_pca(self):
        if self.option('exp_way').lower() == 'fpkm':
            exp = self.quant.option('ref_gene_fpkm').path
        else:
            exp = self.quant.option('ref_gene_tpm').path
        self.exp_pca.set_options({'exp': exp})
        self.exp_pca.on('start', self.set_step, {'start': self.step.exp_pca})
        self.exp_pca.on('end', self.set_step, {'end': self.step.exp_pca})
        self.exp_pca.on('end', self.set_output, 'exp_pca')
        if hasattr(self, 'ellipse'):
            self.exp_pca.on('end', self.run_ellipse)
        self.exp_pca.run()

    def run_ellipse(self):
        self.ellipse.set_options({
            'analysis': 'pca',
            'group_table': self.option('group_table').prop['path'],
            'pc_table': os.path.join(self.exp_pca.output_dir, 'PCA.xls'),
        })
        self.ellipse.on('start', self.set_step, {'start': self.step.ellipse})
        self.ellipse.on('end', self.set_step, {'end': self.step.ellipse})
        self.ellipse.on('end', self.set_output, 'exp_pca')
        self.ellipse.run()

    def run_exp_corr(self):
        if self.option('exp_way').lower() == 'fpkm':
            exp = self.quant.option('ref_gene_fpkm').path
        else:
            exp = self.quant.option('ref_gene_tpm').path
        self.exp_corr.set_options({'exp': exp})
        self.exp_corr.on('start', self.set_step, {'start': self.step.exp_corr})
        self.exp_corr.on('end', self.set_step, {'end': self.step.exp_corr})
        self.exp_corr.on('end', self.set_output, 'exp_corr')
        self.exp_corr.run()

    def run_diffexp(self):
        if self.option('exp_way').lower() == 'fpkm':
            exp = self.quant.option('ref_gene_fpkm').path
        else:
            exp = self.quant.option('ref_gene_tpm').path
        opts = {
            'is_workflow':"yes",
            'count': self.quant.option('ref_gene_count').path,
            'exp': exp,
            'group': self.option('group_table'),
            'cmp': self.option('control_file'),
            'fc': float(self.option('fc')),
            'method': self.option('diff_method'),
            'exp_type': self.option('exp_way'),
        }

        if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
            opts.update({
                'pvalue_padjust': self.option('pvalue_padjust'),
                'pvalue': float(self.option('diff_fdr_ci')),
                'padjust_way': self.option('padjust_way'),
            })
            if self.option('pair_table').is_set:
                opts.update({'is_batch': True, 'has_batch': True, 'batch_matrix': self.option('pair_table').path})
            self.diffexp.set_options(opts)
        if self.option('diff_method').lower() == 'noiseq':
            # self.option('diff_fdr_ci',0.8)
            opts.update({
                "prob": float(self.option('diff_fdr_ci')),
            })
            if self.option('pair_table').is_set:
                opts.update({'is_batch': True, 'has_batch': True, 'batch_matrix': self.option('pair_table').path})
            self.diffexp.set_options(opts)

        self.diffexp.on('start', self.set_step, {'start': self.step.diffexp})
        self.diffexp.on('end', self.set_step, {'end': self.step.diffexp})
        self.diffexp.on('end', self.set_output, 'diffexp')
        self.diffexp.run()

    def run_diff_geneset_analysis(self):
        #modify by 20210927  获取最终差异参数
        diff_option_file = os.path.join(self.diffexp.work_dir,"option_file")
        with open(diff_option_file,"r") as r:
            diff_options = json.load(r)
        self.option("diff_method",diff_options["method"])
        if self.option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            self.option('pvalue_padjust',diff_options["pvalue_padjust"])
            self.option("fc", diff_options["fc"])

        else:
            # 如果是NOIseq,调整次数参数以适应后面导表
            self.option("diff_fdr_ci", diff_options["prob"])
        if self.option('exp_way').lower() == 'fpkm':
            exp = self.quant.option('ref_gene_fpkm').path
        else:
            exp = self.quant.option('ref_gene_tpm').path
        opts = {
            'diff_path': self.diffexp.output_dir,
            'annot_result': self.annot_merge.output_dir,
            'gene_count_file': self.quant.option('ref_gene_count').path,
            'diff_method': self.option('diff_method'),
            'gene_exp_file': exp,
            'group': self.option('group_table').path,
            'level': "G",
            "kegg_version":self.annot_config_dict['kegg']['version'],
            "reactome_version":self.annot_config_dict['reactome']['version'],
            'species': self.organism_name,
        }
        self.diff_geneset_analysis.set_options(opts)
        self.diff_geneset_analysis.on('start', self.set_step, {'start': self.step.diff_geneset_analysis})
        self.diff_geneset_analysis.on('end', self.set_step, {'end': self.step.diff_geneset_analysis})
        self.diff_geneset_analysis.on('end', self.set_output, 'diff_geneset_analysis')
        self.diff_geneset_analysis.run()

    def run_sam_rna(self):
        if self.option('is_assemble'):
            ref_gtf = self.refrna_assemble.option('ref_and_new_gtf').path
        else:
            ref_gtf = self.ref_gtf
        self.sam_rna.set_options({
            'ref_genome': 'customer_mode',
            'ref_genome_custom': self.ref_genome,
            'ref_gtf': ref_gtf,
            # 'bamlist': self.rnaseq_mapping.option('bamlist'),
            'in_bam':self.rnaseq_mapping.option('bam_output'),
            'des': self.des,
            'des_type': self.des_type
        })
        self.sam_rna.on('start', self.set_step, {'start': self.step.sam_rna})
        self.sam_rna.on('end', self.set_step, {'end': self.step.sam_rna})
        self.sam_rna.on('end', self.set_output, 'sam_rna')
        self.sam_rna.run()

    def run_snp_rna(self):
        if self.option('is_assemble'):
            ref_gtf = self.refrna_assemble.option('ref_and_new_gtf').path
        else:
            ref_gtf = self.ref_gtf
        self.snp_rna.set_options({
            'ref_genome_custom': self.ref_genome,
            'ref_gtf': ref_gtf,
            'in_bam': self.rnaseq_mapping.option('bam_output'),
            'des': self.des,
            'des_type': self.des_type
        })
        self.snp_rna.on('start', self.set_step, {'start': self.step.snp_rna})
        self.snp_rna.on('end', self.set_step, {'end': self.step.snp_rna})
        self.snp_rna.on('end', self.set_output, 'snp_rna')
        self.snp_rna.run()




    def run_call_snp_indel(self):
        if self.option('is_assemble'):
            ref_gtf = self.refrna_assemble.option('ref_and_new_gtf').path
        else:
            ref_gtf = self.ref_gtf
        opts = {
            'ref_fasta': self.ref_genome,
            'bam_list': os.path.join(self.rnaseq_mapping.output_dir, "bamlist"),
            'des': self.des,
            'des_type': self.des_type,
            'ref_gtf': ref_gtf,
            "align_method": self.option("align_method").lower()
        }
        self.call_snp_indel.set_options(opts)
        self.call_snp_indel.on('start', self.set_step, {'start': self.step.call_snp_indel})
        self.call_snp_indel.on('end', self.set_step, {'end': self.step.call_snp_indel})
        self.call_snp_indel.on('end', self.set_output, 'call_snp_indel')
        self.call_snp_indel.run()

    def run_somatic_sentieon(self):
        if self.option('is_assemble'):
            ref_gtf = self.refrna_assemble.option('ref_and_new_gtf').path
        else:
            ref_gtf = self.ref_gtf
        opts = {
            'ref_fasta': self.ref_genome,
            'bam_list': os.path.join(self.rnaseq_mapping.output_dir, "bamlist"),
            'des': self.des,
            'des_type': self.des_type,
            'ref_gtf': ref_gtf,
            "align_method": self.option("align_method").lower(),
            'tn_pair' :self.option("tumor_pair_table").prop["path"]
        }
        self.somatic_analysis.set_options(opts)
        self.somatic_analysis.on('start', self.set_step, {'start': self.step.somatic_sentieon})
        self.somatic_analysis.on('end', self.set_step, {'end': self.step.somatic_sentieon})
        self.somatic_analysis.on('end', self.set_output, 'somatic')
        self.somatic_analysis.run()

    def run_rmats(self):
        self.step.add_steps('rmats')
        if self.option('fq_type') == 'PE':
            seq_type = 'paired'
        elif self.option('fq_type') == 'SE':
            seq_type = 'single'
        if self.option('strand_specific'):
            if self.option('strand_dir').startswith('R'):
                lib_type = 'fr-firststrand'
            elif self.option('strand_dir').startswith('F'):
                lib_type = 'fr-secondstrand'
        else:
            lib_type = 'fr-unstranded'
        if self.option('is_assemble'):
            input_gtf = self.refrna_assemble.option('ref_and_new_gtf').prop["path"]
        else:
            input_gtf = self.ref_gtf
        self.as_analysis.set_options({
            'control_table': self.option('control_file').path,
            'group_table': self.option('group_table').path,
            'bam_input': self.rnaseq_mapping.option('bam_output'),
            'input_gtf': input_gtf,
            'seq_type': seq_type,
            'lib_type': lib_type
        })
        self.as_analysis.on('start', self.set_step, {'start': self.step.rmats})
        self.as_analysis.on('end', self.set_step, {'end': self.step.rmats})
        self.as_analysis.on('end', self.set_output, 'rmats')
        self.as_analysis.run()

    def run_as_profile(self):
        ## Stringtie for ASprofile
        self.intermediate_result = os.path.join(self.work_dir, 'intermediate_results')
        if self.option('is_assemble'):
            gtf_dir = os.path.join(self.output_dir, 'assembly/Stringtie/')
            files = glob.glob(gtf_dir + '/' + '*.gtf')
            if os.path.exists(os.path.join(self.intermediate_result, 'Stringtie')):
                shutil.rmtree(os.path.join(self.intermediate_result, 'Stringtie'))
            os.makedirs(os.path.join(self.intermediate_result, 'Stringtie'))
            for file in files:
                os.link(file, os.path.join(self.intermediate_result, 'Stringtie', os.path.basename(file)))
            gtf_path = os.path.join(self.intermediate_result, 'Stringtie/')
            gtfs = os.listdir(gtf_path)
            gtf_sample = os.path.join(gtf_path, 'list')
            with open(gtf_sample, 'w') as file1:
                for gtf in gtfs:
                    sample = re.findall(r'(.*?)_out', gtf)[0]
                    file1.write(gtf + '\t' + sample + '\n')
        self.as_analysis.set_options({
            'gtf_dir':os.path.join(self.intermediate_result, 'Stringtie'),
            'ref_fa': self.ref_genome,
            'ref_gtf': self.ref_gtf,
            'group_table' :self.option('group_table').path
        })
        self.as_analysis.on('start', self.set_step, {'start': self.step.rmats})
        self.as_analysis.on('end', self.set_step, {'end': self.step.rmats})
        self.as_analysis.on('end', self.set_output, 'as_profile')
        self.as_analysis.run()



    def stop(self, reason):
        assert reason in ['rrna', 'mapping']
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.export_task_info()
        self.export_genome_info()
        self.export_sample_info()
        self.export_qc_result()
        if reason == 'rrna':
            msg = 'Workflow will stop because the rRNA ratio is not up to par'
            self.logger.warn(msg)
        elif reason == 'mapping':
            self.export_ref_rna_qc_alignment()
            msg = 'Workflow will stop because the alignment rate is not up to par'
            self.logger.warn(msg)
        receiver = ['caiping.shi@majorbio.com', 'rna_bioinfor@majorbio.com']
        a = SendEmail("897236887@qq.com", "smtp.qq.com", "fhwuvcclstjqbfga", "897236887@qq.com", ','.join(receiver),
                      "WARNING - project_sn ({}), task_id ({})".format(self._sheet.project_sn, self._sheet.id), 465)
        a.send_msg(msg)
        a.send_email()
        # db = Config().get_mongo_client(mtype='project', dydb_forbid=True)[Config().get_mongo_dbname(mtype='project', dydb_forbid=True)]
        # email = db['sg_task_email']
        # email.update({'task_id': self.task_id}, {'$set': {'status': '5'}}, upsert=True)
        super(MedicalTranscriptomeWorkflow, self).end()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] in ['fastp_rna', 'hiseq_qc']:
            self.move2outputdir(obj.output_dir, 'QC_stat')
        if event['data'] == 'hiseq_reads_stat_raw':
            self.move2outputdir(obj.output_dir, 'QC_stat/before_qc')
        if event['data'] == 'hiseq_reads_stat_use':
            self.move2outputdir(obj.output_dir, 'QC_stat/after_qc')
        if event['data'] == 'rnaseq_mapping':
            self.move2outputdir(obj.output_dir, 'mapping')
        if event['data'] == 'map_assessment':
            self.move2outputdir(obj.output_dir, 'map_qc')
        if event['data'] == 'refrna_assemble':
            self.move2outputdir(obj.output_dir, 'assembly')
        if event['data'] in ['sam_rna', 'snp_rna', 'call_snp_indel']:
            self.move2outputdir(obj.output_dir, 'snp')
        if event['data'] == 'quant':
            self.move2outputdir(obj.output_dir, 'express')
        if event['data'] == 'exp_pca':
            self.move2outputdir(obj.output_dir, 'exp_pca')
        if event['data'] == 'exp_corr':
            self.move2outputdir(obj.output_dir, 'exp_corr')
        if event['data'] == 'diffexp':
            self.move2outputdir(obj.output_dir, 'diffexpress')
        if event['data'] == 'rmats':
            self.move2outputdir(obj.output_dir, 'altersplicing')
        if event['data'] == 'as_profile':
            self.move2outputdir(obj.output_dir, 'altersplicing')
        if event['data'] == 'annot_mapdb':
            self.move2outputdir(obj.output_dir, 'annot_mapdb')
        if event['data'] == 'annot_orfpfam':
            self.move2outputdir(obj.output_dir, 'annot_orfpfam')
        if event['data'] == 'annot_merge':
            self.move2outputdir(obj.output_dir, 'annot_merge')
        if event['data'] == 'transcript_abstract':
            self.move2outputdir(obj.output_dir, 'abstract_transcripts')
        if event['data'] == 'annot_class_beta_new':
            self.move2outputdir(obj.output_dir, 'annot_class_new')
        if event['data'] == 'run_annot_class_beta_ref':
            self.move2outputdir(obj.output_dir, 'annot_class_ref')
        if event['data'] == 'annot_filter_new':
            self.move2outputdir(obj.output_dir, 'annot_filter_new')
        if event['data'] == 'annot_filter_ref':
            self.move2outputdir(obj.output_dir, 'annot_filter_ref')
        if event['data'] == 'gene_fusion':
            self.move2outputdir(obj.output_dir, 'gene_fusion')
        if event['data'] == 'diff_geneset_analysis':
            self.move2outputdir(obj.output_dir, 'diff_geneset_analysis')
        if event['data'] == 'somatic':
            self.move2outputdir(obj.output_dir, 'somatic')


    @workfuncdeco
    def move2outputdir(self, olddir, newname):
        if not os.path.isdir(olddir):
            self.set_error('can not find source directory -> (%s)', variables=(olddir), code="13700337")
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])

    def move_file(self, src, dst):
        if os.path.isfile(src):
            os.link(src, dst)
        else:
            os.mkdir(dst)
            for file in os.listdir(src):
                old_path = os.path.join(src, file)
                new_path = os.path.join(dst, file)
                self.move_file(old_path, new_path)
        # self.logger.info('succeed in moving {} to {}'.format(src, dst))

    def set_db(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.stop_timeout_check()
        if os.path.exists(os.path.join(self.work_dir, "ttttt")):
            self.set_upload()
        else:
            if os.path.exists(os.path.join(self.work_dir,"temporary")):
                shutil.rmtree(os.path.join(self.work_dir,"temporary"))
            os.makedirs(os.path.join(self.work_dir,"temporary"))
            self.export_temporary = os.path.join(self.work_dir,"temporary")
            if "mapping" in self.analysis_content:
                self.export_task_info()
                self.export_genome_info()
                self.export_sample_info()
                self.export_qc_result()

                '''
                # if self.option('datatype') == 'rawdata':
                #     self.export_ref_rna_qc_before()
                # self.export_ref_rna_qc_after()
                '''

                self.export_ref_rna_qc_alignment()
                self.export_ref_rna_qc_assessment()
            if "annotation" in self.analysis_content:
                if self.option('is_assemble'):
                    self.export_ref_assembly()
                self.export_annotation()
            if "quantification" in self.analysis_content:
                self.export_all_exp_matrix()
                if self.option('sample_num') == 'multiple':
                    self.export_all_exp_distribution()
                    if len(self.option('group_table').prop['group_dict']) > 1:
                        self.export_add_exp_venn()
                    if self.option('group_table').prop['sample_number'] > 2:
                        self.export_all_exp_pca()
                        self.export_all_exp_corr()
                    self.export_gene_detail()
            if "other" in self.analysis_content:
                if self.option('fq_type') == 'PE':
                    self.export_gene_fusion()
                if self.option('sample_num') == 'multiple':
                    self.export_all_exp_diff()
                    self.export_diff_geneset_analysis()
                    if self.option('is_as') == 'True':
                        self.export_rmats()
                        # self.export_rmats_count()
                if self.option('is_snp') == 'True':
                    self.export_snp()
                if self.option('is_somatic') == 'True':
                    self.export_somatic()

            if self.option("report_img"):
                self.export_report_img()
            self.set_upload()

    # @workfuncdeco
    def set_upload(self):
        self.logger.info('kaishilekaishile1')
        if "quantification" in self.analysis_content:
            self.merge_annotation_exp_matrix()  # 表达量表增加注释信息
        if "other" in self.analysis_content:
            if self.option('sample_num') == 'multiple':
                self.merge_annotation_diffexp_matrix()  # 差异表达量表增加注释信息
        if self.option("upload_offline"):
            # 将质控结果文件和比对结果文件，传输至线下服务器
            self.offline_results = os.path.join(self.work_dir, self.task_id)
            if os.path.isdir(self.offline_results):
                shutil.rmtree(self.offline_results)
            os.mkdir(self.offline_results)
            # os.mkdir(os.path.join(self.offline_results, "QC"))
            # if self.option('qc_soft') == 'fastp':
            #     fastq = self.fastp_rna.option('fq_list').path
            # elif self.option('qc_soft') == 'seqprep':
            #     fastq = self.hiseq_qc.option('fq_list').path
            # with open(fastq, "r") as f:
            #     for line in f:
            #         items = line.strip().split("\t")[1:]
            #         for item in items:
            #             os.link(item, os.path.join(self.offline_results, "QC", os.path.basename(item)))
            os.mkdir(os.path.join(self.offline_results, "Mapping"))
            bam_files = glob.glob(self.rnaseq_mapping.output_dir + "/bam/*.bam")
            for bam in bam_files:
                os.link(bam, os.path.join(self.offline_results, "Mapping", os.path.basename(bam)))

        # 页面交互运行所需中间结果文件，不在结果目录中呈现
        self.intermediate_result = os.path.join(self.work_dir, 'intermediate_results')
        if os.path.isdir(self.intermediate_result):
            shutil.rmtree(self.intermediate_result)
        os.mkdir(self.intermediate_result)

        if "mapping" in self.analysis_content:
            ## AlignBam
            os.makedirs(os.path.join(self.intermediate_result, 'Align/AlignBam'))
            for b in os.listdir(os.path.join(self.output_dir, 'mapping/bam')):
                os.link(os.path.join(self.output_dir, 'mapping/bam', b),
                        os.path.join(self.intermediate_result, 'Align/AlignBam', b))
        if "annotation" in self.analysis_content:
            ## Annotation
            os.makedirs(os.path.join(self.intermediate_result, 'Annotation'))
            CopyFile().linkdir(os.path.join(self.output_dir, 'annot_merge'),
                               os.path.join(self.intermediate_result, 'Annotation'))
            rm_files = glob.glob(os.path.join(self.intermediate_result, 'Annotation/*/*/*/*.KOs.txt')) + \
                       glob.glob(os.path.join(self.intermediate_result, 'Annotation/*/*/*/*.pdf')) + \
                       glob.glob(os.path.join(self.intermediate_result, 'Annotation/*/cog/summary.*.tsv')) + \
                       glob.glob(os.path.join(self.intermediate_result, 'Annotation/*/kegg/kegg_layer_*.xls'))
            for rm_file in rm_files:
                os.remove(rm_file)
            ## Transcripts
            os.makedirs(os.path.join(self.intermediate_result, 'Transcripts'))
            if not self.option('is_assemble'):
                os.link(os.path.join(self.transcript_abstract.work_dir, 'trans2gene'),
                        os.path.join(self.intermediate_result, 'Transcripts/trans2gene.txt'))
                os.link(os.path.join(self.transcript_abstract.work_dir, 'exons.fa'),
                        os.path.join(self.intermediate_result, 'Transcripts/all_transcripts.fa'))
            ## Stringtie for ASprofile
            if self.option('is_assemble'):
                gtf_dir = os.path.join(self.output_dir, 'assembly/Stringtie/')
                files = glob.glob(gtf_dir + '/' + '*.gtf')
                os.makedirs(os.path.join(self.intermediate_result, 'Stringtie'))
                for file in files:
                    os.link(file, os.path.join(self.intermediate_result, 'Stringtie', os.path.basename(file)))
                gtf_path = os.path.join(self.intermediate_result, 'Stringtie/')
                gtfs = os.listdir(gtf_path)
                gtf_sample = os.path.join(gtf_path, 'list')
                with open(gtf_sample, 'w') as file1:
                    for gtf in gtfs:
                        sample = re.findall(r'(.*?)_out', gtf)[0]
                        file1.write(gtf + '\t' + sample + '\n')

        if "quantification" in self.analysis_content:
            ## SequenceDatabase
            os.mkdir(os.path.join(self.intermediate_result, 'SequenceDatabase'))

            seq_db = self.detail.option('database').path
            os.link(seq_db, os.path.join(self.intermediate_result, 'SequenceDatabase', os.path.basename(seq_db)))
            seq_detail = os.path.join(self.detail.output_dir,"detail")
            shutil.copytree(seq_detail,os.path.join(self.intermediate_result, 'SequenceDetail'))


            ## Express
            os.makedirs(os.path.join(self.intermediate_result, 'Express/ExpAnnalysis'))
            if self.option('express_method') == 'RSEM':
                os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.fpkm.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.fpkm.matrix.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.count.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.annot.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.fpkm.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.fpkm.matrix.annot.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.count.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.annot.xls'))
                if self.option('level').lower() == 'transcript':
                    os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix'),
                            os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.tpm.matrix.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.fpkm.matrix'),
                            os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.fpkm.matrix.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.count.matrix'),
                            os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/transcript.count.matrix.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix.annot.xls'),
                            os.path.join(self.intermediate_result,
                                         'Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.fpkm.matrix.annot.xls'),
                            os.path.join(self.intermediate_result,
                                         'Express/ExpAnnalysis/transcript.fpkm.matrix.annot.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.count.matrix.annot.xls'),
                            os.path.join(self.intermediate_result,
                                         'Express/ExpAnnalysis/transcript.count.matrix.annot.xls'))
            else:
                os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.count.matrix'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.tpm.matrix.annot.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.count.matrix.annot.xls'),
                        os.path.join(self.intermediate_result, 'Express/ExpAnnalysis/gene.count.matrix.annot.xls'))
                if self.option('level').lower() == 'transcript':
                    os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix'),
                            os.path.join(self.intermediate_result,
                                         'Express/ExpAnnalysis/transcript.tpm.matrix.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.count.matrix'),
                            os.path.join(self.intermediate_result,
                                         'Express/ExpAnnalysis/transcript.count.matrix.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix.annot.xls'),
                            os.path.join(self.intermediate_result,
                                         'Express/ExpAnnalysis/transcript.tpm.matrix.annot.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.count.matrix.annot.xls'),
                            os.path.join(self.intermediate_result,
                                         'Express/ExpAnnalysis/transcript.count.matrix.annot.xls'))

        if "other" in self.analysis_content:
            if self.option('sample_num') == 'multiple':
                ## Diffexpress
                os.makedirs(os.path.join(self.intermediate_result, 'Diffexpress'))
                diff_output = self.diffexp.output_dir
                files = glob.glob(diff_output + '/' + '*_vs_*.normalize.xls') + glob.glob(diff_output + '/' + '*_vs_*.sizeFactor.xls')
                for file in files:
                    os.link(file, os.path.join(self.intermediate_result, 'Diffexpress', os.path.basename(file)))
                ## SNP
                os.makedirs(os.path.join(self.intermediate_result, 'SNP'))
                if self.option('is_snp') == 'True':
                    if self.option('snp_method').lower() == 'samtools':
                        final_vcf = os.path.join(self.sam_rna.vcf_filter.output_dir,"final.vcf")
                        # final_vcf = os.path.join(self.sam_rna.work_dir, "VcfFilterSamtools/output/final.vcf")
                    elif self.option('snp_method').lower() == 'gatk':
                        final_vcf = os.path.join(self.snp_rna.vcffilter.output_dir, "final.vcf")
                        # final_vcf = os.path.join(self.snp_rna.work_dir, "VcfFilterGatk/output/final.vcf")
                    else:
                        final_vcf = os.path.join(self.call_snp_indel.vcffilter.output_dir, "final.vcf")
                        # final_vcf = os.path.join(self.call_snp_indel.work_dir, "VcfFilterGatk/output/final.vcf")
                    os.link(final_vcf, os.path.join(self.intermediate_result, 'SNP', os.path.basename(final_vcf)))
                ## AS
                if self.option('is_as') == 'True':
                    shutil.copytree(os.path.join(self.output_dir, 'altersplicing'),
                                    os.path.join(self.intermediate_result, 'AS'))

        # 项目结果目录中呈现文件
        self.output_dir = self.output_dir
        self.target_dir = os.path.join(self.work_dir, 'upload')
        if os.path.isdir(self.target_dir):
            shutil.rmtree(self.target_dir)
        os.mkdir(self.target_dir)

        # 04Background 医学版结果目录结构有所变更，01为差异未必每个项目均有，因此先调整后续部分
        os.mkdir(os.path.join(self.target_dir, '04Background'))
        ## get run parameter
        if self.option("get_run_log"):
            self.run_parameter(os.path.join(self.target_dir, '04Background'))
        # ## genome_stat
        # genome_stat = pd.read_table(self.genome_stat, header=0)
        # genome_stat.rename(columns={'Chr': 'Chromosome', 'Size(Mb)': 'Length(MB)'}, inplace=True)
        # genome_stat.to_csv(os.path.join(self.target_dir, '01Background', self.organism_name + ".genome_info.xls"),
        #                    sep="\t", header=True, index=False)
        ## software_info
        software_info = os.path.join(self.target_dir, '04Background', "software_info.xls")
        db = Config().get_mongo_client(mtype='medical_transcriptome', dydb_forbid=True)[Config().get_mongo_dbname(mtype='medical_transcriptome', dydb_forbid=True)]
        my_collection = db['sg_software_database']
        my_results = my_collection.find({})
        with open(software_info, "w") as w:
            w.write("\t".join(["Soft/Database", "Version", "Analysis", "Source"]) + "\n")
            for collection in my_results:
                w.write("\t".join(
                    [str(collection["software_database"]), str(collection["version"]), str(collection["usage"]),
                     str(collection["source"])]) + "\n")
        ## sample_info
        sample_info = os.path.join(self.target_dir, '04Background', "sample_info.xls")
        if self.option('sample_num') == 'multiple':
            group_table = self.option('group_table').prop["path"]
        else:
            group_table = os.path.join(self.work_dir, 'group.txt')
        productive_names = dict()
        mj_names = dict()
        if self.option("productive_table").is_set:
            with open(self.option("productive_table").path, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    items = line.strip().split("\t")
                    if len(items) >= 2:
                        productive_names[items[0]] = items[1]
                    if len(items) >= 3:
                        mj_names[items[0]] = items[2]
        with open(sample_info, "w") as w, open(group_table, "r") as f:
            if self.option("productive_table").is_set:
                if mj_names:
                    w.write(
                        "\t".join(["Species Latin Name", "Sample Productive Name", "MJ_name","Sample Name", "Group Name"]) + "\n")
                else:
                    w.write("\t".join(["Species Latin Name", "Sample Productive Name", "Sample Name", "Group Name"]) + "\n")
            else:
                w.write("\t".join(["Species Latin Name", "Sample Name", "Group Name"]) + "\n")
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    if line.strip().split("\t")[0] in productive_names:
                        if line.strip().split("\t")[0] in mj_names:
                            w.write("\t".join([self.organism_name, productive_names[line.strip().split("\t")[0]],
                                               mj_names[line.strip().split("\t")[0]], line.strip().split("\t")[0],
                                               line.strip().split("\t")[1]]) + "\n")
                        else:
                            w.write("\t".join([self.organism_name, productive_names[line.strip().split("\t")[0]],
                                               line.strip().split("\t")[0], line.strip().split("\t")[1]]) + "\n")
                    else:
                        w.write("\t".join(
                            [self.organism_name, line.strip().split("\t")[0], line.strip().split("\t")[1]]) + "\n")

        # 05Basic_Analysis 基础分析部分数据
        os.mkdir(os.path.join(self.target_dir,'05Basic_Analysis'))
        os.mkdir(os.path.join(self.target_dir,'05Basic_Analysis','01QC'))
        if self.option('datatype') == 'rawdata':
            raw_stat = os.path.join(self.output_dir, 'QC_stat/before_qc/fastq_stat.xls')
            raw_df = pd.read_table(raw_stat)
            raw_df = raw_df[["#Sample_ID", "Total_Reads", "Total_Bases"]]
            raw_df = raw_df.rename(
                columns={"#Sample_ID": "Sample", "Total_Reads": "Raw reads", "Total_Bases": "Raw bases"})
            raw_df = raw_df.set_index("Sample")
        clean_stat = os.path.join(self.output_dir, 'QC_stat/after_qc/final_results')
        clean_df = pd.read_table(clean_stat)
        clean_df = clean_df.rename(
            columns={"#Sample_ID": "Sample", "Total_Reads": "Clean reads", "Total_Bases": "Clean bases",
                     "Error%": "Error rate(%)", "Q20%": "Q20(%)", "Q30%": "Q30(%)", "GC%": "GC content(%)",
                     "read1Dup": "Dup_R1(%)", "read2Dup": "Dup_R2(%)", "PairedDup": "Dup_pair(%)",
                     "rRNA Ratio(%)": "rRNA ratio(%)"})
        clean_df = clean_df.set_index("Sample")
        if self.option('datatype') == 'rawdata':
            qc_df = pd.concat([raw_df, clean_df], axis=1)
        else:
            qc_df = clean_df
        def percet_change(df):
            try:
                df["Dup_R1(%)"] = df["Dup_R1(%)"] * 100
            except:
                pass
            try:
                df["Dup_R2(%)"] = df["Dup_R2(%)"] * 100
            except:
                pass
            try:
                df["Dup_pair(%)"] = df["Dup_pair(%)"] * 100
            except:
                pass
            return df
        qc_df = qc_df.apply(lambda x :percet_change(x),axis=1)
        qc_df.to_csv(os.path.join(self.target_dir, '05Basic_Analysis/01QC/sequencing_data_statistic.xls'),sep="\t")



        if "mapping" in self.analysis_content:
            # 02Align
            os.mkdir(os.path.join(self.target_dir,'05Basic_Analysis', '02Align'))
            align_stat = os.path.join(self.rnaseq_mapping.output_dir, 'Comparison_results')
            os.link(align_stat, os.path.join(self.target_dir, '05Basic_Analysis/02Align/align_stat.xls'))
            ## QualityAssessment
            os.makedirs(os.path.join(self.target_dir, '05Basic_Analysis/02Align/QualityAssessment'))
            chr_distribution = os.path.join(self.target_dir, '05Basic_Analysis/02Align/QualityAssessment', 'chr_distribution.xls')
            chr_stat_files = sorted(
                glob.glob("{}/*.bam_chr_stat.xls".format(self.map_assessment.output_dir + "/chr_stat")))
            dict_b = {}
            samples = {}
            chr = {}
            for file in chr_stat_files:
                sample_name = os.path.basename(file).split(".")[0]
                samples[sample_name] = 1
                with open(file, "r") as f:
                    f.readline()
                    for line in f:
                        items = line.strip().split("\t")
                        if items[0] not in chr:
                            chr[items[0]] = 1
                            dict_b[items[0]] = {}
                        dict_b[items[0]][sample_name] = items[1]
            with open(chr_distribution, "w") as w:
                w.write('Chromosome' + "\t" + "\t".join(sorted(samples.keys())) + "\n")
                for chr in sorted(chr.keys()):
                    w.write(chr)
                    for sample in sorted(samples.keys()):
                        if sample in dict_b[chr]:
                            w.write("\t" + dict_b[chr][sample])
                        else:
                            w.write("\t0")
                    w.write("\n")
            region_distribution = os.path.join(self.target_dir, '05Basic_Analysis/02Align/QualityAssessment', 'region_distribution.xls')
            region_stat_files = sorted(
                glob.glob("{}/*.reads_distribution.txt".format(self.map_assessment.output_dir + "/distribution")))
            dict_c = {}
            samples = {}
            distributions = ["CDS", "5'UTR", "3'UTR", "Introns", "Intergenic"]
            for file in region_stat_files:
                values = []
                sample_name = os.path.basename(file).split(".")[0]
                samples[sample_name] = 1
                with open(file, "r") as f:
                    f.readline()
                    f.next()
                    f.next()
                    f.next()
                    f.next()
                    for line in f:
                        if re.match(r"==", line):
                            continue
                        else:
                            line = line.strip().split()
                            values.append(float(line[2]))
                    values_new = values[:4]
                    values_new.append(sum([values[6], values[9]]))
                    total = sum(values_new)
                    for n, dis in enumerate(distributions):
                        if sample_name not in dict_c:
                            dict_c[sample_name] = {}
                        dict_c[sample_name][dis] = str(values_new[n]) + "(" + str(
                            float("%0.4f" % (values_new[n] / total)) * 100) + "%)"
            with open(region_distribution, "w") as w:
                w.write('Sample' + "\t" + "\t".join(distributions) + "\n")
                for sample in sorted(samples.keys()):
                    w.write(sample)
                    for dis in distributions:
                        w.write("\t" + dict_c[sample][dis])
                    w.write("\n")

        if "annotation" in self.analysis_content:
            # Assemble
            if self.option('is_assemble'):
                set_assemble = Upload()
                set_assemble.set_assembly(os.path.join(self.output_dir, 'assembly'), os.path.join(self.target_dir, '05Basic_Analysis','03Assemble'))
            trans2gene = os.path.join(self.target_dir, '05Basic_Analysis', '03Assemble','Sequence', 'trans2gene.txt')
            if not os.path.exists(os.path.join(self.intermediate_result, 'Transcripts', 'trans2gene.txt')):
                self.move_file(trans2gene, os.path.join(self.intermediate_result, 'Transcripts', 'trans2gene.txt'))
            # if os.path.exists(trans2gene):
            #     os.remove(trans2gene)

            '''
            #医学转录组不再需要单独注释模块
            # Annotation
            os.makedirs(os.path.join(self.target_dir, '05Annotation'))
            annot_file = RefAnnotation()
            merge_annot_v3 = os.path.join(self.output_dir, 'annot_merge')
            if "quantification" in self.analysis_content:
                gene_exp = self.quant.output_dir + "/gene.tpm.matrix"
                trans_exp = self.quant.output_dir + "/transcript.tpm.matrix"
            else:
                gene_exp = None
                trans_exp = None
            annot_dir = self.annot_merge.output_dir
            if not self.option("is_assemble"):
                self.api_annotation.has_new = False
                gene2trans = None
                gene2trans_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
            else:
                gene2trans = annot_dir + "/newannot_class/all_tran2gene.txt"
                gene2trans_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
            annot_file.run(merge_annot_v3, gene2trans, gene2trans_ref, gene_exp=gene_exp, trans_exp=trans_exp,
                           output_dir=os.path.join(self.target_dir, '05Annotation'))
            '''

            if self.option('is_assemble'):
                os.system("cat {} {} > {}".format(
                    self.known_pep,
                    os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.pep'),
                    os.path.join(self.target_dir, '05Basic_Analysis/03Assemble/Sequence/all_pep.fa')
                ))
                os.system("cat {} {} > {}".format(
                    self.known_cds,
                    os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.cds'),
                    os.path.join(self.target_dir, '05Basic_Analysis/03Assemble/Sequence/all_cds.fa')
                ))
                with open(os.path.join(self.target_dir, '05Basic_Analysis/03Assemble/Sequence/all_id.xls'), 'w') as fo, \
                        open(os.path.join(self.annot_merge.output_dir, "allannot_class/all_tran2gene.txt"), 'r') as fin:
                    fo.write("gene_id\ttranscript_id\tprotein_id\n")
                    for line in fin:
                        cols = line.strip("\n").split("\t")
                        fo.write("{}\t{}\t{}\n".format(cols[1], cols[0], cols[4]))

            else:
                os.makedirs(os.path.join(self.target_dir, '05Basic_Analysis','03Assemble'))
                os.makedirs(os.path.join(self.target_dir, '05Basic_Analysis','03Assemble', "Sequence"))
                os.link(self.known_cds,
                        os.path.join(self.target_dir, '05Basic_Analysis/03Assemble/Sequence/all_cds.fa'))
                os.link(self.known_pep,
                        os.path.join(self.target_dir, '05Basic_Analysis/03Assemble/Sequence/all_pep.fa'))

                with open(os.path.join(self.target_dir, '05Basic_Analysis/03Assemble/Sequence/all_id.xls'), 'w') as fo, \
                        open(os.path.join(self.annot_merge.output_dir, "refannot_class/all_tran2gene.txt"), 'r') as fin:
                    fo.write("gene_id\ttranscript_id\tprotein_id\n")
                    for line in fin:
                        cols = line.strip("\n").split("\t")
                        fo.write("{}\t{}\t{}\n".format(cols[1], cols[0], cols[4]))


        #02 Annotion&Express医学转录组02文件夹
        if "quantification" in self.analysis_content:
            # Expression
            os.makedirs(os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis'))
            if self.option('express_method') == 'RSEM':
                os.link(os.path.join(self.quant.output_dir, 'gene.count.matrix.annot.xls'),
                        os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.count.matrix.annot.xls'))
                if self.option('exp_way').lower() == 'tpm':
                    os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix.annot.xls'),
                            os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.TPM.matrix.annot.xls'))
                    #add by fwy 20210207 配合小工具需求
                    os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix'),
                            os.path.join(self.target_dir,
                                         '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.tpm.matrix.xls'))
                if self.option('exp_way').lower() == 'fpkm':
                    os.link(os.path.join(self.quant.output_dir, 'gene.fpkm.matrix.annot.xls'),
                            os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.FPKM.matrix.annot.xls'))
                    # add by fwy 20210207 配合小工具需求
                    os.link(os.path.join(self.quant.output_dir, 'gene.fpkm.matrix'),
                            os.path.join(self.target_dir,
                                         '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.fpkm.matrix.xls'))
                if self.option('level').lower() == 'transcript':
                    os.link(os.path.join(self.quant.output_dir, 'transcript.count.matrix.annot.xls'),
                            os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.count.matrix.annot.xls'))
                    if self.option('exp_way').lower() == 'tpm':
                        os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix.annot.xls'),
                                os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.TPM.matrix.annot.xls'))
                        # add by fwy 20210207 配合小工具需求
                        os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix'),
                                os.path.join(self.target_dir,
                                             '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.tpm.matrix.xls'))
                    if self.option('exp_way').lower() == 'fpkm':
                        os.link(os.path.join(self.quant.output_dir, 'transcript.fpkm.matrix.annot.xls'),
                                os.path.join(self.target_dir,
                                             '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.FPKM.matrix.annot.xls'))
                        # add by fwy 20210207 配合小工具需求
                        os.link(os.path.join(self.quant.output_dir, 'transcript.fpkm.matrix'),
                                os.path.join(self.target_dir,
                                             '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.fpkm.matrix'))
            else:
                os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix.annot.xls'),
                        os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.TPM.matrix.annot.xls'))
                # add by fwy 20210207 配合小工具需求
                os.link(os.path.join(self.quant.output_dir, 'gene.tpm.matrix'),
                        os.path.join(self.target_dir,
                                     '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.tpm.matrix.xls'))
                os.link(os.path.join(self.quant.output_dir, 'gene.count.matrix.annot.xls'),
                        os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.count.matrix.annot.xls'))
                if self.option('level').lower() == 'transcript':
                    # add by fwy 20210207 配合小工具需求
                    os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix'),
                            os.path.join(self.target_dir,
                                         '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.tpm.matrix.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.tpm.matrix.annot.xls'),
                            os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.TPM.matrix.annot.xls'))
                    os.link(os.path.join(self.quant.output_dir, 'transcript.count.matrix.annot.xls'),
                            os.path.join(self.target_dir, '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.count.matrix.annot.xls'))
            if self.option('sample_num') == 'multiple':
                # os.makedirs(os.path.join(self.target_dir, '06Express/ExpCorr'))
                # os.link(os.path.join(self.exp_corr.output_dir, 'sample_correlation.xls'),
                #         os.path.join(self.target_dir, '06Express/ExpCorr/sample_correlation.xls'))
                if self.option('group_table').prop['sample_number'] > 2:
                    os.makedirs(os.path.join(self.target_dir, '02Annotion&Express/02Exp_Corr'))
                    os.link(os.path.join(self.exp_corr.output_dir, 'sample_correlation.xls'),
                            os.path.join(self.target_dir, '02Annotion&Express/02Exp_Corr/sample_correlation.xls'))#modify by fwy 两个样本时不进行相关性分析
                    os.makedirs(os.path.join(self.target_dir, '02Annotion&Express/03Exp_PCA'))
                    # os.link(os.path.join(self.exp_pca.output_dir, 'PCA.xls'),
                    #         os.path.join(self.target_dir, '06Express/ExpPCA/PCA.xls'))
                    with open(os.path.join(self.exp_pca.output_dir, 'Explained_variance_ratio.xls'),"r") as raw:
                         with open(os.path.join(self.target_dir, '02Annotion&Express/03Exp_PCA/Explained_variance_ratio.xls'),"w") as new:
                              data = raw.read()
                              head = " \t Proportion of Variance\n"
                              new.write(head+data)
                    # os.link(os.path.join(self.exp_pca.output_dir, 'Explained_variance_ratio.xls'),
                    #         os.path.join(self.target_dir, '02Annotion&Express/03Exp_PCA/Explained_variance_ratio.xls'))

            os.makedirs(os.path.join(self.target_dir, '05Basic_Analysis', '04Exp_Distribution'))

        if 'other' in self.analysis_content and self.option('sample_num') == 'multiple':
            #表达差异和表达差异基因集相关模块
            #01DiffExpress
            self.gene_name_add = os.path.join(self.annot_merge.output_dir, "allannot_class", "all_annot_tran.xls")
            annot_info = pd.read_table(self.gene_name_add)
            gene2name_df = annot_info[["transcript_id","gene_id", "gene_name"]]
            # gene2name_df["gene_name"].fillna(gene2name_df["gene_id"], inplace=True)
            # self.gene_id2name_dict = dict(zip(gene2name_df.gene_id, gene2name_df.gene_name))
            self.id2namedict ={}
            # self.id2namedict["transcript_id"] = dict(zip(gene2name_df['transcript_id'], [x if x else '-' for x in gene2name_df['gene_name']]))
            self.id2namedict = dict(zip(gene2name_df['gene_id'], [x if x else '-' for x in gene2name_df['gene_name']]))
            add_name = GeneInfoSupple(id2namedict = self.id2namedict ,level ="G")
            if self.option('sample_num') == 'multiple':
                # 01DiffExpress_G
                if os.path.exists(os.path.join(self.target_dir,'01Diff_Express')):
                    shutil.rmtree(os.path.join(self.target_dir,'01Diff_Express'))
                os.makedirs(os.path.join(self.target_dir,'01Diff_Express', '01DiffExpress_G'))
                files = glob.glob(self.diffexp.output_dir + "/*.annot.xls") + glob.glob(
                    self.diffexp.output_dir + "/diff_summary_*.xls") + glob.glob(self.diffexp.output_dir + "/*.pdf")
                for b in files:
                    if not os.path.exists(os.path.join(self.target_dir, '01Diff_Express', '01DiffExpress_G', os.path.basename(b))):
                        os.link(b, os.path.join(self.target_dir, '01Diff_Express', '01DiffExpress_G', os.path.basename(b)))
                #02DiffExpress_ Cluster_Analysis
                os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '02DiffExpress_Cluster_Analysis'))
                selecet_geneset = os.listdir(os.path.join(self.diff_geneset_analysis.output_dir,"cluster"))[0]
                # os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '02DiffExpress_Cluster_Analysis',selecet_geneset))
                shutil.copytree(os.path.join(self.diff_geneset_analysis.output_dir,"cluster",selecet_geneset),
                                os.path.join(self.target_dir, '01Diff_Express','02DiffExpress_Cluster_Analysis',selecet_geneset))
                # 03DiffExpress_ Cluster_Analysis
                #获取基因集信息
                os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion'))
                diff_genesets = []
                for gt in os.listdir(os.path.join(self.diff_geneset_analysis.output_dir)):
                    if gt != "cluster":
                        diff_genesets.append(gt)
                #GO_Annotion
                os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion','01GO'))
                for geneset in diff_genesets:
                    if os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir, geneset,"diff_go_class")):
                        os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                                 '01GO',geneset))
                        # os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_go_class','go_class_table.xls'),
                        #         os.path.join(self.target_dir,'01Diff_Express', '03DiffExpress_Geneset_Annotion',
                        #                          '01GO',geneset,"go_class_stat.xls"))
                        #add by fwy 20210201 新增一列gene_name
                        # add_name = GeneInfoSupple(id2namedict=self.id2namedict, level="G")
                        raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_go_class','go_class_table.xls')
                        file_path = os.path.join(self.target_dir,'01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                                 '01GO',geneset,"go_class_stat.xls")
                        split = ";"
                        add_columns = [geneset+" list"]
                        # add_name.add_gene_name(file_path,split,"go",add_columns,)
                        try:
                            add_name.add_gene_name(raw_path, split, "go", add_columns,file_path )
                        except:
                            os.link(raw_path,file_path)

                #02KEGG_Annotion
                os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion', '02KEGG'))
                for geneset in diff_genesets:
                    if os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_kegg_class")):
                        os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion','02KEGG',geneset))
                        # os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_kegg_class','kegg_stat.xls'),
                        #         os.path.join(self.target_dir,'01Diff_Express', '03DiffExpress_Geneset_Annotion',
                        #                          '02KEGG',geneset,"kegg_class_stat.xls"))

                        # add by fwy 20210201 新增一列gene_name
                        raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_kegg_class','kegg_stat.xls')
                        file_path = os.path.join(self.target_dir,'01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                                 '02KEGG',geneset,"kegg_class_stat.xls")
                        split = ";"
                        kc = pd.read_table(raw_path, sep="\t", index_col=False)
                        add_columns = [x for x in kc.columns if x.endswith("genes")]
                        try:
                            add_name.add_gene_name(raw_path, split, "kegg", add_columns,file_path)
                        except:
                            os.link(raw_path, file_path)

                        # shutil.copytree(os.path.join(self.diff_geneset_analysis.output_dir, geneset,'diff_kegg_class','pathways'),
                        #                 os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                        #                              '02KEGG', geneset ,
                        #                             "pathways"))
                        os.link(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_class', 'pathways.tar.gz'),
                            os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                         '02KEGG', geneset,
                                         "pathways.tar.gz"))


                #03 Reactome _Annotion
                os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion', '03Reactome'))
                for geneset in diff_genesets:
                    if os.path.exists(os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_reactome_class")):
                        os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                                 '03Reactome', geneset))
                        # os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_class',
                        #                      'reactome_class.xls'),
                        #         os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                        #                      '03Reactome', geneset , "Reactome_class_stat.xls"))
                        # add by fwy 20210201 新增一列gene_name
                        raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_class',
                                             'reactome_class.xls')
                        file_path = os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                             '03Reactome', geneset , "Reactome_class_stat.xls")
                        split = ";"
                        rc = pd.read_table(raw_path, sep="\t", index_col=False)
                        add_columns = [x for x in rc.columns if x.endswith("genes")]
                        try:
                            add_name.add_gene_name(raw_path, split, "reactome", add_columns,file_path)
                        except:
                            os.link(raw_path,file_path)

                        try:
                            raw_path = os.path.join(self.diff_geneset_analysis, geneset, 'diff_reactome_class','reactome_path.xls')
                            file_path = os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion', '03Reactome', geneset, "pathway_class_stat.xls")
                        except Exception as e:
                            print(e)


                        # shutil.copytree(
                        #     os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_class', 'svg'),
                        #     os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                        #                  '03Reactome', geneset ,
                        #                   "Reactome_pathways"))
                        os.link(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_class', 'svg.tar.gz'),
                            os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                         '03Reactome', geneset,
                                         "Reactome_pathways.tar.gz"))


                #04DO_Annotion
                if self.organism_name == "Homo_sapiens":
                    os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                             '04DO'))
                    for geneset in diff_genesets:
                        if os.path.exists(
                                os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_do_class")):
                            os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                                     '04DO', geneset))
                            # os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_do_class',
                            #                      'do_level2_class.xls'),
                            #         os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                            #                      '04DO', geneset,
                            #                       "DO_class_stat.xls"))
                            # add by fwy 20210201 新增一列gene_name
                            raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_do_class',
                                                 'do_level2_class.xls')
                            file_path = os.path.join(self.target_dir, '01Diff_Express', '03DiffExpress_Geneset_Annotion',
                                                 '04DO', geneset,
                                                  "DO_class_stat.xls")
                            split = ";"
                            dc = pd.read_table(raw_path, sep="\t", index_col=False)
                            add_columns = [x for x in dc.columns if x.endswith("seqs")]
                            try:
                                add_name.add_gene_name(raw_path, split, "do", add_columns,file_path)
                            except:
                                os.link(raw_path,file_path)

                os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich'))
                # GO_Enrich
                os.makedirs(
                    os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich', '01GO'))
                for geneset in diff_genesets:
                    if os.path.exists(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_go_enrich")):
                        os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                                 '01GO', geneset ))
                        # os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_go_enrich',
                        #                      'go_enrich_geneset_list_gene.xls'),
                        #         os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                        #                      '01GO', geneset,  "go_enrich_stat.xls"))

                        # add by fwy 20210201 新增一列gene_name
                        raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_go_enrich',
                                             'go_enrich_geneset_list_gene.xls')
                        if os.path.exists(raw_path) :
                            file_path = os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                                 '01GO', geneset,  "go_enrich_stat.xls")
                            split = ";"
                            add_columns = ["seq_list"]
                            try:
                                add_name.add_gene_name(raw_path, split, "go", add_columns,file_path)
                            except:
                                    os.link(raw_path,file_path)

                # 02KEGG_Enrich
                os.makedirs(
                    os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich', '02KEGG'))
                for geneset in diff_genesets:
                    if os.path.exists(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_kegg_enrich")):
                        os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                                 '02KEGG', geneset))
                        # os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_enrich','enrich',
                        #                      '{}_gene.list.DE.list.check.kegg_enrichment.xls'.format(geneset)),
                        #         os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                        #                      '02KEGG', geneset , "kegg_enrich_stat.xls"))

                        # add by fwy 20210201 新增一列gene_name
                        raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_enrich','enrich',
                                             '{}_gene.list.DE.list.check.kegg_enrichment.xls'.format(geneset))
                        file_path = os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                             '02KEGG', geneset , "kegg_enrich_stat.xls")
                        split = "|"
                        add_columns = ["Genes"]
                        try:
                            add_name.add_gene_name(raw_path, split, "kegg_enrich", add_columns,file_path)
                        except:
                            os.link(raw_path,file_path)

                    # shutil.copytree(
                    #     os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_enrich','class','pathways'),
                    #     os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                    #                  '02KEGG', geneset,
                    #                    "kegg_enrich_pathways"))

                        os.link(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_kegg_enrich', 'class',
                                         'pathways.tar.gz'),
                            os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                         '02KEGG', geneset,
                                         "kegg_enrich_pathways.tar.gz"))
                #03 Reactome_Enrich
                os.makedirs(
                    os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich', '03Reactome'))
                for geneset in diff_genesets:
                    if os.path.exists(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_reactome_enrich")):
                        os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                                 '03Reactome', geneset ))
                        # os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_enrich',
                        #                      '{}_gene.list.reactome_enrichment.xls'.format(geneset)),
                        #         os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                        #                      '03Reactome', geneset , "Reactome_enrich_stat.xls"))

                        # add by fwy 20210201 新增一列gene_name
                        raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_enrich',
                                             '{}_gene.list.reactome_enrichment.xls'.format(geneset))
                        file_path = os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                             '03Reactome', geneset , "Reactome_enrich_stat.xls")
                        split = "|"
                        add_columns = ["Genes"]
                        try:
                            add_name.add_gene_name(raw_path, split, "reactome", add_columns,file_path)
                        except:
                            os.link(raw_path,file_path)


                        # shutil.copytree(
                        #     os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_enrich', 'svg'),
                        #     os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                        #                  '03Reactome', geneset ,
                        #                   "Reactome_enrich_pathways"))
                        os.link(
                            os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_reactome_enrich', 'svg.tar.gz'),
                            os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                         '03Reactome', geneset,
                                         "Reactome_enrich_pathways.tar.gz"))

                if self.organism_name == "Homo_sapiens":
                    # 04 DO_Enrich
                    os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich','04DO'))
                    for geneset in diff_genesets:
                        if os.path.exists(
                                os.path.join(self.diff_geneset_analysis.output_dir, geneset, "diff_do_enrich")):
                            os.makedirs(os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                                     '04DO', geneset ))
                            os.link(os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_do_enrich', 'do_enrichment.xls'),
                                    os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                                 '04DO', geneset , "do_enrich_stat.xls"))
                            # add by fwy 20210201 新增一列gene_name
                            raw_path = os.path.join(self.diff_geneset_analysis.output_dir, geneset, 'diff_do_enrich', 'do_enrichment.xls')
                            file_path = os.path.join(self.target_dir, '01Diff_Express', '04DiffExpress_Geneset_Enrich',
                                                 '04DO', geneset , "do_enrich_stat.xls")
                            split = "|"
                            add_columns = ["Genes"]
                            try:
                                add_name.add_gene_name(raw_path, split, "do", add_columns,file_path)
                            except:
                                pass
                                # os.link(raw_path,file_path)


            # SNP
            os.makedirs(os.path.join(self.target_dir, '03Gene_structure_analysis'))
            # AS
            self.logger.info('rmats_flag1')
            if self.option('sample_num') == 'multiple' and self.option('is_as') == 'True' and self.option("as_method").lower() == "rmats":
                self.logger.info('rmats_flag2')
                set_rmats = Upload()
                set_rmats.set_rmats(os.path.join(self.output_dir, 'altersplicing'), os.path.join(self.target_dir,'03Gene_structure_analysis', '01AS'),
                          self.task_id)
            elif self.option('is_as') == 'True' and self.option("as_method").lower() == "as_profile":
                shutil.copytree(self.as_analysis.output_dir,
                                os.path.join(self.target_dir, '03Gene_structure_analysis', '01AS'))

            self.logger.info('rmats_flag3')

            #SNP
            if self.option('sample_num') == 'multiple' and self.option('is_snp') == 'True':
                os.makedirs(os.path.join(self.target_dir, '03Gene_structure_analysis','02SNP_InDel_Analysis'))
                CopyFile().linkdir(os.path.join(self.work_dir, 'SnpTmp'), os.path.join(self.target_dir, '03Gene_structure_analysis','02SNP_InDel_Analysis'))
                rm_files = [os.path.join(self.target_dir, '03Gene_structure_analysis','02SNP_InDel_Analysis', 'snp_annotation.xls'),
                            os.path.join(self.target_dir,  '03Gene_structure_analysis','02SNP_InDel_Analysis','data_anno_pre.xls')]
                for rm_file in rm_files:
                    if os.path.exists(rm_file):
                        os.remove(rm_file)
                # SNP_vcf add by fwy 20200509
                os.makedirs(os.path.join(self.target_dir, '03Gene_structure_analysis','02SNP_InDel_Analysis', 'SNP_vcf'))
                if self.option('is_snp') == 'True':
                    if self.option('snp_method').lower() == 'samtools':
                        final_vcf = os.path.join(self.sam_rna.vcf_filter.output_dir, "final.vcf")
                    elif self.option('snp_method').lower() == 'gatk':
                        final_vcf = os.path.join(self.snp_rna.vcffilter.output_dir, "final.vcf")
                    else:
                        final_vcf = os.path.join(self.call_snp_indel.vcffilter.output_dir, "final.vcf")
                    os.link(final_vcf, os.path.join(self.target_dir, '03Gene_structure_analysis','02SNP_InDel_Analysis','SNP_vcf', os.path.basename(final_vcf)))
            # Somatic
            if self.option('sample_num') == 'multiple' and self.option('is_somatic') == 'True':
                # os.makedirs(os.path.join(self.target_dir, '03Gene_structure_analysis', '02Somatic_Analysis'))
                set_somatic = Upload()
                set_somatic.set_somatic(self.somatic_analysis.output_dir,
                                        os.path.join(self.work_dir, 'SomaticTmp'),
                                        os.path.join(self.target_dir, '03Gene_structure_analysis',
                                                     '02Somatic_Analysis'))

                # CopyFile().linkdir(self.somatic_analysis.output_dir,
                #                    os.path.join(self.target_dir, '03Gene_structure_analysis', '02Somatic_Analysis'))

            # Gene_Fusion
            if self.option('fq_type') == 'PE':
                os.makedirs(os.path.join(self.target_dir,'03Gene_structure_analysis', '03GeneFusion'))
                set_gene_fusion = Upload()
                set_gene_fusion.set_gene_fusion(self.gene_fusion.output_dir,os.path.join(self.target_dir,'03Gene_structure_analysis', '03GeneFusion','GeneFusion_Star_Fusion'))
                # shutil.copytree(self.gene_fusion.output_dir,os.path.join(self.target_dir,'03Gene_structure_analysis', '03GeneFusion', 'GeneFusion_Star_Fusion'))
        self.move_chart_file()
        self.def_upload()

    def run_parameter(self, dir_path):
        document = self.database['sg_table_relation'].find_one({})
        targets = document['target']
        parsed_tables = list()
        for target in targets:
            if target[0] in ['sg_status', 'sg_task', 'sg_software_para', 'sg_geneset', 'task_workflow_params',
                             'sg_one_touch', 'sg_file_dowload', 'sg_transcripts', 'sg_result_table_deposit',
                             'sg_splicing_rmats_count', 'sg_mapping', 'sg_qc'] or target[0] in parsed_tables:
                continue
            parsed_tables.append(target[0])
            records = self.database[target[0]].find({'task_id': self.task_id})
            if type(records) == dict:
                record = records
                if 'params' in record and record['params'] and record['params'] != 'none':
                    if type(record['params']) == dict or type(record['params']) == unicode:
                        params = json.loads(record['params'])
                    else:
                        params = record['params']
                    if 'submit_location' in params and params['submit_location']:
                        get_run_log = GetRunLog("medical_transcriptome", table=target[0], main_id=record['main_id'],
                                                dir_path=dir_path, append=True)
                        get_run_log.run()
            else:
                for record in records:
                    if 'params' in record and record['params'] and record['params'] != 'none':
                        if type(record['params']) == dict or type(record['params']) == unicode:
                            params = json.loads(record['params'])
                        else:
                            params = record['params']
                        if 'submit_location' in params and params['submit_location']:
                            get_run_log = GetRunLog("medical_transcriptome", table=target[0], main_id=record['main_id'],
                                                    dir_path=dir_path, append=True)
                            get_run_log.run()

    def move_chart_file(self):
        os.makedirs(os.path.join(self.target_dir, '02Annotion&Express/04Exp_Venn/'))
        os.makedirs(os.path.join(self.target_dir, '01Diff_Express/05DiffExp_Venn/'))
        file2uploads = [
            ("*raw_qc_qual.box*.pdf", '05Basic_Analysis/01QC/'),
            ("*raw_qc_error.line*.pdf", '05Basic_Analysis/01QC/'),
            ("*raw_qc_base.line*.pdf", '05Basic_Analysis/01QC/'),
            ("*clean_qc_qual.box*.pdf", '05Basic_Analysis/01QC/'),
            ("*clean_qc_error.line*.pdf", '05Basic_Analysis/01QC/'),
            ("*clean_qc_base.line*.pdf", '05Basic_Analysis/01QC/'),
            ("*align.satu*.pdf", '05Basic_Analysis/02Align/QualityAssessment/'),
            ("*align_pos_dist*.pdf", '05Basic_Analysis/02Align/QualityAssessment/'),
            ("*align_coverage*.pdf", '05Basic_Analysis/02Align/QualityAssessment/'),
            ("*align_chr_dist*.pdf", '05Basic_Analysis/02Align/QualityAssessment/'),
            # ("*all.assemble*.pdf", '04Assemble/'),
            # ("ref.assemble_relation*.pdf", '04Assemble/'),
            # # ("new.assemble_relation*.pdf", '04Assemble/'),
            # ("*annot_gene_stat*.pdf", '05Annotation/AnnotStatistics/'),
            ("all.exp_relation*.pdf", '02Annotion&Express/03Exp_PCA/'),
            ("all.exp.heat_corr.pdf", '02Annotion&Express/02Exp_Corr/'),
            ("all.exp.upset.pdf", '02Annotion&Express/04Exp_Venn/'),
            ("all.exp.venn.pdf", '02Annotion&Express/04Exp_Venn/'),
            # ("all.diffexp.*.pdf", '07DiffExpress_G/'),
            # ("*.diffexp.volcano.*.pdf", '07DiffExpress_G/'),
            # ("*.diffexp.scatter.*.pdf", '07DiffExpress_G/'),
            # ("all_*splice_stat.*.pdf", '09AS/'),
            ("*exp_distribution.box.pdf", '05Basic_Analysis/04Exp_Distribution/'),
            ("*exp_distribution.violin.pdf", '05Basic_Analysis/04Exp_Distribution/'),
            ("*snp.*stat.*.pdf", "03Gene_structure_analysis/02SNP_InDel_Analysis/"),
            ("*gene_fusion.pdf", "03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/")
        ]

        if self.option("control_file").is_set:
            cmp_list = self.option("control_file").prop["cmp_list"]
            for cmps in cmp_list:
                file2uploads.append((
                    "{}_vs_{}*splice_stat.*.pdf".format(cmps[0], cmps[1]), '03Gene_structure_analysis/01AS/{}_vs_{}_AS/'.format(cmps[0], cmps[1])))
        if self.option('sample_num') == 'multiple':
            file2uploads.append(("diff_genesets.analysis.venn.pdf", '01Diff_Express/05DiffExp_Venn/'))
            file2uploads.append(("*cluster*.pdf", '01Diff_Express/02DiffExpress_Cluster_Analysis/'))
        if self.option('sample_num') == 'multiple':
            file2uploads_geneset = [
                ("{}.go_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "01Diff_Express/03DiffExpress_Geneset_Annotion/01GO/{}".format("{geneset_name}")),
                ("{}.go_enrich.gene_set.*.pdf".format("{geneset_name}"),
                 "01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/{}".format("{geneset_name}")),
                ("{}*kegg_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/{}".format("{geneset_name}")),
                ("{}.kegg_enrich.gene_set.*.pdf".format("{geneset_name}"),
                 "01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/{}".format("{geneset_name}")),
                ("{}*reactome_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/{}".format("{geneset_name}")),
                ("{}.reactome_enrich.gene_set.*.pdf".format("{geneset_name}"),
                 "01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/{}".format("{geneset_name}")),
                ("{}.do_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "01Diff_Express/03DiffExpress_Geneset_Annotion/04DO/{}".format("{geneset_name}")),
                ("{}.do_enrich.gene_set.*.pdf".format("{geneset_name}"),
                 "01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/{}".format("{geneset_name}"))
            ]

        for filefrom, fileto in file2uploads:
            pdf_file = glob.glob(os.path.join(self.chart.work_dir, filefrom))
            for p in pdf_file:
                if os.path.exists(os.path.join(self.target_dir , fileto, os.path.basename(p))):
                    os.remove(os.path.join(self.target_dir, fileto, os.path.basename(p)))
                try:
                    os.link(p, os.path.join(self.target_dir, fileto, os.path.basename(p)))
                except:
                    self.logger.info('cuolacuolaraw:{} new:{}'.format(p,self.target_dir + "/" + fileto + "/" + os.path.basename(p)))
        if 'other' in self.analysis_content and self.option('sample_num') == 'multiple':

            genesets = [os.path.basename(i) for i in glob.glob("{}/*vs*".format(self.diff_geneset_analysis.output_dir))]

            #以基因集为单位,专门为差异一键化分析准备
            for geneset in genesets :
                for filefrom, fileto in file2uploads_geneset:
                    pdf_file = glob.glob((self.chart.work_dir + "/" + filefrom).format(geneset_name = geneset))
                    for p in pdf_file:
                        if os.path.exists(self.target_dir + "/" + fileto.format(geneset_name = geneset) + "/" + os.path.basename(p)):
                            os.remove(self.target_dir + "/" +  fileto.format(geneset_name = geneset) + "/" + os.path.basename(p))
                        try:
                            os.link(p, self.target_dir + "/" +  fileto.format(geneset_name = geneset) + "/" + os.path.basename(p))
                        except:
                            self.logger.info('cuolacuolaraw:{} new:{}'.format(p,
                                                                              self.target_dir + "/" + fileto.format(
                                                                                  geneset_name=geneset) + "/" + os.path.basename(
                                                                                  p)))

        pass




    @workfuncdeco
    def def_upload(self):
        # transfer = MultiFileTransfer()
        intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results')
        if not intermediate_dir.endswith("/"):
            intermediate_dir = intermediate_dir + "/"
        if not self.intermediate_result.endswith("/"):
            self.intermediate_result = self.intermediate_result + "/"
        self.upload_to_s3(self.intermediate_result, intermediate_dir)

        if self.option("report_img"):
            s3 = self._sheet.output.split(":")[0]
            report_img_dir = self.chart.work_dir + '/png/'
            report_img_s3 = s3 + "://commonbucket/files/report_img/medical/" + self.task_id + "/"
            self.upload_to_s3(report_img_dir, report_img_s3)

        # transfer.add_upload(self.intermediate_result, intermediate_dir)
        # transfer.perform()
        if self.option("upload_offline"):
            # 将质控结果文件和比对结果文件传输至线下服务器
            self.logger.info("开始向线下服务器传递比对bam文件，请耐心等待")
            cmd = "scp -r -i ~/.ssh/id_rsa {} dongmei.fu@192.168.10.46:/mnt/ilustre/centos7users/meng.luo/sanger_data".format(
                self.offline_results)
            code = os.system(cmd)
            if code == 0:
                self.logger.info("命令{}执行成功！".format(cmd))
            else:
                self.logger.info("命令{}执行失败！".format(cmd))
        if not self.target_dir.endswith("/"):
            self.target_dir = self.target_dir + "/"
        sdir = self.add_upload_dir(self.target_dir)
        sdir.add_regexp_rules([
            [r'03Align/AlignBam/.*\.bam', 'bam', '样本比对bam文件', 1],
            [r'01Diff_Express/01DiffExpress_G/.*_vs_.*\.xls', 'xls', '差异分析结果表', 0],
            [r'01Diff_Express/01DiffExpress_G/.*\.DE\.list', 'xls', '差异基因列表', 0,],
            [r'01Diff_Express/01DiffExpress_G/.*summary.*\.xls', 'xls', '差异表达基因统计表', 0],
            [r'01Diff_Express/01DiffExpress_G/.*total_diff_stat.*\.xls', 'xls', '差异表达基因详情总表', 0],
            [r'01Diff_Express/02DiffExpress_Cluster_Analysis/.*', '', '差异基因聚类分析', 0],
            [r'01Diff_Express/02DiffExpress_Cluster_Analysis/.*\.expression_matrix\.xls', '', '聚类热图分析表', 0],
            [r'01Diff_Express/02DiffExpress_Cluster_Analysis/.*\.subcluster_1\.xls', '', '子聚类分析表', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/01GO/.*', '', '差异基因集功能注释分析', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/01GO/.*/go_class_stat\.xls', '', 'GO分类统计表 ', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/01GO/.*/go_detail\.xls ', '', '基因/转录本对应GO注释详情表', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*down.*pdf', '', '下调基因KEGG注释统计图', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*up.*pdf', '', '上调基因KEGG注释统计图', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_class_stat\.xls', '', 'KEGG分类统计表', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_pathways', '', 'KEGG通路图', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_pathways.pathways.tar.gz', '', 'KEGG通路图', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_pathways/.*png', '', 'KEGG通路图片png', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_pathways/.*html', '', 'KEGG通路html文件', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*', '', '差异基因集功能注释分析', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_class_stat\.xls', '', 'Reactome分类统计表 ', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/pathway_class_stat\.xls', '', 'Reactome通路统计表 ', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_pathways', '','Reactome通路图 ', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_pathways.tar.gz', '', 'Reactome通路图 ', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_pathways/.*\.png', '', 'Reactome通路图片png', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_pathways/.*\.html', '', 'Reactome通路html文件 ', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/04DO/.*', '', '差异基因集功能注释分析', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/04DO/.*/DO_class_stat\.xls ', '', 'DO分类统计表 ', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/.*', '', '差异基因集功能富集分析', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/.*/go_enrich_stat\.xls', '', 'GO富集分析统计表  ', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/.*/go_detail\.xls ', '', '基因/转录本对应GO富集详情表', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_stat\.xls', '', 'KEGG富集分析统计表 ', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*', '', '差异基因集功能富集分析 ', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_pathways', '', 'KEGG富集通路图', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_pathways/.*png', '', 'KEGG通路图片png', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_pathways/.*html', '', 'KEGG通路html文件', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*', '', '差异基因集功能富集分析', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_enrich_stat\.xls', '',' Reactome富集分析统计表  ', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_pathways', '', 'Reactome通路图 ', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_pathways/.*\.png', '', 'Reactome通路图片png', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_pathways/.*\.html', '','Reactome通路html文件 ', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/.*', '', '差异基因集功能富集分析', 0],
            [r'01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/.*/DO_enrich_stat\.xls ', '', 'DO富集分析统计表 ', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*', '', '差异组别', 0, "211437"],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JC', '', 'JC定量', 0, "211438"],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JCEC', '', 'JCEC定量', 0, "211439"],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JC/SE.detail.xls', 'xls', 'SE可变剪切事件详情表（JC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JCEC/SE.detail.xls', 'xls', 'SE可变剪切事件详情表（JCEC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JC/MXE.detail.xls', 'xls', 'MXE可变剪切事件详情表（JC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JCEC/MXE.detail.xls', 'xls', 'MXE可变剪切事件详情表（JCEC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JC/A3SS.detail.xls', 'xls', 'A3SS可变剪切事件详情表（JC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JCEC/A3SS.detail.xls', 'xls', 'A3SS可变剪切事件详情表（JCEC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JC/A5SS.detail.xls', 'xls', 'A5SS可变剪切事件详情表（JC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JCEC/A5SS.detail.xls', 'xls', 'A5SS可变剪切事件详情表（JCEC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JC/RI.detail.xls', 'xls', 'RI可变剪切事件详情表（JC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/JCEC/RI.detail.xls', 'xls', 'RI可变剪切事件详情表（JCEC）', 0],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/diff_event_stats.xls', 'xls', '差异可变剪切事件统计表', 0 ],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/diff_pattern_stats.JC.xls', 'xls', '差异可变剪切模式变化统计表（JC）', 0 ],
            [r'03Gene_structure_analysis/01AS/.*_vs_.*/diff_pattern_stats.JCEC.xls', 'xls', '差异可变剪切模式变化统计表（JCEC）', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*', '', '单个样本融合信息详情', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*/fusion_inspector', '', 'fusion_inspector文件', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*/star-fusion.fusion_predictions.abridged.tsv', '','基因融合信息详情缩减表', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*/star-fusion.fusion_predictions.tsv','', '基因融合信息详情detail表', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*/fusion_inspector/finspector.fa', '', 'fusion序列文件', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*/fusion_inspector/finspector.bed', '', 'fusion序列文件', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*/fusion_inspector/finspector.junction_reads.bam','', 'fusion bam文件', 0],
            [r'03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion/.*/fusion_inspector/finspector.spanning_reads.bam','', 'fusion bam文件', 0],
            ## 图片描述
            [r"05Basic_Analysis/01QC/.*raw_qc_qual\.box\.pdf", 'pdf', '原始数据碱基质量分布图', 0],
            [r"05Basic_Analysis/01QC/.*raw_qc_error\.line\.pdf", 'pdf', '原始数据碱基错误率分布图', 0],
            [r"05Basic_Analysis/01QC/.*raw_qc_base\.line\.pdf", 'pdf', '原始数据碱基含量分布图', 0],
            [r"05Basic_Analysis/01QC/.*clean_qc_qual\.box\.pdf", 'pdf', '质控数据碱基质量分布图', 0],
            [r"05Basic_Analysis/01QC/.*clean_qc_error\.line\.pdf", 'pdf', '质控数据碱基错误率分布图', 0],
            [r"05Basic_Analysis/01QC/.*clean_qc_base\.line\.pdf", 'pdf', '质控数据碱基含量分布图', 0],
            [r"05Basic_Analysis/02Align/QualityAssessment/.*align\.satu.*\.pdf", 'pdf', '测序饱和度曲线图', 0],
            [r"05Basic_Analysis/02Align/QualityAssessment/.*align_pos_dist.*\.pdf", 'pdf', '不同区域Reads分布统计饼图', 0],
            [r"05Basic_Analysis/02Align/QualityAssessment/.*align_coverage.*\.pdf", 'pdf', '测序覆盖度分布图', 0],
            [r"05Basic_Analysis/02Align/QualityAssessment/.*align_chr_dist.*\.pdf", 'pdf', '不同染色体Reads分布统计柱状图', 0],
            # [r"04Assemble/all\.assemble.*\.column.pdf", 'pdf', '转录本长度分布柱状图', 0],
            # [r"04Assemble/all\.assemble.*\.pie.pdf", 'pdf', '转录本长度分布饼图', 0],
            # [r"04Assemble/.*assemble_relation.*\.columns\.pdf", 'pdf', '基因与转录本关系柱状图', 0],
            # [r"04Assemble/.*assemble_relation.*\.line\.pdf", 'pdf', '外显子与转录本关系曲线图', 0],
            # [r"05Annotation/AnnotStatistics/.*annot_gene_stat.*\.column\.pdf", 'pdf', '功能注释统计柱状图', 0],
            # [r"05Annotation/AnnotStatistics/.*annot_gene_stat.*\.venn\.pdf", 'pdf', '功能注释统计venn图', 0],
            # [r"05Annotation/AnnotStatistics/.*annot_gene_stat.*\.upset\.pdf", 'pdf', '功能注释统计upset图', 0],
            [r"05Basic_Analysis/04Exp_Distribution/.*exp_distribution\.box\.pdf", 'pdf', '表达量分布盒型图', 0],
            [r"05Basic_Analysis/04Exp_Distribution/.*exp_distribution\.density\.pdf", 'pdf', '表达量分布密度图', 0],
            [r"05Basic_Analysis/04Exp_Distribution/.*exp_distribution\.violin\.pdf", 'pdf', '表达量分布小提琴图', 0],
            [r"02Annotion&Express/04Exp_Venn/.*\.venn\.pdf", 'pdf', '样本间venn图', 0],
            [r"02Annotion&Express/04Exp_Venn/.*\.upset\.pdf", 'pdf', '样本间upset图', 0],
            [r"02Annotion&Express/03Exp_PCA/.*all\.exp_relation.*\.pdf", 'pdf', '样本间PCA图', 0],
            [r"02Annotion&Express/02Exp_Corr/.*all\.exp\.heat_corr.*\.pdf", 'pdf', '样本间相关性热图', 0],
            [r"01Diff_Express/01DiffExpress_G/all\.diffexp.*bar\.pdf", 'pdf', '表达量差异统计柱状图', 0],
            [r"01Diff_Express/01DiffExpress_G/all\.diffexp.*bar2\.pdf", 'pdf', '表达量差异统计堆积图', 0],
            [r"01Diff_Express/01DiffExpress_G/.*diffexp\.volcano.*\.pdf", 'pdf', '表达量差异火山图', 0],
            [r"01Diff_Express/01DiffExpress_G/.*scatter.*\.pdf", 'pdf', '表达量差异ma图', 0],
            [r"09AS/.*splice_stat.*\.pdf", 'pdf', '可变剪切事件统计图', 0],
            [r"03Gene_structure_analysis/01AS/.*_vs_.*/.*splice_stat\.*pie\.pdf", 'pdf', '组内差异可变剪切事件统计饼状图', 0],
            [r"03Gene_structure_analysis/01AS/.*_vs_.*/.*splice_stat\.*column\.pdf", 'pdf', '组内差异可变剪切事件统计柱状图', 0],
            [r"03Gene_structure_analysis/02SNP_InDel_Analysis/.*snp\.pos_stat\.pie\.pdf", 'pdf', 'SNP不同区域分布饼图', 0],
            [r"03Gene_structure_analysis/02SNP_InDel_Analysis/.*snp\.type_stat\.column\.pdf", 'pdf', 'SNP类型统计图', 0],
            [r"03Gene_structure_analysis/02SNP_InDel_Analysis/.*snp\.type_stat\.pie\.pdf", 'pdf', 'SNP类型饼图', 0],
            [r"03Gene_structure_analysis/02SNP_InDel_Analysis/.*snp\.depth_stat\.column\.pdf", 'pdf', 'SNP深度统计图', 0],
            [r"03Gene_structure_analysis/02SNP_InDel_Analysis/.*snp\.depth_stat\.pie\.pdf", 'pdf', 'SNP深度饼图', 0],
            [r"03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/.*.pdf", 'pdf', '基因融合circos图', 0],
            #差异基因集模块工作流部分图片
            [r"01Diff_Express/03DiffExpress_Geneset_Annotion/01GO/.*/.*go_annot\.gene_set\.column\.pdf", 'pdf', 'GO注释柱状图', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/.*/.*bar\.pdf", 'pdf','GO富集图片(柱形图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/.*/.*bar_line\.pdf", 'pdf', 'GO富集图片(柱形图-带折线)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/.*/.*buble\.pdf", 'pdf', 'GO富集图片(气泡图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/01GO/.*/.*buble2\.pdf", 'pdf', 'GO富集图片(气泡图-分散型)', 0],
            [r"01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*/.*all_1\.kegg_annot\.{}*genes\.column\.pdf", 'pdf',
             'kegg注释堆积图', 0],
            # [r"01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_annot\.{}*genes\.column\.pdf", 'pdf',
            #  'kegg注释柱状图', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/.*bar\.pdf", 'pdf','kegg富集相关图片(柱形图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/.*bar_line\.pdf", 'pdf','kegg富集相关图片(柱形图-带折线)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/.*buble\.pdf", 'pdf', 'kegg富集相关图片(气泡图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG/.*/.*buble2\.pdf", 'pdf','kegg富集相关图片(气泡图-分散型)', 0],
            # [r"01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/.*reactome_annot\.gene_set\.column.pdf", 'pdf',
            #  'reactome注释柱状图', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*up.*pdf', '', '上调基因reactome注释统计图', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*down.*pdf', '', '下调基因reactome注释统计图', 0],
            [r'01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*gz', '', 'kegg注释通路图', 0],
            [r"01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome/.*/.*all_1\.reactome_annot\.gene_set\.column.pdf",
             'pdf',
             'reactome注释堆积图', 0],
            # [r"01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/.*reactome_enrich\.gene_set.*\.pdf", 'pdf',
            #  'reactome富集相关图片', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/.*bar\.pdf", 'pdf','reactome富集相关图片(柱形图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/.*bar_line\.pdf", 'pdf', 'reactome富集相关图片(柱形图-带折线)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/.*buble\.pdf", 'pdf', 'reactome富集相关图片(气泡图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome/.*/.*buble2\.pdf", 'pdf', 'reactome富集相关图片(气泡图-分散型)', 0],
            [r"01Diff_Express/03DiffExpress_Geneset_Annotion/04DO/.*/.*do_annot\.gene_set\.column\.pdf", 'pdf',
             'DO注释柱状图', 0],
            # [r"01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/.*/.*do_enrich\.gene_set.*\.pdf", 'pdf',
            #  'DO富集相关图片', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/.*/.*bar\.pdf", 'pdf','DO富集相关图片(柱形图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/.*/.*bar_line\.pdf", 'pdf', 'DO富集相关图片(柱形图-带折线)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/.*/.*buble\.pdf", 'pdf', 'DO富集相关图片(气泡图)', 0],
            [r"01Diff_Express/04DiffExpress_Geneset_Enrich/04DO/.*/.*buble2\.pdf", 'pdf', 'DO富集相关图片(气泡图-分散型)', 0],
            [r"01Diff_Express/05DiffExp_Venn/diff_genesets\.analysis\.venn\.pdf", 'pdf',
             '差异基因集venn分析', 0],
            [r"01Diff_Express/02DiffExpress_Cluster_Analysis/.*cluster.*\.pdf", 'pdf',
             '差异基因集聚类分析图片', 0]
        ])
        sdir.add_relpath_rules([
            ['.', '', '流程分析结果目录', 0],
            ['01Diff_Express', '', '差异基因数据挖掘结果目录', 0],
            ['01Diff_Express/01DiffExpress_G', '', '表达差异统计分析', 0],
            ['01Diff_Express/02DiffExpress_Cluster_Analysis', '', '差异基因聚类分析', 0],
            ['01Diff_Express/03DiffExpress_Geneset_Annotion', '', '差异基因集功能注释分析', 0],
            ['01Diff_Express/03DiffExpress_Geneset_Annotion/01GO', '', '差异基因集GO功能注释分析', 0],
            ['01Diff_Express/03DiffExpress_Geneset_Annotion/02KEGG', '', '差异基因集KEGG功能注释分析', 0],
            ['01Diff_Express/03DiffExpress_Geneset_Annotion/03Reactome', '', '差异基因集Reactome功能注释分析', 0],
            ['01Diff_Express/03DiffExpress_Geneset_Annotion/04DO', '', '差异基因集DO功能注释分析', 0],
            ['01Diff_Express/04DiffExpress_Geneset_Enrich', '', '差异基因集功能富集分析', 0],
            ['01Diff_Express/04DiffExpress_Geneset_Enrich/01GO', '', '差异基因集GO功能富集', 0],
            ['01Diff_Express/04DiffExpress_Geneset_Enrich/02KEGG', '', '差异基因集KEGG富集分析分析', 0],
            ['01Diff_Express/04DiffExpress_Geneset_Enrich/03Reactome', '', '差异基因集Reactome富集分析分析', 0],
            ['01Diff_Express/04DiffExpress_Geneset_Enrich/04DO', '', '差异基因集DO功能富集', 0],
            ['01Diff_Express/05DiffExp_Venn', '', 'mRNA差异基因集Venn分析', 0],
            # ['02Annotion&Express', '', '功能注释与表达量分析数据挖掘结果目录', 0],
            ['02Annotion&Express', '', '所有基因数据挖掘结果目录', 0],
            ['02Annotion&Express/01Annot&Exp_Annalysis', '', '功能注释与表达量分析', 0],
            ['02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis', '', '表达量分析文件', 0],
            ['02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.count.matrix.annot.xls', '', '基因count表达定量注释结果表', 0],
            ['02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.count.matrix.annot.xls', '', '转录本count表达定量注释结果表', 0],
            ['02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.TPM.matrix.annot.xls', '', '基因TPM表达定量注释结果表', 0],
            ['02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.TPM.matrix.annot.xls', '', '转录本TPM表达定量注释结果表', 0],
            ['02Annotion&Express/02Exp_Corr', '', '样本间相关性分析', 0],
            ['02Annotion&Express/02Exp_Corr/sample_correlation.xls', '', '样本间相关性系数表', 0],
            ['02Annotion&Express/03Exp_PCA', '', '样本间PCA分析', 0],
            ['02Annotion&Express/03Exp_PCA/Explained_variance_ratio.xls', '', '主成分解释表', 0],
            ['02Annotion&Express/04Exp_Venn', '', '样本间Venn分析', 0],
            ['03Gene_structure_analysis', '', '基因结构分析数据挖掘结果目录', 0],
            ['03Gene_structure_analysis/01AS', '', '可变剪切分析结果目录', 0],
            ['03Gene_structure_analysis/01AS/AS_diff', '', '差异可变剪切分析文件', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis', '', 'SNP/InDel分析结果目录', 0],
            # ['03Gene_structure_analysis/02SNP_InDel_Analysis/SNP', '', 'SNP/InDel分析文件', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/SNP_vcf', '', 'SNP鉴定vcf结果目录', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/SNP_vcf/final.vcf', '', 'SNP鉴定vcf结果目录', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/snp_transition_tranversion_statistics.xls', '', 'SNP频率统计结果表', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/snp_depth_statistics.xls', '',  'SNP深度统计结果表', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/snp_position_distribution.xls', '', 'SNP不同区域分布结果表', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/indel_position_distribution.xls', '', 'InDel不同区域分布结果表', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/indel_anno.xls', '', 'InDel分析结果注释详情表', 0],
            ['03Gene_structure_analysis/02SNP_InDel_Analysis/snp_anno.xls', '', 'SNP分析结果注释详情表',
             0],
            ['03Gene_structure_analysis/02Somatic_Analysis', '', 'Somatic SNV/InDel分析结果目录', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/Somatic_SNV', '', 'Somatic SNV/InDel分析文件', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/SNV_vcf', '', 'SNV鉴定vcf结果目录', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/SNV_vcf/final.vcf', '', '', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/snv_transition_tranversion_statistics.xls', '', 'SNV类型统计结果表', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/snv_depth_statistics.xls', '', 'SNV深度统计结果表', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/snv_position_distribution.xls', '', 'SNV不同区域分布结果表', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/indel_position_distribution.xls', '','InDel不同区域分布结果表', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/indel_anno.xls', '', 'InDel分析结果注释详情表', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/snp_anno.xls', '', 'SNV分析结果注释详情表',   0],
            ['03Gene_structure_analysis/02Somatic_Analysis/Columns_of_dbNSFP_variant.xls', '', 'dbNSFP数据库SNV 保守性预测和致病性分析详情表',0],
            ['03Gene_structure_analysis/02Somatic_Analysis/dbNSFP_variant_stat.xls', '', 'dbNSFP数据库SNV 保守性预测和致病性分析统计表', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/predict.pdf', '', 'SNV的致病性预测分析图', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/predict.png', '','SNV的致病性预测分析图', 0],
            ['03Gene_structure_analysis/02Somatic_Analysis/predict.svg', '','SNV的致病性预测分析图', 0],
            ['03Gene_structure_analysis/03GeneFusion', '', '基因融合结果目录', 0],
            # ['03Gene_structure_analysis/03GeneFusion/', '', '', 0],
            ['03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/Star_fusion', '', '基因融合信息详情', 0],
            ['03Gene_structure_analysis/03GeneFusion/GeneFusion_Star_Fusion/fusion_stat.txt', '', '基因融合信息统计表', 0],
            ['04Background', 'xls', '项目背景结果目录', 0],
            ['04Background/sample_info.xls', 'xls', '样本信息表', 0],
            ['04Background/software_info.xls', 'xls', '软件信息表', 0],
            ['04Background/run_parameter.txt', 'txt', '分析参数日志', 0],
            ['05Basic_Analysis', '', '基础分析结果目录', 0],
            ['05Basic_Analysis/01QC', '', '测序数据质控', 0],
            ['05Basic_Analysis/01QC/sequencing_data_statistic.xls', '', '测序数据统计表', 0],
            ['05Basic_Analysis/02Align', '', '序列比对', 0],
            ['05Basic_Analysis/02Align/align_stat.xls', '', '比对结果统计表', 0],
            ['05Basic_Analysis/02Align/QualityAssessment', '', '比对结果整体评估', 0],
            ['05Basic_Analysis/02Align/QualityAssessment/chr_distribution.xls', '', '不同染色体Reads分布统计表', 0],
            ['05Basic_Analysis/02Align/QualityAssessment/region_distribution.xls', '', '不同区域Reads分布统计表', 0],
            ['05Basic_Analysis/03Assemble', '', '转录本组装', 0],
            ['05Basic_Analysis/03Assemble/classcode_statistics.xls', '', '新转录本类型统计表', 0],
            ['05Basic_Analysis/03Assemble/Sequence', '', '转录本序列文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/all_transcripts.fa', 'fasta', '组装结果序列文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/new_transcripts.fa', 'fasta', '新转录本序列文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/all_cds.fa', 'fasta', 'CDS序列文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/all_pep.fa', 'fasta', '蛋白序列文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/all_id.xls', 'xls', '基因转录本蛋白ID对应关系文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/all_transcripts.gtf', 'gtf', '组装结果GTF文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/new_transcripts.gtf', 'gtf', '新转录本GTF文件', 0],
            ['05Basic_Analysis/03Assemble/Sequence/trans2gene.txt', 'txt', '转录本与基因的对应关系文件', 0],
            ['05Basic_Analysis/04Exp_Distribution', '', '表达量分布', 0]

        ])
        self.update_collections()

    @workfuncdeco
    def update_collections(self):
        if self._sheet.output.endswith("/"):
            intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results')
        else:
            intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results/')
        db = Config().get_mongo_client(mtype='medical_transcriptome')[Config().get_mongo_dbname('medical_transcriptome')]
        col0 = db['sg_task']
        col0.update({'task_id': self.task_id},
                    {'$set': {'refrna_seqdb': os.path.join(intermediate_dir, 'SequenceDatabase/refrna_seqs.db')}},
                    upsert=True)
        col0.update({'task_id': self.task_id},
                    {'$set': {'refrna_seqdetail': os.path.join(intermediate_dir, 'SequenceDetail/')}},
                    upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'ref_gtf': self.ref_gtf}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'ref_genome': self.ref_genome}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'genome_id': self.genome_id}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'genome_id': self.genome_id}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'fq_type': self.option('fq_type')}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'strand_specific': self.option('strand_specific')}},
                    upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'strand_dir': self.option('strand_dir')}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'version': 'v1'}}, upsert=True)
        col0.update({'task_id': self.task_id}, {'$set': {'analysis_strategy': self.option('analysis_strategy')}},
                    upsert=True)
        # if "annotation" in self.analysis_content:
        #     # col0.update({'task_id': self.task_id},
        #     #             {'$set': {'database_version': {"kegg": self.option("kegg_version")}}}, upsert=True)


        if "annotation" in self.analysis_content:
            annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
        else:
            # 异常退出
            annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
        if self.annot_config_dict['kegg']['version'] > "2020":
            if self.option("kegg_specific") not in [None, ""]:
                annot_version_dict['kegg'] += "_spe"
        else:
            del annot_version_dict['kegg']

        col0.update({'task_id': self.task_id},
                    {'$set': {'database_version': annot_version_dict,
                              'annot_group': self.option("annot_group")}}, upsert=True)
        col0.update({'task_id': self.task_id},
                    {'$set': {'organism_name': self.organism_name}}, upsert=True)
        if self.option('is_assemble'):
            col0.update({'task_id': self.task_id},
                        {'$set': {'assemble_t2g': os.path.join(intermediate_dir, 'Transcripts/trans2gene.txt')}},
                        upsert=True)
            col0.update({'task_id': self.task_id},
                        {'$set': {
                            'assemble_fa': os.path.join(self._sheet.output,
                                                        '05Basic_Analysis/03Assemble/Sequence/all_transcripts.fa')}},
                        upsert=True)
            col0.update({'task_id': self.task_id},
                        {'$set': {
                            'as_gtf': os.path.join(self._sheet.output, '05Basic_Analysis/03Assemble/Sequence/all_transcripts.gtf')}},
                        upsert=True)
            col0.update({'task_id': self.task_id},
                        {'$set': {
                            'sample_gtf': os.path.join(intermediate_dir, 'Stringtie/')}},
                        upsert=True)
        else:
            col0.update({'task_id': self.task_id},
                        {'$set': {'assemble_t2g': os.path.join(intermediate_dir, 'Transcripts/trans2gene.txt')}},
                        upsert=True)
            col0.update({'task_id': self.task_id},
                                        {'$set': {'assemble_fa': os.path.join(intermediate_dir, 'Transcripts/all_transcripts.fa')}},
                                        upsert=True)
            col0.update({'task_id': self.task_id}, {'$set': {'as_gtf': self.ref_gtf}}, upsert=True)
        # if self.option('is_assemble'):
        #     col0.update({'task_id': self.task_id},
        #                 {'$set': {'assemble_t2g': os.path.join(intermediate_dir, 'Transcripts/trans2gene.txt')}},
        #                 upsert=True)
        #     col0.update({'task_id': self.task_id},
        #                 {'$set': {
        #                     'assemble_fa': os.path.join(self._sheet.output,
        #                                                 '04Assemble/Sequence/all_transcripts.fa')}},
        #                 upsert=True)
        #     col0.update({'task_id': self.task_id},
        #                 {'$set': {
        #                     'as_gtf': os.path.join(self._sheet.output, '05Basic_Analysis/03Assemble/Sequence/all_transcripts.gtf')}},
        #                 upsert=True)
        #     col0.update({'task_id': self.task_id},
        #                 {'$set': {
        #                     'sample_gtf': os.path.join(intermediate_dir, 'Stringtie/')}},
        #                 upsert=True)
        # else:
        #     col0.update({'task_id': self.task_id},
        #                 {'$set': {'assemble_t2g': os.path.join(intermediate_dir, 'Transcripts/trans2gene.txt')}},
        #                 upsert=True)
        #     col0.update({'task_id': self.task_id},
        #                 {'$set': {'assemble_fa': os.path.join(intermediate_dir, 'Transcripts/all_transcripts.fa')}},
        #                 upsert=True)
        #     col0.update({'task_id': self.task_id}, {'$set': {'as_gtf': self.ref_gtf}}, upsert=True)
        col1 = db['sg_annotation_stat']
        col1.update({'task_id': self.task_id},
                    {'$set': {'result_dir': os.path.join(intermediate_dir, 'Annotation')}}, upsert=True)
        if "quantification" in self.analysis_content:
            col2 = db['sg_exp']
            col2.update({'task_id': self.task_id, 'level': 'G'}, {
                '$set': {'count_file': os.path.join(intermediate_dir,
                                                    'Express/ExpAnnalysis/gene.count.matrix.xls')}}, upsert=True)
            if self.option('exp_way').lower() == 'tpm':
                col2.update({'task_id': self.task_id, 'level': 'G'}, {
                    '$set': {'exp_file': os.path.join(self._sheet.output,
                                                        '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.tpm.matrix.xls')}}, upsert=True)
            else:
                col2.update({'task_id': self.task_id, 'level': 'G'}, {
                    '$set': {'exp_file': os.path.join(self._sheet.output,
                                                      '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/gene.tpm.matrix.xls')}}, upsert=True)
            col2.update({'task_id': self.task_id, 'level': 'T'}, {
                '$set': {'count_file': os.path.join(intermediate_dir,
                                                    'Express/ExpAnnalysis/transcript.count.matrix.xls')}}, upsert=True)
            if self.option('exp_way').lower() == 'tpm':
                col2.update({'task_id': self.task_id, 'level': 'T'}, {
                    '$set': {'exp_file': os.path.join(self._sheet.output,
                                                        '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.tpm.matrix.xls')}}, upsert=True)
            else:
                col2.update({'task_id': self.task_id, 'level': 'T'}, {
                    '$set': {'exp_file': os.path.join(self._sheet.output,
                                                      '02Annotion&Express/01Annot&Exp_Annalysis/ExpAnnalysis/transcript.fpkm.matrix.xls')}}, upsert=True)
        super(MedicalTranscriptomeWorkflow, self).end()

    @workfuncdeco
    def export_task_info(self):
        api = self.api.api('task_info.medical_transcriptome')
        api.add_task_info()

    @workfuncdeco
    def export_genome_info(self):
        api = self.api.api('medical_transcriptome.genome_info')
        api.add_genome_info(file_path=self.genome_stat, species_name=self.option('ref_genome'), species=self.species,
                            ref_anno_version=self.genome_version, hyperlink=self.hyperlink)

    @workfuncdeco
    def export_sample_info(self):
        api = self.api.api('medical_transcriptome.qc')
        if self.option('datatype') == 'rawdata':
            fq_dir =self.option("fastq_dir").path
        else:
            fq_dir = self.option("qc_dir").pat
        sample_list = os.path.join(fq_dir, 'list.txt')
        if self.option('productive_table').is_set:
            api.add_sample_info(sample_list=sample_list, group_file=self.option('group_table').path,
                                bam_path=os.path.join(self.rnaseq_mapping.output_dir, "bam"),
                                productive_table=self.option('productive_table').path)

        else:
            api.add_sample_info(sample_list=sample_list, group_file=self.option('group_table').path,
                                bam_path=os.path.join(self.rnaseq_mapping.output_dir, "bam"))

    @workfuncdeco
    def export_qc_result(self):
        api = self.api.api('medical_transcriptome.qc')
        if self.option("datatype") == "rawdata":
            if self.option('sample_num') == 'multiple':
                if self.option("group_table").is_set:
                    self.group_id, specimen_names, category_names = api.add_sample_group(
                        self.option('group_table').path)
                if self.option('control_file').is_set:
                    self.control_id, compare_detail = api.add_group_compare(self.option('control_file').path,
                                                                            self.group_id)
                qc_id = api.add_qc(fq_type=self.option("fq_type"))
                api.add_qc_detail(qc_id, qc_stat_before = self.hiseq_reads_stat_raw.output_dir,qc_stat_after = self.hiseq_reads_stat_use.output_dir, group=self.option('group_table').path)
                api.add_qc_graph(qc_id, os.path.join(self.hiseq_reads_stat_raw.output_dir), "before")
                api.add_qc_graph(qc_id, os.path.join(self.hiseq_reads_stat_use.output_dir), "after")
            else:
                sp_set = set()
                fastq_dir = self.option('fastq_dir').path
                for line in open(os.path.join(fastq_dir, 'list.txt')):
                    sp_set.add(line.strip().split('\t')[1])
                if len(sp_set) != 1:
                    self.set_error('invalid sample number')
                sample = list(sp_set)[0]
                group_table = os.path.join(self.work_dir, 'group.txt')
                open(group_table, 'w').write('#sample\tgroup\n{}\t{}'.format(sample, sample))
                self.option("group_table").set_path(group_table)
                self.group_id, specimen_names, category_names = api.add_sample_group(group_table)
                qc_id = api.add_qc(fq_type=self.option("fq_type"))
                api.add_qc_detail(qc_id, qc_stat_before=self.hiseq_reads_stat_raw.output_dir,
                                  qc_stat_after=self.hiseq_reads_stat_use.output_dir,group=group_table)
                api.add_qc_graph(qc_id, os.path.join(self.hiseq_reads_stat_raw.output_dir), "before")
                api.add_qc_graph(qc_id, os.path.join(self.hiseq_reads_stat_use.output_dir), "after")
        else:
            if self.option('sample_num') == 'multiple':
                if self.option("group_table").is_set:
                    self.group_id, specimen_names, category_names = api.add_sample_group(
                        self.option('group_table').path)
                if self.option('control_file').is_set:
                    self.control_id, compare_detail = api.add_group_compare(self.option('control_file').path,
                                                                            self.group_id)
                qc_id = api.add_qc(fq_type=self.option("fq_type"))
                api.add_qc_detail(qc_id, qc_stat_after=self.hiseq_reads_stat_use.output_dir)
                api.add_qc_graph(qc_id, os.path.join(self.hiseq_reads_stat_use.output_dir, 'qualityStat'), "after")
            else:
                sp_set = set()
                fastq_dir = self.option('qc').path
                for line in open(os.path.join(fastq_dir, 'list.txt')):
                    sp_set.add(line.strip().split('\t')[1])
                if len(sp_set) != 1:
                    self.set_error('invalid sample number')
                sample = list(sp_set)[0]
                group_table = os.path.join(self.work_dir, 'group.txt')
                open(group_table, 'w').write('#sample\tgroup\n{}\t{}'.format(sample, sample))
                self.option("group_table").set_path(group_table)
                self.group_id, specimen_names, category_names = api.add_sample_group(group_table)
                qc_id = self.add_qc(fq_type=self.option("fq_type"))
                api.add_qc_detail(qc_id, qc_stat_after=self.hiseq_reads_stat_use.output_dir)
                api.add_qc_graph(qc_id, os.path.join(self.hiseq_reads_stat_use.output_dir, 'qualityStat'), "after")


    @workfuncdeco
    def export_ref_rna_qc_after(self):
        api = self.api.api('medical_transcriptome.qc')
        api.add_samples_info(qc_stat=self.hiseq_reads_stat_use.output_dir, fq_type=self.option('fq_type').lower(),
                             about_qc='after')
        api.add_gragh_info(quality_stat=os.path.join(self.hiseq_reads_stat_use.output_dir, 'qualityStat'),
                           about_qc='after')
        if self.option('sample_num') == 'multiple':
            if self.option("group_table").is_set:
                self.group_id, specimen_names, category_names = api.add_sample_group(self.option('group_table').path)
            if self.option('control_file').is_set:
                self.control_id, compare_detail = api.add_group_compare(self.option('control_file').path, self.group_id)

        else:
            sp_set = set()
            if self.option('datatype') == 'rawdata':
                fastq_dir = self.option('fastq_dir').path
            else:
                fastq_dir = self.option('qc_dir').path
            for line in open(os.path.join(fastq_dir, 'list.txt')):
                sp_set.add(line.strip().split('\t')[1])
            if len(sp_set) != 1:
                self.set_error('invalid sample number')
            sample = list(sp_set)[0]
            group_table = os.path.join(self.work_dir, 'group.txt')
            open(group_table, 'w').write('#sample\tgroup\n{}\t{}'.format(sample, sample))
            self.group_id, specimen_names, category_names = api.add_specimen_group(group_table)
        intermediate_dir = self._sheet.output.replace('workflow_results', 'intermediate_results/')
        api.add_bam_path(intermediate_dir)

    @workfuncdeco
    def export_ref_rna_qc_alignment(self):
        api = self.api.api('medical_transcriptome.mapping')
        api.add_mapping_stat(os.path.join(self.rnaseq_mapping.output_dir, 'stat'), method=self.option('align_method').lower(), group=self.option('group_table').path)

    @workfuncdeco
    def export_ref_rna_qc_assessment(self):
        api = self.api.api('medical_transcriptome.mapping')
        if 'saturation' in self.option('map_assess_method'):
            api.add_rpkm_table(os.path.join(self.map_assessment.output_dir, 'saturation'), self.option('group_table').path)
        if 'coverage' in self.option('map_assess_method'):
            api.add_coverage_table(os.path.join(self.map_assessment.output_dir, 'coverage'), self.option('group_table').path)
        if 'distribution' in self.option('map_assess_method'):
            api.add_distribution_table(os.path.join(self.map_assessment.output_dir, 'distribution'), self.option('group_table').path)
        if 'chr_stat' in self.option('map_assess_method'):
            api.add_chrom_distribution_table(os.path.join(self.map_assessment.output_dir, 'chr_stat'), self.option('group_table').path)

    @workfuncdeco
    def export_ref_assembly(self):
        api = self.api.api('medical_transcriptome.ref_assembly')
        params = json.dumps({'task_id': self.task_id, 'submit_location': 'transcripts', 'task_type': 2}, sort_keys=True)
        if self.option('assemble_method').lower() == 'stringtie':
            all_gtf_path = os.path.join(self.refrna_assemble.output_dir, 'Stringtie')
            merged_path = os.path.join(self.refrna_assemble.output_dir, 'StringtieMerge')
        elif self.option('assemble_method').lower() == 'cufflinks':
            all_gtf_path = os.path.join(self.refrna_assemble.output_dir, 'Cufflinks')
            merged_path = os.path.join(self.refrna_assemble.output_dir, 'Cuffmerge')
        statistics_path = os.path.join(self.refrna_assemble.output_dir, 'Statistics')
        api.add_assembly_result(params=params, all_gtf_path=all_gtf_path, merged_path=merged_path,
                                statistics_path=statistics_path)

    @workfuncdeco
    def export_annotation(self):
        self.api_annotation = self.api.api("medical_transcriptome.ref_annotation")
        annot_dir = self.annot_merge.output_dir
        # annot_dir = os.path.join(self.work_dir,"AnnotMerge__2/output")
        if not self.option("is_assemble"):
            self.api_annotation.has_new = False
            trans2gene = None
            trans2gene_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
        else:
            trans2gene = annot_dir + "/newannot_class/all_tran2gene.txt"
            trans2gene_ref = annot_dir + "/refannot_class/all_tran2gene.txt"
        self.api_annotation.species_name = self.option("ref_genome")
        self.api_annotation.has_new = self.option('is_assemble')
        params_dict = {
            "nr_evalue": str(self.option("nr_evalue")),
            "nr_similarity": self.option("nr_similarity"),
            "nr_identity": self.option("nr_identity"),
            "swissprot_evalue": str(self.option("swissprot_evalue")),
            "swissprot_similarity": self.option("swissprot_similarity"),
            "swissprot_identity": self.option("swissprot_identity"),
            "cog_evalue": str(self.option("cog_evalue")),
            "cog_similarity": self.option("cog_similarity"),
            "cog_identity": self.option("cog_identity"),
            "kegg_evalue": str(self.option("kegg_evalue")),
            "kegg_similarity": self.option("kegg_similarity"),
            "kegg_identity": self.option("kegg_identity"),
            "pfam_evalue": str(self.option("pfam_evalue")),
        }
        if "quantification" in self.analysis_content:
            gene_exp = self.quant.output_dir + "/gene.tpm.matrix"
            trans_exp = self.quant.output_dir + "/transcript.tpm.matrix"
        else:
            gene_exp = None
            trans_exp = None
        self.api_annotation.run(
            annot_dir,
            trans2gene,
            trans2gene_ref,
            params_dict=params_dict,
            taxonomy=self.option("kegg_database"),
            exp_level=self.option("level").lower(),
            version="v1",
            gene_exp=gene_exp,
            trans_exp=trans_exp,
        )

        '''
        api = self.api.api('ref_rna_v2.annotation')
        api.species_name = self.option('ref_genome')
        annot_merge_output_dir = self.annot_merge.output_dir
        params_dict = {
            'task_id': self.task_id,
            'submit_location': 'annotationstat',
            'task_type': 2,
            'nr_evalue': self.option('nr_evalue'),
            'nr_identity': self.option('nr_identity'),
            'nr_similarity': self.option('nr_similarity'),
            'swissprot_evalue': self.option('swissprot_evalue'),
            'swissprot_identity': self.option('swissprot_identity'),
            'swissprot_similarity': self.option('swissprot_similarity'),
            'cog_evalue': self.option('cog_evalue'),
            'cog_identity': self.option('cog_identity'),
            'cog_similarity': self.option('cog_similarity'),
            'kegg_evalue': self.option('kegg_evalue'),
            'kegg_identity': self.option('kegg_identity'),
            'kegg_similarity': self.option('kegg_similarity'),
            'pfam_evalue': self.option('pfam_evalue')
        }
        taxonomy = self.option('kegg_database')
        exp_level = self.option('level').lower()
        api.has_new = self.option('is_assemble')
        api.run(annot_merge_output_dir, params_dict, taxonomy, exp_level)
        '''

    @workfuncdeco
    def export_all_exp_matrix(self):
        api = self.api.api('medical_transcriptome.all_exp')
        exp_type = self.option('exp_way').upper()
        exp_matrix = {'T': os.path.join(self.quant.output_dir, 'transcript.{}.matrix'.format(exp_type.lower())),
                      'G': os.path.join(self.quant.output_dir, 'gene.{}.matrix'.format(exp_type.lower()))}
        quant_method = self.option('express_method')
        lib_type = self.option('strand_dir')
        if self.option('sample_num') == 'multiple':
            group_dict = self.option('group_table').prop['group_dict']
        else:
            group_table = os.path.join(self.work_dir, 'group.txt')
            for line in open(group_table):
                if line.strip() and line[0] != '#':
                    sample = line.strip().split('\t')[1]
            group_dict = {sample: [sample]}
        group_id = self.group_id
        project_sn = self.project_sn
        task_id = self.task_id
        def params(exp_level):
            return json.dumps({
                'task_id': task_id,
                'submit_location': 'exp_detail',
                'task_type': 2,
                'method': quant_method,
                'exp_type': exp_type,
                'level' :exp_level
            }, sort_keys=True, separators=(',', ':'))

        self.exp_ids = dict()
        if self.option('level').lower() == 'transcript':
             annot_matrix_tran = os.path.join(self.annot_merge.output_dir,"allannot_class","all_annot_tran.xls")
             # annot_dir = os.path.join(self.work_dir, "AnnotMerge__2/output/allannot_class")
             # annot_matrix_tran = os.path.join(annot_dir,"all_annot_tran.xls")
             self.exp_ids['T'] = api.add_exp(exp_matrix['T'],annot_matrix_tran, quant_method=quant_method, level='T',
                                             lib_type=lib_type, group_dict=group_dict, group_id=group_id,
                                             exp_type=exp_type, add_distribution=False, project_sn=project_sn,
                                             task_id=task_id, params=params("T"))
        annot_matrix_gene = os.path.join(self.annot_merge.output_dir, "all_annot_gene.xls")
        # annot_dir = os.path.join(self.work_dir, "AnnotMerge__2/output/allannot_class")
        annot_matrix_gene = os.path.join(self.annot_merge.output_dir,"allannot_class", "all_annot_gene.xls")
        self.exp_ids['G'] = api.add_exp(exp_matrix['G'],annot_matrix_gene ,quant_method=quant_method, level='G',
                                        lib_type=lib_type, group_dict=group_dict, group_id=group_id,
                                        exp_type=exp_type, add_distribution=False, project_sn=project_sn,
                                        task_id=task_id, params=params("G"))

    @workfuncdeco
    def export_all_exp_distribution(self):
        api = self.api.api('medical_transcriptome.all_exp')
        task_id = self.task_id
        quant = self.quant
        exp_way = self.option('exp_way')

        def exp_matrix(exp_level):
            return quant.option('ref_{}_{}'.format(exp_level, exp_way)).path

        group_dict = self.option('group_table').prop['group_dict']
        exp_ids = self.exp_ids
        group_id = self.group_id

        def params(exp_level):
            return json.dumps({
                'task_id': task_id,
                'submit_location': 'expgraph',
                'task_type': 2,
                'exp_id': str(exp_ids[exp_level]),
                'group_dict': group_dict,
                'group_id': str(group_id),
                'level': exp_level,
                'kind': 'ref'
            }, sort_keys=True, separators=(',', ':'))

        quant_method = self.option('express_method')
        project_sn = self.project_sn
        if self.option('level').lower() == 'transcript':
            api.add_distribution(exp_matrix=exp_matrix('transcript'), group_dict=group_dict, params=params('T'),
                                 level='T', quant_method=quant_method, project_sn=project_sn, task_id=task_id)
        api.add_distribution(exp_matrix=exp_matrix('gene'), group_dict=group_dict, params=params('G'), level='G',
                             quant_method=quant_method, project_sn=project_sn, task_id=task_id)

    @workfuncdeco
    def export_add_exp_venn(self):
        api = self.api.api('medical_transcriptome.all_exp')
        graph_table = os.path.join(self.exp_venn.output_dir, 'venn_graph.xls')
        group_dict = self.option('group_table').prop['group_dict']
        if len(group_dict) > 6:
            group_dict = OrderedDict(group_dict.items()[:6])
        params = json.dumps(dict(
            task_id=self.task_id,
            submit_location='expvenn',
            task_type=2,
            exp_id=str(self.exp_ids['G']),
            group_id=str(self.group_id),
            level='G',
            group_dict=group_dict,
            threshold='1',
            kind='ref',
        ), sort_keys=True, separators=(',', ':'))
        import datetime
        time_now = datetime.datetime.now()
        name = 'ExpVenn_G_{}_{}_{}'.format(
            self.option('express_method'), self.option('exp_way').upper(), time_now.strftime('%Y%m%d_%H%M%S'))
        main_info = dict(
            project_sn=self.project_sn,
            task_id=self.task_id,
            version='v1',
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=str(self.exp_ids['G']),
            desc='Expression venn analysis main table',
            params=params,
            status='start'
        )
        main_id = api.create_db_table('sg_exp_venn', [main_info])
        api.add_exp_venn(graph_table, main_id=main_id)

    @workfuncdeco
    def export_all_exp_corr(self):
        api = self.api.api('medical_transcriptome.all_exp')
        corr_work_dir = self.exp_corr.work_dir
        quant_method = self.option('express_method')
        task_id = self.task_id
        exp_ids = self.exp_ids
        group_dict = self.option('group_table').prop['group_dict']
        group_id = self.group_id

        def params(exp_level):
            return json.dumps({
                'task_id': task_id,
                'submit_location': 'expcorr',
                'task_type': 2,
                'exp_id': str(exp_ids[exp_level]),
                'group_dict': group_dict,
                'group_id': str(group_id),
                'level': exp_level,
                'scm': 'complete',
                'scd': 'euclidean',
                'corr_method': 'pearson',
                'kind': 'ref',
                'Draw_in_groups':'no'
            }, sort_keys=True, separators=(',', ':'))

        project_sn = self.project_sn
        api.add_exp_corr2(corr_work_dir, exp_level='G', quant_method=quant_method, params=params('G'),
                          project_sn=project_sn, task_id=task_id)

    @workfuncdeco
    def export_all_exp_pca(self):
        api = self.api.api('medical_transcriptome.all_exp')
        pca_output_dir = self.exp_pca.output_dir
        quant_method = self.option('express_method')
        task_id = self.task_id
        exp_ids = self.exp_ids
        exp_way = self.option('exp_way')
        group_dict = self.option('group_table').prop['group_dict']
        group_id = self.group_id

        def params(exp_level):
            return json.dumps({
                'task_id': task_id,
                'submit_location': 'exppca',
                'task_type': 2,
                'exp_id': str(exp_ids[exp_level]),
                'group_dict': group_dict,
                'group_id': str(group_id),
                'level': exp_level,
                'kind': 'ref',
                'Draw_in_groups': "no"
            }, sort_keys=True, separators=(',', ':'))

        project_sn = self.project_sn
        main_id = api.add_exp_pca2(pca_output_dir, quant_method=quant_method, exp_id=exp_ids['G'], exp_level='G',
                                   params=params('G'), project_sn=project_sn, task_id=task_id)
        if hasattr(self, 'ellipse'):
            api.insert_ellipse_table(os.path.join(self.ellipse.work_dir, 'ellipse_out.xls'), main_id)

    @workfuncdeco
    def export_all_exp_diff(self):
        api = self.api.api('medical_transcriptome.all_exp')
        diff_output = self.diffexp.output_dir
        exp_ids = self.exp_ids
        group_dict = self.option('group_table').prop['group_dict']
        group_id = self.group_id
        quant_method = self.option('express_method')
        diff_method = self.option('diff_method')
        project_sn = self.project_sn
        task_id = self.task_id
        control_id = self.control_id
        fc = str(float(self.option('fc')))
        # if '.' in fc:
        #     if fc.split('.')[1] == '0':
        #         fc = str(int(float(fc)))
        correct_method = self.option('padjust_way')
        stat_type = self.option('pvalue_padjust')
        stat_cutoff = str(self.option('diff_fdr_ci'))
        # tpm_filter_threshold = str(float(self.option('filter_tpm')))

        # if '.' in tpm_filter_threshold:
        #     if tpm_filter_threshold.split('.')[1] == '0':
        #         tpm_filter_threshold = str(int(float(tpm_filter_threshold)))

        def params(exp_level):
            if diff_method.lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                params_dict = {
                    'task_id': task_id,
                    'submit_location': 'diffgene',
                    'task_type': 2,
                    'exp_id': str(exp_ids[exp_level]),
                    'group_id': str(group_id),
                    'control_id': str(control_id),
                    'level': exp_level,
                    'group_dict': group_dict,
                    'fc': fc,
                    # 'tpm_filter_threshold': tpm_filter_threshold,
                    'filter_method':"no",
                    'stat_type': stat_type,
                    'stat_cutoff': stat_cutoff,
                    'diff_method': diff_method,
                    'kind': 'ref',
                    'is_batch': 'False',
                }
                if stat_type == 'padjust':
                    params_dict.update({'correct_method': correct_method})
            else:
                params_dict = {
                    'task_id': task_id,
                    'submit_location': 'diffgene',
                    'task_type': 2,
                    'exp_id': str(exp_ids[exp_level]),
                    'group_id': str(group_id),
                    'control_id': str(control_id),
                    'level': exp_level,
                    'group_dict': group_dict,
                    'fc': fc,
                    'filter_method': 'no',
                    # 'stat_cutoff': stat_cutoff,
                    'diff_method': diff_method,
                    'kind': 'ref',
                    'is_batch': 'False',
                    # 'prob': "0.8",
                    'prob': stat_cutoff,

                }
            return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

        if diff_method.lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            self.diff_id = api.add_diffexp(diff_output, exp_id=exp_ids['G'], group_dict=group_dict, group_id=group_id, exp_level='G',
                            quant_method=quant_method, diff_method=diff_method, project_sn=project_sn, task_id=task_id,
                            params=params('G'), pvalue_padjust=stat_type)
        else:
            self.diff_id = api.add_diffexp_noiseq(diff_output, exp_id=exp_ids['G'], group_dict=group_dict, group_id=group_id,
                                   exp_level='G',
                                   quant_method=quant_method, diff_method=diff_method, project_sn=project_sn,
                                   task_id=task_id,
                                   params=params('G'))

    @workfuncdeco
    def export_gene_detail(self):
        api = self.api.api('medical_transcriptome.gene_detail')
        refrna_seqdb = self.detail.option('database').path
        if self.option('is_assemble'):
            t2g_file = self.refrna_assemble.option('trans2gene').path
            txpt_fa = self.refrna_assemble.option('all_transcripts_fa').path
            new_cds = os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.cds')
            new_pep = os.path.join(self.annot_orfpfam.output_dir, 'new_transcripts.fa.transdecoder.pep')
        else:
            t2g_file = self.transcript_abstract.option('trans2gene').path
            txpt_fa = self.transcript_abstract.option('trans_fa').path
            new_cds = None
            new_pep = None
        txpt_bed = self.gene_fa.option('transcript_bed').path
        gene_bed = self.gene_fa.option('gene_bed').path
        gene_fa = self.gene_fa.option('gene_fa').path
        biomart_file = self.des
        biomart_type = self.des_type
        species_urls = self.hyperlink
        api.add_gene_detail(refrna_seqdb, t2g_file, txpt_bed, txpt_fa, gene_bed, gene_fa,
                            biomart_file, biomart_type, species_urls, new_cds, new_pep)

        gene_stat  =os.path.join(self.detail.output_dir,"detail","gene_stat")
        trans_stat = os.path.join(self.detail.output_dir, "detail", "trans_stat")
        api_seq = self.api.api('medical_transcriptome.seq_detail')
        api_seq.add_seq_stat(gene_stat, trans_stat)




    @workfuncdeco
    def export_snp(self):
        api = self.api.api('medical_transcriptome.snp')
        task_id = self.task_id
        project_sn = self.project_sn
        new_output = os.path.join(self.work_dir, 'SnpTmp')
        if os.path.exists(new_output):
            shutil.rmtree(new_output)
        os.mkdir(new_output)
        if self.option('snp_method').lower() == 'gatk':
            snp_anno = self.snp_rna.output_dir
            if os.path.exists(snp_anno + "/snp_annotation_statistics.xls"):
                params = dict(
                    task_id=task_id,
                    submit_location='snp',
                    task_type=2,
                    method_type='gatk',
                )
                api.add_snp_main(snp_anno=snp_anno, group=self.option('group_table').path, params=params, task_id=task_id, method_type='gatk',
                                 project_sn=project_sn, new_output=new_output)
        if self.option('snp_method').lower() == 'samtools':
            snp_anno = self.sam_rna.output_dir
            if os.path.exists(snp_anno + "/snp_annotation_statistics.xls"):
                params = dict(
                    task_id=task_id,
                    submit_location='snp',
                    task_type=2,
                    method_type='samtools'
                )
                api.add_snp_main(snp_anno=snp_anno, group=self.option('group_table').path, params=params, task_id=task_id, method_type='samtools',
                                 project_sn=project_sn, new_output=new_output)
        if self.option('snp_method').lower() == 'sentieon':
            snp_anno = self.call_snp_indel.output_dir + '/predeal'
            if os.path.exists(snp_anno + "/snp_annotation_statistics.xls"):
                params = dict(
                    task_id=task_id,
                    submit_location='snp',
                    task_type=2,
                    method_type='sentieon',
                    algorithm="HaplotypeCaller"
                )
                api.add_snp_main(snp_anno=snp_anno, group=self.option('group_table').path, params=params, task_id=task_id, method_type='sentieon',
                                 project_sn=project_sn, new_output=new_output)
        # if os.path.exists(os.path.join(self.work_dir, 'SnpTmp/snp_anno.xls')):
        #     os.remove(os.path.join(self.work_dir, 'SnpTmp/snp_anno.xls'))
        # os.link(os.path.join(snp_anno, 'data_anno_pre.xls'), os.path.join(self.work_dir, 'SnpTmp/snp_anno.xls'))

    @workfuncdeco
    def export_somatic(self):
        api = self.api.api('medical_transcriptome.somatic')
        task_id = self.task_id
        project_sn = self.project_sn
        new_output = os.path.join(self.work_dir, 'SomaticTmp')
        if os.path.exists(new_output):
            shutil.rmtree(new_output)
        os.mkdir(new_output)
        snp_anno = self.somatic_analysis.output_dir + '/predeal'
        if os.path.exists(snp_anno + "/snp_annotation_statistics.xls"):
            params = dict(
                task_id=task_id,
                submit_location='somatic',
                task_type=2,
                method_type='sentieon'
            )
            predict_stat = os.path.join(self.somatic_analysis.somatic_predict.predict_stat.output_dir,"predict_stat")
            api.add_somatic_main(snp_anno=snp_anno, params=params, task_id=task_id, method_type='sentieon',
                             project_sn=project_sn, predict_stat=predict_stat,new_output=new_output)


    @workfuncdeco
    def export_rmats(self):
        api = self.api.api('medical_transcriptome.rmats')
        for p in glob.glob(os.path.join(self.as_analysis.output_dir, '*')):
            if os.path.isdir(p):
                outpath = p
                ctrl, test = os.path.basename(outpath).split('_vs_')
                group_dict = {ctrl: self.option('group_table').prop['group_dict'][ctrl],
                              test: self.option('group_table').prop['group_dict'][test]}
                compare_plan = '{}|{}'.format(ctrl, test)
                params = json.dumps({
                    'task_id': self.task_id,
                    'submit_location': 'splicingrmats',
                    'task_type': 2,
                    'group_id': str(self.group_id),
                    'group_dict': group_dict,
                    'control_id': str(self.control_id),
                    'compare_plan': compare_plan
                }, sort_keys=True, separators=(',', ':'))
                api.add_sg_splicing_rmats(params, outpath)

    @workfuncdeco
    def export_rmats_count(self):
        api = self.api.api('medical_transcriptome.rmats_count')
        api.add_rmats_count(self.as_analysis.output_dir)

    @workfuncdeco
    def export_gene_fusion(self):
        api = self.api.api('medical_transcriptome.gene_fusion')
        project_sn = self.project_sn
        params = json.dumps({
            'task_id': self.task_id,
            'submit_location': 'genefusion',
            'task_type': 2,
            'min_junction_reads':'1',
            'min_sum_frags': '2',
            'min_novel_junction_support': '3',
            'min_spanning_frags_only': '5',
            'min_FFPM': '0.1'
        }, sort_keys=True, separators=(',', ':'))
        fusion_result = self.gene_fusion.output_dir
        extract_all_gene_pos(self.ref_gtf, os.path.join(self.export_temporary, "gene_pos"))
        chr_length_path = os.path.join(os.path.dirname(self.ref_genome), "star_27",
                                            "ctat_genome_lib_build_dir", "ref_genome.fa.star.idx", "chrNameLength.txt")
        assemble_level_file = os.path.join(os.path.dirname(self.ref_gtf), "assembly_level.txt")
        extract_all_chr_length(assemble_level_file, chr_length_path,
                               os.path.join(self.export_temporary, "chr_length"))
        s3_output = os.path.join(self._sheet.output, '03Gene_structure_analysis', "03GeneFusion",
                                 "GeneFusion_Star_Fusion", "Star_fusion")
        api.add_fusion_main(fusion_result=fusion_result,pos_file=os.path.join(self.export_temporary, "gene_pos"),chr_length = os.path.join(self.export_temporary, "chr_length"),circos= 1,task_id=self.task_id,params=params,project_sn=project_sn,s3_output = s3_output, group=self.option('group_table').path)

    def export_as_profile(self):
        api = self.api.api('medical_transcriptome.asprofile')
        params = json.dumps({
            'task_id': self.task_id,
            'submit_location': 'asprofile',
            'task_type': 2,
            'software': 'Asprofile',
        }, sort_keys=True, separators=(',', ':'))
        result_as = os.path.join(self.asprofile.output , 'AS_result_merge.txt')
        as_statistics = os.path.join(self.asprofile.output_dir, 'AS_statistics_merge.txt')
        s3_as_result = os.path.join(self._sheet.output, '03Gene_structure_analysis', '01AS','AS_result_merge.txt')
        main_id =  api.add_asprofile_result(result_as, s3_as_result,params=params)
        group_dict = self.option('group_table').prop['group_dict']
        sample_list = ",".join([j for group in group_dict for j in group_dict[group]])
        api.add_asprofile_statistics(as_statistics,main_id,sample_list)

    def export_diff_geneset_analysis(self):
        diff_geneset_pipline_result = self.diff_geneset_analysis.output_dir
        diff_id = self.diff_id
        task_id = self.task_id
        analysis_names = ["kegg","go","reactome"]
        if self.organism_name == "Homo_sapiens":
            analysis_names.append("do")
        file_json_path = os.path.join(self.diff_geneset_analysis.file_prepare.output_dir, "prepare_json")
        with open(file_json_path,"r") as j:
            file_dict = json.load(j)
        kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
        # all_annot_path = os.path.join(self.annot_merge.output_dir,"allannot_class","all_annot.xls")
        all_annot_path = os.path.join(self.annot_merge.output_dir, "allannot_class", "all_annot_gene.xls")
        all_annot_df  = pd.read_table(all_annot_path)
        annot_df = all_annot_df[["gene_id", "gene_name", "description"]]
        annot_df.to_csv(os.path.join(self.export_temporary,"gene_detail"),sep="\t",index=False)
        gene_detail = os.path.join(self.export_temporary,"gene_detail")
        api = self.api.api('medical_transcriptome.diff_geneset_work_pipline')
        api.add_diff_genest_pipline_table(diff_geneset_pipline_result,diff_id = diff_id, task_id=task_id , analysis_names = analysis_names,
                                          kegg_level_path=kegg_level_path,inter_path= self.export_temporary,exp_id =self.exp_ids['G'])

    def export_report_img(self):
        report_config = os.path.join(self.chart.work_dir, 'report_config.json')
        api =  self.api.api('medical_transcriptome.report_model')
        s3 = self._sheet.output.split(":")[0]
        report_img_s3 = s3 + ":commonbucket/files/report_img/medical/" + self.task_id
        member_id = self._sheet.member_id.split("m_")[-1]
        api.add_report_image(self.task_id, report_config, report_img_s3, member_id)

    def merge_annotation_exp_matrix(self):
        if self.option('sample_num') == 'multiple':
            group_dict = self.option('group_table').prop['group_dict']
        else:
            group_table = os.path.join(self.work_dir, 'group.txt')
            for line in open(group_table):
                if line.strip() and line[0] != '#':
                    sample = line.strip().split('\t')[1]
            group_dict = {sample: [sample]}
        exp_output = self.quant.output_dir
        if self.option('is_assemble'):
            annot = os.path.join(self.output_dir, 'annot_merge/allannot_class/all_annot_tran.xls')
        else:
            annot = os.path.join(self.output_dir, 'annot_merge/refannot_class/all_annot_tran.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        gene_annot_pd = all_annot[all_annot['is_gene'] == 'yes'].drop(
            columns=['transcript_id', 'is_gene', 'gene_name', 'description', 'length'])
        order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"]
        gene_annot_pd = pd.DataFrame(gene_annot_pd, columns=order)
        gene_info_pd = all_annot[all_annot['is_gene'] == 'yes'][['gene_name', 'description', 'length']]
        trans_annot_pd = all_annot.reset_index().drop(
            columns=['gene_id', 'is_gene', 'gene_name', 'description', 'length']).set_index('transcript_id')
        trans_annot_pd = pd.DataFrame(trans_annot_pd, columns=order)
        trans_info_pd = all_annot[['transcript_id', 'gene_name', 'description', 'length']].reset_index().set_index(
            'transcript_id')
        # gene
        ## gene tpm
        gene_tpm_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
        gene_tpm_pd = pd.read_table(gene_tpm_matrix, header=0, index_col=0)
        gene_tpm_dicts = gene_tpm_pd.to_dict('index')
        gene_group_tpm = OrderedDict()
        for seq_id in sorted(gene_tpm_dicts):
            tmp_exp_dict = gene_tpm_dicts[seq_id]
            for group in group_dict:
                if seq_id not in gene_group_tpm:
                    gene_group_tpm[seq_id] = dict()
                gene_group_tpm[seq_id].update(
                    {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
        gene_group_pd = (pd.DataFrame(data=gene_group_tpm, columns=gene_group_tpm.keys())).T
        gene_tpm_result = pd.concat([gene_info_pd, gene_tpm_pd, gene_group_pd, gene_annot_pd], axis=1)
        gene_tpm_out = os.path.join(exp_output, 'gene.tpm.matrix.annot.xls')
        header = ['gene_id']
        header.extend(gene_tpm_result.columns.tolist())
        with open(gene_tpm_out, "w") as w:
            w.write("\t".join(header) + "\n")
        gene_tpm_result.to_csv(gene_tpm_out, header=False, index=True, sep='\t', mode='a')
        ## gene count
        gene_count_matrix = os.path.join(exp_output, 'gene.count.matrix')
        gene_count_pd = pd.read_table(gene_count_matrix, header=0, index_col=0)
        gene_count_dicts = gene_count_pd.to_dict('index')
        gene_group_count = OrderedDict()
        for seq_id in sorted(gene_count_dicts):
            tmp_exp_dict = gene_count_dicts[seq_id]
            for group in group_dict:
                if seq_id not in gene_group_count:
                    gene_group_count[seq_id] = dict()
                gene_group_count[seq_id].update(
                    {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
        gene_group_pd = (pd.DataFrame(data=gene_group_count, columns=gene_group_count.keys())).T
        gene_count_result = pd.concat([gene_info_pd, gene_count_pd, gene_group_pd, gene_annot_pd], axis=1)
        gene_count_out = os.path.join(exp_output, 'gene.count.matrix.annot.xls')
        header = ['gene_id']
        header.extend(gene_count_result.columns.tolist())
        with open(gene_count_out, "w") as w:
            w.write("\t".join(header) + "\n")
        gene_count_result.to_csv(gene_count_out, header=False, index=True, sep='\t', mode='a')
        ## gene fpkm
        if self.option('express_method') == 'RSEM':
            gene_fpkm_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
            gene_fpkm_pd = pd.read_table(gene_fpkm_matrix, header=0, index_col=0)
            gene_fpkm_dicts = gene_fpkm_pd.to_dict('index')
            gene_group_fpkm = OrderedDict()
            for seq_id in sorted(gene_fpkm_dicts):
                tmp_exp_dict = gene_fpkm_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in gene_group_fpkm:
                        gene_group_fpkm[seq_id] = dict()
                    gene_group_fpkm[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            gene_group_pd = (pd.DataFrame(data=gene_group_fpkm, columns=gene_group_fpkm.keys())).T
            gene_fpkm_result = pd.concat([gene_info_pd, gene_fpkm_pd, gene_group_pd, gene_annot_pd], axis=1)
            gene_fpkm_out = os.path.join(exp_output, 'gene.fpkm.matrix.annot.xls')
            header = ['gene_id']
            header.extend(gene_fpkm_result.columns.tolist())
            with open(gene_fpkm_out, "w") as w:
                w.write("\t".join(header) + "\n")
            gene_fpkm_result.to_csv(gene_fpkm_out, header=False, index=True, sep='\t', mode='a')
        if self.option('level').lower() == 'transcript':
            ## transcript tpm
            transcript_tpm_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            transcript_tpm_pd = pd.read_table(transcript_tpm_matrix, header=0, index_col=0)
            transcript_tpm_dicts = transcript_tpm_pd.to_dict('index')
            transcript_group_tpm = OrderedDict()
            for seq_id in sorted(transcript_tpm_dicts):
                tmp_exp_dict = transcript_tpm_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in transcript_group_tpm:
                        transcript_group_tpm[seq_id] = dict()
                    transcript_group_tpm[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            transcript_group_pd = (pd.DataFrame(data=transcript_group_tpm, columns=transcript_group_tpm.keys())).T
            transcript_tpm_result = pd.concat([trans_info_pd, transcript_tpm_pd, transcript_group_pd, trans_annot_pd],
                                              axis=1)
            transcript_tpm_out = os.path.join(exp_output, 'transcript.tpm.matrix.annot.xls')
            header = ['transcript_id']
            header.extend(transcript_tpm_result.columns.tolist())
            with open(transcript_tpm_out, "w") as w:
                w.write("\t".join(header) + "\n")
            transcript_tpm_result.to_csv(transcript_tpm_out, header=False, index=True, sep='\t', mode='a')
            ## transcript fpkm
            if self.option('express_method') == 'RSEM':
                transcript_fpkm_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
                transcript_fpkm_pd = pd.read_table(transcript_fpkm_matrix, header=0, index_col=0)
                transcript_fpkm_dicts = transcript_fpkm_pd.to_dict('index')
                transcript_group_fpkm = OrderedDict()
                for seq_id in sorted(transcript_fpkm_dicts):
                    tmp_exp_dict = transcript_fpkm_dicts[seq_id]
                    for group in group_dict:
                        if seq_id not in transcript_group_fpkm:
                            transcript_group_fpkm[seq_id] = dict()
                        transcript_group_fpkm[seq_id].update({group: round(
                            sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
                transcript_group_pd = (pd.DataFrame(data=transcript_group_fpkm, columns=transcript_group_fpkm.keys())).T
                transcript_fpkm_result = pd.concat(
                    [trans_info_pd, transcript_fpkm_pd, transcript_group_pd, trans_annot_pd], axis=1)
                transcript_fpkm_out = os.path.join(exp_output, 'transcript.fpkm.matrix.annot.xls')
                header = ['transcript_id']
                header.extend(transcript_fpkm_result.columns.tolist())
                with open(transcript_fpkm_out, "w") as w:
                    w.write("\t".join(header) + "\n")
                transcript_fpkm_result.to_csv(transcript_fpkm_out, header=False, index=True, sep='\t', mode='a')
            ## transcript count
            transcript_count_matrix = os.path.join(exp_output, 'transcript.count.matrix')
            transcript_count_pd = pd.read_table(transcript_count_matrix, header=0, index_col=0)
            transcript_count_dicts = transcript_count_pd.to_dict('index')
            transcript_group_count = OrderedDict()
            for seq_id in sorted(transcript_count_dicts):
                tmp_exp_dict = transcript_count_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in transcript_group_count:
                        transcript_group_count[seq_id] = dict()
                    transcript_group_count[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            transcript_group_pd = (pd.DataFrame(data=transcript_group_count, columns=transcript_group_count.keys())).T
            transcript_count_result = pd.concat(
                [trans_info_pd, transcript_count_pd, transcript_group_pd, trans_annot_pd], axis=1)
            transcript_count_out = os.path.join(exp_output, 'transcript.count.matrix.annot.xls')
            header = ['transcript_id']
            header.extend(transcript_count_result.columns.tolist())
            with open(transcript_count_out, "w") as w:
                w.write("\t".join(header) + "\n")
            transcript_count_result.to_csv(transcript_count_out, header=False, index=True, sep='\t', mode='a')

    def merge_annotation_diffexp_matrix(self):
        if self.option("is_assemble"):
            annot = os.path.join(self.output_dir, 'annot_merge/allannot_class/all_annot_tran.xls')
        else:
            annot = os.path.join(self.output_dir, 'annot_merge/refannot_class/all_annot_tran.xls')
        all_annot = pd.read_table(annot, header=0, index_col=0)
        gene_annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(
            columns=['transcript_id', 'is_gene', 'gene_name', 'description', 'length'])
        order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"]
        gene_annot_pd = pd.DataFrame(gene_annot_pd, columns=order)
        gene_info_pd = all_annot[all_annot['is_gene'] == 'yes'][['gene_name', 'description', 'length']]
        diff_output = self.diffexp.output_dir
        duplicate_files = glob.glob(diff_output + '/' + '*.annot.xls') + glob.glob(
            diff_output + '/' + '*_vs_*.normalize.xls') + glob.glob(diff_output + '/' + '*_vs_*.sizeFactor.xls')
        for file in duplicate_files:
            os.remove(os.path.join(diff_output, file))
        target_files = glob.glob(diff_output + "/*.xls")
        for each in target_files:
            gene_pd = pd.read_table(each, header=0, index_col=0)
            gene_result = pd.concat([gene_info_pd, gene_pd, gene_annot_pd], join='inner', axis=1)
            gene_out = each.split('.xls')[0] + '.annot.xls'
            header = ['gene_id']
            header.extend(gene_result.columns.tolist())
            with open(gene_out, "w") as w:
                w.write("\t".join(header) + "\n")
            gene_result.to_csv(gene_out, header=False, index=True, sep='\t', mode='a')

    def get_fq_list(self):
        if self.option('datatype') == 'rawdata':
            if self.option('qc_soft') == 'fastp':
                self.fq_list = self.fastp_rna.option('fq_list').path
            elif self.option('qc_soft') == 'seqprep':
                self.fq_list = self.hiseq_qc.option('fq_list').path
        else:
            dct = dict()
            fq_path = self.option('qc_dir').path
            for line in open(os.path.join(fq_path, 'list.txt')):
                items = line.strip().split('\t')
                if items[1] in dct:
                    dct[items[1]][items[2]] = os.path.join(fq_path, items[0])
                    flag = 1
                else:
                    dct[items[1]] = {items[2]: os.path.join(fq_path, items[0])}
                    flag = 0
            else:
                self.fq_list = os.path.join(fq_path, 'fq_list.txt')
                open(self.fq_list, 'w').writelines(['{}\t{l}\t{r}\n'.format(k, **v) for k, v in dct.items()] if flag \
                                                       else ['{}\t{s}\n'.format(k, **v) for k, v in dct.items()])
        return self.fq_list



class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """


    # def test(self):
    #     from mbio.workflows.medical_transcriptome.medical_transcriptome import MedicalTranscriptomeWorkflow
    #     from biocluster.wsheet import Sheet
    #     import random
    #     id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
    #     data = {
    #         "id": "medical_transcriptome_workflow_" + id,
    #         "type": "workflow",
    #         "name": "medical_transcriptome.medical_transcriptome",
    #         "options": dict(
    #             ref_genome="Homo_sapiens",
    #             genome_id="GM0259",
    #             fastq_dir="/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/remote_input/fastq_dir/raw_data",
    #             group_table='/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/remote_input/group_table/example_group_1528169151.txt',
    #             control_file='/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/remote_input/control_file/example_control_1528169151.txt',
    #         )
    #     }
    #     wsheet = Sheet(data=data)
    #     wf = MedicalTranscriptomeWorkflow(wsheet)
    #     wf.sheet.id = 'test_medical'
    #     wf.sheet.project_sn = 'test_medical'
    #     wf.IMPORT_REPORT_DATA = False
    #     wf.IMPORT_REPORT_AFTER_DATA = False
    #     wf.run()

    def test(self):
        worker = worker_client()
        id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        data = {
            "id": "medical_transcriptome_workflow_" + id,
            "type": "workflow",
            "name": "medical_transcriptome.medical_transcriptome",
            "options": dict(
                ref_genome="Homo_sapiens",
                genome_id="GM0259",
                fastq_dir="/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/remote_input/fastq_dir/raw_data",
                group_table='/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/remote_input/group_table/example_group_1528169151.txt',
                control_file='/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38314/remote_input/control_file/example_control_1528169151.txt',
            )
        }

        info = worker.add_task(data)
        print(info)


if __name__ == '__main__':
    unittest.main()
