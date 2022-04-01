# -*- coding:utf-8 -*-
# __author__ = 'shicaiping,liubinxu,qinjincheng'

import glob
import json
import os
import shutil
import tarfile
import unittest
import pandas as pd
from biocluster.config import Config
from biocluster.workflow import Workflow
from mbio.packages.dna_evolution.send_email import SendEmail
from mbio.packages.ref_rna_v2.functions import tryforgood


class CircrnaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CircrnaWorkflow, self).__init__(wsheet_object)
        options = [
            ## 分析对象
            # 分析文库选择
            # {'name': 'lib_select', 'type': 'string', 'default': 'longRNA,smallRNA'},
            #
            # ## 基础参数设置
            # # 分析对象选择
            # {'name': 'rna_select', 'type': 'string', 'default': 'mRNA,lncRNA,smallRNA,circRNA'},
            # # 数据类型 ['rawdata', 'cleandata']
            {'name': 'datatype', 'type': 'string', 'default': 'rawdata'},
            # # 测序类型 ['PE', 'SE']
            {'name': 'fq_type', 'type': 'string', 'default': 'PE'},
            # # 链特异性 [True, False]
            {'name': 'strand_specific', 'type': 'bool', 'default': True},
            # # 链特异性方向 {'PE': ['RF', 'FR'], 'SE': ['R', 'F']}
            {'name': 'strand_dir', 'type': 'string', 'default': 'RF'},
            # # 质控软件 ['fastp', 'seqprep']
            {'name': 'qc_soft', 'type': 'string', 'default': 'fastp'},
            # # 测序质量 ['phred_33', 'phred_64']
            {'name': 'quality_score_system', 'type': 'string', 'default': 'phred_33'},
            # 长度需求
            {'name': 'length_required', 'type': 'string', 'default': '30'},
            # 接头序列1
            {'name': 'adapter_a', 'type': 'string', 'default': 'AGATCGGAAGAGCACACGTC'},
            # 接头序列2
            {'name': 'adapter_b', 'type': 'string', 'default': 'AGATCGGAAGAGCGTCGTGT'},
            # # 生物学重复 [True, False]
            {'name': 'is_duplicate', 'type': 'bool', 'default': True},
            # # 原始序列文件
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            # # 质控序列文件
            {'name': 'qc_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            # # 分组方案
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            # # 对照组文件
            {'name': 'control_file', 'type': 'infile', 'format': 'sample.control_table'},
            # 配对信息表
            {'name': 'pair_table', 'type': 'infile', 'format': 'sample.group_table'},
            # # 物种类别 ['Animal', 'Plant', 'Protist', 'Fungi']
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            # # 具体物种 sg_genome_db.organism_name
            {'name': 'organism_name', 'type': 'string', 'default': None},
            # # 基因组编号 sg_genome_db.genome_id
            {'name': 'genome_id', 'type': 'string', 'default': None},
            # # 已知lncRNA序列文件 比对参数（选填）
            # {'name': 'knownlnc_fasta', 'type': 'infile', 'format': 'sequence.fasta'},
            # {'name': 'knownlnc_evalue', 'type': 'float', 'default': 1e-5},
            # {'name': 'knownlnc_qcov', 'type': 'float', 'default': 80},
            # {'name': 'knownlnc_scov', 'type': 'float', 'default': 80},
            # # 终止分析阈值条件
            {'name': 'mapping_stop', 'type': 'bool', 'default': True},
            # # 大于 %的样本
            {'name': 'mapping_sample_percent', 'type': 'float', 'default': 50.0},
            # # Mapping Ratio 小于
            {'name': 'mapping_ratio', 'type': 'float', 'default': 60.0},
            # # 终止后续分析
            {'name': 'rrna_stop', 'type': 'bool', 'default': True},
            # # 大于 %的样本
            {'name': 'rrna_sample_percent', 'type': 'float', 'default': 50.0},
            # # rRNA Ratio 小于 %
            {'name': 'rrna_ratio', 'type': 'float', 'default': 15.0},
            #
            # ## 高级参数设置
            # # 比对软件 ['Hisat', 'Tophat', 'STAR']
            {'name': 'align_method', 'type': 'string', 'default': 'Hisat'},
            {'name': 'assemble_method', 'type': 'string', 'default': 'stringtie'},
            # # lncRNA初步筛选length(nt)
            # {'name': 'transcript_len', 'type': 'int', 'default': 200},
            # # lncRNA初步筛选exon（个数）
            # {'name': 'exon_num', 'type': 'int', 'default': 2},
            # # lncRNA初步筛选ORF(nt)
            # {'name': 'orf_len', 'type': 'int', 'default': 300},
            # # 编码能力预测CPC
            # {'name': 'cpc', 'type': 'string', 'default': 'True'},
            # # CPC score
            # {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            # # 编码能力预测CNCI
            # {'name': 'cnci', 'type': 'string', 'default': 'True'},
            # # CNCI score
            # {'name': 'cnci_score', 'type': 'float', 'default': 0.0},
            # # 编码能力预测CPAT
            # {'name': 'cpat', 'type': 'string', 'default': 'True'},
            # # CPAT score
            # {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            # # 编码能力预测PfamScan
            # {'name': 'pfamscan', 'type': 'string', 'default': 'True'},
            # # 候选lncRNA确定标准
            # {'name': 'identify_num', 'type': 'int', 'default': 2},

            # # NR库分类 ['Animal, Plant', 'Protist', 'Fungi', 'All']
            # {'name': 'nr_database', 'type': 'string', 'default': None},
            # # KEGG库分类 ['Animal, Plant', 'Protist', 'Fungi', 'All']
            # {'name': 'kegg_database', 'type': 'string', 'default': None},
            # # NR(GO)-Evalue
            # {'name': 'nrgo_evalue', 'type': 'float', 'default': 1e-5},
            # # SwissProt-Evalue
            # {'name': 'swissprot_evalue', 'type': 'float', 'default': 1e-5},
            # # KEGG-Evalue
            # {'name': 'kegg_evalue', 'type': 'float', 'default': 1e-5},
            # # COG-Evalue
            # {'name': 'cog_evalue', 'type': 'float', 'default': 1e-5},
            # # Pfam-Evalue
            # {'name': 'pfam_evalue', 'type': 'float', 'default': 1e-5},

            # circRNA鉴定软件 ['ciri2+find_circ', 'circ2', 'find_circ', 'circ_finder', 'circexplorer2']
            {'name': 'circ_method', 'type': 'string', 'default': None},
            # min BSJ reads
            {'name': 'junction_reads', 'type': 'int', 'default': 2},
            # max circRNA length
            {'name': 'circrna_length', 'type': 'int', 'default': 100000},

            # 表达定量软件 ['RSEM', 'Salmon', 'Kallisto']
            {'name': 'exp_method', 'type': 'string', 'default': 'Salmon'},
            # 表达定量指标 ['tpm', 'fpkm']
            {'name': 'exp_way', 'type': 'string', 'default': 'tpm'},
            # Filter TPM
            {'name': 'exp_threshold', 'type': 'float', 'default': 0.0},
            # 差异分析软件 ['DESeq2', 'edgeR', 'DEGseq']
            {'name': 'diff_method', 'type': 'string', 'default': 'DESeq2'},
            # 表达量筛选
            {'name': 'diff_filter', 'type': 'string', 'default': None},
            {'name': 'diff_threshold', 'type': 'float', 'default': 0.0},
            # FC
            {'name': 'fc', 'type': 'float', 'default': 2.0},
            # 显著性水平
            {'name': 'pvalue_padjust', 'type': 'string', 'default': 'padjust'},
            {'name': 'diff_fdr_ci', 'type': 'float', 'default': 0.05},
            # 多重验证校正方法 ['BH', 'Bonferroni', 'Holm', 'BY']
            {'name': 'padjust_way', 'type': 'string', 'default': 'BH'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        self.modules = dict()
        self.tools = dict()

        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
        collection = database['sg_genome_db']
        self.genome_doc = collection.find_one({'genome_id': self.option('genome_id')})
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish')

        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        data = os.path.join(self.work_dir, 'data.json')
        if os.path.exists(data):
            with open(data, 'r') as load_f:
                load_dict = json.load(load_f)
                if 'rerun' in load_dict and load_dict['rerun']:
                    self.logger.info("该项目重运行中，先删除mongo库中已有数据")
                    self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        self.program = os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin/python')
        cmd = '{} {}'.format(self.program, self.script)
        cmd += ' {} {}'.format(self.task_id, 'whole_transcriptome')
        code = os.system(cmd)
        if code == 0:
            self.logger.info("命令{}执行成功！".format(cmd))
        else:
            raise Exception("命令{}执行失败！".format(cmd))

    def check_options(self):
        # if self.option('nr_database') == 'All':
        #     self.option('nr_database', 'nr')
        # elif self.option('nr_database') == 'Animal':
        #     self.option('nr_database', 'metazoa')
        # elif self.option('nr_database') == 'Plant':
        #     self.option('nr_database', 'viridiplantae')
        # elif self.option('nr_database') == 'Protist':
        #     self.option('nr_database', 'protist')
        # elif self.option('nr_database') == 'Fungi':
        #     self.option('nr_database', 'fungi')
        # if self.option('kegg_database') == 'All':
        #     self.option('kegg_database', 'All')
        # elif self.option('kegg_database') == 'Animal':
        #     self.option('kegg_database', 'Animals')
        # elif self.option('kegg_database') == 'Plant':
        #     self.option('kegg_database', 'Plants')
        # elif self.option('kegg_database') == 'Fungi':
        #     self.option('kegg_database', 'Fungi')
        # elif self.option('kegg_database') == 'Protist':
        #     self.option('kegg_database', 'Protists')
        for k, v in self.sheet.options().items():
            self.logger.debug('{} -> {}'.format(k, v))
        else:
            return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        self.add_steps()
        super(CircrnaWorkflow, self).run()

    def add_steps(self):
        self.step.add_steps('file_check')
        self.step.add_steps('hiseq_reads_stat_raw')
        self.step.add_steps('fastp_rna')
        self.step.add_steps('hiseq_reads_stat_use')
        self.step.add_steps('rnaseq_mapping')
        self.step.add_steps('map_assessment')
        self.step.add_steps('circ_brush')
        self.step.add_steps('exp_makecirc')
        self.step.add_steps('diff_exp_c')
        self.load_libraries()

    def load_libraries(self):
        self.tools['file_check'] = self.add_tool('whole_transcriptome.longrna.file_check')
        self.modules['hiseq_reads_stat_raw'] = self.add_module('whole_transcriptome.longrna.hiseq_reads_stat')
        if self.option('qc_soft') == 'fastp':
            self.modules['fastp_rna'] = self.add_module('datasplit.fastp_rna')
        else:
            self.modules['hiseq_qc'] = self.add_module('ref_rna_v2.hiseq_qc')
        self.modules['hiseq_reads_stat_use'] = self.add_module('whole_transcriptome.longrna.hiseq_reads_stat')
        self.modules['rnaseq_mapping'] = self.add_module('whole_transcriptome.longrna.rnaseq_mapping')
        self.modules['map_assessment'] = self.add_module('whole_transcriptome.longrna.map_assessment')
        self.modules['circ_brush'] = self.add_module('whole_transcriptome.circ_brush')
        self.tools['exp_makecirc'] = self.add_tool('whole_transcriptome.formation.exp_makecirc')
        self.modules['diff_exp_c'] = self.add_module('whole_transcriptome.diff_exp')
        self.run_file_check()

    def run_file_check(self):
        if self.option('datatype') == 'rawdata':
            fastq_dir = self.option('fastq_dir').path
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir').path
        in_gtf = os.path.join(self.db_path, self.genome_doc['gtf'])
        sample_num = 'multiple'
        opts = {
            'fq_type': self.option('fq_type'),
            'fastq_dir': fastq_dir,
            'in_gtf': in_gtf,
            'sample_num': sample_num,
            'group_table': self.option('group_table'),
            'control_file': self.option('control_file')
        }
        self.tools['file_check'].set_options(opts)
        self.tools['file_check'].on('start', self.set_step, {'start': self.step.file_check})
        self.tools['file_check'].on('end', self.set_step, {'end': self.step.file_check})
        self.tools['file_check'].on('end', self.set_output, 'file_check')
        if self.option('datatype') == 'rawdata':
            self.on_rely([self.modules['hiseq_reads_stat_raw'], self.modules['hiseq_reads_stat_use']], self.check_rrna)
            self.tools['file_check'].on('end', self.run_hiseq_reads_stat_raw)
            if self.option('qc_soft') == 'fastp':
                self.tools['file_check'].on('end', self.run_fastp_rna)
            elif self.option('qc_soft') == 'seqprep':
                self.tools['file_check'].on('end', self.run_hiseq_qc)
        elif self.option('datatype') == 'cleandata':
            self.tools['file_check'].on('end', self.run_hiseq_reads_stat_use)
        self.tools['file_check'].run()

    def run_hiseq_reads_stat_raw(self):
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'quality': quality
        }
        self.modules['hiseq_reads_stat_raw'].set_options(opts)
        self.modules['hiseq_reads_stat_raw'].on('start', self.set_step, {'start': self.step.hiseq_reads_stat_raw})
        self.modules['hiseq_reads_stat_raw'].on('end', self.set_step, {'end': self.step.hiseq_reads_stat_raw})
        self.modules['hiseq_reads_stat_raw'].on('end', self.set_output, 'hiseq_reads_stat_raw')
        self.modules['hiseq_reads_stat_raw'].run()

    def run_hiseq_qc(self):
        if self.option('fq_type') == "PE":
            self.modules['hiseq_qc'].set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_a': self.option('adapter_a'),
                'adapter_b': self.option('adapter_b')
            })
        else:
            self.modules['hiseq_qc'].set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_a': self.option('adapter_a'),
            })
        self.modules['hiseq_qc'].on('start', self.set_step, {'start': self.step.hiseq_qc})
        self.modules['hiseq_qc'].on('end', self.set_step, {'end': self.step.hiseq_qc})
        self.modules['hiseq_qc'].on('end', self.set_output, 'hiseq_qc')
        self.modules['hiseq_qc'].on('end', self.run_hiseq_reads_stat_use)
        self.modules['hiseq_qc'].run()

    def run_fastp_rna(self):
        fq_dir = self.option('fastq_dir').path
        sample_path = os.path.join(fq_dir, 'abs.list.txt')
        open(sample_path, 'w').writelines(
            '{}/{}'.format(fq_dir, line) for line in open(os.path.join(fq_dir, 'list.txt'))
        )
        if self.option('fq_type') == "PE":
            self.modules['fastp_rna'].set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'length_required': self.option('length_required'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_sequence': self.option('adapter_a'),
                'adapter_sequence_r2': self.option('adapter_b')
            })
        else:
            self.modules['fastp_rna'].set_options({
                'sample_path': sample_path,
                'fq_type': self.option('fq_type'),
                'length_required': self.option('length_required'),
                'quality_score_system': self.option('quality_score_system'),
                'adapter_sequence_s': self.option('adapter_a'),
            })
        self.modules['fastp_rna'].on('start', self.set_step, {'start': self.step.fastp_rna})
        self.modules['fastp_rna'].on('end', self.set_step, {'end': self.step.fastp_rna})
        self.modules['fastp_rna'].on('end', self.set_output, 'fastp_rna')
        self.modules['fastp_rna'].on('end', self.run_hiseq_reads_stat_use)
        self.modules['fastp_rna'].run()

    def run_hiseq_reads_stat_use(self):
        if self.option('datatype') == 'rawdata':
            fastq_dir = self.modules['fastp_rna'].option('sickle_dir').path
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir').path
        if self.option('quality_score_system').endswith('33'):
            quality = 33
        elif self.option('quality_score_system').endswith('64'):
            quality = 64
        self.modules['hiseq_reads_stat_use'].set_options({
            'fastq_dir': fastq_dir,
            'fq_type': self.option('fq_type'),
            'quality': quality,
            'dup': True,
            'rfam': True,
            'rrna_ratio': self.option('rrna_ratio'),
        })
        self.modules['hiseq_reads_stat_use'].on('start', self.set_step, {'start': self.step.hiseq_reads_stat_use})
        self.modules['hiseq_reads_stat_use'].on('end', self.set_step, {'end': self.step.hiseq_reads_stat_use})
        self.modules['hiseq_reads_stat_use'].on('end', self.set_output, 'hiseq_reads_stat_use')
        if self.option('datatype') == 'cleandata':
            self.modules['hiseq_reads_stat_use'].on('end', self.check_rrna)
        self.modules['hiseq_reads_stat_use'].run()

    def check_rrna(self):
        if self.option('rrna_stop'):
            if self.modules['hiseq_reads_stat_use'].option('rrna_sample_percent') > self.option('rrna_sample_percent'):
                self.stop('rrna')
            else:
                self.run_rnaseq_mapping()
        else:
            self.run_rnaseq_mapping()

    def run_rnaseq_mapping(self):
        ref_genome = self.option('organism_name')
        genome_version = self.genome_doc['assembly']
        genome_annot_version = self.genome_doc['annot_version']
        mapping_method = self.option('align_method').lower()
        seq_method = self.option('fq_type')
        if self.option('datatype') == 'rawdata':
            fastq_dir = self.modules['fastp_rna'].option('sickle_dir').path
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir').path
        assemble_method = self.option('assemble_method').lower()
        opts = {
            'ref_genome': ref_genome,
            'genome_version': genome_version,
            'genome_annot_version': genome_annot_version,
            'mapping_method': mapping_method,
            'seq_method': seq_method,
            'fastq_dir': fastq_dir,
            'assemble_method': assemble_method,
            'strand_specific': self.option('strand_specific'),
            'mapping_ratio': self.option('mapping_ratio')
        }
        self.modules['rnaseq_mapping'].set_options(opts)
        self.modules['rnaseq_mapping'].on('start', self.set_step, {'start': self.step.rnaseq_mapping})
        self.modules['rnaseq_mapping'].on('end', self.set_step, {'end': self.step.rnaseq_mapping})
        self.modules['rnaseq_mapping'].on('end', self.set_output, 'rnaseq_mapping')
        self.modules['rnaseq_mapping'].on('end', self.check_mapping)
        self.modules['rnaseq_mapping'].run()

    def check_mapping(self):
        if self.option('mapping_stop'):
            if self.modules['rnaseq_mapping'].option('mapping_sample_percent') > self.option('mapping_sample_percent'):
                self.stop('mapping')
            else:
                self.on_rely([self.modules['map_assessment'], self.modules['diff_exp_c']], self.set_db)
                self.run_map_assessment()
                self.run_circ_brush()
        else:
            self.on_rely([self.modules['map_assessment'], self.modules['diff_exp_c']], self.set_db)
            self.run_map_assessment()
            self.run_circ_brush()

    def run_map_assessment(self):
        opts = {
            'bam': self.modules['rnaseq_mapping'].option('bam_output'),
            'bed': self.tools['file_check'].option('bed')
        }
        self.modules['map_assessment'].set_options(opts)
        self.modules['map_assessment'].on('start', self.set_step, {'start': self.step.map_assessment})
        self.modules['map_assessment'].on('end', self.set_step, {'end': self.step.map_assessment})
        self.modules['map_assessment'].on('end', self.set_output, 'map_assessment')
        self.modules['map_assessment'].run()

    def run_circ_brush(self):
        if 'ciri2' in self.option('circ_method') and 'find_circ' in self.option('circ_method'):
            circ_method = 'ciri2,find_circ'
        else:
            circ_method = self.option('circ_method')
        genome = os.path.join(self.db_path, self.genome_doc['dna_fa'])
        if 'lnc_dir' in self.genome_doc:
            annotate = os.path.join(self.db_path, self.genome_doc['lnc_dir'], 'mrna.gtf')
        else:
            annotate = os.path.join(self.db_path, self.genome_doc['gtf'])

        if self.option('datatype') == 'rawdata':
            fastq_dir = self.modules['fastp_rna'].option('sickle_dir').path
        elif self.option('datatype') == 'cleandata':
            fastq_dir = self.option('qc_dir').path
        opts = {
            'circ_method': circ_method,
            'genome': genome,
            'annotate': annotate,
            'fastq_dir': fastq_dir,
            'organism_name': self.option('organism_name'),
            'junction_reads': self.option('junction_reads'),
            'circrna_length': self.option('circrna_length')
        }
        self.modules['circ_brush'].set_options(opts)
        self.modules['circ_brush'].on('start', self.set_step, {'start': self.step.circ_brush})
        self.modules['circ_brush'].on('end', self.set_step, {'end': self.step.circ_brush})
        self.modules['circ_brush'].on('end', self.set_output, 'circ_brush')
        self.modules['circ_brush'].on('end', self.run_exp_makecirc)
        self.modules['circ_brush'].run()

    def run_exp_makecirc(self):

        c_rpm = os.path.join(self.modules['circ_brush'].output_dir, 'RPM.txt')
        c_detail = os.path.join(self.modules['circ_brush'].output_dir, 'detail.txt')
        c_count = os.path.join(self.modules['circ_brush'].output_dir, 'count.txt')
        opts = {

            'c_rpm': c_rpm,
            'c_detail': c_detail,
            'c_count': c_count
        }
        self.tools['exp_makecirc'].set_options(opts)
        self.tools['exp_makecirc'].on('start', self.set_step, {'start': self.step.exp_makecirc})
        self.tools['exp_makecirc'].on('end', self.set_step, {'end': self.step.exp_makecirc})
        self.tools['exp_makecirc'].on('end', self.set_output, 'exp_makecirc')
        self.tools['exp_makecirc'].on('end', self.run_diff_exp_c)
        self.tools['exp_makecirc'].run()

    def run_diff_exp_c(self):
        program = self.option('diff_method')
        count_matrix = os.path.join(self.tools['exp_makecirc'].output_dir, 'count/C.reads.txt')
        group_table = self.option('group_table').path
        control_table = self.option('control_file').path
        exp_matrix = os.path.join(self.tools['exp_makecirc'].output_dir, 'circrna/T.rpm.txt')
        kind_table = os.path.join(self.modules['diff_exp_c'].work_dir, 'kind.txt')
        kind_df = pd.read_table(exp_matrix, usecols=['transcript_id', 'kind'])
        kind_df.to_csv(kind_table, sep='\t', index=False)
        threshold = self.option('diff_threshold')
        method = self.option('padjust_way').lower()
        stat_type = self.option('pvalue_padjust')
        stat_cutoff = self.option('diff_fdr_ci')
        fc = self.option('fc')
        opts = {
            'program': program,
            'count_matrix': count_matrix,
            'group_table': group_table,
            'control_table': control_table,
            'exp_matrix': exp_matrix,
            'kind_table': kind_table,
            'filter': str(self.option('diff_filter')).lower(),
            'threshold': threshold,
            'fc': fc
        }
        if program.lower() in ["degseq", "edger", "deseq2", 'limma']:
            opts.update({
                'method': method,
                'stat_type': stat_type,
                'stat_cutoff': stat_cutoff,
            })
            if self.option('pair_table').is_set:
                opts.update({'is_batch': True, 'has_batch': True, 'batch_matrix': self.option('pair_table').path})
        if program.lower() == 'noiseq':
            opts.update({
                stat_type: float(stat_cutoff),
            })
            if self.option('pair_table').is_set:
                opts.update({'is_batch': True, 'has_batch': True, 'batch_matrix': self.option('pair_table').path})
        self.modules['diff_exp_c'].set_options(opts)
        self.modules['diff_exp_c'].on('start', self.set_step, {'start': self.step.diff_exp_c})
        self.modules['diff_exp_c'].on('end', self.set_step, {'end': self.step.diff_exp_c})
        self.modules['diff_exp_c'].on('end', self.set_output, 'diff_exp_c')
        self.modules['diff_exp_c'].run()

    def set_output(self, event):
        obj = event['bind_object']
        self.move2outputdir(obj.output_dir, event['data'])

    def move2outputdir(self, olddir, newname):
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
        else:
            self.logger.debug('succeed in linking {} to {}'.format(olddir, newdir))

    def move_file(self, src, dst):
        if os.path.isfile(src):
            os.link(src, dst)
        else:
            os.mkdir(dst)
            for file in os.listdir(src):
                old_path = os.path.join(src, file)
                new_path = os.path.join(dst, file)
                self.move_file(old_path, new_path)

    def stop(self, reason):
        assert reason in ('rrna', 'mapping')
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.export_task_info()
        self.export_genome_info()
        if self.option('datatype') == 'rawdata':
            self.export_qc()
        elif self.option('datatype') == 'cleandata':
            self.export_qc_after()
        msg = str()
        if reason == 'rrna':
            msg = 'Workflow will stop because the rRNA ratio is up to par.'
            self.logger.warn(msg)
        elif reason == 'mapping':
            self.export_mapping(False)
            msg = 'Workflow will stop because the mapping ratio is not up to par.'
            self.logger.warn(msg)
        receiver = ['caiping.shi@majorbio.com', 'rna_bioinfor@majorbio.com']
        a = SendEmail("897236887@qq.com", "smtp.qq.com", "fhwuvcclstjqbfga", "897236887@qq.com", ",".join(receiver),
                      "WARNING - project_sn ({}), task_id ({})".format(self._sheet.project_sn, self._sheet.id), 465)
        a.send_msg(msg)
        a.send_email()
        super(CircrnaWorkflow, self).end()

    def set_db(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.export_task_info()
        self.export_genome_info()
        if self.option('datatype') == 'rawdata':
            self.export_qc()
        elif self.option('datatype') == 'cleandata':
            self.export_qc_after()
        self.export_mapping()
        self.export_expression_stat_circ()
        self.export_diff_exp_stat_circ()
        self.export_email()
        self.set_upload()

    def set_upload(self):
        shutil.copy(os.path.join(self.work_dir, 'data.json'), os.path.join(self.output_dir, 'data.json'))
        shutil.copy(self.option('group_table').path, os.path.join(self.output_dir, 'group.txt'))
        shutil.copy(self.option('control_file').path, os.path.join(self.output_dir, 'control.txt'))
        for fname in glob.glob(os.path.join(self.output_dir, 'fastp_rna/fastq/*.fastq')):
            os.remove(fname)
        self.pack_and_compress()
        self.add_upload_dir(self.output_dir)
        self.end()

    def pack_and_compress(self):
        self.logger.warn('start packing folders and compressing them, it will take a long time')
        os.chdir(self.output_dir)
        dirs = os.walk(self.output_dir).next()[1]
        for dirname in dirs:
            tarname = '{}.tar.gz'.format(dirname)
            if dirname in ('file_check', 'rnaseq_mapping'):
                continue
            if not os.path.isfile(os.path.join(self.output_dir, tarname)):
                with tarfile.open(tarname, 'w:gz') as tar:
                    tar.add(dirname)
            shutil.rmtree(dirname)
        os.chdir(self.work_dir)
        self.logger.warn('succeed in packing folders and compressing them')

    def end(self):
        super(CircrnaWorkflow, self).end()

    def export_task_info(self):
        api = self.api.api('whole_transcriptome.task_info')
        api.add_task_info(os.path.join(self.work_dir, 'data.json'))

    def export_genome_info(self):
        api = self.api.api('whole_transcriptome.genome_info')
        file_path = os.path.join(self.db_path, self.genome_doc['gene_stat'])
        species_name = self.option('organism_name')
        ref_anno_version = self.genome_doc['assembly']
        hyperlink = self.genome_doc['ensemble_web']
        api.add_genome_info(file_path, species_name, ref_anno_version, hyperlink)

    def export_qc(self):
        api = self.api.api('whole_transcriptome.qc')
        sample_list = os.path.join(self.option('fastq_dir').path, 'list.txt')
        group_file = self.option('group_table').path
        compare_file = self.option('control_file').path
        fq_type = self.option('fq_type')
        qc_stat_before = self.modules['hiseq_reads_stat_raw'].output_dir
        qc_stat_after = self.modules['hiseq_reads_stat_use'].output_dir
        api.add_sample_info(sample_list=sample_list, library='circle')
        group_id, specimen_names, category_names = api.add_sample_group(group_file=group_file, library='circle')
        control_id, compare_names = api.add_group_compare(compare_file=compare_file, library='circle',
                                                          group_id=group_id)
        qc_id = api.add_qc(fq_type=fq_type, library='circle')
        api.add_qc_detail(qc_id, qc_stat_before, qc_stat_after, 'circle')
        api.add_qc_graph(qc_id, qc_stat_before, 'circle', 'before')
        api.add_qc_graph(qc_id, qc_stat_after, 'circle', 'after')

    def export_qc_after(self):
        api = self.api.api('whole_transcriptome.qc')
        sample_list = os.path.join(self.option('fastq_dir').path, 'list.txt')
        group_file = self.option('group_table').path
        compare_file = self.option('control_file').path
        fq_type = self.option('fq_type')
        stat_output_dir = self.modules['hiseq_reads_stat_use'].output_dir
        api.add_sample_info(sample_list=sample_list, library='circle')
        group_id, specimen_names, category_names = api.add_sample_group(group_file=group_file, library='circle')
        control_id, compare_names = api.add_group_compare(compare_file=compare_file, library='circle',
                                                          group_id=group_id)
        api.add_qc(fq_type=fq_type, library='circle')
        api.add_qc_detail_after(qc_id, stat_output_dir)
        api.add_qc_graph(qc_id, stat_output_dir, 'circle', 'after')

    def export_mapping(self, assess=True):
        api = self.api.api('whole_transcriptome.mapping')
        stat_file = os.path.join(self.modules['rnaseq_mapping'].output_dir, 'stat')
        method = self.option('align_method')
        api.add_mapping_stat(stat_file=stat_file, library='circle', method=method)
        if assess:
            chr_dir = os.path.join(self.modules['map_assessment'].output_dir, 'chr_stat')
            cov_dir = os.path.join(self.modules['map_assessment'].output_dir, 'coverage')
            dis_dir = os.path.join(self.modules['map_assessment'].output_dir, 'distribution')
            sat_dir = os.path.join(self.modules['map_assessment'].output_dir, 'saturation')
            params = json.dumps({'task_id': self.task_id, 'submit_location': 'mapping', 'task_type': 2}, sort_keys=True)
            api.add_chrom_distribution_table(chr_dir, params=params, library='circle')
            api.add_coverage_table(cov_dir, params=params, detail=True, library='circle')
            api.add_distribution_table(dis_dir, params=params, library='circle')
            api.add_rpkm_table(sat_dir, params=params, detail=True, library='circle')

    def export_expression_stat_circ(self):
        api = self.api.api('whole_transcriptome.expression_statcirc')
        stat = os.path.join(self.tools['exp_makecirc'].output_dir, 'count/C.stat.txt')
        api.add_expression_statcirc(stat, self.task_id, self.project_sn)

    def export_diff_exp_stat_circ(self):
        api = self.api.api('whole_transcriptome.diff_exp_statcirc')
        map_dict = {
            'control': self.option('control_file').path,
            'c_output_dir': os.path.join(self.modules['diff_exp_c'].output_dir)
        }
        if self.option('diff_method').lower() in ['degseq', 'deseq2', 'edger', 'limma']:
            arg_dict = {
                'program': self.option('diff_method'),
                'fc': self.option('fc'),
                'qvalue': self.option('diff_fdr_ci'),
                'method': self.option('padjust_way')
            }
        if self.option('diff_method').lower() in ['noiseq']:
            arg_dict = {
                'program': self.option('diff_method'),
                'fc': self.option('fc'),
                'prob': self.option('diff_fdr_ci'),
            }
        api.add_diff_exp_statcirc(map_dict, arg_dict, self.task_id, self.project_sn, 'circ')

    def export_email(self):
        if self.modules['diff_exp_c'].option('email'):
            database = Config().get_mongo_client(mtype='project')[Config().get_mongo_dbname('project')]
            collection = database['sg_task_email']
            collection.update({'task_id': self.task_id}, {'$set': {'status': '5'}}, upsert=True)


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_ath(self):
        from mbio.workflows.whole_transcriptome.circrna import CircrnaWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.longrna',
            'options': {
                'datatype': 'rawdata',
                'fq_type': 'PE',
                'strand_specific': 'True',
                # 'strand_dir': 'RF',
                # 'qc_soft': 'fastp',
                'quality_score_system': 'phred_33',
                # 'is_duplicate': 'True',
                'fastq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190813056/fastq_dir/raw_data',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190813056/group_table/group.txt',
                'control_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190813056/control_file/control.txt',
                'taxonmy': 'Plant',
                'organism_name': 'Arabidopsis_thaliana',
                'genome_id': 'GM0348',
                'circ_method': 'ciri2'
                # 'nr_database': 'Plant',
                # 'kegg_database': 'Plant',
            }
        }
        wsheet = Sheet(data=data)
        wf = CircrnaWorkflow(wsheet)
        wf.sheet.id = 'longrna'
        wf.sheet.project_sn = 'longrna'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ath')])
    unittest.TextTestRunner(verbosity=2).run(suite)
