#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == 'gaohao'

import os
import re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
from biocluster.config import Config
import json
from collections import defaultdict
import time
import datetime
import shutil


class FungiModule(Module):
    """
    真菌基因组单个样品的工作流
    """
    def __init__(self, work_id):
        super(FungiModule, self).__init__(work_id)
        options = [
            {"name": "raw_dir", "type": "infile", "format": "bacgenome.raw_dir"},###rawdata的文件目录
            {"name": "asse_dir", "type": "infile", "format": "sequence.fasta_dir"},  ###assemble的文件目录
            {"name": "sample_name", "type": "string"},###样品名称
            {"name": "analysis", "type": "string", "default": "uncomplete"},###流程分析模式complete，uncomplete
            {"name": "busco_database", "type": "string"},
            {"name": "ref_protein", "type": "string"},
            {'name': 'qc_tool', 'type': 'string', 'default': 'fastp'}  # 质控软件选择
        ]
        self.add_option(options)
        self.sequence = self.add_module('fungi_genome.fungi_genome') #将数据解压、合并统计
        self.fungi_qc = self.add_module('fungi_genome.fungi_genome_qc')  #将二代数据质控
        self.assemble = self.add_module('fungi_genome.bac_assemble')  # 二代数据组装
        self.assemble_com = self.add_tool('fungi_genome.canu_assemble') #三代组装
        self.assemble_assess = self.add_module('fungi_genome.assemble_assess')  # 组装数据评估
        self.assemble_busco = self.add_tool('fungi_genome.busco')  # 组装数据busco评估
        self.assemble_cegma = self.add_tool('fungi_genome.cegma')  # 组装数据cegma评估
        self.genome_assess =self.add_module('fungi_genome.genome_assess')  #基因组大小评估
        self.repeat = self.add_tool('predict.repeatmasker_and_trf')  # 重复序列预测
        self.predict_gene =self.add_module('fungi_genome.fungi_predict') #基因预测，非编码预测
        self.anno =self.add_module('fungi_genome.annotation') #六种基因数据库注释
        self.cyps_anno = self.add_module('annotation.cyps_anno')#p450
        self.promote = self.add_tool('gene_structure.promote')#启动子生成
        self.summary =self.add_tool('bacgenome.bac_summary') #项目总览统计
        self.pathogenic_system = self.add_module('fungi_genome.pathogenic_system') #致病系统注释
        self.gbk_file =self.add_tool('fungi_genome.bac_gbk') #生成gbk文件
        self.step.add_steps('sequence', 'fungi_qc','assemble','assemble_busco','assemble_cegma','assemble_assess','genome_assess','repeat','predict_gene','annotation','promote','pathogenic_system','summary','gbk','cyps_anno')
        self.seq_type = ''
        self.list = [self.repeat,self.assemble_cegma,self.assemble_busco]
        self.list2 = [self.predict_gene,self.genome_assess]
        self.anno_list = [self.anno,self.cyps_anno,self.pathogenic_system]
        self.sum_list = [self.promote,self.summary,self.gbk_file]

    def check_options(self):
        """
        检查参数
        """
        if not self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
            raise OptionError("请提供raw数据或assemble数据！不能同时为空！", code="22100601")
        if  not self.option("sample_name"):
            raise OptionError("必须输入样本名!", code="22100602")

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

    def run(self):
        super(FungiModule, self).run()
        self.get_specimen()
        if not os.path.exists(self.output_dir + '/' + self.option("sample_name")):
            os.mkdir(self.output_dir + '/' + self.option("sample_name"))
        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/specimen.data'):
            os.remove(self.output_dir + '/' + self.option("sample_name") + '/specimen.data')
        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/gene.data'):
            os.remove(self.output_dir + '/' + self.option("sample_name") + '/gene.data')
        os.link(self.work_dir + '/specimen.data',self.output_dir + '/' + self.option("sample_name") + '/specimen.data')
        os.link(self.work_dir + '/gene.data',self.output_dir + '/' + self.option("sample_name") + '/gene.data')
        if self.option("raw_dir").is_set:
            raw_path = self.option("raw_dir").prop['path'] + '/' + 'list.txt'
            self.seq_type = self.get_sequence_type(raw_path)
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                if self.seq_type in ['PE','PE,MP','MP,PE']:
                    self.on_rely(self.sum_list, self.end)
                    self.on_rely(self.anno_list,self.run_su)
                    self.on_rely(self.list2,self.run_anno)
                    self.on_rely(self.list, self.run_gene)
                    self.assemble_assess.on("end", self.run_assemble_busco)
                    self.assemble_assess.on("end", self.run_assemble_cegma)
                    self.assemble_assess.on("end", self.run_repeat)
                    self.assemble.on("end", self.run_assess)
                    self.fungi_qc.on("end", self.run_assemble)
                    self.sequence.on("end", self.run_fungi_qc)
                    self.run_sequence()

                elif self.seq_type in ['pacbio','Pacbio','PACBIO']:
                    self.on_rely(self.sum_list, self.end)
                    self.on_rely(self.anno_list, self.run_su)
                    self.predict_gene.on("end", self.run_anno)
                    self.on_rely(self.list, self.run_predict)
                    self.assemble_assess.on("end", self.run_assemble_busco)
                    self.assemble_assess.on("end", self.run_assemble_cegma)
                    self.assemble_assess.on("end", self.run_repeat)
                    self.assemble_com.on("end", self.run_assess)
                    self.fungi_qc.on("end", self.run_assemble)
                    self.sequence.on("end", self.run_fungi_qc)
                    self.run_sequence()

                elif self.seq_type in ['PE,pacbio','pacbio,PE','PE,Pacbio','Pacbio,PE']:
                    self.on_rely(self.sum_list, self.end)
                    self.on_rely(self.anno_list,self.run_su)
                    self.on_rely(self.list2,self.run_anno)
                    self.on_rely(self.list, self.run_gene)
                    self.assemble_assess.on("end", self.run_assemble_busco)
                    self.assemble_assess.on("end", self.run_assemble_cegma)
                    self.assemble_assess.on("end", self.run_repeat)
                    self.assemble_com.on("end", self.run_assess)
                    self.fungi_qc.on("end", self.run_assemble)
                    self.sequence.on("end", self.run_fungi_qc)
                    self.run_sequence()

            elif not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                self.on_rely(self.sum_list, self.end)
                self.on_rely(self.anno_list, self.run_su)
                self.predict_gene.on("end", self.run_anno)
                self.on_rely(self.list, self.run_predict)
                self.assemble_assess.on("end", self.run_assemble_busco)
                self.assemble_assess.on("end", self.run_assemble_cegma)
                self.assemble_assess.on("end", self.run_repeat)
                self.run_assess()

            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                if self.seq_type in ['PE','PE,MP','MP,PE']:
                    self.on_rely(self.sum_list, self.end)
                    self.on_rely(self.anno_list,self.run_su)
                    self.on_rely(self.list2,self.run_anno)
                    self.on_rely(self.list, self.run_gene)
                    self.assemble_assess.on("end", self.run_assemble_busco)
                    self.assemble_assess.on("end", self.run_assemble_cegma)
                    self.assemble_assess.on("end", self.run_repeat)
                    self.fungi_qc.on("end", self.run_assess)
                    self.sequence.on("end", self.run_fungi_qc)
                    self.run_sequence()

                elif self.seq_type in ['pacbio','Pacbio','PACBIO']:
                    self.on_rely(self.sum_list, self.end)
                    self.on_rely(self.anno_list, self.run_su)
                    self.predict_gene.on("end", self.run_anno)
                    self.on_rely(self.list, self.run_predict)
                    self.assemble_assess.on("end", self.run_assemble_busco)
                    self.assemble_assess.on("end", self.run_assemble_cegma)
                    self.assemble_assess.on("end", self.run_repeat)
                    self.fungi_qc.on("end", self.run_assess)
                    self.sequence.on("end", self.run_fungi_qc)
                    self.run_sequence()

                elif self.seq_type in ['PE,pacbio','pacbio,PE','PE,Pacbio','Pacbio,PE']:
                    self.on_rely(self.sum_list, self.end)
                    self.on_rely(self.anno_list,self.run_su)
                    self.on_rely(self.list2,self.run_anno)
                    self.on_rely(self.list, self.run_gene)
                    self.assemble_assess.on("end", self.run_assemble_busco)
                    self.assemble_assess.on("end", self.run_assemble_cegma)
                    self.assemble_assess.on("end", self.run_repeat)
                    self.fungi_qc.on("end", self.run_assess)
                    self.sequence.on("end", self.run_fungi_qc)
                    self.run_sequence()

    def run_anno(self):
        self.run_annotation()
        self.run_pathogenic_system()
        self.run_cyps()

    def run_su(self):
        self.run_promote()
        self.run_summary()
        self.run_gbk()

    def run_gene(self):
        self.run_predict()
        self.run_genome_assess()

    def run_sequence(self):
        opts = ''
        if self.option('raw_dir').is_set:
            opts = {
                'raw_dir': self.option('raw_dir'),
                'sample_name': self.option('sample_name'),
                'sequence_type': self.seq_type,
            }
        self.set_run(opts, self.sequence, 'sequence', self.step.sequence)

    def run_fungi_qc(self):
        #pacbio_dir_path = self.sequence.option('pacbio_dir').prop['path']
        #dir_path = self.sequence.option('dir').prop['path']
        opts = ''
        if  self.sequence.option('dir').is_set and not self.sequence.option('pacbio_dir').is_set:
            opts = {
                'analysis': self.option('analysis'),
                'fastq_dir': self.sequence.option('dir'),
                'sample_name': self.option('sample_name'),
                'qc_tool': self.option('qc_tool')
            }
        #elif os.path.exists(pacbio_dir_path) and not os.path.exists(dir_path):
        elif self.sequence.option('pacbio_dir').is_set and not self.sequence.option('dir').is_set:
            opts = {
                'analysis': self.option('analysis'),
                'pacbio_dir': self.sequence.option('pacbio_dir'),
                'sample_name': self.option('sample_name'),
                'qc_tool': self.option('qc_tool')
            }
        #elif os.path.exists(pacbio_dir_path) and os.path.exists(dir_path):
        elif self.sequence.option('pacbio_dir').is_set and self.sequence.option('dir').is_set:
            opts = {
                'analysis': self.option('analysis'),
                'fastq_dir': self.sequence.option('dir'),
                'sample_name': self.option('sample_name'),
                'pacbio_dir':self.sequence.option('pacbio_dir'),
                'qc_tool': self.option('qc_tool')
            }
        self.set_run(opts, self.fungi_qc, 'fungi_qc', self.step.fungi_qc)

    def run_assemble(self):
        seq_type =self.get_assemble_type(self.option("raw_dir").prop['path'] + '/' + 'list.txt')
        self.logger.info(seq_type)
        if self.option('analysis') in ['uncomplete'] and not re.search(r'pacbio',self.seq_type):
            clean_dir = self.fungi_qc.output_dir + '/cleandata'
            sample_info = self.sequence.output_dir + '/data/' + 'sample_info'
            opts = {
                'fq_dir': clean_dir,
                'seq_type': seq_type,
                'sample_name': self.option('sample_name'),
                'sample_info': sample_info,
            }
            self.set_run(opts, self.assemble, 'assemble', self.step.assemble)
        if self.option('analysis') in ['uncomplete'] and re.search(r'pacbio',self.seq_type):
            clean_fq = self.sequence.option('pacbio_dir').prop['path'] + '/all.pacbio.fq'
            size = self.get_size(self.sequence.option('pacbio_dir').prop['path'] + '/pacbio.rawdata.list')
            self.logger.info(self.sequence.option('pacbio_dir').prop['path'] + '/pacbio.rawdata.list')
            self.logger.info(size)
            opts = {
                'subread_fq': clean_fq,
                'genomeSize': size,
                'sample_name': self.option('sample_name'),
            }
            self.set_run(opts, self.assemble_com, 'assemble', self.step.assemble)

    def run_assess(self):
        seq_path = ''
        if self.option('asse_dir').is_set:
            seq_path =self.get_seq()
        elif self.option('raw_dir').is_set  and not re.search(r'pacbio',self.seq_type):
            seq_path =self.assemble.option('scaffold')
        elif self.option('raw_dir').is_set and re.search(r'pacbio',self.seq_type):
            seq_path =self.assemble_com.output_dir + '/' + self.option('sample_name') + '.contigs.fasta'
        opts = {
            'seq_scaf': seq_path,
            'sample_name': self.option('sample_name'),
        }
        self.set_run(opts, self.assemble_assess, 'assemble_assess', self.step.assemble_assess)

    def run_genome_assess(self):
        bases = self.get_bases()
        opts = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess.option('scaffold')
            opts = {
                'seq': scaf_seq,
                'fastq_dir': self.sequence.option('dir'),
                'bases': bases,
                'sample_name': self.option('sample_name'),
            }
        self.set_run(opts, self.genome_assess, 'genome_assess', self.step.genome_assess)

    def run_repeat(self):
        scaf_seq = self.assemble_assess.option('scaffold')
        opts = {
            'input_genome': scaf_seq,
            'analysis_type': '2',
        }
        self.set_run(opts, self.repeat, 'repeat', self.step.repeat)

    def run_assemble_busco(self):
        scaf_seq = self.assemble_assess.option('scaffold')
        opts = {
            'scaf_fa': scaf_seq,
            'sample_name': self.option('sample_name'),
            'database':self.option('busco_database'),
        }
        self.set_run(opts, self.assemble_busco, 'assemble_busco', self.step.assemble_busco)

    def run_assemble_cegma(self):
        scaf_seq = self.assemble_assess.option('scaffold')
        opts = {
            'scaf_fa':scaf_seq,
            'sample_name': self.option('sample_name'),
        }
        self.set_run(opts, self.assemble_cegma, 'assemble_cegma', self.step.assemble_cegma)

    def run_predict(self):
        scaf_seq = self.assemble_assess.option('scaffold')
        opts = {
            'masked':self.repeat.work_dir + '/' + self.option('sample_name') + '_scaf.fna.masked',
            'fasta': scaf_seq,
            'sample': self.option('sample_name'),
            'cegma_gff':self.assemble_cegma.work_dir + '/' + self.option('sample_name') + '.cegma.gff',
            'ref_protein':self.option('ref_protein'),
        }
        self.set_run(opts, self.predict_gene, 'predict_gene', self.step.predict_gene)

    def run_annotation(self):
        opts = {
            'gene_seq': self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS.faa',
            'gene_gff': self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS.gff',
            'database_list': 'nr_v20200604,swissprot_v20200617,pfam_v33.1,eggnog,kegg_v94.2,cazy_v8',
            'sample': self.option('sample_name'),
        }
        self.set_run(opts, self.anno, 'annotation', self.step.annotation)

    def run_cyps(self):
        opts = {
            'query': self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS.faa',
            'sample': self.option('sample_name'),
        }
        self.set_run(opts, self.cyps_anno, 'cyps_anno', self.step.cyps_anno)

    def run_summary(self):
        gene_stat = ''
        asse = ''
        if self.option('analysis') in ['uncomplete']:
            gene_stat =self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS_statistics.xls'
            asse = self.assemble_assess.output_dir + '/assembly/' + self.option('sample_name') + '_assembly_summary.xls'
        opts = {
            'gene_statistics': gene_stat,
            'rrna_gff': self.predict_gene.output_dir + '/rRNA/' + self.option('sample_name') + '_rRNA.gff',
            'trna_gff': self.predict_gene.output_dir + '/tRNA/' + self.option('sample_name') + '_tRNA.gff',
            'assemble': asse,
            'cog': self.anno.output_dir + '/COG/' + self.option('sample_name') + '_cog_anno.xls',
            'kegg': self.anno.output_dir + '/KEGG/' + self.option('sample_name') + '_kegg_anno.xls',
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.summary, 'summary', self.step.summary)

    def run_gbk(self):
        """
        细菌v3增加此功能，生成合并的gbk文件和gff文件
        :return:
        """
        self.logger.info(">>>start run_gbk")
        gene_gff = self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS.gff'
        rrna_gff = self.predict_gene.output_dir + '/rRNA/' + self.option('sample_name') + '_rRNA.gff'
        trna_gff = self.predict_gene.output_dir + '/tRNA/' + self.option('sample_name') + '_tRNA.gff'
        protein_seq = self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS.faa'
        anno_summary = self.anno.output_dir + '/Summary/' + self.option('sample_name') + '_anno_summary.xls'
        scaf_seq = self.assemble_assess.output_dir + '/assembly/' + self.option('sample_name') + '_scaf.fna'
        opts = {
            'gen_gff': gene_gff,
            'rrna_gff': rrna_gff,
            'trna_gff': trna_gff,
            'pro_fa': protein_seq,
            'genome_fa': scaf_seq,
            'anno': anno_summary,
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.gbk_file, 'gbk', self.step.gbk)

    def run_promote(self):
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq = self.assemble_assess.work_dir + '/all.fasta'
        opts = {
            'sequence': self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS.fnn',
            'assemble': scaf_seq,
            'sample': self.option('sample_name'),
            'pro_tidy_type': 1   #zouguanqing 20180713
        }
        self.set_run(opts, self.promote, 'promote', self.step.promote)

    def run_pathogenic_system(self):
        opts = {
            'query':  self.predict_gene.output_dir + '/CDS_predict/' + self.option('sample_name') + '_CDS.faa',
            'sample': self.option('sample_name'),
        }
        self.set_run(opts, self.pathogenic_system, 'pathogenic_system', self.step.pathogenic_system)

    def set_output(self,event):
        """
        将各个模块的结果输出至output
        """
        if event['data'] == 'fungi_qc':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/fastx'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/fastx')
            shutil.copytree(self.fungi_qc.output_dir + '/fastx',
                            self.output_dir + '/' + self.option("sample_name") + '/fastx')
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/data_QC'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/data_QC')
            shutil.copytree(self.fungi_qc.output_dir + '/data_QC',
                            self.output_dir + '/' + self.option("sample_name") + '/data_QC')
        if event['data'] == 'genome_assess':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
            shutil.copytree(self.genome_assess.output_dir,
                            self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
        if event['data'] == 'assemble_assess':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
            shutil.copytree(self.assemble_assess.output_dir,
                            self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
        if not os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment'):
            os.makedirs(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment')
        if event['data'] == 'assemble_busco':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment/' + self.option("sample_name") + '_busco.xls'):
                os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment/' + self.option("sample_name") + '_busco.xls')
            os.link(self.assemble_busco.work_dir + '/' + self.option("sample_name") + '_busco.xls',
                    self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment/' + self.option("sample_name") + '_busco.xls')
        if event['data'] == 'assemble_cegma':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment/' + self.option("sample_name") + '_cegma.xls'):
                os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment/' + self.option("sample_name") + '_cegma.xls')
            os.link(self.assemble_cegma.work_dir + '/' + self.option("sample_name") + '_cegma.xls',
                    self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly_assessment/' + self.option("sample_name") + '_cegma.xls')
        if event['data'] == 'predict_gene':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/predict'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/predict')
            shutil.copytree(self.predict_gene.output_dir,
                            self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/predict')
        if event['data'] == 'repeat':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/repeats'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/repeats')
            shutil.copytree(self.repeat.output_dir ,
                            self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/repeats')
        if event['data'] == 'annotation':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/annotation'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/annotation')
            shutil.copytree(self.anno.output_dir, self.output_dir + '/' + self.option("sample_name") + '/annotation')
        if event['data'] == 'cyps_anno':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/P450'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/P450')
            shutil.copytree(self.cyps_anno.output_dir + '/P450',self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/P450')
        if event['data'] == 'promote':
            if os.path.exists(self.output_dir + '/' + self.option(
                    "sample_name") + '/structral_genome/promoter_predict'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
            shutil.copytree(self.promote.output_dir,self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
        if event['data'] == 'summary':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/project_overview'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/project_overview')
            shutil.copytree(self.summary.output_dir,self.output_dir + '/' + self.option("sample_name") + '/project_overview')
            os.rename(self.output_dir + '/' + self.option("sample_name") + '/project_overview/' + self.option("sample_name") + '.project.summary',
                      self.output_dir + '/' + self.option("sample_name") + '/project_overview/genome_overview.xls')
        if event['data'] == 'pathogenic_system':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
            shutil.copytree(self.pathogenic_system.output_dir,self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
        ## gbk文件
        if event['data'] == 'gbk':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/project_overview/All_predict_info'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/project_overview/All_predict_info')
            shutil.copytree(self.gbk_file.output_dir,self.output_dir + '/' + self.option("sample_name") + '/project_overview/All_predict_info')

    def end(self):
        super(FungiModule, self).end()

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

    def move_dir(self, olddir, newname):  # 原函数名move2outputdir
        """
        移动一个目录下所有文件/文件夹到workflow输出路径下，供set_output调用
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="22100601")
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

    def get_sequence_type(self,file):
        dict ={}
        type =[]
        with open(file, "rb") as l:
            raw_lines = l.readlines()
            for line in raw_lines[1:]:
                line2 = line.strip('\r\n').split("\t")
                dict[line2[5]] = line2[5]
        for k in dict.iterkeys():
            type.append(k)
        if len(type) ==1:
            sequence_type=type[0]
        else:
            sequence_type = ','.join(type)
        return sequence_type

    def get_assemble_type(self,file):
        dict = {}
        type = []
        with open(file, "rb") as l:
            raw_lines = l.readlines()
            for line in raw_lines[1:]:
                line2 = line.strip('\r\n').split("\t")
                if line2[5] not in ['pacbio']:
                    dict[line2[5]] = line2[5]
        for k in dict.iterkeys():
            type.append(k)
        if len(type) == 1:
            sequence_type = type[0]
        else:
            sequence_type = ','.join(type)
        return sequence_type

    def get_size(self,file):
        with open(file,'r') as f:
            lines = f.readlines()
            size =lines[0].rstrip('\r\n').split('\t')[4]
        return size

    def get_seq(self):
        path=''
        if self.option('analysis') in ['uncomplete'] and self.option('asse_dir').is_set:
            list=self.option('asse_dir').prop['path'] + '/list.txt'
            with open(list,'r') as f:
                lines=f.readlines()
                for line in lines[1:]:
                    line = line.rstrip('\n\r').split('\t')
                    path=self.option('asse_dir').prop['path'] + '/' + line[1]
                    print(path)
        return path

    def get_gene_prefix(self):
        de =''
        prefix ={}
        if self.option('analysis') in ['complete']:
            file =self.assemble_assess.work_dir + '/plasmid.type.xls'
            with open(file,'r') as f:
                lines =f.readlines()
                for line in lines[0:]:
                    line =line.rstrip('\r\n').split('\t')
                    prefix[line[0]]=line[1]
            f.close()
            de =str(prefix)
        return de

    def get_bases(self):
        tatol =0
        base_file =self.fungi_qc.output_dir + '/data_QC/' + self.option('sample_name') + '_Illumina_statistics.xls'
        self.logger.info(base_file)
        if os.path.exists(base_file):
            with open(base_file,'r') as f:
                lines =f.readlines()
                for line in lines[1:]:
                    line2 =line.rstrip('\r\n').split('\t')
                    if int(eval(line2[1])) < 1000:
                        self.logger.info(line2[4])
                        tatol += int(eval(line2[4]))
        self.logger.info(tatol)
        return tatol

    def get_seq_type(self,file_dir):
        list =[]
        files =os.listdir(file_dir)
        for file in files:
            lst=file.split('.')
            list.append(lst[0])
        if len(list) == 1:
            seq_type =list[0]
        else:
            seq_type = ','.join(list)
        return seq_type

    def get_specimen(self):
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                output3 = self.work_dir + "/" + 'specimen.data'
                output4 = self.work_dir + "/" + 'gene.data'
                with open(raw_list_path, "rb") as l, open(output3, 'w') as file3, open(
                        output4, 'w') as file4:
                    lines = l.readlines()
                    lib_list = []
                    raw_list = []
                    for line in lines[1:]:
                        line2 = line.strip().split()
                        if len(line2) == 6:
                            lib_type = line2[5] + line2[2]
                            lib_list.append(lib_type)
                            raw_list.append(line2[1])
                    if len(lib_list) == 1:
                        lib = lib_list[0]
                    else:
                        lib = ','.join(lib_list)
                    if len(raw_list) == 1:
                        raw_file = raw_list[0]
                    else:
                        raw_file = ';'.join(raw_list)
                    file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                    file3.write(self.option('sample_name') + '\t' + raw_file + '\t' + lib + '\n')
                    file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                    file4.write(self.option('sample_name') + '\t' + '-' + '\t' + '-' + '\t' + 'gene' + '\n')
            elif not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                ass_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                output3 = self.work_dir + "/" + 'specimen.data'
                output4 = self.work_dir + "/" + 'gene.data'
                with open(ass_list_path, "rb") as l, open(output3, 'w') as file3, open(output4, 'w') as file4:
                    file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                    file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                    lines = l.readlines()
                    for line in lines[1:]:
                        line2 = line.strip('\r\n').split()
                        if len(line2) == 2:
                            file3.write(self.option('sample_name') + '\t' + '-' + '\t' + '-' + '\n')
                            file4.write(self.option('sample_name') + '\t' + line2[1] + '\t' + '-' + '\t' + 'gene' + '\n')

            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                ass_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                output3 = self.work_dir + "/" + 'specimen.data'
                output4 = self.work_dir + "/" + 'gene.data'
                with open(raw_list_path, "rb") as l,open(ass_list_path, "rb") as f, open(output3, 'w') as file3, open(
                        output4, 'w') as file4:
                    lines = l.readlines()
                    lib_list = []
                    raw_list = []
                    for line in lines[1:]:
                        line2 = line.strip().split()
                        if len(line2) == 6:
                            lib_type = line2[5] + line2[2]
                            lib_list.append(lib_type)
                            raw_list.append(line2[1])
                    lines2 =f.readlines()
                    assem = ''
                    for line in lines2[1:]:
                        line2 = line.strip().split()
                        if len(line2) == 2:
                            assem =line2[1]
                    if len(lib_list) == 1:
                        lib = lib_list[0]
                    else:
                        lib = ','.join(lib_list)
                    if len(raw_list) == 1:
                        raw_file = raw_list[0]
                    else:
                        raw_file = ';'.join(raw_list)
                    file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                    file3.write(self.option('sample_name') + '\t' + raw_file + '\t' + lib + '\n')
                    file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                    file4.write(self.option('sample_name') + '\t' + assem + '\t' + '-' + '\t' + 'gene' + '\n')