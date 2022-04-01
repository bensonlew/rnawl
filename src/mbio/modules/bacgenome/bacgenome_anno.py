#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == 'zouguanqing'

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


class BacgenomeAnnoModule(Module):
    """
    微生物基因组单个样品的注释工作流
    """
    def __init__(self, work_id):
        super(BacgenomeAnnoModule, self).__init__(work_id)
        options = [
            #{"name": "anno_dir", "type": "infile", "format":"bacgenome.anno_dir"},
            {"name": "gff_file", "type": "string", "default":""},
            {"name": "all_fasta", "type": "string", "default":""},
            {"name": "seq_dir", "type": "string"},
            {"name": "sample_name", "type": "string"}, ###样品名称
            {"name": "analysis", "type": "string", "default": "uncomplete"}, ###流程分析模式complete，uncomplete
            {"name": "nr_evalue","type":"string","default":"1e-5"},
            {"name": "swissprot_evalue","type":"string","default":"1e-5"},
            {"name": "cog_evalue","type":"string","default":"1e-5"},
            {"name": "kegg_evalue","type":"string","default":"1e-5"},
            {"name": "go_evalue","type": "string", "default":"1e-5"},
            {"name": "txt_info","type":"string","default":""}  # 传过来的额外的信息，例如"{plasmid_num:2,chr_num:1}"
        ]
        self.add_option(options)

        self.pre_anno = self.add_tool('bacgenome.gff_deal')
        self.anno =self.add_module('bacgenome.annotation') #六种基因数据库注释
        self.gbk_file =self.add_tool('bacgenome.bac_gbk')  #gbk文件生成
        self.promote = self.add_tool('gene_structure.promote')#启动子生成
        self.island = self.add_module('bacgenome.island') #基因组岛预测
        self.crispr = self.add_tool('bacgenome.crisprcas') #CRISPR的预测
        self.prephage =self.add_tool('bacgenome.prephage') #prephage的预测
        self.cgview = self.add_tool('bacgenome.cgview') #cgview画图
        self.tree = self.add_module('bacgenome.bac_tree') #进化树
        self.circos = self.add_module('bacgenome.combine_circos')#circos画图
        self.pathogenic_system = self.add_module('annotation.pathogenic_system') #致病系统注释
        self.antismash = self.add_tool('annotation.antismash_v5')
        self.antismash_com = self.add_module('bacgenome.antismash')
        self.overview = self.add_tool('bacgenome.anno_summary')  #zouguanqing
        self.is_predict = self.add_module("bacgenome.is_predict")  #qingchen.zhang ## is可移动元件预测
        self.integron = self.add_module("bacgenome.integron_predict")  #qingchen.zhang ## 整合子可移动元件预测
        self.repeatmasker = self.add_module('bacgenome.repeatmasker')#qingchen.zhang ## 散在重复序列预测
        self.gbk_real_file =self.add_tool('bacgenome.bac_real_gbk')  #gbk文件生成

        self.step.add_steps('pre_anno','annotation','promote','gbk','cgview','circos','antismash','crispr','prephage','island','pathogenic_system','overview','tree', 'is_predict', 'integron', "repeatmasker", 'real_gbk') #,
        self.seq_type = ''

        self.complete_anno =[self.pathogenic_system,self.anno,self.tree] #



    def check_options(self):
        """
        检查参数
        """
        if not self.option("analysis"):
            raise OptionError("必须输入流程运行类型！")
        if self.option("analysis") not in ['complete','uncomplete']:
            raise OptionError("请提分析流程参数！")


        if  not self.option("sample_name"):
            raise OptionError("必须输入样本名！")

        if self.option('analysis') in ['uncomplete']:
            self.ana2 = [self.circos, self.prephage, self.crispr,  self.cgview, self.promote, self.island, self.antismash, self.is_predict, self.integron, self.repeatmasker, self.gbk_real_file]
        else:
            self.ana2 = [self.circos, self.prephage, self.crispr,  self.cgview, self.promote, self.island, self.antismash_com, self.is_predict, self.integron, self.repeatmasker, self.gbk_real_file]

        self.fna = self.option('all_fasta')
        self.ori_gff = self.option('gff_file')
        self.seq_dir =  self.option('seq_dir')

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_run(self, opts, module, event, step, start=True):
        self.logger.info( "%s" % module._options)
        module.set_options(opts)
        module.on('start', self.set_step, {'start': step})
        module.on('end', self.set_step, {'end': step})
        module.on('end', self.set_output, event)
        if start:
            module.run()

    def run(self):
        super(BacgenomeAnnoModule, self).run()
        #self.get_specimen()
        # if not os.path.exists(self.output_dir + '/' + self.option("sample_name")):
        #     os.mkdir(self.output_dir + '/' + self.option("sample_name"))
        # if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/specimen.data'):
        #     os.remove(self.output_dir + '/' + self.option("sample_name") + '/specimen.data')
        # if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/gene.data'):
        #     os.remove(self.output_dir + '/' + self.option("sample_name") + '/gene.data')
        #os.link(self.work_dir + '/specimen.data',self.output_dir + '/' + self.option("sample_name") + '/specimen.data')
        #os.link(self.work_dir + '/gene.data',self.output_dir + '/' + self.option("sample_name") + '/gene.data')

        self.gbk_file.on("end", self.run_antismash)
        self.gbk_file.on("end", self.run_circos)
        self.gbk_file.on("end", self.run_prephage)
        self.gbk_file.on("end", self.run_real_gbk)
        self.gbk_file.on("end", self.run_repeatmasker)
        self.gbk_file.on("end", self.run_is_predict)
        self.gbk_file.on("end", self.run_integron)
        self.gbk_file.on("end", self.run_island)
        self.gbk_file.on("end", self.run_cgview)
        self.gbk_file.on("end", self.run_crispr)
        self.gbk_file.on("end", self.run_promote)
        self.on_rely(self.complete_anno,self.run_gbk)  #并行跑完注释，致病，进化树后，跑gbk，然后并行跑ana2列表分析
        self.on_rely(self.ana2,self.run_overview)
        self.overview.on("end",self.end)
        self.pre_anno.on("end", self.run_annotation)
        self.pre_anno.on("end", self.run_pathogenic_system)
        self.pre_anno.on("end", self.run_tree)
        self.run_pre_anno()

    def run_pre_anno(self):

        opts = {
            'fna' : self.fna,
            'all_gff' : self.ori_gff,
            'sample_name': self.option('sample_name'),
            'txt_info':self.option('txt_info'),
            'analysis' : self.option("analysis")
        }
        self.set_run(opts, self.pre_anno,'pre_anno', self.step.pre_anno)


    def run_tree(self):
        if self.option('analysis') in ['uncomplete']:
            opts = {
                'seq_faa': self.pre_anno.option('faa'),
                'sample_name': self.option('sample_name'),
                'analysis': self.option('analysis'),
            }
            self.set_run(opts, self.tree, 'tree', self.step.tree)
        elif self.option('analysis') in ['complete']:
            seq_fa = self.pre_anno.option('faa')  #faa 序列
            scaf_seq = self.fna
            opts = {
                'seq_faa': seq_fa,
                'sample_name': self.option('sample_name'),
                'analysis': self.option('analysis'),
                'gene_gff': self.pre_anno.option('gene_gff'),
                'rrna_gff': self.pre_anno.option('rrna_gff'),
                'genome_fa':scaf_seq,
            }
            self.set_run(opts, self.tree, 'tree', self.step.tree)

    def run_annotation(self):
        opts = {
            'gene_seq': self.pre_anno.option('faa'),
            'gene_gff': self.pre_anno.option('gene_gff'),
            'database_list': 'nr,swissprot,pfam,eggnog,kegg,cazy',
            'sample': self.option('sample_name'),
            'analysis': self.option('analysis'),
            'has_two_component' : 'T',
            'produce_mark' : True
        }

        opts['nr_evalue'] = self.option('nr_evalue')
        opts['go_evalue'] = self.option('go_evalue')
        opts['cog_evalue'] = self.option('cog_evalue')
        opts['kegg_evalue'] = self.option('kegg_evalue')
        opts['swissprot_evalue'] = self.option('swissprot_evalue')

        self.set_run(opts, self.anno, 'annotation', self.step.annotation)



    def run_gbk(self):
        self.logger.info(">>>in run_gbk")
        scaf_seq = self.fna
        opts = {
            'gen_gff':    self.pre_anno.option('gene_gff'),
            'rrna_gff':   self.pre_anno.option('rrna_gff'),
            'trna_gff':   self.pre_anno.option('trna_gff'),
            'pro_fa': self.pre_anno.option('faa'),
            'genome_fa': scaf_seq,
            'anno': self.anno.option('tidy_summary'),
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.gbk_file, 'gbk', self.step.gbk)

    def run_real_gbk(self):
        """
        生成gbk文件，基因的名称为真实的基因的名称
        :return:
        """
        self.logger.info(">>>in run_gbk")
        scaf_seq = self.fna
        opts = {
            'gen_gff': self.pre_anno.option('gene_gff'),
            'rrna_gff': self.pre_anno.option('rrna_gff'),
            'trna_gff': self.pre_anno.option('trna_gff'),
            'pro_fa': self.pre_anno.option('faa'),
            'genome_fa': scaf_seq,
            'anno': self.anno.option('tidy_summary'),
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.gbk_real_file, 'real_gbk', self.step.real_gbk)

    def run_promote(self):
        scaf_seq = self.fna

        opts = {
            'sequence': self.pre_anno.option('ffn') ,
            'assemble': scaf_seq,
            'sample': self.option('sample_name'),
        }
        self.set_run(opts, self.promote, 'promote', self.step.promote)

    def run_crispr(self):
        scaf_seq = self.fna

        opts = {
            'seq': scaf_seq,
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.crispr, 'crispr', self.step.crispr)

    def run_prephage(self):
        scaf_seq = self.fna

        opts = {
            'prot_seq': self.pre_anno.option('faa'),
            'scaf_seq': scaf_seq,
            'gene_gff':   self.pre_anno.option('gene_gff'),
            'anno': self.anno.option('tidy_summary'),
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.prephage, 'prephage', self.step.prephage)

    def run_island(self):
        seq_dir = self.seq_dir

        opts = {
            'fa_dir': seq_dir,
            'gbk_dir': self.gbk_file.option('gbk_dir'),
            'anno': self.anno.option('tidy_summary'),
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.island, 'island', self.step.island)

    def run_is_predict(self):
        """
        is预测
        :return:
        """
        scaf_seq = self.fna
        opts = {
            'genome_fa': scaf_seq,
            'gene_faa': self.pre_anno.option('faa'),
            'gene_fna': self.pre_anno.option('ffn'),
            'gene_gff': self.pre_anno.option('gene_gff').prop['path'],
            'sample': self.option('sample_name')
        }
        self.set_run(opts, self.is_predict, 'is_predict', self.step.is_predict)

    def run_integron(self):
        """
        整合子预测
        :return:
        """
        scaf_seq = self.fna
        opts = {
            'genome_fa': scaf_seq,
            'gene_faa': self.pre_anno.option('faa'),
            'gene_gff': self.pre_anno.option('gene_gff').prop['path'],
            'sample': self.option('sample_name')
        }
        self.set_run(opts, self.integron, 'integron', self.step.integron)

    def run_repeatmasker(self):
        """
        散在重复预测
        :return:
        """
        scaf_seq = self.fna
        opts = {
            "genome_fa": scaf_seq,
            "sample": self.option('sample_name'),
            "analysis": self.option('analysis'),
        }
        self.set_run(opts, self.repeatmasker, 'repeatmasker', self.step.repeatmasker)

    def run_cgview(self):
        scaf_seq = self.fna

        opts = {
            'gen_gff':  self.pre_anno.option('gene_gff'),
            'rrna_gff': self.pre_anno.option('rrna_gff'),
            'trna_gff': self.pre_anno.option('trna_gff'),
            'pro_fa':  self.pre_anno.option('faa'),
            'genome_fa': scaf_seq,
            'anno': self.anno.option('tidy_summary'),
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.cgview, 'cgview', self.step.cgview)

    def run_circos(self):
        scaf_seq = self.fna
        if self.option('analysis') in ['uncomplete']:
            anno_cog =self.anno.output_dir + '/COG/' + self.option('sample_name') + '_cog_anno.xls'
            opts = {
                'gene':  self.pre_anno.option('gene_gff'),
                'rrna': self.pre_anno.option('rrna_gff'),
                'trna': self.pre_anno.option('trna_gff'),
                'location': 'Scaffold',
                'assemble': scaf_seq,
                'anno_cog': anno_cog,
                'specimen_id': self.option('sample_name'),
            }
            self.set_run(opts, self.circos, 'circos', self.step.circos)
        elif self.option('analysis') in ['complete']:

            seq_dir = self.seq_dir
            seq_type = self.get_seq_type(seq_dir)
            self.logger.info(seq_type)

            anno_cog =self.anno.output_dir + '/COG/' + self.option('sample_name') + '_cog_anno.xls'
            opts = {
                'gene':  self.pre_anno.option('gene_gff'),
                'rrna': self.pre_anno.option('rrna_gff'),
                'trna': self.pre_anno.option('trna_gff'),
                'location':  seq_type,
                'assemble': scaf_seq,
                'anno_cog': anno_cog,
                'specimen_id': self.option('sample_name'),
            }
            self.set_run(opts, self.circos, 'circos', self.step.circos)

    def run_pathogenic_system(self):
        opts = {
            'query':  self.pre_anno.option('faa'),
            'sample': self.option('sample_name'),
            "analysis_type" : self.option("analysis"),
            "analysis" : "card,vfdb,tcdb,phi,tmhmm"
        }
        self.set_run(opts, self.pathogenic_system, 'pathogenic_system', self.step.pathogenic_system)

    def run_antismash(self):
        if self.option('analysis') in ['uncomplete']:
            opts = {
                'genome_gbk': self.gbk_file.option('all_gbk'),
            }
            self.set_run(opts, self.antismash, 'antismash', self.step.antismash)
        elif self.option('analysis') in ['complete']:
            self.logger.info(self.gbk_file.option('gbk_dir').prop['path'])
            gbk =self.gbk_file.option('gbk_dir').prop['path']
            opts = {
                'gbk_dir': gbk,
                'sample_name': self.option('sample_name'),
                'analysis': self.option('analysis'),
            }
            self.logger.info(opts)
            self.set_run(opts, self.antismash_com, 'antismash', self.step.antismash)

    def run_overview(self):
        ## 区分完成图和扫描图
        if self.option('analysis') in ['complete']:
            tcdb = self.pathogenic_system.output_dir + '/TCDB/{}_whole_genome_tcdb_anno.xls'.format(self.option('sample_name'))
            tmhmm = self.pathogenic_system.output_dir + '/TMHMM/{}_whole_genome_tmhmm_anno.xls'.format(self.option('sample_name'))
        else:
            tcdb = self.pathogenic_system.output_dir + '/TCDB/{}_tcdb_anno.xls'.format(self.option('sample_name'))
            tmhmm = self.pathogenic_system.output_dir + '/TMHMM/{}_tmhmm_anno.xls'.format(self.option('sample_name'))


        sum = self.anno.output_dir + '/Summary/{}_anno_summary.xls'.format(self.option('sample_name'))
        phi = self.pathogenic_system.output_dir + '/PHI/{}_phi_anno.xls'.format(self.option('sample_name'))
        card = self.pathogenic_system.output_dir + '/CARD/{}_card_anno.xls'.format(self.option('sample_name'))

        vfdb = self.pathogenic_system.output_dir + '/VFDB/{}_vfdb_anno.xls'.format(self.option('sample_name'))
        cazy = self.anno.output_dir + '/CAZy/{}_anno_cazy.xls'.format(self.option('sample_name'))
        promote = self.promote.output_dir + '/{}_promoter_result.xls'.format(self.option('sample_name'))
        island = self.island.output_dir + '/Genomic_Islands/{}_GI_detail.xls'.format(self.option('sample_name'))
        prephage = self.prephage.output_dir + '/prephage/{}_prephage_detail.xls'.format(self.option('sample_name'))
        two_component = self.anno.output_dir + '/Two_component/{}.senser_regulator.xls'.format(self.option('sample_name'))
        if self.option("analysis") in ["complete"]:
            antismash = self.antismash_com.output_dir + '/{}_gene_antismash.xls'.format(self.option('sample_name'))
        else:
            antismash = self.antismash.output_dir + '/gene_antismash.xls'
        signal_n = self.pathogenic_system.output_dir + '/SIGNALP/{}_Gram-_SignalP.txt'.format(self.option('sample_name'))
        signal_p =  self.pathogenic_system.output_dir + '/SIGNALP/{}_Gram+_SignalP.txt'.format(self.option('sample_name'))
        gene_ori_info = self.pre_anno.work_dir + '/gene_desc.txt'
        files_list = ';'.join([sum,phi,card,tcdb,tmhmm,vfdb,cazy,promote,island,prephage,two_component,antismash,signal_p,signal_n,gene_ori_info])   ##原流程没有 gene_ori_info
        files_label = "sum;phi;card;tcdb;tmhmm;vfdb;cazy;promoter;island;prephage;two_sys;antismash;signal_p;signal_n;gene_ori_info"   ##原流程没有 gene_ori_info
        opts = {
            'files_list': files_list,
            'files_label': files_label ,
            'merge_column': 'Gene ID',
            'project_name': 'bacgenome'
        }
        self.set_run(opts, self.overview, 'overview',self.step.overview)

    def set_output(self,event):
        """
        将各个模块的结果输出至output
        """
        if not os.path.exists(self.output_dir + '/' + self.option("sample_name")):
            os.mkdir(self.output_dir + '/' + self.option("sample_name"))

        if event['data'] == 'overview':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/'+ self.option('sample_name') + '_summary.xls'):
                os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option('sample_name') + '_summary.xls')
            os.link(self.overview.output_dir + '/summary.xls', self.output_dir + '/' + self.option("sample_name") + '/' + self.option('sample_name') + '_summary.xls')

        if event['data'] == 'pre_anno':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary'):
                os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
            os.link(self.pre_anno.output_dir + '/' + self.option("sample_name") + '.project.summary',
                    self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')

            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.cds.ffn'):
                os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.cds.ffn')
            os.link(self.pre_anno.output_dir + '/' + self.option("sample_name") + '.cds.ffn',
                    self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.cds.ffn')
            self.move_dir(self.pre_anno.output_dir, self.option("sample_name")+'/gene')  #20190620


        if event['data'] == 'annotation':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/annotation'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/annotation')
            shutil.copytree(self.anno.output_dir, self.output_dir + '/' + self.option("sample_name") + '/annotation')
        if event['data'] == 'tree':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/tree'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/tree')
            shutil.copytree(self.tree.output_dir,self.output_dir + '/' + self.option("sample_name") + '/tree')
        if event['data'] == 'promote':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
            shutil.copytree(self.promote.output_dir, self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
        if event['data'] == 'gbk':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/gbk'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/gbk')
            shutil.copytree(self.gbk_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/gbk')
        if event['data'] == 'cgview':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/cgview'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/cgview')
            shutil.copytree(self.cgview.output_dir, self.output_dir + '/' + self.option("sample_name") + '/cgview')
        if event['data'] == 'circos':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/circos'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/circos')
            shutil.copytree(self.circos.output_dir,self.output_dir + '/' + self.option("sample_name") + '/circos')
        if event['data'] == 'antismash':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
            if self.option('analysis') in ['uncomplete']:
                shutil.copytree(self.antismash.output_dir, self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
            else:
                shutil.copytree(self.antismash_com.output_dir, self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
        if event['data'] == 'island':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands')
            shutil.copytree(self.island.output_dir + '/Genomic_Islands', self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands')
        if event['data'] == 'crispr':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
            shutil.copytree(self.crispr.output_dir + '/CRISPR_Cas', self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
        if event['data'] == 'prephage':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
            shutil.copytree(self.prephage.output_dir + '/prephage' , self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
        if event['data'] == 'pathogenic_system':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
            shutil.copytree(self.pathogenic_system.output_dir, self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
        if event['data'] == 'is_predict':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Is_Predict'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Is_Predict')
            shutil.copytree(self.is_predict.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Is_Predict')
        if event['data'] == 'integron':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Integron'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Integron')
            shutil.copytree(self.integron.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Integron')
        if event['data'] == 'repeatmasker':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat')
            shutil.copytree(self.repeatmasker.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat')
            if os.path.exists(self.fna):
                if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                    os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                os.link(self.fna, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
        if event['data'] == 'real_gbk':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
            shutil.copytree(self.gbk_real_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/real_gbk')

    def end(self):
        self.logger.info('BacgenomeAnno end writed by hand')
        super(BacgenomeAnnoModule, self).end()

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
            self.set_error('需要移动到output目录的文件夹不存在。', code="21400901")
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



    def get_seq_type(self,file_dir):
        list =[]
        files =os.listdir(file_dir)
        for file in files:
            lst=file.split('.')
            list.append('.'.join(lst[0:-1]))
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
        elif self.option("analysis") in ["complete"]:
            if not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                ass_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                output3 = self.work_dir + "/" + 'specimen.data'
                output4 = self.work_dir + "/" + 'gene.data'
                with open(ass_list_path, "rb") as l, open(output3, 'w') as file3, open(output4, 'w') as file4:
                    lines = l.readlines()
                    size = {}
                    pla_num = 0
                    file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                    file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                    for line in lines[1:]:
                        line2 = line.strip('\r\n').split()
                        if len(line2) == 3:
                            if re.search(r'plasmid',line2[2]):
                                pla_num +=1
                            file =os.path.join(self.option("asse_dir").prop['path'],line2[1])
                            num =os.path.getsize(file)
                            des = line2[1] + '\t' +  line2[2]
                            size[des]= num
                    dict_b = sorted(size.items(), key=lambda size: size[1], reverse=True)
                    num = 0
                    for b in dict_b:
                        b = b[0].split('\t')
                        if b[1] == 'chromosome':
                            file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'gene' + '\n')
                        if pla_num == 1:
                            if b[1] == 'plasmid':
                                file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'p_gene' + '\n')
                        elif pla_num >1:
                            if b[1] == 'plasmid':
                                num += 1
                                if num == 1:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pA_gene' + '\n')
                                elif num == 2:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pB_gene' + '\n')
                                elif num == 3:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pC_gene' + '\n')
                                elif num == 4:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pD_gene' + '\n')
                                elif num == 5:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pE_gene' + '\n')
                                elif num == 6:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pF_gene' + '\n')
                                elif num == 7:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pG_gene' + '\n')
                                elif num == 8:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pH_gene' + '\n')
                                elif num == 9:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pI_gene' + '\n')
                                elif num == 10:
                                    file4.write(self.option('sample_name') + '\t' + b[0] + '\t' + b[1] + '\t' + 'pJ_gene' + '\n')
                    file3.write(self.option('sample_name') + '\t' + '-' + '\t' + '-' + '\n')
            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                ass_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                output3 = self.work_dir + "/" + 'specimen.data'
                output4 = self.work_dir + "/" + 'gene.data'
                with open(raw_list_path, "rb") as l,open(ass_list_path, "rb") as f, open(output3, 'w') as file3, open(output4, 'w') as file4:
                    lib_list = []
                    raw_list = []
                    lines = f.readlines()
                    size = {}
                    pla_num = 0
                    file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                    file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                    for line in lines[1:]:
                        line2 = line.strip('\r\n').split()
                        if len(line2) == 3:
                            if re.search(r'plasmid', line2[2]):
                                pla_num += 1
                            file = os.path.join(self.option("asse_dir").prop['path'], line2[1])
                            num = os.path.getsize(file)
                            des = line2[1] + '\t' + line2[2]
                            size[des] = num
                    dict_b = sorted(size.iteritems(), key=lambda size: size[1], reverse=True)
                    num = 0
                    for key in dict_b:
                        key = key[0].split('\t')
                        if key[1] == 'chromosome':
                            file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'gene' + '\n')
                        if pla_num == 1:
                            if key[1] == 'plasmid':
                                file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'p_gene' + '\n')
                        elif pla_num > 1:
                            if key[1] == 'plasmid':
                                num += 1
                                if num == 1:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pA_gene' + '\n')
                                elif num == 2:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pB_gene' + '\n')
                                elif num == 3:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pC_gene' + '\n')
                                elif num == 4:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pD_gene' + '\n')
                                elif num == 5:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pE_gene' + '\n')
                                elif num == 6:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pF_gene' + '\n')
                                elif num == 7:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pG_gene' + '\n')
                                elif num == 8:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pH_gene' + '\n')
                                elif num == 9:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pI_gene' + '\n')
                                elif num == 10:
                                    file4.write(self.option('sample_name') + '\t' + key[0] + '\t' + key[1] + '\t' + 'pJ_gene' + '\n')
                    lines2 = l.readlines()
                    for line in lines2[1:]:
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
                    file3.write(self.option('sample_name') + '\t' + raw_file + '\t' + lib + '\n')