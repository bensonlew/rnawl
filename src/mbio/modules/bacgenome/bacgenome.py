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
from Bio import SeqIO
import shutil
from mbio.files.bacgenome.methy_dir import MethyDirFile
from mbio.packages.metagbin.common_function import link_dir,link_file


class BacgenomeModule(Module):
    """
    微生物基因组单个样品的工作流
    """
    def __init__(self, work_id):
        super(BacgenomeModule, self).__init__(work_id)
        options = [
            {"name": "raw_dir", "type": "infile", "format": "bacgenome.raw_dir"},###rawdata的文件目录
            {"name": "asse_dir", "type": "infile", "format": "sequence.fasta_dir"},  ###assemble的文件目录
            {"name": "sample_name", "type": "string"},###样品名称
            {"name": "analysis", "type": "string", "default": "uncomplete"}, ###流程分析模式complete，uncomplete
            {"name": "software_list", "type":"string", "default":""} , #使用的软件，逗号分割}
            {"name": "trans_code", "type" : "string","default":"11"},
            {"name": "nr_evalue","type":"string","default":"1e-5"},
            {"name": "swissprot_evalue","type":"string","default":"1e-5"},
            {"name": "cog_evalue","type":"string","default":"1e-5"},
            {"name": "kegg_evalue","type":"string","default":"1e-5"},
            {"name": "go_evalue","type": "string", "default":"1e-5"},
            {"name": "p_trans_code","type": "string","default":"11"}, #
            {"name": "p_software_list", "type": "string","default":"genemark"},  #质粒的基因预测使用的软件，逗号分割
            {'name': 'qc_tool', 'type': 'string', 'default': 'fastp'},
            {"name": "coverage", "type": "float", "default": 0.6},# 最小的coverage
            {"name": "identity", "type": "float", "default": 0.8},## 最小的identity
            {'name': 'assemble_tool', 'type': 'string'},  # 组装软件选择
            {'name': 'kmer', 'type': 'string', 'default': '21-47'},  # 质控软件选择
        ]
        self.add_option(options)
        self.sequence = self.add_module('sequence.bac_genome') #将数据解压、合并统计
        self.bac_qc = self.add_module('bacgenome.bacgenome_qc')  #将二代/三代数据质控
        self.assemble = self.add_module('bacgenome.bac_assemble')  # 二代数据组装
        self.complete_assemble = self.add_module('bacgenome.unicycler2')  # 二代数据组装
        self.assemble_assess_un = self.add_module('bacgenome.assemble_assess')  # 二代组装数据评估
        self.assemble_assess =self.add_tool('bacgenome.complete_asse_stat') #完成图组装数据评估
        self.genome_assess =self.add_module('bacgenome.genome_assess')  #基因组评估
        self.predict_gene =self.add_module('bacgenome.dnabac_predict') #基因预测，非编码预测，重复序列预测

        self.anno =self.add_module('bacgenome.annotation') #六种基因数据库注释
        self.gbk_file =self.add_tool('bacgenome.bac_gbk')  #gbk文件生成
        self.promote = self.add_tool('gene_structure.promote')#启动子生成
        self.island = self.add_module('bacgenome.island') #基因组岛预测
        self.crispr = self.add_tool('bacgenome.crisprcas') #CRISPR的预测
        #self.prephage =self.add_tool('bacgenome.prephage') #prephage的预测
        self.prephage = self.add_tool('bacgenome.phigaro')  # prephage的预测
        self.cgview = self.add_tool('bacgenome.cgview') #cgview画图
        self.tree = self.add_module('bacgenome.bac_tree') #进化树
        self.summary =self.add_tool('bacgenome.bac_summary') #项目总览统计
        self.circos = self.add_module('bacgenome.combine_circos')#circos画图
        self.pathogenic_system = self.add_module('annotation.pathogenic_system') #致病系统注释
        self.antismash = self.add_tool('annotation.antismash_v5')
        self.antismash_com = self.add_module('bacgenome.antismash')
        self.overview = self.add_tool('bacgenome.anno_summary')  #zouguanqing
        self.is_predict = self.add_module("bacgenome.is_predict")  #qingchen.zhang ## is可移动元件预测
        self.integron = self.add_module("bacgenome.integron_predict")  #qingchen.zhang ## 整合子可移动元件预测
        self.repeatmasker = self.add_module('bacgenome.repeatmasker')#qingchen.zhang ## 散在重复序列预测
        self.gbk_real_file =self.add_tool('bacgenome.bac_real_gbk')  #gbk文件生成
        self.resfinder_predict = self.add_tool('bacgenome.resfinder') #耐药基因resfinder预测
        self.srna_predict = self.add_tool('bacgenome.srna_predict')  # srna预测
        self.draf_pla = self.add_tool('bacgenome.draft_plasmid_predict')  # 扫描图质粒预测
        self.com_pla = self.add_tool('bacgenome.complete_plasmid_predict')  # 完成图图质粒预测
        self.busco = self.add_tool('bacgenome.busco')  #完成图时用于质量评估

        self.step.add_steps('busco','srna_predict', 'sequence', 'bac_qc','assemble','assemble_assess','genome_assess','predict_gene','tree','annotation','promote',
                            'gbk','cgview','circos','antismash','island','crispr','prephage','pathogenic_system','summary','overview', 'methylation', 'is_predict', 'integron', "repeatmasker", "real_gbk", "resfinder","plasmid")
        self.seq_type = ''

        self.annoo=[self.pathogenic_system,self.anno,self.genome_assess,self.tree, self.srna_predict]
        self.complete_anno =[self.pathogenic_system,self.anno,self.tree,self.srna_predict]  #完成图
        self.annoo2 = [self.pathogenic_system,self.anno,self.genome_assess,self.tree, self.srna_predict]
        self.complete_anno2 =[self.pathogenic_system,self.anno,self.tree, self.srna_predict]

        self.ana =[self.circos,self.prephage,self.crispr,self.island,self.cgview,self.promote,self.summary, self.is_predict, self.integron, self.repeatmasker, self.gbk_real_file, self.resfinder_predict]
        self.ana2 = [self.circos, self.prephage, self.crispr, self.island, self.cgview, self.promote, self.summary,self.antismash, self.is_predict, self.integron, self.repeatmasker, self.gbk_real_file, self.resfinder_predict]

    def check_options(self):
        """
        检查参数
        """
        if not self.option("analysis"):
            raise OptionError("必须输入流程运行类型！", code="21400901")
        if self.option("analysis") not in ['complete','uncomplete']:
            raise OptionError("请提分析流程参数！", code="21400902")
        else:
            if self.option("analysis") in ['uncomplete']:
                if not self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                    raise OptionError("请提供raw数据或assemble数据！不能同时为空！", code="21400903")
            if self.option("analysis") in ['complete']:
                if not self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                    raise OptionError("请提供raw数据或assemble数据！不能同时为空！", code="21400903")
        if  not self.option("sample_name"):
            raise OptionError("必须输入样本名！", code="21400905")

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
        super(BacgenomeModule, self).run()
        if self.option("raw_dir").is_set:
            raw_path = self.option("raw_dir").prop['path'] + '/' + 'list.txt'
            self.seq_type = self.get_sequence_type(raw_path)
            md = MethyDirFile()  ## add by qingchen.zhang 11 lines 完成图 上传原始序列时才有此逻辑
            md.set_path(self.option('raw_dir').prop['path'])
            md.get_info()
            if self.option("sample_name") in md.prop['bam_dict'] and md.prop['bam_dict'][self.option("sample_name")] != {}: ## fix by # qingchen.zhang
                self.methylation = self.add_tool("bacgenome.motifmaker")  # motifmaker甲基化
                self.annoo2.append(self.methylation)
                self.complete_anno2.append(self.methylation)
            elif self.option("sample_name") in md.prop["filedict"] and md.prop["filedict"][self.option("sample_name")] != {}:
                self.methylation = self.add_module("bacgenome.methylation")  #zouguanqing
                self.annoo2.append(self.methylation)
                self.complete_anno2.append(self.methylation)
        self.logger.info(self.seq_type)
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                self.overview.on('end',self.end)
                self.on_rely(self.ana2,self.run_overview)  # 总览
                self.gbk_file.on("end", self.run_resfinder)
                self.gbk_file.on("end", self.run_real_gbk)
                self.gbk_file.on("end", self.run_repeatmasker)
                self.gbk_file.on("end", self.run_is_predict)
                self.gbk_file.on("end", self.run_integron)
                self.gbk_file.on("end", self.run_antismash)
                self.gbk_file.on("end", self.run_circos)
                self.gbk_file.on("end", self.run_prephage)
                self.gbk_file.on("end", self.run_island)
                self.gbk_file.on("end", self.run_cgview)
                self.gbk_file.on("end", self.run_crispr)
                self.gbk_file.on("end", self.run_promote)
                self.gbk_file.on("end",self.run_summary)
                self.on_rely(self.annoo,self.run_gbk)
                self.predict_gene.on("end", self.run_annotation)
                self.predict_gene.on("end", self.run_pathogenic_system)
                self.predict_gene.on("end", self.run_genome_assess)
                self.predict_gene.on("end", self.run_tree)
                self.predict_gene.on("end", self.run_srna)
                self.draf_pla.on("end", self.run_predict)
                self.assemble_assess_un.on("end", self.run_plasmid)
                self.assemble.on("end", self.run_assess_uncomplete)
                self.bac_qc.on("end", self.run_assemble)
                self.sequence.on("end", self.run_bac_qc)
                self.run_sequence()
                
            elif not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                self.overview.on('end',self.end)
                self.on_rely(self.ana2,self.run_overview)
                self.gbk_file.on("end", self.run_resfinder)
                self.gbk_file.on("end", self.run_real_gbk)
                self.gbk_file.on("end", self.run_repeatmasker)
                self.gbk_file.on("end", self.run_is_predict)
                self.gbk_file.on("end", self.run_integron)
                self.gbk_file.on("end", self.run_antismash)
                self.gbk_file.on("end", self.run_circos)
                self.gbk_file.on("end", self.run_prephage)
                self.gbk_file.on("end", self.run_island)
                self.gbk_file.on("end", self.run_cgview)
                self.gbk_file.on("end", self.run_crispr)
                self.gbk_file.on("end", self.run_promote)
                self.gbk_file.on("end", self.run_summary)
                self.on_rely(self.complete_anno, self.run_gbk)
                self.predict_gene.on("end", self.run_annotation)
                self.predict_gene.on("end", self.run_pathogenic_system)
                self.predict_gene.on("end", self.run_tree)
                self.predict_gene.on("end", self.run_srna)
                self.draf_pla.on("end", self.run_predict)
                self.assemble_assess_un.on("end", self.run_plasmid)
                self.run_assess_uncomplete()
            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                self.overview.on('end',self.end)
                self.on_rely(self.ana2,self.run_overview)
                self.gbk_file.on("end", self.run_resfinder)
                self.gbk_file.on("end", self.run_real_gbk)
                self.gbk_file.on("end", self.run_repeatmasker)
                self.gbk_file.on("end", self.run_is_predict)
                self.gbk_file.on("end", self.run_integron)
                self.gbk_file.on("end", self.run_antismash)
                self.gbk_file.on("end", self.run_circos)
                self.gbk_file.on("end", self.run_crispr)
                self.gbk_file.on("end", self.run_prephage)
                self.gbk_file.on("end", self.run_island)
                self.gbk_file.on("end", self.run_cgview)
                self.gbk_file.on("end", self.run_promote)
                self.gbk_file.on("end", self.run_summary)
                self.on_rely(self.annoo,self.run_gbk)
                self.predict_gene.on("end", self.run_annotation)
                self.predict_gene.on("end", self.run_pathogenic_system)
                self.predict_gene.on("end", self.run_genome_assess)
                self.predict_gene.on("end", self.run_srna)
                self.predict_gene.on("end", self.run_tree)
                self.draf_pla.on("end", self.run_predict)
                self.assemble_assess_un.on("end", self.run_plasmid)
                self.bac_qc.on("end", self.run_assess_uncomplete)
                self.sequence.on("end", self.run_bac_qc)
                self.run_sequence()

        elif self.option("analysis") in ["complete"]:
            if not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                self.overview.on("end",self.end)
                self.antismash_com.on("end", self.run_overview)
                self.on_rely(self.ana, self.run_antismash)
                self.gbk_file.on("end", self.run_resfinder)
                self.gbk_file.on("end", self.run_real_gbk)
                self.gbk_file.on("end", self.run_repeatmasker)
                self.gbk_file.on("end", self.run_is_predict)
                self.gbk_file.on("end", self.run_integron)
                self.gbk_file.on("end", self.run_circos)
                self.gbk_file.on("end", self.run_prephage)
                self.gbk_file.on("end", self.run_island)
                self.gbk_file.on("end", self.run_cgview)
                self.gbk_file.on("end", self.run_crispr)
                self.gbk_file.on("end", self.run_promote)
                self.gbk_file.on("end", self.run_summary)
                self.on_rely(self.complete_anno,self.run_gbk)
                self.predict_gene.on("end", self.run_annotation)
                self.predict_gene.on("end", self.run_pathogenic_system)
                self.predict_gene.on("end", self.run_tree)
                self.predict_gene.on("end", self.run_srna)
                self.busco.on("end", self.run_predict)
                self.assemble_assess.on("end", self.run_buso)
                self.com_pla.on("end", self.run_assess_complete)
                self.run_plasmid()

            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                if re.search(r'PE',self.seq_type):
                    self.overview.on("end",self.end)
                    self.antismash_com.on("end", self.run_overview)
                    self.on_rely(self.ana, self.run_antismash)
                    self.gbk_file.on("end", self.run_resfinder)
                    self.gbk_file.on("end", self.run_real_gbk)
                    self.gbk_file.on("end", self.run_repeatmasker)
                    self.gbk_file.on("end", self.run_is_predict)
                    self.gbk_file.on("end", self.run_integron)
                    self.gbk_file.on("end", self.run_circos)
                    self.gbk_file.on("end", self.run_prephage)
                    self.gbk_file.on("end", self.run_island)
                    self.gbk_file.on("end", self.run_cgview)
                    self.gbk_file.on("end", self.run_crispr)
                    self.gbk_file.on("end", self.run_promote)
                    self.gbk_file.on("end", self.run_summary)
                    self.on_rely(self.annoo2, self.run_gbk)
                    self.predict_gene.on("end", self.run_annotation)
                    self.predict_gene.on("end", self.run_pathogenic_system)
                    self.predict_gene.on("end", self.run_genome_assess)
                    self.predict_gene.on("end", self.run_methylation)
                    self.predict_gene.on("end", self.run_tree)
                    self.predict_gene.on("end", self.run_srna)
                    self.busco.on("end", self.run_predict)
                    self.assemble_assess.on("end", self.run_buso)
                    self.com_pla.on("end", self.run_assess_complete)
                    self.bac_qc.on("end", self.run_plasmid)
                    self.sequence.on("end", self.run_bac_qc)
                    self.run_sequence()

            elif self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                if re.search(r'Pacbio', self.seq_type) and re.search(r'PE',self.seq_type):
                    self.overview.on("end", self.end)
                    self.antismash_com.on("end", self.run_overview)
                    self.on_rely(self.ana, self.run_antismash) # 总览
                    self.gbk_file.on("end", self.run_resfinder)
                    self.gbk_file.on("end", self.run_real_gbk)
                    self.gbk_file.on("end", self.run_repeatmasker)
                    self.gbk_file.on("end", self.run_is_predict)
                    self.gbk_file.on("end", self.run_integron)
                    self.gbk_file.on("end", self.run_circos)
                    self.gbk_file.on("end", self.run_prephage)
                    self.gbk_file.on("end", self.run_island)
                    self.gbk_file.on("end", self.run_cgview)
                    self.gbk_file.on("end", self.run_crispr)
                    self.gbk_file.on("end", self.run_promote)
                    self.gbk_file.on("end", self.run_summary)
                    self.on_rely(self.annoo2, self.run_gbk)
                    self.predict_gene.on("end", self.run_annotation)
                    self.predict_gene.on("end", self.run_pathogenic_system)
                    self.predict_gene.on("end", self.run_genome_assess)
                    self.predict_gene.on("end", self.run_methylation)
                    self.predict_gene.on("end", self.run_srna)
                    self.predict_gene.on("end", self.run_tree)
                    self.busco.on("end", self.run_predict)
                    self.assemble_assess.on("end", self.run_buso)
                    self.com_pla.on("end", self.run_assess_complete)
                    self.complete_assemble.on("end", self.run_plasmid)
                    self.bac_qc.on("end", self.run_assemble)
                    self.sequence.on("end", self.run_bac_qc)
                    self.run_sequence()
                elif re.search(r'Nanopore', self.seq_type) and re.search(r'PE',self.seq_type):
                    self.overview.on("end", self.end)
                    self.antismash_com.on("end", self.run_overview)
                    self.on_rely(self.ana, self.run_antismash)  # 总览
                    self.gbk_file.on("end", self.run_resfinder)
                    self.gbk_file.on("end", self.run_real_gbk)
                    self.gbk_file.on("end", self.run_repeatmasker)
                    self.gbk_file.on("end", self.run_is_predict)
                    self.gbk_file.on("end", self.run_integron)
                    self.gbk_file.on("end", self.run_circos)
                    self.gbk_file.on("end", self.run_prephage)
                    self.gbk_file.on("end", self.run_island)
                    self.gbk_file.on("end", self.run_cgview)
                    self.gbk_file.on("end", self.run_crispr)
                    self.gbk_file.on("end", self.run_promote)
                    self.gbk_file.on("end", self.run_summary)
                    self.on_rely(self.annoo, self.run_gbk)
                    self.predict_gene.on("end", self.run_annotation)
                    self.predict_gene.on("end", self.run_pathogenic_system)
                    self.predict_gene.on("end", self.run_genome_assess)
                    self.predict_gene.on("end", self.run_srna)
                    self.predict_gene.on("end", self.run_tree)
                    self.busco.on("end", self.run_predict)
                    self.assemble_assess.on("end", self.run_buso)
                    self.com_pla.on("end", self.run_assess_complete)
                    self.complete_assemble.on("end", self.run_plasmid)
                    self.bac_qc.on("end", self.run_assemble)
                    self.sequence.on("end", self.run_bac_qc)
                    self.run_sequence()
                else:
                    raise OptionError("完成图的rawdata必须是PE+Pacbio或PE+Nanopore！", code="21400903")


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

        files_list = [sum]
        files_label = ['sum']
        for f,l in zip([phi,card,tcdb,tmhmm,vfdb,cazy,promote],"phi;card;tcdb;tmhmm;vfdb;cazy;promoter".split(';')):
            if os.path.exists(f):
                files_list.append(f)
                files_label.append(l)
        files_list = ';'.join(files_list)
        files_label = ';'.join(files_label)

        if os.path.exists(prephage):
            files_list += ';'+prephage
            files_label += ';prephage'
        if os.path.exists(island):
            files_list += ';'+island
            files_label += ';island'
        if os.path.exists(antismash):
            files_list += ';'+antismash
            files_label += ';antismash'
        if os.path.exists(two_component):
            files_list += ';'+two_component
            files_label += ';two_sys'
        if os.path.exists(signal_p):
            files_list += ';'+signal_p
            files_label += ';signal_p'

        if os.path.exists(signal_n):
            files_list += ';'+signal_n
            files_label += ';signal_n'


        opts = {
            'files_list': files_list,
            'files_label' : files_label ,
            'merge_column' : 'Gene ID',
            'project_name' : 'bacgenome'
        }

        self.set_run(opts, self.overview, 'overview',self.step.overview)


    def run_sequence(self):
        opts = ''
        if self.option('raw_dir').is_set:
            opts = {
                'analysis': self.option('analysis'),
                'raw_dir': self.option('raw_dir'),
                'sample_name': self.option('sample_name'),
                'sequence_type': self.seq_type,
            }
        self.set_run(opts, self.sequence, 'sequence', self.step.sequence)

    def run_plasmid(self):
        if self.option('analysis') in ['uncomplete']:
            seq_path = ''
            if not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                seq_path = self.get_seq()
            elif self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                seq_path = self.assemble_assess_un.option('scaffold')
            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                seq_path = self.get_seq()
            opts = {
                'fasta': seq_path,
                'sample_name': self.option('sample_name'),
            }
            self.set_run(opts, self.draf_pla, 'plasmid', self.step.plasmid)
        elif self.option('analysis') in ['complete']:
            seq_path = ''
            status_path = ''
            if not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                seq_path = self.get_seq()
                self.get_status(seq_path, self.work_dir + "/all.status.txt")
                status_path = self.work_dir + "/all.status.txt"
            elif self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                seq_path = self.complete_assemble.output_dir +"/Unicycler/" + self.option('sample_name') +".scaffold.fna"
                status_path = self.complete_assemble.output_dir +"/Unicycler/all.status.txt"
            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                seq_path = self.get_seq()
                self.get_status(seq_path,self.work_dir+"/all.status.txt")
                status_path = self.work_dir+"/all.status.txt"
            opts = {
                'fasta': seq_path,
                'sample_name': self.option('sample_name'),
                'status': status_path,
            }
            self.set_run(opts, self.com_pla, 'plasmid', self.step.plasmid)


    def run_bac_qc(self):
        opts = ''
        if self.option('qc_tool') == 'no':
            opts = {
                'skip_module': 'T'
            }
        else:
            if self.sequence.option('dir').is_set and self.option('analysis') in ['uncomplete']:
                opts = {
                    'analysis': self.option('analysis'),
                    'fastq_dir': self.sequence.option('dir'),
                    'sample_name': self.option('sample_name'),
                    'qc_tool': self.option('qc_tool'),
                }
            elif self.sequence.option('pacbio_dir').is_set and self.sequence.option('dir').is_set and self.option('analysis') in ['complete']:
                opts = {
                    'analysis': self.option('analysis'),
                    'fastq_dir': self.sequence.option('dir'),
                    'sample_name': self.option('sample_name'),
                    'pacbio_dir': self.sequence.option('pacbio_dir'),
                    'qc_tool': self.option('qc_tool'),
                }
            elif self.sequence.option('nanopore_fq').is_set and self.sequence.option('dir').is_set and self.option('analysis') in ['complete']:
                opts = {
                    'analysis': self.option('analysis'),
                    'fastq_dir': self.sequence.option('dir'),
                    'sample_name': self.option('sample_name'),
                    'nanopore_fq': self.sequence.option('nanopore_fq'),
                    'qc_tool': self.option('qc_tool'),
                }
            elif (self.sequence.option('dir').is_set and self.option('analysis') in ['complete'] and not self.sequence.option('nanopore_fq').is_set) or (not self.sequence.option('pacbio_dir').is_set and self.sequence.option('dir').is_set and self.option('analysis') in ['complete']):
                opts = {
                    'analysis': self.option('analysis'),
                    'fastq_dir': self.sequence.option('dir'),
                    'sample_name': self.option('sample_name'),
                    'qc_tool': self.option('qc_tool'),
                }
            else:
                opts = {
                    'skip_module': 'T'
                }
        self.set_run(opts, self.bac_qc, 'bac_qc', self.step.bac_qc)

    def run_assemble(self):
        if self.option('analysis') in ['uncomplete'] and self.option('raw_dir').is_set:
            seq_type = self.get_assemble_type(self.option("raw_dir").prop['path'] + '/' + 'list.txt')
            self.logger.info(seq_type)
            clean_dir = self.bac_qc.output_dir + '/cleandata'
            sample_info = self.sequence.output_dir + '/data/' + 'sample_info'
            opts = {
                'fq_dir': clean_dir,
                'seq_type': seq_type,
                'sample_name': self.option('sample_name'),
                'sample_info': sample_info,
                'genome_size' : self.genome_size,
                'kmer':self.option("kmer")
            }
            self.set_run(opts, self.assemble, 'assemble', self.step.assemble)
        elif self.option('analysis') in ['complete'] and self.option('raw_dir').is_set:
            fq1 = ''
            fq2 = ''
            subread_fq = ''
            type = ''
            for i in os.listdir(self.bac_qc.output_dir + '/cleandata'):
                if re.search(".clean.1.fastq",i):
                    fq1 = self.bac_qc.output_dir + '/cleandata/' + i
                elif re.search(".clean.2.fastq",i):
                    fq2 = self.bac_qc.output_dir + '/cleandata/' + i
            if self.sequence.option('nanopore_fq').is_set:
                subread_fq = self.sequence.option('nanopore_fq').prop['path']
                type = 'nanopore'
            elif self.sequence.option('pacbio_dir').is_set:
                subread_fq = self.sequence.option("pacbio_dir").prop["path"] + '/all.pacbio.fq'
                type = 'pacbio'
            opts = {
                'read1': fq1,
                'read2': fq2,
                'sample_name': self.option('sample_name'),
                'data_type': type,
                'subread_fq' : subread_fq
            }
            self.set_run(opts, self.complete_assemble, 'assemble', self.step.assemble)

    def run_assess_uncomplete(self):
        if self.option('asse_dir').is_set :
            seq_path =self.get_seq()
        else:
            seq_path =self.assemble.option('scaffold')
        opts =''
        if self.option('raw_dir').is_set and self.option('asse_dir').is_set:
            opts = {
                'seq_scaf': seq_path,
                'sample_name': self.option('sample_name'),
                'type':"busco",
            }
        elif self.option('raw_dir').is_set and not self.option('asse_dir').is_set:
            opts = {
                'seq_scaf': seq_path,
                'sample_name': self.option('sample_name'),
                'type': "busco",
            }
        else:
            opts = {
                'seq_scaf': seq_path,
                'sample_name': self.option('sample_name'),
            }
        self.set_run(opts, self.assemble_assess_un, 'assemble_assess', self.step.assemble_assess)

    def run_assess_complete(self):
        self.seq_dir1 = self.get_split_fasta(self.com_pla.option('scf').prop['path'])
        fa = ''
        table = ''
        if self.option('analysis') in ['complete']:
            if self.option('raw_dir').is_set and self.option('asse_dir').is_set:
                fa = self.com_pla.option('scf').prop['path']
                table = self.com_pla.option('table').prop['path']
            elif not self.option('raw_dir').is_set and self.option('asse_dir').is_set:
                fa = self.com_pla.option('scf').prop['path']
                table = self.com_pla.option('table').prop['path']
            elif self.option('raw_dir').is_set and not self.option('asse_dir').is_set:
                fa = self.com_pla.option('scf').prop['path']
                table = self.com_pla.option('table').prop['path']
        opts = {
            'fa': fa,
            'table': table,
            'sample_name': self.option('sample_name'),
        }
        self.set_run(opts, self.assemble_assess, 'assemble_assess', self.step.assemble_assess)

    def run_buso(self):
        opts = {
            "scaf_fa": self.com_pla.option('scf').prop['path'],
            "database": "bacteria",
            "sample_name": self.option('sample_name'),
        }
        self.set_run(opts, self.busco, 'busco', self.step.busco)

    def run_genome_assess(self):
        bases = self.get_bases()
        opts = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
            opts = {
                'seq': scaf_seq,
                'fastq_dir': self.sequence.option('dir'),
                'bases': bases,
                'sample_name': self.option('sample_name'),
            }
        elif self.option('analysis') in ['complete']:
            scaf_seq =self.com_pla.option('scf')
            opts = {
                'seq': scaf_seq,
                'fastq_dir': self.sequence.option('dir'),
                'bases': bases,
                'sample_name': self.option('sample_name'),
            }
        self.set_run(opts, self.genome_assess, 'genome_assess', self.step.genome_assess)

    def run_srna(self):
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
            opts = {
                'fasta': scaf_seq,
                'sample_name': self.option('sample_name'),
            }
        elif self.option('analysis') in ['complete']:
            scaf_seq =self.com_pla.option('scf')
            opts = {
                'fasta': scaf_seq,
                'sample_name': self.option('sample_name'),
            }
        self.set_run(opts, self.srna_predict, 'srna_predict', self.step.srna_predict)

    def run_predict(self):
        opts = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq =self.assemble_assess_un.option('scaffold')
            opts = {
                'genome': scaf_seq,
                'sample':self.option('sample_name'),
                'software_list':self.option('software_list'),  #zouguanqing
                'trans_code' : self.option('trans_code')
            }
        if self.option('analysis') in ['complete']:
            chr_seq = self.assemble_assess.option('chr').prop['path']
            if self.assemble_assess.option('pla').is_set:
                pla_seq = self.assemble_assess.option('pla').prop['path']
                de =self.get_gene_prefix()
                self.logger.info(de)
                opts = {
                    'genome': chr_seq,
                    'sample': self.option('sample_name'),
                    'genome_plas':pla_seq,
                    'plasmid_prefix':de,
                    'software_list':self.option('software_list'),
                    'trans_code' : self.option('trans_code'),
                    'p_software_list' : self.option('p_software_list'),
                    'p_trans_code' :  self.option('p_trans_code')
                }
            else:
                opts = {
                    'genome': chr_seq,
                    'sample': self.option('sample_name'),
                    'software_list':self.option('software_list'),
                    'trans_code' : self.option('trans_code')
                }
        self.set_run(opts, self.predict_gene, 'predict_gene', self.step.predict_gene)

    def run_tree(self):
        if self.option('analysis') in ['uncomplete']:
            opts = {
                'seq_faa': self.predict_gene.option('sample_gene_faa'),
                'seq_fa': self.predict_gene.option('sample_gene_fnn'),
                'gene_gff': self.predict_gene.option('sample_gene_gff'),
                'sample_name': self.option('sample_name'),
                'analysis': self.option('analysis'),
            }
            self.set_run(opts, self.tree, 'tree', self.step.tree)
        elif self.option('analysis') in ['complete']:
            seq_faa = ''
            seq_fa = self.predict_gene.output_dir + '/predict/CDS_predict/' + self.option('sample_name') + '_whole_genome_CDS.fnn'
            if os.path.exists(self.predict_gene.output_dir + '/predict/CDS_predict/' + self.option('sample_name') + '_Chromosome_CDS.faa'):
                seq_faa =self.predict_gene.output_dir + '/predict/CDS_predict/' + self.option('sample_name') + '_Chromosome_CDS.faa'
            else:
                seq_faa = self.predict_gene.work_dir + '/chromosome_CDS.faa'
            scaf_seq = self.com_pla.option('scf')
            opts = {
                'seq_faa': seq_faa,
                'seq_fa': seq_fa,
                'sample_name': self.option('sample_name'),
                'analysis': self.option('analysis'),
                'gene_gff': self.predict_gene.option('sample_gene_gff'),
                'rrna_gff':self.predict_gene.option('sample_rrna_gff'),
                'genome_fa':scaf_seq,
            }
            self.set_run(opts, self.tree, 'tree', self.step.tree)

    def run_methylation(self):  #zouguanqing
        md = MethyDirFile()
        md.set_path(self.option('raw_dir').prop['path'])
        md.get_info()
        scaf_seq = self.com_pla.option('scf')
        if self.option("sample_name") in md.prop["bam_dict"] and md.prop["bam_dict"][self.option("sample_name")] != {}:## fix by qingchen.zhang
            #self.methylation = self.add_tool("bacgenome.motifmaker")  # motifmaker甲基化
            opts = {
                "input": md.prop["methy_files"][self.option("sample_name")],
                "ref_database": scaf_seq,
                "sample": self.option("sample_name")
            }
            self.set_run(opts, self.methylation, 'methylation', self.step.methylation)
        elif self.option("sample_name") in md.prop["filedict"] and md.prop["filedict"][self.option("sample_name")] != {}:
            #self.methylation = self.add_module("bacgenome.methylation")  #zouguanqing
            opts = {
                "methy_dir" : md,
                "ref_fa" :scaf_seq,
            }
            self.set_run(opts, self.methylation, 'methylation', self.step.methylation)

    def run_annotation(self):
        opts = {
            'gene_seq': self.predict_gene.option('sample_gene_faa'),
            'gene_gff': self.predict_gene.option('sample_gene_gff'),
            'database_list': 'nr_v20200604,swissprot_v20200617,pfam_v33.1,eggnog,kegg_v94.2,cazy_v8',
            'sample': self.option('sample_name'),
            'analysis': self.option('analysis'),
            "has_two_component" : 'T',
            'produce_mark' : True
        }
        opts['nr_evalue'] = self.option('nr_evalue')
        opts['cog_evalue'] = self.option('cog_evalue')
        opts['kegg_evalue'] = self.option('kegg_evalue')
        opts['swissprot_evalue'] = self.option('swissprot_evalue')

        self.set_run(opts, self.anno, 'annotation', self.step.annotation)

    def run_summary(self):
        gene_stat = ''
        asse = ''
        if self.option('analysis') in ['uncomplete']:
            gene_stat =self.predict_gene.output_dir + '/predict/CDS_predict/' + self.option('sample_name') + '_CDS_statistics.xls'
            asse = self.assemble_assess_un.output_dir + '/assembly/' + self.option('sample_name') + '_assembly_summary.xls'
        elif self.option('analysis') in ['complete']:
            gene_stat = self.predict_gene.output_dir + '/predict/CDS_predict/' + self.option('sample_name') + '_whole_genome_CDS_statistics.xls'
            asse = self.assemble_assess.output_dir + '/' + self.option('sample_name') + '_assembly_summary.xls'
        opts = {
            'gene_statistics': gene_stat,
            'rrna_gff': self.predict_gene.option('sample_rrna_gff'),
            'trna_gff': self.predict_gene.option('sample_trna_gff'),
            'assemble': asse,
            'cog': self.anno.output_dir + '/COG/' + self.option('sample_name') + '_cog_anno.xls',
            'kegg': self.anno.output_dir + '/KEGG/' + self.option('sample_name') + '_kegg_anno.xls',
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.summary, 'summary', self.step.summary)

    def run_gbk(self):
        self.logger.info(">>>in run_gbk")
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq =self.com_pla.option('scf')
        opts = {
            'gen_gff': self.predict_gene.option('sample_gene_gff'),
            'rrna_gff': self.predict_gene.option('sample_rrna_gff'),
            'trna_gff': self.predict_gene.option('sample_trna_gff'),
            'pro_fa': self.predict_gene.option('sample_gene_faa'),
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
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq =self.com_pla.option('scf')
        opts = {
            'gen_gff': self.predict_gene.option('sample_gene_gff'),
            'rrna_gff': self.predict_gene.option('sample_rrna_gff'),
            'trna_gff': self.predict_gene.option('sample_trna_gff'),
            'pro_fa': self.predict_gene.option('sample_gene_faa'),
            'genome_fa': scaf_seq,
            'anno': self.anno.option('tidy_summary'),
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.gbk_real_file, 'real_gbk', self.step.real_gbk)

    def run_promote(self):
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq = self.com_pla.option('scf')
        opts = {
            'sequence': self.predict_gene.option('sample_gene_fnn'),
            'assemble': scaf_seq,
            'sample': self.option('sample_name'),
        }
        self.set_run(opts, self.promote, 'promote', self.step.promote)

    def run_crispr(self):
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq = self.com_pla.option('scf')
        opts = {
            'seq': scaf_seq,
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.crispr, 'crispr', self.step.crispr)

    def run_prephage(self):
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq = self.com_pla.option('scf')
        opts = {
            'fasta': scaf_seq,
            'sample_name': self.option('sample_name'),
            'anno': self.anno.option('tidy_summary'),
        }
        self.set_run(opts, self.prephage, 'prephage', self.step.prephage)

    def run_island(self):
        seq_dir = ''
        if self.option('analysis') in ['uncomplete']:
            seq_dir =self.assemble_assess_un.output_dir + '/scaffold'
        elif self.option('analysis') in ['complete']:
            seq_dir = self.seq_dir1
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
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq =self.com_pla.option('scf')
        opts = {
            'genome_fa': scaf_seq,
            'gene_faa': self.predict_gene.option('sample_gene_faa'),
            'gene_fna': self.predict_gene.option('sample_gene_fnn'),
            'gene_gff': self.predict_gene.option('sample_gene_gff').prop['path'],
            'sample': self.option('sample_name')
        }
        self.set_run(opts, self.is_predict, 'is_predict', self.step.is_predict)

    def run_integron(self):
        """
        integron整合子预测
        :return:
        """
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq = self.com_pla.option('scf')
        opts = {
            'genome_fa': scaf_seq,
            'gene_faa': self.predict_gene.option('sample_gene_faa'),
            'gene_gff': self.predict_gene.option('sample_gene_gff').prop['path'],
            'sample': self.option('sample_name')
        }
        self.set_run(opts, self.integron, 'integron', self.step.integron)

    def run_repeatmasker(self):
        """
        散在重复预测
        :return:
        """
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            self.scaf_seq = self.assemble_assess_un.option('scaffold').prop['path']
        elif self.option('analysis') in ['complete']:
            self.scaf_seq = self.com_pla.option('scf')
        opts = {
            "genome_fa": self.scaf_seq,
            "sample": self.option('sample_name'),
            "analysis": self.option('analysis'),
        }
        self.set_run(opts, self.repeatmasker, 'repeatmasker', self.step.repeatmasker)

    def run_resfinder(self):
        """
        耐药基因resfinder预测
        :return:
        """
        options = {
            "gene_fa": self.predict_gene.option('sample_gene_fnn'),
            "gene_gff": self.predict_gene.option('sample_gene_gff'),
            "sample_name": self.option('sample_name'),
            "min_cov": self.option("coverage"),
            "min_iden": self.option("identity"),
            "resfinder_database": "true",
            "disinfinder_database": "true",
            "species_name": ""
            }
        self.set_run(options, self.resfinder_predict, 'resfinder', self.step.resfinder)

    def run_cgview(self):
        scaf_seq = ''
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
        elif self.option('analysis') in ['complete']:
            scaf_seq =self.com_pla.option('scf')
        opts = {
            'gen_gff':  self.predict_gene.option('sample_gene_gff'),
            'rrna_gff': self.predict_gene.option('sample_rrna_gff'),
            'trna_gff': self.predict_gene.option('sample_trna_gff'),
            'pro_fa':  self.predict_gene.option('sample_gene_faa'),
            'genome_fa': scaf_seq,
            'anno': self.anno.option('tidy_summary'),
            'sample_name': self.option('sample_name'),
            'analysis': self.option('analysis'),
        }
        self.set_run(opts, self.cgview, 'cgview', self.step.cgview)

    def run_circos(self):
        if self.option('analysis') in ['uncomplete']:
            scaf_seq = self.assemble_assess_un.option('scaffold')
            anno_cog =self.anno.output_dir + '/COG/' + self.option('sample_name') + '_cog_anno.xls'
            opts = {
                'gene':  self.predict_gene.option('sample_gene_gff'),
                'rrna': self.predict_gene.option('sample_rrna_gff'),
                'trna': self.predict_gene.option('sample_trna_gff'),
                'location': 'Scaffold',
                'assemble': scaf_seq,
                'anno_cog': anno_cog,
                'specimen_id': self.option('sample_name'),
            }
            self.set_run(opts, self.circos, 'circos', self.step.circos)
        elif self.option('analysis') in ['complete']:
            scaf_seq = self.com_pla.option('scf')
            seq_dir =self.seq_dir1
            seq_type = self.get_seq_type(seq_dir)
            self.logger.info(seq_type)
            anno_cog =self.anno.output_dir + '/COG/' + self.option('sample_name') + '_cog_anno.xls'
            opts = {
                'gene': self.predict_gene.option('sample_gene_gff'),
                'rrna': self.predict_gene.option('sample_rrna_gff'),
                'trna': self.predict_gene.option('sample_trna_gff'),
                'location': seq_type,
                'assemble': scaf_seq,
                'anno_cog': anno_cog,
                'specimen_id': self.option('sample_name'),
            }
            self.set_run(opts, self.circos, 'circos', self.step.circos)

    def run_pathogenic_system(self):
        opts = {
            'query':  self.predict_gene.option('sample_gene_faa'),
            'sample': self.option('sample_name'),
            "analysis_type" : self.option("analysis"),
            "analysis" : "card,vfdb,tcdb,phi,tmhmm"

        }
        self.set_run(opts, self.pathogenic_system, 'pathogenic_system', self.step.pathogenic_system)

    def run_antismash(self):
        if self.option('analysis') in ['uncomplete']:
            opts = {
                'genome_gbk': self.gbk_file.option('antismash_gbk'),
            }
            self.set_run(opts, self.antismash, 'antismash', self.step.antismash)
        elif self.option('analysis') in ['complete']:
            self.logger.info(self.gbk_file.option('antismash_gbk_co').prop['path'])
            gbk =self.get_antismash_gbk(self.gbk_file.option('antismash_gbk_co').prop['path']) ## modify by qingchen.zhang # 原因antismash5.0 升级后原来处理的gbk_dir文件不再适用
            opts = {
                'gbk_dir': gbk,
                'sample_name': self.option('sample_name'),
                'analysis': self.option('analysis'),
            }
            self.logger.info(opts)
            self.set_run(opts, self.antismash_com, 'antismash', self.step.antismash)
            self.logger.info('aaaaa')

    def set_output(self,event):
        """
        将各个模块的结果输出至output
        """
        obj = event['bind_object']
        if event['data'] == 'overview':
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/'+ self.option('sample_name') + '_summary.xls'):
                os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option('sample_name') + '_summary.xls')
            os.link(self.overview.output_dir + '/summary.xls', self.output_dir + '/' + self.option("sample_name") + '/' + self.option('sample_name') + '_summary.xls')
            ## 下面的seq和gff用于打通小工具使用，上传至s3nb，并且将初始化的基因组对象换成路径
            if self.option('analysis') in ['uncomplete']:
                self.scaf_seq = self.assemble_assess_un.option('scaffold').prop['path']
            elif self.option('analysis') in ['complete']:
                self.scaf_seq =self.com_pla.option('scf').prop['path']
            if os.path.exists(self.scaf_seq):
                shutil.copyfile(self.scaf_seq, self.output_dir + '/' + self.option("sample_name") + '/' +self.option("sample_name") +'.all.fna')
            if os.path.exists(self.gbk_real_file.output_dir + "/" + self.option("sample_name") +'.all.gff'):
                shutil.copyfile(self.gbk_real_file.output_dir + "/" + self.option("sample_name") +'.all.gff', self.output_dir + '/' + self.option("sample_name") + '/' +self.option("sample_name") +'.all.gff')

        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                if event['data'] == 'bac_qc':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/fastx'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/fastx')
                    shutil.copytree(self.bac_qc.output_dir + '/fastx',self.output_dir + '/' + self.option("sample_name") + '/fastx')
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/data_QC'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                    shutil.copytree(self.bac_qc.output_dir +'/data_QC' ,self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                if event['data'] == 'genome_assess':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                    shutil.copytree(self.genome_assess.output_dir, self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                if event['data'] == 'assemble_assess':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                    shutil.copytree(self.assemble_assess_un.output_dir, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                if event['data'] == 'predict_gene':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                    shutil.copytree(self.predict_gene.output_dir + '/predict', self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
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
                if event['data'] == 'summary':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary'):
                        os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                    os.link(self.summary.output_dir + '/' + self.option("sample_name") + '.project.summary',self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
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
                    shutil.copytree(self.antismash.output_dir, self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
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
                    shutil.copytree(self.prephage.output_dir , self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
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
                    self.scaf_seq = self.assemble_assess_un.option('scaffold').prop['path']
                    if os.path.exists(self.scaf_seq):
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                            os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                        os.link(self.scaf_seq, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                if event['data'] == 'real_gbk':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                    shutil.copytree(self.gbk_real_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                if event['data'] == 'resfinder':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                    shutil.copytree(self.resfinder_predict.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                if event['data'] == 'srna_predict':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/srna_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                    shutil.copytree(self.srna_predict.output_dir, self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                if event['data'] == 'plasmid':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                    shutil.copytree(self.draf_pla.output_dir, self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')

            if not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                if event['data'] == 'assemble_assess':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                    shutil.copytree(self.assemble_assess_un.output_dir, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                if event['data'] == 'predict_gene':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                    shutil.copytree(self.predict_gene.output_dir + '/predict', self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
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
                if event['data'] == 'summary':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary'):
                        os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                    os.link(self.summary.output_dir + '/' + self.option("sample_name") + '.project.summary',self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
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
                    shutil.copytree(self.antismash.output_dir, self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
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
                    shutil.copytree(self.prephage.output_dir, self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
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
                    self.scaf_seq = self.assemble_assess_un.option('scaffold').prop['path']
                    if os.path.exists(self.scaf_seq):
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                            os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                        os.link(self.scaf_seq, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                if event['data'] == 'real_gbk':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                    shutil.copytree(self.gbk_real_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                if event['data'] == 'resfinder':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                    shutil.copytree(self.resfinder_predict.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                if event['data'] == 'srna_predict':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/srna_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                    shutil.copytree(self.srna_predict.output_dir, self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                if event['data'] == 'plasmid':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                    shutil.copytree(self.draf_pla.output_dir, self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')

            if self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                if event['data'] == 'bac_qc':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/fastx'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/fastx')
                    shutil.copytree(self.bac_qc.output_dir + '/fastx',self.output_dir + '/' + self.option("sample_name") + '/fastx')
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/data_QC'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                    shutil.copytree(self.bac_qc.output_dir + '/data_QC',self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                if event['data'] == 'genome_assess':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                    shutil.copytree(self.genome_assess.output_dir, self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                if event['data'] == 'assemble_assess':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                    shutil.copytree(self.assemble_assess_un.output_dir, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                if event['data'] == 'predict_gene':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                    shutil.copytree(self.predict_gene.output_dir + '/predict', self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
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
                if event['data'] == 'summary':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary'):
                        os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                    os.link(self.summary.output_dir + '/' + self.option("sample_name") + '.project.summary',self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
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
                    shutil.copytree(self.antismash.output_dir, self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
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
                    shutil.copytree(self.prephage.output_dir , self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
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
                    self.scaf_seq = self.assemble_assess_un.option('scaffold').prop['path']
                    if os.path.exists(self.scaf_seq):
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                            os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                        os.link(self.scaf_seq, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                if event['data'] == 'real_gbk':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                    shutil.copytree(self.gbk_real_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                if event['data'] == 'resfinder':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                    shutil.copytree(self.resfinder_predict.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                if event['data'] == 'srna_predict':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/srna_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                    shutil.copytree(self.srna_predict.output_dir, self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                if event['data'] == 'plasmid':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                    shutil.copytree(self.draf_pla.output_dir, self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')

        elif self.option("analysis") in ["complete"]:
            if not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                if event['data'] == 'assemble_assess':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                    shutil.copytree(self.assemble_assess.output_dir, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                if event['data'] == 'busco':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly/busco'):
                        shutil.rmtree(
                                self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly/busco')
                    shutil.copytree(self.busco.output_dir, self.output_dir + '/' + self.option(
                            "sample_name") + '/assembly_predict/assembly/busco')
                if event['data'] == 'predict_gene':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                    shutil.copytree(self.predict_gene.output_dir + '/predict', self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
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
                if event['data'] == 'summary':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary'):
                        os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                    os.link(self.summary.output_dir + '/' + self.option("sample_name") + '.project.summary',self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
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
                    shutil.copytree(self.prephage.output_dir, self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
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
                    self.scaf_seq =self.assemble_assess.work_dir + '/all.fasta'
                    if os.path.exists(self.scaf_seq):
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                            os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                        os.link(self.scaf_seq, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                if event['data'] == 'real_gbk':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                    shutil.copytree(self.gbk_real_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                if event['data'] == 'resfinder':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                    shutil.copytree(self.resfinder_predict.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                if event['data'] == 'srna_predict':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/srna_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                    shutil.copytree(self.srna_predict.output_dir, self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                if event['data'] == 'plasmid':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                    shutil.copytree(self.com_pla.output_dir, self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')

            elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                if event['data'] == 'methylation':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/methylation'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/methylation')
                    if os.path.exists(self.methylation.output_dir + '/' + self.option("sample_name")):
                        shutil.copytree(self.methylation.output_dir + '/' + self.option("sample_name"),
                                        self.output_dir + '/' + self.option("sample_name") + '/methylation')
                    else:
                        shutil.copytree(self.methylation.output_dir,
                                        self.output_dir + '/' + self.option("sample_name") + '/methylation')

                if re.search(r'PE', self.seq_type):
                    if event['data'] == 'bac_qc':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/fastx'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/fastx')
                        shutil.copytree(self.bac_qc.output_dir + '/fastx',self.output_dir + '/' + self.option("sample_name") + '/fastx')
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/data_QC'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                        shutil.copytree(self.bac_qc.output_dir + '/data_QC',self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                    if event['data'] == 'assemble_assess':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                        shutil.copytree(self.assemble_assess.output_dir, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                    if event['data'] == 'busco':
                        if os.path.exists(
                                self.output_dir + '/' + self.option(
                                    "sample_name") + '/assembly_predict/assembly/busco'):
                            shutil.rmtree(
                                self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly/busco')
                        shutil.copytree(self.busco.output_dir, self.output_dir + '/' + self.option(
                            "sample_name") + '/assembly_predict/assembly/busco')
                    if event['data'] == 'genome_assess':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                        shutil.copytree(self.genome_assess.output_dir, self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                    if event['data'] == 'predict_gene':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                        shutil.copytree(self.predict_gene.output_dir + '/predict', self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                    if event['data'] == 'annotation':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/annotation'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/annotation')
                        shutil.copytree(self.anno.output_dir,self.output_dir + '/' + self.option("sample_name") + '/annotation')
                    if event['data'] == 'tree':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/tree'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/tree')
                        shutil.copytree(self.tree.output_dir,self.output_dir + '/' + self.option("sample_name") + '/tree')
                    if event['data'] == 'promote':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
                        shutil.copytree(self.promote.output_dir, self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
                    if event['data'] == 'summary':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary'):
                            os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                        os.link(self.summary.output_dir + '/' + self.option("sample_name") + '.project.summary',self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                    if event['data'] == 'gbk':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/gbk'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/gbk')
                        shutil.copytree(self.gbk_file.output_dir,self.output_dir + '/' + self.option("sample_name") + '/gbk')
                    if event['data'] == 'cgview':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/cgview'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/cgview')
                        shutil.copytree(self.cgview.output_dir,self.output_dir + '/' + self.option("sample_name") + '/cgview')
                    if event['data'] == 'circos':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/circos'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/circos')
                        shutil.copytree(self.circos.output_dir,self.output_dir + '/' + self.option("sample_name") + '/circos')
                    if event['data'] == 'antismash':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
                        shutil.copytree(self.antismash_com.output_dir, self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
                    if event['data'] == 'island':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands')
                        shutil.copytree(self.island.output_dir + '/Genomic_Islands',self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands')
                    if event['data'] == 'crispr':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
                        shutil.copytree(self.crispr.output_dir + '/CRISPR_Cas',self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
                    if event['data'] == 'prephage':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
                        shutil.copytree(self.prephage.output_dir, self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
                    if event['data'] == 'pathogenic_system':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
                        shutil.copytree(self.pathogenic_system.output_dir,self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
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
                        self.scaf_seq =self.com_pla.option('scf').prop['path']
                        if os.path.exists(self.scaf_seq):
                            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                                os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                            os.link(self.scaf_seq, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                    if event['data'] == 'real_gbk':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                        shutil.copytree(self.gbk_real_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                    if event['data'] == 'resfinder':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                        shutil.copytree(self.resfinder_predict.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                    if event['data'] == 'srna_predict':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/srna_predict'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                        shutil.copytree(self.srna_predict.output_dir,
                                            self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                    if event['data'] == 'plasmid':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                        shutil.copytree(self.com_pla.output_dir,
                                        self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                else:
                    if event['data'] == 'bac_qc':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/data_QC'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/data_QC'):
                            shutil.copytree(self.bac_qc.output_dir + '/data_QC', self.output_dir + '/' + self.option("sample_name") + '/data_QC')   #guanqing.zou 20180912
                    if event['data'] == 'assemble_assess':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                        shutil.copytree(self.assemble_assess.output_dir, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                    if event['data'] == 'busco':
                        if os.path.exists(
                                self.output_dir + '/' + self.option(
                                    "sample_name") + '/assembly_predict/assembly/busco'):
                            shutil.rmtree(
                                self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly/busco')
                        shutil.copytree(self.busco.output_dir, self.output_dir + '/' + self.option(
                            "sample_name") + '/assembly_predict/assembly/busco')
                    if event['data'] == 'predict_gene':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                        shutil.copytree(self.predict_gene.output_dir + '/predict', self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                    if event['data'] == 'annotation':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/annotation'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/annotation')
                        shutil.copytree(self.anno.output_dir,self.output_dir + '/' + self.option("sample_name") + '/annotation')
                    if event['data'] == 'tree':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/tree'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/tree')
                        shutil.copytree(self.tree.output_dir,self.output_dir + '/' + self.option("sample_name") + '/tree')
                    if event['data'] == 'promote':
                        if os.path.exists( self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict'):
                            shutil.rmtree( self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
                        shutil.copytree(self.promote.output_dir, self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
                    if event['data'] == 'summary':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary'):
                            os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                        os.link(self.summary.output_dir + '/' + self.option("sample_name") + '.project.summary',self.output_dir + '/' + self.option("sample_name") + '/' + self.option("sample_name") + '.project.summary')
                    if event['data'] == 'gbk':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/gbk'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/gbk')
                        shutil.copytree(self.gbk_file.output_dir,self.output_dir + '/' + self.option("sample_name") + '/gbk')
                    if event['data'] == 'cgview':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/cgview'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/cgview')
                        shutil.copytree(self.cgview.output_dir,self.output_dir + '/' + self.option("sample_name") + '/cgview')
                    if event['data'] == 'circos':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/circos'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/circos')
                        shutil.copytree(self.circos.output_dir,self.output_dir + '/' + self.option("sample_name") + '/circos')
                    if event['data'] == 'antismash':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
                        shutil.copytree(self.antismash_com.output_dir, self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
                    if event['data'] == 'island':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands')
                        shutil.copytree(self.island.output_dir + '/Genomic_Islands',self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands')
                    if event['data'] == 'crispr':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
                        shutil.copytree(self.crispr.output_dir + '/CRISPR_Cas',self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
                    if event['data'] == 'prephage':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
                        shutil.copytree(self.prephage.output_dir, self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
                    if event['data'] == 'pathogenic_system':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
                        shutil.copytree(self.pathogenic_system.output_dir,self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
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
                        self.scaf_seq =self.com_pla.option('scf').prop['path']
                        if os.path.exists(self.scaf_seq):
                            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                                os.remove(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                            os.link(self.scaf_seq, self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                    if event['data'] == 'real_gbk':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                        shutil.copytree(self.gbk_real_file.output_dir, self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                    if event['data'] == 'resfinder':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                        shutil.copytree(self.resfinder_predict.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                    if event['data'] == 'srna_predict':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/srna_predict'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                        shutil.copytree(self.srna_predict.output_dir,
                                            self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                    if event['data'] == 'plasmid':
                        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict'):
                            shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                        shutil.copytree(self.com_pla.output_dir,
                                        self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
            elif self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                if event['data'] == 'methylation':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/methylation'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/methylation')
                    if os.path.exists(self.methylation.output_dir + '/' + self.option("sample_name")):
                        shutil.copytree(self.methylation.output_dir + '/' + self.option("sample_name"), self.output_dir + '/' + self.option("sample_name") + '/methylation')
                    else:
                        shutil.copytree(self.methylation.output_dir,  self.output_dir + '/' + self.option("sample_name") + '/methylation')
                if event['data'] == 'bac_qc':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/fastx'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/fastx')
                    shutil.copytree(self.bac_qc.output_dir + '/fastx',
                                    self.output_dir + '/' + self.option("sample_name") + '/fastx')
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/data_QC'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                    shutil.copytree(self.bac_qc.output_dir + '/data_QC',
                                    self.output_dir + '/' + self.option("sample_name") + '/data_QC')
                if event['data'] == 'assemble':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/Unicycler'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/Unicycler')
                    shutil.copytree(self.complete_assemble.output_dir + '/Unicycler', self.output_dir + '/' + self.option("sample_name") + '/Unicycler')
                if event['data'] == 'assemble_assess':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                    shutil.copytree(self.assemble_assess.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly')
                if event['data'] == 'busco':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly/busco'):
                        shutil.rmtree(
                                self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/assembly/busco')
                    shutil.copytree(self.busco.output_dir, self.output_dir + '/' + self.option(
                            "sample_name") + '/assembly_predict/assembly/busco')
                if event['data'] == 'genome_assess':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                    shutil.copytree(self.genome_assess.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/genomic_assessment')
                if event['data'] == 'predict_gene':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                    shutil.copytree(self.predict_gene.output_dir + '/predict',
                                    self.output_dir + '/' + self.option("sample_name") + '/assembly_predict/predict')
                if event['data'] == 'annotation':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/annotation'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/annotation')
                    shutil.copytree(self.anno.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/annotation')
                if event['data'] == 'tree':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/tree'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/tree')
                    shutil.copytree(self.tree.output_dir, self.output_dir + '/' + self.option("sample_name") + '/tree')
                if event['data'] == 'promote':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict'):
                        shutil.rmtree(
                            self.output_dir + '/' + self.option("sample_name") + '/structral_genome/promoter_predict')
                    shutil.copytree(self.promote.output_dir, self.output_dir + '/' + self.option(
                        "sample_name") + '/structral_genome/promoter_predict')
                if event['data'] == 'summary':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/' + self.option(
                            "sample_name") + '.project.summary'):
                        os.remove(self.output_dir + '/' + self.option("sample_name") + '/' + self.option(
                            "sample_name") + '.project.summary')
                    os.link(self.summary.output_dir + '/' + self.option("sample_name") + '.project.summary',
                            self.output_dir + '/' + self.option("sample_name") + '/' + self.option(
                                "sample_name") + '.project.summary')
                if event['data'] == 'gbk':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/gbk'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/gbk')
                    shutil.copytree(self.gbk_file.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/gbk')
                if event['data'] == 'cgview':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/cgview'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/cgview')
                    shutil.copytree(self.cgview.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/cgview')
                if event['data'] == 'circos':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/circos'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/circos')
                    shutil.copytree(self.circos.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/circos')
                if event['data'] == 'antismash':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH'):
                        shutil.rmtree(
                            self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
                    shutil.copytree(self.antismash_com.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/metabolic_system/antiSMASH')
                if event['data'] == 'island':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands'):
                        shutil.rmtree(
                            self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Genomic_Islands')
                    shutil.copytree(self.island.output_dir + '/Genomic_Islands', self.output_dir + '/' + self.option(
                        "sample_name") + '/mobile_elements/Genomic_Islands')
                if event['data'] == 'crispr':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas'):
                        shutil.rmtree(
                            self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
                    shutil.copytree(self.crispr.output_dir + '/CRISPR_Cas',
                                    self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/CRISPR_Cas')
                if event['data'] == 'prephage':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
                    shutil.copytree(self.prephage.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/prephage')
                if event['data'] == 'pathogenic_system':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
                    shutil.copytree(self.pathogenic_system.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system')
                if event['data'] == 'is_predict':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Is_Predict'):
                        shutil.rmtree(
                            self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Is_Predict')
                    shutil.copytree(self.is_predict.output_dir + '/' + self.option("sample_name"),
                                    self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Is_Predict')
                if event['data'] == 'integron':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Integron'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Integron')
                    shutil.copytree(self.integron.output_dir + '/' + self.option("sample_name"),
                                    self.output_dir + '/' + self.option("sample_name") + '/mobile_elements/Integron')
                if event['data'] == 'repeatmasker':
                    if os.path.exists(self.output_dir + '/' + self.option(
                            "sample_name") + '/assembly_predict/predict/Interpersed_repeat'):
                        shutil.rmtree(self.output_dir + '/' + self.option(
                            "sample_name") + '/assembly_predict/predict/Interpersed_repeat')
                    shutil.copytree(self.repeatmasker.output_dir + '/' + self.option("sample_name"),
                                    self.output_dir + '/' + self.option(
                                        "sample_name") + '/assembly_predict/predict/Interpersed_repeat')
                    self.scaf_seq = self.com_pla.option('scf').prop['path']
                    if os.path.exists(self.scaf_seq):
                        if os.path.exists(self.output_dir + '/' + self.option(
                                "sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta'):
                            os.remove(self.output_dir + '/' + self.option(
                                "sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                        os.link(self.scaf_seq, self.output_dir + '/' + self.option(
                            "sample_name") + '/assembly_predict/predict/Interpersed_repeat/all.fasta')
                if event['data'] == 'real_gbk':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/real_gbk'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                    shutil.copytree(self.gbk_real_file.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/real_gbk')
                if event['data'] == 'resfinder':
                    if os.path.exists(
                            self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder'):
                        shutil.rmtree(
                            self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                    shutil.copytree(self.resfinder_predict.output_dir + '/' + self.option("sample_name"),
                                    self.output_dir + '/' + self.option("sample_name") + '/pathogenic_system/Resfinder')
                if event['data'] == 'srna_predict':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/srna_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                    shutil.copytree(self.srna_predict.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/srna_predict')
                if event['data'] == 'plasmid':
                    if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict'):
                        shutil.rmtree(self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')
                    shutil.copytree(self.com_pla.output_dir,
                                    self.output_dir + '/' + self.option("sample_name") + '/plasmid_predict')


    def end(self):
        if self.option('analysis') in ['uncomplete']:
            self.get_specimen()
        elif self.option('analysis') in ['complete']:
            self.get_complete_specimen()
        if not os.path.exists(self.output_dir + '/' + self.option("sample_name")):
            os.mkdir(self.output_dir + '/' + self.option("sample_name"))
        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/specimen.data'):
            os.remove(self.output_dir + '/' + self.option("sample_name") + '/specimen.data')
        if os.path.exists(self.output_dir + '/' + self.option("sample_name") + '/gene.data'):
            os.remove(self.output_dir + '/' + self.option("sample_name") + '/gene.data')
        self.logger.info(self.work_dir + '/specimen.data')
        os.link(self.work_dir + '/specimen.data', self.output_dir + '/' + self.option("sample_name") + '/specimen.data')
        os.link(self.work_dir + '/gene.data', self.output_dir + '/' + self.option("sample_name") + '/gene.data')
        super(BacgenomeModule, self).end()

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
                    
    def get_status(self, file, fli2):
        with open(fli2, "w") as g:
            for fa_iterator in SeqIO.parse(file, "fasta"):
                id = fa_iterator.id
                g.write("{}\t{}\n".format(str(id), "circle"))

    def get_sequence_type(self,file):
        list1 = []
        with open(file, "rb") as l:
            raw_lines = l.readlines()
            for line in raw_lines[1:]:
                line2 = line.strip('\r\n').split("\t")
                if line2[5] in ['pacbio', 'Pacbio']:
                    list1.append("Pacbio")
                elif line2[5] in ['nanopore', 'Nanopore']:
                    list1.append("Nanopore")
                else:
                    list1.append(line2[5])
        if len(set(list1)) == 1:
            sequence_type = list1[0]
        else:
            sequence_type = ','.join(set(list1))
        return sequence_type

    def get_assemble_type(self,file):
        list1 = []
        geneome_size = []
        with open(file, "rb") as l:
            raw_lines = l.readlines()
            for line in raw_lines[1:]:
                line2 = line.strip('\r\n').split("\t")
                if line2[5] in ['pacbio','Pacbio']:
                    list1.append("Pacbio")
                elif line2[5] in ['nanopore','Nanopore']:
                    list1.append("Nanopore")
                else:
                    list1.append(line2[5])
                geneome_size.append(line2[4])
        if len(set(list1)) == 1:
            sequence_type = list1[0]
        else:
            sequence_type = ','.join(set(list1))

        ## 获取基因组指定的大小，如果没有指定则定为4m
        new_size = []
        for size in geneome_size:
            if re.match('^[\d\.]+$',size):
                size = 1000000 * float(size)
            elif re.match('^[\d\.]+[kK]$',size):
                size = 1000 * float(size[:-1])
            elif re.match('^[\d\.]+[mM]$',size):
                size = 1000000 * float(size[:-1])
            else:
                continue
            new_size.append(size)
        if new_size:
            self.genome_size = min(new_size)
        else:
            self.genome_size = 5000000
        return sequence_type

    def get_seq(self):
        path=''
        if self.option('asse_dir').is_set:
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
            file =self.assemble_assess.option("output").prop['path']
            with open(file,'r') as f:
                lines =f.readlines()
                for line in lines[0:]:
                    line =line.rstrip('\r\n').split('\t')
                    prefix[line[0]]=line[1]
            de =str(prefix)
        return de

    def get_bases(self):
        tatol =0
        base_file =self.bac_qc.output_dir + '/data_QC/' + self.option('sample_name') + '_Illumina_statistics.xls'
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
                    if len(line2) == 6 or len(line2) == 7:
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
                    if len(line2) == 2 or len(line2) == 3:
                        file3.write(self.option('sample_name') + '\t' + '-' + '\t' + '-' + '\n')
                        file4.write(self.option('sample_name') + '\t' + line2[1] + '\t' + '-' + '\t' + 'gene' + '\n')

        elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
            raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
            ass_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
            output3 = self.work_dir + "/" + 'specimen.data'
            output4 = self.work_dir + "/" + 'gene.data'
            with open(raw_list_path, "rb") as l, open(ass_list_path, "rb") as f, open(output3, 'w') as file3, open(
                    output4, 'w') as file4:
                lines = l.readlines()
                lib_list = []
                raw_list = []
                for line in lines[1:]:
                    line2 = line.strip().split()
                    if len(line2) == 6 or len(line2) > 6:
                        lib_type = line2[5] + line2[2]
                        lib_list.append(lib_type)
                        raw_list.append(line2[1])
                lines2 = f.readlines()
                assem = ''
                for line in lines2[1:]:
                    line2 = line.strip().split()
                    if len(line2) == 2 or len(line2) == 3:
                        assem = line2[1]
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


    def get_complete_specimen(self):
        if not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
            asemble_fa =self.com_pla.option('scf').prop['path']
            base_name = os.path.basename(asemble_fa)
            output3 = self.work_dir + "/" + 'specimen.data'
            output4 = self.work_dir + "/" + 'gene.data'
            list1 =[]
            for fa_iterator in SeqIO.parse(asemble_fa, "fasta"):
                list1.append(fa_iterator.id)
            with open(output3, 'w') as file3, open(output4, 'w') as file4:
                file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                for i in list1:
                    if re.search("chromosome", i.lower()):
                        file4.write(self.option('sample_name') + '\t' + base_name + '\tchromosome\t' + 'gene' + '\n')
                    elif re.search(r'plasmid', i.lower()):
                        gene_pre = i[0].lower()+i[-1]+"_gene"
                        file4.write(self.option('sample_name') + '\t' + base_name + '\tplasmid\t' + gene_pre + '\n')
                file3.write(self.option('sample_name') + '\t' + '-' + '\t' + '-' + '\n')
        elif self.option("raw_dir").is_set and self.option("asse_dir").is_set:
            asemble_fa = self.com_pla.option('scf').prop['path']
            base_name = os.path.basename(asemble_fa)
            raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
            output3 = self.work_dir + "/" + 'specimen.data'
            output4 = self.work_dir + "/" + 'gene.data'
            list1 = []
            for fa_iterator in SeqIO.parse(asemble_fa, "fasta"):
                list1.append(fa_iterator.id)
            with open(raw_list_path, "rb") as l,  open(output3, 'w') as file3, open(
                    output4, 'w') as file4:
                lib_list = []
                raw_list = []
                file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                lines2 = l.readlines()
                for line in lines2[1:]:
                    line2 = line.strip().split()
                    if len(line2) == 6 or len(line2) > 6:
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
                for i in list1:
                    if re.search("chromosome", i.lower()):
                        file4.write(self.option('sample_name') + '\t' + base_name + '\tchromosome\t' + 'gene' + '\n')
                    elif re.search(r'plasmid', i.lower()):
                        gene_pre = i[0].lower()+i[-1]+"_gene"
                        file4.write(self.option('sample_name') + '\t' + base_name + '\tplasmid\t' + gene_pre + '\n')
        elif self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
            asemble_fa = self.com_pla.option('scf').prop['path']
            base_name = os.path.basename(asemble_fa)
            raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
            output3 = self.work_dir + "/" + 'specimen.data'
            output4 = self.work_dir + "/" + 'gene.data'
            list1 = []
            for fa_iterator in SeqIO.parse(asemble_fa, "fasta"):
                list1.append(fa_iterator.id)
            with open(raw_list_path, "rb") as l,  open(output3, 'w') as file3, open(
                    output4, 'w') as file4:
                lib_list = []
                raw_list = []
                file3.write('Sample Initial Name' + '\t' + 'File Name' + '\t' + 'Library' + '\n')
                file4.write('Sample Name' + '\t' + 'File Name' + '\t' + 'Genome Type' + '\t' + 'prefix_gene' + '\n')
                lines2 = l.readlines()
                for line in lines2[1:]:
                    line2 = line.strip().split()
                    if len(line2) == 6 or len(line2) > 6:
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
                for i in list1:
                    if re.search("chromosome", i.lower()):
                        file4.write(self.option('sample_name') + '\t' + base_name + '\tchromosome\t' + 'gene' + '\n')
                    elif re.search(r'plasmid', i.lower()):
                        gene_pre = i[0].lower()+i[-1]+"_gene"
                        file4.write(self.option('sample_name') + '\t' + base_name + '\tplasmid\t' + gene_pre + '\n')


    def get_antismash_gbk(self, gbk_dir):
        """
        功能是将完成图gbk_dir中间的一层location文件夹名称去掉，并去掉文件中的样本名称
        :param gbk_dir:
        :return:
        """
        new_gbk_dir = os.path.join(self.work_dir, "analysis_gbk")
        if os.path.exists(new_gbk_dir):
            shutil.rmtree(new_gbk_dir)
        os.mkdir(new_gbk_dir)
        if os.path.isdir(gbk_dir):
            alllists = os.listdir(gbk_dir)
            for gbk_d in alllists:
                dir_path = os.path.join(gbk_dir, gbk_d)
                if os.path.isdir(dir_path):
                    gbk_lists = os.listdir(dir_path)
                    # link_dir(dir_path, new_gbk_dir)
                    for gbk_d in gbk_lists:
                        new_gbk_file_name = gbk_d.split("_")[-1] ## 去掉样本名称
                        gbk_file_path = os.path.join(dir_path, gbk_d)
                        new_gbk_file_path = os.path.join(new_gbk_dir, new_gbk_file_name)
                        link_file(gbk_file_path, new_gbk_file_path)## 此功能函数允许多个文件夹链接到同一个文件夹下，前提是文件的名称不能是相同的，否则覆盖
                elif os.path.isfile(dir_path):
                    new_gbk_file_path = os.path.join(new_gbk_dir, gbk_d)
                    link_file(dir_path, new_gbk_file_path)
                else:
                    raise Exception("错误的文件夹类型")
        return new_gbk_dir

    def get_split_fasta(self, fasta):
        if os.path.exists(self.work_dir+"/seq_dir"):
            shutil.rmtree(self.work_dir+"/seq_dir")
        os.mkdir(self.work_dir+"/seq_dir")
        for fa_iterator in SeqIO.parse(fasta, "fasta"):
            id = fa_iterator.id
            SeqIO.write([fa_iterator], self.work_dir+"/seq_dir/"+id+".fasta", "fasta")
        return self.work_dir+"/seq_dir"