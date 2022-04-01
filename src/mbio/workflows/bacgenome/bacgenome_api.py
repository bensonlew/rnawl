# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""微生物基因组分析工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
import os,re
import json
import shutil
from biocluster.config import Config
import time
import datetime
import gevent
import functools
import shutil
import types
from collections import defaultdict


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

class BacgenomeApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        微生物基因组workflow option参数设置
        """
        self._sheet = wsheet_object
        super(BacgenomeApiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "raw_dir", "type": "infile", "format": "bacgenome.raw_dir"},  ###rawdata的文件目录
            {"name": "asse_dir", "type": "infile", "format": "bacgenome.asse_dir"},  ###组装的文件目录
            {"name": "analysis", "type": "string", "default": "uncomplete"}, ###流程分析模式complete，uncomplete
            {"name": "data_type", "type": "string", }
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.samples ={}
        self.modules = []
        self.trees = []
        #self.hgene_tree = self.add_batch('graph.phy_tree', ignore_error=True,batch_type="tool")
        self.hgene_tree = self.add_tool('graph.phy_tree')
        #self.tree = self.add_batch('graph.phy_tree', ignore_error=True, batch_type="tool")
        self.tree = self.add_tool('graph.phy_tree')
        self.list =[self.hgene_tree,self.tree]
        self.logger.info(self._sheet.output)
        self.step.add_steps('hgene_tree', 'tree')
        self.remote_dir = self._sheet.output + '/'

    def check_options(self):
        """
        检查参数
        """
        if not self.option("analysis"):
            raise OptionError("请提供流程分析模式！", code="11400101")
        if not self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
            raise OptionError("必须输入原始序列文件夹或组装序列文件夹其中一个！", code="11400102")


    def run(self):
        """
        运行:genome_workflow
        :return:
        """
        self.get_list()
        #self.run_api()
        gevent.spawn_later(5, self.end)
        super(BacgenomeApiWorkflow, self).run()


    def end(self):
        self.send_files()
        super(BacgenomeApiWorkflow, self).end()

    def run_api(self, test=False):
        print "aaaa"
        task_id =self._sheet.id
        project_sn =self._sheet.project_sn
        self.gbk = self.api.api("bacgenome.gbk")
        self.circos = self.api.api("bacgenome.circos")
        self.crispr = self.api.api("bacgenome.crispr")
        self.island = self.api.api("bacgenome.island")
        self.prephage = self.api.api("bacgenome.prephage")
        self.datastat = self.api.api("bacgenome.genome_qc")
        self.genome_size = self.api.api("bacgenome.genome_size")
        self.assess_gc = self.api.api("bacgenome.assess_gc")
        self.assemble = self.api.api("bacgenome.assemble")
        self.assess_kmer = self.api.api("bacgenome.assess_kmer")
        self.anno_antismash = self.api.api("bacgenome.anno_antismash")
        self.anno_card = self.api.api("bacgenome.anno_card")
        self.anno_cazy = self.api.api("bacgenome.anno_cazy")
        self.anno_cog = self.api.api("bacgenome.anno_cog")
        self.anno_go = self.api.api("bacgenome.anno_go")
        self.anno_kegg = self.api.api("bacgenome.anno_kegg")
        self.anno_nr = self.api.api("bacgenome.anno_nr")
        self.anno_pfam = self.api.api("bacgenome.anno_pfam")
        self.anno_phi = self.api.api("bacgenome.anno_phi")
        self.anno_ref = self.api.api("bacgenome.anno_ref")
        self.anno_summary = self.api.api("bacgenome.anno_summary")
        self.anno_swissprot = self.api.api("bacgenome.anno_swissprot")
        self.anno_tcdb = self.api.api("bacgenome.anno_tcdb")
        self.anno_tmhmm = self.api.api("bacgenome.anno_tmhmm")
        self.anno_vfdb = self.api.api("bacgenome.anno_vfdb")
        self.gene_graph = self.api.api("bacgenome.gene_graph")
        self.gene_predict = self.api.api("bacgenome.gene_predict")
        self.promote = self.api.api("bacgenome.promote")
        self.repeat_predict = self.api.api("bacgenome.repeat_predict")
        self.rrna_predict = self.api.api("bacgenome.rrna_predict")
        self.secretory = self.api.api("bacgenome.secretory")
        self.summary_map = self.api.api("bacgenome.summary_map")
        self.trna_predict = self.api.api("bacgenome.trna_predict")
        if self.option('raw_dir').is_set:
            seq_type = self.get_assemble_type(self.option('raw_dir').prop['path'] + '/list.txt')
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set:
                datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                         "Trimmomatic, SeqPrep, Sickle, FastqTotalHighQualityBase.jar",
                                                         "30", "20")
            elif self.option("asse_dir").is_set:
                datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', "", "assemble", "", "",
                                                         "")
        elif self.option("analysis") in ["complete"]:
            if self.option("raw_dir").is_set:
                if re.search(r'PE', seq_type):
                    datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                             "Trimmomatic, SeqPrep, Sickle, FastqTotalHighQualityBase.jar;smrtanalysis",
                                                             "30", "20")
                else:
                    datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                             "smrtanalysis", "", "")
            elif self.option("asse_dir").is_set:
                datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', "", "assemble", "", "",
                                                         "")
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set:
                assemble_id = self.assemble.add_assemble(self.option('analysis'), seq_type, 'rawdata', '组装评估',
                                                         'kmer (21-41)', 'SOAPdenovo v2.04,GapCloser v1.12')
            elif self.option("asse_dir").is_set:
                assemble_id = self.assemble.add_assemble(self.option('analysis'), '', 'assemble', '组装评估', '', '')
        elif self.option("analysis") in ["complete"]:
            if self.option("raw_dir").is_set:
                assemble_id = self.assemble.add_assemble(self.option('analysis'), seq_type, 'rawdata', '组装评估', '', '')
            elif self.option("asse_dir").is_set:
                assemble_id = self.assemble.add_assemble(self.option('analysis'), '', 'assemble', '组装评估', '', '')
        if self.option("analysis") in ["uncomplete"]:
            gene_id = self.gene_predict.add_gene_predict(self.remote_dir , '/assembly_predict/predict/CDS_predict/',
                                                         '_CDS.', params='{"soft": "Glimmer"}')
        elif self.option("analysis") in ["complete"]:
            gene_id = self.gene_predict.add_gene_predict(self.remote_dir , '/assembly_predict/predict/CDS_predict/',
                                                         '_whole_genome_CDS.', params='{"soft": "Glimmer and GeneMarkS"}')
        trna_id = self.trna_predict.add_trna_predict(self.remote_dir, '/assembly_predict/predict/tRNA/', '_tRNA.',
                                                     params='soft:tRNAscan-SE')
        rrna_id = self.rrna_predict.add_rrna_predict(self.remote_dir, '/assembly_predict/predict/rRNA/', '_rRNA.',
                                                     params='soft:barrnap')
        repeat_id = self.repeat_predict.add_repeat_predict(params='soft:TRF')
        nr_id = self.anno_nr.add_anno_nr(params='{"NR":"Diamond"}')
        cog_id = self.anno_cog.add_anno_cog(params="{COG:Diamond}")
        gbk_id = self.gbk.add_gbk(project_sn=project_sn, task_id=task_id)
        go_id = self.anno_go.add_anno_go(params="{GO:blast2go}")
        antismash_id =self.anno_antismash.add_anno_antismash(params="{antismash:antismash}")
        cazy_id = self.anno_cazy.add_anno_cazy(params="{CAZy:hmmscan}")
        kegg_id = self.anno_kegg.add_anno_kegg(self.remote_dir, "/annotation/KEGG/", "_kegg_pathway_img", params="{KEGG:Diamond}")
        pfam_id = self.anno_pfam.add_anno_pfam(params="{Pfam:hmmer3}")
        swissport_id = self.anno_swissprot.add_anno_swissprot(params="{Swiss-prot:blast}")
        summary_id = self.anno_summary.add_anno_summary(params="{NR,KEGG,COG,GO,Swissprot,Pfam}")
        island_id = self.island.add_island("isand")
        prephage_id = self.prephage.add_prephage("prephage")
        crispr_id = self.crispr.add_crispr("crispr")
        gc_id = self.assess_gc.add_assess_gc("GC_Depth",link_path =self.remote_dir)
        size_id = self.genome_size.add_assess_size("genome size")
        n = 1
        print 'bbbbb'
        for sample in self.samples.keys():
            print 'cccc'
            print sample
            self.logger.info(">>>start wait_file,file route: %s" % (self.output_dir + '/' + sample))
            self.logger.info(">>>wait_file end")
            if self.option("analysis") in ["uncomplete"]:
                self.datastat.add_datastat_specimen(datastat_id,
                                                    self.output_dir + '/' + sample + '/' + 'specimen.data')
                self.datastat.add_datastat_uncomplete_gene(datastat_id,
                                                           self.output_dir + '/' + sample + '/' + 'gene.data')
            elif self.option("analysis") in ["complete"]:
                self.datastat.add_datastat_specimen(datastat_id,
                                                    self.output_dir + '/' + sample + '/' + 'specimen.data')
                self.datastat.add_datastat_complete_gene(datastat_id,
                                                         self.output_dir + '/' + sample + '/' + 'gene.data')
            if self.option("analysis") in ["uncomplete"]:
                gbk_fi = os.path.join(self.remote_dir,sample + "/project_overview/gbk/")
                self.gbk.add_gbk_detail(path = self.output_dir + '/' + sample + '/gbk/gbk',main_id =gbk_id,sample_name= sample,gbk_path= gbk_fi,task_id=task_id)
            elif self.option("analysis") in ["complete"]:
                gbk_fi = os.path.join(self.remote_dir,sample + "/project_overview/gbk/")
                self.gbk.add_gbk_detail(path=self.output_dir + '/' + sample + '/gbk/seq_gbk', main_id=gbk_id,sample_name=sample, gbk_path=gbk_fi, task_id=task_id)
            if self.option("analysis") in ["uncomplete"]:
                self.datastat.add_datastat_uncomplete_summary(datastat_id,
                                                              self.output_dir + '/' + sample + '/' + sample + '.project.summary')
            elif self.option("analysis") in ["complete"]:
                self.datastat.add_datastat_complete_summary(datastat_id,
                                                            self.output_dir + '/' + sample + '/' + sample + '.project.summary')
            if self.option("analysis") in ["uncomplete"]:
                for type in ['hgene']:
                    self.datastat.add_datastat_tree_detail(datastat_id,
                                                           self.output_dir + '/' + sample + '/tree/' + type + '/' + sample + '.phylo_tree.nwk',
                                                           sample, type, 'single')
            elif self.option("analysis") in ["complete"]:
                for type in ['16s', 'hgene']:
                    self.datastat.add_datastat_tree_detail(datastat_id,
                                                           self.output_dir + '/' + sample + '/tree/' + type + '/' + sample + '.phylo_tree.nwk',
                                                           sample, type, 'single')
            if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                my_seq = self.get_lib_type(self.option('raw_dir').prop['path'] + '/list.txt')
                self.logger.info(my_seq)
                if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls'):
                    self.datastat.add_qc_stat_uncomplete(datastat_id,
                                                         self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls')
                if sample in my_seq.keys():
                    self.logger.info(my_seq[sample].keys())
                    for type in my_seq[sample].keys():
                        self.logger.info(type)
                        for s in ['raw', 'clean']:
                            self.datastat.add_datastat_graphic(datastat_id,
                                                               self.output_dir + '/' + sample + '/fastx/' + type + '_l.' + s + '_fastxstat',
                                                               sample, s, "left", type)
                            self.datastat.add_datastat_graphic(datastat_id,
                                                               self.output_dir + '/' + sample + '/fastx/' + type + '_r.' + s + '_fastxstat',
                                                               sample, s, "right", type)

            elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                my_seq = self.get_lib_type(self.option('raw_dir').prop['path'] + '/list.txt')
                if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls'):
                    self.datastat.add_qc_stat_uncomplete(datastat_id,
                                                              self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls')
                    if sample in my_seq.keys():
                        for type in my_seq[sample].keys():
                            for s in ['raw', 'clean']:
                                self.datastat.add_datastat_graphic(datastat_id,
                                                                   self.output_dir + '/' + sample + '/fastx/' + type + '_l.' + s + '_fastxstat',
                                                                   sample, s, "left", type)
                                self.datastat.add_datastat_graphic(datastat_id,
                                                                   self.output_dir + '/' + sample + '/fastx/' + type + '_r.' + s + '_fastxstat',
                                                                   sample, s, "right", type)
                if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '.PacBio_statistics.xls'):
                    self.datastat.add_qc_stat_complete(datastat_id,
                                                            self.output_dir + '/' + sample + '/data_QC/' + sample + '.PacBio_statistics.xls',
                                                            'pacbio',sample)
                    for s in ['raw', 'clean']:
                        if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '.' + s + '.len.xls'):
                            self.datastat.add_datastat_pacbio_graphic(datastat_id,self.output_dir + '/' + sample + '/data_QC/' +sample+ '.' + s + '.len.xls',sample, s, "len")
                        if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + 'pacbio.' + s + '.qv.xls'):
                            self.datastat.add_datastat_pacbio_graphic(datastat_id, self.output_dir + '/' + sample + '/data_QC/' + 'pacbio.' + s + '.qv.xls',sample, s, "qv")
            if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "1k",self.remote_dir + sample + "/genomic_assessment/depth_gc_1000/")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "3k",self.remote_dir + sample + "/genomic_assessment/depth_gc_3000/")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "5k",self.remote_dir + sample + "/genomic_assessment/depth_gc_5000/")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "8k",self.remote_dir + sample + "/genomic_assessment/depth_gc_8000/")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "10k",self.remote_dir + sample + "/genomic_assessment/depth_gc_10000/")
                kmer_id = self.assess_kmer.add_assess_kmer(sample, "kmer")
                self.assess_kmer.add_assess_kmer_detail(kmer_id,
                                                        self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + ".frequency.xls")
                self.genome_size.add_assess_size_detail(size_id,
                                                        self.output_dir + '/' + sample + "/genomic_assessment/genome_size/" + sample + ".summary.xls",
                                                        sample)
            elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                if re.search(r'PE', seq_type):
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "1k",self.remote_dir + sample + "/genomic_assessment/depth_gc_1000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "3k",self.remote_dir + sample + "/genomic_assessment/depth_gc_3000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "5k",self.remote_dir + sample + "/genomic_assessment/depth_gc_5000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "8k",self.remote_dir + sample + "/genomic_assessment/depth_gc_8000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "10k",self.remote_dir + sample + "/genomic_assessment/depth_gc_10000/")
                    kmer_id = self.assess_kmer.add_assess_kmer(sample, "kmer")
                    self.assess_kmer.add_assess_kmer_detail(kmer_id,
                                                            self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + ".frequency.xls")
                    self.genome_size.add_assess_size_detail(size_id,
                                                            self.output_dir + '/' + sample + "/genomic_assessment/genome_size/" + sample + ".summary.xls",
                                                            sample)
            if self.option("analysis") in ["uncomplete"]:
                self.assemble.add_assemble_stat_uncomplete(assemble_id,
                                                           self.output_dir + '/' + sample + '/assembly_predict/assembly/assembly/' + sample + '_assembly_summary.xls',
                                                           sample)
                for type in ['scaffold', 'contig']:
                    self.assemble.add_assemble_seq(assemble_id,
                                                   self.output_dir + '/' + sample + '/assembly_predict/assembly/assembly/' + sample + '_assembly_' + type + '_details.xls',
                                                   sample, type,self.remote_dir)
                    for window in ['1k', '2k', '5k']:
                        if window == '1k':
                            self.assemble.add_assemble_graphic(assemble_id,
                                                               self.output_dir + '/' + sample + '/assembly_predict/assembly/len/' + sample + '.1000.' + type + 's.len.xls',
                                                               sample, type, '1k')
                        elif window == '2k':
                            self.assemble.add_assemble_graphic(assemble_id,
                                                               self.output_dir + '/' + sample + '/assembly_predict/assembly/len/' + sample + '.2000.' + type + 's.len.xls',
                                                               sample, type, '2k')
                        elif window == '5k':
                            self.assemble.add_assemble_graphic(assemble_id,
                                                               self.output_dir + '/' + sample + '/assembly_predict/assembly/len/' + sample + '.5000.' + type + 's.len.xls',
                                                               sample, type, '5k')
            elif self.option("analysis") in ["complete"]:
                self.assemble.add_assemble_stat_complete(assemble_id,
                                                         self.output_dir + '/' + sample + '/assembly_predict/assembly/' + sample + '_assembly_summary.xls',
                                                         sample)
                self.assemble.add_assemble_complete_seq(assemble_id,
                                                        self.output_dir + '/' + sample + '/assembly_predict/assembly/' + sample + '_assembly_details.xls',
                                                        sample, self.remote_dir)
            if self.option("analysis") in ["uncomplete"]:
                self.gene_predict.add_gene_predict_detail(gene_id, sample,
                                                          self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS.gff")
                self.gene_predict.add_gene_predict_specimen(gene_id, sample,
                                                            self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS_statistics.xls")
                self.gene_predict.add_gene_predict_seq(gene_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS.fnn",
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS.faa")
                self.gene_predict.add_gene_predict_bar(gene_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/' + "length_distribute.txt")
            elif self.option("analysis") in ["complete"]:
                self.gene_predict.add_gene_predict_detail(gene_id, sample,
                                                          self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS.gff")
                self.gene_predict.add_gene_predict_specimen(gene_id, sample,
                                                            self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS_statistics.xls")
                self.gene_predict.add_gene_predict_seq(gene_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS.fnn",
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS.faa")
                self.gene_predict.add_gene_predict_bar(gene_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/' + "length_distribute.txt")
            self.trna_predict.add_trna_predict_detail(trna_id, sample,
                                                      self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.gff")
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.fnn"):
                self.trna_predict.add_trna_predict_seq(trna_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.fnn",
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.struc")
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/assembly_predict/predict/rRNA/' + sample + "_rRNA.fnn"):
                self.rrna_predict.add_rrna_predict_seq(rrna_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/rRNA/' + sample + "_rRNA.fnn")
            self.rrna_predict.add_rrna_predict_detail(rrna_id, sample,
                                                      self.output_dir + '/' + sample + '/assembly_predict/predict/rRNA/' + sample + "_rRNA.gff")
            if self.option("analysis") in ["uncomplete"]:
                self.repeat_predict.add_repeat_predict_detail(repeat_id, sample,
                                                              self.output_dir + '/' + sample + '/assembly_predict/predict/repeats/' + sample + "_TRF.gff")
            elif self.option("analysis") in ["complete"]:
                self.repeat_predict.add_repeat_predict_detail(repeat_id, sample,
                                                              self.output_dir + '/' + sample + '/assembly_predict/predict/repeats/' + sample + "_whole_genome_TRF.gff")
            if self.option("analysis") in ["uncomplete"]:
                self.anno_nr.add_anno_nr_detail(nr_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/NR/" + sample + "_anno_nr.xls")
            elif self.option("analysis") in ["complete"]:
                self.anno_nr.add_anno_nr_detail(nr_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/NR/" + sample + "_whole_genome_anno_nr.xls")
            self.anno_cog.add_anno_cog_detail(cog_id, sample,
                                              self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls")
            self.anno_go.add_anno_go_specimen(go_id, sample,
                                              self.output_dir + '/' + sample + "/annotation/GO/" + sample + "_go_level2_statistics.xls")
            self.anno_go.add_anno_go_detail(go_id, sample,
                                              self.output_dir + '/' + sample + "/annotation/GO/" + sample + "_go_anno.xls")
            self.anno_cazy.add_anno_cazy_detail(cazy_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/CAZy/" + sample + "_anno_cazy.xls")
            self.anno_kegg.add_anno_kegg_detail(kegg_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls")
            self.anno_kegg.add_anno_kegg_level(kegg_id, sample,
                                               self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_level_stat.xls")
            if self.option("analysis") in ["uncomplete"]:
                self.anno_pfam.add_anno_pfam_detail(pfam_id, sample,
                                                    self.output_dir + '/' + sample + "/annotation/Pfam/" + sample + "_anno_pfam.xls")
            elif self.option("analysis") in ["complete"]:
                self.anno_pfam.add_anno_pfam_detail(pfam_id, sample,
                                                    self.output_dir + '/' + sample + "/annotation/Pfam/" + sample + "_whole_genome_anno_pfam.xls")
            if self.option("analysis") in ["uncomplete"]:
                self.anno_swissprot.add_anno_swissprot_detail(swissport_id, sample,
                                                              self.output_dir + '/' + sample + "/annotation/Swissprot/" + sample + "_anno_swissprot.xls")
            elif self.option("analysis") in ["complete"]:
                self.anno_swissprot.add_anno_swissprot_detail(swissport_id, sample,
                                                              self.output_dir + '/' + sample + "/annotation/Swissprot/" + sample + "_whole_genome_anno_swissprot.xls")
            self.anno_summary.add_anno_summary_detail(summary_id, sample,
                                                      self.output_dir + '/' + sample + "/annotation/Summary/" + sample + "_anno_summary.xls")
            if self.option("analysis") in ["complete"]:
                if os.path.exists(self.output_dir + '/' + sample + "/metabolic_system/antiSMASH"):
                    files = os.listdir(self.output_dir + '/' + sample + "/metabolic_system/antiSMASH")
                    if len(files) != 0:
                        for file in files:
                            self.anno_antismash.add_anno_antismash_stat(antismash_id, sample,self.output_dir + '/' + sample + '/metabolic_system/antiSMASH/' + file)

            elif self.option("analysis") in ["uncomplete"]:
                if os.path.exists(self.output_dir + '/' + sample + "/metabolic_system/antiSMASH/antismash_anno.xls"):
                    self.anno_antismash.add_anno_antismash_stat(antismash_id, sample,self.output_dir + '/' + sample + "/metabolic_system/antiSMASH/antismash_anno.xls")
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_summary.xls"):
                self.island.add_island_detail(island_id,
                                              self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_summary.xls",
                                              sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_detail.xls"):
                self.island.add_island_gene(island_id,
                                            self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_detail.xls",
                                            sample)
            if os.path.exists(self.output_dir + '/' + sample + '/mobile_elements/prophage/' + sample + ".stat.xls"):
                self.prephage.add_prephage_stat(prephage_id,
                                                self.output_dir + '/' + sample + '/mobile_elements/prophage/' + sample + ".stat.xls",
                                                sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/prophage/' + sample + "_prophage_summary.xls"):
                self.prephage.add_prephage_detail(prephage_id,
                                                  self.output_dir + '/' + sample + '/mobile_elements/prophage/' + sample + "_prophage_summary.xls",
                                                  sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/prophage/' + sample + "_prophage_detail.xls"):
                self.prephage.add_prephage_gene(prephage_id,
                                                self.output_dir + '/' + sample + '/mobile_elements/prophage/' + sample + "_prophage_detail.xls",
                                                sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_summary.xls"):
                self.crispr.add_crispr_stat(crispr_id,
                                            self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_summary.xls",
                                            sample)
                self.crispr.add_crispr_detail(crispr_id,
                                              self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_summary.xls",
                                              sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_detail.xls"):
                self.crispr.add_crispr_psa(crispr_id,
                                           self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_detail.xls",
                                           sample)
            if n == 1:
                prom_id = self.promote.add_promote(
                    self.output_dir + '/' + sample + "/structral_genome/promoter_predict",
                    main=True, specimen_id=sample, update_id=summary_id)
                self.logger.info("promoter end")
                self.api_gene_graph = self.api.api("bacgenome.gene_graph")
                if self.option("analysis") in ["uncomplete"]:
                    gene_graph_id = self.gene_graph.add_gene_graph(
                        self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_CDS.fnn",
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls", main=True,
                        specimen_id=sample)
                    self.logger.info("gene_graph end")
                elif self.option("analysis") in ["complete"]:
                    gene_graph_id = self.gene_graph.add_gene_graph(
                        self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_whole_genome_CDS.fnn",
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls",
                        main=True, specimen_id=sample)
                    self.logger.info("gene_graph end")
                vfdb_id = self.anno_vfdb.add_anno_vfdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/VFDB", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("vfdb end")
                phi_id = self.anno_phi.add_anno_phi(
                    self.output_dir + '/' + sample + "/pathogenic_system/PHI", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("phi end")
                tcdb_id = self.anno_tcdb.add_anno_tcdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/TCDB", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("tcdb end")
                self.api_card = self.api.api("bacgenome.anno_card")
                card_id = self.anno_card.add_anno_card(
                    self.output_dir + '/' + sample + "/pathogenic_system/CARD", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("card end")
                tmhmm_id = self.anno_tmhmm.add_anno_tmhmm(
                    self.output_dir + '/' + sample + "/pathogenic_system/TMHMM", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("tmhmm end")
                secretory_id = self.secretory.add_secretory(
                    self.output_dir + '/' + sample + "/pathogenic_system/secretion_system",
                    anno_summary_id=summary_id, kegg_id=kegg_id, main=True, specimen_id=sample)
                n +=1
            else:
                self.promote.add_promote(
                    self.output_dir + '/' + sample + "/structral_genome/promoter_predict",
                    main_id=prom_id, specimen_id=sample, update_id=summary_id)
                self.logger.info("promoter end")
                if self.option("analysis") in ["uncomplete"]:
                    self.gene_graph.add_gene_graph(
                        self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_CDS.fnn",
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls", main_id=gene_graph_id,
                        specimen_id= sample)
                    self.logger.info("gene_graph end")
                elif self.option("analysis") in ["complete"]:
                    self.gene_graph.add_gene_graph(
                        self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_whole_genome_CDS.fnn",
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls",
                        main_id=gene_graph_id, specimen_id=sample)
                    self.logger.info("gene_graph end")
                self.anno_vfdb.add_anno_vfdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/VFDB", main_id=vfdb_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("vfdb end")
                self.anno_phi.add_anno_phi(
                    self.output_dir + '/' + sample + "/pathogenic_system/PHI", main_id=phi_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("phi end")
                self.anno_tcdb.add_anno_tcdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/TCDB", main_id=tcdb_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("tcdb end")
                self.anno_card.add_anno_card(
                    self.output_dir + '/' + sample + "/pathogenic_system/CARD", main_id=card_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("card end")
                self.anno_tmhmm.add_anno_tmhmm(
                    self.output_dir + '/' + sample + "/pathogenic_system/TMHMM", main_id=tmhmm_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("tmhmm end")
                self.secretory.add_secretory(self.output_dir + '/' + sample + "/pathogenic_system/secretion_system",
                                             anno_summary_id=summary_id, kegg_id=kegg_id, main_id=secretory_id,
                                             specimen_id=sample)
        self.logger.info("secretory end")
        self.api_circos_table = self.api.api("bacgenome.circos_table")
        self.api_circos_table.add_circos_table(project_sn=project_sn, task_id=task_id)
        self.logger.info("circos_table end")
        self.summary_map = self.api.api("bacgenome.summary_map")
        self.summary_map.add_map(task_id=task_id)
        self.circos = self.api.api("bacgenome.circos")
        if self.option("analysis") in ["complete"]:
            seq = self.get_seq_type()
            for key, value in seq.iteritems():
                locatin_list = value
                for location in self.re_location(locatin_list):
                    params = {'specimen_id': key, 'task_type': 2, 'submit_location': 'circos',
                              'task_id': task_id,
                              'location': location}
                    name = key + '_Circos_' + '_' + location + '_' + datetime.datetime.now().strftime(
                        "%Y%m%d_%H%M%S%f")[:-3]
                    self.circos.add_circos(
                        self.output_dir + '/' + key + '/circos/' + location, project_sn=project_sn, name=name,
                        location=location, task_id=task_id, main=True, specimen_id=key, params=params,
                        link_path=self.remote_dir + sample + '/circos/'+ location + '/')
        if self.option("analysis") in ["uncomplete"]:
            for sample in self.samples.keys():
                params = {'specimen_id': sample, 'task_type': 2, 'submit_location': 'circos', 'task_id': task_id}
                name = sample + '_Circos_' + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
                self.circos.add_circos(
                    self.output_dir + '/' + sample + '/circos/Scaffold', project_sn=project_sn, name=name,
                    task_id=task_id, main=True, specimen_id=sample,
                    params=params, link_path=self.remote_dir + sample + '/circos/Scaffold/')
        self.cgview = self.api.api("bacgenome.cgview")
        if self.option("analysis") in ["uncomplete"]:
            for sample in self.samples.keys():
                params = {'genome_type': 'Scaffold', 'species_name': '', 'specimen_id': sample,
                          "submit_location": 'cgview',
                          "task_type": 2}
                file = self.output_dir + '/' + sample + '/cgview'
                name = sample + '_Cgview_' + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
                self.cgview.add_cgview(file, params=params, project_sn=project_sn,specimen_id=sample,genome_type= 'scaffold',
                                       task_id=task_id, main='true', name=name, link_path=self.remote_dir + sample + '/cgview/')
        if self.option("analysis") in ["complete"]:
            for sample in self.samples.keys():
                file = self.output_dir + '/' + sample + '/cgview'
                files = os.listdir(file)
                for fi in files:
                    tmp = ''
                    if fi.startswith('Chromosome'):
                        name = 'Chr' + fi[10:]
                    elif fi.startswith('Plasmid'):
                        name = 'p' + fi[7:]
                    else:
                        name = fi
                    fil = self.output_dir + '/' + sample + '/cgview/' + fi
                    params = {'genome_type': name, 'species_name': '', 'specimen_id': sample,
                              "submit_location": 'cgview',
                              "task_type": 2}
                    name2 = sample + '_Cgview_' + '_' + fi + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[
                                                                   :-3]
                    self.cgview.add_cgview(fil, genome_type=name, specimen_id=sample, name=name2, params=params,
                                           project_sn=project_sn,
                                           task_id=task_id, main='true',
                                           link_path=self.remote_dir + sample + '/cgview/')
        if self.option("analysis") in ["uncomplete"]:
            if os.path.exists(self.hgene_tree.work_dir + '/phylo_tree.nwk'):
                    self.datastat.add_datastat_tree_detail(datastat_id,self.hgene_tree.work_dir + '/phylo_tree.nwk',
                                                       '', 'hgene', 'all')
        if self.option("analysis") in ["complete"]:
            if os.path.exists(self.hgene_tree.work_dir + '/phylo_tree.nwk'):
                self.datastat.add_datastat_tree_detail(datastat_id, self.hgene_tree.work_dir + '/phylo_tree.nwk',
                                                       '', 'hgene', 'all')
            if os.path.exists(self.tree.work_dir + '/phylo_tree.nwk'):
                self.datastat.add_datastat_tree_detail(datastat_id,self.tree.work_dir + '/phylo_tree.nwk',
                                                       '', '16s', 'all')

    def re_location(self,location_list):
        rename = []
        for one in location_list.split(','):
            if one.startswith('Chromosome'):
                name = 'Chr' + one[10:]
            elif one.startswith('Plasmid'):
                name = 'p' + one[7:]
            else:
                name = one
            rename.append(name)
        return rename


    def send_files(self):
        """
        结果放置到/upload_results
        """
        dir_o = self.output_dir
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        repaths = []
        regexps = []
        for sample in self.samples.keys():
            if self.option("analysis") in ["uncomplete"]:
                files = os.listdir(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly/"))
                for file in files:
                    if re.search(r'.fna.index.', file):
                        os.remove(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly/" + file))
            if self.option("analysis") in ["uncomplete"]:
                if os.path.exists(os.path.join(dir_o, sample + "/gbk")):
                    self.move_dir(os.path.join(dir_o, sample + "/gbk/gbk/Scaffold"), os.path.join(dir_up, sample + "/project_overview/gbk"))
            elif self.option("analysis") in ["complete"]:
                if os.path.exists(os.path.join(dir_o, sample + "/gbk")):
                    files = os.listdir(os.path.join(dir_o, sample + "/gbk/seq_gbk"))
                    for file in files:
                        self.move_file(os.path.join(dir_o, sample + "/gbk/seq_gbk/" + file+ '/' + sample + '_' + file + '.gbk'), os.path.join(dir_up, sample + "/project_overview/gbk/"  + sample + '_' + file + '.gbk'))
            if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                if os.path.exists(os.path.join(dir_o, sample + "/data_QC")):
                    self.move_dir(os.path.join(dir_o, sample + "/data_QC"), os.path.join(dir_up, sample + "/data_QC"))
            elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                if os.path.exists(os.path.join(dir_o, sample + "/data_QC")):
                    self.move_dir(os.path.join(dir_o, sample + "/data_QC"), os.path.join(dir_up, sample + "/data_QC"))
            if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency")):
                    self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency"),os.path.join(dir_up, sample + "/genomic_assessment/kmer_frequency"))
                    for i in [1000,3000,5000,8000,10000]:
                        self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/depth_gc_" + str(i)),os.path.join(dir_up, sample + "/genomic_assessment/depth_gc_" + str(i)))
            elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency")):
                    self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency"),os.path.join(dir_up, sample + "/genomic_assessment/kmer_frequency"))
                    for i in [1000,3000,5000,8000,10000]:
                        self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/depth_gc_" + str(i)),os.path.join(dir_up, sample + "/genomic_assessment/depth_gc_" + str(i)))
            if self.option("analysis") in ["uncomplete"]:
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly")):
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly"),os.path.join(dir_up, sample + "/assembly_predict/assembly"))
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/scaffold"),os.path.join(dir_up,sample + "/assembly_predict/assembly/scaffold/"))
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/contig"),os.path.join(dir_up, sample + "/assembly_predict/assembly/contig"))
            elif self.option("analysis") in ["complete"]:
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/")):
                    self.move_file(os.path.join(dir_o, sample + "/assembly_predict/assembly/" + sample + '_assembly_summary.xls'),os.path.join(dir_up, sample + "/assembly_predict/assembly/"  + sample + '_assembly_summary.xls'))
                    self.move_file(os.path.join(dir_o, sample + "/assembly_predict/assembly/" + sample + '_assembly_details.xls'),os.path.join(dir_up, sample + "/assembly_predict/assembly/" + sample + '_assembly_details.xls'))
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/seq_dir")):
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/seq_dir"),os.path.join(dir_up, sample + "/assembly_predict/assembly/seq_dir"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/CDS_predict")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/CDS_predict"),os.path.join(dir_up, sample + "/assembly_predict/predict/CDS_predict"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/tRNA")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/tRNA"),os.path.join(dir_up, sample + "/assembly_predict/predict/tRNA"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/rRNA")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/rRNA"),os.path.join(dir_up, sample + "/assembly_predict/predict/rRNA"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/repeats")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/repeats"),os.path.join(dir_up, sample + "/assembly_predict/predict/repeats"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Swissprot")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Swissprot"),os.path.join(dir_up, sample + "/annotation/Swissprot"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Summary")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Summary"),os.path.join(dir_up, sample + "/annotation/Summary"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Pfam")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Pfam"),os.path.join(dir_up, sample + "/annotation/Pfam"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/KEGG")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/KEGG"),os.path.join(dir_up, sample + "/annotation/KEGG"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/KEGG/gene_cazy_class_stat.xls")):
                os.link(os.path.join(dir_o, sample + "/annotation/KEGG/kegg_pathway_img.tar.gz"),
                        os.path.join(dir_up, sample + "/annotation/KEGG/kegg_pathway_img.tar.gz"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/GO")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/GO"),os.path.join(dir_up, sample + "/annotation/GO"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/NR")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/NR"),os.path.join(dir_up, sample + "/annotation/NR"))  #guanqing.zou 20180903
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/COG")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/COG"),os.path.join(dir_up, sample + "/annotation/COG"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/CAZy")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/CAZy"),os.path.join(dir_up, sample + "/metabolic_system/CAZy"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/gene_cazy_family_stat.xls")):
                os.link(os.path.join(dir_o, sample + "/annotation/gene_cazy_family_stat.xls"),
                              os.path.join(dir_up, sample + "/metabolic_system/CAZy/" + sample + '_cazy_family_stat.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/gene_cazy_class_stat.xls")):
                os.link(os.path.join(dir_o, sample + "/annotation/gene_cazy_class_stat.xls"),
                              os.path.join(dir_up, sample + "/metabolic_system/CAZy/" + sample + '_cazy_class_stat.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/structral_genome/promoter_predict")):
                self.move_dir(os.path.join(dir_o, sample + "/structral_genome/promoter_predict"), os.path.join(dir_up, sample + "/structral_genome/promoter_predict"))
            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/Genomic_Islands")):
                self.move_dir(os.path.join(dir_o, sample + "/mobile_elements/Genomic_Islands"), os.path.join(dir_up, sample + "/mobile_elements/Genomic_Islands"))
            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/CRISPR_Cas")):
                self.move_dir(os.path.join(dir_o, sample + "/mobile_elements/CRISPR_Cas"), os.path.join(dir_up, sample + "/mobile_elements/CRISPR_Cas"))
            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prophage")):
                if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prophage/" + sample + '_prophage_summary.xls')):
                    self.move_file(os.path.join(dir_o, sample + "/mobile_elements/prophage/" + sample + '_prophage_summary.xls'),
                      os.path.join(dir_up, sample + "/mobile_elements/prophage/" + sample + '_prophage_summary.xls'))
                if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prophage/" + sample + '_prophage_detail.xls')):
                    self.move_file(os.path.join(dir_o, sample + "/mobile_elements/prophage/" + sample + '_prophage_detail.xls'),
                    os.path.join(dir_up, sample + "/mobile_elements/prophage/" + sample + '_prophage_detail.xls'))
                if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prophage/" + sample + '_prophage.fna')):
                    self.move_file(os.path.join(dir_o, sample + "/mobile_elements/prophage/" + sample + '_prophage.fna'),
                    os.path.join(dir_up, sample + "/mobile_elements/prophage/" + sample + '_prophage.fna'))
            if os.path.exists(os.path.join(dir_o, sample + "/metabolic_system/antiSMASH")):
                self.move_file(os.path.join(dir_o, sample + "/metabolic_system/antiSMASH/antismash_anno.xls"), os.path.join(dir_up, sample + "/metabolic_system/antiSMASH/antismash_anno.xls"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/VFDB")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/VFDB"), os.path.join(dir_up, sample + "/pathogenic_system/VFDB"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/TMHMM")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/TMHMM"), os.path.join(dir_up, sample + "/pathogenic_system/TMHMM"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/TCDB")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/TCDB"), os.path.join(dir_up, sample + "/pathogenic_system/TCDB"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/PHI")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/PHI"), os.path.join(dir_up, sample + "/pathogenic_system/PHI"))
            if os.path.exists(os.path.join(dir_o, sample + "/cgview")):
                self.move_dir(os.path.join(dir_o, sample + "/cgview"), os.path.join(dir_up, sample + "/cgview"))
            if os.path.exists(os.path.join(dir_o, sample + "/circos")):
                self.move_dir(os.path.join(dir_o, sample + "/circos"),os.path.join(dir_up, sample + "/circos"))
            repaths += [
                ["%s" % sample, "", "样品%s流程分析结果目录" % sample,0,'130001'],
                ["%s/project_overview" % sample, "", "样品%s的基因组总览目录" % sample,0,'130002'],
                ["%s/project_overview/gbk" % sample, "", "样品%s的基因组的gbk文件目录" % sample,0,'130003'],
                ["%s/data_QC" % sample, "", "样品%s的测序数据质控统计目录" % sample,0,'130004'],
                ["%s/genomic_assessment" % sample, "", "样品%s的基因组评估结果目录" % sample,0,'130005'],
                ["%s/genomic_assessment/depth_gc_.*" % sample, "", "样品%s的基因组评估gc_depth结果目录" % sample, 1],
                ["%s/genomic_assessment/kmer_frequency" % sample, "", "样品%s的基因组kmer评估目录" % sample,0,'130006'],
                ["%s/assembly_predict" % sample, "", "样品%s的基因组组装与预测结果目录" % sample,0,'130007'],
                ["%s/assembly_predict/assembly" % sample, "", "样品%s的基因组组装结果目录" % sample,0,'130008'],
                ["%s/assembly_predict/assembly/scaffold" % sample, "", "样品%s的基因组组装scaffold结果目录" % sample, 1],
                ["%s/assembly_predict/assembly/contig" % sample, "", "样品%s的基因组组装contig结果目录" % sample, 1],
                ["%s/assembly_predict/assembly/seq_dir" % sample, "", "样品%s的基因组完成图seq_dir结果目录" % sample, 1],
                ["%s/assembly_predict/predict" % sample, "", "样品%s的基因组基因预测目录" % sample,0,'130009'],
                ["%s/assembly_predict/predict/CDS_predict" % sample, "", "样品%s的编码基因预测目录" % sample,0,'130010'],
                ["%s/assembly_predict/predict/repeats" % sample, "", "样品%s的基因组重复序列预测结果目录" % sample,0,'130011'],
                ["%s/assembly_predict/predict/rRNA" % sample, "", "样品%s的rRNA预测结果目录" % sample,0,'130012'],
                ["%s/assembly_predict/predict/tRNA" % sample, "", "样品%s的tRNA预测结果目录" % sample,0,'130013'],
                ["%s/annotation" % sample, "", "样品%s的基因注释结果目录" % sample,0,'130014'],
                ["%s/annotation/NR" % sample, "", "样品%s的基因NR注释结果目录" % sample,0,'130015'],
                ["%s/annotation/Swissprot" % sample, "", "样品%s的基因Swiss-Prot注释结果目录" % sample,0,'130016'],
                ["%s/annotation/Pfam" % sample, "", "样品%s的基因Pfam注释结果目录" % sample,0,'130017'],
                ["%s/annotation/COG" % sample, "", "样品%s的基因COG注释结果目录" % sample,0,'130018'],
                ["%s/annotation/GO" % sample, "", "样品%s的基因GO注释结果目录" % sample,0,'130019'],
                ["%s/annotation/KEGG" % sample, "", "样品%s的基因KEGG注释结果目录" % sample,0,'130020'],
                ["%s/annotation/Summary" % sample, "", "样品%s的基因组注释结果汇总目录" % sample,0,'130021'],
                ["%s/structral_genome" % sample, "", "样品%s的结构基因组分析目录" % sample,0,'130022'],
                ["%s/structral_genome/promoter_predict" % sample, "", "样品%s的启动子预测结果目录" % sample,0,'130023'],
                ["%s/mobile_elements" % sample, "", "样品%s的可移动元件分析结果目录" % sample,0,'130024'],
                ["%s/mobile_elements/Genomic_Islands" % sample, "", "样品%s的基因组岛预测结果目录" % sample,0,'130025'],
                ["%s/mobile_elements/prophage" % sample, "", "样品%s的前噬菌体预测结果目录" % sample,0,'130026'],
                ["%s/mobile_elements/CRISPR_Cas" % sample, "", "样品%s的CRISPR_Cas系统预测结果目录" % sample,0,'130027'],
                ["%s/metabolic_system" % sample, "", "样品%s的代谢系统分析目录" % sample,0,'130028'],
                ["%s/metabolic_system/CAZy" % sample, "", "样品%s的碳水化合物活性酶注释结果目录" % sample,0,'130029'],
                ["%s/metabolic_system/antiSMASH" % sample, "", "样品%s的次级代谢产物合成基因簇分析结果目录" % sample,0,'130030'],
                ["%s/pathogenic_system" % sample, "", "样品%s的致病系统分析目录" % sample,0,'130031'],
                ["%s/pathogenic_system/VFDB" % sample, "", "样品%s的毒力基因预测结果目录" % sample,0,'130032'],
                ["%s/pathogenic_system/CARD" % sample, "", "样品%s的耐药基因预测结果目录" % sample,0,'130033'],
                ["%s/pathogenic_system/PHI" % sample, "", "样品%s的病原菌与宿主互作分析结果目录" % sample,0,'130034'],
                ["%s/pathogenic_system/TCDB" % sample, "", "样品%s的转运蛋白分析结果目录" % sample,0,'130035'],
                ["%s/pathogenic_system/TMHMM" % sample, "", "样品%s的跨膜蛋白分析结果目录" % sample,0,'130036'],
                ["%s/pathogenic_system/secretion_system" % sample, "", "样品%s的分泌系统分析结果目录" % sample,0,'130037'],
                ["%s/structral_genome" % sample, "", "样品%s的结构基因组分析目录" % sample, 0, '130022'],
                ["%s/circos" % sample, "", "样品%s的基因圈图分析circos目录" % sample, 1],
                ["%s/cgview" % sample, "", "样品%s的基因圈图分析cgview目录" % sample, 1],
            ]
            if self.option("analysis") in ["uncomplete"]:
                regexps += [
                    [r"%s/project_overview/gbk/.+\.gbk" % sample, "", "基因组的gbk文件",0,'130038'],
                    [r"%s/data_QC/.+_Illumina_statistics.xls" % sample, "xls", "二代测序数据质控统计表",0,'130039'],
                    [r"%s/data_QC/.+_PacBio_statistics.xls" % sample, "xls", "三代测序数据质控统计表",0,'130040'],
                    [r"%s/genomic_assessment/kmer_frequency/.+_Kmer_frequency.xls" % sample, "xls", "基因组评估kmer频率表",0,'130041'],
                    [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'130042'],
                    [r"%s/assembly_predict/assembly/.+_assembly_details.xls" % sample, "xls", "组装结果详情表",0,'130043'],
                    [r"%s/assembly_predict/assembly/.+_assembly_scaffold_details.xls" % sample, "xls","组装结果scaffolds详情表",0,'130044'],
                    [r"%s/assembly_predict/assembly/.+_assembly_contig_details.xls" % sample, "xls", "组装结果contigs详情表",0,'130045'],
                    [r"%s/assembly_predict/assembly/.+.agp" % sample, "", "基因组组装的scaffold与contig对应关系文件",0,'130046'],
                    [r"%s/assembly_predict/assembly/.+_scaf.fna" % sample, "", "基因组组装的scaffold文件",0,'130047'],
                    [r"%s/assembly_predict/assembly/.+_ctg.fna" % sample, "", "基因组组装的contigs文件",0,'130048'],
                    [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'130042'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.faa" % sample, "", "编码基因预测氨基酸序列",0,'130049'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.fnn" % sample, "", "编码基因预测核苷酸序列",0,'130050'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.gff" % sample, "", "编码基因预测gff格式统计表",0,'130051'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_CDS_statistics.xls" % sample, "xls", "编码基因预测统计表",0,'130052'],
                    [r"%s/assembly_predict/predict/repeats/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0,'130053'],
                    [r"%s/assembly_predict/predict/repeats/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0,'130054'],
                    [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0,'130055'],
                    [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0,'130056'],
                    [r"%s/assembly_predict/predict/rRNA/.+_rRNA.fnn" % sample, "", "rRNA预测序列文件",0,'130057'],
                    [r"%s/assembly_predict/predict/rRNA/.+_rRNA.gff" % sample, "", "rRNA预测gff统计文件",0,'130058'],
                    [r"%s/assembly_predict/predict/tRNA/.+_tRNA.fnn" % sample, "", "tRNA预测序列文件",0,'130059'],
                    [r"%s/assembly_predict/predict/tRNA/.+_tRNA.gff" % sample, "", "tRNA预测结果文件",0,'130060'],
                    [r"%s/assembly_predict/predict/tRNA/.+_tRNA.struc" % sample, "", "tRNA预测二级结构文件",0,'130061'],
                    [r"%s/annotation/NR/.+_anno_nr.xls" % sample, "xls", "基因NR注释结果详情表",0,'130062'],
                    [r"%s/annotation/Swissprot/.+_anno_swissprot.xls" % sample, "xls", "基因Swiss-Prot注释结果详情表",0,'130063'],
                    [r"%s/annotation/Pfam/.+_anno_pfam.xls" % sample, "xls", "基因Pfam注释结果详情表",0,'130064'],
                    [r"%s/annotation/COG/.+_cog_anno.xls" % sample, "xls", "基因COG注释结果详情表",0,'130065'],
                    [r"%s/annotation/COG/.+_cog_summary.xls" % sample, "xls", "COG注释结果汇总统计表",0,'130066'],
                    [r"%s/annotation/GO/.+_go_anno.xls" % sample, "xls", "基因GO注释结果详情表",0,'130067'],
                    [r"%s/annotation/GO/.+_go_level2_statistics.xls" % sample, "xls", "GO注释level2统计文件",0,'130068'],
                    [r"%s/annotation/GO/.+_go_list.xls" % sample, "xls", "基因与GO注释ID对应表",0,'130069'],
                    [r"%s/annotation/KEGG/.+_kegg_anno.xls" % sample, "xls", "KEGG注释结果文件",0,'130070'],
                    [r"%s/annotation/KEGG/.+_kegg_level_stat.xls" % sample, "xls", "KEGG注释level统计文件",0,'130071'],
                    [r"%s/annotation/KEGG/kegg_pathway_img.tar.gz" % sample, "", "KEGG注释pathway通路图片的压缩文件",0,'130072'],
                    [r"%s/annotation/Summary/.+_anno_summary.xls" % sample, "xls", "基因组注释结果汇总表",0,'130073'],
                    [r"%s/structral_genome/promoter_predict/.+_promoter_result.xls" % sample, "xls", "全基因组水平上启动子预测结果表",0,'130074'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_summary.xls" % sample, "xls", "全基因组水平上基因组岛预测结果文件",0,'130075'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_detail.xls" % sample, "xls", "全基因组水平上基因组岛预测结果详情表",0,'130076'],
                    [r"%s/mobile_elements/prophage/.+_prophage_summary.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果统计表",0,'130077'],
                    [r"%s/mobile_elements/prophage/.+_prophage_detail.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果详情表",0,'130078'],
                    [r"%s/mobile_elements/prophage/.+_prophage.fna" % sample, "", "全基因组水平上前噬菌体预测结果序列文件",0,'130079'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_summary.xls" % sample, "xls", "CRISPR_Cas系统预测结果统计表",0,'130080'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_detail.xls" % sample, "xls", "CRISPR_Cas系统预测结果详情表",0,'130081'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_anno.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释结果表",0,'130082'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_family_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释family统计表",0,'130083'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_class_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释class统计表",0,'130084'],
                    [r"%s/metabolic_system/antiSMASH/.+_antismash_anno.xls" % sample, "xls", "次级代谢产物合成基因簇预测结果",0,'130085'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_align.xls" % sample, "xls", "全基因组水平上毒力基因预测表",0,'130086'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_anno.xls" % sample, "xls", "全基因组水平上毒力基因注释表",0,'130087'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_level.xls" % sample, "xls", "全基因组水平上毒力基因分级统计表",0,'130088'],
                    [r"%s/pathogenic_system/CARD/.+_card_align.xls" % sample, "xls", "全基因组水平上耐药基因预测表",0,'130089'],
                    [r"%s/pathogenic_system/CARD/.+_card_anno.xls" % sample, "xls", "全基因组水平上耐药基因注释表",0,'130090'],
                    [r"%s/pathogenic_system/CARD/.+_card_category.xls" % sample, "xls", "全基因组水平上耐药基因分类统计表",0,'130091'],
                    [r"%s/pathogenic_system/PHI/.+_phi_align.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因预测表",0,'130092'],
                    [r"%s/pathogenic_system/PHI/.+_phi_anno.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因注释表",0,'130093'],
                    [r"%s/pathogenic_system/PHI/.+_phi_stat.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析PHI ID统计表",0,'130094'],
                    [r"%s/pathogenic_system/PHI/.+_phi_phenotype.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析表型统计表",0,'130095'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align.xls" % sample, "xls", "转运蛋白预测表",0,'130096'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align_top1.xls" % sample, "xls", "转运蛋白预测表(TOP1)",0,'130097'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_anno.xls" % sample, "xls", "转运蛋白注释表",0,'130098'],
                    [r"%s/pathogenic_system/TCDB/.+_whole_genome_tcdb_anno.xls" % sample, "xls", "全基因组水平上转运蛋白注释表",0,'130099'],
                    [r"%s/pathogenic_system/TMHMM/.+_tmhmm_anno.xls" % sample, "xls", "跨膜蛋白分析结果详情表",0,'130100'],
                    [r"%s/pathogenic_system/TMHMM/.+_whole_genome_tmhmm_anno.xls" % sample, "xls", "全基因组水平上跨膜蛋白注释表",0,'130101'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_genes.xls" % sample, "xls","分泌系统相关基因统计表",0,'130102'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_type.xls" % sample, "xls","分泌系统的类型统计表",0,'130103'],
                ]
            elif self.option("analysis") in ["complete"]:
                regexps += [
                    [r"%s/project_overview/gbk/.+_chromosome[0-9]\.gbk" % sample, "", "基因组染色体的gbk文件",0,'130104'],
                    [r"%s/project_overview/gbk/.+_plasmid[A-Z]\.gbk" % sample, "", "基因组质粒的gbk文件",0,'130105'],
                    [r"%s/data_QC/.+_Illumina_statistics.xls" % sample, "xls", "二代测序数据质控统计表",0,'130039'],
                    [r"%s/data_QC/.+_PacBio_statistics.xls" % sample, "xls", "三代测序数据质控统计表",0,'130040'],
                    [r"%s/genomic_assessment/kmer_frequency/.+_Kmer_frequency.xls" % sample, "xls", "基因组评估kmer频率表",0,'130041'],
                    [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'130042'],
                    [r"%s/assembly_predict/assembly/.+_assembly_details.xls" % sample, "xls", "组装结果详情表",0,'130043'],
                    [r"%s/assembly_predict/assembly/.+_chromosome[0-9]\.fna" % sample, "", "基因组完成图组装的染色体文件",0,'130106'],
                    [r"%s/assembly_predict/assembly/.+_plasmid[A-Z]\.fna" % sample, "", "基因组完成图组装的质粒文件",0,'130107'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS_statistics.xls" % sample, "xls","全基因组水平编码基因预测统计表",0,'130108'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS_statistics.xls" % sample, "xls","核染色体水平编码基因预测统计表",0,'130109'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS_statistics.xls" % sample, "xls","质粒水平编码基因预测统计表",0,'130110'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS.gff" % sample, "","全基因组水平编码基因预测gff格式统计表",0,'130111'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS.gff" % sample, "","核染色体水平编码基因预测gff格式统计表",0,'130112'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS.gff" % sample, "","质粒水平编码基因预测gff格式统计表",0,'130113'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS.faa" % sample, "","全基因组水平编码基因预测氨基酸序列",0,'130114'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS.faa" % sample, "","核染色体水平编码基因预测氨基酸序列",0,'130115'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS.faa" % sample, "","质粒水平编码基因预测氨基酸序列",0,'130116'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS.fnn" % sample, "","全基因组水平编码基因预测核苷酸序列",0,'130117'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS.fnn" % sample, "","染色体水平编码基因预测核苷酸序列",0,'130118'],
                    [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS.fnn" % sample, "","质粒水平编码基因预测核苷酸序列",0,'130119'],
                    [r"%s/assembly_predict/predict/repeats/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0,'130053'],
                    [r"%s/assembly_predict/predict/repeats/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0,'130054'],
                    [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0,'130055'],
                    [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0,'130056'],
                    [r"%s/assembly_predict/predict/rRNA/.+_rRNA.fnn" % sample, "", "rRNA预测序列文件",0,'130057'],
                    [r"%s/assembly_predict/predict/rRNA/.+_rRNA.gff" % sample, "", "rRNA预测gff统计文件",0,'130058'],
                    [r"%s/assembly_predict/predict/tRNA/.+_tRNA.fnn" % sample, "", "tRNA预测序列文件",0,'130059'],
                    [r"%s/assembly_predict/predict/tRNA/.+_tRNA.gff" % sample, "", "tRNA预测结果文件",0,'130060'],
                    [r"%s/assembly_predict/predict/tRNA/.+_tRNA.struc" % sample, "", "tRNA预测二级结构文件",0,'130061'],
                    [r"%s/annotation/NR/.+_whole_genome_anno_nr.xls" % sample, "xls", "全基因组水平上基因NR注释结果详情表",0,'130120'],
                    [r"%s/annotation/NR/.+_chromosome[0-9]_anno_nr.xls" % sample, "xls", "核染色体水平上基因NR注释结果详情表",0,'130121'],
                    [r"%s/annotation/NR/.+_plasmid[A-Z]_anno_nr.xls" % sample, "xls", "质粒水平上基因NR注释结果详情表",0,'130122'],
                    [r"%s/annotation/Swissprot/.+_chromosome[0-9]_anno_swissprot.xls" % sample, "xls","核染色体水平上基因Swiss-Prot注释结果详情表",0,'130123'],
                    [r"%s/annotation/Swissprot/.+_plasmid[A-Z]_anno_swissprot.xls" % sample, "xls","质粒水平上基因Swiss-Prot注释结果详情表",0,'130124'],
                    [r"%s/annotation/Pfam/.+_whole_genome_anno_pfam.xls" % sample, "xls", "全基因组水平上基因Pfam注释结果详情表",0,'130125'],
                    [r"%s/annotation/Pfam/.+_chromosome[0-9]_anno_pfam.xls" % sample, "xls", "核染色体水平上基因Pfam注释结果详情表",0,'130126'],
                    [r"%s/annotation/Pfam/.+_plasmid[A-Z]_anno_pfam.xls" % sample, "xls", "质粒水平上基因Pfam注释结果详情表",0,'130127'],
                    [r"%s/annotation/COG/.+_cog_anno.xls" % sample, "xls", "基因COG注释结果详情表",0,'130065'],
                    [r"%s/annotation/COG/.+_cog_summary.xls" % sample, "xls", "COG注释结果汇总统计表",0,'130066'],
                    [r"%s/annotation/GO/.+_go_anno.xls" % sample, "xls", "基因GO注释结果详情表",0,'130067'],
                    [r"%s/annotation/GO/.+_go_level2_statistics.xls" % sample, "xls", "GO注释level2统计文件",0,'130068'],
                    [r"%s/annotation/GO/.+_go_list.xls" % sample, "xls", "基因与GO注释ID对应表",0,'130069'],
                    [r"%s/annotation/KEGG/.+_kegg_anno.xls" % sample, "xls", "KEGG注释结果文件",0,'130070'],
                    [r"%s/annotation/KEGG/.+_kegg_level_stat.xls" % sample, "xls", "KEGG注释level统计文件",0,'130071'],
                    [r"%s/annotation/KEGG/kegg_pathway_img.tar.gz" % sample, "", "KEGG注释pathway通路图片的压缩文件",0,'130072'],
                    [r"%s/annotation/Summary/.+_anno_summary.xls" % sample, "xls", "基因组注释结果汇总表",0,'130073'],
                    [r"%s/structral_genome/promoter_predict/.+_promoter_result.xls" % sample, "xls", "全基因组水平上启动子预测结果表",0,'130074'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_summary.xls" % sample, "xls", "全基因组水平上基因组岛预测结果文件",0,'130075'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_detail.xls" % sample, "xls", "全基因组水平上基因组岛预测结果详情表",0,'130076'],
                    [r"%s/mobile_elements/prophage/.+_prophage_summary.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果统计表",0,'130077'],
                    [r"%s/mobile_elements/prophage/.+_prophage_detail.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果详情表",0,'130078'],
                    [r"%s/mobile_elements/prophage/.+_prophage.fna" % sample, "", "全基因组水平上前噬菌体预测结果序列文件",0,'130079'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_summary.xls" % sample, "xls", "CRISPR_Cas系统预测结果统计表",0,'130080'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_detail.xls" % sample, "xls", "CRISPR_Cas系统预测结果详情表",0,'130081'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_anno.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释结果表",0,'130082'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_family_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释family统计表",0,'130083'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_class_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释class统计表",0,'130084'],
                    [r"%s/metabolic_system/antiSMASH/.+_antismash_anno.xls" % sample, "xls", "次级代谢产物合成基因簇预测结果",0,'130085'],
                    [r"%s/metabolic_system/antiSMASH/.+_chromosome[0-9]_antismash_anno.xls" % sample, "xls","核染色体水平上次级代谢产物合成基因簇预测结果",0,'130128'],
                    [r"%s/metabolic_system/antiSMASH/.+_plasmid[A-Z]_antismash_anno.xls" % sample, "xls","质粒水平上次级代谢产物合成基因簇预测结果",0,'130129'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_align.xls" % sample, "xls", "全基因组水平上毒力基因预测表",0,'130086'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_anno.xls" % sample, "xls", "全基因组水平上毒力基因注释表",0,'130087'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_level.xls" % sample, "xls", "全基因组水平上毒力基因分级统计表",0,'130088'],
                    [r"%s/pathogenic_system/CARD/.+_card_align.xls" % sample, "xls", "全基因组水平上耐药基因预测表",0,'130089'],
                    [r"%s/pathogenic_system/CARD/.+_card_anno.xls" % sample, "xls", "全基因组水平上耐药基因注释表",0,'130090'],
                    [r"%s/pathogenic_system/CARD/.+_card_category.xls" % sample, "xls", "全基因组水平上耐药基因分类统计表",0,'130091'],
                    [r"%s/pathogenic_system/PHI/.+_phi_align.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因预测表",0,'130092'],
                    [r"%s/pathogenic_system/PHI/.+_phi_anno.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因注释表",0,'130093'],
                    [r"%s/pathogenic_system/PHI/.+_phi_stat.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析PHI ID统计表",0,'130094'],
                    [r"%s/pathogenic_system/PHI/.+_phi_phenotype.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析表型统计表",0,'130095'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align.xls" % sample, "xls", "转运蛋白预测表",0,'130096'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align_top1.xls" % sample, "xls", "转运蛋白预测表(TOP1)",0,'130097'],
                    [r"%s/pathogenic_system/TCDB/.+_whole_genome_tcdb_anno.xls" % sample, "xls", "全基因组水平上转运蛋白注释表",0,'130099'],
                    [r"%s/pathogenic_system/TCDB/.+_chromosome[0-9]_tcdb_anno.xls" % sample, "xls", "核染色体水平上转运蛋白注释表",0,'130130'],
                    [r"%s/pathogenic_system/TCDB/.+_plasmid[A-Z]_tcdb_anno.xls" % sample, "xls", "质粒水平上转运蛋白注释表",0,'130131'],
                    [r"%s/pathogenic_system/TMHMM/.+_whole_genome_tmhmm_anno.xls" % sample, "xls", "全基因组水平上跨膜蛋白注释表",0,'130101'],
                    [r"%s/pathogenic_system/TMHMM/.+_chromosome[0-9]_tmhmm_anno.xls" % sample, "xls", "核染色体水平上跨膜蛋白注释表",0,'130132'],
                    [r"%s/pathogenic_system/TMHMM/.+_plasmid[A-Z]_tmhmm_anno.xls" % sample, "xls", "质粒水平上跨膜蛋白注释表",0,'130133'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_genes.xls" % sample, "xls","分泌系统相关基因统计表",0,'130102'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_type.xls" % sample, "xls","分泌系统的类型统计表",0,'130103'],
                ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

    def get_lib_type(self,file):
        lib_type = {}
        with open(file,'r') as f:
            lines =f.readlines()
            for line in lines[1:]:
                lin =line.rstrip('\r\n').split('\t')
                if re.search(r'PE',lin[5]):
                    des = lin[5] + str(lin[2])
                    if lin[0] in lib_type.keys():
                        lib_type[line[0]][des] = lin[5]
                    else:
                        lib_type[lin[0]] = {des: lin[5]}

        return lib_type

    def get_assemble_type(self,file):
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



    def get_list(self):
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                if os.path.exists(list_path):
                    self.logger.info(list_path)
                with open(list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 6:
                            self.samples[line[0]] = line[0]
                        else:
                            #raise Exception('raw_dir的list.txt文件格式有误')
                            self.set_error('raw_dir的list.txt文件格式有误', code="11400101")
            if self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                if os.path.exists(list_path):
                    self.logger.info(list_path)
                with open(list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 2:
                            self.samples[line[0]] = line[0]
                        else:
                            #raise Exception('asse_dir的list.txt文件格式有误')
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400102")
            if self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                raw_sample ={}
                assemble_sample ={}
                with open(raw_list_path, "rb") as f:
                    lines=f.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 6:
                            raw_sample[line[0]] = line[0]
                        else:
                            #raise Exception('raw_dir的list.txt文件格式有误')
                            self.set_error('raw_dir的list.txt文件格式有误', code="11400103")
                with open(assemble_list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 2:
                            assemble_sample[line[0]] = line[0]
                        else:
                            #raise Exception('asse_dir的list.txt文件格式有误')
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400104")
                for key in raw_sample.keys():
                    if key in assemble_sample.keys():
                        self.samples[key]=key
                    else:
                        #raise Exception('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致')
                        self.set_error('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致', code="11400105")

        elif self.option("analysis") in ["complete"]:
            if self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                if os.path.exists(list_path):
                    self.logger.info(list_path)
                with open(list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 3:
                            self.samples[line[0]] = line[0]
                        else:
                            #raise Exception('asse_dir的list.txt文件格式有误')
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400106")
                        if line[2] not in ['chromosome', 'plasmid']:
                            #raise Exception('asse_dir的list.txt文件的Genome Type有误！')
                            self.set_error('asse_dir的list.txt文件的Genome Type有误！', code="11400107")

            if self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                raw_sample ={}
                assemble_sample ={}
                with open(raw_list_path, "rb") as f:
                    lines=f.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 6:
                            raw_sample[line[0]] = line[0]
                        else:
                            #raise Exception('raw_dir的list.txt文件格式有误')
                            self.set_error('raw_dir的list.txt文件格式有误', code="11400108")
                with open(assemble_list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 3:
                            assemble_sample[line[0]] = line[0]
                        else:
                            #raise Exception('asse_dir的list.txt文件格式有误')
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400109")
                        if line[2] not in ['chromosome', 'plasmid']:
                            #raise Exception('asse_dir的list.txt文件的Genome Type有误！')
                            self.set_error('asse_dir的list.txt文件的Genome Type有误！', code="11400110")
                for key in raw_sample.keys():
                    if key in assemble_sample.keys():
                        self.samples[key]=key
                    else:
                        #raise Exception('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致')
                        self.set_error('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致', code="11400111")

    def run_split_dir(self):
        self.get_list()
        samples = self.samples
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path = os.path.join(self.work_dir,sample)
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    reslut_path = os.path.join(reslut_path, 'data')
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    output = reslut_path + "/list.txt"
                    file = open(output, 'w')
                    with open(list_path, "rb") as l:
                        lines = l.readlines()
                        file.write(lines[0])
                        for line in lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file.write(line )
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files_path = raw.split(',')
                                        for file_path in files_path:
                                            if os.path.exists(reslut_path + "/" + file_path):
                                                os.remove(reslut_path + "/" + file_path)
                                            os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                         reslut_path + "/" + file_path)
                                else:
                                    files_path = line2[1].split(',')
                                    for file_path in files_path:
                                        if os.path.exists(reslut_path + "/" + file_path):
                                            os.remove(reslut_path + "/" + file_path)
                                        os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                     reslut_path + "/" + file_path)

            if self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    reslut_path = os.path.join(reslut_path, 'assemble')
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    output = reslut_path + "/list.txt"
                    file = open(output, 'w')
                    with open(list_path, "rb") as l:
                        lines = l.readlines()
                        file.write(lines[0])
                        for line in lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file.write(line)
                                if os.path.exists(reslut_path + "/" + line2[1]):
                                    os.remove(reslut_path + "/" + line2[1])
                                os.link(self.option("asse_dir").prop['path'] + "/" + line2[1],
                                             reslut_path + "/" + line2[1])

            if self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path1 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path1 = os.path.join(reslut_path1, 'data')
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path2 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    reslut_path2 = os.path.join(reslut_path2, 'assemble')
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    output1 = reslut_path1 + "/list.txt"
                    file1 = open(output1, 'w')
                    output2 = reslut_path2 + "/list.txt"
                    file2 = open(output2, 'w')
                    with open(raw_list_path, "rb") as l, open(assemble_list_path, "rb") as f:
                        raw_lines = l.readlines()
                        ass_lines = f.readlines()
                        file1.write(raw_lines[0] )
                        file2.write(ass_lines[0])
                        for line in raw_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file1.write(line)
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files_path = raw.split(',')
                                        for file_path in files_path:
                                            if os.path.exists(reslut_path1 + "/" + file_path):
                                                os.remove(reslut_path1 + "/" + file_path)
                                            os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                         reslut_path1 + "/" + file_path)
                                else:
                                    files_path = line2[1].split(',')
                                    for file_path in files_path:
                                        if os.path.exists(reslut_path1 + "/" + file_path):
                                            os.remove(reslut_path1 + "/" + file_path)
                                        os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                     reslut_path1 + "/" + file_path)
                        for line in ass_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file2.write(line)
                                if os.path.exists(reslut_path2 + "/" + line2[1]):
                                    os.remove(reslut_path2 + "/" + line2[1])
                                os.link(self.option("asse_dir").prop['path'] + "/" + line2[1],
                                             reslut_path2 + "/" + line2[1])

        elif self.option("analysis") in ["complete"]:
            if self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    reslut_path = os.path.join(reslut_path, 'assemble')
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    output = reslut_path + "/list.txt"
                    file = open(output, 'w')
                    with open(list_path, "rb") as l:
                        lines = l.readlines()
                        file.write(lines[0])
                        for line in lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file.write(line)
                                if os.path.exists(reslut_path + "/" + line2[1]):
                                    os.remove(reslut_path + "/" + line2[1])
                                os.link(self.option("asse_dir").prop['path'] + "/" + line2[1],reslut_path + "/" + line2[1])

            if self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path1 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path1 = os.path.join(reslut_path1, 'data')
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path2 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    reslut_path2 = os.path.join(reslut_path2, 'assemble')
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    output1 = reslut_path1 + "/list.txt"
                    file1 = open(output1, 'w')
                    output2 = reslut_path2 + "/list.txt"
                    file2 = open(output2, 'w')
                    with open(raw_list_path, "rb") as l, open(assemble_list_path, "rb") as f:
                        raw_lines = l.readlines()
                        ass_lines = f.readlines()
                        file1.write(raw_lines[0])
                        file2.write(ass_lines[0])
                        for line in raw_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file1.write(line)
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files_path = raw.split(',')
                                        for file_path in files_path:
                                            if os.path.exists(reslut_path1 + "/" + file_path):
                                                os.remove(reslut_path1 + "/" + file_path)
                                            os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                        reslut_path1 + "/" + file_path)
                                else:
                                    files_path = line2[1].split(',')
                                    for file_path in files_path:
                                        if os.path.exists(reslut_path1 + "/" + file_path):
                                            os.remove(reslut_path1 + "/" + file_path)
                                        os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                    reslut_path1 + "/" + file_path)
                        for line in ass_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file2.write(line)
                                if os.path.exists(reslut_path2 + "/" + line2[1]):
                                    os.remove(reslut_path2 + "/" + line2[1])
                                os.link(self.option("asse_dir").prop['path'] + "/" + line2[1],
                                            reslut_path2 + "/" + line2[1])

    def get_seq_type(self):
        dict = {}
        if self.option("analysis") in ["complete"]:
            ass_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
            for sample in self.samples.keys():
                list = []
                chr_num = 0
                pla_num = 0
                seq_type = ''
                with open(ass_list_path, 'r') as f:
                    lines = f.readlines()
                    for line in lines[1:]:
                        line2 = line.rstrip('\r\n').split('\t')
                        if sample == line2[0]:
                            if line2[2] == 'chromosome' or line2[2] == 'Chromosome':
                               chr_num +=1
                            elif line2[2] == 'plasmid' or line2[2] == 'Plasmid':
                               pla_num +=1
                if chr_num ==1:
                    list.append('Chromosome')
                if chr_num >1:
                    for i in range(1,chr_num +1):
                        list.append('Chromosome' + str(i))
                if pla_num ==1:
                    list.append('Plasmid')
                if pla_num >1:
                    for i in range(1,pla_num +1):
                        if i ==1:
                            list.append('PlasmidA')
                        if i == 2:
                            list.append('PlasmidB')
                        if i ==3:
                            list.append('PlasmidC')
                        if i ==4:
                            list.append('PlasmidD')
                        if i ==5:
                            list.append('PlasmidE')
                        if i ==6:
                            list.append('PlasmidF')
                        if i ==7:
                            list.append('PlasmidG')
                if len(list) ==1:
                    seq_type =list[0]
                else:
                    seq_type =','.join(list)
                dict[sample] = seq_type
        return dict

    def move_dir(self, olddir, newdir):  # 原函数名move2outputdir
        """
        移动一个目录下所有文件/文件夹到workflow输出路径下，供set_output调用
        """
        start = time.time()
        if not os.path.isdir(olddir):
            #raise Exception('需要移动到output目录的文件夹不存在。')
            self.set_error('需要移动到output目录的文件夹不存在。', code="11400112")
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
