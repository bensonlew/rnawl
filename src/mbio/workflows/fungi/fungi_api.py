# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""真菌基因组分析工作流"""

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

class FungiApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        真菌基因组workflow option参数设置
        """
        self._sheet = wsheet_object
        super(FungiApiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "raw_dir", "type": "infile", "format": "fungi_genome.raw_dir"},  ###rawdata的文件目录
            {"name": "asse_dir", "type": "infile", "format": "fungi_genome.asse_dir"},  ###组装的文件目录
            {"name": "analysis", "type": "string", "default": "uncomplete"}, ###流程分析模式complete，uncomplete
            {"name": "data_type", "type": "string"},
            {"name": "faa_file", "type": "infile","format": "sequence.fasta"},#训练集蛋白文件
            {"name": "busco_database", "type": "string"},#busco数据库参数
            {"name": "ncbi", "type": "string"},  #
            {"name": "genbank", "type": "string"},  #
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        '''初始化module/tool'''
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.samples ={}
        self.modules = []
        self.logger.info(self._sheet.output)
        self.remote_dir = self._sheet.output + '/'
        

    def check_options(self):
        """
        检查参数
        """
        if not self.option("analysis"):
            raise OptionError("请提供流程分析模式！", code="12100101")
        if not self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
            raise OptionError("必须输入原始序列文件夹或组装序列文件夹其中一个！", code="12100102")


    def run(self):
        """
        运行:genome_workflow
        :return:
        """
        self.get_list()
        self.run_api()
        gevent.spawn_later(5, self.end)
        super(FungiApiWorkflow, self).run()

    def wait_file(self, path, wait_times=1):
        '''
        增加等待文件结果方法
        :param path: 结果文件路径
        :param wait_times: 等待次数
        :return: 文件路径
        :time: 20180425
        '''
        while wait_times < 11:
            if not os.path.exists(path):
                time.sleep(10)
                wait_times += 1
                self.wait_file(path, wait_times=wait_times)
            return path
        self.logger.info("超过文件等待次数，需检查文件%s" % path)
        return

    def end(self):
        self.run_api()
        self.send_files()
        super(FungiApiWorkflow, self).end()

    def run_api(self, test=False):
        task_id = self._sheet.id
        self.datastat = self.api.api("fungi_genome.genome_qc")
        self.genome_size = self.api.api("fungi_genome.genome_size")
        self.assess_gc = self.api.api("fungi_genome.assess_gc")
        self.assemble = self.api.api("fungi_genome.assemble")
        self.assess_kmer = self.api.api("fungi_genome.assess_kmer")
        self.anno_card = self.api.api("fungi_genome.anno_card")
        self.anno_cazy = self.api.api("fungi_genome.anno_cazy")
        self.anno_cog = self.api.api("fungi_genome.anno_cog")
        self.anno_go = self.api.api("fungi_genome.anno_go")
        self.anno_kegg = self.api.api("fungi_genome.anno_kegg")
        self.anno_nr = self.api.api("fungi_genome.anno_nr")
        self.anno_pfam = self.api.api("fungi_genome.anno_pfam")
        self.anno_phi = self.api.api("fungi_genome.anno_phi")
        self.anno_summary = self.api.api("fungi_genome.anno_summary")
        self.anno_swissprot = self.api.api("fungi_genome.anno_swissprot")
        self.anno_cyps = self.api.api("fungi_genome.anno_cyps")
        self.anno_dfvf = self.api.api("fungi_genome.anno_dfvf")
        self.anno_signalp = self.api.api("fungi_genome.anno_signalp")
        self.anno_tcdb = self.api.api("fungi_genome.anno_tcdb")
        self.anno_tmhmm = self.api.api("fungi_genome.anno_tmhmm")
        self.gene_predict = self.api.api("fungi_genome.gene_predict")
        self.promote = self.api.api("fungi_genome.promote")
        self.repeat_predict = self.api.api("fungi_genome.repeat_predict")
        self.rrna_predict = self.api.api("fungi_genome.rrna_predict")
        self.secretory = self.api.api("fungi_genome.secretory")
        self.summary_map = self.api.api("fungi_genome.summary_map")
        self.trna_predict = self.api.api("fungi_genome.trna_predict")
        seq_type = ''
        if self.option('raw_dir').is_set:
            seq_type = self.get_assemble_type(self.option('raw_dir').prop['path'] + '/list.txt')
        datastat_id = ''
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and seq_type in ['PE', 'PE,MP', 'MP,PE']:
                datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                         "Trimmomatic, SeqPrep, Sickle, FastqTotalHighQualityBase.jar",
                                                         "30", "20")
            elif self.option("raw_dir").is_set and seq_type in ['pacbio', 'Pacbio', 'PACBIO']:
                datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                         "bamtools", "", "")
            elif self.option("raw_dir").is_set and seq_type in ['PE,pacbio', 'pacbio,PE', 'PE,Pacbio', 'Pacbio,PE']:
                datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                         "Trimmomatic, SeqPrep, Sickle, FastqTotalHighQualityBase.jar;bamtools",
                                                         "30", "20")
            elif self.option("asse_dir").is_set:
                datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', "", "assemble", "", "",
                                                         "")
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not re.search(r'pacbio', seq_type):
                assemble_id = self.assemble.add_assemble(self.option('analysis'), seq_type, 'rawdata', '组装评估',
                                                         'SOAPdenovo v2.04,GapCloser v1.12', 'kmer (21-41)')
            elif self.option("raw_dir").is_set and re.search(r'pacbio', seq_type):
                assemble_id = self.assemble.add_assemble(self.option('analysis'), seq_type, 'rawdata', '组装评估',
                                                         'canu V1.7', '')
            elif self.option("asse_dir").is_set:
                assemble_id = self.assemble.add_assemble(self.option('analysis'), '', 'assemble', '组装评估', '', '')
        if self.option("analysis") in ["uncomplete"]:
            gene_id = self.gene_predict.add_gene_predict(self.remote_dir, '/assembly_predict/predict/CDS_predict/',
                                                         '_CDS.', params='{"soft": "Maker2"}')
        trna_id = self.trna_predict.add_trna_predict(self.remote_dir, '/assembly_predict/predict/tRNA/', '_tRNA.',
                                                     params='soft:tRNAscan-SE')
        rrna_id = self.rrna_predict.add_rrna_predict(self.remote_dir, '/assembly_predict/predict/rRNA/', '_rRNA.',
                                                     params='soft:barrnap')
        repeat_id = self.repeat_predict.add_repeat_predict(params='soft:TRF')
        nr_id = self.anno_nr.add_anno_nr(params='{"NR":"Diamond"}')
        cog_id = self.anno_cog.add_anno_cog(params="{COG:Diamond}")
        go_id = self.anno_go.add_anno_go(params="{GO:blast2go}")
        cazy_id = self.anno_cazy.add_anno_cazy(params="{CAZy:hmmscan}")
        resu = self.remote_dir
        kegg_id = self.anno_kegg.add_anno_kegg(resu, "/annotation/KEGG/", "_kegg_pathway_img", params="{KEGG:Diamond}")
        pfam_id = self.anno_pfam.add_anno_pfam(params="{Pfam:hmmer3}")
        swissport_id = self.anno_swissprot.add_anno_swissprot(params="{Swiss-prot:blast}")
        summary_id = self.anno_summary.add_anno_summary(params="{NR,KEGG,COG,GO,Swissprot,Pfam}")
        gc_id = self.assess_gc.add_assess_gc("GC_Depth")
        size_id = self.genome_size.add_assess_size("genome size")
        n = 1
        for sample in self.samples.keys():
            if self.option("analysis") in ["uncomplete"]:
                self.datastat.add_datastat_specimen(datastat_id,
                                                    self.output_dir + '/' + sample + '/' + 'specimen.data')
                self.logger.info('specimen.data')
                self.datastat.add_datastat_uncomplete_gene(datastat_id,
                                                           self.output_dir + '/' + sample + '/' + 'gene.data')
                self.logger.info('gene.data')
            if self.option("analysis") in ["uncomplete"]:
                self.datastat.add_datastat_uncomplete_summary(datastat_id,
                                                              self.output_dir + '/' + sample + '/project_overview/genome_overview.xls')
                self.logger.info('genome_overview.xls')
            if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                my_seq = self.get_lib_type(self.option('raw_dir').prop['path'] + '/list.txt')
                self.logger.info(my_seq)
                if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls'):
                    self.datastat.add_qc_stat_uncomplete(datastat_id,
                                                         self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls')
                    self.logger.info('Illumina_statistics.xls')
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

                if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '_PacBio_statistics.xls'):
                    self.datastat.add_qc_stat_complete(datastat_id,
                                                            self.output_dir + '/' + sample + '/data_QC/' + sample + '_PacBio_statistics.xls',
                                                            'pacbio',sample)
                    self.logger.info('PacBio_statistics.xls')
                    self.datastat.add_datastat_pacbio_graphic(datastat_id,
                                                              self.output_dir + '/' + sample + '/fastx/' + sample + '.clean.len.xls',sample,'', "len")
                    self.logger.info('pacbio.len.xls')
            if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "10k",
                                                    "/fungi/" + task_id + "/data/gc_depth/" + sample + "/depth_gc_10000")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "20k",
                                                    "/fungi/" + task_id + "/data/gc_depth/" + sample + "/depth_gc_20000")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "30k",
                                                    "/fungi/" + task_id + "/data/gc_depth/" + sample + "/depth_gc_30000")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "50k",
                                                    "/fungi/" + task_id + "/data/gc_depth/" + sample + "/depth_gc_50000")
                self.assess_gc.add_assess_gc_detail(gc_id, sample, "100k",
                                                    "/fungi/" + task_id + "/data/gc_depth/" + sample + "/depth_gc_100000")
                kmer_id = self.assess_kmer.add_assess_kmer(sample, "kmer")
                if os.path.exists(self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + ".frequency.xls"):
                    self.assess_kmer.add_assess_kmer_detail(kmer_id,
                                                        self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + ".frequency.xls")
                if os.path.exists(self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + "_Kmer_frequency.xls"):
                    self.genome_size.add_heter(size_id,
                                           self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + "_Kmer_frequency.xls",sample)
                if os.path.exists(self.output_dir + '/' + sample + "/genomic_assessment/genome_size/" + sample + ".summary.xls"):
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
                                                   sample, type,
                                                   '/fungi/' + task_id + '/data/assemble/' + sample + '/' + type)
                    self.assemble.add_assemble_graphic(assemble_id,
                                                       self.output_dir + '/' + sample + '/assembly_predict/assembly/len/' + sample + '.' + type + 's.len.xls',
                                                       sample, type)
                self.assemble.add_assemble_assess_busco(assemble_id,
                                                        self.output_dir + '/' + sample + '/assembly_predict/assembly_assessment/' + sample + '_busco.xls',
                                                        sample)
                self.assemble.add_assemble_assess_cegma(assemble_id,
                                                        self.output_dir + '/' + sample + '/assembly_predict/assembly_assessment/' + sample + '_cegma.xls',
                                                        sample)

            if self.option("analysis") in ["uncomplete"]:
                self.gene_predict.add_gene_predict_detail(gene_id, sample,
                                                          self.output_dir + '/' + sample + '/assembly_predict/predict/predict/CDS_predict/' + sample + "_CDS.gff")
                self.gene_predict.add_gene_predict_specimen(gene_id, sample,
                                                            self.output_dir + '/' + sample + '/assembly_predict/predict/predict/CDS_predict/' + sample + "_CDS_statistics.xls")
                self.gene_predict.add_gene_predict_seq(gene_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/predict/CDS_predict/' + sample + "_CDS.fnn",
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/predict/CDS_predict/' + sample + "_CDS.faa")
                self.gene_predict.add_gene_predict_bar(gene_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/predict/CDS_predict/' + sample + "_CDS_length.xls")
            self.trna_predict.add_trna_predict_detail(trna_id, sample,
                                                      self.output_dir + '/' + sample + '/assembly_predict/predict/predict/tRNA/' + sample + "_tRNA.gff")
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/assembly_predict/predict/predict/tRNA/' + sample + "_tRNA.fnn"):
                self.trna_predict.add_trna_predict_seq(trna_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/predict/tRNA/' + sample + "_tRNA.fnn",
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/predict/tRNA/' + sample + "_tRNA.struc")
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/assembly_predict/predict/predict/rRNA/' + sample + "_rRNA.fnn"):
                self.rrna_predict.add_rrna_predict_seq(rrna_id, sample,
                                                       self.output_dir + '/' + sample + '/assembly_predict/predict/predict/rRNA/' + sample + "_rRNA.fnn")
            self.rrna_predict.add_rrna_predict_detail(rrna_id, sample,
                                                      self.output_dir + '/' + sample + '/assembly_predict/predict/predict/rRNA/' + sample + "_rRNA.gff")
            if self.option("analysis") in ["uncomplete"]:
                self.repeat_predict.add_repeat_predict_detail(repeat_id, sample,
                                                              self.output_dir + '/' + sample + '/assembly_predict/predict/repeats/' + sample + "_Rep.gff",
                                                              self.output_dir + '/' + sample + '/assembly_predict/predict/repeats/' + sample + "_Rep.tbl")
            if self.option("analysis") in ["uncomplete"]:
                self.anno_nr.add_anno_nr_detail(nr_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/NR/" + sample + "_anno_nr.xls")
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
            if self.option("analysis") in ["uncomplete"]:
                self.anno_swissprot.add_anno_swissprot_detail(swissport_id, sample,
                                                              self.output_dir + '/' + sample + "/annotation/Swissprot/" + sample + "_anno_swissprot.xls")
            self.anno_summary.add_anno_summary_detail(summary_id, sample,
                                                      self.output_dir + '/' + sample + "/annotation/Summary/" + sample + "_anno_summary.xls")
            if n == 1:
                prom_id = self.promote.add_promote(
                    self.output_dir + '/' + sample + "/structral_genome/promoter_predict",
                    main=True, specimen_id=sample, update_id=summary_id)
                self.logger.info("promoter end")
                self.gene_graph = self.api.api("fungi_genome.gene_graph")
                if self.option("analysis") in ["uncomplete"]:
                    gene_graph_id = self.gene_graph.add_gene_graph(
                        self.output_dir + '/' + sample + "/assembly_predict/predict/predict/CDS_predict/" + sample + "_CDS.gff",
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls", main=True,
                        specimen_id=sample)
                    self.logger.info("gene_graph end")
                dfvf_id = self.anno_dfvf.add_anno_dfvf(
                    self.output_dir + '/' + sample + "/pathogenic_system/DFVF", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("dfvf end")
                signalp_id = self.anno_signalp.add_anno_signalp(
                    self.output_dir + '/' + sample + "/pathogenic_system/SIGNALP", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("signalp end")
                cyps_id = self.anno_cyps.add_anno_cyps(
                    self.output_dir + '/' + sample + "/metabolic_system/P450", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("cyps end")
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
                n += 1
            else:
                self.promote.add_promote(
                    self.output_dir + '/' + sample + "/structral_genome/promoter_predict",
                    main_id=prom_id, specimen_id=sample, update_id=summary_id)
                self.logger.info("promoter end")
                if self.option("analysis") in ["uncomplete"]:
                    self.gene_graph.add_gene_graph(
                        self.output_dir + '/' + sample + "/assembly_predict/predict/predict/CDS_predict/" + sample + "_CDS.gff",
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls",
                        main_id=gene_graph_id, specimen_id=sample)
                    self.logger.info("gene_graph end")
                self.anno_dfvf.add_anno_dfvf(
                    self.output_dir + '/' + sample + "/pathogenic_system/DFVF", main_id=dfvf_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("vfdb end")
                self.anno_signalp.add_anno_signalp(
                    self.output_dir + '/' + sample + "/pathogenic_system/SIGNALP", main_id=signalp_id,
                    specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("signalp end")
                self.anno_cyps.add_anno_cyps(
                    self.output_dir + '/' + sample + "/metabolic_system/P450", main_id=cyps_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("cyps end")
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
        self.summary_map.add_map(task_id=task_id)

    def send_files(self):
        """
        结果放置到/upload_results
        """
        dir_o = self.output_dir
        #dir_o = './output'
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        repaths = []
        regexps = []
        for sample in self.samples.keys():
            files =os.listdir(os.path.join(dir_o,sample + "/assembly_predict/assembly/assembly/"))
            for file in files:
                if re.search(r'.fna.index.',file):
                    os.remove(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly/" +file))
            if os.path.exists(os.path.join(dir_o, sample + "/project_overview")):

                self.move_dir(os.path.join(dir_o, sample + "/project_overview"), os.path.join(dir_up, sample + "/project_overview"))
            if os.path.exists(os.path.join(dir_o, sample + "/data_QC")):
                self.move_dir(os.path.join(dir_o, sample + "/data_QC"), os.path.join(dir_up, sample + "/data_QC"))
            if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/" + sample + '_gc_depth.xls')):
                self.move_file(os.path.join(dir_o, sample + "/genomic_assessment/" + sample + '_gc_depth.xls'),
                              os.path.join(dir_up, sample + "/genomic_assessment/GC_depth/" + sample + '_gc_depth.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/genome_size/" + sample + '.summary.xls')):
                self.move_file(os.path.join(dir_o, sample + "/genomic_assessment/genome_size/" + sample + '.summary.xls'),
                              os.path.join(dir_up, sample + "/genomic_assessment/genome_size/" + sample + '_genome_size_assessment.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency/" + sample + '_Kmer_frequency.xls')):
                self.move_file(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency/" + sample + '_Kmer_frequency.xls'),
                              os.path.join(dir_up, sample + "/genomic_assessment/kmer_frequency/"  + sample + '_Kmer_frequency.xls'))
                for i in ['depth_gc_10000', 'depth_gc_20000', 'depth_gc_30000', 'depth_gc_50000', 'depth_gc_100000']:
                     self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/" + i),
                       os.path.join(dir_up,  sample + "/genomic_assessment/GC_depth_png/" + i))   #guanqing.zou 20181010
                    #    os.path.join(self.file_dir, "gc_depth" + '/' + sample + '/' + i))
                #self.upload_to_s3(os.path.join(dir_o, sample + "/genomic_assessment/*"),os.path.join(self.file_dir, "gc_depth" + '/' + sample + '/'))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly"),
                              os.path.join(dir_up, sample + "/assembly_predict/assembly"))
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/scaffold"),
                              os.path.join(dir_up, sample + "/assembly_predict/assembly/scaffold"))
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/contig"),
                              os.path.join(dir_up, sample + "/assembly_predict/assembly/contig"))
                              #os.path.join(self.file_dir, "assemble/" + sample + '/contig'))

                #self.upload_to_s3(os.path.join(dir_o, sample + "/assembly_predict/assembly/scaffold/*"),
                #             os.path.join(self.file_dir, "assemble/" + sample + '/scaffold/'))
                #self.upload_to_s3(os.path.join(dir_o, sample + "/assembly_predict/assembly/contig/*"),
                #             os.path.join(self.file_dir, "assemble/" + sample + '/contig/'))
            for type in ['scaffolds','contigs']:
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/len/" + sample + '.' +type + '.len.xls')):
                    os.link(os.path.join(dir_o, sample + "/assembly_predict/assembly/len/" + sample + '.' +type + '.len.xls'),
                                   os.path.join(dir_up,sample + "/assembly_predict/assembly/" + sample + '_' + type + '_length_statistics.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly_assessment")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly_assessment"),
                              os.path.join(dir_up, sample + "/assembly_predict/assembly_assessment"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/predict/CDS_predict")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/predict/CDS_predict"),
                              os.path.join(dir_up, sample + "/assembly_predict/predict/CDS_predict"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/predict/tRNA")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/predict/tRNA"),
                              os.path.join(dir_up, sample + "/assembly_predict/predict/tRNA"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/predict/rRNA")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/predict/rRNA"),
                              os.path.join(dir_up, sample + "/assembly_predict/predict/rRNA"))
            if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/repeats")):
                self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/repeats"),
                              os.path.join(dir_up, sample + "/assembly_predict/predict/repeats"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Swissprot")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Swissprot"),
                              os.path.join(dir_up, sample + "/annotation/Swissprot"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Summary")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Summary"),
                              os.path.join(dir_up, sample + "/annotation/Summary"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Pfam")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Pfam"),
                              os.path.join(dir_up, sample + "/annotation/Pfam"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/NR")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/NR"),
                              os.path.join(dir_up, sample + "/annotation/NR"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/KEGG")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/KEGG"),
                              os.path.join(dir_up, sample + "/annotation/KEGG"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/GO")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/GO"),
                              os.path.join(dir_up, sample + "/annotation/GO"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/COG")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/COG"),
                              os.path.join(dir_up, sample + "/annotation/COG"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/CAZy/" + sample + '_anno_cazy.xls')):
                self.move_file(os.path.join(dir_o, sample + "/annotation/CAZy/" + sample + '_anno_cazy.xls'),
                              os.path.join(dir_up, sample + "/metabolic_system/CAZy/" + sample + '_anno_cazy.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/gene_cazy_family_stat.xls")):
                os.link(os.path.join(dir_o, sample + "/annotation/gene_cazy_family_stat.xls"),
                              os.path.join(dir_up, sample + "/metabolic_system/CAZy/" + sample + '_cazy_family_stat.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/gene_cazy_class_stat.xls")):
                os.link(os.path.join(dir_o, sample + "/annotation/gene_cazy_class_stat.xls"),
                              os.path.join(dir_up, sample + "/metabolic_system/CAZy/" + sample + '_cazy_class_stat.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/metabolic_system/P450")):
                self.move_dir(os.path.join(dir_o, sample + "/metabolic_system/P450"),
                              os.path.join(dir_up, sample + "/metabolic_system/P450"))
            if os.path.exists(os.path.join(dir_o, sample + "/structral_genome/promoter_predict")):
                self.move_dir(os.path.join(dir_o, sample + "/structral_genome/promoter_predict"),
                              os.path.join(dir_up, sample + "/structral_genome/promoter_predict"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/DFVF")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/DFVF"),
                              os.path.join(dir_up, sample + "/pathogenic_system/DFVF"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/SignalP")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/SignalP"),
                              os.path.join(dir_up, sample + "/pathogenic_system/SignalP"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/TMHMM")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/TMHMM"),
                              os.path.join(dir_up, sample + "/pathogenic_system/TMHMM"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/TCDB")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/TCDB"),
                              os.path.join(dir_up, sample + "/pathogenic_system/TCDB"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/PHI")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/PHI"),
                              os.path.join(dir_up, sample + "/pathogenic_system/PHI"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/CARD")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/CARD"),
                              os.path.join(dir_up, sample + "/pathogenic_system/CARD"))
            repaths += [
                [".", "", "基础分析结果文件夹", 0, "140000"],
                ["%s" % sample, "", "样品%s流程分析结果目录" % sample,0,'140001'],
                ["%s/project_overview" % sample, "", "样品%s的基因组总览目录" % sample,0,'140002'],
                ["%s/data_QC" % sample, "", "样品%s的测序数据质控统计目录" % sample,0,'140003'],
                ["%s/genomic_assessment" % sample, "", "样品%s的基因组评估结果目录" % sample,0,'140004'],
                ["%s/genomic_assessment/kmer_frequency" % sample, "", "样品%s的基因组kmer评估目录" % sample,0,'140005'],
                ["%s/genomic_assessment/GC_depth" % sample, "", "样品%s的基因组GC_depth评估目录" % sample,0,'140006'],
                ["%s/genomic_assessment/GC_depth_png" % sample, "", "" ,1],  #guanqing.zou 20181010
                ["%s/genomic_assessment/genome_size" % sample, "", "样品%s的基因组大小评估目录" % sample,0,'140007'],
                ["%s/assembly_predict" % sample, "", "样品%s的基因组组装与预测结果目录" % sample,0,'140008'],
                ["%s/assembly_predict/assembly" % sample, "", "样品%s的基因组组装结果目录" % sample,0,'140009'],
                ["%s/assembly_predict/assembly/scaffold" % sample, "", "" ,1],   ##guanqing.zou 20181010
                ["%s/assembly_predict/assembly/contig" % sample, "", "" ,1],   ##guanqing.zou 20181010
                ["%s/assembly_predict/assembly_assessment" % sample, "", "样品%s的基因组组装结果评估目录" % sample,0,'140010'],
                ["%s/assembly_predict/predict" % sample, "", "样品%s的基因组基因预测目录" % sample,0,'140011'],
                ["%s/assembly_predict/predict/CDS_predict" % sample, "", "样品%s的编码基因预测目录" % sample,0,'140012'],
                ["%s/assembly_predict/predict/repeats" % sample, "", "样品%s的基因组重复序列预测结果目录" % sample,0,'140013'],
                ["%s/assembly_predict/predict/rRNA" % sample, "", "样品%s的rRNA预测结果目录" % sample,0,'140014'],
                ["%s/assembly_predict/predict/tRNA" % sample, "", "样品%s的tRNA预测结果目录" % sample,0,'140015'],
                ["%s/assembly_predict/predict/pseudogene" % sample, "", "样品%s的假基因预测结果目录" % sample,0,'140016'],
                ["%s/annotation" % sample, "", "样品%s的基因注释结果目录" % sample,0,'140017'],
                ["%s/annotation/NR" % sample, "", "样品%s的基因NR注释结果目录" % sample,0,'140018'],
                ["%s/annotation/Swissprot" % sample, "", "样品%s的基因Swiss-Prot注释结果目录" % sample,0,'140019'],
                ["%s/annotation/Pfam" % sample, "", "样品%s的基因Pfam注释结果目录" % sample,0,'140020'],
                ["%s/annotation/COG" % sample, "", "样品%s的基因COG注释结果目录" % sample,0,'140021'],
                ["%s/annotation/GO" % sample, "", "样品%s的基因GO注释结果目录" % sample,0,'140022'],
                ["%s/annotation/KEGG" % sample, "", "样品%s的基因KEGG注释结果目录" % sample,0,'140023'],
                ["%s/annotation/KEGG/.+_kegg_pathway_img" % sample, "", "KEGG注释pathway通路图片目录",0,'140024'],
                ["%s/annotation/Summary" % sample, "", "样品%s的基因组注释结果汇总目录" % sample,0,'140025'],
                ["%s/structral_genome" % sample, "", "样品%s的结构基因组分析目录" % sample,0,'140026'],
                ["%s/structral_genome/promoter_predict" % sample, "", "样品%s的启动子预测结果目录" % sample,0,'140027'],
                ["%s/metabolic_system" % sample, "", "样品%s的代谢系统分析目录" % sample,0,'140028'],
                ["%s/metabolic_system/CAZy" % sample, "", "样品%s的碳水化合物活性酶注释结果目录" % sample,0,'140029'],
                ["%s/metabolic_system/P450" % sample, "", "样品%s的细胞色素P450分析结果目录" % sample,0,'140030'],
                ["%s/pathogenic_system" % sample, "", "样品%s的致病系统分析目录" % sample,0,'140031'],
                ["%s/pathogenic_system/DFVF" % sample, "", "样品%s的毒力基因预测结果目录" % sample,0,'140032'],
                ["%s/pathogenic_system/SignalP" % sample, "", "样品%s的分泌蛋白预测结果目录" % sample,0,'140033'],
                ["%s/pathogenic_system/CARD" % sample, "", "样品%s的耐药基因预测结果目录" % sample,0,'140034'],
                ["%s/pathogenic_system/PHI" % sample, "", "样品%s的病原菌与宿主互作分析结果目录" % sample,0,'140035'],
                ["%s/pathogenic_system/TCDB" % sample, "", "样品%s的转运蛋白分析结果目录" % sample,0,'140036'],
                ["%s/pathogenic_system/TMHMM" % sample, "", "样品%s的跨膜蛋白分析结果目录" % sample,0,'140037'],
                ["%s/pathogenic_system/secretion_system" % sample, "", "样品%s的分泌系统分析结果目录" % sample,0,'140038'],
            ]
            regexps += [
                [r"%s/project_overview/genome_overview.xls" % sample, "", "基因组概况统计",0,'140039'],
                [r"%s/data_QC/.+_Illumina_statistics.xls" % sample, "xls", "二代测序数据质控统计表",0,'140040'],
                [r"%s/data_QC/.+_PacBio_statistics.xls" % sample, "xls", "三代测序数据质控统计表",0,'140041'],
                [r"%s/genomic_assessment/kmer_frequency/.+_Kmer_frequency.xls" % sample, "xls", "基因组评估kmer频率表",0,'140042'],
                [r"%s/genomic_assessment/GC_depth/.+_gc_depth.xls" % sample, "xls", "基因组评估GC_depth表",0,'140043'],
                [r"%s/genomic_assessment/genome_size/.+_genome_size_assessment.xls" % sample, "xls", "基因组大小评估结果表",0,'140044'],
                [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'140045'],
                [r"%s/assembly_predict/assembly/.+_assembly_details.xls" % sample, "xls", "组装结果详情表",0,'140046'],
                [r"%s/assembly_predict/assembly/.+_assembly_scaffold_details.xls" % sample, "xls","组装结果scaffolds详情表",0,'140047'],
                [r"%s/assembly_predict/assembly/.+_assembly_contig_details.xls" % sample, "xls", "组装结果contigs详情表",0,'140048'],
                [r"%s/assembly_predict/assembly/.+_scaffolds_length_statistics.xls" % sample, "xls", "组装结果scaffolds长度统计详情表",0,'140049'],
                [r"%s/assembly_predict/assembly/.+_contigs_length_statistics.xls" % sample, "xls", "组装结果contigs长度统计详情表",0,'140050'],
                [r"%s/assembly_predict/assembly/.+.agp" % sample, "", "基因组组装的scaffold与contig对应关系文件",0,'140051'],
                [r"%s/assembly_predict/assembly/.+_scaf.fna" % sample, "", "基因组组装的scaffold文件",0,'140052'],
                [r"%s/assembly_predict/assembly_assessment/.+_cegma.xls" % sample, "", "基因组组装结果CEGMA评估表",0,'140053'],
                [r"%s/assembly_predict/assembly_assessment/.+_busco.xls" % sample, "", "基因组组装结果BUSCO评估表",0,'140054'],
                [r"%s/assembly_predict/assembly/.+_ctg.fna" % sample, "", "基因组组装的contigs文件",0,'140055'],
                [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'140045'],
                [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.faa" % sample, "", "编码基因预测氨基酸序列",0,'140056'],
                [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.fnn" % sample, "", "编码基因预测核苷酸序列",0,'140057'],
                [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.gff" % sample, "", "编码基因预测gff格式统计表",0,'140058'],
                [r"%s/assembly_predict/predict/CDS_predict/.+_CDS_statistics.xls" % sample, "xls", "编码基因预测统计表",0,'140059'],
                [r"%s/assembly_predict/predict/CDS_predict/.+_CDS_length.xls" % sample, "xls", "编码基因长度统计表",0,'140060'],
                [r"%s/assembly_predict/predict/repeats/.+_Rep.tbl" % sample, "", "散在重复序列预测tbl结果文件",0,'140061'],
                [r"%s/assembly_predict/predict/repeats/.+_Rep.out" % sample, "", "散在重复序列预测out结果文件",0,'140062'],
                [r"%s/assembly_predict/predict/repeats/.+_Rep.gff" % sample, "", "散在重复序列预测gff结果文件",0,'140063'],
                [r"%s/assembly_predict/predict/rRNA/.+_rRNA.fnn" % sample, "", "rRNA预测序列文件",0,'140064'],
                [r"%s/assembly_predict/predict/rRNA/.+_rRNA.gff" % sample, "", "rRNA预测gff统计文件",0,'140065'],
                [r"%s/assembly_predict/predict/tRNA/.+_tRNA.fnn" % sample, "", "tRNA预测序列文件",0,'140066'],
                [r"%s/assembly_predict/predict/tRNA/.+_tRNA.gff" % sample, "", "tRNA预测结果文件",0,'140067'],
                [r"%s/assembly_predict/predict/tRNA/.+_tRNA.struc" % sample, "", "tRNA预测二级结构文件",0,'140068'],
                [r"%s/assembly_predict/predict/pseudogene/.+_pseudogene.xls" % sample, "", "假基因预测结果文件",0,'140069'],
                [r"%s/annotation/NR/.+_anno_nr.xls" % sample, "xls", "基因NR注释结果详情表",0,'140070'],
                [r"%s/annotation/Swissprot/.+_anno_swissprot.xls" % sample, "xls", "基因Swiss-Prot注释结果详情表",0,'140071'],
                [r"%s/annotation/Pfam/.+_anno_pfam.xls" % sample, "xls", "基因Pfam注释结果详情表",0,'140072'],
                [r"%s/annotation/COG/.+_cog_anno.xls" % sample, "xls", "基因COG注释结果详情表",0,'140073'],
                [r"%s/annotation/COG/.+_cog_summary.xls" % sample, "xls", "COG注释结果汇总统计表",0,'140074'],
                [r"%s/annotation/GO/.+_go_anno.xls" % sample, "xls", "基因GO注释结果详情表",0,'140075'],
                [r"%s/annotation/GO/.+_go_level2_statistics.xls" % sample, "xls", "GO注释level2统计文件",0,'140076'],
                [r"%s/annotation/KEGG/.+_kegg_anno.xls" % sample, "xls", "KEGG注释结果文件",0,'140077'],
                [r"%s/annotation/KEGG/.+_kegg_level_stat.xls" % sample, "xls", "KEGG注释level统计文件",0,'140078'],
                [r"%s/annotation/Summary/.+_anno_summary.xls" % sample, "xls", "基因组注释结果汇总表",0,'140079'],
                [r"%s/structral_genome/promoter_predict/.+_promoter_result.xls" % sample, "xls", "启动子预测结果表",0,'140080'],
                [r"%s/metabolic_system/CAZy/.+_anno_cazy.xls" % sample, "xls", "碳水化合物活性酶注释结果表",0,'140081'],
                [r"%s/metabolic_system/CAZy/.+_cazy_family_stat.xls" % sample, "xls", "碳水化合物活性酶注释family统计表",0,'140082'],
                [r"%s/metabolic_system/CAZy/.+_cazy_class_stat.xls" % sample, "xls", "碳水化合物活性酶注释class统计表",0,'140083'],
                [r"%s/metabolic_system/P450/.+_anno_p450.xls" % sample, "xls", "细胞色素P450分析结果文件",0,'140084'],
                [r"%s/pathogenic_system/DFVF/.+_dfvf_align_table.xls" % sample, "xls", "Diamond比对结果表",0,'140085'],
                [r"%s/pathogenic_system/DFVF/.+_dfvf.anno.xls" % sample, "xls", "毒力基因注释结果表",0,'140086'],
                [r"%s/pathogenic_system/SignalP/.+_SignalP.xls" % sample, "xls", "分泌蛋白预测结果文件",0,'140087'],
                [r"%s/pathogenic_system/CARD/.+_card_align.xls" % sample, "xls", "耐药基因预测表",0,'140088'],
                [r"%s/pathogenic_system/CARD/.+_card_anno.xls" % sample, "xls", "耐药基因注释表",0,'140089'],
                [r"%s/pathogenic_system/CARD/.+_card_category.xls" % sample, "xls", "耐药基因分类统计表",0,'140090'],
                [r"%s/pathogenic_system/PHI/.+_phi_align.xls" % sample, "xls", "病原菌与宿主互作相关基因预测表",0,'140091'],
                [r"%s/pathogenic_system/PHI/.+_phi_anno.xls" % sample, "xls", "病原菌与宿主互作相关基因注释表",0,'140092'],
                [r"%s/pathogenic_system/PHI/.+_phi_stat.xls" % sample, "xls", "病原菌与宿主互作分析PHI ID统计表",0,'140093'],
                [r"%s/pathogenic_system/PHI/.+_phi_phenotype.xls" % sample, "xls", "病原菌与宿主互作分析表型统计表",0,'140094'],
                [r"%s/pathogenic_system/TCDB/.+_tcdb_align.xls" % sample, "xls", "转运蛋白预测表",0,'140095'],
                [r"%s/pathogenic_system/TCDB/.+_tcdb_align_top1.xls" % sample, "xls", "转运蛋白预测表(TOP1)",0,'140096'],
                [r"%s/pathogenic_system/TCDB/.+_tcdb_anno.xls" % sample, "xls", "转运蛋白注释表",0,'140097'],
                [r"%s/pathogenic_system/TCDB/.+_whole_genome_tcdb_anno.xls" % sample, "xls", "转运蛋白注释表",0,'140098'],
                [r"%s/pathogenic_system/TMHMM/.+_tmhmm_anno.xls" % sample, "xls", "跨膜蛋白分析结果详情表",0,'140099'],
                [r"%s/pathogenic_system/TMHMM/.+_genome_tmhmm_anno.xls" % sample, "xls", "跨膜蛋白注释表",0,'140100'],
                [r"%s/pathogenic_system/secretion_system/.+_secretion_system_genes.xls" % sample, "xls","分泌系统相关基因统计表",0,'140101'],
                [r"%s/pathogenic_system/secretion_system/.+_secretion_system_type.xls" % sample, "xls","分泌系统的类型统计表",0,'140102'],
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
                        lib_type[lin[0]][des] = lin[5]
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
                            self.set_error('raw_dir的list.txt文件格式有误', code="12100101")
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
                            self.set_error('asse_dir的list.txt文件格式有误', code="12100102")
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
                            self.set_error('raw_dir的list.txt文件格式有误', code="12100103")
                with open(assemble_list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 2:
                            assemble_sample[line[0]] = line[0]
                        else:
                            self.set_error('asse_dir的list.txt文件格式有误', code="12100104")
                for key in raw_sample.keys():
                    if key in assemble_sample.keys():
                        self.samples[key]=key
                    else:
                        self.set_error('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致', code="12100105")

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
                            self.set_error('asse_dir的list.txt文件格式有误', code="12100106")
                        if line[2] not in ['chromosome', 'plasmid']:
                            self.set_error('asse_dir的list.txt文件的Genome Type有误！', code="12100107")

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
                            self.set_error('raw_dir的list.txt文件格式有误', code="12100108")
                with open(assemble_list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 3:
                            assemble_sample[line[0]] = line[0]
                        else:
                            self.set_error('asse_dir的list.txt文件格式有误', code="12100109")
                        if line[2] not in ['chromosome', 'plasmid']:
                            self.set_error('asse_dir的list.txt文件的Genome Type有误！', code="12100110")
                for key in raw_sample.keys():
                    if key in assemble_sample.keys():
                        self.samples[key]=key
                    else:
                        self.set_error('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致', code="12100111")

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
                    for i in range(1,chr_num +1):
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
            self.set_error('需要移动到output目录的文件夹不存在。', code="12100112")
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
