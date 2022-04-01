# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""宏基因组binning分析工作流"""

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
from mbio.packages.metagbin.common_function import link_dir


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

class MetagbinApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        宏基因组binning workflow option参数设置
        """
        self._sheet = wsheet_object
        super(MetagbinApiWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {'name': 'speciman_info', 'type': 'infile', 'format': 'meta_genomic.specimen_info'},  # 样本集信息表
            {'name': 'qc', 'type': 'bool', 'default': False},  # 是否需要质控
            {'name': 'qc_quality', 'type': 'int', 'default': 20},  # 质控质量值标准
            {'name': 'qc_length', 'type': 'int', 'default': 30},  # 质控最短序列长度
            {'name': 'rm_host', 'type': 'bool', 'default': False},  # 是否需要去除宿主
            {'name': 'ref_database', 'type': 'string', 'default': ''},  # 宿主参考序列库中对应的物种名，eg：E.coli ,B.taurus
            {'name': 'ref_undefined', "type": 'infile', 'format': 'sequence.fasta_dir'},
            {'name': 'ref_undefined_name', 'type': 'string', 'default': 'undefined'},  # 自定义参考宿主名称，适应页面参数
            {'name': 'cdhit_identity', 'type': 'float', 'default': 0.95},  # 基因序列聚类相似度
            {'name': 'cdhit_coverage', 'type': 'float', 'default': 0.9},  # 基因序列聚类覆盖度
            {'name': 'min_contig', 'type': 'int', 'default': 1000},  # 基因序列聚类覆盖度
            {'name': 'sofware_bin', 'type': 'string', 'default': 'metabat'},  # bin运行软件方法
            {"name": "assemble_dir", "type": "infile", "format": "metagbin.assemble_dir"},  # 组装文件fasta文件夹
            {'name': 'clean_fq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.bin_name = ''
        self.remote_dir = self._sheet.output + '/'

    def check_options(self):
        """
        检查参数
        """
        pass

    def set_output(self, event):
        """
        将各个模块的结果输出至output
        """
        pass

    def get_binname(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            bin = lines[1].rstrip('\n\r\t').split('\t')[0]
        return bin

    def get_taxon(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            bin_name = lines[1].rstrip('\n\r\t').split('\t')[-1]
        return bin_name

    def get_sof(self, file):
        with open(file, 'rb') as f2:
            lines = f2.readlines()
            line = lines[1].rstrip('\r\n\t').split('\t')
            result = line[0]
        return result

    def get_in(self, file):
        with open (file, 'r') as f:
            lines = f.readlines()
            line = lines[1].rstrip('\r\n\t').split('\t')
            insert = line[0]
            max = line[1]
        return insert,max

    def end(self):
        self.run_move()
        self.send_files()
        super(MetagbinApiWorkflow, self).end()

    def run_api(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.output_dir = "/mnt/lustre/users/sanger/sg-users/gh/bin/workflow/output"
        self.bin_name = self.get_binname(self.output_dir + '/Binning_result/all.bin.summary.xls')
        self.api_bin_stat()
        self.api_assemble()
        self.api_predict()
        self.api_anno()
        self.api_16s_taxon()

    def run_move(self):
        bin_name = 'G_' + self.bin_name
        if os.path.exists(self.output_dir + '/' + bin_name + "/Assembly_summary"):
            shutil.rmtree(self.output_dir + '/' + bin_name + "/Assembly_summary")

    def api_bin_stat(self):
        min_contig=self.option('min_contig')
        sofware_bin = self.option('sofware_bin')
        pi_path = self.api.api("metagbin.common_api")
        sofware_bin.replace("metabat","metabat2")
        bin_params = {"soft": sofware_bin, "cuto_ff": min_contig}
        bin_id = pi_path.add_main('bin', name="Bin_origin", params= bin_params, desc ="binning的结果内容!")
        mongo_key = 'bin_id,len,scf_num,lon_len,n50,aver_scf,comp,cont,strain_heter,domain'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/all.bin.summary.xls', 'bin_detail', bin_id,
                                mongo_key, has_head=True, main_name = 'bi_id',convert_dic={'comp':float})
        mongo_key2 = 'bin_id,total_marker,uniq_marker,a,b,c,d,e,f'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/all.marker.xls', 'bin_marker', bin_id,
                                mongo_key2, has_head=True, main_name = 'bi_id')
        mongo_key3 = 'bin_id,taxon'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/summary.anno.xls', 'bin_taxon', bin_id,
                                mongo_key3, has_head=True, main_name='bi_id')
        mongo_key4 = 'sca_id,bin_id,start,end,gc,coverage'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/all_bins.scf.xls', 'bin_scaf', bin_id,
                                mongo_key4, has_head=True, main_name='bi_id')
        mongo_key4 = 'sca_id,bin_id,start,end,strand,seq'
        pi_path.add_main_detail(self.output_dir + '/Binning_result/all_bins.16s.xls', 'bin_s16', bin_id,
                                mongo_key4, has_head = True, main_name='bi_id',main_table='bin', update_dic=
                                {'main_id': bin_id, 'seq_path': self.remote_dir + 'Binning_result/merge.bam.gz',
                                 'ref_seq': self.remote_dir + 'Binning_result/all.uniContigs.fa.gz', 'path': self.remote_dir + 'Binning_result/bin/'})

    def api_assemble(self):
        bin_name = 'G_' + self.bin_name
        api_assembly = self.api.api('metagbin.assembly')
        soft = self.get_sof(self.output_dir + '/' + bin_name + "/Assembly_summary/Assembly_soft.txt")
        insert,max = self.get_in(self.output_dir + '/' + bin_name + "/Assembly_summary/config.txt")
        assembly_params = {'genome_id': bin_name, 'soft': soft, 'task_type': 2, 'submit_location': "assemble"}
        assembly_id = api_assembly.add_assembly(bin_name, params=assembly_params)
        file_path = self.output_dir + '/' + bin_name + "/Assembly/Assembly_summary/" + bin_name + \
                    "_assembly_summary.xls"
        bin_cover_path1 = self.output_dir + '/' + bin_name + "/Assembly_summary/Bin_coverage.xls"
        assem_cover_path2 = self.output_dir + '/' + bin_name + "/Assembly/Assembly_summary/Genome_coverage.xls"
        genome_path = self.remote_dir + bin_name + "/Assembly/Assembly_summary/" + bin_name + "_scaffold.fa"
        api_assembly.add_assembly_detail(assembly_id, bin_name, insert, max, file_path=file_path, bin_cover_path1=bin_cover_path1,
                                         assem_cover_path2=assem_cover_path2, genome_path=genome_path)
        self.logger.info('组装结果统计写入mongo成功')

        api_assembly_assess = self.api.api('metagbin.assembly_assess')
        assess_path = self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/' + bin_name + '_assess.xls'
        soft = self.get_sof(self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/Assess_soft.txt')
        assess_params = { 'genome_id': bin_name, 'soft': soft, 'task_type': 2, 'submit_location': "assemble_assess"}
        assess_id = api_assembly_assess.add_assembly_assess(bin_name, params=assess_params)
        api_assembly_assess.add_assembly_assess_detail(assess_id, soft, bin_name, file_path=assess_path)
        self.logger.info('组装评估结果写入mongo成功')

        api_gc_depth = self.api.api("metagbin.assess_gc")
        gc_params = {"genome_id": bin_name, 'soft': 'Bowtie2', 'task_type': 2, 'submit_location': "assemble_gc"}
        gc_path = self.remote_dir + bin_name + "/Assembly/GC_depth/"
        api_gc_depth.add_assembly_gc(bin_name, gc_path, params=gc_params)
        self.logger.info('组装Gc_depth结果写入mongo成功')

    def api_predict(self):
        bin_name = 'G_' + self.bin_name
        api_predict_cds = self.api.api('metagbin.predict_cds')
        cds_params = {"genome_id": bin_name, 'soft': 'Glimmer', 'task_type': 2, 'submit_location': "predict_cds"}
        pre_id = api_predict_cds.add_predict_cds(bin_name, params=cds_params)
        sample_stat = self.output_dir + '/' + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS_statistics.xls"
        api_predict_cds.add_predict_cds_stat(pre_id, bin_name, sample_stat)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.gff"
        fnn_path = self.remote_dir + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.fnn"
        faa_path = self.remote_dir + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.faa"
        gff_path = self.remote_dir + bin_name + "/Gene_predict/CDS_predict/" + bin_name + "_CDS.gff"
        api_predict_cds.add_predict_cds_detail(pre_id, bin_name, predict_gff = predict_gff, fnn_path =fnn_path, faa_path = faa_path, gff_path = gff_path)
        distribute = self.output_dir + '/' + bin_name + "/Gene_predict/length_distribute.txt"
        api_predict_cds.add_predict_cds_bar(pre_id, bin_name, distribute)
        self.logger.info('基因预测cds结果写入mongo成功')

        api_predict_repeat = self.api.api('metagbin.predict_repeat')
        repeat_params = {"genome_id": bin_name, 'soft': 'Tandem Repeats Finder', "parameters":"Match=2, Mismatch=7, Delta=7, PM=80, PI=10, Minscore=80, MaxPeriod=500", 'task_type': 2, 'submit_location': "predict_repeat"}
        repeat_id = api_predict_repeat.add_predict_repeat(bin_name, params=repeat_params)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/repeats/" + bin_name + "_TRF.gff"
        if os.path.exists(predict_gff):
            api_predict_repeat.add_predict_repeat_detail(repeat_id, bin_name, predict_gff)
        self.logger.info('重复序列预测结果写入mongo成功')

        api_predict_rrna = self.api.api('metagbin.predict_rrna')
        rrna_params = {"genome_id": bin_name, 'soft': 'Barrnap', 'task_type': 2, 'submit_location': "predict_rrna"}
        rrna_id = api_predict_rrna.add_predict_rrna(bin_name, params=rrna_params)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/rRNA/" + bin_name + "_rRNA.gff"
        rna_path = self.remote_dir + bin_name + "/Gene_predict/rRNA/" + bin_name + '-16S_rRNA.fa'
        fnn = self.output_dir + '/' + bin_name + "/Gene_predict/rRNA/" + bin_name + "_rRNA.fnn"
        api_predict_rrna.add_predict_rrna_detail(rrna_id, bin_name, predict_gff)
        if os.path.exists(predict_gff):
            api_predict_rrna.add_predict_rrna_detail(rrna_id, bin_name, predict_gff = predict_gff)
            if os.path.exists(self.output_dir + '/' + bin_name + "/Gene_predict/rRNA/" + bin_name + '-16S_rRNA.fa'):
                api_predict_rrna.add_predict_rrna_seq(rrna_id, bin_name, fnn=fnn, predict_gff=rna_path)
            else:
                api_predict_rrna.add_predict_rrna_seq(rrna_id, bin_name, fnn=fnn)
        self.logger.info('rRNA预测结果写入mongo成功')

        api_predict_trna = self.api.api('metagbin.predict_trna')
        trna_params = {"genome_id": bin_name, 'soft': 'tRNAscan-SE', 'task_type': 2, 'submit_location': "predict_trna"}
        trna_id = api_predict_trna.add_predict_trna(bin_name, params=trna_params)
        predict_gff = self.output_dir + '/' + bin_name + "/Gene_predict/tRNA/" + bin_name + "_tRNA.gff"
        if os.path.exists(predict_gff):
            api_predict_trna.add_predict_trna_detail(trna_id, bin_name, predict_gff)
        self.logger.info('tRNA预测结果写入mongo成功')

        api_genome_path = self.api.api("metagbin.common_api")
        soft = self.get_sof(self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/Assess_soft.txt')
        assess_path = self.output_dir + '/' + bin_name + '/Assembly/Genome_assessment/' + bin_name + '_assess.xls'
        complete = ''
        if soft == 'BUSCO':
            with open(assess_path, 'rb') as f:
                line = f.readlines()
                line = line[1].strip().split('\t')
                complete = line[0]
            assess_taxon = 'Bacteria'
        else:
            with open(assess_path, 'rb') as f1:
                line = f1.readlines()
                line = line[1].strip().split('\t')
                complete = line[0]
            assess_taxon = 'Archaea'
        genome_params = {"genome_id": bin_name, "complete": complete}
        api_genome_path.add_main("genome", name="Genome_origin", params=genome_params, others={"bin_id": self.bin_name, 'genome_id': bin_name, "taxon": assess_taxon})
        self.logger.info('基因组管理表导入mongo成功')

    def api_anno(self):
        bin_name = 'G_' + self.bin_name
        anno_nr = self.api.api("metagbin.anno_nr")
        nr_params = {"genome_id": bin_name, "nr": "Diamond", 'task_type': 2, 'submit_location': "anno_nr"}
        nr_id = anno_nr.add_anno_nr(params=nr_params)
        anno_nr.add_anno_nr_detail(nr_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/NR/" + bin_name + "_anno_nr.xls")
        anno_cog = self.api.api("metagbin.anno_cog")
        cog_params = {"genome_id": bin_name, "cog": "Diamond", 'task_type': 2, 'submit_location': "anno_cog"}
        cog_id = anno_cog.add_anno_cog(params=cog_params)
        anno_cog.add_anno_cog_detail(cog_id, bin_name,
                                          self.output_dir + '/' + bin_name + "/annotation/COG/" + bin_name + "_cog_anno.xls")
        anno_kegg = self.api.api("metagbin.anno_kegg")
        kegg_params = {"genome_id": bin_name, "kegg": "Diamond", 'task_type': 2, 'submit_location': "anno_kegg"}
        kegg_id = anno_kegg.add_anno_kegg(self.remote_dir, "/annotation/KEGG/", "_kegg_pathway_img",
                                               params = kegg_params)
        anno_kegg.add_anno_kegg_detail(kegg_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/KEGG/" + bin_name + "_kegg_anno.xls")
        anno_kegg.add_anno_kegg_level(kegg_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/KEGG/" + bin_name + "_kegg_level_stat.xls")
        anno_cazy = self.api.api("metagbin.anno_cazy")
        cazy_params = {"genome_id": bin_name, "cazy": "hmmscan", 'task_type': 2, 'submit_location': "annocazy"}
        cazy_id = anno_cazy.add_anno_cazy(params=cazy_params)
        anno_cazy.add_anno_cazy_detail(cazy_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/CAZY/" + bin_name + "_anno_cazy.xls")
        api_card = self.api.api("metagbin.anno_card")
        card_params = {"genome_id": bin_name, "card": "Diamond", 'task_type': 2, 'submit_location': "annocard"}
        api_card.add_anno_card(self.output_dir + '/' + bin_name + "/annotation/CARD", main=True, genome_id=bin_name, params=card_params)
        api_sum = self.api.api("metagbin.anno_summary")
        api_stat = self.api.api("metagbin.anno_stat")
        sum_params = {"genome_id": bin_name, 'task_type': 2, 'submit_location': "annosummary"}
        stat_params = {"genome_id": bin_name, 'task_type': 2, 'submit_location': "anno_stat"}
        sum_id = api_sum.add_anno_summary(params=sum_params)
        stat_id = api_stat.add_anno_summary(params=stat_params)
        api_sum.add_anno_summary_detail(sum_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/Summary/" + bin_name + "_anno_summary.xls")
        api_stat.add_anno_stat_detail(stat_id, bin_name, self.output_dir + '/' + bin_name + "/annotation/Summary/" + bin_name + "_anno_stat.xls")

    def api_16s_taxon(self):
        bin_name = 'G_' + self.bin_name
        if os.path.exists(self.output_dir + '/' + bin_name + "/Taxon_identify/16S_identify/blasr.xls"):
            pi_path = self.api.api("metagbin.common_api")
            taxonomy_params = {"query": bin_name, "db": "silva128", 'task_type': 2, 'submit_location': "identif_sixteen"}
            taxonomy_params = json.dumps(taxonomy_params, sort_keys=True, separators=(',', ':'))
            s16_id = pi_path.add_main("identif_16s", name="16S_origin", params=taxonomy_params)
            pi_path.add_main_detail(self.output_dir + '/' + bin_name + "/Taxon_identify/16S_identify/blasr.xls", "identif_16s_detail", s16_id, "genome_id,taxon,identify", has_head=True,
                                main_name="s16_id", main_table='identif_16s', update_dic={'main_id': s16_id})
            pi_path.updat_genome_taxon(self._sheet.id, bin_name, self.output_dir + '/' + bin_name + "/Taxon_identify/16S_identify/blasr.xls")

    def send_files(self):
        bin_name = 'G_' + self.bin_name
        dir_o = self.output_dir
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        repaths = []
        regexps = []
        if os.path.exists(os.path.join(dir_o, "Binning_result")):
            link_dir(os.path.join(dir_o, "Binning_result"), os.path.join(dir_up,"Binning_result"))
        if os.path.exists(os.path.join(dir_o, bin_name)):
            link_dir(os.path.join(dir_o, bin_name), os.path.join(dir_up, bin_name))
        repaths += [
            [".", "", "工作流分析结果目录", 0],
            ["Binning_result", "", "Binning的相关结果目录", 0],
            ["Binning_result/all.scaffolds.xls", "xls", "所有scaffold的统计信息", 0],
            ["Binning_result/all.marker.xls", "xls", "所有bin的marker基因统计信息", 0],
            ["Binning_result/all.bin.summary.xls", "xls", "每个bin的统计信息", 0],
            ["Binning_result/all.bin.rRNA.xls", "xls", "存在16sRNA的bin的统计信息", 0],
            ["Binning_result/all.bin.anno.xls", "xls", "每个bin的物种注释信息", 0],
            ["Binning_result/bin", "", "bin的序列文件目录", 0],
            ["%s" %bin_name, "", "%s的结果目录" %bin_name, 0],
            ["%s/Assembly" % bin_name, "", "%s的基因组组装结果目录" % bin_name, 0],
            ["%s/Assembly/Assembly_summary" % bin_name, "", "%s的基因组组装统计结果目录" % bin_name, 0],
            ["%s/Assembly/Assembly_summary/Bin_coverage.xls" % bin_name, "", "%s的基因组binning覆盖度结果表" % bin_name, 0],
            ["%s/Assembly/Assembly_summary/Genome_coverage.xls" % bin_name, "", "%s的基因组基因组Genome覆盖度结果表" % bin_name, 0],
            ["%s/Assembly/Genomic_assessment" % bin_name, "", "%s的基因组评估结果目录" % bin_name, 0],
            ["%s/Assembly/GC_depth" % bin_name, "", "%s的基因组GC_depth" % bin_name, 0],
            ["%s/Gene_predict" % bin_name, "", "%s的基因组的基因预测目录" % bin_name, 0],
            ["%s/Gene_predict/CDS_predict" % bin_name, "", "%s基因组编码基因预测目录" % bin_name, 0],
            ["%s/Gene_predict/rRNA" % bin_name, "", "%s基因组rRNA预测结果目录" % bin_name, 0],
            ["%s/Gene_predict/tRNA" % bin_name, "", "%s基因组tRNA预测结果目录" % bin_name, 0],
            ["%s/Gene_predict/Repeats" % bin_name, "", "%s基因组重复序列预测结果目录" % bin_name, 0],
            ["%s/annotation" % bin_name, "", "%s基因组的功能注释结果目录" % bin_name, 0],
            ["%s/annotation/NR" % bin_name, "", "%s基因组的NR注释结果目录" % bin_name, 0],
            ["%s/annotation/COG" % bin_name, "", "%s基因组的COG注释结果目录" % bin_name, 0],
            ["%s/annotation/KEGG" % bin_name, "", "%s基因组的KEGG注释结果目录" % bin_name, 0],
            ["%s/annotation/CAZY" % bin_name, "", "%s基因组的CAZY注释结果目录" % bin_name, 0],
            ["%s/annotation/CARD" % bin_name, "", "%s基因组的CARD注释结果目录" % bin_name, 0],
            ["%s/annotation/Summary" % bin_name, "", "%s基因组的注释结果汇总目录" % bin_name, 0],
            ["%s/Taxon_identify" % bin_name, "", "%s基因组物种鉴定模块目录" % bin_name, 0],
            ["%s/Taxon_identify/16S_identify" % bin_name, "", "%s基因组的16S物种分类结果目录" % bin_name, 0],
            ["%s/Taxon_identify/16S_identify/blasr_result.xls" % bin_name, "xls", "%s基因组的16S物种分类结果表" % bin_name, 0],
        ]
        regexps += [
            [r"Binning_result/bin/bin[0-9]+\.fa", "", "bin的序列文件", 0],
            [r"%s/Assembly/Assembly_summary/.+_assembly_scaffold_details.xls" % bin_name, "", "%s的基因组组装scaffold统计结果" % bin_name, 0],
            [r"%s/Assembly/Assembly_summary/.+_scafflod.fa" % bin_name, "", "%s的基因组组装序列结果文件" % bin_name, 0],
            [r"%s/Assembly/Genomic_assessment/.+_assess.xls" % bin_name, "xls", "%s的基因组评估结果文件" % bin_name, 0],
            [r"%s/Assembly/GC_depth/.+_1000.png" % bin_name, "", "%s的基因组GC_depth_1000统计图(PNG)" % bin_name, 0],
            [r"%s/Assembly/GC_depth/.+_1000.svg" % bin_name, "", "%s的基因组GC_depth_1000统计图(SVG)" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS.faa" % bin_name, "", "%s基因组编码基因预测蛋白序列" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS.fnn" % bin_name, "", "%s基因组编码基因预测核苷酸序列" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS.gff" % bin_name, "", "%s基因组编码基因预测统计表(gff格式)" % bin_name, 0],
            [r"%s/Gene_predict/CDS_predict/.+_CDS_statistics.xls" % bin_name, "", "%s基因组编码基因预测统计表" % bin_name, 0],
            [r"%s/Gene_predict/rRNA/.+_rRNA.fnn" % bin_name, "", "%s基因组rRNA预测核苷酸序列文件" % bin_name, 0],
            [r"%s/Gene_predict/rRNA/.+_rRNA.gff" % bin_name, "", "%s基因组rRNA预测统计表(gff格式)" % bin_name, 0],
            [r"%s/Gene_predict/rRNA/.+-16s_rRNA.fnn" % bin_name, "", "%s基因组rRNA预测16S核苷酸序列文件" % bin_name, 0],
            [r"%s/Gene_predict/tRNA/.+_tRNA.gff" % bin_name, "", "%s基因组tRNA预测统计表(gff格式)" % bin_name, 0],
            [r"%s/Gene_predict/tRNA/.+_tRNA.fnn" % bin_name, "", "%s基因组tRNA预测核苷酸序列文件" % bin_name, 0],
            [r"%s/Gene_predict/Repeats/.+_TRF.gff" % bin_name, "", "%s基因组重复序列预测统计表(gff格式)" % bin_name, 0],
            [r"%s/annotation/NR/.+_anno_nr.xls" % bin_name, "", "%s基因组的NR注释结果表" % bin_name, 0],
            [r"%s/annotation/COG/.+_cog_summary.xls" % bin_name, "", "%s基因组的COG注释结果分类统计表" % bin_name, 0],
            [r"%s/annotation/COG/.+_cog_anno.xls" % bin_name, "", "%s基因组的COG注释结果表" % bin_name, 0],
            [r"%s/annotation/KEGG/" % bin_name, "", "%s基因组的KEGG注释结果目录" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_anno_cazy.xls" % bin_name, "", "%s基因组的CAZY注释分类结果表" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_cazy_class_stat.xls" % bin_name, "", "%s基因组的CAZY注释class分类结果表" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_cazy_family_stat.xls" % bin_name, "", "%s基因组的CAZY注释family分类结果表" % bin_name, 0],
            [r"%s/annotation/CAZY/.+_cazy_parse_anno.xls" % bin_name, "", "%s基因组的CAZY注释结果表" % bin_name, 0],
            [r"%s/annotation/CARD/.+_card_align.xls" % bin_name, "", "%s基因组的CARD比对结果表" % bin_name, 0],
            [r"%s/annotation/CARD/.+_card_anno.xls" % bin_name, "", "%s基因组的CARD注释结果表" % bin_name, 0],
            [r"%s/annotation/CARD/.+_card_category.xls" % bin_name, "", "%s基因组的CARD注释分类结果表" % bin_name, 0],
            [r"%s/annotation/Summary/.+_anno_summary.xls" % bin_name, "", "%s基因组的注释结果汇总表" % bin_name, 0],
        ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

    def run(self):
        """
        运行:genome_workflow
        :return:
        """
        task_info = self.api.api('task_info.metagbin_task_info')
        task_info.add_task_info()
        self.run_api()
        self.logger.info(self.work_dir)
        gevent.spawn_later(5, self.end)
        super(MetagbinApiWorkflow, self).run()




