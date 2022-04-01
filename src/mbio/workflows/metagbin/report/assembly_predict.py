# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'@20190130#last_modifild:20200420@gaohao

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
import datetime
import gevent
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir


class AssemblyPredictWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        宏基因组binning组装和预测交互分析
        :param wsheet_object:
        :return:
        """
        self._sheet = wsheet_object
        super(AssemblyPredictWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "insert_size", "type": "int"},  #样本的插入片段大小
            {"name": "max_length", "type": "int"}, #样本的最大读长
            {"name": "bin_id", "type": "string"},    #挑选的bin的名称
            {"name": "bam_file", "type": "infile", "format": "metagbin.file_gz, metagbin.file_spr_dir"},  # 交互分析输入cleandata的fastq压缩文件，质控后的文件
            {"name": "ref_fa", "type": "infile","format": "sequence.fasta_dir"},   #binning模块的bin结果目录--文件夹，根据此参数和bin_id去找到参考序列
            {"name": "taxon", "type": "string", "default": "Bacteria"}, #古菌还是细菌
            {"name": "workflow", "type": "string", "default": "interaction"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "int", "default": 2},
            {"name": "main_id", "type": "string", "default": ""},
            {"name": "update_info", "type": "string", "default": ""},
            {"name": "task_id", "type": "string"},#任务，用于导主表和更新状态
            {"name": "project_sn", "type": "string"},#项目sn号，用于导表和更新状态
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.spring = self.add_module("metagbin.spring_files")
        self.assembly = self.add_module("metagbin.bin_assembly")
        self.predict = self.add_module("metagbin.stat_predict")
        self.step.add_steps('assembly','predict')
        self.remote_dir = self._sheet.output

    def check_options(self):
        """
        参数检查
        :return:
        """
        self.logger.info('正在进行参数检查')
        if not self.option("ref_fa").is_set:
            raise OptionError("请提供参考序列的文件夹！")
        if not self.option("bam_file").is_set:
            raise OptionError("请提供组装的序列文件！")
        if not self.option("bin_id"):
            raise OptionError("请提供bin的名称！")
        if not self.option("taxon"):
            raise OptionError("请提供参考的bin为细菌还是古菌！")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_spring(self):
        self.spring.set_options({
            "file_dir": self.option("bam_file").prop['path'],
            "type": "decompress"
        })
        self.spring.on("end", self.run_assembly)
        self.spring.run()

    def run_assembly(self):
        """
        对上传序列进行组装
        :return:
        """
        self.logger.info('开始对序列进行组装')
        bin_id = '_'.join(self.option('bin_id').split("_")[1:])
        if re.search(r'.gz', self.option('bam_file').prop['path']):
            opts = ({
                'ref_fa': self.option('ref_fa'),#参考bin的文件夹
                'bin_id': str(bin_id), #去掉前面的G_
                'bam_file': self.option('bam_file'),#一般为压缩后的文件夹
                'insert_size': str(self.option('insert_size')),
                'max_length': str(self.option('max_length')),
                'workflow': "interaction"
            })
        else:
            opts = ({
                'ref_fa': self.option('ref_fa'),#参考bin的文件夹
                'bin_id': str(bin_id), #去掉前面的G_
                'fastq_dir': self.spring.output_dir,#cleandata_dir
                'insert_size': str(self.option('insert_size')),
                'max_length': str(self.option('max_length')),
                'workflow': "workflow"
            })

        self.assembly.set_options(opts)
        self.assembly.on('start', self.set_step, {'start': self.step.assembly})
        self.assembly.on('end', self.set_step, {'end': self.step.assembly})
        self.assembly.on('end', self.set_step, {'start': self.step.predict})
        self.assembly.on('end', self.run_predict)
        self.assembly.run()

    def run_predict(self):
        """
        对组装结果进行预测
        :return:
        """
        self.logger.info('开始对组装结果进行基因预测')
        bin_id = '_'.join(self.option('bin_id').split("_")[1:])
        opts = ({
            'genome': self.assembly.option('scaffold'),
            'taxon' : self.option('taxon'),
            'bin_id': str(bin_id),
            'fastq1': self.assembly.option('fastq1'),
            'fastq2': self.assembly.option('fastq2'),
            'fastqs': self.assembly.option('fastqs'),
        })
        self.predict.set_options(opts)
        self.predict.on('start', self.set_step, {'start': self.step.predict})
        self.predict.on('start', self.set_step, {'end': self.step.predict})
        self.predict.on('end', self.set_output)
        self.predict.run()

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.genome_id = str(self.option('bin_id')).replace(".", "_")
        if os.path.exists(self.output_dir + '/Assembly_predict'):
            shutil.rmtree(self.output_dir + '/Assembly_predict')
        os.mkdir(self.output_dir + '/Assembly_predict')
        result_path = self.output_dir + '/Assembly_predict'
        link_dir(self.predict.output_dir, result_path)
        if not os.path.exists(result_path + "/Assembly/Assembly_summary/Bin_coverage.xls"):
            os.link(self.assembly.output_dir + '/Bin_coverage.xls', result_path + "/Assembly/Assembly_summary/Bin_coverage.xls")
        self.set_db()

    def remove_file(self):
        """
        删除文件
        :return:
        """
        if os.path.exists(self.output_dir + '/Assembly_predict/Assembly/Assembly_summary/Bin_coverage.xls'):
            os.remove(self.output_dir + '/Assembly_predict/Assembly/Assembly_summary/Bin_coverage.xls')
        if os.path.exists(self.output_dir + '/Assembly_predict/Assembly/Assembly_summary/Genome_coverage.xls'):
            os.remove(self.output_dir + '/Assembly_predict/Assembly/Assembly_summary/Genome_coverage.xls')
        if os.path.exists(self.output_dir + '/Assembly_predict/Gene_predict/length_distribute.txt'):
            os.remove(self.output_dir + '/Assembly_predict/Gene_predict/length_distribute.txt')

    def set_db(self):
        """
        将结果导入mongo库
        :return:
        """
        self.logger.info('正在写入mongo数据库')
        #self.genome_id = str(self.option('bin_id'))
        genome_id = self.genome_id
        old_genome_id = str(self.option('bin_id'))
        api_assembly = self.api.api('metagbin.assembly')
        api_genome_path = self.api.api("metagbin.common_api")
        assembly_id = self.option("main_id")
        soft_path = self.assembly.output_dir + "/Assembly_soft.txt"
        with open(soft_path, 'rb') as f:
            lines = f.readlines()
            assem_soft = lines[1].strip("\r\n").split('\t')[0]
        assembly_params = {'genome_id': genome_id,
                           "soft": assem_soft,
                            "submit_location": self.option("submit_location"),
                            "task_type": int(self.option("task_type"))}
        assembly_name = "Assembly_" + genome_id+ "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        #api_assembly.add_assembly(genome_id,  params=assembly_params, name=assembly_name)

        file_path = self.output_dir + "/Assembly_predict/Assembly/Assembly_summary/"+ old_genome_id +"_assembly_summary.xls"
        bin_cover_path1 = self.assembly.output_dir + "/Bin_coverage.xls"
        #insert_path = self.assembly.output_dir + "/config.txt"
        assem_cover_path2 = self.output_dir + "/Assembly_predict/Assembly/Assembly_summary/Genome_coverage.xls"
        genome_path = self.remote_dir + "Assembly/Assembly_summary/" + old_genome_id +"_scaffold.fa"
        anno_path = self.output_dir + '/Assembly_predict/Assembly/Assembly_summary/assembly.anno.xls'
        api_assembly.add_assembly_detail(assembly_id, genome_id, self.option('insert_size'), self.option('max_length'), file_path=file_path,bin_cover_path1=bin_cover_path1, assem_cover_path2=assem_cover_path2,
                                                      genome_path=genome_path, table_name="assembly", name=assembly_name, params=assembly_params, stat=anno_path)
        self.logger.info('组装结果统计写入mongo成功')

        api_assembly_assess = self.api.api('metagbin.assembly_assess')
        soft_path = self.predict.output_dir + '/Assembly/Genome_assessment/Assess_soft.txt'
        with open(soft_path, 'rb') as f2:
            lines = f2.readlines()
            assess_soft = lines[1].strip("\r\n").split('\t')[0]
        taxon_table = self.output_dir + '/Assembly_predict/Assembly/Assembly_summary/assembly.anno.xls'
        with open(taxon_table, 'r') as f:
            lines = f.readlines()
            assess_taxon = lines[1].strip('\r\n').split('\t')[-1]
        assess_path = self.output_dir + '/Assembly_predict/Assembly/Genome_assessment/'+old_genome_id+'_assess.xls'
        assess_params = {'genome_id': genome_id,
                         "soft": assess_soft,
                         "submit_location": "assemble_assess",
                         "task_type":self.option("task_type")}
        assess_name = "Assembly_name_" + genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        assess_id = api_assembly_assess.add_assembly_assess(genome_id, params=assess_params, name= assess_name)
        api_genome_path.add_sg_status("sg_status", submit_location="assemble_assess", params=assess_params,
                                      table_id=assess_id, table_name=assess_name, genome_id=genome_id, type_name="assembly_assess")
        api_assembly_assess.add_assembly_assess_detail(assess_id, assess_soft, genome_id, file_path=assess_path)
        self.logger.info('组装评估结果写入mongo成功')

        api_gc_depth = self.api.api("metagbin.assess_gc")
        gc_params = {"genome_id": genome_id,
                    "soft": "Bowtie2",
                    "submit_location": "assemble_gc",
                    "task_type":self.option("task_type")}
        gc_path = self.remote_dir + "Assembly/GC_depth/"
        gc_name = "GC_depth_" + genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        gc_id = api_gc_depth.add_assembly_gc(genome_id, gc_path=gc_path, params=gc_params, name=gc_name)
        api_genome_path.add_sg_status("sg_status", submit_location="assemble_gc", params=gc_params,
                                      table_id=gc_id, table_name=gc_name, genome_id=genome_id, type_name="assembly_gc")
        self.logger.info('组装Gc_depth结果写入mongo成功')

        api_predict_cds = self.api.api('metagbin.predict_cds')
        cds_soft = "Glimmer"
        cds_params = {"genome_id": genome_id,
                    "soft": cds_soft,
                    "submit_location": "predict_cds",
                    "task_type":self.option("task_type")}
        cds_name = "Predict_cds_" + genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        pre_id = api_predict_cds.add_predict_cds(genome_id, params=cds_params, name=cds_name)
        api_genome_path.add_sg_status("sg_status", submit_location="predict_cds", params=cds_params,
                                      table_id=pre_id, table_name=cds_name, genome_id=genome_id, type_name="predict_cds")
        sample_stat = self.output_dir + "/Assembly_predict/Gene_predict/CDS_predict/"+old_genome_id+"_CDS_statistics.xls"
        api_predict_cds.add_predict_cds_stat(pre_id, genome_id, sample_stat)
        predict_cds_gff = self.output_dir + "/Assembly_predict/Gene_predict/CDS_predict/"+old_genome_id + "_CDS.gff"
        fnn_path = self.remote_dir + "Gene_predict/CDS_predict/"+old_genome_id + "_CDS.fnn"
        faa_path = self.remote_dir + "Gene_predict/CDS_predict/"+old_genome_id + "_CDS.faa"
        gff_path = self.remote_dir + "Gene_predict/CDS_predict/"+old_genome_id + "_CDS.gff"
        api_predict_cds.add_predict_cds_detail(pre_id, genome_id, predict_gff=predict_cds_gff, fnn_path=fnn_path, faa_path=faa_path, gff_path=gff_path)
        distribute = self.output_dir +"/Assembly_predict/Gene_predict/length_distribute.txt"
        api_predict_cds.add_predict_cds_bar(pre_id, genome_id, distribute)
        self.logger.info('基因预测cds结果写入mongo成功')

        api_predict_repeat = self.api.api('metagbin.predict_repeat')
        repeat_soft = "Tandem Repeats Finder"
        repeat_params ={"genome_id": genome_id,
                        "parameters": "Match=2, Mismatch=7, Delta=7, PM=80, PI=10, Minscore=80, MaxPeriod=500",
                        "soft": repeat_soft,
                        "submit_location": "predict_repeat",
                        "task_type":self.option("task_type")}
        repeat_name = "Predict_repeat_" + genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        repeat_id = api_predict_repeat.add_predict_repeat(genome_id, params=repeat_params, name=repeat_name)
        api_genome_path.add_sg_status("sg_status", submit_location="predict_repeat", params=repeat_params,
                                      table_id=repeat_id, table_name=repeat_name, genome_id=genome_id, type_name="predict_repeat")
        predict_repeat_gff = self.output_dir + "/Assembly_predict/Gene_predict/repeats/"+old_genome_id + "_TRF.gff"
        if os.path.exists(predict_repeat_gff):
            api_predict_repeat.add_predict_repeat_detail(repeat_id, genome_id, predict_gff=predict_repeat_gff)
            self.logger.info('重复序列预测结果写入mongo成功')
        else:
            self.logger.error('未能成功预测出重复序列')

        api_predict_rrna = self.api.api('metagbin.predict_rrna')
        rrna_soft = "Barrnap"
        rrna_params = {"genome_id": genome_id,
                       "soft": rrna_soft,
                        "submit_location": "predict_rrna",
                        "task_type":self.option("task_type")}
        rrna_name = "Predict_rrna" + genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        rrna_id = api_predict_rrna.add_predict_rrna(genome_id, params=rrna_params, name=rrna_name)
        api_genome_path.add_sg_status("sg_status", submit_location="predict_rrna", params=rrna_params,
                                      table_id=rrna_id, table_name=rrna_name, genome_id=genome_id, type_name="predict_rrna")
        predict_rna_gff = self.output_dir + "/Assembly_predict/Gene_predict/rRNA/"+old_genome_id + "_rRNA.gff"
        rna_path = self.output_dir + "/Assembly_predict/Gene_predict/rRNA/"+ old_genome_id +'-16S_rRNA.fa'
        new_rna_path = self.remote_dir + "Gene_predict/rRNA/" + old_genome_id +'-16S_rRNA.fa'
        fnn = self.output_dir + "/Assembly_predict/Gene_predict/rRNA/"+old_genome_id + "_rRNA.fnn"
        if os.path.exists(predict_rna_gff):
            api_predict_rrna.add_predict_rrna_detail(rrna_id, genome_id, predict_gff=predict_rna_gff)
            if os.path.exists(rna_path):
                api_predict_rrna.add_predict_rrna_seq(rrna_id, genome_id, fnn=fnn, predict_gff=new_rna_path)
            else:
                api_predict_rrna.add_predict_rrna_seq(rrna_id, genome_id, fnn=fnn)
                self.logger.error("未能成功预测出16s_rRNA")
        else:
            self.logger.error('未能成功预测出rRNA')
        self.logger.info('rRNA预测结果写入mongo成功')

        api_predict_trna = self.api.api('metagbin.predict_trna')
        trna_soft = "tRNAscan-SE"
        trna_params = {"genome_id": genome_id,
                       "soft": trna_soft,
                        "submit_location": "predict_trna",
                        "task_type":self.option("task_type")}
        trna_name = "Predict_trna_" + genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        trna_id = api_predict_trna.add_predict_trna(genome_id, params=trna_params, name=trna_name)
        api_genome_path.add_sg_status("sg_status", submit_location="predict_trna", params=trna_params,
                                      table_id=trna_id, table_name=trna_name, genome_id=genome_id, type_name="predict_trna")
        predict_trna_gff = self.output_dir + "/Assembly_predict/Gene_predict/tRNA/"+old_genome_id + "_tRNA.gff"
        if os.path.exists(predict_trna_gff):
            api_predict_trna.add_predict_trna_detail(trna_id, genome_id, predict_gff=predict_trna_gff)
        else:
            self.logger.error('未能成功预测出tRNA')
        self.logger.info('tRNA预测结果写入mongo成功')

        if self.option('taxon') == 'Bacteria':
            with open(assess_path, 'rb') as f8:
                lines = f8.readlines()
                complete = float(lines[1].strip("\r\n").split('\t')[0])
        else:
            with open(assess_path, 'rb') as f9:
                lines = f9.readlines()
                complete = float(lines[1].strip("\r\n").split('\t')[0])
        genome_params={"genome_id": genome_id, "complete": complete}
        old_task_id = str(self._sheet.id).split('_')[0]
        num = str(self._sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        genome_name = "Genome_" + genome_id + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        api_genome_path.add_main("genome", name=genome_name, params=genome_params,
                                 others={"bin_id": '_'.join(self.option('bin_id').split("_")[1:]), "genome_id":genome_id, "taxon": assess_taxon},
                                 task_id=task_id)
        #api_genome_path.add_sg_status("genome", submit_location=self.option("submit_location"), params=genome_params,table_id=genome_manage_id, table_name="genome", genome_id=genome_id)
        self.logger.info('基因组管理表导入mongo成功')
        self.end()

    def end(self):
        self.logger.info('正在上传分析结果')
        result_dir = self.output_dir + "/Assembly_predict"
        self.logger.info('正在删除文件')
        self.remove_file()
        repaths = [
            [".", "", "基因组交互分析结果目录", 0, ""],
            ["Assembly", "", "基因组组装结果目录", 0, ""],
            ["Gene_predict", "", "基因组的基因预测目录", 0, ""],
            ["Assembly/Assembly_summary", "", "%s基因组组装统计结果目录" % self.genome_id,0,''],
            ["Assembly/GC_depth", "", "%s基因组GC_depth" % self.genome_id,0, ''],
            ["Assembly/Genome_assessment", "", "%s基因组评估结果目录" % self.genome_id,0,''],
            ["Gene_predict/CDS_predict", "", "%s基因组的编码基因预测目录" % self.genome_id,0,''],
            ["Gene_predict/tRNA", "", "%s基因组的tRNA预测结果目录" % self.genome_id,0,''],
            ["Gene_predict/rRNA", "", "%s基因组的rRNA预测结果目录" % self.genome_id,0,''],
            ["Gene_predict/repeats", "", "%s基因组的重复序列预测结果目录" % self.genome_id,0,''],
        ]
        regexps = [
            [r"/Assembly/.+_assembly_scaffold_details.xls", "xls", "%s基因组组装scaffold统计结果"% self.genome_id,0,''],
            [r"/Assembly/.+_scaffold.fa", "", "%s基因组组装结果"% self.genome_id,0,''],
            #[r"/Assembly/Bin_coverage.xls", "xls", "%s基因组组装的binning覆盖度结果"% self.genome_id,0,''],
            #[r"/Assembly/Genome_coverage.xls", "xls", "%s基因组组装的genome覆盖度结果"% self.genome_id,0,''],
            [r"/Assembly/Genome_assessment/.+_assess.xls", "xls", "%s基因组评估结果表"% self.genome_id,0,''],
            [r"/Assembly/GC_depth/gc_depth.png", "", "%s基因组GC_depth_1000PNG统计图"% self.genome_id,0,''],
            [r"/Assembly/GC_depth/gc_depth.svg", "", "%s基因组GC_depth_1000SVG统计图"% self.genome_id,0,''],
            [r"/Genomic_assessment/predict/CDS_predict/.+_CDS.faa", "", "%s基因组的编码基因预测蛋白序列"% self.genome_id,0,''],
            [r"/Gene_predict/CDS_predict/.+_CDS.fnn", "", "%s基因组的编码基因预测核苷酸序列"% self.genome_id,0,''],
            [r"/Gene_predict/CDS_predict/.+_CDS.gff", "", "%s基因组的编码基因预测gff格式统计结果"% self.genome_id,0,''],
            [r"/Gene_predict/CDS_predict/.+_CDS_statistics.xls", "xls", "%s基因组的编码基因预测统计表"% self.genome_id,0,''],
            [r"/Gene_predict/CDS_predict/length_distribute.txt", "", "%s编码基因预测长度分布表"% self.genome_id,0,''],
            [r"/Gene_predict/repeats/.+/.+_TRF.gff", "", "%s基因组的重复序列预测gff格式统计结果"% self.genome_id,0,''],
            [r"/Gene_predict/rRNA/.+.16S_rRNA.fa", "", "%s基因组的rRNA预测16S核苷酸序列"% self.genome_id,0,''],
            [r"/Gene_predict/rRNA/.+_rRNA.gff", "", "%srRNA预测gff统计文件"% self.genome_id,0,''],
            [r"/Gene_predict/rRNA/.+_rRNA.fnn", "", "%s基因组的rRNA预测核苷酸序列"% self.genome_id,0,''],
            [r"/Gene_predict/rRNA/.+_rRNA.gff", "", "%s基因组的rRNA预测gff格式统计结果"% self.genome_id,0,''],
            [r"/Gene_predict/tRNA/.+_tRNA.gff", "", "%s基因组的tRNA预测gff格式统计结果"% self.genome_id,0,''],
            [r"/Gene_predict/tRNA/.+_tRNA.fnn", "", "%s基因组的tRNA预测核苷酸序列"% self.genome_id,0,''],
        ]
        resultdir = self.add_upload_dir(result_dir)
        resultdir.add_relpath_rules(repaths)
        resultdir.add_relpath_rules(regexps)
        super(AssemblyPredictWorkflow, self).end()

    def run(self):
        if re.search(r'.gz', self.option('bam_file').prop['path']):
            self.run_assembly()
        else:
            self.run_spring()
        super(AssemblyPredictWorkflow, self).run()




