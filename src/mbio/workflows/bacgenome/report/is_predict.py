# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
import types
import shutil
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
import gevent
from mbio.packages.bacgenome.common import link_dir


class IsPredictWorkflow(Workflow):
    """
    细菌基因组插入序列分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(IsPredictWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimen_id", "type": "string"},  # 样品名称
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.common_api = self.api.api('bacgenome.common_api')
        self.modules = []

    def download_file(self):
        """
        download file from s3
        :return:
        """
        self.logger.info("开始下载组装序列文件！")
        assemble_dir = os.path.join(self.work_dir, "assemble_dir")
        if os.path.exists(assemble_dir):
            shutil.rmtree(assemble_dir)
        (self.assemble_dir, self.analysis_type)= self.common_api.down_seq_files(assemble_dir, self.option("task_id"), self.samples)
        self.logger.info("开始下载基因序列文件！")
        gene_dir = os.path.join(self.work_dir, "gene_dir")
        if os.path.exists(gene_dir):
            shutil.rmtree(gene_dir)
        self.gene_dir = self.common_api.down_gene_files(gene_dir, self.option("task_id"), self.samples)
        if self.analysis_type in ['complete']:
            self.cat_reads()

    def cat_reads(self):
        """
        合并reads 完成图合并fasta文件，形成单个样本文件一个
        :return:
        """
        assemble_dir2 = os.path.join(self.work_dir, "assemble_dir2")
        if os.path.exists(assemble_dir2):
            shutil.rmtree(assemble_dir2)
        os.mkdir(assemble_dir2)
        for sample in self.samples:
            sample_path = os.path.join(assemble_dir2, sample + "_scaf.fna")
            assemble_dir = os.path.join(self.assemble_dir, sample)
            dir_list = os.listdir(assemble_dir)
            for file2 in dir_list:
                os.system("cat {} >> {}".format(os.path.join(assemble_dir, file2), sample_path))

    def run_is(self):
        """
        整合子预测
        :return:
        """
        self.logger.info("开始用is软件进行预测！")
        self.samples = str(self.option("specimen_id")).split(",")
        self.download_file()
        for sample in self.samples:
            self.is_predict = self.add_module('bacgenome.is_predict')
            if self.analysis_type in ['uncomplete']:
                options = {
                    "genome_fa": os.path.join(self.assemble_dir, sample +"_scaf.fna"),
                    "gene_faa": os.path.join(self.gene_dir, sample, sample + ".faa"),
                    "gene_fna": os.path.join(self.gene_dir, sample, sample + ".fnn"),
                    "gene_gff": os.path.join(self.gene_dir, sample, sample + ".gff"),
                    "sample": sample
                }
            else:
                options = {
                    "genome_fa": os.path.join(self.work_dir,"assemble_dir2", sample +"_scaf.fna"),
                    "gene_faa": os.path.join(self.gene_dir, sample, sample + ".faa"),
                    "gene_fna": os.path.join(self.gene_dir, sample, sample + ".fnn"),
                    "gene_gff": os.path.join(self.gene_dir, sample, sample + ".gff"),
                    "sample": sample
                }
            self.is_predict.set_options(options)
            self.modules.append(self.is_predict)
        if len(self.modules) > 1:
            self.on_rely(self.modules, self.set_output)
        else:
            self.modules[0].on('end', self.set_output)
        for module in self.modules:
            module.run()
            gevent.sleep(0)

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_is()
        super(IsPredictWorkflow, self).run()

    def set_output(self):
        """
        设置结果目录文件
        :return:
        """
        self.logger.info("开始设置结果文件目录！")
        for module in self.modules:
            for sample_dir in os.listdir(module.output_dir):
                sample_dir_path = module.output_dir + "/" +sample_dir
                output_sample = self.output_dir + "/" +sample_dir
                if os.path.exists(output_sample):
                    shutil.rmtree(output_sample)
                if os.path.exists(sample_dir_path):
                    if len(os.listdir(sample_dir_path)) != 0:
                        self.logger.info(sample_dir_path)
                        link_dir(sample_dir_path, output_sample)
        self.logger.info("设置结果文件目录完成!")
        self.set_db()

    def set_db(self):
        """
        导入MongoDB数据
        :return:
        """
        self.logger.info("start mongo>>>>>>>>>>>>")
        api_is = self.api.api('bacgenome.anno_is')
        main_id = self.option("main_id")
        for sample in self.samples:
            sample_dir = os.path.join(self.output_dir, sample)
            sample_path = os.path.join(sample_dir, sample + ".stat.xls")
            if os.path.exists(sample_path):
                api_is.add_anno_is_stat(sample_path, main_id=main_id, sample=sample)
            else:
                api_is.add_main_id(main_id=main_id)
            stat_file = os.path.join(sample_dir, sample + ".raw")
            seq_file = os.path.join(sample_dir, sample + "_sequence.fna")
            gff_file = os.path.join(sample_dir, sample + ".is.xls")
            enzyme_file = os.path.join(sample_dir, sample + ".enzyme.xls")
            gene_file = os.path.join(sample_dir, sample + ".gene.ffn")
            blast_file = os.path.join(sample_dir, sample + ".blast.xls")
            if os.path.exists(gff_file):
                if os.path.exists(enzyme_file) and os.path.exists(gene_file):
                    api_is.add_anno_is_detail(stat_file,seq_file, main_id=main_id, sample=sample, gff_file=gff_file, transposon_file=enzyme_file,blast_file=blast_file)
                    api_is.add_anno_is_summary(gff_file, seq_file, sample, main_id=main_id,transposon_file=enzyme_file, gene_file=gene_file)
                else:
                    api_is.add_anno_is_detail(stat_file,seq_file, main_id=main_id, sample=sample, blast_file=blast_file)
                    api_is.add_anno_is_summary(gff_file, seq_file, sample, main_id=main_id)
                os.remove(gene_file)
                os.remove(stat_file)
                os.remove(os.path.join(sample_dir, sample + ".sum"))
                os.remove(os.path.join(sample_dir, sample +  ".gff"))
                os.remove(os.path.join(sample_dir, sample +  ".orf.faa"))
                os.remove(os.path.join(sample_dir, sample +  ".orf.fna"))
                os.remove(os.path.join(sample_dir, sample +  "_sequence.fna"))
        self.logger.info("end MongoDB<<<<<<<<<<<<<")
        self.end()

    def end(self):
        repaths = [
            [".", "", "插入序列预测结果目录",0,"130199"],
        ]
        regexps = [
            [r'.*_stat.xls', 'xls', '插入序列预测统计表',0,"130200"],
            [r'.*.blast.xls', 'xls', '插入序列预测blast结果表',0,"130201"],
            [r'.*.is.xls', 'xls', '插入序列预测结果文件',0,"130202"],
            [r'.*.out', 'xls', '插入序列预测详情表',0,"130203"],
            [r'.*.enzyme.xls', 'xls', '插入序列预测转座酶详情表',0,"130204"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(IsPredictWorkflow, self).end()