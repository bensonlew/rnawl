# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
import gevent
import shutil
from mbio.packages.bacgenome.common import link_dir


class IntegronWorkflow(Workflow):
    """
    细菌基因组整合子预测
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(IntegronWorkflow, self).__init__(wsheet_object)
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


    def run_integron(self):
        """
        整合子预测
        :return:
        """
        self.logger.info("开始用integron_finder软件进行预测！")
        self.samples = str(self.option("specimen_id")).split(",")
        self.logger.info(self.samples)
        self.download_file()
        for sample in self.samples:
            self.integron_predict = self.add_module('bacgenome.integron_predict')
            if self.analysis_type in ['uncomplete']:
                options = {
                    "genome_fa": os.path.join(self.assemble_dir, sample +"_scaf.fna"),
                    "gene_faa": os.path.join(self.gene_dir, sample, sample + ".faa"),
                    "gene_gff": os.path.join(self.gene_dir, sample, sample + ".gff"),
                    "sample": sample
                }
            else:
                options = {
                    "genome_fa": os.path.join(self.work_dir, "assemble_dir2", sample +"_scaf.fna"),
                    "gene_faa": os.path.join(self.gene_dir, sample, sample + ".faa"),
                    "gene_gff": os.path.join(self.gene_dir, sample,sample + ".gff"),
                    "sample": sample
                }
            self.integron_predict.set_options(options)
            self.modules.append(self.integron_predict)
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
        self.run_integron()
        super(IntegronWorkflow, self).run()

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
                        link_dir(sample_dir_path, output_sample)
        self.logger.info("设置结果文件目录完成!")
        self.set_db()


    def set_db(self):
        self.logger.info("start mongo>>>>>>>>>>>>")
        api_integron = self.api.api('bacgenome.anno_integron')
        for sample in self.samples:
            self.logger.info(sample)
            main_id = self.option("main_id")
            sample_dir = os.path.join(self.output_dir, sample)
            sample_path = os.path.join(sample_dir, sample + ".sample.xls")
            # params = {"software":"Integron Finder version 1.5"}
            # main_id = api_integron.add_anno_integron(task_id="tsg_37312", project_sn="188_5eb378c64a3f6",params=params,specimen_id="GCF_000016585.1")
            stat_file = os.path.join(sample_dir, sample + ".stat.xls")
            seq_file = os.path.join(sample_dir, sample + ".integron.fna")
            if os.path.exists(sample_path):
                api_integron.add_anno_integron_stat(sample_path, main_id=main_id)
            else:
                api_integron.add_main_id(main_id=main_id)
            if os.path.exists(stat_file) and os.path.exists(seq_file):
                api_integron.add_anno_integron_detail(stat_file,seq_file, main_id=main_id)
            summary_path = os.path.join(sample_dir, sample + ".integrons")
            sequence_path = os.path.join(sample_dir, sample + "_sequence.fna")
            if os.path.exists(summary_path) and os.path.exists(sequence_path):
                api_integron.add_anno_integron_summary(summary_path, sequence_path, sample, main_id=main_id)
            if os.path.exists(stat_file):
                os.rename(stat_file, os.path.join(sample_dir, sample + "_Integron_detail.xls"))
            if os.path.exists(sample_path):
                os.rename(sample_path, os.path.join(sample_dir, sample + "_stat.xls"))
        self.end()

    def end(self):
        repaths = [
            [".", "", "整合子预测结果目录",0,"130192"],
        ]
        regexps = [
            [r'.*_stat.xls', 'xls', '整合子分析统计表',0,"130193"],
            [r'.*_Integron_detail.xls', 'xls', '整合子分析详情表',0,"130194"],
            [r'.*.integrons', 'txt', '整合子元素表',0,"130195"],
            [r'.*.summary', 'txt', '整合子分类统计表',0,"130196"],
            [r'.*.integron.fna', 'fna', '整合子序列文件',0,"130197"],
            [r'.*_sequence.fna', 'fna', '整合子元素序列文件',0,"130198"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(IntegronWorkflow, self).end()
