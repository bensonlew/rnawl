# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'

""" 从文件中获取数据统计信息"""
import os
import datetime
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId


class DataStatWorkflow(Workflow):
    """
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DataStatWorkflow,self).__init__(wsheet_object)
        options = [
            #{"name": "data_table", "type": "infile", "format": "sequence.profile_table"}
            {'name': 'test', 'type': 'bool', 'default': False}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        

    def check_options(self):

        pass

    def run(self):
        self.set_db()
        super(DataStatWorkflow, self).run()

    def set_db(self):
        """
        将运行结果保存到mongo数据库中
        :return:
        """
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.api_stat = self.api.api("bac_comp_genome.gene_stat")
        self.api_gene = self.api.api("bac_comp_genome.gene_predict")
        gene_predict = "/mnt/ilustre/users/sanger-dev/home/zhangqingchen/compare_genome/compare_result/compare_test/module/total_sample.xls"
        self.seq_type = {}
        with open(gene_predict, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('#'):
                    pass
                else:
                    line = line.strip().split("\t")
                    sample = line[0]
                    type = line[1]
                    self.seq_type[sample] = type
        main_id = "5d24d7d717b2bf6bb03cf058"
        params = {
            "soft":{
                "cds": "Glimmer, GeneMark, Prodigal",
                "rrna": "Barrnap",
                "trna": "tRNAscan-SE",
            }
        }
        dir_path = "/mnt/ilustre/users/sanger-dev/workspace/20191028/Single_gene_predict_total_1238/GenePredict/output"
        stat_id = self.api_stat.add_gene_stat(params=params)
        predict_id = self.api_gene.add_gene_predict(params=params)
        data_id = "5db68d3717b2bf627d8f7b6d" ###用于更新样本信息表的信息
        genomes_id = "5db68d3717b2bf627d8f7b6d" ###用于更新基因组基本特征信息表
        for sample in os.listdir(dir_path):
            type = self.seq_type[sample]
            sample_dir_path = os.path.join(dir_path, sample)
            sample_path = os.path.join(sample_dir_path, sample + '_sample_stat.xls')
            self.api_stat.add_gene_stat_detail(stat_id, sample, sample_path)
            self.logger.info("导入统计表格成功！")

            self.api_stat.add_genomes_detail(data_id, sample, sample_path)
            if type in ["seq", "gff"]:##参考库中的数据不再更新，seq或者seq+gff这两种格式的数据需要再重新更新表
                sample_stat_path = os.path.join(sample_dir_path, sample + '_stat.xls') #需要对其进行进一步的处理
                self.api_stat.add_data_detail(genomes_id, sample, sample_stat_path)

            self.logger.info("更新样本信息表和基因组特征表成功！")

            cds_path = os.path.join(sample_dir_path, sample + '_CDS.gff')
            self.api_gene.add_gene_cds_detail(predict_id, sample, cds_path)
            self.logger.info("导入cds预测结果成功！")

            rrna_path = os.path.join(sample_dir_path, sample + '_rRNA.gff')
            s16_path = os.path.join(sample_dir_path, sample + '_16S.fna')
            if os.path.exists(s16_path):
                self.api_gene.add_gene_rrna_detail(predict_id, sample, rrna_path, fnn=s16_path)
            else:
                self.api_gene.add_gene_rrna_detail(predict_id, sample, rrna_path)
            self.logger.info("导入rrna预测结果成功！")

            trna_path = os.path.join(sample_dir_path, sample + '_tRNA.gff')
            self.api_gene.add_gene_trna_detail(predict_id, sample, trna_path, type)
            self.logger.info("导入trna预测结果成功！")


        self.end()

    def end(self):
        super(DataStatWorkflow, self).end()


