#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'fwy'

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna_v3.functions import workfuncdeco
from  mbio.packages.ref_rna_v3.large.wsheet import Sheet
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
import json
import glob
import os
import pandas as pd
import collections
from collections import OrderedDict
import re
import dask
import unittest
import datetime
import gevent
import functools
import time


# 定义用于统计导表时间的装饰器
def time_count(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run ' + func_name + ' at ' + start_time)
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End ' + func_name + ' at ' + end_time)
        print("{}函数执行时间约为{}s".format(func.__name__, end - start))

    return wrapper


class SetDbAgent(Agent):
    """
    SetDb:用于导表,将导表函数加入
    """

    def __init__(self, parent):
        super(SetDbAgent, self).__init__(parent)
        options = [
            {"name": "sheet_data_json", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件
            {"name": "option_data_json", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件
            {'name': 'analysis_content', 'type': 'string'}, #输入
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'is_assemble', 'type': 'bool', 'default': True},
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'productive_table', 'type': 'infile', 'format': 'sample.group_table'},
            # 质控序列文件
            # {'name': 'qc_dir', 'type': 'string', 'default': None},
            {'name': 'qc_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},
            {'name': 'group', 'type': 'infile', 'format': 'denovo_rna_v2.group_table'},
            {'name': 'function_json', "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'control', 'type': 'infile', 'format': 'denovo_rna_v2.compare_table'},


        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        # if not self.option("merge_file").is_set:
        #     raise OptionError("请输入VCF格式文件", code="35600802")
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        super(SetDbAgent, self).end()


class SetDbTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(SetDbTool, self).__init__(config)
        self.sheet = self.get_parent_sheet_data()
        self.workflow_options  = self.get_option_data()
        self.function_infos = self.get_function_json()
        self.group_id = ""
        self.group_detail = ""
        self.group_category = ""
        self.control_id = ""
        self.compare_detail = ''
        self.task_id = self.option("task_id")
        self.project_sn =  self.sheet.project_sn
        self.exp_ids = dict()
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        
        #整理接口






    def get_parent_sheet_data(self):
        with open(self.option("sheet_data_json").prop["path"], 'r') as f:
            a = json.loads(f.read())
        b = Sheet(data = a)
        return b

    def get_option_data(self):
        with open(self.option("option_data_json").prop["path"], 'r') as f:
            a = json.loads(f.read())
        return a

    def get_function_json(self):
        with open(self.option("function_json").prop["path"], 'r') as f:
            a = json.loads(f.read())
        return a

    def workflow_option(self, name, value=None):
        if name not in self.workflow_options.keys():
            raise Exception("参数%s不存在，请先添加参数" % name)
        if value is None:
            return self.workflow_options[name]
        else:
            self.workflow_options[name] = value

    def build_seq_database(self):
        self.export_seq = self.api.api("denovo_rna_v2.seq_detail")
        cds = self.function_infos["build_seq_database"]["cds"]
        pep = self.function_infos["build_seq_database"]["pep"]
        fasta = self.function_infos["build_seq_database"]["fasta"]
        trans2unigene = self.function_infos["build_seq_database"]["trans2unigene"]
        seq_db = self.function_infos["build_seq_database"]["seq_db"]
        self.export_seq.build_seq_database(seq_db, cds, pep, fasta, trans2unigene, task_id=self.task_id)
        self.seq_stat_seq = self.api.api("denovo_rna_v2.seq_stat_detail")
        self.logger.info("开始进行sg_seq_stat的导表")
        gene_stat = self.function_infos["build_seq_database"]["gene_stat"]
        self.logger.info("gene_stat的路径是{}".format(gene_stat))
        trans_stat = self.function_infos["build_seq_database"]["trans_stat"]
        self.logger.info("trans_stat的路径是{}".format(trans_stat))
        self.seq_stat_seq.add_seq_stat(gene_stat, trans_stat)
        self.logger.info("开始进行sg_seq_stat的结束")

    def run_api(self, test=False):
        greenlets_list_first = []
        greenlets_list_sec = []
        # greenlets_list_third = []
        task_info = self.api.api('task_info.denovo_task_info')
        task_info.add_task_info()
        self.logger.info("进行第一阶段导表")
        greenlets_list_first.append(gevent.spawn(self.export_qc))
        greenlets_list_first.append(gevent.spawn(self.export_denovo_assembly))
        greenlets_list_first.append(gevent.spawn(self.export_denovo_align))
        gevent.joinall(greenlets_list_first)
        self.logger.info("进行第二阶段导表")
        greenlets_list_sec.append(gevent.spawn(self.export_denovo_annotation))
        if self.workflow_option("sample_num") == "multiple":
            if self.workflow_option("is_snp") == "True":
                greenlets_list_sec.append(gevent.spawn(self.export_snp))
        greenlets_list_sec.append(gevent.spawn(self.export_ssr))
        if self.workflow_option("tf_database") != "Other":
            greenlets_list_sec.append(gevent.spawn(self.export_tf))
        greenlets_list_sec.append(gevent.spawn(self.export_cds))
        gevent.joinall(greenlets_list_sec)
        self.logger.info("进行第三阶段导表")
        # greenlets_list_third.append(gevent.spawn(self.export_expression))
        # gevent.joinall(greenlets_list_third)
        self.export_expression()
        self.logger.info("导表完成")

    @time_count
    def export_qc(self):
        gevent.sleep()
        api_qc = self.api.api("denovo_rna_v2.denovo_rna_qc")
        qc_stat_before = self.function_infos["export_qc"]["qc_stat_before"]
        fq_type = self.workflow_option("fq_type").lower()
        if self.workflow_option("datatype") == "rawdata":
            api_qc.add_samples_info(qc_stat_before, fq_type=fq_type, about_qc="before", group=self.option('group').path)
        # api_qc.add_samples_alias(self.option("alias_table").prop["path"],about_qc="before")
        quality_stat_after = self.function_infos["export_qc"]["quality_stat_after"]  # 质控数据结果统计
        quality_stat_before = self.function_infos["export_qc"]["quality_stat_before"]  # 原始数据结果统计
        if self.workflow_option("datatype") == "rawdata":
            api_qc.add_gragh_info(quality_stat_before, "before")
        qc_stat_after = self.function_infos["export_qc"]["qc_stat_after"]
        api_qc.add_samples_info(qc_stat_after, fq_type=fq_type, about_qc="after", group=self.option('group').path)
        # api_qc.add_samples_alias(self.option("alias_table").prop["path"],about_qc="after")
        api_qc.add_gragh_info(quality_stat_after, "after")
        if self.option("group").is_set:
            self.group_id, self.group_detail, self.group_category = api_qc.add_specimen_group(
                self.option("group").prop["path"])
            self.logger.info("group_detail为：" + str(self.group_detail))
            if self.option('productive_table').is_set:
                api_qc.add_productive_name(samples=self.option('group').prop['sample'],
                                   productive_table=self.option('productive_table').path)
        if self.option("control").is_set:
            self.control_id, compare_detail = api_qc.add_control_group(self.option("control").prop["path"],
                                                                            self.group_id)
            self.compare_detail = compare_detail
        # if self.option("is_snp") == True:
        intermediate_dir = self.sheet.output.replace('workflow_results', 'intermediate_results/')
        api_qc.add_bam_path(intermediate_dir)

    @time_count
    def export_denovo_assembly(self):
        '''
        导入组装结果表格 liubinxu
        '''
        gevent.sleep()
        api_denovo_assembly = self.api.api("denovo_rna_v2.denovoass_assemble2")
        assemble_filter_dir =  self.function_infos["export_denovo_assembly"]["assemble_filter_dir"]
        filter_evolution = self.function_infos["export_denovo_assembly"]["filter_evolution"]
        unigene_evaluation = self.function_infos["export_denovo_assembly"]["unigene_evaluation"]
        filter_unigene_evaluation = self.function_infos["export_denovo_assembly"]["filter_unigene_evaluation"]
        # result_dir = self.assemble.output_dir
        self.logger.info("{}".format(self.workflow_option("optimize")))
        if self.workflow_option("optimize") == "True" or self.workflow_option("optimize") == True:
            api_denovo_assembly.run(assemble_filter_dir, unigene_evaluation, filter_evolution,
                                         filter_unigene_evaluation)
        else:
            api_denovo_assembly.run(assemble_filter_dir, unigene_evaluation)

    @time_count
    def export_denovo_align(self):
        '''
        导入mapping率统计表格 liubinxu
        '''
        gevent.sleep()
        self.api_denovo_align = self.api.api("denovo_rna_v2.denovo_align")
        result_dir = self.function_infos["export_denovo_align"]["result_dir"]
        self.api_denovo_align.run(result_dir, group=self.option('group').path)

    @time_count
    def export_denovo_annotation(self):
        '''
        导入注释结果表格 liubinxu
        '''
        gevent.sleep()
        api_denovo_annotation = self.api.api("denovo_rna_v2.denovo_annotation")
        result_dir = self.function_infos["export_denovo_annotation"]["result_dir"]
        trans2gene = self.function_infos["export_denovo_annotation"]["trans2gene"]
        exp_output = self.function_infos["export_denovo_annotation"]["exp_output"]
        gene_exp = os.path.join(exp_output, 'gene.tpm.matrix')
        trans_exp = os.path.join(exp_output, 'transcript.tpm.matrix')
        params = {
            "nr_evalue": str(self.workflow_option("nr_evalue")),
            "nr_similarity": str(0),
            "nr_identity": str(0),
            "swissprot_evalue": str(self.workflow_option("swissprot_evalue")),
            "swissprot_similarity": str(0),
            "swissprot_identity": str(0),
            "cog_evalue": str(self.workflow_option("string_evalue")),
            "cog_similarity": str(0),
            "cog_identity": str(0),
            "kegg_evalue": str(self.workflow_option("kegg_evalue")),
            "kegg_similarity": str(0),
            "kegg_identity": str(0),
            "pfam_evalue": str(self.workflow_option("pfam_evalue")),
            "submit_location": "annotationstat",
            "task_id": self.task_id,
            "task_type": 2,
        }
        api_denovo_annotation.anno_type = 'origin'
        api_denovo_annotation.run(result_dir, trans2gene, params, taxon=self.workflow_option("kegg_database"), version="v2",
                                       exp_level="T", gene_exp=gene_exp, trans_exp=trans_exp)
        # self.api_denovo_annotation.run(result_dir, trans2gene, params, taxon=self.option("kegg_database"))

    @time_count
    def export_snp(self):
        gevent.sleep()
        api_snpfinal = self.api.api("denovo_rna_v3.snp_api")
        self.logger.info("开始进行Snpfinal的导表")
        task_id = self.task_id
        project_sn = self.project_sn
        # call_vcf_path = self.snp.work_dir + '/Snp/' + 'call.vcf'
        params = dict(
            task_id=task_id,
            submit_location="snp_detail",
            task_type=2,
            method=self.workflow_option("snp_method")
        )
        snpfinal_work_dir = self.function_infos["export_snp"]["snpfinal_work_dir"]



        if snpfinal_work_dir + '/new_snp_rewrite' is None and snpfinal_work_dir + '/new_indel_rewrite' is None:
            self.logger.info("此次分析没有call出snp和indel")

        if snpfinal_work_dir + '/new_snp_rewrite' is not None and snpfinal_work_dir + '/new_indel_rewrite' is not None:
            new_snp_rewrite = snpfinal_work_dir + '/new_snp_rewrite'
            new_indel_rewrite = snpfinal_work_dir + '/new_indel_rewrite'
            depth_path = snpfinal_work_dir + '/depth_new_per'
            hh_path = snpfinal_work_dir + '/statis_hh'
            tt_new_per_path = snpfinal_work_dir + '/tt_new_per'
            cds_path = snpfinal_work_dir + "/statis_cds"
            anno_path = snpfinal_work_dir + "/snp_anno_stat"
            api_snpfinal.add_snp(new_snp_rewrite, new_indel_rewrite, depth_path, hh_path, tt_new_per_path, cds_path,
                                 anno_path, project_sn=project_sn, task_id=task_id, params=params,
                                 group=self.option('group').path)
            self.logger.info("Snpfinal的导表成功，此次分析call出snp和indel")

        if snpfinal_work_dir + '/new_snp_rewrite' is not None and snpfinal_work_dir + '/new_indel_rewrite' is None:
            new_snp_rewrite = snpfinal_work_dir + '/new_snp_rewrite'
            depth_path = snpfinal_work_dir + '/depth_new_per'
            hh_path = snpfinal_work_dir + '/statis_hh'
            tt_new_per_path = snpfinal_work_dir + '/tt_new_per'
            cds_path = snpfinal_work_dir + "/statis_cds"
            anno_path = snpfinal_work_dir + "/snp_anno_stat"
            api_snpfinal.add_snp(new_snp_rewrite=new_snp_rewrite, new_indel_rewrite=None, depth_path=depth_path,
                                 hh_path=hh_path, tt_new_per_path=tt_new_per_path, cds_path=cds_path,
                                 anno_path=anno_path, project_sn=project_sn, task_id=task_id, params=params,
                                 group=self.option('group').path)
            self.logger.info("Snpfinal的导表成功，此次分析call出snp,但是没有call出indel")

        if snpfinal_work_dir + '/new_snp_rewrite' is None and snpfinal_work_dir + '/new_indel_rewrite' is not None:
            new_indel_rewrite = snpfinal_work_dir + '/new_indel_rewrite'
            api_snpfinal.add_snp(new_snp_rewrite=None, new_indel_rewrite=new_indel_rewrite, project_sn=project_sn,
                                 task_id=task_id, params=params, group=self.option('group').path)
            self.logger.info("Snpfinal的导表成功，此次分析call出indel,但是没有call出snp")

    @time_count
    def export_ssr(self):
        gevent.sleep()
        task_id = self.task_id
        project_sn = self.project_sn
        params = dict(
            task_id=task_id,
            submit_location="ssr",
            task_type=2,
            rept_1=10,
            rept_2=6,
            rept_3=5,
            rept_4=5,
            rept_5=5,
            rept_6=5,
            ssr_distance=100,
        )
        api_ssr = self.api.api("denovo_rna_v2.ssr")
        self.logger.info("开始进行ssr的导表")
        ssr_work_dir = self.function_infos["export_ssr"]["ssr_work_dir"]
        ssr_statistic_path = ssr_work_dir + '/' + 'ssr_type.txt'
        ssr_detail_path = ssr_work_dir + '/' + 'tmp.txt'
        ssr_class_path = ssr_work_dir + '/' + 'ssr_repeats_class.txt'
        api_ssr.add_ssr(ssr_statistic_path, ssr_detail_path, ssr_class_path, name=None, params=params,
                        project_sn=project_sn, task_id=task_id)

    @time_count
    def export_tf(self):
        gevent.sleep()
        task_id = self.task_id
        project_sn = self.project_sn
        params = dict(
            task_id=task_id,
            task_type=2,
            search_pfam='True',
            p_length=50,
            Markov_length=3000,
            cpu=20,
            hmmcan1='noali',
            hmmcan2='acc',
            hmmcan3='notextw',
            E=self.workflow_option("tf_evalue")
        )
        # self.fasta_name = self.assemble.option("filter_fa").prop["path"].split("/")[-1]
        # bed = '{}.transdecoder.bed'.format(self.fasta_name)
        bed = "all_predicted.bed"
        cds_predict_work_dir = self.function_infos["export_tf"]["cds_predict_work_dir"]
        bedpath = cds_predict_work_dir + '/' + bed
        if self.workflow_option("tf_database").lower() == "animal":
            tf_unigene_path = os.path.join(cds_predict_work_dir, 'Predict', 'merge_only_unigene_animal')
            tf_transcript_path = os.path.join(cds_predict_work_dir, 'Predict', 'merge_only_transcript_animal')
            api_tf_api = self.api.api("denovo_rna_v2.tf_api")
            self.logger.info("开始进行动物tf的导表")
            api_tf_api.add_tf_unigene(tf_unigene_path, bedpath=bedpath, name=None, params=params,
                                      project_sn=project_sn, task_id=task_id)
            self.logger.info("动物的tf基因水平的导表完成")
            api_tf_api.add_tf_transcript(tf_transcript_path, bedpath=bedpath, name=None, params=params,
                                         project_sn=project_sn, task_id=task_id)
            self.logger.info("动物的tf基因，转录水平的导表完成")

        if self.workflow_option("tf_database").lower() == "plant":
            tf_unigene_path = os.path.join(cds_predict_work_dir, 'Predict', 'merge_only_unigene_plant')
            tf_transcript_path = os.path.join(cds_predict_work_dir, 'Predict', 'merge_only_transcript_plant')
            api_tf_api = self.api.api("denovo_rna_v2.tf_api")
            self.logger.info("开始进行植物tf的导表")
            api_tf_api.add_tf_unigene(tf_unigene_path, bedpath=bedpath, name=None, params=params,
                                      project_sn=project_sn, task_id=task_id)
            self.logger.info("植物的tf基因水平的导表完成")
            api_tf_api.add_tf_transcript(tf_transcript_path, bedpath=bedpath, name=None, params=params,
                                         project_sn=project_sn, task_id=task_id)
            self.logger.info("植物的tf基因，转录水平的导表完成")


    @time_count
    def export_cds(self):
        gevent.sleep()
        task_id = self.task_id
        project_sn = self.project_sn
        params = dict(
            task_id=task_id,
            task_type=2,
            search_pfam='True',
            p_length=50,
            Markov_length=3000,
            cpu=20,
            hmmcan1='noali',
            hmmcan2='acc',
            hmmcan3='notextw',
            E=self.workflow_option("tf_evalue")
        )
        api_cds = self.api.api("denovo_rna_v2.cdslen")
        self.logger.info("开始进行cds的导表")
        cds_predict_output_dir = self.function_infos["export_cds"]["cds_predict_output_dir"]
        cds_unigene_length = cds_predict_output_dir + '/' + 'cds_len_unigene.txt'
        cds_transcript_length = cds_predict_output_dir + '/' + 'cds_len_transcript.txt'
        all_predicted = cds_predict_output_dir + '/' + 'all_predicted.xls'
        api_cds.cds(all_predicted, project_sn=self.project_sn, task_id=self.task_id,
                    cds_unigene_length=cds_unigene_length, cds_transcript_length=cds_transcript_length, params=params)

    @time_count
    def export_expression(self):
        gevent.sleep()
        all_exp = self.api.api("denovo_rna_v2.all_exp")
        if self.option("group").is_set:
            group_dict = self.option('group').prop['group_dict']
            group_id = self.group_id
        else:
            group_dict = None
            group_id = None
        if self.option("control").is_set:
            control_id = self.control_id
        else:
            control_id = None
        quant_method = self.workflow_option('express_method')
        task_id = self.task_id
        project_sn = self.project_sn

        # add exp matrix
        ## 还需确认路径信息，因为后续交互需要用到express.output_dir
        exp_output = self.function_infos["export_denovo_annotation"]["exp_output"]
        # if self.workflow_option("express_method") == "RSEM":
        #     exp_output = self.align.output_dir
        # else:
        #     exp_output = self.express.output_dir
        if self.workflow_option("express_method") == "RSEM" and self.workflow_option("exp_way") == "fpkm":
            params = dict(
                task_id=task_id,
                submit_location="exp_detail",
                task_type=2,
                method=quant_method,
                exp_type='FPKM',
            )
            libtype = self.function_infos["export_expression"]["libtype"]
            exp_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
            if self.workflow_option("sample_num") == "multiple":
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=True,
                                               exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=False,
                                               exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)

            exp_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
            if self.workflow_option("sample_num") == "multiple":
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=True,
                                              exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=False,
                                              exp_type='FPKM', project_sn=project_sn, task_id=task_id, params=params)

        else:
            params = dict(
                task_id=task_id,
                submit_location="exp_detail",
                task_type=2,
                method=quant_method,
                exp_type='TPM',
            )
            libtype = self.function_infos["export_expression"]["libtype"]

            exp_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            if self.workflow_option("sample_num") == "multiple":
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=True,
                                               exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                trans_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='T',
                                               lib_type=libtype,
                                               group_dict=group_dict, group_id=group_id, add_distribution=False,
                                               exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)

            exp_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
            if self.workflow_option("sample_num") == "multiple":
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=True,
                                              exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)
            else:
                gene_exp_id = all_exp.add_exp(exp_matrix, quant_method=quant_method, exp_level='G',
                                              lib_type=libtype,
                                              group_dict=group_dict, group_id=group_id, add_distribution=False,
                                              exp_type='TPM', project_sn=project_sn, task_id=task_id, params=params)


        if self.workflow_option("sample_num") == "multiple":
            # add gene corr
            corr_output = self.function_infos["export_expression"]["corr_output"]
            params = dict(
                task_id=task_id,
                submit_location='expcorr',
                task_type=2,
                exp_id=str(gene_exp_id),
                group_id=str(group_id),
                exp_level="G",
                group_dict=group_dict,
                scm="complete",
                scd="euclidean",
                # quant_method=quant_method,
                corr_method="pearson",
                draw_in_groups="no"
            )
            all_exp.add_exp_corr2(corr_output, exp_level='G', quant_method=quant_method, params=params,
                                  project_sn=project_sn, task_id=task_id)

            # add gene venn
            if len(self.option('group').prop['group_dict']) > 1:
                graph_table = self.function_infos["export_expression"]["graph_table"]
                group_dict = self.option('group').prop['group_dict']
                if len(group_dict) > 6:
                    group_dict = OrderedDict(group_dict.items()[:6])
                params = json.dumps(dict(
                    task_id=self.task_id,
                    submit_location='expvenn',
                    task_type=2,
                    exp_id=str(gene_exp_id),
                    group_id=str(group_id),
                    exp_level='G',
                    group_dict=group_dict,
                    threshold=1,
                ), sort_keys=True, separators=(',', ':'))
                import datetime
                time_now = datetime.datetime.now()
                name = 'ExpVenn_G_{}_{}_{}'.format(
                    self.workflow_option('express_method'), self.workflow_option('exp_way').upper(), time_now.strftime('%Y%m%d_%H%M%S'))
                main_info = dict(
                    project_sn=self.project_sn,
                    task_id=self.task_id,
                    version='v2',
                    name=name,
                    created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                    exp_id=str(gene_exp_id),
                    desc='Expression venn analysis main table',
                    params=params,
                    status='start'
                )
                main_id = all_exp.create_db_table('sg_exp_venn', [main_info])
                all_exp.add_exp_venn(graph_table, main_id=main_id)

            # add gene pca
            if self.option("group").prop["sample_number"] > 2:
                pca_output = self.function_infos["export_expression"]["pca_output"]
                params = dict(
                    task_id=task_id,
                    submit_location="exppca",
                    task_type=2,
                    exp_id=str(gene_exp_id),
                    group_id=str(group_id),
                    exp_level="G",
                    group_dict=self.option('group').prop['group_dict'],
                    draw_in_groups="no"
                    # quant_method=quant_method,
                )
                main_id = all_exp.add_exp_pca2(pca_output, quant_method=quant_method, exp_id=gene_exp_id, exp_level="G",
                                               params=params, project_sn=project_sn, task_id=task_id)
                if 'ellipse' in self.function_infos["export_expression"]:
                    all_exp.insert_ellipse_table(self.function_infos["export_expression"]["ellipse"], main_id)

            # add transcript diff
            if self.workflow_option("level").lower() == "transcript":
                diff_output = self.function_infos["export_expression"]["trans_diff_output"]
                uniform_output = self.function_infos["export_expression"]["trans_uniform_output"]
                exp_id, exp_level = trans_exp_id, 'T'
                diff_method = self.workflow_option('diff_method')
                stat_type = self.workflow_option('pvalue_padjust')
                params = dict(
                    task_id=task_id,
                    submit_location="diff_detail",
                    task_type=2,
                    exp_id=str(exp_id),
                    group_id=str(group_id),
                    control_id=str(control_id),
                    exp_level=exp_level,
                    group_dict=self.option('group').prop['group_dict'],
                    fc=str(float(self.workflow_option('fc'))),
                    stat_type=stat_type,
                    stat_cutoff=self.workflow_option('diff_fdr_ci'),
                    # correct_method=self.option('padjust_way'),
                    # quant_method=quant_method,
                    # filter_method="no",
                    is_batch="False",
                    diff_method=diff_method,
                )
                if self.workflow_option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                    params.update({"correct_method": self.workflow_option('padjust_way')})
                all_exp.add_diffexp_all(uniform_output, diff_output, exp_id=exp_id,
                                        group_dict=group_dict, group_id=group_id,
                                        exp_level=exp_level, quant_method=quant_method,
                                        diff_method=diff_method,
                                        project_sn=project_sn, task_id=task_id, params=params,
                                        pvalue_padjust=stat_type
                                        )

            # add gene diff
            diff_output = self.function_infos["export_expression"]["gene_diff_output"]
            uniform_output = self.function_infos["export_expression"]["gene_uniform_output"]
            exp_id, exp_level = gene_exp_id, 'G'
            diff_method = self.workflow_option('diff_method')
            stat_type = self.workflow_option('pvalue_padjust')
            params = dict(
                task_id=task_id,
                submit_location="diff_detail",
                task_type=2,
                exp_id=str(exp_id),
                group_id=str(group_id),
                control_id=str(control_id),
                exp_level=exp_level,
                group_dict=self.option('group').prop['group_dict'],
                fc=str(float(self.workflow_option('fc'))),
                # correct_method=self.option('padjust_way'),
                stat_type=stat_type,
                stat_cutoff=self.workflow_option('diff_fdr_ci'),
                # quant_method=quant_method,
                # filter_method="no",
                is_batch="False",
                diff_method=diff_method,
            )
            if self.workflow_option("diff_method").lower() in ["degseq", "edger", "deseq2", 'limma']:
                params.update({"correct_method": self.workflow_option('padjust_way')})

            all_exp.add_diffexp_all(uniform_output, diff_output, exp_id=exp_id,
                                    group_dict=group_dict, group_id=group_id,
                                    exp_level=exp_level, quant_method=quant_method,
                                    diff_method=diff_method,
                                    project_sn=project_sn, task_id=task_id, params=params,
                                    pvalue_padjust=stat_type
                                    )


    def run(self):
        super(SetDbTool, self).run()
        self.build_seq_database()
        self.run_api()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.large.set_db",
            "instant": False,
            "options": dict(
                sheet_data_json = "/mnt/ilustre/users/sanger-dev/workspace/20210506/Refrna_ngj8_he50cvn7kga74plmlp7usb/workflow_sheet_data.json",
                option_data_json = "/mnt/ilustre/users/sanger-dev/workspace/20210506/Refrna_ngj8_he50cvn7kga74plmlp7usb/option_data.json",
                analysis_content = "[\"mapping\",\"annotation\",\"quantification\",\"other\"]"
                # merge_file = "/mnt/ilustre/users/isanger/workspace/20210318/Refrna_st73_2iu0i6el0bq566603jk30t/Quant/quant_large.json",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()