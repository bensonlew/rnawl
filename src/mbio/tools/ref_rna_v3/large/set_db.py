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
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'function_json', "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'control_file', 'type': 'infile', 'format': 'sample.control_table'},
            {'name': 'report_img', 'type': 'bool', 'default': False},

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
        self.analysis_content = self.option("analysis_content")
        self.workflow_options  = self.get_option_data()
        self.function_infos = self.get_function_json()
        self.group_id = ""
        self.control_id = ""
        self.task_id = self.option("task_id")
        self.project_sn =  self.sheet.project_sn
        self.exp_ids = dict()
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.diff_id = ''
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


    @workfuncdeco
    def export_task_info(self):
        api = self.api.api('task_info.ref_rna_v2')
        api.add_task_info()

    @workfuncdeco
    def export_genome_info(self):
        api = self.api.api('ref_rna_v2.genome_info')
        genome_stat = self.function_infos["export_genome_info"]["file_path"]
        species_name = self.function_infos ["export_genome_info"]["species_name"]
        species = self.function_infos["export_genome_info"]["species"]
        ref_anno_version = self.function_infos["export_genome_info"]["ref_anno_version"]
        hyperlink = self.function_infos["export_genome_info"]["hyperlink"]
        api.add_genome_info(file_path=genome_stat, species_name=species_name, species=species,
                            ref_anno_version=ref_anno_version, hyperlink=hyperlink)

    @workfuncdeco
    def export_ref_rna_qc_before(self):
        api = self.api.api('ref_rna_v2.ref_rna_qc')
        qc_stat = self.function_infos["export_ref_rna_qc_before"]["qc_stat"]
        fq_type = self.function_infos["export_ref_rna_qc_before"]["fq_type"]
        quality_stat = self.function_infos["export_ref_rna_qc_before"]["quality_stat"]
        api.add_samples_info(qc_stat = qc_stat, fq_type = fq_type,
                             about_qc = 'before',group=self.option('group_table').path)
        api.add_gragh_info(quality_stat= quality_stat,about_qc='before')

    @workfuncdeco
    def export_productive_name(self):
        api = self.api.api('ref_rna_v2.ref_rna_qc')
        if self.option('productive_table').is_set:
            api.add_productive_name(samples=self.option('group_table').prop['sample'],
                                    productive_table=self.option('productive_table').path)
        else:
            pass

    @workfuncdeco
    def export_ref_rna_qc_after(self):
        api = self.api.api('ref_rna_v2.ref_rna_qc')
        qc_stat = self.function_infos["export_ref_rna_qc_after"]["qc_stat"]
        fq_type = self.function_infos["export_ref_rna_qc_after"]["fq_type"]
        quality_stat = self.function_infos["export_ref_rna_qc_after"]["quality_stat"]
        api.add_samples_info(qc_stat=qc_stat, fq_type=fq_type,
                             about_qc='after', group=self.option('group_table').path)
        api.add_gragh_info(quality_stat=quality_stat, about_qc='after')
        if self.workflow_option('sample_num') == 'multiple':
            if self.option("group_table").is_set:
                self.group_id, specimen_names, category_names = api.add_specimen_group(self.option('group_table').path)
            if self.option('control_file').is_set:
                self.control_id, compare_detail = api.add_control_group(self.option('control_file').path, self.group_id)
        else:
            sp_set = set()
            if self.workflow_options('datatype') == 'rawdata':
                fastq_dir = self.option('fastq_dir').path
            else:
                fastq_dir = self.option('qc_dir').path
            for line in open(os.path.join(fastq_dir, 'list.txt')):
                sp_set.add(line.strip().split('\t')[1])
            if len(sp_set) != 1:
                self.set_error('invalid sample number')
            sample = list(sp_set)[0]
            group_table = os.path.join(self.work_dir, 'group.txt')
            open(group_table, 'w').write('#sample\tgroup\n{}\t{}'.format(sample, sample))
            self.group_id, specimen_names, category_names = api.add_specimen_group(group_table)
        intermediate_dir = self.sheet.output.replace('workflow_results', 'intermediate_results/')
        api.add_bam_path(intermediate_dir)

    @workfuncdeco
    def export_ref_rna_qc_alignment(self):
        api = self.api.api('ref_rna_v2.ref_rna_qc')
        stat_file = self.function_infos["export_ref_rna_qc_alignment"]["alignment"]
        if self.workflow_option('align_method').lower() == 'hisat':
            api.add_hisat_mapping_stat(stat_file,self.option('group_table').path)
        elif self.workflow_option('align_method').lower() == 'tophat':
            api.add_tophat_mapping_stat(stat_file,self.option('group_table').path)
        elif self.workflow_option('align_method').lower() == 'star':
            api.add_star_mapping_stat(stat_file,self.option('group_table').path)


    @workfuncdeco
    def export_ref_rna_qc_assessment(self):
        api = self.api.api('ref_rna_v2.ref_rna_qc')
        if 'saturation' in self.workflow_option('map_assess_method'):
            api.add_rpkm_table(file_path=self.function_infos["export_ref_rna_qc_assessment"]["saturation"], group=self.option('group_table').path)
        if 'coverage' in self.workflow_option('map_assess_method'):
            api.add_coverage_table(coverage=self.function_infos["export_ref_rna_qc_assessment"]["coverage"], group=self.option('group_table').path)
        if 'distribution' in self.workflow_option('map_assess_method'):
            api.add_distribution_table(distribution=self.function_infos["export_ref_rna_qc_assessment"]["distribution"], group=self.option('group_table').path)
        if 'chr_stat' in self.workflow_option('map_assess_method'):
            api.add_chrom_distribution_table(distribution=self.function_infos["export_ref_rna_qc_assessment"]["chr_stat"], group=self.option('group_table').path)

    @workfuncdeco
    def export_ref_assembly(self):
        api = self.api.api('ref_rna_v2.ref_assembly')
        params = json.dumps({'task_id': self.task_id, 'submit_location': 'transcripts', 'task_type': 2}, sort_keys=True)
        all_gtf_path = self.function_infos["export_ref_assembly"]["all_gtf_path"]
        merged_path =  self.function_infos["export_ref_assembly"]["merged_path"]
        statistics_path = self.function_infos["export_ref_assembly"]["statistics_path"]
        api.add_assembly_result(params=params, all_gtf_path=all_gtf_path, merged_path=merged_path,
                                statistics_path=statistics_path)

    @workfuncdeco
    def export_annotation(self):
        api_annotation = self.api.api("ref_rna_v2.ref_annotation")
        annot_dir = self.function_infos["export_annotation"]["annot_dir"]
        if not self.option("is_assemble"):
            api_annotation.has_new = False
            trans2gene = None
            trans2gene_ref = self.function_infos["export_annotation"]["trans2gene_ref"]
        else:
            trans2gene = self.function_infos["export_annotation"]["trans2gene"]
            trans2gene_ref = self.function_infos["export_annotation"]["trans2gene_ref"]
        api_annotation.species_name = self.function_infos["export_annotation"]["ref_genome"]
        api_annotation.has_new = self.option('is_assemble')
        params_dict = {
            "nr_evalue": str(self.workflow_option("nr_evalue")),
            "nr_similarity": self.workflow_option("nr_similarity"),
            "nr_identity": self.workflow_option("nr_identity"),
            "swissprot_evalue": str(self.workflow_option("swissprot_evalue")),
            "swissprot_similarity": self.workflow_option("swissprot_similarity"),
            "swissprot_identity": self.workflow_option("swissprot_identity"),
            "cog_evalue": str(self.workflow_option("cog_evalue")),
            "cog_similarity": self.workflow_option("cog_similarity"),
            "cog_identity": self.workflow_option("cog_identity"),
            "kegg_evalue": str(self.workflow_option("kegg_evalue")),
            "kegg_similarity": self.workflow_option("kegg_similarity"),
            "kegg_identity": self.workflow_option("kegg_identity"),
            "pfam_evalue": str(self.workflow_option("pfam_evalue")),
        }
        if "quantification" in self.analysis_content:
            gene_exp = self.function_infos["export_annotation"]["gene_exp"]
            trans_exp = self.function_infos["export_annotation"]["trans_exp"]
        else:
            gene_exp = None
            trans_exp = None
        api_annotation.run(
            annot_dir,
            trans2gene,
            trans2gene_ref,
            params_dict=params_dict,
            taxonomy=self.workflow_option("kegg_database"),
            exp_level=self.workflow_option("level").lower(),
            version="v3",
            gene_exp=gene_exp,
            trans_exp=trans_exp,
        )

    @workfuncdeco
    def export_all_exp_matrix(self):
        api = self.api.api('ref_rna_v2.all_exp')
        exp_type = self.function_infos["export_all_exp_matrix"]["exp_type"]
        exp_matrix = self.function_infos["export_all_exp_matrix"]["exp_matrix"]
        quant_method = self.function_infos["export_all_exp_matrix"]["quant_method"]
        lib_type = self.function_infos["export_all_exp_matrix"]["lib_type"]
        if self.workflow_option('sample_num') == 'multiple':
            group_dict = self.option('group_table').prop['group_dict']
        else:
            group_table = os.path.join(self.work_dir, 'group.txt')
            for line in open(group_table):
                if line.strip() and line[0] != '#':
                    sample = line.strip().split('\t')[1]
            group_dict = {sample: [sample]}
        group_id = self.group_id
        project_sn = self.project_sn
        task_id = self.task_id
        params = json.dumps({
            'task_id': task_id,
            'submit_location': 'exp_detail',
            'task_type': 2,
            'method': quant_method,
            'exp_type': exp_type
        }, sort_keys=True, separators=(',', ':'))
        self.exp_ids = dict()
        if self.workflow_option('level').lower() == 'transcript':
            self.exp_ids['T'] = api.add_exp(exp_matrix=exp_matrix['T'], quant_method=quant_method, exp_level='T',
                                            lib_type=lib_type, group_dict=group_dict, group_id=group_id,
                                            exp_type=exp_type, add_distribution=False, project_sn=project_sn,
                                            task_id=task_id, params=params)
        self.exp_ids['G'] = api.add_exp(exp_matrix=exp_matrix['G'], quant_method=quant_method, exp_level='G',
                                        lib_type=lib_type, group_dict=group_dict, group_id=group_id,
                                        exp_type=exp_type, add_distribution=False, project_sn=project_sn,
                                        task_id=task_id, params=params)

    @workfuncdeco
    def export_all_exp_distribution(self):
        api = self.api.api('ref_rna_v2.all_exp')
        task_id = self.task_id
        exp_way = self.function_infos["export_all_exp_distribution"]["exp_way"]
        group_dict = self.option('group_table').prop['group_dict']
        exp_ids = self.exp_ids
        group_id = self.group_id

        def params(exp_level):
            return json.dumps({
                'task_id': task_id,
                'submit_location': 'expgraph',
                'task_type': 2,
                'exp_id': str(exp_ids[exp_level]),
                'group_dict': group_dict,
                'group_id': str(group_id),
                'exp_level': exp_level,
                'type': 'ref'
            }, sort_keys=True, separators=(',', ':'))

        quant_method = self.workflow_option('express_method')
        project_sn = self.project_sn
        if self.workflow_option('level').lower() == 'transcript':
            api.add_distribution(exp_matrix=self.function_infos["export_all_exp_distribution"]["trans_exp_matrix"], group_dict=group_dict, params=params('T'),
                                 exp_level='T', quant_method=quant_method, project_sn=project_sn, task_id=task_id)
        api.add_distribution(exp_matrix = self.function_infos["export_all_exp_distribution"]["gene_exp_matrix"], group_dict=group_dict, params=params('G'), exp_level='G',
                             quant_method=quant_method, project_sn=project_sn, task_id=task_id)

    @workfuncdeco
    def export_add_exp_venn(self):
        api = self.api.api('ref_rna_v2.all_exp')
        graph_table = self.function_infos["export_add_exp_venn"]["graph_table"]
        group_dict = self.option('group_table').prop['group_dict']
        if len(group_dict) > 6:
            group_dict = OrderedDict(group_dict.items()[:6])
        params = json.dumps(dict(
            task_id=self.task_id,
            submit_location='expvenn',
            task_type=2,
            exp_id=str(self.exp_ids['G']),
            group_id=str(self.group_id),
            exp_level='G',
            group_dict=group_dict,
            threshold='1',
            type='ref',
        ), sort_keys=True, separators=(',', ':'))
        import datetime
        time_now = datetime.datetime.now()
        name = 'ExpVenn_G_{}_{}_{}'.format(
            self.workflow_option('express_method'), self.workflow_option('exp_way').upper(), time_now.strftime('%Y%m%d_%H%M%S'))
        main_info = dict(
            project_sn=self.project_sn,
            task_id=self.task_id,
            version='v3',
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            exp_id=str(self.exp_ids['G']),
            desc='Expression venn analysis main table',
            params=params,
            status='start'
        )
        main_id = api.create_db_table('sg_exp_venn', [main_info])
        api.add_exp_venn(graph_table, main_id=main_id)

    @workfuncdeco
    def export_all_exp_pca(self):
        api = self.api.api('ref_rna_v2.all_exp')
        pca_output_dir = self.function_infos["export_all_exp_pca"]["pca_output_dir"]
        quant_method = self.workflow_option('express_method')
        task_id = self.task_id
        exp_ids = self.exp_ids
        exp_way = self.workflow_option('exp_way')
        group_dict = self.option('group_table').prop['group_dict']
        group_id = self.group_id

        def params(exp_level):
            return json.dumps({
                'task_id': task_id,
                'submit_location': 'exppca',
                'task_type': 2,
                'exp_id': str(exp_ids[exp_level]),
                'group_dict': group_dict,
                'group_id': str(group_id),
                'exp_level': exp_level,
                'type': 'ref',
                'Draw_in_groups': "no"
            }, sort_keys=True, separators=(',', ':'))

        project_sn = self.project_sn
        main_id = api.add_exp_pca2(pca_output_dir, quant_method=quant_method, exp_id=exp_ids['G'], exp_level='G',
                                   params=params('G'), project_sn=project_sn, task_id=task_id)
        if 'ellipse' in self.function_infos["export_all_exp_pca"]:

            api.insert_ellipse_table(self.function_infos["export_all_exp_pca"]["ellipse"], main_id)

    @workfuncdeco
    def export_all_exp_corr(self):
        api = self.api.api('ref_rna_v2.all_exp')
        corr_work_dir = self.function_infos["export_all_exp_corr"]["corr_work_dir"]
        quant_method = self.workflow_option('express_method')
        task_id = self.task_id
        exp_ids = self.exp_ids
        group_dict = self.option('group_table').prop['group_dict']
        group_id = self.group_id

        def params(exp_level):
            return json.dumps({
                'task_id': task_id,
                'submit_location': 'expcorr',
                'task_type': 2,
                'exp_id': str(exp_ids[exp_level]),
                'group_dict': group_dict,
                'group_id': str(group_id),
                'exp_level': exp_level,
                'scm': 'complete',
                'scd': 'euclidean',
                'corr_method': 'pearson',
                'type': 'ref'
            }, sort_keys=True, separators=(',', ':'))

        project_sn = self.project_sn
        api.add_exp_corr2(corr_work_dir, exp_level='G', quant_method=quant_method, params=params('G'),
                          project_sn=project_sn, task_id=task_id)

    @workfuncdeco
    def export_gene_detail(self):
        api = self.api.api('ref_rna_v3.gene_detail')
        refrna_seqdb = self.function_infos["export_gene_detail"]["refrna_seqdb"]
        t2g_file = self.function_infos["export_gene_detail"]["t2g_file"]
        txpt_fa = self.function_infos["export_gene_detail"]["txpt_fa"]
        new_cds = self.function_infos["export_gene_detail"]["new_cds"]
        new_pep = self.function_infos["export_gene_detail"]["new_pep"]
        txpt_bed  = self.function_infos["export_gene_detail"]["txpt_bed"]
        gene_bed = self.function_infos["export_gene_detail"]["gene_bed"]
        gene_fa = self.function_infos["export_gene_detail"]["gene_fa"]
        biomart_file = self.function_infos["export_gene_detail"]["biomart_file"]
        biomart_type = self.function_infos["export_gene_detail"]["biomart_type"]
        species_urls = self.function_infos["export_gene_detail"]["species_urls"]

        api.add_gene_detail(refrna_seqdb, t2g_file, txpt_bed, txpt_fa, gene_bed, gene_fa,
                            biomart_file, biomart_type, species_urls, new_cds, new_pep)

        gene_stat = self.function_infos["export_gene_detail"]["gene_stat"]
        trans_stat = self.function_infos["export_gene_detail"]["trans_stat"]
        api_seq = self.api.api('ref_rna_v2.seq_detail')
        api_seq.add_seq_stat(gene_stat, trans_stat)

    @workfuncdeco
    def export_all_exp_diff(self):
        api = self.api.api('ref_rna_v2.all_exp')
        diff_output = self.function_infos["export_all_exp_diff"]["diff_output"]
        exp_ids = self.exp_ids
        group_dict = self.option('group_table').prop['group_dict']
        group_id = self.group_id
        quant_method = self.workflow_option('express_method')
        diff_method = self.workflow_option('diff_method')
        project_sn = self.project_sn
        task_id = self.task_id
        control_id = self.control_id
        fc = str(float(self.workflow_option('fc')))
        # if '.' in fc:
        #     if fc.split('.')[1] == '0':
        #         fc = str(int(float(fc)))
        correct_method = self.workflow_option('padjust_way')
        stat_type = self.workflow_option('pvalue_padjust')
        stat_cutoff = str(self.workflow_option('diff_fdr_ci'))
        tpm_filter_threshold = str(float(self.workflow_option('filter_tpm')))

        # if '.' in tpm_filter_threshold:
        #     if tpm_filter_threshold.split('.')[1] == '0':
        #         tpm_filter_threshold = str(int(float(tpm_filter_threshold)))

        def params(exp_level):
            if diff_method.lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                params_dict = {
                    'task_id': task_id,
                    'submit_location': 'diff_detail',
                    'task_type': 2,
                    'exp_id': str(exp_ids[exp_level]),
                    'group_id': str(group_id),
                    'control_id': str(control_id),
                    'exp_level': exp_level,
                    'group_dict': group_dict,
                    'fc': fc,
                    'tpm_filter_threshold': tpm_filter_threshold,
                    'stat_type': stat_type,
                    'stat_cutoff': stat_cutoff,
                    'diff_method': diff_method,
                    'type': 'ref',
                    'is_batch': 'False',
                }
                if stat_type == 'padjust':
                    params_dict.update({'correct_method': correct_method})
            else:
                params_dict = {
                    'task_id': task_id,
                    'submit_location': 'diff_detail',
                    'task_type': 2,
                    'exp_id': str(exp_ids[exp_level]),
                    'group_id': str(group_id),
                    'control_id': str(control_id),
                    'exp_level': exp_level,
                    'group_dict': group_dict,
                    'fc': fc,
                    'tpm_filter_threshold': tpm_filter_threshold,
                    'stat_cutoff': stat_cutoff,
                    'diff_method': diff_method,
                    'type': 'ref',
                    'is_batch': 'False',
                    'prob': float(self.workflow_option('diff_fdr_ci'))
                }
            return json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

        s3_output_dir = os.path.join(self.sheet.output, "07DiffExpress_G")
        if diff_method.lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            self.diff_id = api.add_diffexp(diff_output, exp_id=exp_ids['G'], group_dict=group_dict, group_id=group_id, exp_level='G',
                            quant_method=quant_method, diff_method=diff_method, project_sn=project_sn, task_id=task_id,
                            params=params('G'), pvalue_padjust=stat_type,s3_output_dir = s3_output_dir)
        else:
            self.diff_id = api.add_diffexp_noiseq(diff_output, exp_id=exp_ids['G'], group_dict=group_dict, group_id=group_id,
                                   exp_level='G',
                                   quant_method=quant_method, diff_method=diff_method, project_sn=project_sn,
                                   task_id=task_id,
                                   params=params('G'),s3_output_dir = s3_output_dir)

    @workfuncdeco
    def export_rmats(self):
        api = self.api.api('ref_rna_v3.rmats')
        rmats_output = self.function_infos["export_rmats"]["rmats_output"]
        for p in glob.glob(os.path.join(rmats_output, '*')):
            if os.path.isdir(p):
                outpath = p
                ctrl, test = os.path.basename(outpath).split('_vs_')
                group_dict = {ctrl: self.option('group_table').prop['group_dict'][ctrl],
                              test: self.option('group_table').prop['group_dict'][test]}
                compare_plan = '{}|{}'.format(ctrl, test)
                params = json.dumps({
                    'task_id': self.task_id,
                    'submit_location': 'splicingrmats',
                    'task_type': 2,
                    'group_id': str(self.group_id),
                    'group_dict': group_dict,
                    'control_id': str(self.control_id),
                    'compare_plan': compare_plan
                }, sort_keys=True, separators=(',', ':'))
                api.add_sg_splicing_rmats(params, outpath)

    @workfuncdeco
    def export_rmats_count(self):
        api = self.api.api('ref_rna_v3.rmats_count')
        rmats_output = self.function_infos["export_rmats"]["rmats_output"]
        api.add_rmats_count(rmats_output, group=self.option('group_table').path)

    @workfuncdeco
    def export_snp(self):
        api = self.api.api('ref_rna_v3.ref_snp')
        task_id = self.task_id
        project_sn = self.project_sn
        new_output = self.function_infos["export_snp"]["new_output"]
        if os.path.exists(new_output):
            shutil.rmtree(new_output)
        os.mkdir(new_output)
        if self.workflow_option('snp_method').lower() == 'gatk':
            snp_anno = self.function_infos["export_snp"]["snp_anno"]
            if os.path.exists(snp_anno + "/snp_annotation_statistics.xls"):
                params = dict(
                    task_id=task_id,
                    submit_location='snp',
                    task_type=2,
                    method_type='gatk'
                )
                api.add_snp_main(snp_anno=snp_anno, group=self.option('group_table').path, params=params,
                                 task_id=task_id, method_type='gatk',
                                 project_sn=project_sn, new_output=new_output)
        if self.workflow_option('snp_method').lower() == 'samtools':
            snp_anno =  self.function_infos["export_snp"]["snp_anno"]
            if os.path.exists(snp_anno + "/snp_annotation_statistics.xls"):
                params = dict(
                    task_id=task_id,
                    submit_location='snp',
                    task_type=2,
                    method_type='samtools'
                )
                api.add_snp_main(snp_anno=snp_anno, group=self.option('group_table').path, params=params,
                                 task_id=task_id, method_type='samtools',
                                 project_sn=project_sn, new_output=new_output)
        if self.workflow_option('snp_method').lower() == 'sentieon':
            snp_anno =  self.function_infos["export_snp"]["snp_anno"]
            if os.path.exists(snp_anno + "/snp_annotation_statistics.xls"):
                params = dict(
                    task_id=task_id,
                    submit_location='snp',
                    task_type=2,
                    method_type='sentieon'
                )
                api.add_snp_main(snp_anno=snp_anno, group=self.option('group_table').path, params=params,
                                 task_id=task_id, method_type='sentieon',
                                 project_sn=project_sn, new_output=new_output)
        # if os.path.exists(os.path.join(self.work_dir, 'SnpTmp/snp_anno.xls')):
        #     os.remove(os.path.join(self.work_dir, 'SnpTmp/snp_anno.xls'))
        # os.link(os.path.join(snp_anno, 'data_anno_pre.xls'), os.path.join(self.work_dir, 'SnpTmp/snp_anno.xls'))

    @workfuncdeco
    def export_diff_geneset_analysis(self):
        export_temporary = self.function_infos["export_diff_geneset_analysis"]["export_temporary"]
        diff_geneset_pipline_result =  self.function_infos["export_diff_geneset_analysis"]["diff_geneset_pipline_result"]
        diff_id = self.diff_id
        task_id = self.task_id
        analysis_names = self.function_infos["export_diff_geneset_analysis"]["analysis_names"]
        kegg_level_path = self.function_infos["export_diff_geneset_analysis"]["kegg_level_path"]
        gene_detail = self.function_infos["export_diff_geneset_analysis"]["gene_detail"]
        api = self.api.api('ref_rna_v2.diff_geneset_work_pipline')
        api.add_diff_genest_pipline_table(diff_geneset_pipline_result, diff_id=diff_id, task_id=task_id,
                                          analysis_names=analysis_names,
                                          kegg_level_path=kegg_level_path, inter_path=export_temporary,
                                          exp_id=self.exp_ids['G'])

    @workfuncdeco
    def export_project_overview(self):
        api = self.api.api('ref_rna_v3.project_overview')
        group_dict = self.option('group_table').prop['group_dict']
        sample_num = 0
        for g in group_dict:
            sample_num += len(group_dict[g])
        api.add_project_overview(task_id=self.task_id, group_dict=group_dict, sample_num=sample_num,
                                 exp_level=self.workflow_option("level").lower())

    @workfuncdeco
    def export_report_img(self):
        report_config = self.function_infos['export_report_img']['report_config']
        api = self.api.api('ref_rna_v2.report_model')
        # s3 = self._sheet.output.split(":")[0]
        report_img_s3 = self.function_infos["export_report_img"]["report_img_s3"]
        api.add_report_image(self.task_id, report_config, report_img_s3)

    def run(self):
        super(SetDbTool, self).run()
        if "mapping" in self.analysis_content:
            self.export_task_info()
            self.export_genome_info()
            if self.workflow_option('datatype') == 'rawdata':
                self.export_ref_rna_qc_before()
            self.export_ref_rna_qc_after()
            self.export_productive_name()
            self.export_ref_rna_qc_alignment()
            self.export_ref_rna_qc_assessment()
        if "annotation" in self.analysis_content:
            if self.option('is_assemble'):
                self.export_ref_assembly()
            self.export_annotation()
        if "quantification" in self.analysis_content:
            self.export_all_exp_matrix()
            if self.workflow_option('sample_num') == 'multiple':
                self.export_all_exp_distribution()
                if len(self.option('group_table').prop['group_dict']) > 1:
                    self.export_add_exp_venn()
                if self.option('group_table').prop['sample_number'] > 2:
                    self.export_all_exp_pca()
                    self.export_all_exp_corr()
                self.export_gene_detail()
        if "other" in self.analysis_content:
            if self.workflow_option('sample_num') == 'multiple':
                self.export_all_exp_diff()
                self.export_diff_geneset_analysis()
                if self.workflow_option('is_as') == 'True':
                    self.export_rmats()
                    self.export_rmats_count()
            if self.workflow_option('is_snp') == 'True':
                self.export_snp()
        self.export_project_overview()
        if self.option("report_img"):
            self.export_report_img()
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