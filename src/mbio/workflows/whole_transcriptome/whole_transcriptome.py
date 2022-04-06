# -*- coding:utf-8 -*-
# __author__ = 'shicaiping,liubinxu,qinjincheng'

import glob
import json
import os
import shutil
import unittest
import pandas as pd
from Bio import SeqIO
from biocluster.config import Config
from biocluster.workflow import Workflow
from mbio.packages.whole_transcriptome.utils import check_map_dict
from mbio.packages.ref_rna_v2.functions import tryforgood
from mbio.packages.project_demo.delete_demo import DeleteDemoMongo
import time
from collections import OrderedDict


class WholeTranscriptomeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WholeTranscriptomeWorkflow, self).__init__(wsheet_object)
        options = [
            ## 分析对象
            # 分析文库选择 ['long', 'small', 'circle']
            {'name': 'lib_select', 'type': 'string', 'default': 'longRNA-seq,smallRNA-seq'},
            # longRNA-Seq任务
            {'name': 'long_task_id', 'type': 'string', 'default': None},
            # smallRNA-Seq任务
            {'name': 'small_task_id', 'type': 'string', 'default': None},
            # circRNA-Seq任务
            {'name': 'circ_task_id', 'type': 'string', 'default': None},

            ## 关联分析 mRNA
            # SNP分析
            {'name': 'is_snp', 'type': 'bool', 'default': True},
            # SNP分析方法 ['sentieon', 'gatk', 'samtools']
            {'name': 'snp_method', 'type': 'string', 'default': 'sentieon'},
            # 可变剪切分析
            {'name': 'is_as', 'type': 'bool', 'default': True},

            ## 关联分析 lncRNA
            # 相关性算法 ['pearson', 'spearman', 'kendall']
            {'name': 'lt_cor_way', 'type': 'string', 'default': 'spearman'},
            # 相关性系数
            {'name': 'lt_cor_cut', 'type': 'string', 'default': '0.9'},
            # 显著性水平 ['padjust', 'pvalue']
            {'name': 'lt_sig_way', 'type': 'string', 'default': 'padjust'},
            # 显著性水平
            {'name': 'lt_sig_cut', 'type': 'string', 'default': '0.05'},
            # 多重检验校正方法 ['BH', 'Bonferroni', 'Holm', 'BY']
            {'name': 'lt_adj_way', 'type': 'string', 'default': 'BH'},
            # cis作用距离上游
            {'name': 'up_dis', 'type': 'string', 'default': '10'},
            # cis作用距离下游
            {'name': 'down_dis', 'type': 'string', 'default': '10'},
            # mirbase一级分类 ['Animal', 'Plant']
            {'name': 'mirbase_category', 'type': 'string', 'default': None},
            # mirbase具体物种
            {'name': 'mirbase_specie', 'type': 'string', 'default': None},

            ## 关联分析 smallRNA
            # 靶mRNA软件 ['miranda', 'rnahybird', 'targetscan', 'psrobot', 'targetfinder']
            {'name': 'm_method', 'type': 'string'},
            # miRanda score
            {'name': 'm_miranda_score', 'type': 'string', 'default': '160'},
            # miRanda energy
            {'name': 'm_miranda_energy', 'type': 'string', 'default': '-20'},
            # miRanda strict
            {'name': 'm_miranda_strict', 'type': 'string', 'default': 'on'},
            # RNAhybird num
            {'name': 'm_rnahybird_num', 'type': 'string', 'default': '100'},
            # RNAhybird energy
            {'name': 'm_rnahybird_energy', 'type': 'string', 'default': '-20'},
            # RNAhybird p-value
            {'name': 'm_rnahybird_pvalue', 'type': 'string', 'default': '0.01'},
            # psRobot score
            {'name': 'm_ps_robot_score', 'type': 'string', 'default': '2.5'},
            # Targetfinder score
            {'name': 'm_targetfinder_score', 'type': 'string', 'default': '4'},
            # 靶mRNA候选标准
            {'name': 'm_min_support', 'type': 'string', 'default': '1'},
            # 靶lncRNA软件 ['miranda', 'rnahybird', 'psrobot', 'targetfinder']
            {'name': 'l_method', 'type': 'string'},
            # miRanda score
            {'name': 'l_miranda_score', 'type': 'string', 'default': '160'},
            # miRanda energy
            {'name': 'l_miranda_energy', 'type': 'string', 'default': '-20'},
            # miRanda strict
            {'name': 'l_miranda_strict', 'type': 'string', 'default': 'on'},
            # RNAhybird num
            {'name': 'l_rnahybird_num', 'type': 'string', 'default': '100'},
            # RNAhybird energy
            {'name': 'l_rnahybird_energy', 'type': 'string', 'default': '-20'},
            # RNAhybird p-value
            {'name': 'l_rnahybird_pvalue', 'type': 'string', 'default': '0.01'},
            # psRobot score
            {'name': 'l_ps_robot_score', 'type': 'string', 'default': '2.5'},
            # Targetfinder score
            {'name': 'l_targetfinder_score', 'type': 'string', 'default': '4'},
            # 靶mRNA候选标准
            {'name': 'l_min_support', 'type': 'string', 'default': '1'},
            # 靶circRNA软件 ['miranda', 'rnahybird', 'psrobot', 'targetfinder']
            {'name': 'c_method', 'type': 'string'},
            # miRanda score
            {'name': 'c_miranda_score', 'type': 'string', 'default': '160'},
            # miRanda energy
            {'name': 'c_miranda_energy', 'type': 'string', 'default': '-20'},
            # miRanda strict
            {'name': 'c_miranda_strict', 'type': 'string', 'default': 'on'},
            # RNAhybird num
            {'name': 'c_rnahybird_num', 'type': 'string', 'default': '100'},
            # RNAhybird energy
            {'name': 'c_rnahybird_energy', 'type': 'string', 'default': '-20'},
            # RNAhybird p-value
            {'name': 'c_rnahybird_pvalue', 'type': 'string', 'default': '0.01'},
            # psRobot score
            {'name': 'c_ps_robot_score', 'type': 'string', 'default': '2.5'},
            # Targetfinder score
            {'name': 'c_targetfinder_score', 'type': 'string', 'default': '4'},
            # 靶mRNA候选标准
            {'name': 'c_min_support', 'type': 'string', 'default': '1'},
            {'name': 'report_img', 'type': 'bool', 'default': True},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.project_type = 'whole_transcriptome'
        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        self.modules = dict()
        self.tools = dict()

        # 用于在重运行时，删除已经导入到mongo库的表，避免数据重复
        # time.sleep(10)
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()
        # data = os.path.join(self.work_dir, 'data.json')
        # if os.path.exists(data):
        #     with open(data, 'r') as load_f:
        #         load_dict = json.load(load_f)
        #         if 'rerun' in load_dict and load_dict['rerun']:
        #             self.logger.info("该项目重运行中，先删除mongo库中已有数据")
        #             self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        # self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        # self.program = os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin/python')
        # cmd = '{} {}'.format(self.program, self.script)
        # cmd += ' {} {}'.format(self.task_id, 'whole_transcriptome')
        # code = os.system(cmd)
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))
        delete = DeleteDemoMongo(self.task_id, 'whole_transcriptome')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        if self.option("mirbase_category"):
            if self.option("mirbase_category").lower() == "plant":
                if not self.option("c_method"):
                    self.option("c_method", "psrobot")
                if not self.option("l_method"):
                    self.option("l_method", "psrobot")
                if not self.option("m_method"):
                    self.option("m_method", "psrobot")
            else:
                if not self.option("c_method"):
                    self.option("c_method", "miranda")
                if not self.option("l_method"):
                    self.option("l_method", "miranda")
                if not self.option("m_method"):
                    self.option("m_method", "miranda")
        else:
            if not self.option("c_method"):
                self.option("c_method", "miranda")
            if not self.option("l_method"):
                self.option("l_method", "miranda")
            if not self.option("m_method"):
                self.option("m_method", "miranda")
        for k, v in self.sheet.options().items():
            self.logger.debug('{} -> {}'.format(k, v))
        else:
            return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        self.add_steps()
        super(WholeTranscriptomeWorkflow, self).run()

    def add_steps(self):
        self.step.add_steps('transfer_l')
        self.step.add_steps('transfer_s')
        self.step.add_steps('transfer_c')
        self.step.add_steps('whole_snp')
        self.step.add_steps('rmats')
        self.step.add_steps('formation_mrna')
        self.step.add_steps('target_cistrans')
        self.step.add_steps('lncrna_family')
        self.step.add_steps('mirna_precursor')
        self.step.add_steps('formation_lncrna')
        self.step.add_steps('target_mirna')
        self.step.add_steps('atcg_bias')
        self.step.add_steps('mirna_edit')
        self.step.add_steps('smallrna_family_analyse')
        self.step.add_steps('formation_mirna')
        self.step.add_steps('formation_circrna')
        self.step.add_steps('formation_longrna')
        self.step.add_steps('diff_split')
        self.step.add_steps('diff_split_g')
        self.step.add_steps('gene_detail')
        self.step.add_steps('genesets_analysis')
        self.load_libraries()

    def load_libraries(self):
        self.tools['transfer_l'] = self.add_tool('whole_transcriptome.transfer')
        self.tools['transfer_s'] = self.add_tool('whole_transcriptome.transfer')
        self.tools['transfer_c'] = self.add_tool('whole_transcriptome.transfer')

        self.modules['whole_snp'] = self.add_module('whole_transcriptome.whole_snp')
        self.modules['rmats'] = self.add_module('whole_transcriptome.rmats')
        self.modules['formation_mrna'] = self.add_module('whole_transcriptome.formation')

        self.modules['target_cistrans'] = self.add_module('whole_transcriptome.target_cistrans')
        self.modules['lncrna_family'] = self.add_module('lnc_rna.lncrna_family')
        self.modules['mirna_precursor'] = self.add_module('lnc_rna.mirna_precursor')
        self.modules['formation_lncrna'] = self.add_module('whole_transcriptome.formation')

        self.modules['target_mirna'] = self.add_module('whole_transcriptome.target_mirna')
        self.tools['atcg_bias'] = self.add_tool('small_rna.atcg_bias')
        self.modules['mirna_edit'] = self.add_module('small_rna.mirna_edit')
        self.tools['smallrna_family_analyse'] = self.add_tool('small_rna.smallrna_family_analyse')
        self.modules['formation_mirna'] = self.add_module('whole_transcriptome.formation')

        self.modules['formation_circrna'] = self.add_module('whole_transcriptome.formation')

        self.modules['formation_longrna'] = self.add_module('whole_transcriptome.formation')
        self.tools['diff_split'] = self.add_tool('whole_transcriptome.formation.diff_split')
        self.tools['diff_split_g'] = self.add_tool('whole_transcriptome.formation.diff_split')
        self.modules['gene_detail'] = self.add_module('whole_transcriptome.gene_detail')
        self.modules['genesets_analysis'] = self.add_module('whole_transcriptome_v1_2.workflow_diffgt.diff_geneset_all_pipline')
        self.download_data()

    def download_data(self):
        database = Config().get_mongo_client(mtype=self.project_type)[Config().get_mongo_dbname(self.project_type)]
        collection = database['task']
        self.lib_dict = {'long': True, 'small': False, 'circle': False}
        if 'small' in self.option('lib_select'):
            self.lib_dict['small'] = True
        if 'circ' in self.option('lib_select'):
            self.lib_dict['circle'] = True
        relies = list()
        if self.lib_dict['long']:
            self.long_task_info = collection.find_one({'task_id': self.option('long_task_id')})
            relies.append(self.tools['transfer_l'])
        if self.lib_dict['small']:
            self.small_task_info = collection.find_one({'task_id': self.option('small_task_id')})
            relies.append(self.tools['transfer_s'])
        if self.lib_dict['circle']:
            self.circle_task_info = collection.find_one({'task_id': self.option('circ_task_id')})
            relies.append(self.tools['transfer_c'])
        self.on_rely(relies, self.enter_pipeline)
        for lib, func in (
                ('long', self.run_transfer_l), ('small', self.run_transfer_s), ('circle', self.run_transfer_c)):
            if self.lib_dict[lib]:
                func()

    def run_transfer_l(self):
        if self.long_task_info['output'].endswith("/"):
            indir = self.long_task_info['output']
        else:
            indir = self.long_task_info['output'] + "/"
        opts = {
            'indir': indir
        }
        self.tools['transfer_l'].set_options(opts)
        self.tools['transfer_l'].on('start', self.set_step, {'start': self.step.transfer_l})
        self.tools['transfer_l'].on('end', self.set_step, {'end': self.step.transfer_l})
        self.tools['transfer_l'].run()

    def run_transfer_s(self):
        if self.small_task_info['output'].endswith("/"):
            indir = self.small_task_info['output']
        else:
            indir = self.small_task_info['output'] + "/"
        self.logger.debug('indir {}'.format(indir))
        opts = {
            'indir': indir
        }
        self.tools['transfer_s'].set_options(opts)
        self.tools['transfer_s'].on('start', self.set_step, {'start': self.step.transfer_s})
        self.tools['transfer_s'].on('end', self.set_step, {'end': self.step.transfer_s})
        self.tools['transfer_s'].run()

    def run_transfer_c(self):
        if self.circle_task_info['output'].endswith("/"):
            indir = self.circle_task_info['output']
        else:
            indir = self.circle_task_info['output'] + "/"
        opts = {
            'indir': indir
        }
        self.tools['transfer_c'].set_options(opts)
        self.tools['transfer_c'].on('start', self.set_step, {'start': self.step.transfer_c})
        self.tools['transfer_c'].on('end', self.set_step, {'end': self.step.transfer_c})
        self.tools['transfer_c'].run()

    def enter_pipeline(self):
        if self.lib_dict['circle']:
            if not os.path.isdir(os.path.join(self.tools['transfer_c'].output_dir, 'exp_make')):
                shutil.copytree(os.path.join(self.tools['transfer_c'].output_dir, 'exp_makecirc'),
                                os.path.join(self.tools['transfer_c'].output_dir, 'exp_make'))
            if not os.path.isdir(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/circrna')):
                os.mkdir(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/circrna'))
            for fpath in ('exp_make/circrna/T.rpm.txt', 'exp_make/count/C.reads.txt'):
                shutil.copy(os.path.join(self.tools['transfer_c'].output_dir, fpath),
                            os.path.join(self.tools['transfer_l'].output_dir, fpath))
            for dname in ('circ_brush', 'diff_exp_c'):
                if os.path.isdir(os.path.join(self.tools['transfer_l'].output_dir, dname)):
                    shutil.rmtree(os.path.join(self.tools['transfer_l'].output_dir, dname))
                shutil.copytree(os.path.join(self.tools['transfer_c'].output_dir, dname),
                                os.path.join(self.tools['transfer_l'].output_dir, dname))

        # database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2')]
        database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[
            Config().get_mongo_dbname(mtype='ref_rna_v2', dydb_forbid=True)]
        collection = database['sg_genome_db']
        self.genome_doc = collection.find_one({'genome_id': self.long_task_info['options']['genome_id']})
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish')

        database = Config().get_mongo_client(mtype='small_rna', dydb_forbid=True)[
            Config().get_mongo_dbname(mtype='small_rna', dydb_forbid=True)]
        # database = Config().get_mongo_client(mtype='small_rna')[Config().get_mongo_dbname('small_rna')]
        collection = database['mirbase']
        if ',' in str(self.option('mirbase_specie')):
            mirbase_specie = self.option('mirbase_specie').strip(",").split(',')[0]
        else:
            mirbase_specie = self.option('mirbase_specie')
        self.mirbase_doc = collection.find_one({'organism': mirbase_specie})

        self.relies = list()
        self.funcs = list()

        if 'mRNA' in self.long_task_info['options']['rna_select']:
            self.set_mrna()
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            self.set_lncrna()
        if self.lib_dict['small']:
            self.set_smallrna()
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            self.set_circrna()
        self.set_common()

        # self.on_rely(self.relies, self.set_db)
        self.on_rely(self.relies, self.run_chart)
        for func in self.funcs:
            func()

    def set_mrna(self):
        if self.option('is_snp'):
            self.relies.append(self.modules['whole_snp'])
            self.funcs.append(self.run_whole_snp)
        if self.option('is_as'):
            self.relies.append(self.modules['rmats'])
            self.funcs.append(self.run_rmats)
        self.relies.append(self.modules['formation_mrna'])
        self.funcs.append(self.run_formation_mrna)
        self.relies.append(self.modules['formation_longrna'])
        self.funcs.append(self.run_formation_longrna)
        self.relies.append(self.tools['diff_split'])
        self.relies.append(self.tools['diff_split_g'])
        self.relies.append(self.modules['genesets_analysis'])
        self.funcs.append(self.run_diff_split)
        self.funcs.append(self.run_diff_split_gene)
        # self.funcs.append(self.run_geneset_analysis)

    def set_lncrna(self):
        self.relies.append(self.modules['target_cistrans'])
        # self.funcs.append(self.run_target_cistrans)
        self.relies.append(self.modules['lncrna_family'])
        self.funcs.append(self.run_lncrna_family)
        self.relies.append(self.modules['mirna_precursor'])
        self.funcs.append(self.run_mirna_precursor)
        self.relies.append(self.modules['formation_lncrna'])
        self.funcs.append(self.run_formation_lncrna)

    def set_smallrna(self):
        self.relies.append(self.modules['target_mirna'])
        # self.funcs.append(self.run_target_mirna)
        self.relies.append(self.tools['atcg_bias'])
        self.funcs.append(self.run_atcg_bias)
        self.relies.append(self.modules['mirna_edit'])
        self.funcs.append(self.run_mirna_edit)
        self.relies.append(self.tools['smallrna_family_analyse'])
        self.funcs.append(self.run_smallrna_family_analyse)
        self.relies.append(self.modules['formation_mirna'])
        self.funcs.append(self.run_formation_mirna)

    def set_circrna(self):
        self.relies.append(self.modules['formation_circrna'])
        self.funcs.append(self.run_formation_circrna)

    def set_common(self):
        self.relies.append(self.modules['gene_detail'])
        self.funcs.append(self.run_gene_detail)

    def run_gene_detail(self):
        rnas = list()
        opts = dict()
        if 'mRNA' in self.long_task_info['options']['rna_select']:
            rnas.append('mrna')
            gene_type = os.path.join(self.tools['transfer_l'].output_dir,
                                    'large_gush/filter_by_express/filtered_file/gene_type.xls')
            trans_type = os.path.join(self.tools['transfer_l'].output_dir,
                                      'large_gush/filter_by_express/filtered_file/trans_type.xls')
            mrna_gtf = os.path.join(self.tools['transfer_l'].output_dir,
                                    'large_gush/filter_by_express/filtered_file/all_mrna.gtf')
            all_gtf = os.path.join(self.tools['transfer_l'].output_dir,
                                    'large_gush/filter_by_express/filtered_file/all.gtf')
            dna_fa = os.path.join(self.db_path, self.genome_doc['dna_fa'])
            relation_file = os.path.join(self.tools['transfer_l'].output_dir,
                                         'annotation/allannot_class/all_tran2gene.txt')
            biomart_file = os.path.join(self.db_path, self.genome_doc['bio_mart_annot'])
            biomart_type = self.genome_doc['biomart_gene_annotype']
            ref_cds = os.path.join(self.db_path, self.genome_doc['cds'])
            ref_pep = os.path.join(self.db_path, self.genome_doc['pep'])
            new_cds = os.path.join(self.tools['transfer_l'].output_dir,
                                   'annotation/newannot_orfpfam/novel_mrna.fa.transdecoder.cds')
            new_pep = os.path.join(self.tools['transfer_l'].output_dir,
                                   'annotation/newannot_orfpfam/novel_mrna.fa.transdecoder.pep')
            genome_id = self.genome_doc['genome_id']
            opts.update({
                'trans_type': trans_type,
                'gene_type': gene_type,
                'mrna_gtf': mrna_gtf,
                'all_gtf': all_gtf,
                'dna_fa': dna_fa,
                'relation_file': relation_file,
                'biomart_file': biomart_file,
                'biomart_type': biomart_type,
                'ref_cds': ref_cds,
                'ref_pep': ref_pep,
                'new_cds': new_cds,
                'new_pep': new_pep,
                'genome_id': genome_id})
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            rnas.append('lncrna')
            lncrna_gtf = os.path.join(self.tools['transfer_l'].output_dir,
                                      'large_gush/filter_by_express/filtered_file/all_lncrna.gtf')
            known_lncrna_detail = os.path.join(self.tools['transfer_l'].output_dir,
                                               'large_gush/known_lnc_identify/known_lncrna_detail.xls')
            novel_lncrna_detail = os.path.join(self.tools['transfer_l'].output_dir,
                                               'large_gush/new_lncrna_predict/novel_lncrna_predict_detail.xls')
            lnc_relation_file = os.path.join(self.tools['transfer_l'].output_dir,
                                             'large_gush/filter_by_express/filtered_file/trans_type.xls')
            opts.update({
                'lncrna_gtf': lncrna_gtf,
                'known_lncrna_detail': known_lncrna_detail,
                'novel_lncrna_detail': novel_lncrna_detail,
                'lnc_relation_file': lnc_relation_file,
            })
        if self.lib_dict['small']:
            rnas.append('mirna')
            known_mirna_detail = os.path.join(self.tools['transfer_s'].output_dir,
                                              'srna/known_mirna/known_mirna_detail.xls')
            novel_mirna_detail = os.path.join(self.tools['transfer_s'].output_dir,
                                              'srna/novel_mirna/novel_mirna_detail.xls')
            opts.update({
                'known_mirna_detail': known_mirna_detail,
                'novel_mirna_detail': novel_mirna_detail
            })
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            rnas.append('circrna')
            if self.lib_dict['circle']:
                circrna_detail = os.path.join(self.tools['transfer_c'].output_dir, 'circ_brush/detail.txt')
            else:
                circrna_detail = os.path.join(self.tools['transfer_l'].output_dir, 'circ_brush/detail.txt')
            opts.update({'circrna_detail': circrna_detail})
        opts.update({'rna_type': ','.join(rnas)})
        self.modules['gene_detail'].set_options(opts)
        self.modules['gene_detail'].on('start', self.set_step, {'start': self.step.gene_detail})
        self.modules['gene_detail'].on('end', self.set_step, {'end': self.step.gene_detail})
        self.modules['gene_detail'].on('end', self.set_output, 'gene_detail')
        self.modules['gene_detail'].run()

    def run_formation_circrna(self):
        if self.lib_dict['circle']:
            key = 'transfer_c'
        else:
            key = 'transfer_l'
        exp_matrix = os.path.join(self.tools[key].output_dir, 'exp_make/circrna/T.rpm.txt')
        group_table = os.path.join(self.tools[key].output_dir, 'group.txt')
        calls = 'graph,venn,corr,pca'
        graph_log = False
        threshold = 0.0
        opts = {
            'exp_matrix': exp_matrix,
            'group_table': group_table,
            'calls': calls,
            'graph_log': graph_log,
            'threshold': threshold,
        }
        self.modules['formation_circrna'].set_options(opts)
        self.modules['formation_circrna'].on('start', self.set_step, {'start': self.step.formation_circrna})
        self.modules['formation_circrna'].on('end', self.set_step, {'end': self.step.formation_circrna})
        self.modules['formation_circrna'].on('end', self.set_output, 'formation_circrna')
        self.modules['formation_circrna'].run()

    def run_target_mirna(self):
        s_id_set = set(
            pd.read_table(os.path.join(self.tools['transfer_s'].output_dir, 'diff_exp/summary.txt'))['seq_id'])
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            l_id_set = set(
                pd.read_table(os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_t/summary.txt'))['seq_id'])
        else:
            l_id_set = set(
                pd.read_table(os.path.join(self.tools['diff_split'].output_dir, 'mrna/summary.txt'))['seq_id'])
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            c_id_set = set(
                pd.read_table(os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_c/summary.txt'))['seq_id'])
        else:
            c_id_set = set()

        target_dir = os.path.join(self.modules['target_mirna'].work_dir, 'target')
        novol = os.path.join(target_dir, 'new.mirna.fasta')
        known = os.path.join(target_dir, 'ref.mirna.fasta')
        m_novol = os.path.join(target_dir, 'new.mrna.fasta')
        m_known = os.path.join(target_dir, 'ref.mrna.fasta')
        l_novol = os.path.join(target_dir, 'new.lncrna.fasta')
        l_known = os.path.join(target_dir, 'ref.lncrna.fasta')
        c_novol = os.path.join(target_dir, 'new.circrna.fasta')
        c_known = os.path.join(target_dir, 'ref.circrna.fasta')

        if not os.path.isdir(target_dir):
            os.mkdir(target_dir)
            SeqIO.write([record for record in SeqIO.parse(
                os.path.join(self.tools['transfer_s'].output_dir, 'srna/novel_mirna/novel_mature_seq.fa'), 'fasta') if
                         record.id in s_id_set], novol, 'fasta')
            SeqIO.write([record for record in SeqIO.parse(
                os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna/mature.fa'), 'fasta') if
                         record.id in s_id_set], known, 'fasta')
            SeqIO.write([record for record in SeqIO.parse(
                os.path.join(self.tools['transfer_l'].output_dir,
                             'large_gush/filter_by_express/filtered_file/novel_mrna.fa'), 'fasta') if
                         record.id in l_id_set], m_novol, 'fasta')
            SeqIO.write([record for record in SeqIO.parse(
                os.path.join(self.tools['transfer_l'].output_dir,
                             'large_gush/filter_by_express/filtered_file/known_mrna.fa'), 'fasta') if
                         record.id in l_id_set], m_known, 'fasta')
            SeqIO.write([record for record in SeqIO.parse(
                os.path.join(self.tools['transfer_l'].output_dir,
                             'large_gush/filter_by_express/filtered_file/novel_lncrna.fa'), 'fasta') if
                         record.id in l_id_set], l_novol, 'fasta')
            SeqIO.write([record for record in SeqIO.parse(
                os.path.join(self.tools['transfer_l'].output_dir,
                             'large_gush/filter_by_express/filtered_file/known_lncrna.fa'), 'fasta') if
                         record.id in l_id_set], l_known, 'fasta')
            if c_id_set:
                SeqIO.write([record for record in SeqIO.parse(
                    os.path.join(self.tools['transfer_l'].output_dir, 'circ_brush/new.fasta'), 'fasta') if
                             record.id in c_id_set], c_novol, 'fasta')
                SeqIO.write([record for record in SeqIO.parse(
                    os.path.join(self.tools['transfer_l'].output_dir, 'circ_brush/ref.fasta'), 'fasta') if
                             record.id in c_id_set], c_known, 'fasta')
            else:
                open(c_novol, 'w').close()
                open(c_known, 'w').close()

        anno_detail = os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/all_annot.xls')
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            circ_detail = os.path.join(self.tools['transfer_l'].output_dir, 'circ_brush/detail.txt')
        else:
            circ_detail = os.path.join(self.modules['target_mirna'].work_dir, 'circ_detail.txt')
            open(circ_detail, 'w').close()
        lnc_ref_detail = os.path.join(self.tools['transfer_l'].output_dir,
                                      'large_gush/known_lnc_identify/known_lncrna_detail.xls')
        lnc_new_detail = os.path.join(self.tools['transfer_l'].output_dir,
                                      'large_gush/new_lncrna_predict/novel_lncrna_predict_detail.xls')

        taxonomy = self.small_task_info['options']['taxonmy']

        m_method = self.option('m_method')
        m_miranda_score = self.option('m_miranda_score')
        m_miranda_energy = self.option('m_miranda_energy')
        m_miranda_strict = self.option('m_miranda_strict')
        m_rnahybird_num = self.option('m_rnahybird_num')
        m_rnahybird_energy = self.option('m_rnahybird_energy')
        m_rnahybird_pvalue = self.option('m_rnahybird_pvalue')
        m_ps_robot_score = self.option('m_ps_robot_score')
        m_targetfinder_score = self.option('m_targetfinder_score')
        m_min_support = self.option('m_min_support')

        l_method = self.option('l_method')
        l_miranda_score = self.option('l_miranda_score')
        l_miranda_energy = self.option('l_miranda_energy')
        l_miranda_strict = self.option('l_miranda_strict')
        l_rnahybird_num = self.option('l_rnahybird_num')
        l_rnahybird_energy = self.option('l_rnahybird_energy')
        l_rnahybird_pvalue = self.option('l_rnahybird_pvalue')
        l_ps_robot_score = self.option('l_ps_robot_score')
        l_targetfinder_score = self.option('l_targetfinder_score')
        l_min_support = self.option('l_min_support')

        c_method = self.option('c_method')
        c_miranda_score = self.option('c_miranda_score')
        c_miranda_energy = self.option('c_miranda_energy')
        c_miranda_strict = self.option('c_miranda_strict')
        c_rnahybird_num = self.option('c_rnahybird_num')
        c_rnahybird_energy = self.option('c_rnahybird_energy')
        c_rnahybird_pvalue = self.option('c_rnahybird_pvalue')
        c_ps_robot_score = self.option('c_ps_robot_score')
        c_targetfinder_score = self.option('c_targetfinder_score')
        c_min_support = self.option('c_min_support')

        opts = {
            'novol': novol,
            'known': known,
            'm_novol': m_novol,
            'm_known': m_known,
            'l_novol': l_novol,
            'l_known': l_known,
            'c_novol': c_novol,
            'c_known': c_known,
            'taxonomy': taxonomy,
            'anno_detail': anno_detail,
            'circ_detail': circ_detail,
            'lnc_ref_detail': lnc_ref_detail,
            'lnc_new_detail': lnc_new_detail,

            'm_method': m_method,
            'm_miranda_score': m_miranda_score,
            'm_miranda_energy': m_miranda_energy,
            'm_miranda_strict': m_miranda_strict,
            'm_rnahybird_num': m_rnahybird_num,
            'm_rnahybird_energy': m_rnahybird_energy,
            'm_rnahybird_pvalue': m_rnahybird_pvalue,
            'm_ps_robot_score': m_ps_robot_score,
            'm_targetfinder_score': m_targetfinder_score,
            'm_min_support': m_min_support,

            'l_method': l_method,
            'l_miranda_score': l_miranda_score,
            'l_miranda_energy': l_miranda_energy,
            'l_miranda_strict': l_miranda_strict,
            'l_rnahybird_num': l_rnahybird_num,
            'l_rnahybird_energy': l_rnahybird_energy,
            'l_rnahybird_pvalue': l_rnahybird_pvalue,
            'l_ps_robot_score': l_ps_robot_score,
            'l_targetfinder_score': l_targetfinder_score,
            'l_min_support': l_min_support,

            'c_method': c_method,
            'c_miranda_score': c_miranda_score,
            'c_miranda_energy': c_miranda_energy,
            'c_miranda_strict': c_miranda_strict,
            'c_rnahybird_num': c_rnahybird_num,
            'c_rnahybird_energy': c_rnahybird_energy,
            'c_rnahybird_pvalue': c_rnahybird_pvalue,
            'c_ps_robot_score': c_ps_robot_score,
            'c_targetfinder_score': c_targetfinder_score,
            'c_min_support': c_min_support
        }
        opts = {k: str(v) for k, v in opts.items()}
        self.modules['target_mirna'].set_options(opts)
        self.modules['target_mirna'].on('end', self.set_output, 'target_mirna')
        self.modules['target_mirna'].on('start', self.set_step, {'start': self.step.target_mirna})
        self.modules['target_mirna'].on('end', self.set_step, {'end': self.step.target_mirna})
        self.modules['target_mirna'].run()

    def run_atcg_bias(self):
        known = os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna/mature.fa')
        novel = os.path.join(self.tools['transfer_s'].output_dir, 'srna/novel_mirna/novel_mature_seq.fa')
        opts = {'known': known, 'novel': novel}
        self.tools['atcg_bias'].set_options(opts)
        self.tools['atcg_bias'].on('end', self.set_output, 'atcg_bias')
        self.tools['atcg_bias'].on('start', self.set_step, {'start': self.step.atcg_bias})
        self.tools['atcg_bias'].on('end', self.set_step, {'end': self.step.atcg_bias})
        self.tools['atcg_bias'].run()

    def run_mirna_edit(self):
        list_file = os.path.join(self.tools['transfer_s'].output_dir, 'mirna_qc/clean_data/list.txt')
        lines = list()
        for line in open(list_file):
            sample = line.strip().split('\t')[1]
            lines.append('{}\t{}\n'.format(os.path.join(self.tools['transfer_s'].output_dir,
                                                        'mirna_qc/clean_data/{}_clip_s.fastq.trimmed'.format(sample)),
                                           sample))
        open(list_file, 'w').writelines(lines)
        if self.small_task_info['version'] >= 'v1.2':
            species = self.small_task_info['options']['mirna_specie']
        else:
            species = self.small_task_info['options']['mirbase_specie']
        if ',' in species:
            species = species.split(',')[0]
        hairpin_fa = os.path.join(self.config.SOFTWARE_DIR, 'database/mirbase/hairpin.fa')
        mature_fa = os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna/mature.fa')
        index = os.path.join(self.db_path, self.genome_doc['dna_index'])
        opts = {
            'list_file': list_file,
            'species': species,
            'hairpin_fa': hairpin_fa,
            'mature_fa': mature_fa,
            'index': index
        }
        self.modules['mirna_edit'].set_options(opts)
        self.modules['mirna_edit'].on('end', self.set_output, 'mirna_edit')
        self.modules['mirna_edit'].on('start', self.set_step, {'start': self.step.mirna_edit})
        self.modules['mirna_edit'].on('end', self.set_step, {'end': self.step.mirna_edit})
        self.modules['mirna_edit'].run()

    def run_smallrna_family_analyse(self):
        mir = os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna_count.xls')
        matfa = os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna/mature.fa')
        novofa = os.path.join(self.tools['transfer_s'].output_dir, 'srna/novel_mirna/novel_mature_seq.fa')
        opts = {'mir': mir, 'matfa': matfa, 'novofa': novofa}
        self.tools['smallrna_family_analyse'].set_options(opts)
        self.tools['smallrna_family_analyse'].on('end', self.set_output, 'smallrna_family_analyse')
        self.tools['smallrna_family_analyse'].on('start', self.set_step, {'start': self.step.smallrna_family_analyse})
        self.tools['smallrna_family_analyse'].on('end', self.set_step, {'end': self.step.smallrna_family_analyse})
        self.tools['smallrna_family_analyse'].run()

    def run_formation_mirna(self):
        if os.path.isdir(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mirna')):
            shutil.rmtree(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mirna'))
        os.makedirs(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mirna'))
        count_df = pd.read_table(os.path.join(self.tools['transfer_s'].output_dir, 'srna/total_mirna_count.xls'),
                                 index_col=0)
        count_df.index.name = 'seq_id'
        count_df.to_csv(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/S.reads.txt'), sep='\t')
        ref_df = pd.read_table(os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna_norm.xls'),
                               index_col=0)
        ref_df['level'] = 'T'
        ref_df['category'] = 'miRNA'
        ref_df['kind'] = 'ref'
        new_df = pd.read_table(os.path.join(self.tools['transfer_s'].output_dir, 'srna/novel_mirna_norm.xls'),
                               index_col=0)
        new_df['level'] = 'T'
        new_df['category'] = 'miRNA'
        new_df['kind'] = 'new'
        all_df = pd.concat([ref_df, new_df])
        all_df.index.name = 'transcript_id'
        all_df.to_csv(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mirna/T.tpm.txt'), sep='\t')
        exp_matrix = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mirna/T.tpm.txt')
        group_table = os.path.join(self.tools['transfer_s'].output_dir, 'group.txt')
        calls = 'graph,venn,corr,pca'
        opts = {
            'exp_matrix': exp_matrix,
            'group_table': group_table,
            'calls': calls,
        }
        self.modules['formation_mirna'].set_options(opts)
        self.modules['formation_mirna'].on('start', self.set_step, {'start': self.step.formation_mirna})
        self.modules['formation_mirna'].on('end', self.set_step, {'end': self.step.formation_mirna})
        self.modules['formation_mirna'].on('end', self.set_output, 'formation_mirna')
        self.modules['formation_mirna'].run()

    def run_mirna_precursor(self):
        lncrna_fa = os.path.join(self.tools['transfer_l'].output_dir,
                                 'large_gush/filter_by_express/filtered_file/all_lncrna.fa')
        species = self.mirbase_doc['organism']
        known_list = os.path.join(self.tools['transfer_l'].output_dir,
                                  'large_gush/filter_by_express/filtered_file/known_lncrna_ids.list')
        novel_list = os.path.join(self.tools['transfer_l'].output_dir,
                                  'large_gush/filter_by_express/filtered_file/novel_lncrna_ids.list')
        t2g = os.path.join(self.tools['transfer_l'].output_dir,
                           'large_gush/filter_by_express/filtered_file/t2g.txt')
        opts = {
            'lncrna_fa': lncrna_fa,
            'species': species,
            'known_list': known_list,
            'novel_list': novel_list,
            't2g': t2g
        }
        self.modules['mirna_precursor'].set_options(opts)
        self.modules['mirna_precursor'].on('end', self.set_output, 'mirna_precursor')
        self.modules['mirna_precursor'].on('start', self.set_step, {'start': self.step.mirna_precursor})
        self.modules['mirna_precursor'].on('end', self.set_step, {'end': self.step.mirna_precursor})
        self.modules['mirna_precursor'].run()

    def run_lncrna_family(self):
        lncrna_fa = os.path.join(self.tools['transfer_l'].output_dir,
                                 'large_gush/filter_by_express/filtered_file/all_lncrna.fa')
        known_list = os.path.join(self.tools['transfer_l'].output_dir,
                                  'large_gush/filter_by_express/filtered_file/known_lncrna_ids.list')
        novel_list = os.path.join(self.tools['transfer_l'].output_dir,
                                  'large_gush/filter_by_express/filtered_file/novel_lncrna_ids.list')
        t2g = os.path.join(self.tools['transfer_l'].output_dir,
                           'large_gush/filter_by_express/filtered_file/t2g.txt')
        opts = {'lncrna_fa': lncrna_fa, 'known_list': known_list, 'novel_list': novel_list, 't2g': t2g}
        self.modules['lncrna_family'].set_options(opts)
        self.modules['lncrna_family'].on('end', self.set_output, 'lncrna_family')
        self.modules['lncrna_family'].on('start', self.set_step, {'start': self.step.lncrna_family})
        self.modules['lncrna_family'].on('end', self.set_step, {'end': self.step.lncrna_family})
        self.modules['lncrna_family'].run()

    def run_target_cistrans(self):
        known = os.path.join(self.tools['transfer_l'].output_dir,
                             'large_gush/filter_by_express/filtered_file/known_lncrna.gtf')
        novol = os.path.join(self.tools['transfer_l'].output_dir,
                             'large_gush/filter_by_express/filtered_file/novel_lncrna.gtf')
        mrna_gtf = os.path.join(self.tools['transfer_l'].output_dir,
                                'large_gush/filter_by_express/filtered_file/all_mrna.gtf')
        annotation = os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/all_annot.xls')
        up_dis = self.option('up_dis')
        down_dis = self.option('down_dis')
        exp_way = self.long_task_info['options']['exp_way']

        exp_matrix_lnc = os.path.join(self.tools['transfer_l'].output_dir,
                                      'large_gush/filter_by_express/classifyquant/T.{}.txt'.format(exp_way))
        exp_matrix_target = os.path.join(self.tools['transfer_l'].output_dir,
                                         'large_gush/filter_by_express/classifyquant/G.{}.txt'.format(exp_way))
        padjust_way = self.option('lt_adj_way').lower()
        filter_lnc_file = os.path.join(self.tools['diff_split'].output_dir, 'lncrna/summary.txt')
        filter_target_file = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_g/summary.txt')
        opts = {
            'known': known,
            'novol': novol,
            'mrna_gtf': mrna_gtf,
            'annotation': annotation,
            'up_dis': up_dis,
            'down_dis': down_dis,
            'exp_matrix_lnc': exp_matrix_lnc,
            'exp_matrix_target': exp_matrix_target,
            'cor_cutoff': self.option('lt_cor_cut'),
            'corr_way': self.option('lt_cor_way'),
            'padjust_way': padjust_way,
            'filter_lnc_file': filter_lnc_file,
            'filter_target_file': filter_target_file
        }
        if self.option('lt_sig_way') == 'padjust':
            opts['qvalue_cutoff'] = 0.05
        elif self.option('lt_sig_way') == 'pvalue':
            opts['pvalue_cutoff'] = 0.05
        self.modules['target_cistrans'].set_options(opts)
        self.modules['target_cistrans'].on('end', self.set_output, 'target_cistrans')
        self.modules['target_cistrans'].on('start', self.set_step, {'start': self.step.target_cistrans})
        self.modules['target_cistrans'].on('end', self.set_step, {'end': self.step.target_cistrans})
        self.modules['target_cistrans'].run()

    def run_formation_lncrna(self):
        exp_way = self.long_task_info['options']['exp_way']
        exp_matrix = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/lncrna/T.{}.txt'.format(exp_way))
        group_table = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        calls = 'graph,venn'
        opts = {
            'exp_matrix': exp_matrix,
            'group_table': group_table,
            'calls': calls,
        }
        self.modules['formation_lncrna'].set_options(opts)
        self.modules['formation_lncrna'].on('start', self.set_step, {'start': self.step.formation_lncrna})
        self.modules['formation_lncrna'].on('end', self.set_step, {'end': self.step.formation_lncrna})
        self.modules['formation_lncrna'].on('end', self.set_output, 'formation_lncrna')
        self.modules['formation_lncrna'].run()

    def run_whole_snp(self):
        organism_name = self.long_task_info['options']['organism_name']
        # genome_version = self.genome_doc['assembly']
        annot_version = self.genome_doc['annot_version']
        call_type = self.option('snp_method').lower()
        in_bam = os.path.join(self.tools['transfer_l'].output_dir, 'rnaseq_mapping/bam')
        bam_list = os.path.join(self.tools['transfer_l'].output_dir, 'rnaseq_mapping/bamlist')
        ref_genome_custom = os.path.join(self.db_path, self.genome_doc['dna_fa'])
        des = os.path.join(self.db_path, self.genome_doc['bio_mart_annot'])
        des_type = self.genome_doc['biomart_gene_annotype']
        ref_gtf = os.path.join(self.tools['transfer_l'].output_dir,
                               'large_gush/filter_by_express/filtered_file/all_mrna.gtf')
        opts = {
            'in_bam': in_bam,
            'ref_genome_custom': ref_genome_custom,
            'des': des,
            'des_type': des_type,
            'ref_gtf': ref_gtf,
            'align_method': self.long_task_info['options']['align_method'],
            'bam_list': bam_list
        }
        large_chr = {'Ginkgo_biloba': ['v1.0'], 'Triticum_turgidum': ['ensembl_45', 'iwgsc_refseq'],'Triticum_aestivum':['ensembl_48']}
        if organism_name in large_chr and annot_version in large_chr[organism_name]:
            opts.update({'analysis_format': 'cram'})
        self.modules['whole_snp'].set_options(opts)
        self.modules['whole_snp'].on('start', self.set_step, {'start': self.step.whole_snp})
        self.modules['whole_snp'].on('end', self.set_step, {'end': self.step.whole_snp})
        self.modules['whole_snp'].on('end', self.set_output, 'whole_snp')
        self.modules['whole_snp'].run()

    def run_rmats(self):
        control_table = os.path.join(self.tools['transfer_l'].output_dir, 'control.txt')
        group_table = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        bam_input = os.path.join(self.tools['transfer_l'].output_dir, 'rnaseq_mapping/bam')
        input_gtf = os.path.join(self.tools['transfer_l'].output_dir,
                                 'large_gush/filter_by_express/filtered_file/all_mrna.gtf')
        if self.long_task_info['options']['fq_type'] == 'PE':
            seq_type = 'paired'
        elif self.long_task_info['options']['fq_type'] == 'SE':
            seq_type = 'single'
        if self.long_task_info['options']['strand_specific'] == 'True':
            if self.long_task_info['options']['strand_dir'].startswith('R'):
                lib_type = 'fr-firststrand'
            elif self.long_task_info['options']['strand_dir'].startswith('F'):
                lib_type = 'fr-secondstrand'
        else:
            lib_type = 'fr-unstranded'
        opts = {
            'control_table': control_table,
            'group_table': group_table,
            'bam_input': bam_input,
            'input_gtf': input_gtf,
            'seq_type': seq_type,
            'lib_type': lib_type
        }
        self.modules['rmats'].set_options(opts)
        self.modules['rmats'].on('start', self.set_step, {'start': self.step.rmats})
        self.modules['rmats'].on('end', self.set_step, {'end': self.step.rmats})
        self.modules['rmats'].on('end', self.set_output, 'rmats')
        self.modules['rmats'].run()

    def run_formation_mrna(self):
        exp_way = self.long_task_info['options']['exp_way']
        exp_matrix = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/T.{}.txt'.format(exp_way))
        group_table = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        calls = 'graph,venn'
        opts = {
            'exp_matrix': exp_matrix,
            'group_table': group_table,
            'calls': calls,
        }
        self.modules['formation_mrna'].set_options(opts)
        self.modules['formation_mrna'].on('start', self.set_step, {'start': self.step.formation_mrna})
        self.modules['formation_mrna'].on('end', self.set_step, {'end': self.step.formation_mrna})
        self.modules['formation_mrna'].on('end', self.set_output, 'formation_mrna')
        self.modules['formation_mrna'].run()

    def run_formation_longrna(self):
        exp_way = self.long_task_info['options']['exp_way']
        exp_matrix = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/longrna/T.{}.txt'.format(exp_way))
        group_table = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        calls = 'corr,pca'
        opts = {
            'exp_matrix': exp_matrix,
            'group_table': group_table,
            'calls': calls,
        }
        self.modules['formation_longrna'].set_options(opts)
        self.modules['formation_longrna'].on('start', self.set_step, {'start': self.step.formation_longrna})
        self.modules['formation_longrna'].on('end', self.set_step, {'end': self.step.formation_longrna})
        self.modules['formation_longrna'].on('end', self.set_output, 'formation_longrna')
        self.modules['formation_longrna'].run()

    def run_diff_split(self):
        diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_t')
        relation = os.path.join(self.tools['transfer_l'].output_dir,
                                'large_gush/filter_by_express/filtered_file/trans_type.xls')
        opts = {'diff_dir': diff_dir, 'relation': relation}
        self.tools['diff_split'].set_options(opts)
        self.tools['diff_split'].on('start', self.set_step, {'start': self.step.diff_split})
        self.tools['diff_split'].on('end', self.set_step, {'end': self.step.diff_split})
        self.tools['diff_split'].on('end', self.set_output, 'diff_split')
        self.tools['diff_split'].on('end', self.run_geneset_analysis)
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            self.tools['diff_split'].on('end', self.run_target_cistrans)
        if self.lib_dict['small']:
            self.tools['diff_split'].on('end', self.run_target_mirna)
        self.tools['diff_split'].run()

    def run_diff_split_gene(self):
        diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_g')
        relation = os.path.join(self.tools['transfer_l'].output_dir,
                                'large_gush/filter_by_express/filtered_file/gene_type.xls')
        opts = {'diff_dir': diff_dir, 'relation': relation}
        self.tools['diff_split_g'].set_options(opts)
        self.tools['diff_split_g'].on('start', self.set_step, {'start': self.step.diff_split_g})
        self.tools['diff_split_g'].on('end', self.set_step, {'end': self.step.diff_split_g})
        self.tools['diff_split_g'].on('end', self.set_output, 'diff_split_g')
        self.tools['diff_split_g'].run()

    def run_geneset_analysis(self):
        diff_dir = os.path.join(self.tools['diff_split'].output_dir, 'mrna')
        exp = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/T.tpm.txt')
        annot_result = os.path.join(self.tools['transfer_l'].output_dir, 'annotation')
        group_file = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        level ="T"
        if "database_version" in self.long_task_info :
            kegg_version = self.long_task_info["database_version"].get("kegg", "2018")
        else:
            kegg_version =None
        species = self.long_task_info['options']['organism_name']
        opts = {
            'diff_path': diff_dir,
            'annot_result': annot_result,
            'diff_method': self.option('diff_method'),
            'transcript_exp_file': exp,
            'group': group_file,
            'level': level,
            "kegg_version": kegg_version,
            'species': species,
        }
        self.modules['genesets_analysis'].set_options(opts)
        self.modules['genesets_analysis'].on('start', self.set_step, {'start': self.step.genesets_analysis})
        self.modules['genesets_analysis'].on('end', self.set_step, {'end': self.step.genesets_analysis})
        self.modules['genesets_analysis'].on('end', self.set_output, 'genesets_analysis')
        self.modules['genesets_analysis'].run()


    def set_output(self, event):
        obj = event['bind_object']
        self.move2outputdir(obj.output_dir, event['data'])

    def move2outputdir(self, olddir, newname):
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])

    def move_file(self, src, dst):
        if os.path.isfile(src):
            os.link(src, dst)
        else:
            os.mkdir(dst)
            for file in os.listdir(src):
                old_path = os.path.join(src, file)
                new_path = os.path.join(dst, file)
                self.move_file(old_path, new_path)
        self.logger.debug('succeed in linking {} to {}'.format(src, dst))


    def set_db(self):
        self.pre_export()
        self.export_geneset()
        self.export_genome_info()
        self.export_qc_l()
        self.export_qc_s()
        self.export_qc_c()
        self.export_mapping_l()
        self.export_mapping_s()
        self.export_mapping_c()
        self.export_assembly()
        self.export_annotation()
        self.export_exp_detail()
        self.export_exp_graph()
        self.export_exp_venn()
        self.export_exp_corr()
        self.export_exp_pca()
        self.export_gensets_analysis()
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            self.export_lnc_identify()
            self.export_lnc_target()
            self.export_lncrna_family()
            self.export_mirna_precursor()
        self.export_srna()
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            self.export_circrna()
        self.export_diff()
        self.export_small_target()
        if self.option('is_as'):
            self.export_rmats()
            self.export_rmats_count()
        if self.option('is_snp'):
            self.export_snp()
        self.export_atcg_bias()
        self.export_mirna_edit()
        self.export_smallrna_family_analyse()
        self.export_gene_detail()
        self.export_task()
        # self.option("report_img",True)
        if self.option("report_img"):
            self.export_report_img()
        self.set_upload()

    def pre_export(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False

    def export_genome_info(self):
        api = self.api.api('whole_transcriptome.genome_info')
        file_path = os.path.join(self.db_path, self.genome_doc['gene_stat'])
        species_name = self.long_task_info['options']['organism_name']
        ref_anno_version = self.genome_doc['assembly']
        hyperlink = self.genome_doc['ensemble_web']
        api.add_genome_info(file_path, species_name, ref_anno_version, hyperlink)

    def export_qc_l(self):
        api = self.api.api('whole_transcriptome.qc')
        sample_list = os.path.join(self.tools['transfer_l'].output_dir, 'fastp_rna/fastq/list.txt')
        group_file = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        compare_file = os.path.join(self.tools['transfer_l'].output_dir, 'control.txt')
        if os.path.exists(os.path.join(self.tools['transfer_l'].output_dir, 'productive_table.txt')):
            productive_table = os.path.join(self.tools['transfer_l'].output_dir, 'productive_table.txt')
        else:
            productive_table = None
        fq_type = self.long_task_info['options']['fq_type']
        qc_stat_before = os.path.join(self.tools['transfer_l'].output_dir, 'hiseq_reads_stat_raw')
        qc_stat_after = os.path.join(self.tools['transfer_l'].output_dir, 'hiseq_reads_stat_use')
        api.add_sample_info(sample_list=sample_list, library='long', group_file=group_file, productive_table=productive_table)
        group_id, specimen_names, category_names = api.add_sample_group(group_file=group_file, library='long')
        control_id, compare_names = api.add_group_compare(compare_file=compare_file, library='long', group_id=group_id)
        qc_id = api.add_qc(fq_type=fq_type, library='long')
        api.add_qc_detail(qc_id, qc_stat_before, qc_stat_after, 'long', group=group_file)
        api.add_qc_graph(qc_id, qc_stat_before, 'long', 'before')
        api.add_qc_graph(qc_id, qc_stat_after, 'long', 'after')
        self.long_group_id = group_id
        self.long_control_id = control_id

    def export_qc_s(self):
        if self.lib_dict['small']:
            api = self.api.api('whole_transcriptome.qc')
            sample_list = os.path.join(self.tools['transfer_s'].output_dir, 'mirna_qc/clean_data/list.txt')
            group_file = os.path.join(self.tools['transfer_s'].output_dir, 'group.txt')
            compare_file = os.path.join(self.tools['transfer_s'].output_dir, 'control.txt')
            if os.path.exists(os.path.join(self.tools['transfer_s'].output_dir, 'productive_table.txt')):
                productive_table = os.path.join(self.tools['transfer_s'].output_dir, 'productive_table.txt')
            else:
                productive_table = None
            fq_type = self.small_task_info['options']['fq_type']
            qc_stat_before = os.path.join(self.tools['transfer_s'].output_dir, 'hiseq_reads_stat_raw')
            qc_stat_after = os.path.join(self.tools['transfer_s'].output_dir, 'hiseq_reads_stat_use')
            qc_result_dir = os.path.join(self.tools['transfer_s'].output_dir, 'mirna_qc/clean_data')
            api.add_sample_info(sample_list=sample_list, library='small', group_file=group_file, productive_table=productive_table)
            group_id, specimen_names, category_names = api.add_sample_group(group_file=group_file, library='small')
            control_id, compare_names = api.add_group_compare(compare_file=compare_file, library='small',
                                                              group_id=group_id)
            qc_id = api.add_qc(fq_type=fq_type, library='small')
            api.add_qc_detail(qc_id, qc_stat_before, qc_stat_after, 'small', qc_result_dir, group=group_file)
            api.add_qc_graph(qc_id, qc_stat_before, 'small', 'before')
            api.add_qc_graph(qc_id, qc_stat_after, 'small', 'after', qc_result_dir)
            self.small_group_id = group_id
            self.small_control_id = control_id

    def export_qc_c(self):
        if self.lib_dict['circle']:
            api = self.api.api('whole_transcriptome.qc')
            sample_list = os.path.join(self.tools['transfer_c'].output_dir, 'fastp_rna/fastq/list.txt')
            group_file = os.path.join(self.tools['transfer_c'].output_dir, 'group.txt')
            compare_file = os.path.join(self.tools['transfer_c'].output_dir, 'control.txt')
            if os.path.exists(os.path.join(self.tools['transfer_c'].output_dir, 'productive_table.txt')):
                productive_table = os.path.join(self.tools['transfer_c'].output_dir, 'productive_table.txt')
            else:
                productive_table = None
            fq_type = self.circle_task_info['options']['fq_type']
            qc_stat_before = os.path.join(self.tools['transfer_c'].output_dir, 'hiseq_reads_stat_raw')
            qc_stat_after = os.path.join(self.tools['transfer_c'].output_dir, 'hiseq_reads_stat_use')
            api.add_sample_info(sample_list=sample_list, library='circle', group_file=group_file, productive_table=productive_table)
            group_id, specimen_names, category_names = api.add_sample_group(group_file=group_file, library='circle')
            control_id, compare_names = api.add_group_compare(compare_file=compare_file, library='circle',
                                                              group_id=group_id)
            qc_id = api.add_qc(fq_type=fq_type, library='circle')
            api.add_qc_detail(qc_id, qc_stat_before, qc_stat_after, 'circle', group=group_file)
            api.add_qc_graph(qc_id, qc_stat_before, 'circle', 'before')
            api.add_qc_graph(qc_id, qc_stat_after, 'circle', 'after')
            self.circle_group_id = group_id
            self.circle_control_id = control_id

    def export_mapping_l(self, assess=True):
        api = self.api.api('whole_transcriptome.mapping')
        stat_file = os.path.join(self.tools['transfer_l'].output_dir, 'rnaseq_mapping/stat')
        method = self.long_task_info['options']['align_method']
        group_file = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        api.add_mapping_stat(stat_file=stat_file, library='long', method=method, group=group_file)
        if assess:
            chr_dir = os.path.join(self.tools['transfer_l'].output_dir, 'map_assessment/chr_stat')
            cov_dir = os.path.join(self.tools['transfer_l'].output_dir, 'map_assessment/coverage')
            dis_dir = os.path.join(self.tools['transfer_l'].output_dir, 'map_assessment/distribution')
            sat_dir = os.path.join(self.tools['transfer_l'].output_dir, 'map_assessment/saturation')
            params = json.dumps({'task_id': self.task_id, 'submit_location': 'mapping', 'task_type': 2}, sort_keys=True,
                                separators=(',', ':'))
            api.add_chrom_distribution_table(chr_dir, params=params, library='long', group=group_file)
            api.add_coverage_table(cov_dir, params=params, detail=True, library='long')
            api.add_distribution_table(dis_dir, params=params, library='long', group=group_file)
            api.add_rpkm_table(sat_dir, params=params, detail=True, library='long')

    def export_mapping_s(self):
        if self.lib_dict['small']:
            api = self.api.api('whole_transcriptome.mapping')
            os.path.join(self.tools['transfer_s'].output_dir, 'mapper_and_stat')
            stat_file = os.path.join(self.tools['transfer_s'].output_dir, 'mapper_and_stat')
            distribution = os.path.join(self.tools['transfer_s'].output_dir, 'mapper_and_stat')
            sample_list_file = os.path.join(self.tools['transfer_s'].output_dir, 'mirna_qc/clean_data/list.txt')
            group_file = os.path.join(self.tools['transfer_s'].output_dir, 'group.txt')
            method = 'bowtie'
            params = json.dumps({'task_id': self.task_id, 'submit_location': 'mapping', 'task_type': 2}, sort_keys=True,
                                separators=(',', ':'))
            api.add_mapping_stat(stat_file=stat_file, library='small', method=method, sample_list_file=sample_list_file, group=group_file)
            api.add_chrom_distribution_table(distribution=distribution, params=params, library='small',
                                             sample_list_file=sample_list_file, group=group_file)

    def export_mapping_c(self, assess=True):
        if self.lib_dict['circle']:
            group_file = os.path.join(self.tools['transfer_c'].output_dir, 'group.txt')
            api = self.api.api('whole_transcriptome.mapping')
            stat_file = os.path.join(self.tools['transfer_c'].output_dir, 'rnaseq_mapping/stat')
            method = 'Hisat'
            api.add_mapping_stat(stat_file=stat_file, library='circle', method=method, group=group_file)
            if assess:
                chr_dir = os.path.join(self.tools['transfer_c'].output_dir, 'map_assessment/chr_stat')
                cov_dir = os.path.join(self.tools['transfer_c'].output_dir, 'map_assessment/coverage')
                dis_dir = os.path.join(self.tools['transfer_c'].output_dir, 'map_assessment/distribution')
                sat_dir = os.path.join(self.tools['transfer_c'].output_dir, 'map_assessment/saturation')
                params = json.dumps({'task_id': self.task_id, 'submit_location': 'mapping', 'task_type': 2},
                                    sort_keys=True,
                                    separators=(',', ':'))
                api.add_chrom_distribution_table(chr_dir, params=params, library='circle', group=group_file)
                api.add_coverage_table(cov_dir, params=params, detail=True, library='circle')
                api.add_distribution_table(dis_dir, params=params, library='circle', group=group_file)
                api.add_rpkm_table(sat_dir, params=params, detail=True, library='circle')

    def export_assembly(self):
        api = self.api.api('whole_transcriptome.assembly')
        map_dict = {
            'step': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/step.pk'),
            'code': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/code.pk')
        }
        check_map_dict(map_dict)
        api.add_assembly(map_dict, self.task_id, self.project_sn)

    def export_annotation(self):
        api = self.api.api('whole_transcriptome.annotation')

        way = self.long_task_info['options']['exp_way']
        map_dict = {
            'all_t2g': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/all_tran2gene.txt'),
            'ref_t2g': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/refannot_class/all_tran2gene.txt'),
            'new_t2g': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/newannot_class/all_tran2gene.txt'),
            'T_go': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/go/go_venn_tran.txt'),
            'T_kegg': os.path.join(self.tools['transfer_l'].output_dir,
                                   'annotation/allannot_class/kegg/kegg_venn_tran.txt'),
            'T_cog': os.path.join(self.tools['transfer_l'].output_dir,
                                  'annotation/allannot_class/cog/cog_venn_tran.txt'),
            'T_nr': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/nr/nr_venn_tran.txt'),
            'T_swissprot': os.path.join(self.tools['transfer_l'].output_dir,
                                        'annotation/allannot_class/swissprot/swissprot_venn_tran.txt'),
            'T_pfam': os.path.join(self.tools['transfer_l'].output_dir,
                                   'annotation/allannot_class/pfam/pfam_venn_tran.txt'),
            'G_go': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/go/go_venn_gene.txt'),
            'G_kegg': os.path.join(self.tools['transfer_l'].output_dir,
                                   'annotation/allannot_class/kegg/kegg_venn_gene.txt'),
            'G_cog': os.path.join(self.tools['transfer_l'].output_dir,
                                  'annotation/allannot_class/cog/cog_venn_gene.txt'),
            'G_nr': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/nr/nr_venn_gene.txt'),
            'G_swissprot': os.path.join(self.tools['transfer_l'].output_dir,
                                        'annotation/allannot_class/swissprot/swissprot_venn_gene.txt'),
            'G_pfam': os.path.join(self.tools['transfer_l'].output_dir,
                                   'annotation/allannot_class/pfam/pfam_venn_gene.txt'),
            'T_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/T.reads.txt'),
            'G_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/G.reads.txt'),
            'T_exp': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/T.{}.txt'.format(way)),
            'G_exp': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/G.{}.txt'.format(way)),
        }
        check_map_dict(map_dict)
        api.add_annotation_stat(map_dict, self.task_id, self.project_sn)

        map_dict = {
            'T_new_2': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/go/go_lev2_tran.stat.xls'),
            'T_new_3': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/go/go_lev3_tran.stat.xls'),
            'T_new_4': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/go/go_lev4_tran.stat.xls'),
            'T_new_gos': os.path.join(self.tools['transfer_l'].output_dir,
                                      'annotation/newannot_class/go/go_list_tran.xls'),
            'G_new_2': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/go/go_lev2_gene.stat.xls'),
            'G_new_3': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/go/go_lev3_gene.stat.xls'),
            'G_new_4': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/go/go_lev4_gene.stat.xls'),
            'G_new_gos': os.path.join(self.tools['transfer_l'].output_dir,
                                      'annotation/newannot_class/go/go_list_gene.xls'),

            'T_ref_2': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/go/go_lev2_tran.stat.xls'),
            'T_ref_3': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/go/go_lev3_tran.stat.xls'),
            'T_ref_4': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/go/go_lev4_tran.stat.xls'),
            'T_ref_gos': os.path.join(self.tools['transfer_l'].output_dir,
                                      'annotation/refannot_class/go/go_list_tran.xls'),
            'G_ref_2': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/go/go_lev2_gene.stat.xls'),
            'G_ref_3': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/go/go_lev3_gene.stat.xls'),
            'G_ref_4': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/go/go_lev4_gene.stat.xls'),
            'G_ref_gos': os.path.join(self.tools['transfer_l'].output_dir,
                                      'annotation/refannot_class/go/go_list_gene.xls'),

            'T_all_2': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/go/go_lev2_tran.stat.xls'),
            'T_all_3': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/go/go_lev3_tran.stat.xls'),
            'T_all_4': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/go/go_lev4_tran.stat.xls'),
            'G_all_2': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/go/go_lev2_gene.stat.xls'),
            'G_all_3': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/go/go_lev3_gene.stat.xls'),
            'G_all_4': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/go/go_lev4_gene.stat.xls'),
        }
        check_map_dict(map_dict)
        api.add_annotation_go(map_dict=map_dict, task_id=self.task_id, project_sn=self.project_sn)

        map_dict = {
            'T_new_c': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_layer_tran.xls'),
            'T_new_l': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_pathway_tran.xls'),
            'T_new_p': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_pathway_tran_dir'),
            'T_new_t': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_gene_tran.xls'),
            'G_new_c': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_layer_gene.xls'),
            'G_new_l': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_pathway_gene.xls'),
            'G_new_p': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_pathway_gene_dir'),
            'G_new_t': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_class/kegg/kegg_gene_gene.xls'),

            'T_ref_c': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_layer_tran.xls'),
            'T_ref_l': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_pathway_tran.xls'),
            'T_ref_p': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_pathway_tran_dir'),
            'T_ref_t': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_gene_tran.xls'),
            'G_ref_c': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_layer_gene.xls'),
            'G_ref_l': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_pathway_gene.xls'),
            'G_ref_p': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_pathway_gene_dir'),
            'G_ref_t': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/refannot_class/kegg/kegg_gene_gene.xls'),

            'T_all_l': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/kegg/kegg_pathway_tran.xls'),
            'T_all_p': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/kegg/kegg_pathway_tran_dir'),
            'G_all_l': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/kegg/kegg_pathway_gene.xls'),
            'G_all_p': os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/allannot_class/kegg/kegg_pathway_gene_dir'),
        }
        api.add_annotation_kegg(map_dict=map_dict, task_id=self.task_id, project_sn=self.project_sn)

        map_dict = {
            'T_new': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/newannot_class/cog/summary.T.tsv'),
            'G_new': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/newannot_class/cog/summary.G.tsv'),
            'T_ref': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/refannot_class/cog/summary.T.tsv'),
            'G_ref': os.path.join(self.tools['transfer_l'].output_dir, 'annotation/refannot_class/cog/summary.G.tsv'),
        }
        check_map_dict(map_dict)
        api.add_annotation_cog(map_dict=map_dict, task_id=self.task_id, project_sn=self.project_sn)

    def export_lnc_identify(self):
        api = self.api.api('whole_transcriptome.lnc_identify')
        params = dict(
            cpc=self.long_task_info['options']['cpc'] == 'True',
            cnci=self.long_task_info['options']['cnci'] == 'True',
            cpat=self.long_task_info['options']['cpat'] == 'True',
            pfamscan=self.long_task_info['options']['pfamscan'] == 'True',
            fasta_file=str(),
            gtf_file=str(),
            identify_num=self.long_task_info['options']['identify_num'],
            transcript_len=self.long_task_info['options']['transcript_len'],
            exon_num=self.long_task_info['options']['exon_num'],
            orf_len=self.long_task_info['options']['orf_len'],
            cpc_score=self.long_task_info['options']['cpc_score'],
            taxonmy=self.long_task_info['options']['taxonmy'],
            cnci_score=self.long_task_info['options']['cnci_score'],
            cpat_score=self.long_task_info['options']['cpat_score'])

        predictions_detail_path = os.path.join(self.tools['transfer_l'].output_dir,
                                               'large_gush/known_lnc_identify/known_lncrna_detail.xls')
        api.known_lncrna_info(predictions_detail_path)

        soft_list = [lnc_soft for lnc_soft in ['cpc', 'cnci', 'cpat', 'pfamscan'] if params[lnc_soft]]

        lnc_annot_path = os.path.join(self.tools['transfer_l'].output_dir,
                                      'large_gush/lncannot/novel_lncrna_vs_lncrna.xls')
        if not os.path.isfile(lnc_annot_path):
            lnc_annot_path = None

        predictions_detail_path = os.path.join(self.tools['transfer_l'].output_dir,
                                               'large_gush/filter_by_express/filtered_lncnovel/novel_lncrna_predict_detail.xls')
        predictions_stat_path = os.path.join(self.tools['transfer_l'].output_dir,
                                             'large_gush/filter_by_express/filtered_lncnovel/novel_lncrna_stat.json')
        api.new_lncrna_predict(
            predictions_detail_path=predictions_detail_path,
            predictions_stat_path=predictions_stat_path,
            tools=','.join(soft_list).replace('pfamscan', 'pfam'),
            params=params,
            lnc_annot_path=lnc_annot_path)

        lncrna_stat_in_sample = os.path.join(self.tools['transfer_l'].output_dir,
                                             'large_gush/lncrna_stat/lncrna_stat_in_sample.xls')
        lncrna_stat_in_category = os.path.join(self.tools['transfer_l'].output_dir,
                                               'large_gush/lncrna_stat/lncrna_stat_in_category.xls')
        api.lncrna_stat(lncrna_stat_in_sample, lncrna_stat_in_category)

    def export_srna(self):
        if self.lib_dict['small']:
            api = self.api.api('whole_transcriptome.srna')
            group_table = os.path.join(self.tools['transfer_s'].output_dir, 'group.txt')
            known_mirna = os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna/known_mirna_detail.xls')
            novel_mirna = os.path.join(self.tools['transfer_s'].output_dir, 'srna/novel_mirna/novel_mirna_detail.xls')
            pdfs_known = os.path.join(self.tools['transfer_s'].output_dir, 'srna/known_mirna/structure_pdf')
            pdfs_novel = os.path.join(self.tools['transfer_s'].output_dir, 'srna/novel_mirna/structure_pdf')
            ncrna_stat = os.path.join(self.tools['transfer_s'].output_dir, 'srna/srna_stat/ncrna_stat.xls')
            srna_stat = os.path.join(self.tools['transfer_s'].output_dir, 'srna/srna_stat/srna_stat.xls')
            mirna_stat = os.path.join(self.tools['transfer_s'].output_dir, 'srna/srna_stat/mirna_stat.xls')
            srna_stat_for_graph = os.path.join(self.tools['transfer_s'].output_dir,
                                               'srna/srna_stat/srna_stat_for_graph.xls')
            if self.small_task_info['version'] >= 'v1.2':
                category = str(self.small_task_info['options']['taxonmy']).lower()
            else:
                category = str(self.small_task_info['options']['mirbase_category']).lower()

            params = json.dumps(
                {'task_id': self.task_id, 'submit_location': 'srna', 'task_type': 2, 'method': 'quantifier'},
                sort_keys=True, separators=(',', ':'))
            api.add_known_mirna(known_mirna, project_sn=self.project_sn, task_id=self.task_id, params=params,
                                pdfs=pdfs_known)
            params = json.dumps(
                {'task_id': self.task_id, 'submit_location': 'srna', 'task_type': 2, 'method': 'mirdeep'},
                sort_keys=True, separators=(',', ':'))
            api.add_novel_mirna(novel_mirna, project_sn=self.project_sn, task_id=self.task_id, params=params,
                                category=category, pdfs=pdfs_novel)
            params = json.dumps(
                {'task_id': self.task_id, 'submit_location': 'srna', 'task_type': 2, 'method': 'mirdeep'},
                sort_keys=True, separators=(',', ':'))
            api.add_mirna_stat(mirna_stat, project_sn=self.project_sn, task_id=self.task_id, params=params, group=group_table)
            params = json.dumps(
                {'task_id': self.task_id, 'submit_location': 'srna', 'task_type': 2, 'method': 'mirdeep'},
                sort_keys=True, separators=(',', ':'))
            api.add_srna_stat(srna_stat, srna_stat_for_graph, project_sn=self.project_sn, task_id=self.task_id,
                              params=params, group=group_table)

    def export_circrna(self):
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            api = self.api.api('whole_transcriptome.circrna')
            if self.lib_dict['circle']:
                taxonomy = self.circle_task_info['options']['taxonmy']
                circ_detail = os.path.join(self.tools['transfer_c'].output_dir, 'circ_brush/detail.txt')
            else:
                taxonomy = self.long_task_info['options']['taxonmy']
                circ_detail = os.path.join(self.tools['transfer_l'].output_dir, 'circ_brush/detail.txt')
            api.add_circrna(task_id=self.task_id, project_sn=self.project_sn, circ_detail=circ_detail,
                            taxonomy=taxonomy)

    def export_exp_detail(self):
        api = self.api.api('whole_transcriptome.expression')
        map_dict = {
            'mRNA_fpkm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/G.fpkm.txt'),
            'mRNA_tpm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/G.tpm.txt'),
            'lncRNA_fpkm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/lncrna/G.fpkm.txt'),
            'lncRNA_tpm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/lncrna/G.tpm.txt')
        }
        check_map_dict(map_dict)
        rna_list = list()
        if 'mRNA' in self.long_task_info['options']['rna_select']:
            rna_list.append('mRNA')
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            rna_list.append('lncRNA')
        categories = ','.join(rna_list)
        way = self.long_task_info['options']['exp_way']
        api.add_exp_g(map_dict, categories, way, self.task_id, self.project_sn)

        if self.lib_dict['circle']:
            circRNA_rpm = os.path.join(self.tools['transfer_c'].output_dir, 'exp_make/circrna/T.rpm.txt')
        else:
            circRNA_rpm = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/circrna/T.rpm.txt')
        map_dict = {
            'mRNA_fpkm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/T.fpkm.txt'),
            'mRNA_tpm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/T.tpm.txt'),
            'lncRNA_fpkm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/lncrna/T.fpkm.txt'),
            'lncRNA_tpm': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/lncrna/T.tpm.txt'),
        }
        if self.lib_dict['small']:
            map_dict['miRNA_tpm'] = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mirna/T.tpm.txt')
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            map_dict['circRNA_rpm'] = circRNA_rpm
        check_map_dict(map_dict)
        rna_list = list()
        if 'mRNA' in self.long_task_info['options']['rna_select']:
            rna_list.append('mRNA')
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            rna_list.append('lncRNA')
        if self.lib_dict['small']:
            rna_list.append('miRNA')
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            rna_list.append('circRNA')
        categories = ','.join(rna_list)
        way = self.long_task_info['options']['exp_way']
        api.add_exp_t(map_dict, categories, way, self.task_id, self.project_sn)

    def export_exp_graph(self):
        api = self.api.api('whole_transcriptome.expression_new')

        map_dict = {
            'sample_box': os.path.join(self.modules['formation_mrna'].output_dir, 'exp_graph/sample_box_data.pk'),
            'group_box': os.path.join(self.modules['formation_mrna'].output_dir, 'exp_graph/group_box_data.pk'),
            'sample_density': os.path.join(self.modules['formation_mrna'].output_dir,
                                           'exp_graph/sample_density_data.pk'),
            'group_density': os.path.join(self.modules['formation_mrna'].output_dir, 'exp_graph/group_density_data.pk'),
            'sample_volin': os.path.join(self.modules['formation_mrna'].output_dir, 'exp_graph/sample_volin_data.pk'),
            'group_volin': os.path.join(self.modules['formation_mrna'].output_dir, 'exp_graph/group_volin_data.pk'),
        }
        check_map_dict(map_dict)
        category = 'mRNA'
        exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
        level = 'T'
        kind = 'all'
        group_id = self.long_group_id
        group_dict = api.get_group_dict(group_id)
        api.add_exp_graph(map_dict, category, exp_id, level, kind, group_id, group_dict, self.task_id, self.project_sn)

        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            map_dict = {
                'sample_box': os.path.join(self.modules['formation_lncrna'].output_dir, 'exp_graph/sample_box_data.pk'),
                'group_box': os.path.join(self.modules['formation_lncrna'].output_dir, 'exp_graph/group_box_data.pk'),
                'sample_density': os.path.join(self.modules['formation_lncrna'].output_dir,
                                               'exp_graph/sample_density_data.pk'),
                'group_density': os.path.join(self.modules['formation_lncrna'].output_dir,
                                              'exp_graph/group_density_data.pk'),
                'sample_volin': os.path.join(self.modules['formation_lncrna'].output_dir,
                                             'exp_graph/sample_volin_data.pk'),
                'group_volin': os.path.join(self.modules['formation_lncrna'].output_dir,
                                            'exp_graph/group_volin_data.pk'),
            }
            check_map_dict(map_dict)
            category = 'lncRNA'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.long_group_id
            group_dict = api.get_group_dict(group_id)
            api.add_exp_graph(map_dict, category, exp_id, level, kind, group_id, group_dict, self.task_id,
                              self.project_sn)

        if self.lib_dict['small']:
            map_dict = {
                'sample_box': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_graph/sample_box_data.pk'),
                'group_box': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_graph/group_box_data.pk'),
                'sample_density': os.path.join(self.modules['formation_mirna'].output_dir,
                                               'exp_graph/sample_density_data.pk'),
                'group_density': os.path.join(self.modules['formation_mirna'].output_dir,
                                              'exp_graph/group_density_data.pk'),
                'sample_volin': os.path.join(self.modules['formation_mirna'].output_dir,
                                             'exp_graph/sample_volin_data.pk'),
                'group_volin': os.path.join(self.modules['formation_mirna'].output_dir,
                                            'exp_graph/group_volin_data.pk'),
            }
            check_map_dict(map_dict)
            category = 'miRNA'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.small_group_id
            group_dict = api.get_group_dict(group_id)
            api.add_exp_graph(map_dict, category, exp_id, level, kind, group_id, group_dict, self.task_id,
                              self.project_sn)

        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            map_dict = {
                'sample_box': os.path.join(self.modules['formation_circrna'].output_dir,
                                           'exp_graph/sample_box_data.pk'),
                'group_box': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_graph/group_box_data.pk'),
                'sample_density': os.path.join(self.modules['formation_circrna'].output_dir,
                                               'exp_graph/sample_density_data.pk'),
                'group_density': os.path.join(self.modules['formation_circrna'].output_dir,
                                              'exp_graph/group_density_data.pk'),
                'sample_volin': os.path.join(self.modules['formation_circrna'].output_dir,
                                             'exp_graph/sample_volin_data.pk'),
                'group_volin': os.path.join(self.modules['formation_circrna'].output_dir,
                                            'exp_graph/group_volin_data.pk'),
            }
            # check_map_dict(map_dict)
            category = 'circRNA'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            if self.lib_dict['circle']:
                group_id = self.circle_group_id
            else:
                group_id = self.long_group_id
            group_dict = api.get_group_dict(group_id)
            api.add_exp_graph(map_dict, category, exp_id, level, kind, group_id, group_dict, self.task_id,
                              self.project_sn)

    def export_exp_venn(self):
        api = self.api.api('whole_transcriptome.expression_new')

        map_dict = {
            'venn': os.path.join(self.modules['formation_mrna'].output_dir, 'exp_venn/venn.txt')
        }
        check_map_dict(map_dict)
        category = 'mRNA'
        exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
        level = 'T'
        kind = 'all'
        group_id = self.long_group_id
        group_dict = api.get_group_dict(group_id)
        threshold = 1.0
        api.add_exp_venn(map_dict, category, exp_id, level, kind, group_id, group_dict, threshold, self.task_id,
                         self.project_sn)

        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            map_dict = {
                'venn': os.path.join(self.modules['formation_lncrna'].output_dir, 'exp_venn/venn.txt')
            }
            check_map_dict(map_dict)
            category = 'lncRNA'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.long_group_id
            group_dict = api.get_group_dict(group_id)
            threshold = 1.0
            api.add_exp_venn(map_dict, category, exp_id, level, kind, group_id, group_dict, threshold, self.task_id,
                             self.project_sn)

        if self.lib_dict['small']:
            map_dict = {
                'venn': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_venn/venn.txt')
            }
            check_map_dict(map_dict)
            category = 'miRNA'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.small_group_id
            group_dict = api.get_group_dict(group_id)
            threshold = 1.0
            api.add_exp_venn(map_dict, category, exp_id, level, kind, group_id, group_dict, threshold, self.task_id,
                             self.project_sn)

        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            map_dict = {
                'venn': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_venn/venn.txt')
            }
            check_map_dict(map_dict)
            category = 'circRNA'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            if self.lib_dict['circle']:
                group_id = self.circle_group_id
            else:
                group_id = self.long_group_id
            group_dict = api.get_group_dict(group_id)
            threshold = 0.0
            api.add_exp_venn(map_dict, category, exp_id, level, kind, group_id, group_dict, threshold, self.task_id,
                             self.project_sn)

    def export_exp_corr(self):
        api = self.api.api('whole_transcriptome.expression_new')

        if self.lib_dict['circle']:
            map_dict = {
                'corr': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_corr/corr.txt'),
                'tree': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_corr/tree.txt')
            }
            check_map_dict(map_dict)
            library = 'circle'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.circle_group_id
            group_dict = api.get_group_dict(group_id)
            way = 'rpm'
            take_mean = 'no'
            corr_method = 'pearson'
            take_log = 'no'
            dist_method = 'euclidean'
            clus_method = 'complete'
            api.add_exp_corr(map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, corr_method,
                             take_log, dist_method, clus_method, self.task_id, self.project_sn)

        if self.lib_dict['small']:
            map_dict = {
                'corr': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_corr/corr.txt'),
                'tree': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_corr/tree.txt')
            }
            check_map_dict(map_dict)
            library = 'small'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.small_group_id
            group_dict = api.get_group_dict(group_id)
            way = 'tpm'
            take_mean = 'no'
            corr_method = 'pearson'
            take_log = 'no'
            dist_method = 'euclidean'
            clus_method = 'complete'
            api.add_exp_corr(map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, corr_method,
                             take_log, dist_method, clus_method, self.task_id, self.project_sn)

        map_dict = {
            'corr': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_corr/corr.txt'),
            'tree': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_corr/tree.txt')
        }
        check_map_dict(map_dict)
        library = 'long'
        exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
        level = 'T'
        kind = 'all'
        group_id = self.long_group_id
        group_dict = api.get_group_dict(group_id)
        way = self.long_task_info['options']['exp_way']
        take_mean = 'no'
        corr_method = 'pearson'
        take_log = 'no'
        dist_method = 'euclidean'
        clus_method = 'complete'
        api.add_exp_corr(map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, corr_method, take_log,
                         dist_method, clus_method, self.task_id, self.project_sn)

    def export_exp_pca(self):
        api = self.api.api('whole_transcriptome.expression_new')

        if self.lib_dict['circle']:
            map_dict = {
                'evr': os.path.join(self.modules['formation_circrna'].output_dir,
                                    'exp_pca/explained_variance_ratio.txt'),
                'pca': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_pca/pca.txt'),
                'ellipse': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_pca/ellipse.txt')
            }
            if not os.path.isfile(map_dict['ellipse']):
                map_dict.pop('ellipse')
            library = 'circle'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.circle_group_id
            group_dict = api.get_group_dict(group_id)
            way = 'rpm'
            take_mean = 'no'
            api.add_exp_pca(map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, self.task_id,
                            self.project_sn)

        if self.lib_dict['small']:
            map_dict = {
                'evr': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_pca/explained_variance_ratio.txt'),
                'pca': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_pca/pca.txt'),
                'ellipse': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_pca/ellipse.txt')
            }
            if not os.path.isfile(map_dict['ellipse']):
                map_dict.pop('ellipse')
            library = 'small'
            exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
            level = 'T'
            kind = 'all'
            group_id = self.small_group_id
            group_dict = api.get_group_dict(group_id)
            way = 'tpm'
            take_mean = 'no'
            api.add_exp_pca(map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, self.task_id,
                            self.project_sn)

        map_dict = {
            'evr': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_pca/explained_variance_ratio.txt'),
            'pca': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_pca/pca.txt'),
            'ellipse': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_pca/ellipse.txt')
        }
        if not os.path.isfile(map_dict['ellipse']):
            map_dict.pop('ellipse')
        library = 'long'
        exp_id = api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id']
        level = 'T'
        kind = 'all'
        group_id = self.long_group_id
        group_dict = api.get_group_dict(group_id)
        way = self.long_task_info['options']['exp_way']
        take_mean = 'no'
        api.add_exp_pca(map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, self.task_id,
                        self.project_sn)

    def export_gensets_analysis(self):
        if os.path.exists(os.path.join(self.modules["genesets_analysis"].output_dir, "results_info")):
            pass
        else:
            self.run_diff_geneset_analysis = True
            if os.path.exists(os.path.join(self.work_dir, "temporary")):
                shutil.rmtree(os.path.join(self.work_dir, "temporary"))
            os.makedirs(os.path.join(self.work_dir, "temporary"))
            self.export_temporary = os.path.join(self.work_dir, "temporary")
            api = self.api.api('whole_transcriptome.diff_geneset_work_pipline')
            diff_geneset_pipline_result = self.modules["genesets_analysis"].output_dir
            diff_id = None
            task_id = self.task_id
            analysis_names = ["kegg", "go", "cog"]
            file_json_path = os.path.join(self.modules["genesets_analysis"].file_prepare.output_dir, "prepare_json")
            with open(file_json_path, "r") as j:
                file_dict = json.load(j)
            kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
            task_infos = {"task_id":self.task_id,"project_sn":self.project_sn}
            api.add_diff_genest_pipline_table(diff_geneset_pipline_result, diff_id=diff_id, task_id=task_id,
                                              analysis_names=analysis_names,
                                              kegg_level_path=kegg_level_path, inter_path=self.export_temporary,
                                              group_id=self.long_group_id,task_infos = task_infos)



    def export_diff(self):
        api = self.api.api('whole_transcriptome.expression_new')

        map_dict = dict()
        # diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_g')
        diff_dir = os.path.join(self.tools['diff_split_g'].output_dir, 'mrna')
        for fname in os.listdir(diff_dir):
            if fname.endswith('.detail.txt'):
                map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
            elif 'summary' in fname:
                map_dict['summary'] = os.path.join(diff_dir, fname)
            elif 'volcano' in fname:
                map_dict['volcano'] = os.path.join(diff_dir, fname)
            elif 'scatter' in fname:
                map_dict['scatter'] = os.path.join(diff_dir, fname)
        check_map_dict(map_dict)
        self.logger.debug(map_dict)
        self.exp_id_G = str(api.db['exp'].find_one({'task_id': self.task_id, 'level': 'G'})['main_id'])
        self.exp_id_T = str(api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id'])
        category = 'mRNA'
        level = 'G'
        kind = 'all'
        background = 'mRNA,lncRNA'
        group_id = self.long_group_id
        group_dict = api.get_group_dict(group_id)
        control_id = self.long_control_id
        # try:
        #     diff_method = self.long_task_info['options']['long_diff_method']
        # except:
        #     diff_method = self.long_task_info['options']['diff_method']
        diff_method = self.long_task_info['options']['diff_method']
        try:
            stat_type = self.long_task_info['options']['pvalue_padjust']
            stat_cutoff = self.long_task_info['options']['diff_fdr_ci']
            fc = self.long_task_info['options']['fc']
        except:
            stat_type = "padjust"
            stat_cutoff = 0.05
            fc = 2.0
        way = self.long_task_info['options']['exp_way']
        api.add_diff(
            map_dict=map_dict,
            category=category,
            level=level,
            kind=kind,
            background=background,
            filter='none',
            threshold=0.0,
            group_id=group_id,
            group_dict=group_dict,
            control_id=control_id,
            stat_type=stat_type,
            stat_cutoff=float(stat_cutoff),
            fc=float(fc),
            diff_method=diff_method,
            correct_method='BH',
            way=way,
            task_id=self.task_id,
            project_sn=self.project_sn,
            exp_id=self.exp_id_G
        )

        map_dict = dict()
        diff_dir = os.path.join(self.tools['diff_split'].output_dir, 'mrna')
        for fname in os.listdir(diff_dir):
            if fname.endswith('.detail.txt'):
                map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
            elif 'summary' in fname:
                map_dict['summary'] = os.path.join(diff_dir, fname)
            elif 'volcano' in fname:
                map_dict['volcano'] = os.path.join(diff_dir, fname)
            elif 'scatter' in fname:
                map_dict['scatter'] = os.path.join(diff_dir, fname)
        check_map_dict(map_dict)
        self.logger.debug(map_dict)
        category = 'mRNA'
        level = 'T'
        kind = 'all'
        background = 'mRNA,lncRNA'
        group_id = self.long_group_id
        group_dict = api.get_group_dict(group_id)
        control_id = self.long_control_id
        # try:
        #     diff_method = self.long_task_info['options']['long_diff_method']
        # except:
        #     diff_method = self.long_task_info['options']['diff_method']
        way = self.long_task_info['options']['exp_way']
        api.add_diff(
            map_dict=map_dict,
            category=category,
            level=level,
            kind=kind,
            background=background,
            filter='none',
            threshold=0.0,
            group_id=group_id,
            group_dict=group_dict,
            control_id=control_id,
            stat_type=stat_type,
            stat_cutoff=float(stat_cutoff),
            fc=float(fc),
            diff_method=diff_method,
            correct_method='BH',
            way=way,
            task_id=self.task_id,
            project_sn=self.project_sn,
            exp_id=self.exp_id_T
        )

        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            map_dict = dict()
            diff_dir = os.path.join(self.tools['diff_split'].output_dir, 'lncrna')
            for fname in os.listdir(diff_dir):
                if fname.endswith('.detail.txt'):
                    map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
                elif 'summary' in fname:
                    map_dict['summary'] = os.path.join(diff_dir, fname)
                elif 'volcano' in fname:
                    map_dict['volcano'] = os.path.join(diff_dir, fname)
                elif 'scatter' in fname:
                    map_dict['scatter'] = os.path.join(diff_dir, fname)
            check_map_dict(map_dict)
            self.logger.debug(map_dict)
            category = 'lncRNA'
            level = 'T'
            kind = 'all'
            background = 'mRNA,lncRNA'
            group_id = self.long_group_id
            group_dict = api.get_group_dict(group_id)
            control_id = self.long_control_id
            # try:
            #     diff_method = self.long_task_info['options']['long_diff_method']
            # except:
            #     diff_method = self.long_task_info['options']['diff_method']
            way = self.long_task_info['options']['exp_way']
            api.add_diff(
                map_dict=map_dict,
                category=category,
                level=level,
                kind=kind,
                background=background,
                filter='none',
                threshold=0.0,
                group_id=group_id,
                group_dict=group_dict,
                control_id=control_id,
                stat_type=stat_type,
                stat_cutoff=float(stat_cutoff),
                fc=float(fc),
                diff_method=diff_method,
                correct_method='BH',
                way=way,
                task_id=self.task_id,
                project_sn=self.project_sn,
                exp_id=self.exp_id_T
            )

        if self.lib_dict['small']:
            map_dict = dict()
            diff_dir = os.path.join(self.tools['transfer_s'].output_dir, 'diff_exp')
            for fname in os.listdir(diff_dir):
                if fname.endswith('.detail.txt'):
                    map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
                elif 'summary' in fname:
                    map_dict['summary'] = os.path.join(diff_dir, fname)
                elif 'volcano' in fname:
                    map_dict['volcano'] = os.path.join(diff_dir, fname)
                elif 'scatter' in fname:
                    map_dict['scatter'] = os.path.join(diff_dir, fname)
            check_map_dict(map_dict)
            self.logger.debug(map_dict)
            category = 'miRNA'
            level = 'T'
            kind = 'all'
            background = None
            group_id = self.small_group_id
            group_dict = api.get_group_dict(group_id)
            control_id = self.small_control_id
            try:
                stat_type = self.long_task_info['options']['pvalue_padjust']
                stat_cutoff = self.long_task_info['options']['diff_fdr_ci']
                fc = self.long_task_info['options']['fc']
            except:
                stat_type = "padjust"
                stat_cutoff = 0.05
                fc = 2.0
            diff_method = self.small_task_info['options']['diff_method']
            # stat_type = self.small_task_info['options']['pvalue_padjust']
            # stat_cutoff = self.small_task_info['option']['diff_fdr_ci']
            # fc = self.small_task_info['option']['fc']

            way = 'tpm'
            api.add_diff(
                map_dict=map_dict,
                category=category,
                level=level,
                kind=kind,
                background=background,
                filter='none',
                threshold=0.0,
                group_id=group_id,
                group_dict=group_dict,
                control_id=control_id,
                stat_type=stat_type,
                stat_cutoff=float(stat_cutoff),
                fc=float(fc),
                diff_method=diff_method,
                correct_method='BH',
                way=way,
                task_id=self.task_id,
                project_sn=self.project_sn,
                exp_id=self.exp_id_T
            )

        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            map_dict = dict()
            if self.lib_dict['circle']:
                diff_dir = os.path.join(self.tools['transfer_c'].output_dir, 'diff_exp_c')
            else:
                diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_c')
            for fname in os.listdir(diff_dir):
                if fname.endswith('.detail.txt'):
                    map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
                elif 'summary' in fname:
                    map_dict['summary'] = os.path.join(diff_dir, fname)
                elif 'volcano' in fname:
                    map_dict['volcano'] = os.path.join(diff_dir, fname)
                elif 'scatter' in fname:
                    map_dict['scatter'] = os.path.join(diff_dir, fname)
            check_map_dict(map_dict)
            self.logger.debug(map_dict)
            category = 'circRNA'
            level = 'T'
            kind = 'all'
            background = None
            if self.lib_dict['circle']:
                group_id = self.circle_group_id
                control_id = self.circle_control_id
                try:
                    diff_method = self.circle_task_info['options']['diff_method']
                    stat_type = self.circle_task_info['options']['pvalue_padjust']
                    stat_cutoff = self.circle_task_info['options']['diff_fdr_ci']
                    fc = self.circle_task_info['options']['fc']
                except:
                    diff_method = 'Deseq2'
                    stat_type = 'padjust'
                    stat_cutoff = 0.05
                    fc = 2.0

            else:
                group_id = self.long_group_id
                control_id = self.long_control_id
                # diff_method = self.long_task_info['options']['diff_method']
                # stat_type = self.long_task_info['options']['pvalue_padjust']
                # stat_cutoff = self.long_task_info['option']['diff_fdr_ci']
                # fc = self.long_task_info['option']['fc']
                try:
                    diff_method = self.long_task_info['options']['diff_method_circ']
                    stat_type = self.long_task_info['options']['pvalue_padjust_circ']
                    stat_cutoff = self.long_task_info['options']['diff_fdr_ci_circ']
                    fc = self.long_task_info['options']['fc_circ']
                except:
                    diff_method = 'Deseq2'
                    stat_type = 'padjust'
                    stat_cutoff = 0.05
                    fc = 2.0


            group_dict = api.get_group_dict(group_id)
            way = 'rpm'
            api.add_diff(
                map_dict=map_dict,
                category=category,
                level=level,
                kind=kind,
                background=background,
                filter='none',
                threshold=0.0,
                group_id=group_id,
                group_dict=group_dict,
                control_id=control_id,
                stat_type=stat_type,
                stat_cutoff=float(stat_cutoff),
                fc=float(fc),
                diff_method=diff_method,
                correct_method='BH',
                way=way,
                task_id=self.task_id,
                project_sn=self.project_sn,
                exp_id=self.exp_id_T
            )

    def export_small_target(self):
        if self.lib_dict['small']:
            api = self.api.api('whole_transcriptome.small_target')
            target_file = self.modules['target_mirna'].output_dir
            params_dict = {'task_id': self.task_id, 'submit_location': 'target_detail', 'task_type': 2}
            for method in ('m_miranda', 'm_targetscan', 'm_psrobot', 'm_targetfinder', 'm_rnahybrid'):
                params_dict[method] = 'yes' if method[2:] in self.option('m_method') else 'no'
            for method in ('l_miranda', 'l_targetscan', 'l_psrobot', 'l_targetfinder', 'l_rnahybrid'):
                params_dict[method] = 'yes' if method[2:] in self.option('l_method') else 'no'
            for method in ('c_miranda', 'c_targetscan', 'c_psrobot', 'c_targetfinder', 'c_rnahybrid'):
                params_dict[method] = 'yes' if method[2:] in self.option('c_method') else 'no'
            for name in ['m_miranda_score', 'm_miranda_energy', 'm_miranda_strict', 'm_rnahybird_num',
                         'm_rnahybird_energy',
                         'm_rnahybird_pvalue', 'm_ps_robot_score', 'm_targetfinder_score', 'm_min_support',
                         'l_miranda_score', 'l_miranda_energy', 'l_miranda_strict', 'l_rnahybird_num',
                         'l_rnahybird_energy',
                         'l_rnahybird_pvalue', 'l_ps_robot_score', 'l_targetfinder_score', 'l_min_support',
                         'c_miranda_score', 'c_miranda_energy', 'c_miranda_strict', 'c_rnahybird_num',
                         'c_rnahybird_energy',
                         'c_rnahybird_pvalue', 'c_ps_robot_score', 'c_targetfinder_score', 'c_min_support']:
                params_dict[name] = self.option(name)
            self.target_mirna_arg_dict = params_dict
            target_dir = os.path.join(self.modules['target_mirna'].work_dir, 'target')
            new_seq = os.path.join(target_dir, 'new.mirna.fasta')
            known_seq = os.path.join(target_dir, 'ref.mirna.fasta')
            anno_type = 'origin'
            species_name = self.small_task_info['options']['organism_name']
            api.import_target_detail(target_file, params_dict, new_seq, known_seq, anno_type, species_name)

    def export_lnc_target(self):
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            api = self.api.api('whole_transcriptome.lnc_target')
            exp_id_G = str(api.db['exp'].find_one({'task_id': self.task_id, 'level': 'G'})['main_id'])
            exp_id_T = str(api.db['exp'].find_one({'task_id': self.task_id, 'level': 'T'})['main_id'])
            group_id = self.long_group_id
            group_dict = api.get_group_dict(group_id)
            params_dict = {
                'task_id': self.task_id,
                'submit_location': 'target_cistrans',
                'task_type': 2,
                'up_dis': self.option('up_dis'),
                'down_dis': self.option('down_dis'),
                'group_id': str(group_id),
                'group_dict': group_dict,
                'corr_cutoff': str(self.option('lt_cor_cut')),
                'corr_way': self.option('lt_cor_way'),
                'pvalue_type': self.option('lt_sig_way'),
                'pvalue_cutoff': str(self.option('lt_sig_cut')),
                'padjust_way': self.option('lt_adj_way'),
                'exp_id_G': exp_id_G,
                'exp_id_T': exp_id_T
            }
            cis_target = os.path.join(self.modules['target_cistrans'].output_dir, 'cis_annot.xls')
            trans_target = os.path.join(self.modules['target_cistrans'].output_dir, 'trans_annot.xls')
            filter_lnc_file = os.path.join(self.tools['diff_split'].output_dir, 'lncrna/summary.txt')
            diff_exp_g = os.path.join(self.tools['diff_split_g'].output_dir, 'mrna')
            # filter_target_file = os.path.join(, 'diff_exp_g/summary.txt')
            filter_target_file = os.path.join(diff_exp_g, 'summary.txt')
            api.import_target_cistrans(params_dict, cis_target, trans_target, filter_lnc_file, filter_target_file)

    def export_rmats(self):
        api = self.api.api('whole_transcriptome.rmats')
        group_id = self.long_group_id
        group_dict = api.get_group_dict(group_id)
        control_id = self.long_control_id
        compare_names = api.get_compare_names(control_id)
        self.logger.info("开始导表rmats相关")
        for compare_plan in compare_names:
            groups = compare_plan.split('|')
            my_group_dict = {group: samples for group, samples in group_dict.items() if group in groups}
            params = json.dumps({
                'task_id': self.task_id,
                'submit_location': 'splicingrmats',
                'task_type': 2,
                'group_id': str(group_id),
                'group_dict': my_group_dict,
                'control_id': str(control_id),
                'compare_plan': compare_plan
            }, sort_keys=True, separators=(',', ':'))
            self.logger.debug(params)
            outpath = os.path.join(self.modules['rmats'].output_dir, '{}_vs_{}'.format(*groups))
            api.add_splicing_rmats(params, outpath)

    def export_rmats_count(self):
        api = self.api.api('whole_transcriptome.rmats_count')
        group_file = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        outpath = self.modules['rmats'].output_dir
        api.add_rmats_count(outpath, group=group_file)

    def export_snp(self):
        api = self.api.api('whole_transcriptome.snp')
        snp_anno = self.modules['whole_snp'].output_dir
        group_file = os.path.join(self.tools['transfer_l'].output_dir, 'group.txt')
        method_type = self.option('snp_method')
        task_id = self.task_id
        project_sn = self.project_sn
        new_output = os.path.join(self.modules['whole_snp'].work_dir, 'upload')
        if not os.path.isdir(new_output):
            os.makedirs(new_output)
        params = dict(task_id=task_id, submit_location='snp', task_type=2, method_type=self.option('snp_method'))
        api.add_snp_main(snp_anno=snp_anno, group=group_file, params=params, method_type=method_type, task_id=task_id,
                         project_sn=project_sn, new_output=new_output)

    def export_lncrna_family(self):
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            api = self.api.api('whole_transcriptome.lncrna_family')
            tabular = self.modules['lncrna_family'].option('tabular').path
            api.add_lncrna_family(tabular)

    def export_mirna_precursor(self):
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            api = self.api.api('whole_transcriptome.mirna_precursor')
            tabular = self.modules['mirna_precursor'].option('tabular').path
            api.add_mirna_precursor(tabular)

    def export_atcg_bias(self):
        if self.lib_dict['small']:
            api = self.api.api('whole_transcriptome.atcg_bias')
            family_path = self.tools['atcg_bias'].output_dir
            params = json.dumps({'task_id': self.task_id, 'submit_location': 'atcg_analyse', 'task_type': 2},
                                sort_keys=True, separators=(',', ':'))
            api.run(family_path, params)

    def export_mirna_edit(self):
        if self.lib_dict['small']:
            api = self.api.api('whole_transcriptome.mirna_edit')
            result_dir = self.modules['mirna_edit'].output_dir
            params = json.dumps({'task_id': self.task_id, 'submit_location': 'mirna_edit', 'task_type': 2},
                                sort_keys=True, separators=(',', ':'))
            api.add_mirna_edit(result_dir, self.task_id, self.project_sn, params)

    def export_smallrna_family_analyse(self):
        if self.lib_dict['small']:
            api = self.api.api('whole_transcriptome.smallrna_family_analyse')
            family_path = self.tools['smallrna_family_analyse'].output_dir
            params = json.dumps({'task_id': self.task_id, 'submit_location': 'family_analyse', 'task_type': 2},
                                sort_keys=True, separators=(',', ':'))
            api.run(family_path, params)

    def export_gene_detail(self):
        api = self.api.api('whole_transcriptome.gene_detail')
        rna_types = self.modules['gene_detail'].option('rna_type')
        result_dir = self.modules['gene_detail'].output_dir
        api.add_gene_detail(rna_types, result_dir)
        api.add_seq_stat(rna_types, result_dir)


    def export_geneset(self):
        api = self.api.api('whole_transcriptome.geneset')

        map_dict = dict()
        level = 'G'
        category = 'mRNA'
        source = 'DE_mR_G_detail'
        category_df = pd.read_table(os.path.join(self.tools['transfer_l'].output_dir,
                                                 'large_gush/filter_by_express/filtered_file/gene_type.xls'),
                                    names=['seq_id', 'category'], index_col=0, usecols=[0, 2])
        mrna_diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_g_mrna')
        if not os.path.isdir(mrna_diff_dir):
            os.mkdir(mrna_diff_dir)
        # gene_diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_g')
        gene_diff_dir = os.path.join(self.tools['diff_split_g'].output_dir, 'mrna')
        for fname in os.listdir(gene_diff_dir):
            if fname.endswith('.detail.txt'):
                df = pd.read_table(os.path.join(gene_diff_dir, fname), index_col='seq_id')
                df = df.join(category_df)
                # df = df[df['category'] == category]
                df.to_csv(os.path.join(mrna_diff_dir, fname), sep='\t')
        diff_dir = mrna_diff_dir
        for fname in os.listdir(diff_dir):
            if fname.endswith('.detail.txt'):
                map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
        check_map_dict(map_dict)
        api.add_geneset_diff(
            map_dict=map_dict,
            level=level,
            category=category,
            source=source,
            task_id=self.task_id,
            project_sn=self.project_sn
        )

        map_dict = dict()
        level = 'T'
        category = 'mRNA'
        source = 'DE_mR_T_detail'
        diff_dir = os.path.join(self.tools['diff_split'].output_dir, 'mrna')
        for fname in os.listdir(diff_dir):
            if fname.endswith('.detail.txt'):
                map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
        check_map_dict(map_dict)
        api.add_geneset_diff(
            map_dict=map_dict,
            level=level,
            category=category,
            source=source,
            task_id=self.task_id,
            project_sn=self.project_sn
        )

        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            map_dict = dict()
            level = 'T'
            category = 'lncRNA'
            source = 'DE_lncR_T_detail'
            diff_dir = os.path.join(self.tools['diff_split'].output_dir, 'lncrna')
            for fname in os.listdir(diff_dir):
                if fname.endswith('.detail.txt'):
                    map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
            check_map_dict(map_dict)
            api.add_geneset_diff(
                map_dict=map_dict,
                level=level,
                category=category,
                source=source,
                task_id=self.task_id,
                project_sn=self.project_sn
            )

        if self.lib_dict['small']:
            map_dict = dict()
            level = 'T'
            category = 'miRNA'
            source = 'DE_miR_T_detail'
            diff_dir = os.path.join(self.tools['transfer_s'].output_dir, 'diff_exp')
            for fname in os.listdir(diff_dir):
                if fname.endswith('.detail.txt'):
                    map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
            check_map_dict(map_dict)
            api.add_geneset_diff(
                map_dict=map_dict,
                level=level,
                category=category,
                source=source,
                task_id=self.task_id,
                project_sn=self.project_sn
            )

        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            map_dict = dict()
            level = 'T'
            category = 'circRNA'
            source = 'DE_circR_T_detail'
            if self.lib_dict['circle']:

                diff_dir = os.path.join(self.tools['transfer_c'].output_dir, 'diff_exp_c')
            else:
                diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_c')
            for fname in os.listdir(diff_dir):
                if fname.endswith('.detail.txt'):
                    map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
            check_map_dict(map_dict)
            api.add_geneset_diff(
                map_dict=map_dict,
                level=level,
                category=category,
                source=source,
                task_id=self.task_id,
                project_sn=self.project_sn
            )

    def export_task(self):
        api = self.api.api('whole_transcriptome.task_info')
        api.add_task_info(os.path.join(self.work_dir, 'data.json'), annot_group=self.long_task_info.get('annot_group'))
        api.add_lib_rna(self.task_id, self.lib_dict)
        api.add_sub_task(self.task_id, self.option('long_task_id'), self.option('small_task_id'),
                         self.option('circ_task_id'))
        api.add_pdf_dir(self.task_id)


    def export_report_img(self):
        report_config = os.path.join(self.chart.work_dir, 'report_config.json')
        api = self.api.api('whole_transcriptome.report_model')
        s3 = self._sheet.output.split(":")[0]
        report_img_s3 = s3 + ":commonbucket/files/report_img/wholerna/" + self.task_id
        api.add_report_image(self.task_id, report_config, report_img_s3)

    def get_group_dict(self,group_table):
        sample_list = set()
        group_dict = OrderedDict()
        with open(group_table) as f:
            header_line = f.readline()
            for line in f:
                if not line.strip():
                    continue
                tmp_list = line.strip().split()
                sample_list.add(tmp_list[0])
                for g in tmp_list[1:]:
                    group_dict.setdefault(g, list())
                    group_dict[g].append(tmp_list[0])
            for g in group_dict.keys():
                group_dict[g] = sorted(list(set(group_dict[g])))
        return  group_dict, sorted(sample_list)

    def run_chart(self):
        self.chart = self.add_tool("whole_transcriptome.chart")
        chart_dict = {
            "type": "workflow"
        }

        if self.lib_dict['long']:
            long_group_dict,longRNA_seq_samples = self.get_group_dict(os.path.join(self.tools['transfer_l'].output_dir, 'group.txt'))
            chart_dict.update({
                "long_samples": longRNA_seq_samples,
                'long_group_dict': long_group_dict
            })
        if self.lib_dict['small']:
            small_group_dict,smallRNA_seq_samples = self.get_group_dict(os.path.join(self.tools['transfer_s'].output_dir, 'group.txt'))
            chart_dict.update({
                "small_samples": smallRNA_seq_samples,
                'small_group_dict': small_group_dict
            })
        if self.lib_dict['circle']:
              circ_group_dict ,circRNA_seq_samples = self.get_group_dict(os.path.join(self.tools['transfer_c'].output_dir, 'group.txt'))
              chart_dict.update({
                  "circ_samples": circRNA_seq_samples,
                  'circ_group_dict': circ_group_dict
              })
        qc_dict ={}
        if self.lib_dict['long']:
            qc_dict.update({
                    "long":{
                        "qc_file_raw":"{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                table_dir=os.path.join(self.tools["transfer_l"].output_dir,"hiseq_reads_stat_raw"), sample_name='{sample_name}'),
                        "qc_file_use" :"{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                table_dir=os.path.join(self.tools["transfer_l"].output_dir,"hiseq_reads_stat_use"), sample_name='{sample_name}')
                    }
            })
        if self.lib_dict['small']:
            qc_dict.update({
                "small": {
                    "qc_file_raw": "{table_dir}/qualityStat/{sample_name}.qual_stat".format(
                        table_dir=os.path.join(self.tools["transfer_s"].output_dir, "hiseq_reads_stat_raw"),
                        sample_name='{sample_name}'),
                    "qc_file_use": "{table_dir}/{sample_name}_clean.length.txt".format(
                        table_dir=os.path.join(self.tools["transfer_s"].output_dir, "mirna_qc", "clean_data"),
                        sample_name='{sample_name}')
                }
            })
        if self.lib_dict['circle']:
            qc_dict.update({
                "circle": {
                    "qc_file_raw": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                        table_dir=os.path.join(self.tools["transfer_c"].output_dir, "hiseq_reads_stat_raw"),
                        sample_name='{sample_name}'),
                    "qc_file_use": "{table_dir}/qualityStat/{sample_name}.l.qual_stat,{table_dir}/qualityStat/{sample_name}.r.qual_stat".format(
                        table_dir=os.path.join(self.tools["transfer_c"].output_dir, "hiseq_reads_stat_use"),
                        sample_name='{sample_name}')
                }
            })
        chart_dict.update({
            "qc" : qc_dict
        })
        map_assess_dict = {
                    "long": {
                        'map_saturation': "{table_dir}/saturation/satur_{sample_name}.eRPKM.xls".format(
                            table_dir=os.path.join(self.tools["transfer_l"].output_dir, "map_assessment"),
                            sample_name='{long_samples}'),
                        'map_coverage': '{table_dir}/coverage/{sample_name}.geneBodyCoverage.txt'.format(
                            table_dir=os.path.join(self.tools["transfer_l"].output_dir, "map_assessment"),
                            sample_name='{long_samples}'),
                        "map_region_stat":'{table_dir}/distribution/{sample_name}.reads_distribution.txt'.format(
                            table_dir=os.path.join(self.tools["transfer_l"].output_dir, "map_assessment"),
                            sample_name='{long_samples}'),
                        "chr_reads_stat": '{table_dir}/chr_stat/{sample_name}.bam_chr_stat.xls'.format(
                            table_dir=os.path.join(self.tools["transfer_l"].output_dir, "map_assessment"),
                            sample_name='{long_samples}')
                    }
                }
        if self.lib_dict['small']:
            map_assess_dict.update({
                "small":  os.path.join(self.tools['transfer_s'].output_dir, 'mapper_and_stat')
            })
        if self.lib_dict['circle']:
            map_assess_dict.update({
                "circle": {
                    'map_saturation': "{table_dir}/saturation/satur_{sample_name}.eRPKM.xls".format(
                        table_dir=os.path.join(self.tools["transfer_c"].output_dir, "map_assessment"),
                        sample_name='{long_samples}'),
                    'map_coverage': '{table_dir}/coverage/{sample_name}.geneBodyCoverage.txt'.format(
                        table_dir=os.path.join(self.tools["transfer_c"].output_dir, "map_assessment"),
                        sample_name='{long_samples}'),
                    "map_region_stat": '{table_dir}/distribution/{sample_name}.reads_distribution.txt'.format(
                        table_dir=os.path.join(self.tools["transfer_c"].output_dir, "map_assessment"),
                        sample_name='{long_samples}'),
                    "chr_reads_stat": '{table_dir}/chr_stat/{sample_name}.bam_chr_stat.xls'.format(
                        table_dir=os.path.join(self.tools["transfer_c"].output_dir, "map_assessment"),
                        sample_name='{long_samples}')
                }
            })


        chart_dict.update({
                "map_assess":map_assess_dict
        })

        assemble_dict = {}
        assemble_dict.update({
            "assemble_step":os.path.join(self.tools['transfer_l'].output_dir, 'assembly/step.pk'),
            "assemble_code": os.path.join(self.tools['transfer_l'].output_dir, 'assembly/code.pk'),
        })
        chart_dict.update({
            "assemble":assemble_dict
        })
        express_dict = {}
        express_venn = {}
        express_corr_dict = {}
        express_pca_dict = {}

        #express_dis &express_venn
        if 'mRNA' in self.long_task_info['options']['rna_select']:
                express_dict.update({
                        "mRNA":{
                            "gene_exp":os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/G.tpm.txt'),
                            "trans_exp":os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/T.tpm.txt')
                        }
                })
                express_venn.update({
                        "mRNA": os.path.join(self.modules['formation_mrna'].output_dir, 'exp_venn/venn.txt')
                })
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
                express_dict.update({
                        "lncRNA":{
                            "gene_exp":os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/lncrna/G.tpm.txt'),
                            "trans_exp":os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/lncrna/T.tpm.txt')
                        }
                })
                express_venn.update({
                        "lncRNA": os.path.join(self.modules['formation_lncrna'].output_dir, 'exp_venn/venn.txt')
                })
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
                if self.lib_dict['circle']:
                    circRNA_rpm = os.path.join(self.tools['transfer_c'].output_dir, 'exp_make/circrna/T.rpm.txt')
                else:
                    circRNA_rpm = os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/circrna/T.rpm.txt')
                express_dict.update({
                        "circRNA":{
                            "trans_exp":circRNA_rpm
                        }
                })
                express_venn.update({
                        "circRNA": os.path.join(self.modules['formation_circrna'].output_dir, 'exp_venn/venn.txt')
                })
        if self.lib_dict['small']:
                express_dict.update({
                        "smallRNA": {
                            "trans_exp": os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mirna/T.tpm.txt')
                        }
                })
                express_venn.update({
                        "smallRNA": os.path.join(self.modules['formation_mirna'].output_dir, 'exp_venn/venn.txt')
                })
        chart_dict.update({
            "express_dis":express_dict,
            "express_venn":express_venn
        })

        #exp_corr
        map_dict = {
            'corr': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_corr/corr.txt'),
            'tree': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_corr/tree.txt')
        }
        express_corr_dict.update({
            "long":map_dict
        })
        if self.lib_dict['circle']:
            circle_map_dict = {
                'corr': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_corr/corr.txt'),
                'tree': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_corr/tree.txt')
            }
            express_corr_dict.update({
                "circle": circle_map_dict
            })
        if self.lib_dict['small']:
            small_map_dict = {
                'corr': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_corr/corr.txt'),
                'tree': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_corr/tree.txt')
            }
            express_corr_dict.update({
                "small": small_map_dict
            })
        chart_dict.update({
            "express_corr" :express_corr_dict
        })

        # exp_pca
        map_dict = {
            'evr': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_pca/explained_variance_ratio.txt'),
            'pca': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_pca/pca.txt'),
            'ellipse': os.path.join(self.modules['formation_longrna'].output_dir, 'exp_pca/ellipse.txt')
        }
        if not os.path.isfile(map_dict['ellipse']):
            map_dict.pop('ellipse')
        express_pca_dict.update({
            "long":map_dict
        })
        if self.lib_dict['small']:
            small_map_dict = {
                'evr': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_pca/explained_variance_ratio.txt'),
                'pca': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_pca/pca.txt'),
                'ellipse': os.path.join(self.modules['formation_mirna'].output_dir, 'exp_pca/ellipse.txt')
            }
            if not os.path.isfile(small_map_dict['ellipse']):
                small_map_dict.pop('ellipse')
            express_pca_dict.update({
                "small" :small_map_dict
            })
        if self.lib_dict['circle']:
            circle_map_dict = {
                'evr': os.path.join(self.modules['formation_circrna'].output_dir,
                                    'exp_pca/explained_variance_ratio.txt'),
                'pca': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_pca/pca.txt'),
                'ellipse': os.path.join(self.modules['formation_circrna'].output_dir, 'exp_pca/ellipse.txt')
            }
            if not os.path.isfile(circle_map_dict['ellipse']):
                circle_map_dict.pop('ellipse')
            express_pca_dict.update({
                "circle": circle_map_dict
            })
        chart_dict.update({
            "express_pca" : express_pca_dict
        })

        #Diff_exp
        diff_exp_dict ={}
        #mRNA 这个确定存在
        g_mrna_diff_dir = os.path.join(self.tools['diff_split_g'].output_dir, 'mrna')
        g_map_dict = dict()
        g_map_dict['summary'] = os.path.join(g_mrna_diff_dir, "summary.txt")
        g_map_dict['volcano'] = os.path.join(g_mrna_diff_dir, "volcano.txt")
        g_map_dict['scatter'] = os.path.join(g_mrna_diff_dir, "scatter.txt")
        t_mrna_diff_dir = os.path.join(self.tools['diff_split'].output_dir, 'mrna')
        t_map_dict = dict()
        t_map_dict['summary'] = os.path.join(t_mrna_diff_dir, "summary.txt")
        t_map_dict['volcano'] = os.path.join(t_mrna_diff_dir, "volcano.txt")
        t_map_dict['scatter'] = os.path.join(t_mrna_diff_dir, "scatter.txt")
        diff_method = self.long_task_info['options']['diff_method']
        try:
            stat_type = self.long_task_info['options']['pvalue_padjust']
            stat_cutoff = self.long_task_info['options']['diff_fdr_ci']
            fc = self.long_task_info['options']['fc']
        except:
            stat_type = "padjust"
            stat_cutoff = 0.05
            fc = 2.0
        long_diff_dict={}
        long_diff_dict.update({
                "g_mrna_dict" :g_map_dict,
                "t_mrna_dict" : t_map_dict,
                "options":{
                    "diff_method":diff_method,
                    "stat_type":stat_type,
                    "stat_cutoff":stat_cutoff,
                    "fc":fc
                }
        })
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            t_lncrna_diff_dir = os.path.join(self.tools['diff_split'].output_dir, 'lncrna')
            t_map_dict = dict()
            t_map_dict['summary'] = os.path.join(t_lncrna_diff_dir, "summary.txt")
            t_map_dict['volcano'] = os.path.join(t_lncrna_diff_dir, "volcano.txt")
            t_map_dict['scatter'] = os.path.join(t_lncrna_diff_dir, "scatter.txt")
            long_diff_dict.update({
                "t_lncrna_dict": t_map_dict,
            })
        chart_dict.update({
            "diff_long" :long_diff_dict
        })
        if self.lib_dict['small']:
            t_map_dict = dict()
            t_small_diff_dir = os.path.join(self.tools['transfer_s'].output_dir, 'diff_exp')
            t_map_dict['summary'] = os.path.join(t_small_diff_dir, "summary.txt")
            t_map_dict['volcano'] = os.path.join(t_small_diff_dir, "volcano.txt")
            t_map_dict['scatter'] = os.path.join(t_small_diff_dir, "scatter.txt")
            diff_method = self.small_task_info['options']['diff_method']
            try:
                stat_type = self.long_task_info['options']['pvalue_padjust']
                stat_cutoff = self.long_task_info['options']['diff_fdr_ci']
                fc = self.long_task_info['options']['fc']
            except:
                stat_type = "padjust"
                stat_cutoff = 0.05
                fc = 2.0
            small_diff_dict = {}
            small_diff_dict.update({
                "t_smallrna_dict": t_map_dict,
                "options": {
                    "diff_method": diff_method,
                    "stat_type": stat_type,
                    "stat_cutoff": stat_cutoff,
                    "fc": fc
                }
            })
            chart_dict.update({
                "diff_small": small_diff_dict
            })
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            t_map_dict = dict()
            if self.lib_dict['circle']:
                diff_dir = os.path.join(self.tools['transfer_c'].output_dir, 'diff_exp_c')
            else:
                diff_dir = os.path.join(self.tools['transfer_l'].output_dir, 'diff_exp_c')
            t_map_dict['summary'] = os.path.join(diff_dir, "summary.txt")
            t_map_dict['volcano'] = os.path.join(diff_dir, "volcano.txt")
            t_map_dict['scatter'] = os.path.join(diff_dir, "scatter.txt")
            try:
                diff_method = self.long_task_info['options']['diff_method_circ']
                stat_type = self.long_task_info['options']['pvalue_padjust_circ']
                stat_cutoff = self.long_task_info['options']['diff_fdr_ci_circ']
                fc = self.long_task_info['options']['fc_circ']
            except:
                diff_method = 'Deseq2'
                stat_type = 'padjust'
                stat_cutoff = 0.05
                fc = 2.0
            circ_diff_dict = {}
            circ_diff_dict.update({
                "t_circlerna_dict": t_map_dict,
                "options": {
                    "diff_method": diff_method,
                    "stat_type": stat_type,
                    "stat_cutoff": stat_cutoff,
                    "fc": fc
                }
            })
            chart_dict.update({
                "diff_circle": circ_diff_dict
            })
        way = self.long_task_info['options']['exp_way']

        chart_dict.update({
            "all_annot_stat": "{table_dir}/allannot_class/all_stat.xls".format(
                table_dir=os.path.join(self.tools['transfer_l'].output_dir, 'annotation'))
            # 'annot_stat_T_exp': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/T.{}.txt'.format(way)),
            # 'annot_stat_G_exp': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/mrna/G.{}.txt'.format(way))
        })
        if os.path.exists( "{table_dir}/refannot_class/all_stat.xls".format(
                table_dir=os.path.join(self.tools['transfer_l'].output_dir, 'annotation'))):
            chart_dict.update({
                "ref_annot_stat":"{table_dir}/refannot_class/all_stat.xls".format(
                table_dir=os.path.join(self.tools['transfer_l'].output_dir, 'annotation'))
            })
        if os.path.exists( "{table_dir}/refannot_class/all_stat.xls".format(
                table_dir=os.path.join(self.tools['transfer_l'].output_dir, 'annotation'))):
            chart_dict.update({
                "new_annot_stat":"{table_dir}/newannot_class/all_stat.xls".format(
                table_dir=os.path.join(self.tools['transfer_l'].output_dir, 'annotation'))
            })


        #lncRNA_predict
        lncRNA_predict_dict = {}
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            predictions_stat_path = os.path.join(self.tools['transfer_l'].output_dir,
                                                 'large_gush/filter_by_express/filtered_lncnovel/novel_lncrna_stat.json')
            lncRNA_predict_dict.update({"predictions_stat_path":predictions_stat_path})
            lncrna_stat_in_sample = os.path.join(self.tools['transfer_l'].output_dir,
                                                 'large_gush/lncrna_stat/lncrna_stat_in_sample.xls')
            lncrna_stat_in_category = os.path.join(self.tools['transfer_l'].output_dir,
                                                   'large_gush/lncrna_stat/lncrna_stat_in_category.xls')

            lncRNA_predict_dict.update({
                "lncrna_stat_in_sample": lncrna_stat_in_sample,
                "lncrna_stat_in_category" : lncrna_stat_in_category
            })
            chart_dict.update({
                "lncRNA_predict_dict" : lncRNA_predict_dict
            })

        # miRNA_predict
        miRNA_predict_dict ={}
        if self.lib_dict['small']:
            mirna_stat = os.path.join(self.tools['transfer_s'].output_dir, 'srna/srna_stat/mirna_stat.xls')
            miRNA_predict_dict.update({
                "mirna_stat" :mirna_stat
            })
            miRNA_predict_dict.update({
                "srna_stat": os.path.join(self.tools['transfer_s'].output_dir, 'srna/srna_stat/srna_stat_for_graph.xls')
            })
            chart_dict.update({
                "miRNA_predict_dict" : miRNA_predict_dict
            })

        #circRNA_predict
        circRNA_predict_dict = {}
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            if self.lib_dict['circle']:
                circ_detail = os.path.join(self.tools['transfer_c'].output_dir, 'circ_brush/detail.txt')
            else:
                circ_detail = os.path.join(self.tools['transfer_l'].output_dir, 'circ_brush/detail.txt')
            circRNA_predict_dict.update({
                "circ_detail" :circ_detail
            })
            chart_dict.update({
                "circRNA_predict_dict": circRNA_predict_dict
            })

        if os.path.exists(os.path.join(self.modules["genesets_analysis"].output_dir, "results_info")):
            self.run_diff_geneset_analysis = False
            pass
        else:
            self.run_diff_geneset_analysis = True
        #diff_geneset_analysis
        if self.run_diff_geneset_analysis:
            cluster_geneset_name = os.listdir(os.path.join(self.modules['genesets_analysis'].output_dir, "cluster"))[0]
            cluster_dir = os.path.join(self.modules['genesets_analysis'].output_dir, "cluster", cluster_geneset_name)
            geneset_analysis_dict = {}
            genesets = sorted(
                [os.path.basename(i) for i in glob.glob("{}/*vs*".format(self.modules["genesets_analysis"].output_dir))])
            geneset_analysis_dict.update({"genesets": genesets})
            all_geneset_file_path  = os.path.join(self.modules['genesets_analysis'].file_prepare.output_dir)
            geneset_analysis_dict.update({
                "diff_geneset_venn":all_geneset_file_path,
                "cluster_geneset_name": cluster_geneset_name,
                "cluster_exp": os.path.join(cluster_dir, "expression_matrix.xls"),
                "cluster_tree": os.path.join(cluster_dir, "seq.cluster_tree.txt"),
                "group_dict" :long_group_dict ,
                "sample_tree": os.path.join(cluster_dir, "sample.cluster_tree.txt"),
                "subcluster_list": glob.glob(cluster_dir + "/*subcluster_*.xls"),
                "gene_annot_file": os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/all_annot.xls')
            })
            file_json_path = os.path.join(self.modules["genesets_analysis"].file_prepare.output_dir, "prepare_json")
            with open(file_json_path, "r") as j:
                file_dict = json.load(j)
            kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]
            geneset_analysis_dict.update({
                "cog_class" : "{table_dir}/{geneset_name}/diff_cog_class/cog_class_table.xls".format(
                        table_dir=self.modules["genesets_analysis"].output_dir, geneset_name="{geneset_name}"),
                "go_class": "{table_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(
                        table_dir=self.modules["genesets_analysis"].output_dir, geneset_name="{geneset_name}"),
                "go_enrich": "{table_dir}/{geneset_name}/diff_go_enrich/go_enrich_geneset_list_gene.xls".format(
                    table_dir=self.modules["genesets_analysis"].output_dir, geneset_name="{geneset_name}"),
                "kegg_class": "{table_dir}/{geneset_name}/diff_kegg_class/kegg_stat.xls".format(
                    table_dir=self.modules["genesets_analysis"].output_dir, geneset_name="{geneset_name}"),
                "kegg_enrich": "{table_dir}/{geneset_name}/diff_kegg_enrich/enrich/{geneset_name}_gene.list.DE.list.check.kegg_enrichment.xls".format(
                    table_dir=self.modules["genesets_analysis"].output_dir, geneset_name="{geneset_name}"),
                "kegg_level": kegg_level_path,
            })
            chart_dict.update({"diff_geneset_analysis":geneset_analysis_dict})

        # rmats
        if self.option('is_as'):
            chart_dict.update({
                "splice_stat": "{table_dir}/sample.event.count.{splicestat}.txt".format(table_dir=self.modules['rmats'].output_dir,
                                                                                        splicestat='{splicestat}'),
                "splice_diff": "{table_dir}/{control}_vs_{test}/event_stats.file.txt".format(
                    table_dir=self.modules['rmats'].output_dir, control='{control}', test='{test}'),
                "splice_psi": "{table_dir}/{control}_vs_{test}/psi_stats.file.txt".format(
                    table_dir=self.modules['rmats'].output_dir, control='{control}', test='{test}')
            })
            cmp_list = self.modules['rmats'].option("control_table").prop["cmp_list"]
            chart_dict.update({
                    "long_cmp_list": cmp_list
            })

        if self.option('is_snp'):
            snp_result = os.path.join(self.modules['whole_snp'].output_dir)
            chart_dict.update({
                "snp_distribution": "{table_dir}/snp_position_distribution.xls".format(
                    table_dir=snp_result),
                "indel_distribution": "{table_dir}/indel_position_distribution.xls".format(
                    table_dir=snp_result),
                "snp_stat": "{table_dir}/snp_transition_tranversion_statistics.xls".format(
                    table_dir=snp_result),
                "snp_depth": "{table_dir}/snp_depth_statistics.xls".format(table_dir=snp_result)
            })

        if self.lib_dict['small']:
            miRNA_struction_dict ={}
            family_path = self.tools['atcg_bias'].output_dir
            miRNA_struction_dict.update({
                "all_first_bias_per":os.path.join(family_path,"all_first_bias_per.xls"),
                "known_first_bias_per":os.path.join(family_path,"known_first_bias_per.xls"),
                "novel_first_bias_per": os.path.join(family_path, "novel_first_bias_per.xls"),
                "all_loc_bias_per":os.path.join(family_path,"all_loc_bias_per.xls"),
                "known_loc_bias_per":os.path.join(family_path,"known_loc_bias_per.xls"),
                "novel_loc_bias_per":os.path.join(family_path,"novel_loc_bias_per.xls"),
            })
            miRNA_edit_dir = self.modules['mirna_edit'].output_dir
            miRNA_struction_dict.update({
                "miRNA_edit_dir" : self.modules['mirna_edit'].output_dir
            })
            chart_dict.update({
                "miRNA_struction":miRNA_struction_dict
            })


        with open(self.work_dir + "/chart_workflow.json", 'w') as json_f:
            json.dump(chart_dict, json_f, sort_keys=True, indent=4)
        self.chart.set_options({
            "file_json": self.work_dir + "/chart_workflow.json"
        })
        # self.chart.on('end', self.set_db)
        self.chart.on('end', self.set_db)
        self.chart.run()





    def set_upload(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        self.set_upload_mrna()
        self.set_upload_lncrna()
        self.set_upload_mirna()
        self.set_upload_circrna()
        self.set_upload_other()
        self.move_chart_file()
        self.def_upload_dir()
        self.end()

    def set_upload_mrna(self):
        from mbio.packages.whole_transcriptome.catalogue import mrna
        mrna.database = Config().get_mongo_client(mtype='whole_transcriptome')[
            Config().get_mongo_dbname('whole_transcriptome')]
        mrna_dir = os.path.join(self.output_dir, 'mrna')
        if os.path.isdir(mrna_dir):
            shutil.rmtree(mrna_dir)
        os.mkdir(mrna_dir)
        mrna.set_background(self.task_id, os.path.join(mrna_dir, '01_Background'))
        map_dict = {
            'new_gtf': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/new.gtf'),
            'all_gtf': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/all.gtf'),
            'new_fasta': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/new.fasta'),
            'all_fasta': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/all.fasta'),
            "ref_cds" : os.path.join(self.db_path, self.genome_doc['cds']),
            "ref_pep": os.path.join(self.db_path, self.genome_doc['pep']),
            "new_cds": os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_orfpfam/novel_mrna.fa.transdecoder.cds'),
            "new_pep": os.path.join(self.tools['transfer_l'].output_dir,
                                    'annotation/newannot_orfpfam/novel_mrna.fa.transdecoder.pep'),
            "all_id": os.path.join(self.tools['transfer_l'].output_dir, 'annotation/allannot_class/all_tran2gene.txt'),
        }
        mrna.set_basic_analysis(map_dict, self.task_id, os.path.join(mrna_dir, '02_Basic_Analysis'))
        mrna.set_annotation(self.task_id, os.path.join(mrna_dir, '03_Annotation'))
        map_dict = {
            'g_anno': os.path.join(mrna_dir, '03_Annotation/01_Anno_Detail/all_gene_anno_detail.xls'),
            't_anno': os.path.join(mrna_dir, '03_Annotation/01_Anno_Detail/all_transcript_anno_detail.xls'),
            'g_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/G.reads.txt'),
            't_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/T.reads.txt')}
        mrna.set_express(map_dict, self.task_id, os.path.join(mrna_dir, '04_Express'))
        way = self.long_task_info['options']['exp_way']
        map_dict = {
            't_count': os.path.join(mrna_dir, '04_Express/01_Exp_Annalysis/transcript_count_anno.xls'),
            't_exp': os.path.join(mrna_dir, '04_Express/01_Exp_Annalysis/transcript_{}_anno.xls'.format(way)),
            't_anno': os.path.join(mrna_dir, '03_Annotation/01_Anno_Detail/all_transcript_anno_detail.xls'),
            'g_count': os.path.join(mrna_dir, '04_Express/01_Exp_Annalysis/gene_count_anno.xls'),
            'g_exp': os.path.join(mrna_dir, '04_Express/01_Exp_Annalysis/gene_{}_anno.xls'.format(way)),
            'g_anno': os.path.join(mrna_dir, '03_Annotation/01_Anno_Detail/all_gene_anno_detail.xls')}
        mrna.set_diff_express(map_dict, self.task_id, os.path.join(mrna_dir, '05_Diff_Express'))
        if self.option('is_as'):
            map_dict = {'input_dir': self.modules['rmats'].output_dir,
                        'g_anno': os.path.join(mrna_dir, '03_Annotation/01_Anno_Detail/all_gene_anno_detail.xls')}
            mrna.set_as(map_dict, self.task_id, os.path.join(mrna_dir, '06_AS'))
        if self.option('is_snp'):
            map_dict = {'upload': os.path.join(self.modules['whole_snp'].work_dir, 'upload')}
            mrna.set_snp_indel(map_dict, os.path.join(mrna_dir, '07_SNP_InDel'))
        if os.path.exists(os.path.join(self.modules["genesets_analysis"].output_dir, "results_info")):
            pass
        else:
            mrna.set_genesets_analysis(self.modules["genesets_analysis"].output_dir,os.path.join(mrna_dir, '08_Diff_geneset_analysis'))


        # shutil.copytree(self.modules["genesets_analysis"].output_dir,os.path.join(mrna_dir,"diff_geneset_analysis"))


    def set_upload_lncrna(self):
        if 'lncRNA' in self.long_task_info['options']['rna_select']:
            from mbio.packages.whole_transcriptome.catalogue import lncrna
            lncrna.database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]
            lncrna_dir = os.path.join(self.output_dir, 'lncrna')
            if os.path.isdir(lncrna_dir):
                shutil.rmtree(lncrna_dir)
            os.mkdir(lncrna_dir)
            map_dict = {'large_gush_dir': os.path.join(self.tools['transfer_l'].output_dir, 'large_gush')}

            def add_lnc_ref( genome_id):
                database = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
                collection = database['sg_genome_db']
                genome_doc = collection.find_one({'genome_id': genome_id})
                lnc_ref = True if 'lnc_dir' in genome_doc else False
                if 'lncRNA' not in self.long_task_info['options']['rna_select']:
                    lnc_ref = False
                return lnc_ref
            lnc_ref = add_lnc_ref(self.long_task_info['options']['genome_id'])
            lncrna.set_lncrna_analysis(map_dict, self.task_id, os.path.join(lncrna_dir, '01_lncRNA_Analysis'),lnc_ref=lnc_ref)
            map_dict = {'t_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/T.reads.txt')}
            lncrna.set_express(map_dict, self.task_id, os.path.join(lncrna_dir, '02_Express'))
            way = self.long_task_info['options']['exp_way']
            map_dict = {'t_count': os.path.join(lncrna_dir, '02_Express/01_Exp_Annalysis/lncRNA_count.xls'),
                        't_exp': os.path.join(lncrna_dir, '02_Express/01_Exp_Annalysis/lncRNA_{}.xls'.format(way))}
            lncrna.set_diff_express(map_dict, self.task_id, os.path.join(lncrna_dir, '03_Diff_Express'))
            map_dict = {'cistrans_dir': self.modules['target_cistrans'].output_dir}
            lncrna.set_lncrna_target(map_dict, self.task_id, os.path.join(lncrna_dir, '04_LncRNA_Target'))
            lncrna.set_lncrna_structure(self.task_id, os.path.join(lncrna_dir, '05_LncRNA_Structure'))

    def set_upload_mirna(self):
        if self.lib_dict['small']:
            from mbio.packages.whole_transcriptome.catalogue import mirna
            mirna.database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]
            mirna_dir = os.path.join(self.output_dir, 'mirna')
            if os.path.isdir(mirna_dir):
                shutil.rmtree(mirna_dir)
            os.mkdir(mirna_dir)
            mirna.set_background(self.task_id, os.path.join(mirna_dir, '01_Background'))
            mirna.set_basic_analysis(self.task_id, os.path.join(mirna_dir, '02_Basic_Analysis'))
            map_dict = {'ref_structure_dir': os.path.join(self.tools['transfer_s'].output_dir,
                                                          'srna/known_mirna/structure_pdf'),
                        'ref_detail_table': os.path.join(self.tools['transfer_s'].output_dir,
                                                         'srna/known_mirna/known_mirna_detail.xls'),
                        'ref_hairpin_fasta': os.path.join(self.tools['transfer_s'].output_dir,
                                                          'srna/known_mirna/hairpin.fa'),
                        'ref_mature_fasta': os.path.join(self.tools['transfer_s'].output_dir,
                                                         'srna/known_mirna/mature.fa'),
                        'new_structure_dir': os.path.join(self.tools['transfer_s'].output_dir,
                                                          'srna/novel_mirna/structure_pdf'),
                        'new_detail_table': os.path.join(self.tools['transfer_s'].output_dir,
                                                         'srna/novel_mirna/novel_mirna_detail.xls'),
                        'new_hairpin_fasta': os.path.join(self.tools['transfer_s'].output_dir,
                                                          'srna/novel_mirna/novel_precursor_seq.fa'),
                        'new_mature_fasta': os.path.join(self.tools['transfer_s'].output_dir,
                                                         'srna/novel_mirna/novel_mature_seq.fa'),
                        'mirna_stat_table': os.path.join(self.tools['transfer_s'].output_dir,
                                                         'srna/srna_stat/mirna_stat.xls'),
                        'srna_stat_table': os.path.join(self.tools['transfer_s'].output_dir,
                                                        'srna/srna_stat/srna_stat.xls')}
            mirna.set_srna_analysis(map_dict, os.path.join(mirna_dir, '03_sRNA_Analysis'))
            map_dict = {'s_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/S.reads.txt')}
            mirna.set_express(map_dict, self.task_id, os.path.join(mirna_dir, '04_Express'))
            map_dict = {'s_count': os.path.join(mirna_dir, '04_Express/01_Exp_Annalysis/miRNA_count.xls'),
                        's_exp': os.path.join(mirna_dir, '04_Express/01_Exp_Annalysis/miRNA_tpm.xls')}
            mirna.set_diff_express(map_dict, self.task_id, os.path.join(mirna_dir, '05_Diff_Express'))
            map_dict = {'target_dir': self.modules['target_mirna'].output_dir}
            arg_dict = self.target_mirna_arg_dict
            mirna.set_mirna_target(map_dict, arg_dict, os.path.join(mirna_dir, '06_miRNA_Target'))
            mirna.set_mirna_structure(self.task_id, os.path.join(mirna_dir, '07_miRNA_Structure'))

    def set_upload_circrna(self):
        if 'circRNA' in self.long_task_info['options']['rna_select'] or self.lib_dict['circle']:
            from mbio.packages.whole_transcriptome.catalogue import circrna
            circrna.database = Config().get_mongo_client(mtype='whole_transcriptome')[
                Config().get_mongo_dbname('whole_transcriptome')]
            circrna_dir = os.path.join(self.output_dir, 'circrna')
            if os.path.isdir(circrna_dir):
                shutil.rmtree(circrna_dir)
            os.mkdir(circrna_dir)
            if self.lib_dict['circle']:
                circrna.set_background(self.task_id, os.path.join(circrna_dir, '01_Background'))
                circrna.set_basic_analysis(self.task_id, os.path.join(circrna_dir, '02_Basic_Analysis'))
            map_dict = {'fasta': os.path.join(self.tools['transfer_l'].output_dir, 'circ_brush/circrna.fasta')}
            circrna.set_circrna_predict(map_dict, self.task_id, os.path.join(circrna_dir, '03_circRNA_Predict'))
            map_dict = {'c_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/C.reads.txt')}
            circrna.set_express(map_dict, self.task_id, os.path.join(circrna_dir, '04_Express'))
            map_dict = {'c_count': os.path.join(circrna_dir, '04_Express/01_Exp_Annalysis/circRNA_count.xls'),
                        'c_exp': os.path.join(circrna_dir, '04_Express/01_Exp_Annalysis/circRNA_rpm.xls')}
            circrna.set_diff_express(map_dict, self.task_id, os.path.join(circrna_dir, '05_Diff_Express'))

    def set_upload_other(self):
        my_dir = os.path.join(self.output_dir, 'other')
        if os.path.isdir(my_dir):
            shutil.rmtree(my_dir)
        os.mkdir(my_dir)
        shutil.copytree(os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count'),
                        os.path.join(my_dir, 'count'))
        bam_dir = os.path.join(my_dir, 'mapping')
        os.mkdir(bam_dir)
        lines = list()
        for source in glob.glob(os.path.join(self.tools['transfer_l'].output_dir, 'rnaseq_mapping/bam/*.bam')):
            fname = os.path.basename(source)
            link_name = os.path.join(bam_dir, fname)
            os.link(source, link_name)
            lines.append(os.path.join(self._sheet.output, 'other/mapping', fname) + '\n')
        open(os.path.join(my_dir, 'bam.list'), 'w').writelines(lines)
        anno_dir = os.path.join(my_dir, 'annotation')
        os.mkdir(anno_dir)
        for source in glob.glob(
                os.path.join(self.tools['transfer_l'].output_dir, 'large_gush/filter_by_express/filtered_file/*.gtf')):
            fname = os.path.basename(source)
            link_name = os.path.join(anno_dir, fname)
            os.link(source, link_name)
        rmats_dir = os.path.join(my_dir, 'rmats')
        shutil.copytree(self.modules['rmats'].output_dir, rmats_dir)
        gene_detail_dir = os.path.join(my_dir, 'gene_detail')
        shutil.copytree(self.modules['gene_detail'].output_dir, gene_detail_dir)
        if self.option("is_snp"):
            snp_dir=os.path.join(my_dir,"snp_vcf")
            os.mkdir(snp_dir)
            snp_vcf_dir = os.path.join(snp_dir,"final.vcf")
            if self.option("snp_method").lower() == "gatk":
                vcf_path = os.path.join(self.modules['whole_snp'].snp.vcffilter.output_dir,"final.vcf")
            elif self.option("snp_method").lower() == "sentieon":
                vcf_path = os.path.join(self.modules['whole_snp'].snp.vcffilter.output_dir, "final.vcf")
            else:
                vcf_path = os.path.join(self.modules['whole_snp'].snp.vcf_filter.output_dir, "final.vcf")
            if os.path.isfile(vcf_path):
                os.link(vcf_path, snp_vcf_dir)

    def move_chart_file(self):
        self.target_dir = self.output_dir
        os.makedirs(os.path.join(self.target_dir, '06Express/ExpVenn/'))

        file2uploads = [
            ("*long_*qc_*.pdf", 'mrna/02_Basic_Analysis/01_QC/'),
            ("*small_*qc_*.pdf", 'mirna/02_Basic_Analysis/01_QC/'),
            ("*small_clean_length_stat*.pdf", 'mirna/02_Basic_Analysis/01_QC/'),
            ("*circle_*qc_*.pdf", 'circrna/02_Basic_Analysis/01_QC/'),
            ("*long*map_satur*.pdf", 'mrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*long.map_coverage.line.pdf", 'mrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*long.seq_region_distribution*pdf", 'mrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*long.map_reads_stat*.pdf", 'mrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*small_seq_reads_distribution*pdf",'mirna/02_Basic_Analysis/02_Align/'),
            ("*circle*map_satur*.pdf*", 'circrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*circle.map_coverage.line.pdf", 'circrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*circle.seq_region_distribution*pdf", 'circrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*circle.map_reads_stat*.pdf", 'circrna/02_Basic_Analysis/02_Align/Quality_Assessment/'),
            ("*assemble_length_distribution*pdf", 'mrna/02_Basic_Analysis/03_Assemble/'),
            ("*new_transcriptome*assemble_distribution*pdf", 'mrna/02_Basic_Analysis/03_Assemble/'),
            ("*annot_*_stat*pdf", 'mrna/03_Annotation/02_Anno_Statistics/'),
            # ("*annot_*_stat*pdf", 'mrna/03_Annotation/02_Anno_Statistics/'),
            ("*lncRNA_predict*pdf", 'lncrna/01_lncRNA_Analysis/03_LncRNA_stat/'),
            ("*category_all*lnc_predict_stat*pdf", 'lncrna/01_lncRNA_Analysis/03_LncRNA_stat/'),
            ("*sample*mirna_predict_stat*pdf", 'mirna/03_sRNA_Analysis/03_sRNA_stat/'),
            ("*sRNAs_distribution*pdf", 'mirna/03_sRNA_Analysis/03_sRNA_stat/'),
            ("*allcircrna_predict_stat.column.pdf","circrna/03_circRNA_Predict/"),
            ("*mRNA*exp_distribution*pdf", 'mrna/04_Express/01_Exp_Annalysis/'),
            ("*mRNA_all*exp*venn*pdf", 'mrna/04_Express/01_Exp_Annalysis/'),
            ("*lncRNA*exp_distribution*pdf", 'lncrna/02_Express/01_Exp_Annalysis/'),
            ("*lncRNA_all*exp*venn*pdf", 'lncrna/02_Express/01_Exp_Annalysis/'),
            ("*circRNA*exp_distribution*pdf", 'circrna/04_Express/01_Exp_Annalysis/'),
            ("*circRNA_all*exp*venn*pdf", 'circrna/04_Express/01_Exp_Annalysis/'),
            ("*smallRNA*exp_distribution*pdf", 'mirna/04_Express/01_Exp_Annalysis/'),
            ("*smallRNA_all*exp*venn*pdf", 'mirna/04_Express/01_Exp_Annalysis/'),
            ("*long_all*exp*heat_corr*pdf", 'mrna/04_Express/02_Exp_Corr/'),
            ("*circle_all*exp*heat_corr*pdf", 'circrna/04_Express/02_Exp_Corr/'),
            ("*small_all*exp*heat_corr*pdf", 'mirna/04_Express/02_Exp_Corr/'),
            ("*long_all*exp_relation_pca*pdf", 'mrna/04_Express/03_Exp_PCA/'),
            ("*circle_all**exp_relation_pca*pdf", 'circrna/04_Express/03_Exp_PCA/'),
            ("*small_all*exp_relation_pca*pdf", 'mirna/04_Express/03_Exp_PCA/'),
            ("*mRNA*differential.summary*pdf",'mrna/05_Diff_Express/'),
            ("*mRNA*diff.volcano*pdf", 'mrna/05_Diff_Express/'),
            ("*mRNA*diff.scatter*pdf", 'mrna/05_Diff_Express/'),
            ("*lncRNA*differential.summary*pdf", 'lncrna/03_Diff_Express/'),
            ("*lncRNA*diff.volcano*pdf", 'lncrna/03_Diff_Express/'),
            ("*lncRNA*diff.scatter*pdf", 'lncrna/03_Diff_Express/'),
            ("*smallRNA*differential.summary*pdf", 'mirna/05_Diff_Express/'),
            ("*smallRNA*diff.volcano*pdf", 'mirna/05_Diff_Express/'),
            ("*smallRNA*diff.scatter*pdf", 'mirna/05_Diff_Express/'),
            ("*circRNA*differential.summary*pdf", 'circrna/05_Diff_Express/'),
            ("*circRNA*diff.volcano*pdf", 'circrna/05_Diff_Express/'),
            ("*circRNA*diff.scatter*pdf", 'circrna/05_Diff_Express/'),
            ("all_*splice_stat.*.pdf", 'mrna/06_AS'),
            ("*snp.*stat.*.pdf", "mrna/07_SNP_InDel/"),
            ("*first_bias_per*.pdf", "mirna/07_miRNA_Structure/01_miRNA_bias"),
            ("*all_loc_bias_per*.pdf", "mirna/07_miRNA_Structure/01_miRNA_bias"),
            ("*miRNAedit_distribution*.pdf", "mirna/07_miRNA_Structure/02_miRNA_edit")
        ]

        if self.option('is_as'):
            cmp_list = self.modules['rmats'].option("control_table").prop["cmp_list"]
            for cmps in cmp_list:
                file2uploads.append((
                    "{}_vs_{}*splice_stat.*.pdf".format(cmps[0], cmps[1]), 'mrna/06_AS/{}_vs_{}/'.format(cmps[0], cmps[1])))

        for filefrom, fileto in file2uploads:
            pdf_file = glob.glob(self.chart.work_dir + "/" + filefrom)
            for p in pdf_file:
                if os.path.exists(self.target_dir + "/" + fileto + "/" + os.path.basename(p)):
                    os.remove(self.target_dir + "/" + fileto + "/" + os.path.basename(p))
                try:
                    os.link(p, self.target_dir + "/" + fileto + "/" + os.path.basename(p))
                except:
                    self.logger.info("move_chart_error: ")
                    self.logger.info("pdf_path:{}".format(p))
                    self.logger.info("target_path : {}".format(self.target_dir + "/" + fileto + "/" + os.path.basename(p)))

        if self.run_diff_geneset_analysis:
            file2uploads_geneset = [
                ("geneset.cluster.heat_corr.pdf".format("{geneset_name}"),
                 "mrna/08_Diff_geneset_analysis/01_{}_Cluster_Analysis".format("{geneset_name}")),
                ("{}.cog_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "mrna/08_Diff_geneset_analysis/02_{}_COG_Annotation".format("{geneset_name}")),
                ("{}.go_annot.gene_set.column.pdf".format("{geneset_name}"),
                 "mrna/08_Diff_geneset_analysis/03_{}_GO_Annotation".format("{geneset_name}")),
                ("{}*.kegg_annot.*.column.pdf".format("{geneset_name}"),
                 "mrna/08_Diff_geneset_analysis/04_{}_KEGG_Annotation".format("{geneset_name}")),
                ("{}.go_enrich.gene_set.*.pdf".format("{geneset_name}"),
                 "mrna/08_Diff_geneset_analysis/05_{}_GO_Enrich".format("{geneset_name}")),
                ("{}*kegg_enrich.gene_set.*.pdf".format("{geneset_name}"),
                 "mrna/08_Diff_geneset_analysis/06_{}_KEGG_Enrich".format("{geneset_name}")),
                ("{}.*venn*.pdf".format("diff_genesets"),
                 "mrna/08_Diff_geneset_analysis/07_GenesetVenn")

            ]
            genesets = [os.path.basename(i) for i in glob.glob("{}/*vs*".format(self.modules['genesets_analysis'].output_dir))]

            # 以基因集为单位,专门为差异一键化分析准备
            for geneset in genesets:
                for filefrom, fileto in file2uploads_geneset:
                    pdf_file = glob.glob((self.chart.work_dir + "/" + filefrom).format(geneset_name=geneset))
                    for p in pdf_file:
                        if os.path.exists(
                                self.target_dir + "/" + fileto.format(geneset_name=geneset) + "/" + os.path.basename(
                                        p)):
                            os.remove(
                                self.target_dir + "/" + fileto.format(geneset_name=geneset) + "/" + os.path.basename(p))
                        try:
                            os.link(p, self.target_dir + "/" + fileto.format(
                                geneset_name=geneset) + "/" + os.path.basename(p))
                        except:
                            self.logger.info('cuolacuolaraw:{} new:{}'.format(p,
                                                                              self.target_dir + "/" + fileto.format(
                                                                                  geneset_name=geneset) + "/" + os.path.basename(
                                                                                  p)))

        pass

    def def_upload_dir(self):
        sdir = self.add_upload_dir(self.output_dir)
        way = self.long_task_info['options']['exp_way']
        try:
            l_method = self.long_task_info['options']['long_diff_method'].lower()
        except:
            l_method = self.long_task_info['options']['diff_method'].lower()
        if self.lib_dict['small']:
            s_method = self.small_task_info['options']['diff_method'].lower()
        else:
            s_method = l_method
        if self.lib_dict['circle']:
            c_method = self.circle_task_info['options']['diff_method'].lower()
        else:
            c_method = l_method

        if self.option("report_img"):
            s3 = self._sheet.output.split(":")[0]
            report_img_dir = self.chart.work_dir + '/png/'
            report_img_s3 = s3 + "://commonbucket/files/report_img/wholerna/" + self.task_id + "/"
            self.upload_to_s3(report_img_dir, report_img_s3)

        sdir.add_regexp_rules([
            [r'mrna/05_Diff_Express/.*_vs_.*_{}_G_anno.xls'.format(l_method), 'xls', '差异分析结果表', 0],
            [r'mrna/05_Diff_Express/.*_vs_.*_{}_T_anno.xls'.format(l_method), 'xls', '差异分析结果表', 0],
            [r'mrna/06_AS/.*_vs_.*', '', '差异组别', 0],
            [r'mrna/06_AS/.*_vs_.*/JC', '', 'JC定量', 0],
            [r'mrna/06_AS/.*_vs_.*/JC/A3SS_detail.xls', 'xls', 'A3SS可变剪切事件详情表（JC）', 0],
            [r'mrna/06_AS/.*_vs_.*/JC/SE_detail.xls', 'xls', 'SE可变剪切事件详情表（JC） ', 0],
            [r'mrna/06_AS/.*_vs_.*/JC/A5SS_detail.xls', 'xls', 'A5SS可变剪切事件详情表（JC） ', 0],
            [r'mrna/06_AS/.*_vs_.*/JC/MXE_detail.xls', 'xls', 'MXE可变剪切事件详情表（JC）', 0],
            [r'mrna/06_AS/.*_vs_.*/JC/RI_detail.xls', 'xls', 'RI可变剪切事件详情表（JC）', 0],
            [r'mrna/06_AS/.*_vs_.*/JCEC', '', 'JCEC定量', 0],
            [r'mrna/06_AS/.*_vs_.*/JCEC/A3SS_detail.xls', 'xls', 'A3SS可变剪切事件详情表（JC）', 0],
            [r'mrna/06_AS/.*_vs_.*/JCEC/SE_detail.xls', 'xls', 'SE可变剪切事件详情表（JC） ', 0],
            [r'mrna/06_AS/.*_vs_.*/JCEC/A5SS_detail.xls', 'xls', 'A5SS可变剪切事件详情表（JC） ', 0],
            [r'mrna/06_AS/.*_vs_.*/JCEC/MXE_detail.xls', 'xls', 'MXE可变剪切事件详情表（JC）', 0],
            [r'mrna/06_AS/.*_vs_.*/JCEC/RI_detail.xls', 'xls', 'RI可变剪切事件详情表（JC）', 0],
            [r'mrna/06_AS/.*_vs_.*/diff_event_stat.xls', '', '差异可变剪切事件统计表', 0],
            [r'mrna/06_AS/.*_vs_.*/diff_pattern_stat.JC.xls', '', '差异可变剪切模式变化统计表（JC）', 0],
            [r'mrna/06_AS/.*_vs_.*/diff_pattern_stat.JCEC.xls', '', '差异可变剪切模式变化统计表（JCEC）', 0],
            [r'lncrna/03_Diff_Express/.*_vs_.*_{}.xls'.format(l_method), 'xls', '差异分析结果表', 0],
            [r'mirna/06_miRNA_Target/target_align_detail/.*.txt.gz', 'gz', '靶基因比对详细信息', 0],
            [r'mirna/03_sRNA_Analysis/01_Known_miRNA/known_pre_structure/*.pdf', 'pdf', '已知miRNA前体结构图片', 0],
            [r'mirna/03_sRNA_Analysis/02_Novel_miRNA/novel_pre_structure/*.pdf', 'pdf', '新miRNA前体结构图片', 0],
            [r'mirna/05_Diff_Express/.*_vs_.*_{}.xls'.format(s_method), 'xls', '差异分析结果表', 0],
            [r'mirna/06_miRNA_Target/target_align_detail/*.txt.gz', 'gz', '中间文件', 0],
              ## 图片描述
            [r"mrna/02_Basic_Analysis/01_QC/.*raw_qc_qual\.box\.pdf", 'pdf', 'longRNA原始数据碱基质量分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*raw_qc_error\.line\.pdf", 'pdf', 'longRNA原始数据碱基错误率分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*raw_qc_base\.line\.pdf", 'pdf', 'longRNA原始数据碱基含量分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*clean_qc_qual\.box\.pdf", 'pdf', 'longRNA质控数据碱基质量分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*clean_qc_error\.line\.pdf", 'pdf', 'longRNA质控数据碱基错误率分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*clean_qc_base\.line\.pdf", 'pdf', 'longRNA质控数据碱基含量分布图', 0],
            [r"mirna/02_Basic_Analysis/01_QC/.*raw_qc_qual\.box\.pdf", 'pdf', 'smallRNA原始数据碱基质量分布图', 0],
            [r"mirna/02_Basic_Analysis/01_QC/.*raw_qc_error\.line\.pdf", 'pdf', 'smallRNA原始数据碱基错误率分布图', 0],
            [r"mirna/02_Basic_Analysis/01_QC/.*raw_qc_base\.line\.pdf", 'pdf', 'smallRNA原始数据碱基含量分布图', 0],
            [r"mirna/02_Basic_Analysis/01_QC/.*small_clean_length_stat\.columns\.pdf", 'pdf', 'smallRNA质控数据reads长度分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*raw_qc_qual\.box\.pdf", 'pdf', 'longRNA原始数据碱基质量分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*raw_qc_error\.line\.pdf", 'pdf', 'longRNA原始数据碱基错误率分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*raw_qc_base\.line\.pdf", 'pdf', 'longRNA原始数据碱基含量分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*clean_qc_qual\.box\.pdf", 'pdf', 'longRNA质控数据碱基质量分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*clean_qc_error\.line\.pdf", 'pdf', 'longRNA质控数据碱基错误率分布图', 0],
            [r"mrna/02_Basic_Analysis/01_QC/.*clean_qc_base\.line\.pdf", 'pdf', 'longRNA质控数据碱基含量分布图', 0],
            [r"circrna/02_Basic_Analysis/01_QC/.*raw_qc_qual\.box\.pdf", 'pdf', 'circleRNA原始数据碱基质量分布图', 0],
            [r"circrna/02_Basic_Analysis/01_QC/.*raw_qc_error\.line\.pdf", 'pdf', 'circleRNA原始数据碱基错误率分布图', 0],
            [r"circrna/02_Basic_Analysis/01_QC/.*raw_qc_base\.line\.pdf", 'pdf', 'circleRNA原始数据碱基含量分布图', 0],
            [r"circrna/02_Basic_Analysis/01_QC/.*clean_qc_qual\.box\.pdf", 'pdf', 'circleRNA质控数据碱基质量分布图', 0],
            [r"circrna/02_Basic_Analysis/01_QC/.*clean_qc_error\.line\.pdf", 'pdf', 'circleRNA质控数据碱基错误率分布图', 0],
            [r"circrna/02_Basic_Analysis/01_QC/.*clean_qc_base\.line\.pdf", 'pdf', 'circleRNA质控数据碱基含量分布图', 0],
            [r"mrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*map_satu.*\.pdf", 'pdf', 'longRNA测序饱和度曲线图', 0],
            [r"mrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*map_coverage.*\.pdf", 'pdf', 'longRNA测序覆盖度分布图', 0],
            [r"mrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*seq_region_distribution.*\.pdf", 'pdf', 'longRNA不同区域Reads分布统计饼图',0],
            [r"mrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*map_reads_stat.*\.pdf", 'pdf', 'longRNA不同染色体Reads分布统计柱状图',0],
            [r"mirna/02_Basic_Analysis/02_Align/.*small_seq_reads_distribution.*\.pdf", 'pdf', 'smallRNA不同染色体Reads分布统计柱状图', 0],
            [r"circrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*map_satu.*\.pdf", 'pdf', 'circleRNA测序饱和度曲线图', 0],
            [r"circrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*map_coverage.*\.pdf", 'pdf', 'circleRNA测序覆盖度分布图', 0],
            [r"circrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*seq_region_distribution.*\.pdf", 'pdf',
             'circleRNA不同区域Reads分布统计饼图', 0],
            [r"circrna/02_Basic_Analysis/02_Align/Quality_Assessment/.*map_reads_stat.*\.pdf", 'pdf',
             'circleRNA不同染色体Reads分布统计柱状图', 0],
            [r"mrna/02_Basic_Analysis/03_Assemble/.*\.assemble_length_distribution.*\.column\.pdf", 'pdf', '转录本长度分布柱状图', 0],
            [r"mrna/02_Basic_Analysis/03_Assemble/new_transcriptome\.assemble_distribution.*\.pdf", 'pdf', '新转录本类型分布饼图',0],
            [r"mrna/03_Annotation/02_Anno_Statistics/.*annot_.*_stat.*.\pdf", 'pdf', '功能注释统计柱状图', 0],
            [r"lncrna/01_lncRNA_Analysis/03_LncRNA_stat/.*new.*lncRNA_predict.*venn.\pdf", 'pdf', 'lncRNA预测Venn图', 0],
            [r"lncrna/01_lncRNA_Analysis/03_LncRNA_stat/.*new.*lncRNA_predict.*upset.\pdf", 'pdf', 'lncRNA预测Upset图', 0],
            [r"lncrna/01_lncRNA_Analysis/03_LncRNA_stat/.*lnc_predict_stat\.column\.pdf", 'pdf', 'lncRNA统计柱状图', 0],
            [r"mirna/03_sRNA_Analysis/03_sRNA_stat/.*mirna_predict_stat\.column\.pdf", 'pdf', 'miRNA统计图', 0],
            [r"mirna/03_sRNA_Analysis/03_sRNA_stat/.*sRNAs_distribution\.pie\.pdf", 'pdf', 'sRNAs统计图', 0],
            [r"circrna/03_circRNA_Predict/.*allcircrna_predict_stat\.column\.pdf", 'pdf', 'circRNA统计图', 0],
            [r"mrna/04_Express/01_Exp_Annalysis/.*exp_distribution\.box\.pdf", 'pdf', 'mRNA表达量分布盒型图', 0],
            [r"mrna/04_Express/01_Exp_Annalysis/.*exp_distribution\.density\.pdf", 'pdf', 'mRNA表达量分布密度图', 0],
            [r"mrna/04_Express/01_Exp_Annalysis/.*exp_distribution\.violin\.pdf", 'pdf', 'mRNA表达量分布小提琴图', 0],
            [r"mrna/04_Express/01_Exp_Annalysis/.*\.venn\.pdf", 'pdf', 'mRNAvenn图', 0],
            [r"lncrna/02_Express/01_Exp_Annalysis/.*exp_distribution\.box\.pdf", 'pdf', 'lncRNA表达量分布盒型图', 0],
            [r"lncrna/02_Express/01_Exp_Annalysis/.*exp_distribution\.density\.pdf", 'pdf', 'lncRNA表达量分布密度图', 0],
            [r"lncrna/02_Express/01_Exp_Annalysis/.*exp_distribution\.violin\.pdf", 'pdf', 'lncRNA表达量分布小提琴图', 0],
            [r"lncrna/02_Express/01_Exp_Annalysis/.*\.venn\.pdf", 'pdf', 'lncRNAvenn图', 0],
            [r"circrna/04_Express/01_Exp_Annalysis/.*exp_distribution\.box\.pdf", 'pdf', 'circleRNA表达量分布盒型图', 0],
            [r"circrna/04_Express/01_Exp_Annalysis/.*exp_distribution\.density\.pdf", 'pdf', 'circleRNA表达量分布密度图', 0],
            [r"circrna/04_Express/01_Exp_Annalysis/.*exp_distribution\.violin\.pdf", 'pdf', 'circleRNA表达量分布小提琴图', 0],
            [r"circrna/04_Express/01_Exp_Annalysis/.*\.venn\.pdf", 'pdf', 'circleRNAvenn图', 0],
            [r"mirna/04_Express/01_Exp_Annalysis/.*exp_distribution\.box\.pdf", 'pdf', 'mirna表达量分布盒型图', 0],
            [r"mirna/04_Express/01_Exp_Annalysis/.*exp_distribution\.density\.pdf", 'pdf', 'mirna表达量分布密度图', 0],
            [r"mirna/04_Express/01_Exp_Annalysis/.*exp_distribution\.violin\.pdf", 'pdf', 'mirna表达量分布小提琴图', 0],
            [r"mirna/04_Express/01_Exp_Annalysis/.*\.venn\.pdf", 'pdf', 'mirnavenn图', 0],
            [r"mrna/04_Express/03_Exp_PCA/.*all\.exp_relation.*\.pdf", 'pdf', 'longRNA样本间PCA图', 0],
            [r"mrna/04_Express/02_Exp_Corr/.*all\.exp\.heat_corr.*\.pdf", 'pdf', 'longRNA样本间相关性热图', 0],
            [r"circrna/04_Express/03_Exp_PCA/.*all\.exp_relation.*\.pdf", 'pdf', 'circleRNA样本间PCA图', 0],
            [r"circrna/04_Express/02_Exp_Corr/.*all\.exp\.heat_corr.*\.pdf", 'pdf', 'circleRNA样本间相关性热图', 0],
            [r"mirna/04_Express/03_Exp_PCA/.*all\.exp_relation.*\.pdf", 'pdf', 'miRNA样本间PCA图', 0],
            [r"mirna/04_Express/02_Exp_Corr/.*all\.exp\.heat_corr.*\.pdf", 'pdf', 'miRNA样本间相关性热图', 0],
            [r"mrna/05_Diff_Express/.*\.summary.*\.pdf", 'pdf', 'mRNA表达量差异统计图', 0],
            [r"mrna/05_Diff_Express/.*diff\.volcano.*\.pdf", 'pdf', 'mRNA表达量差异火山图', 0],
            [r"mrna/05_Diff_Express/.*diff\.scatter.*\.pdf", 'pdf', 'mRNA表达量差异散点图', 0],
            [r"lncrna/03_Diff_Express/.*\.summary.*\.pdf", 'pdf', 'lncRNA表达量差异统计图', 0],
            [r"lncrna/03_Diff_Express/.*diff\.volcano.*\.pdf", 'pdf', 'lncRNA表达量差异火山图', 0],
            [r"lncrna/03_Diff_Express/.*diff\.scatter.*\.pdf", 'pdf', 'lncRNA表达量差异散点图', 0],
            [r"lncrna/03_Diff_Express/.*\.summary.*\.pdf", 'pdf', 'lncRNA表达量差异统计图', 0],
            [r"lncrna/03_Diff_Express/.*diff\.volcano.*\.pdf", 'pdf', 'lncRNA表达量差异火山图', 0],
            [r"lncrna/03_Diff_Express/.*diff\.scatter.*\.pdf", 'pdf', 'lncRNA表达量差异散点图', 0],
            [r"mirna/05_Diff_Express/.*\.summary.*\.pdf", 'pdf', 'miRNA表达量差异统计图', 0],
            [r"mirna/05_Diff_Express/.*diff\.volcano.*\.pdf", 'pdf', 'miRNA表达量差异火山图', 0],
            [r"mirna/05_Diff_Express/.*diff\.scatter.*\.pdf", 'pdf', 'miRNA表达量差异散点图', 0],
            [r"circrna/05_Diff_Express/.*\.summary.*\.pdf", 'pdf', 'circleRNA表达量差异统计图', 0],
            [r"circrna/05_Diff_Express/.*diff\.volcano.*\.pdf", 'pdf', 'circleRNA表达量差异火山图', 0],
            [r"circrna/05_Diff_Express/.*diff\.scatter.*\.pdf", 'pdf', 'circleRNA表达量差异散点图', 0],
            [r"mrna/06_AS/.*splice_stat.*\.pdf", 'pdf', '可变剪切事件统计图', 0],
            [r"mrna/06_AS/.*_vs_.*/.*splice_stat\.*pie\.pdf", 'pdf', '组内差异可变剪切事件统计饼状图', 0],
            [r"mrna/06_AS/.*_vs_.*/.*splice_stat\.*column\.pdf", 'pdf', '组内差异可变剪切事件统计柱状图', 0],
            [r"mrna/07_SNP_InDel/.*snp\.pos_stat\.pie\.pdf", 'pdf', 'SNP不同区域分布饼图', 0],
            [r"mrna/07_SNP_InDel/.*snp\.type_stat\.column\.pdf", 'pdf', 'SNP类型统计图', 0],
            [r"mrna/07_SNP_InDel/.*snp\.type_stat\.pie\.pdf", 'pdf', 'SNP类型饼图', 0],
            [r"mrna/07_SNP_InDel/.*snp\.depth_stat\.column\.pdf", 'pdf', 'SNP深度统计图', 0],
            [r"mrna/07_SNP_InDel/.*snp\.depth_stat\.pie\.pdf", 'pdf', 'SNP深度饼图', 0],
            [r"mirna/07_miRNA_Structure/02_miRNA_edit/.*miRNAedit_distribution.*\.pdf", 'pdf', 'miRNA碱基编辑类型分布图', 0],


        ])

        sdir.add_relpath_rules([
            ['.', '', '流程分析结果目录', 0],
            ['mrna', '', 'mRNA结果目录', 0],
            ['mrna/01_Background', '', '项目背景结果目录', 0],
            ['mrna/01_Background/genome_info.xls', 'xls', '基因组注释信息表', 0],
            ['mrna/01_Background/sample_info.xls', 'xls', '样本信息表', 0],
            ['mrna/01_Background/software_info.xls', 'xls', '软件信息表', 0],
            ['mrna/02_Basic_Analysis', '', '基础分析结果目录', 0],
            ['mrna/02_Basic_Analysis/01_QC', '', '测序数据质控', 0],
            ['mrna/02_Basic_Analysis/01_QC/QC_stat.xls', 'xls', '测序数据统计表', 0],
            ['mrna/02_Basic_Analysis/02_Align', '', '序列比对分析', 0],
            ['mrna/02_Basic_Analysis/02_Align/align_stat.xls', 'xls', '比对结果统计表', 0],
            ['mrna/02_Basic_Analysis/02_Align/Quality_Assessment', '', '转录组质量评估', 0],
            ['mrna/02_Basic_Analysis/02_Align/Quality_Assessment/region_distribution.xls', 'xls', '不同区域Reads分布统计表', 0],
            ['mrna/02_Basic_Analysis/02_Align/Quality_Assessment/chr_distribution.xls', 'xls', '不同染色体Reads分布统计表', 0],
            ['mrna/02_Basic_Analysis/03_Assemble', '', '转录本组装', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence', '', '转录本序列文件', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence/new_transcript.gtf', 'gtf', '新转录本GTF文件', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence/all_transcript.gtf', 'gtf', '组装结果GTF文件', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence/all_transcript.fa', 'fasta', '组装结果序列', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence/new_transcript.fa', 'fasta', '新转录本序列', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence/all_cds.fa', 'fasta', 'CDS序列文件', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence/all_pep.fa', 'fasta', '蛋白序列文件', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/Sequence/all_id.xls', 'xls', '基因转录本蛋白ID对应关系文件', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/length_distribution_200.xls', 'xls', '转录本长度分布表-步长200', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/length_distribution_300.xls', 'xls', '转录本长度分布表-步长300', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/length_distribution_600.xls', 'xls', '转录本长度分布表-步长600', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/length_distribution_1000.xls', 'xls', '转录本长度分布表-步长1000', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/200.assemble_length_distribution.columns.pdf', 'pdf', '转录本长度分布柱状图-步长200', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/300.assemble_length_distribution.columns.pdf', 'pdf', '转录本长度分布柱状图-步长300', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/600.assemble_length_distribution.columns.pdf', 'pdf', '转录本长度分布柱状图-步长600', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/1000.assemble_length_distribution.columns.pdf', 'pdf', '转录本长度分布柱状图-步长1000', 0],
            ['mrna/02_Basic_Analysis/03_Assemble/classcode_stat.xls', 'xls', '新转录本类型统计表', 0],
            ['mrna/03_Annotation', '', '功能注释结果目录', 0],
            ['mrna/03_Annotation/01_Anno_Detail', '', '注释结果详情', 0],
            ['mrna/03_Annotation/01_Anno_Detail/ref_gene_anno_detail.xls', 'xls', '已知基因注释结果详情表', 0],
            ['mrna/03_Annotation/01_Anno_Detail/new_gene_anno_detail.xls', 'xls', '新基因注释结果详情表', 0],
            ['mrna/03_Annotation/01_Anno_Detail/all_gene_anno_detail.xls', 'xls', '所有基因注释结果详情表', 0],
            ['mrna/03_Annotation/01_Anno_Detail/ref_transcript_anno_detail.xls', 'xls', '已知转录本注释结果详情表', 0],
            ['mrna/03_Annotation/01_Anno_Detail/new_transcript_anno_detail.xls', 'xls', '新转录本注释结果详情表', 0],
            ['mrna/03_Annotation/01_Anno_Detail/all_transcript_anno_detail.xls', 'xls', '所有转录本注释结果详情表', 0],
            ['mrna/03_Annotation/02_Anno_Statistics', '', '注释结果统计', 0],
            ['mrna/03_Annotation/02_Anno_Statistics/ref_anno_stat.xls', 'xls', '已知基因/转录本结果统计表', 0],
            ['mrna/03_Annotation/02_Anno_Statistics/new_anno_stat.xls', 'xls', '新基因/转录本结果统计表', 0],
            ['mrna/03_Annotation/02_Anno_Statistics/all_anno_stat.xls', 'xls', '所有基因/转录本结果统计表', 0],
            ['mrna/04_Express', '', '表达量分析结果目录', 0],
            ['mrna/04_Express/01_Exp_Annalysis', '', '表达定量结果', 0],
            ['mrna/04_Express/01_Exp_Annalysis/gene_count_anno.xls', '', '基因count表达定量与功能注释结果表', 0],
            ['mrna/04_Express/01_Exp_Annalysis/transcript_count_anno.xls', '', '转录本count表达定量与功能注释结果表', 0],
            ['mrna/04_Express/01_Exp_Annalysis/gene_{}_anno.xls'.format(way), '',
             '基因{}表达定量与功能注释结果表'.format(way.upper()), 0],
            ['mrna/04_Express/01_Exp_Annalysis/transcript_{}_anno.xls'.format(way), '',
             '转录本{}表达定量与功能注释结果表'.format(way.upper()), 0],
            ['mrna/04_Express/02_Exp_Corr', '', '样本间相关性分析', 0],
            ['mrna/04_Express/02_Exp_Corr/sample_correlation.xls', 'xls', '样本间相关性系数表', 0],
            ['mrna/04_Express/03_Exp_PCA', '', '样本间PCA分析', 0],
            ['mrna/04_Express/03_Exp_PCA/explained_variance_ratio.xls', '', '主成分解释表', 0],
            ['mrna/05_Diff_Express', '', '表达量差异分析结果目录', 0],
            ['mrna/05_Diff_Express/total_diff_stat.G.{}_anno.xls'.format(l_method), 'xls', '基因差异表达与功能注释详情总表', 0],
            ['mrna/05_Diff_Express/total_diff_stat.T.{}_anno.xls'.format(l_method), 'xls', '基因差异表达与功能注释详情总表', 0],
            ['mrna/05_Diff_Express/diff_summary_G_{}_anno.xls'.format(l_method), 'xls', '基因差异表达与功能注释详情总表', 0],
            ['mrna/05_Diff_Express/diff_summary_T_{}_anno.xls'.format(l_method), 'xls', '基因差异表达与功能注释详情总表', 0],
            ['mrna/06_AS', '', '可变剪切分析结果目录', 0],
            ['mrna/06_AS/event_type.JC.xls', 'xls', '可变剪切事件统计表（JC）', 0],
            ['mrna/06_AS/event_type.JCEC.xls', 'xls', '可变剪切事件统计表（JCEC）', 0],
            ['mrna/07_SNP_InDel', '', 'SNP/InDel分析结果目录', 0],
            ['mrna/07_SNP_InDel/snp_anno.xls', 'xls', 'SNP分析结果注释详情表', 0],
            ['mrna/07_SNP_InDel/snp_position_distribution.xls', 'xls', 'SNP不同区域分布结果表', 0],
            ['mrna/07_SNP_InDel/snp_transition_tranversion_statistics.xls', 'xls', 'SNP类型统计结果表', 0],
            ['mrna/07_SNP_InDel/snp_depth_statistics.xls', 'xls', 'SNP深度统计结果表', 0],
            ['mrna/07_SNP_InDel/indel_anno.xls', 'xls', 'InDel分析结果注释详情表', 0],
            ['mrna/07_SNP_InDel/indel_position_distribution.xls', 'xls', 'InDel不同区域分布结果表', 0],
            ['mrna/08_Diff_geneset_analysis', '', '基因集一键化分析结果目录', 0],

            ['lncrna', '', 'lncRNA结果目录', 0],
            ['lncrna/01_lncRNA_Analysis', '', 'lncRNA分析结果目录', 0],
            ['lncrna/01_lncRNA_Analysis/01_Known_lncRNA', '', '已知lncRNA鉴定', 0],
            ['lncrna/01_lncRNA_Analysis/01_Known_lncRNA/known_lncrna.gtf', 'gtf', '已知lncRNA的GTF文件', 0],
            ['lncrna/01_lncRNA_Analysis/01_Known_lncRNA/known_lncrna.fa', 'fasta', '已知lncRNA序列', 0],
            ['lncrna/01_lncRNA_Analysis/01_Known_lncRNA/known_lncrna_detail.xls', 'xls', '已知lncRNA详情表', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA', '', '新lncRNA预测', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_mrna.gtf', 'gtf', '新mRNA的GTF文件', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_mrna.fa', 'fasta', '新mRNA序列', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_lncrna.gtf', 'gtf', '新lncRNA的GTF文件', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_lncrna.fa', 'fasta', '新lncRNA序列', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_lncrna_predict_detail.xls', 'xls', '新lncRNA预测详情表', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_lncrna_predict_stat.xls', 'xls', '新lncRNA预测统计表', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_mrna_ids.list', 'list', '中间文件', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/novel_lncrna_ids.list', 'list', '中间文件', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/cpc_output.xls', 'xls', 'CPC预测结果表', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/pfam_output.xls', 'xls', 'Pfam预测结果表', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/cnci_output.xls', 'xls', 'CNCI预测结果表', 0],
            ['lncrna/01_lncRNA_Analysis/02_Novel_lncRNA/cpat_output.xls', 'xls', 'CPAT预测结果表', 0],
            ['lncrna/01_lncRNA_Analysis/03_LncRNA_stat', '', 'lncRNA统计', 0],
            ['lncrna/01_lncRNA_Analysis/03_LncRNA_stat/lncrna_stat_in_sample.xls', '', 'lncRNA统计表', 0],
            ['lncrna/01_lncRNA_Analysis/03_LncRNA_stat/lncrna_stat_in_category.xls', '', 'lncRNA分类统计表', 0],
            ['lncrna/02_Express', '', '表达量分析结果目录', 0],
            ['lncrna/02_Express/01_Exp_Annalysis', '', '表达定量结果', 0],
            ['lncrna/02_Express/01_Exp_Annalysis/lncRNA_count.xls', '', 'count表达量分析结果表', 0],
            ['lncrna/02_Express/01_Exp_Annalysis/lncRNA_{}.xls'.format(way), '', '{}表达量分析结果表'.format(way.upper()), 0],
            ['lncrna/02_Express/02_Exp_Corr', '', '样本间相关性分析', 0],
            ['lncrna/02_Express/02_Exp_Corr/sample_correlation.xls', 'xls', '样本间相关性分析', 0],
            ['lncrna/02_Express/03_Exp_PCA', '', '样本间PCA分析', 0],
            ['lncrna/02_Express/03_Exp_PCA/explained_variance_ratio.xls', '', '主成分解释表', 0],
            ['lncrna/03_Diff_Express', '', '表达量差异分析结果目录', 0],
            ['lncrna/03_Diff_Express/total_diff_stat.{}.xls'.format(l_method), 'xls', '差异表达基因详情总表', 0],
            ['lncrna/03_Diff_Express/diff_summary_{}.xls'.format(l_method), 'xls', '差异表达基因统计详情表', 0],
            ['lncrna/04_LncRNA_Target', '', 'lncRNA靶基因预测结果目录', 0],
            ['lncrna/04_LncRNA_Target/tar_predict_detail.xls', 'xls', '靶基因预测详情表', 0],
            ['lncrna/05_LncRNA_Structure', '', 'lncRNA结构分析结果目录', 0],
            ['lncrna/05_LncRNA_Structure/01_LncRNA_Family', '', 'lncRNA家族分析', 0],
            ['lncrna/05_LncRNA_Structure/01_LncRNA_Family/lncRNA_family.xls', '', 'lncRNA家族信息表', 0],
            ['lncrna/05_LncRNA_Structure/02_miRNA_Pre', '', 'miRNA前体预测', 0],
            ['lncrna/05_LncRNA_Structure/02_miRNA_Pre/miRNA_precursor.xls', '', 'miRNA前体预测详情表', 0],

            ['mirna', '', 'miRNA结果目录', 0],
            ['mirna/01_Background', '', '项目背景结果目录', 0],
            ['mirna/01_Background/genome_info.xls', 'xls', '基因组注释信息表', 0],
            ['mirna/01_Background/sample_info.xls', 'xls', '样本信息表', 0],
            ['mirna/01_Background/software_info.xls', 'xls', '软件信息表', 0],
            ['mirna/02_Basic_Analysis', '', '基础分析结果目录', 0],
            ['mirna/02_Basic_Analysis/01_QC', '', '测序数据质控', 0],
            ['mirna/02_Basic_Analysis/01_QC/QC_stat.xls', 'xls', '测序数据统计表', 0],
            ['mirna/02_Basic_Analysis/02_Align', '', '序列比对分析', 0],
            ['mirna/02_Basic_Analysis/02_Align/align_stat.xls', 'xls', '比对结果统计表', 0],
            ['mirna/02_Basic_Analysis/02_Align/chr_distribution.xls', 'xls', '不同染色体Reads分布统计表', 0],
            ['mirna/03_sRNA_Analysis', '', 'sRNA分析结果目录', 0],
            ['mirna/03_sRNA_Analysis/01_Known_miRNA', '', '已知miRNA分析', 0],
            ['mirna/03_sRNA_Analysis/01_Known_miRNA/known_pre_structure', '', '已知miRNA前体结构图', 0],
            ['mirna/03_sRNA_Analysis/01_Known_miRNA/known_miRNA_detail.xls', 'xls', '已知miRNA详情表', 0],
            ['mirna/03_sRNA_Analysis/01_Known_miRNA/known_hairpin.fa', 'fasta', '已知miRNA前体序列', 0],
            ['mirna/03_sRNA_Analysis/01_Known_miRNA/known_mature.fa', 'fasta', '已知miRNA成熟体序列', 0],
            ['mirna/03_sRNA_Analysis/02_Novel_miRNA', '', '新miRNA预测', 0],
            ['mirna/03_sRNA_Analysis/02_Novel_miRNA/novel_pre_structure', '', '新miRNA前体结构图', 0],
            ['mirna/03_sRNA_Analysis/02_Novel_miRNA/novel_miRNA_detail.xls', 'xls', '新miRNA预测详情表', 0],
            ['mirna/03_sRNA_Analysis/02_Novel_miRNA/novel_hairpin.fa', 'fasta', '新miRNA前体序列', 0],
            ['mirna/03_sRNA_Analysis/02_Novel_miRNA/novel_mature.fa', 'fasta', '新miRNA成熟体序列', 0],
            ['mirna/03_sRNA_Analysis/03_sRNA_stat', '', 'sRNA统计', 0],
            ['mirna/03_sRNA_Analysis/03_sRNA_stat/miRNA_stat.xls', '', 'miRNA统计表', 0],
            ['mirna/03_sRNA_Analysis/03_sRNA_stat/sRNA_stat.xls', '', 'sRNA统计表', 0],
            ['mirna/04_Express', '', '表达量分析结果目录', 0],
            ['mirna/04_Express/01_Exp_Annalysis', '', '表达定量结果', 0],
            ['mirna/04_Express/01_Exp_Annalysis/miRNA_count.xls', '', 'count表达量分析结果表', 0],
            ['mirna/04_Express/01_Exp_Annalysis/miRNA_tpm.xls', '', 'TPM表达量分析结果表', 0],
            ['mirna/04_Express/02_Exp_Corr', '', '样本间相关性分析', 0],
            ['mirna/04_Express/02_Exp_Corr/sample_correlation.xls', 'xls', '样本间相关系数表', 0],
            ['mirna/04_Express/03_Exp_PCA', '', '样本间PCA分析', 0],
            ['mirna/04_Express/03_Exp_PCA/explained_variance_ratio.xls', '', '主成分解释表', 0],
            ['mirna/05_Diff_Express', '', '表达量差异分析结果目录', 0],
            ['mirna/05_Diff_Express/total_diff_stat_{}.xls'.format(s_method), 'xls', '差异表达基因详情总表', 0],
            ['mirna/05_Diff_Express/diff_summary_{}.xls'.format(s_method), 'xls', '差异表达基因统计详情表', 0],
            ['mirna/06_miRNA_Target', '', 'miRNA靶基因预测结果目录', 0],
            ['mirna/06_miRNA_Target/tar_predict_detail.xls', 'xls', '靶基因预测详情表', 0],
            ['mirna/06_miRNA_Target/target_seqs', '', '靶基因序列', 0],
            ['mirna/06_miRNA_Target/target_seqs/mRNA_known_target.fa', 'fasta', '已知mRNA靶基因序列', 0],
            ['mirna/06_miRNA_Target/target_seqs/mRNA_novel_target.fa', 'fasta', '新mRNA靶基因序列', 0],
            ['mirna/06_miRNA_Target/target_seqs/lncRNA_known_target.fa', 'fasta', '已知lncRNA靶基因序列', 0],
            ['mirna/06_miRNA_Target/target_seqs/lncRNA_novel_target.fa', 'fasta', '新lncRNA靶基因序列', 0],
            ['mirna/06_miRNA_Target/target_seqs/circRNA_known_target.fa', 'fasta', '已知circRNA靶基因序列', 0],
            ['mirna/06_miRNA_Target/target_seqs/circRNA_novel_target.fa', 'fasta', '新circRNA靶基因序列', 0],
            ['mirna/06_miRNA_Target/target_align_detail', '', '靶基因比对详细信息', 0],
            ['mirna/07_miRNA_Structure', '', 'miRNA结构分析结果目录', 0],
            ['mirna/07_miRNA_Structure/01_miRNA_bias', '', 'miRNA碱基偏好性分析', 0],
            ['mirna/07_miRNA_Structure/01_miRNA_bias/first_bias_per.xls', 'xls', '不同长度miRNA首位碱基偏好性统计表', 0],
            ['mirna/07_miRNA_Structure/01_miRNA_bias/loc_bias.xls', 'xls', 'miRNA不同位点碱基的偏好性统计表', 0],
            ['mirna/07_miRNA_Structure/02_miRNA_edit', '', 'miRNA碱基编辑分析', 0],
            ['mirna/07_miRNA_Structure/02_miRNA_edit/edit_detail.xls', 'xls', 'miRNA 碱基编辑详情表', 0],
            ['mirna/07_miRNA_Structure/03_miRNA_family', '', 'miRNA家族分析', 0],
            ['mirna/07_miRNA_Structure/03_miRNA_family/family_species.xls', 'xls', 'miRNA的家族信息表', 0],
            [r"mirna/07_miRNA_Structure/01_miRNA_bias/allall_loc_bias_per.column.pdf", 'pdf', 'all_miRNA不同位点碱基偏好性分布图', 0],
            [r"mirna/07_miRNA_Structure/01_miRNA_bias/allfirst_bias_per.column.pdf", 'pdf', 'all_miRNA首位碱基偏好性分布图', 0],
            [r"mirna/07_miRNA_Structure/01_miRNA_bias/knownall_loc_bias_per.column.pdf", 'pdf', 'known_miRNA不同位点碱基偏好性分布图', 0],
            [r"mirna/07_miRNA_Structure/01_miRNA_bias/knownfirst_bias_per.column.pdf", 'pdf', 'known_miRNA首位碱基偏好性分布图', 0],
            [r"mirna/07_miRNA_Structure/01_miRNA_bias/novelall_loc_bias_per.column.pdf", 'pdf', 'novel_miRNA不同位点碱基偏好性分布图', 0],
            [r"mirna/07_miRNA_Structure/01_miRNA_bias/novelfirst_bias_per.column.pdf", 'pdf', 'novel_miRNA首位碱基偏好性分布图', 0],

            ['circrna', '', 'circRNA结果目录', 0],
            ['circrna/01_Background', '', '项目背景结果目录', 0],
            ['circrna/01_Background/genome_info.xls', 'xls', '基因组注释信息表', 0],
            ['circrna/01_Background/sample_info.xls', 'xls', '样本信息表', 0],
            ['circrna/01_Background/software_info.xls', 'xls', '软件信息表', 0],
            ['circrna/02_Basic_Analysis', '', '基础分析结果目录', 0],
            ['circrna/02_Basic_Analysis/01_QC', '', '测序数据质控', 0],
            ['circrna/02_Basic_Analysis/01_QC/QC_stat.xls', 'xls', '测序数据统计表', 0],
            ['circrna/02_Basic_Analysis/02_Align', '', '序列比对分析', 0],
            ['circrna/02_Basic_Analysis/02_Align/align_stat.xls', 'xls', '比对结果统计表', 0],
            ['circrna/02_Basic_Analysis/02_Align/Quality_Assessment', '', '转录组质量评估', 0],
            ['circrna/02_Basic_Analysis/02_Align/Quality_Assessment/region_distribution.xls', 'xls', '不同区域Reads分布统计表',
             0],
            ['circrna/02_Basic_Analysis/02_Align/Quality_Assessment/chr_distribution.xls', 'xls', '不同染色体Reads分布统计表', 0],
            ['circrna/03_circRNA_Predict', '', 'circRNA预测结果目录', 0],
            ['circrna/03_circRNA_Predict/circRNA_predict_detail.xls', 'xls', 'circRNA预测详情表', 0],
            ['circrna/03_circRNA_Predict/circRNA_stat.xls', 'xls', 'circRNA分类统计表', 0],
            ['circrna/03_circRNA_Predict/circRNA.fa', 'fasta', 'circRNA预测序列', 0],
            ['circrna/04_Express', '', '表达量分析结果目录', 0],
            ['circrna/04_Express/01_Exp_Annalysis', '', '表达定量结果', 0],
            ['circrna/04_Express/01_Exp_Annalysis/circRNA_count.xls', '', 'count表达量分析结果表', 0],
            ['circrna/04_Express/01_Exp_Annalysis/circRNA_rpm.xls', '', 'RPM表达量分析结果表', 0],
            ['circrna/04_Express/02_Exp_Corr', '', '样本间相关性分析', 0],
            ['circrna/04_Express/02_Exp_Corr/sample_correlation.xls', 'xls', '样本间相关系数表', 0],
            ['circrna/04_Express/03_Exp_PCA', '', '样本间PCA分析', 0],
            ['circrna/04_Express/03_Exp_PCA/explained_variance_ratio.xls', '', '主成分解释表', 0],
            ['circrna/05_Diff_Express', '', '表达量差异分析结果目录', 0],
            ['circrna/05_Diff_Express/total_diff_stat_{}.xls'.format(c_method), 'xls', '差异表达基因详情总表', 0],
            ['circrna/05_Diff_Express/diff_summary_{}.xls'.format(c_method), 'xls', '差异表达基因统计详情表', 0],

            ['other', '', '其他信息结果目录', 0],
            #一键化
            # modify by fwy 20210927
            [r'mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis', '', '聚类分析', 0],
            [r'mrna/08_Diff_geneset_analysis/01_*/expression_matrix.xls', 'xls', '聚类热图分析表', 0],
            [r'mrna/08_Diff_geneset_analysis/01_*/subcluster.*', 'xls', '子聚类分析图', 0],
            [r'mrna/08_Diff_geneset_analysis/02_*COG_Annotation', '', 'COG功能注释', 0],
            [r'mrna/08_Diff_geneset_analysis/02_*COG.*/cog_class_table.xls', 'xls', 'COG分类统计表', 0],
            [r'mrna/08_Diff_geneset_analysis/03_*GO_Annotation', '', 'GO功能注释', 0],
            [r'mrna/08_Diff_geneset_analysis/03_*GO_Annotation/.*.xls', 'xls', 'GO分类统计表', 0],
            [r'mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation', '', 'KEGG功能注释', 0],
            [r'mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation/kegg_stat.xls', '', 'KEGG分类统计表', 0],
            [r'mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation/pathways.tar.gz', '', 'KEGG通路图', 0],
            [r'mrna/08_Diff_geneset_analysis/05_*GO_Enrich', '', 'GO功能富集', 0],
            [r'mrna/08_Diff_geneset_analysis/05_*GO_Enrich/go_enrich_stat.xls', '', 'GO富集分析统计表', 0],
            [r'mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich', '', 'KEGG富集分析', 0],
            [r'mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/kegg_erich_stat.xls', '', 'KEGG富集分析统计表', 0],
            [r'mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/pathways.tar.gz', '', 'KEGG富集通路图', 0],

            #图片
            ['mrna/08_Diff_geneset_analysis', '', '基因集分析结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/07_GenesetVenn', '', '差异基因集venn分析结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/07_GenesetVenn/*.pdf', '', '差异基因集venn图', 0],
            ['mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis', '', '基因集聚类分析结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/', '', '基因集聚类分析结果目录', 0],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/seq.cluster_tree.txt", "txt", "基因/转录本聚类树文件", 0, "211531"],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/seq.kmeans_cluster.txt", "txt", "基因/转录本聚类树文件", 0],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/sample.cluster_tree.txt", "txt", "样本聚类树文件", 0, "211532"],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/expression_matrix.xls", "xls", "聚类热图分析表", 0, "211533"],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/seq.subcluster_*.xls", "xls", "子聚类分析表", 0],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/heatmap.pdf", 'pdf', "聚类热图", 0],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/*heat_corr.pdf", 'pdf', "聚类热图", 0],
            ["mrna/08_Diff_geneset_analysis/01_*Cluster_Analysis/*.line.pdf", 'pdf', "子聚类折线图", 0],
            ['mrna/08_Diff_geneset_analysis/02_*COG_Annotation', '', '基因集聚COG分类结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/02_*COG_Annotation/', '', '基因集聚COG分类结果目录', 0],
            ["mrna/08_Diff_geneset_analysis/02_*COG_Annotation/cog_class_stat.xls", "", "COG分类统计表", 0, "211529"],
            ["mrna/08_Diff_geneset_analysis/02_*COG_Annotation/*.pdf", "pdf", "COG分类统计图", 0],
            ['mrna/08_Diff_geneset_analysis/03_*GO_Annotation', '', '基因集GO分类结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/03_*GO_Annotation', '', '基因集GO分类结果目录', 0],
            ["mrna/08_Diff_geneset_analysis/03_*GO_Annotation/go_class_stat.xls", "", "GO分类统计表", 0, "211527"],
            ["mrna/08_Diff_geneset_analysis/03_*GO_Annotation/*.pdf", "pdf", "GO分类统计图", 0],
            ['mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation', '', '基因集KEGG分类结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation', '', '基因集KEGG分类结果目录', 0],
            ["mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation/pathways.tar.gz", "", "KEGG通路图", 0, "211555"],

            ["mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation/kegg_class_stat.xls", "", "Pathway分类统计表", 0, "211557"],
            ["mrna/08_Diff_geneset_analysis/04_*KEGG_Annotation/*.pdf", "pdf", "KEGG分类统计图", 0],
            ['mrna/08_Diff_geneset_analysis/05_GO_Enrich', '', '基因集GO富集分析结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/05_*GO_Enrich', '', '基因集GO富集分析结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/05_*GO_Enrich/go_enrich_stat.xls', 'xls', 'GO富集分析统计表', 0, "211536"],
            ["mrna/08_Diff_geneset_analysis/05_*GO_Enrich/*bar.pdf", "pdf", "GO富集分析柱形图", 0],
            ["mrna/08_Diff_geneset_analysis/05_*GO_Enrich/*bar_line.pdf", "pdf", "GO富集分析柱形图(带折线)", 0],
            ["mrna/08_Diff_geneset_analysis/05_*GO_Enrich/*buble.pdf", "pdf", "GO富集分析气泡图(单基因集)", 0],
            ["mrna/08_Diff_geneset_analysis/05_*GO_Enrich/*buble2.pdf", "pdf", "GO富集分析气泡图(分散型)", 0],
            ['mrna/08_Diff_geneset_analysis/06_KEGG_Enrich', '', '基因集KEGG富集分析结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich', '', '基因集KEGG富集分析结果目录', 0],
            ['mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/kegg_enrich_stat.xls', '', 'KEGG富集分析统计表', 0, "211538"],
            ['mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/pathways.tar.gz', ' ', 'KEGG富集通路图', 0, "211539"],
            ["mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/*bar.pdf", "pdf", "KEGG富集分析柱形图", 0],
            ["mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/*bar_line.pdf", "pdf", "KEGG富集分析柱形图(带折线)", 0],
            ["mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/*buble.pdf", "pdf", "KEGG富集分析气泡图(单基因集)", 0],
            ["mrna/08_Diff_geneset_analysis/06_*KEGG_Enrich/*buble2.pdf", "pdf", "KEGG富集分析气泡图(分散型)", 0],

        ])

    def end(self):
        super(WholeTranscriptomeWorkflow, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the workflow. Just run this script to do test.
    """

    def test(self):
        from mbio.workflows.whole_transcriptome.whole_transcriptome import WholeTranscriptomeWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome',
            'options': {
                'lib_select': 'longRNA-seq,smallRNA-seq',
                'long_task_id': 'tsg_36023',
                'small_task_id': 'tsg_36024',

                'is_snp': 'True',
                'snp_method': 'sentieon',
                'is_as': 'True',

                'lt_cor_way': 'spearman',
                'lt_cor_cut': '0.9',
                'lt_sig_way': 'padjust',
                'lt_sig_cut': '0.05',
                'lt_adj_way': 'BH',
                'up_dis': '10',
                'down_dis': '10',
                'mirbase_category': 'Animal',
                'm_method': 'miranda',
                'm_min_support': '1',
                'm_miranda_energy': '-20.0',
                'm_miranda_score': '160.0',
                'm_miranda_strict': 'on',
                'm_ps_robot_score': '2.5',
                'm_rnahybird_energy': '-20.0',
                'm_rnahybird_num': '100',
                'm_rnahybird_pvalue': '0.01',
                'm_targetfinder_score': '4.0',

                'l_method': 'miranda',
                'l_min_support': '1',
                'l_miranda_energy': '-20.0',
                'l_miranda_score': '160.0',
                'l_miranda_strict': 'on',
                'l_ps_robot_score': '2.5',
                'l_rnahybird_energy': '-20.0',
                'l_rnahybird_num': '100',
                'l_rnahybird_pvalue': '0.01',
                'l_targetfinder_score': '4.0',

                'c_method': 'miranda',
                'c_min_support': '1',
                'c_miranda_energy': '-20.0',
                'c_miranda_score': '160.0',
                'c_miranda_strict': 'on',
                'c_ps_robot_score': '2.5',
                'c_rnahybird_energy': '-20.0',
                'c_rnahybird_num': '100',
                'c_rnahybird_pvalue': '0.01',
                'c_targetfinder_score': '4.0',
            }
        }
        wsheet = Sheet(data=data)
        wf = WholeTranscriptomeWorkflow(wsheet)
        wf.sheet.id = 'whole_transcriptome'
        wf.sheet.project_sn = 'whole_transcriptome'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
