# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import shutil
import unittest

import pandas as pd
from biocluster.config import Config
from biocluster.workflow import Workflow
import re
from mbio.packages.whole_transcriptome.utils import check_map_dict
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.whole_transcriptome.chart.chart import Chart
import glob
from shutil import copyfile

database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]


class DiffExpWorkflow(Workflow):
    '''
    last_modify: 2019.11.14
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffExpWorkflow, self).__init__(wsheet_object)
        LEVEL = ('G', 'T')
        CATEGORY = ('mRNA', 'lncRNA', 'miRNA', 'circRNA')
        KIND = ('all', 'ref', 'new')
        PROGRAM = ('DESeq2', 'DEGseq', 'edgeR', 'limma', 'Noiseq')
        FILTER = ('none', 'mean', 'max', 'min', 'sum', 'median')
        METHOD = ('Bonferroni', 'Holm', 'BH', 'BY')
        STAT_TYPE = ('pvalue', 'padjust')
        options = [
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': LEVEL[1]},
            {'name': 'category', 'type': 'string', 'default': CATEGORY[0]},
            {'name': 'kind', 'type': 'string', 'default': KIND[0]},
            {'name': 'background', 'type': 'string', 'default': None},
            {'name': 'group_id', 'type': 'string', 'default': None},
            {'name': 'group_dict', 'type': 'string', 'default': None},
            {'name': 'control_id', 'type': 'string', 'default': None},
            {'name': 'program', 'type': 'string', 'default': PROGRAM[0]},
            {'name': 'filter', 'type': 'string', 'default': FILTER[0]},
            {'name': 'threshold', 'type': 'float', 'default': 0.0},
            {'name': 'method', 'type': 'string', 'default': METHOD[2]},
            {'name': 'stat_type', 'type': 'string', 'default': STAT_TYPE[1]},
            {'name': 'stat_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'fc', 'type': 'float', 'default': 2.0},
            {'name': 'prob', 'type': 'float', 'default': 0.8}, #新增参数，NOIseq时用来判断是否限制的参数
            {'name': 'is_batch', 'type': 'bool', 'default': False},
            {'name': 'has_batch', 'type': 'bool', 'default': True},
            {'name': 'batch_matrix', 'type': 'infile', 'format': "whole_transcriptome.common"},
            {'name': 'is_rmbe', 'type': 'string', 'default': 'false'},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None},
            {'name': 'exp_id', 'type': 'string', 'default': None}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.database = Config().get_mongo_client(mtype='whole_transcriptome')[
            Config().get_mongo_dbname('whole_transcriptome')]
        self.task_info = self.database['task'].find_one({'task_id': self.option('task_id')})
        self.inter_dirs = []

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(DiffExpWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        self.run_transfer()
        # self.run_transfer()
        super(DiffExpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="diff", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_transfer(self):
        self.transfer = self.add_tool('whole_transcriptome.transfer')
        indir = os.path.join(self.task_info['output'], self.option('category').lower())
        if not indir.endswith("/"):
            indir = indir + "/"
        self.logger.info(indir)
        self.transfer.set_options({'indir': indir})
        self.transfer.on('end', self.run_exp_build)
        self.transfer.run()

    def run_exp_build(self):
        self.exp_build = self.add_tool('whole_transcriptome_v1_1.formation.exp_build')
        opts = {
            'task_id': self.option('task_id'),
            'level': self.option('level'),
            'kind': self.option('kind'),
            'group_dict': self.option('group_dict'),
            'control_id': self.option('control_id'),
            'is_rmbe': self.option('is_rmbe')
        }
        if self.option('background') == 'mRNA,lncRNA':
            opts['library'] = 'long'
        else:
            opts['category'] = self.option('category')
        self.exp_build.set_options(opts)
        self.exp_build.on('end', self.run_diff_prepare)
        self.exp_build.run()

    def run_diff_prepare(self):
        self.diff_prepare = self.add_tool('whole_transcriptome.formation.diff_prepare')
        self.diff_prepare.set_options({
            'task_id': self.option('task_id'),
            'level': self.option('level'),
            'category': self.option('category'),
            'kind': self.option('kind'),
            'background': self.option('background')
        })
        self.diff_prepare.on('end', self.run_diff_exp)
        self.diff_prepare.run()

    def run_diff_exp(self):
        self.diff_exp = self.add_module('whole_transcriptome_v1_1.diff_exp')
        opts = {
            'program': self.option('program'),
            'count_matrix': self.diff_prepare.option('count_matrix'),
            'group_table': self.exp_build.option('group_table'),
            'control_table': self.exp_build.option('control_table'),
            'exp_matrix': self.exp_build.option('exp_matrix'),
            'kind_table': self.diff_prepare.option('kind_table'),
            'filter': self.option('filter'),
            'threshold': self.option('threshold'),
            'fc': self.option('fc'),
            'is_batch': self.option('is_batch')
        }
        if self.option('is_batch') == False:
            if self.option("program").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                opts.update({
                    'method': self.option('method').lower(),
                    'stat_type': self.option('stat_type'),
                    'stat_cutoff': self.option('stat_cutoff'),
                })
            if self.option('program').lower() in ['noiseq']:
                opts.update(
                    prob=self.option('prob')
                )

        else:
            if self.option("program").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                opts.update({
                    'method': self.option('method').lower(),
                    'stat_type': self.option('stat_type'),
                    'stat_cutoff': self.option('stat_cutoff'),
                    'is_batch': self.option('is_batch'),
                    'has_batch': self.option('has_batch')
                })
            if self.option('program').lower() in ['noiseq']:
                opts.update(
                    prob=self.option('prob'),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                )
            if self.option('has_batch') == True:
                opts.update(
                    batch_matrix=self.option('batch_matrix')
                )

        self.diff_exp.set_options(opts)
        self.diff_exp.on('end', self.set_output)
        self.diff_exp.run()

    def set_output(self):
        for fname in os.listdir(self.diff_exp.output_dir):
            src = os.path.join(self.diff_exp.output_dir, fname)
            dst = os.path.join(self.output_dir, fname)
            shutil.copy(src, dst)
        self.split_diff()
        self.set_db()

    def split_diff(self):
        if self.option('background') == 'mRNA,lncRNA' and self.option('level') == 'T':
            exp_id = self.database['exp'].find_one({'task_id': self.option('task_id'), 'level': 'T'})['main_id']
            cursor = self.database['exp_detail'].find({'exp_id': exp_id, 'category': {'$in': ['mRNA', 'lncRNA']}})
            data = [{'seq_id': document['transcript_id'], 'category': document['category']} for document in cursor]
            relate_df = pd.DataFrame(data).set_index('seq_id')
            for fname in os.listdir(self.output_dir):
                df = pd.read_table(os.path.join(self.output_dir, fname), index_col=0)
                df = df.join(relate_df)
                df.index.name = 'seq_id'
                if self.option('category') == 'mRNA':
                    mrna_df = df[df['category'] == 'mRNA']
                    mrna_df = mrna_df.drop(['category'], axis=1)
                    mrna_df.to_csv(os.path.join(self.output_dir, fname), sep='\t')
                elif self.option('category') == 'lncRNA':
                    lncrna_df = df[df['category'] == 'lncRNA']
                    lncrna_df = lncrna_df.drop(['category'], axis=1)
                    lncrna_df.to_csv(os.path.join(self.output_dir, fname), sep='\t')

        if self.option('background') == 'mRNA,lncRNA' and self.option('level') == 'G':
            exp_id = self.database['exp'].find_one({'task_id': self.option('task_id'), 'level': 'G'})['main_id']
            cursor = self.database['exp_detail'].find({'exp_id': exp_id, 'category': {'$in': ['mRNA', 'lncRNA']}})
            data = [{'seq_id': document['gene_id'], 'category': document['category']} for document in cursor]
            relate_df = pd.DataFrame(data).set_index('seq_id')
            for fname in os.listdir(self.output_dir):
                df = pd.read_table(os.path.join(self.output_dir, fname), index_col=0)
                df = df.join(relate_df)
                df.index.name = 'seq_id'
                if self.option('category') == 'mRNA':
                    mrna_df = df[df['category'] == 'mRNA']
                    mrna_df = mrna_df.drop(['category'], axis=1)
                    mrna_df.to_csv(os.path.join(self.output_dir, fname), sep='\t')
                elif self.option('category') == 'lncRNA':
                    lncrna_df = df[df['category'] == 'lncRNA']
                    lncrna_df = lncrna_df.drop(['category'], axis=1)
                    lncrna_df.to_csv(os.path.join(self.output_dir, fname), sep='\t')

    def set_db(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        api = self.api.api('whole_transcriptome.expression_new')
        map_dict = dict()
        diff_dir = self.output_dir
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
        task_id = self.option('task_id')
        category = self.option('category')
        level = self.option('level')
        kind = self.option('kind')
        background = self.option('background')
        threshold = 'NA' if self.option('filter') == 'none' else self.option('threshold')
        group_id = self.option('group_id')
        group_dict = api.get_group_dict(group_id)
        control_id = self.option('control_id')
        stat_type = self.option('stat_type')
        stat_cutoff = self.option('stat_cutoff')
        fc = self.option('fc')
        diff_method = self.option('program')
        correct_method = self.option('method')
        way = api.db['exp'].find_one({'task_id': task_id, 'level': level})['way']
        project_sn = self.sheet.project_sn
        diff_id = self.option('main_id')
        api.add_diff(
            map_dict=map_dict,
            category=category,
            level=level,
            kind=kind,
            background=background,
            filter=self.option('filter'),
            threshold=self.option('threshold'),
            group_id=group_id,
            group_dict=group_dict,
            control_id=control_id,
            stat_type=stat_type,
            stat_cutoff=stat_cutoff,
            fc=fc,
            diff_method=diff_method,
            correct_method=correct_method,
            way=way,
            task_id=task_id,
            project_sn=project_sn,
            exp_id=self.option('exp_id'),
            diff_id=diff_id)
        self.set_upload()

    def Chart(self,upload_dir):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        level = self.option('level')
        group_dict = json.loads(self.option("group_dict"))
        exp = self.diff_prepare.option('count_matrix').path
        map_dict = dict()
        category = self.option('category')
        diff_dir = self.output_dir
        for fname in os.listdir(diff_dir):
            if fname.endswith('.detail.txt'):
                map_dict[fname[:-11]] = os.path.join(diff_dir, fname)
            elif 'summary' in fname:
                map_dict['summary'] = os.path.join(diff_dir, fname)
            elif 'volcano' in fname:
                map_dict['volcano'] = os.path.join(diff_dir, fname)
            elif 'scatter' in fname:
                map_dict['scatter'] = os.path.join(diff_dir, fname)
        summary = map_dict["summary"]
        exp_count = pd.read_table(exp).shape[0]
        volcano = map_dict["volcano"]
        scatter = map_dict["scatter"]
        cmp_list = list(set(pd.read_table(scatter)["compare"]))
        chart.chart_diff_summary(volcano, summary, level, rna_type=category,
                                 group_dict=group_dict, cmp_list=cmp_list,
                                 exp_count=exp_count)
        chart.chart_diff_volcano(volcano, level, rna_type=category)
        chart.chart_diff_scatter(scatter, level, rna_type=category)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            copyfile(p, os.path.join(upload_dir, os.path.basename(p)))
            #
            # if os.path.exists(
            #         self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution")):
            #     os.remove(
            #         self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))
            # copyfile(p, self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))

    def set_upload(self):

        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/02 Diff_Express')
        # rna_dir = os.path.join(self.transfer.output_dir, '{}'.format(self.option('category').lower()))
        rna_dir = self.transfer.output_dir
        way = self.task_info['long_task']['options']['exp_way']
        if self.option('category') == 'mRNA':
            from mbio.packages.whole_transcriptome.catalogue import mrna as rna
            if self.option('level') == 'T':
                map_dict = {
                    't_count': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/transcript_count_anno.xls'),
                    't_exp': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/transcript_{}_anno.xls'.format(way)),
                    't_anno': os.path.join(rna_dir, '03_Annotation/01_Anno_Detail/all_transcript_anno_detail.xls')}
            elif self.option('level') == 'G':
                map_dict = {
                    'g_count': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/gene_count_anno.xls'),
                    'g_exp': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/gene_{}_anno.xls'.format(way)),
                    'g_anno': os.path.join(rna_dir, '03_Annotation/01_Anno_Detail/all_gene_anno_detail.xls')}
        elif self.option('category') == 'lncRNA':
            from mbio.packages.whole_transcriptome.catalogue import lncrna as rna
            map_dict = {'t_count': os.path.join(rna_dir, '02_Express/01_Exp_Annalysis/lncRNA_count.xls'),
                        't_exp': os.path.join(rna_dir, '02_Express/01_Exp_Annalysis/lncRNA_{}.xls'.format(way))}
        elif self.option('category') == 'miRNA':
            from mbio.packages.whole_transcriptome.catalogue import mirna as rna
            map_dict = {'s_count': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/miRNA_count.xls'),
                        's_exp': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/miRNA_tpm.xls')}
        elif self.option('category') == 'circRNA':
            from mbio.packages.whole_transcriptome.catalogue import circrna as rna
            map_dict = {'c_count': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/circRNA_count.xls'),
                        'c_exp': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/circRNA_rpm.xls')}
        upload_dir = os.path.join(self.work_dir, 'upload')
        if not os.path.isdir(upload_dir):
            os.mkdir(upload_dir)
        self.Chart(upload_dir)
        rna.set_diff_express(map_dict, self.option('task_id'), upload_dir, self.option('main_id'))
        if os.path.exists(os.path.join(upload_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(upload_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(upload_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(upload_dir)
        self.inter_dirs = [
            ["02 Diff_Express", "", "表达量差异分析结果目录",0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '{}{}表达量差异分析文件'.format(
                self.option('category'), {'G': '基因', 'T': '转录本'}[self.option('level')]), 0],
            ['total_diff_stat.*.xls', 'xls', '表达量差异详情总表',0],
            ['*_vs_*.xls', 'xls', '表达量差异结果表', 0],
            ['*summary*.xls', 'xls', '表达量差异统计表', 0],
            # ['diff_stat_*.xls', 'xls', '表达量差异统计表', 0],
            # ['diff_stat_*_anno.xls', 'xls', '表达量差异统计详情表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'\S+_vs_\S+_{}_anno.xls'.format(self.option('program').lower()), 'xls', '表达量差异结果表', 0],
            ['.*scatter.pdf', 'pdf', '表达量差异散点图', 0],
            ['.*volcano.pdf', 'pdf', '表达量差异火山图', 0],
            ['.*bar.pdf', 'pdf', '表达量差异统计柱形图', 0],
            ['.*bar2.pdf', 'pdf', '表达量差异统计堆叠图', 0]
        ])
        self.end()

    def end(self):
        super(DiffExpWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.whole_transcriptome_v1_1.report.diff_exp import DiffExpWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome_v1_1.report.diff_exp',
            'options': {
                'task_id': 'tsg_36088',
                'level': 'T',
                'category': 'mRNA',
                'kind': 'ref',
                'background': 'mRNA,lncRNA',
                'group_id': '5dd5152817b2bf7328f25108',
                'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}),
                'control_id': '5dd5152817b2bf7328f25109',
                'program': 'edgeR',
                'filter': 'none',
                'threshold': '0.0',
                'method': 'BH',
                'stat_type': 'padjust',
                'stat_cutoff': '0.05',
                'fc': '2.0',
                'is_batch': False,
                # 'has_batch': True,
                # 'batch_matrix': "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/whole_transcriptome_v2/batch/batch"
            }
        }
        wsheet = Sheet(data=data)
        wf = DiffExpWorkflow(wsheet)
        wf.sheet.project_sn = '188_5dba6f542345b'
        wf.sheet.task_id = 'tsg_36088'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
