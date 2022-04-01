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
        PROGRAM = ('DESeq2', 'DEGseq', 'edgeR')
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
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.task_info = database['task'].find_one({'task_id': self.option('task_id')})
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
        super(DiffExpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="diff", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_transfer(self):
        self.transfer = self.add_tool('whole_transcriptome.transfer')
        indir = self.task_info['output'] + '/{}/'.format(self.option('category').lower())
        self.logger.info(indir)
        self.transfer.set_options({'indir': indir})
        self.transfer.on('end', self.run_exp_build)
        self.transfer.run()

    def run_exp_build(self):
        self.exp_build = self.add_tool('whole_transcriptome.formation.exp_build')
        opts = {
            'task_id': self.option('task_id'),
            'level': self.option('level'),
            'kind': self.option('kind'),
            'group_dict': self.option('group_dict'),
            'control_id': self.option('control_id')
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
        self.diff_exp = self.add_module('whole_transcriptome.diff_exp')
        self.diff_exp.set_options({
            'program': self.option('program'),
            'count_matrix': self.diff_prepare.option('count_matrix'),
            'group_table': self.exp_build.option('group_table'),
            'control_table': self.exp_build.option('control_table'),
            'exp_matrix': self.exp_build.option('exp_matrix'),
            'kind_table': self.diff_prepare.option('kind_table'),
            'filter': self.option('filter'),
            'threshold': self.option('threshold'),
            'method': self.option('method').lower(),
            'stat_type': self.option('stat_type'),
            'stat_cutoff': self.option('stat_cutoff'),
            'fc': self.option('fc'),
        })
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
            exp_id = database['exp'].find_one({'task_id': self.option('task_id'), 'level': 'T'})['main_id']
            cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': {'$in': ['mRNA', 'lncRNA']}})
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
            exp_id = database['exp'].find_one({'task_id': self.option('task_id'), 'level': 'G'})['main_id']
            cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': {'$in': ['mRNA', 'lncRNA']}})
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
        api = self.api.api('whole_transcriptome.expression')
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
            diff_id=diff_id)
        self.set_upload()

    def Chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        level = self.option('level')
        group_dict = json.loads(self.option("group_dict"))
        exp = self.diff_prepare.option('count_matrix')
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
            copyfile(p,os.path.join(self.output_dir,os.path.basename(p)))
            #
            # if os.path.exists(
            #         self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution")):
            #     os.remove(
            #         self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))
            # copyfile(p, self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))



    def set_upload(self):
        self.Chart()
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/02_Diff_Express')
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
            if os.path.exists(map_dict['s_count']):
                pass
            else:
                map_dict['s_count'] = os.path.join(self.diff_prepare.output_dir, "T.miRNA.all.txt")

            if os.path.exists(map_dict['s_exp']):
                pass
            else:
                map_dict['s_exp'] = os.path.join(self.exp_build.output_dir, "T.miRNA.all.txt")
        elif self.option('category') == 'circRNA':
            from mbio.packages.whole_transcriptome.catalogue import circrna as rna
            map_dict = {'c_count': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/circRNA_count.xls'),
                        'c_exp': os.path.join(rna_dir, '04_Express/01_Exp_Annalysis/circRNA_rpm.xls')}
        upload_dir = os.path.join(self.work_dir, 'upload')
        if not os.path.isdir(upload_dir):
            os.mkdir(upload_dir)
        rna.set_diff_express(map_dict, self.option('task_id'), upload_dir, self.option('main_id'))
        if os.path.exists(os.path.join(upload_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(upload_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(upload_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(upload_dir)
        self.inter_dirs = [
            ["02_Diff_Express", "", "表达量差异分析结果目录",0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '{}{}表达量差异分析文件'.format(
                self.option('category'), {'G': '基因', 'T': '转录本'}[self.option('level')]), 0],
            ['total_diff_stat_{}_anno.xls'.format(self.option('program')), 'xls', '表达量差异详情总表', 0],
            ['diff_stat_{}.xls'.format(self.option('program')), 'xls', '表达量差异统计表', 0],
            ['diff_stat_{}_anno.xls'.format(self.option('program')), 'xls', '表达量差异统计详情表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'\S+_vs_\S+_{}_anno.xls'.format(self.option('program').lower()), 'xls', '表达量差异结果表', 0]
        ])
        self.end()

    def end(self):
        super(DiffExpWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.whole_transcriptome.report.diff_exp import DiffExpWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.report.diff_exp',
            'options': {
                'task_id': 'tsg_36088',
                'level': 'T',
                'category': 'mRNA',
                'kind': 'ref',
                'background': 'mRNA,lncRNA',
                'group_id': '5dccfc9b17b2bf67f7f25ced',
                'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}),
                'control_id': '5dccfc9b17b2bf67f7f25cee',
                'program': 'DESeq2',
                'filter': 'none',
                'threshold': '0.0',
                'method': 'BH',
                'stat_type': 'padjust',
                'stat_cutoff': '0.05',
                'fc': '2.0',
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
