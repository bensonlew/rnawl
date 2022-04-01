# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import os
import shutil
import sys

import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId

from mbio.packages.lnc_rna.copy_file import CopyFile

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('lncRNA')
logger.setLevel(logging.DEBUG)

database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]


def set_lncrna_analysis(map_dict, task_id, output_dir,lnc_ref=True):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    for sub_dir in ['01_Known_lncRNA', '02_Novel_lncRNA', '03_LncRNA_stat']:
        if not os.path.isdir(os.path.join(output_dir, sub_dir)):
            os.mkdir(os.path.join(output_dir, sub_dir))
    for fname in ['known_lncrna.fa', 'known_lncrna.gtf', 'known_lncrna_detail.xls']:
        CopyFile().linkfile(os.path.join(map_dict['large_gush_dir'], 'known_lnc_identify', fname),
                            os.path.join(output_dir, '01_Known_lncRNA', fname), mode="cp")
    for fname in ['cnci_output.txt', 'cpat_output.txt', 'cpc_output.txt', 'novel_lncrna.fa', 'novel_lncrna.gtf',
                  'novel_lncrna_ids.list', 'novel_lncrna_predict_detail.xls', 'novel_mrna.fa', 'novel_mrna.gtf',
                  'novel_mrna_ids.list', 'pfam_output.txt']:
        os.path.join(map_dict['large_gush_dir'], 'filter_by_express/filtered_lncnovel', fname)
        CopyFile().linkfile(os.path.join(map_dict['large_gush_dir'], 'filter_by_express/filtered_lncnovel', fname),
                            os.path.join(output_dir, '02_Novel_lncRNA', fname), mode="cp")
    for prefix in ('cnci_output', 'cpat_output', 'cpc_output', 'pfam_output'):
        if os.path.isfile(os.path.join(output_dir, '02_Novel_lncRNA', '{}.txt'.format(prefix))):
            os.rename(os.path.join(output_dir, '02_Novel_lncRNA', '{}.txt'.format(prefix)),
                      os.path.join(output_dir, '02_Novel_lncRNA', '{}.xls'.format(prefix)))
    for fname in ['lncrna_stat_in_sample.xls', 'lncrna_stat_in_category.xls']:
        os.path.join(map_dict['large_gush_dir'], 'lncrna_stat', fname)
        CopyFile().linkfile(os.path.join(map_dict['large_gush_dir'], 'lncrna_stat', fname),
                            os.path.join(output_dir, '03_LncRNA_stat', fname), mode="cp")
    sdf = pd.read_table(os.path.join(output_dir, '03_LncRNA_stat/lncrna_stat_in_sample.xls'))
    if 'type' in sdf:
        sdf = sdf[sdf['type'] == 'LT'].reindex(['sample_name', 'known_num', 'new_num', 'total_num'], axis=1)
        sdf = sdf.rename({'sample_name': 'Sample', 'known_num': 'Known lncRNA num', 'new_num': 'Novel lncRNA num',
                          'total_num': 'Total'}, axis=1)
    sdf.to_csv(os.path.join(output_dir, '03_LncRNA_stat/lncrna_stat_in_sample.xls'), sep='\t', index=False)
    cdf = pd.read_table(os.path.join(output_dir, '03_LncRNA_stat/lncrna_stat_in_category.xls'))
    if 'type' in sdf:
        cdf = cdf[cdf['type'] == 'LT'].reindex(['biotype', 'known', 'novel', 'total'], axis=1)
        cdf = cdf.rename({'biotype': 'Type', 'known': 'Known lncRNA num', 'novel': 'Novel lncRNA num',
                          'total': 'Total'}, axis=1)
    cdf.to_csv(os.path.join(output_dir, '03_LncRNA_stat/lncrna_stat_in_category.xls'), sep='\t', index=False)
    if not lnc_ref:
        shutil.rmtree(os.path.join(output_dir, '01_Known_lncRNA'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_express(map_dict, task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_exp_annalysis(map_dict, task_id, os.path.join(output_dir, '01_Exp_Annalysis'))
    set_exp_corr(task_id, os.path.join(output_dir, '02_Exp_Corr'))
    set_exp_pca(task_id, os.path.join(output_dir, '03_Exp_PCA'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_exp_annalysis(map_dict, task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    samples = [document['old_name'] for document in database['specimen'].find({'task_id': task_id, 'library': 'long'})]
    main_dict = database['exp'].find_one({'task_id': task_id, 'level': 'T'})
    exp_id, way = main_dict['main_id'], main_dict['way']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'lncRNA'})
    t_detail_df = pd.DataFrame(list(cursor)).set_index('transcript_id')
    t_detail_df = t_detail_df.reindex(['gene_id'] + samples, axis=1)
    t_exp_df = t_detail_df
    t_exp_df.to_csv(os.path.join(output_dir, 'lncRNA_{}.xls'.format(way)), sep='\t')
    row_index = t_exp_df.index
    t_count_df = pd.read_table(map_dict['t_count'], index_col=0)
    t_count_df.index.name = 'transcript_id'
    t_count_df = t_count_df.reindex(row_index)
    t_count_df.to_csv(os.path.join(output_dir, 'lncRNA_count.xls'), sep='\t')


def set_exp_corr(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    main_dict = database['exp_corr'].find_one({'task_id': task_id, 'library': 'long'})
    corr_id, samples = main_dict['main_id'], main_dict['samples']
    cursor = database['exp_corr_detail'].find({'corr_id': corr_id})
    df = pd.DataFrame(list(cursor)).set_index('sample')
    df = df.reindex(samples, axis=1)
    df.to_csv(os.path.join(output_dir, 'sample_correlation.xls'), sep='\t')


def set_exp_pca(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    ratio_dict = database['exp_pca'].find_one({'task_id': task_id, 'library': 'long'})['ratio_dict']
    names = ['PC{}'.format(i) for i in range(1, len(ratio_dict) + 1)]
    df = pd.DataFrame([(name, ratio_dict[name]) for name in names], columns=['Name', 'Proportion of Variance'])
    df.to_csv(os.path.join(output_dir, 'explained_variance_ratio.xls'), sep='\t', index=False)


def set_diff_express(map_dict, task_id, output_dir, main_id=None):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    samples = [document['old_name'] for document in database['specimen'].find({'task_id': task_id, 'library': 'long'})]
    way = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['way']
    t_count_df = pd.read_table(map_dict['t_count'], index_col='transcript_id').reindex(samples, axis=1)
    t_count_df = t_count_df.rename({sample: '{}_count'.format(sample) for sample in samples}, axis=1)
    t_exp_df = pd.read_table(map_dict['t_exp'], index_col='transcript_id').reindex(samples, axis=1)
    t_exp_df = t_exp_df.rename({sample: '{}_{}'.format(sample, way) for sample in samples}, axis=1)
    exp_id = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['main_id']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'lncRNA'})
    t_anno_df = pd.DataFrame(list(cursor))
    t_anno_df = t_anno_df.set_index('transcript_id').reindex(['gene_id', 'gene_name', 'description'], axis=1)
    if main_id:
        main_dict = database['diff'].find_one({'main_id': ObjectId(main_id)})
    else:
        main_dict = database['diff'].find_one({'task_id': task_id, 'category': 'lncRNA'})
    diff_id, params, cmp_list = main_dict['main_id'], main_dict['params'], main_dict['cmp_list']
    diff_method = json.loads(params)['diff_method'].lower()
    total_df = t_anno_df.copy()
    total_columns = list()
    df_dict = dict()
    for compare in cmp_list:
        ctrl, case = compare.split('|')
        vs_pair = '{}_vs_{}'.format(ctrl, case)
        div_pair = '{}/{}'.format(case, ctrl)
        cursor = database['diff_detail'].find({'diff_id': diff_id, 'compare': compare})
        diff_df = pd.DataFrame(list(cursor)).set_index('seq_id')
        diff_df.index.name = 'transcript_id'
        diff_df = diff_df.join(t_count_df)
        diff_df = diff_df.join(t_exp_df)
        columns = list(t_count_df.columns) + list(t_exp_df.columns) + ['group1', 'group2', 'fc', 'log2fc', 'pvalue',
                                                                       'padjust', 'significant', 'regulate']
        diff_df = diff_df.reindex(columns, axis=1)
        diff_df.rename({'group1': '{}_mean_{}'.format(ctrl, way),
                        'group2': '{}_mean_{}'.format(case, way)}, axis=1, inplace=True)
        diff_df = diff_df.join(t_anno_df)
        diff_df.to_csv(os.path.join(output_dir, '{}_vs_{}_{}.xls'.format(ctrl, case, diff_method)), sep='\t')
        columns = ['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
        diff_df = diff_df.reindex(columns, axis=1)
        mapper = {'fc': '{}_fc({})'.format(vs_pair, div_pair), 'log2fc': '{}_log2fc({})'.format(vs_pair, div_pair),
                  'pvalue': '{}_pvalue'.format(vs_pair), 'padjust': '{}_padjust'.format(vs_pair),
                  'significant': '{}_significant'.format(vs_pair), 'regulate': '{}_regulate'.format(vs_pair)}
        diff_df = diff_df.rename(mapper, axis=1)
        total_df = total_df.join(diff_df)
        total_columns.extend(diff_df.columns)
    anno_columns = list(t_anno_df.columns)
    total_columns = anno_columns + total_columns
    total_df = total_df.reindex(total_columns, axis=1)
    total_df.to_csv(os.path.join(output_dir, 'total_diff_stat.{}.xls'.format(diff_method)), sep='\t')
    cursor = database['diff_summary'].find({'diff_id': diff_id})
    summary = list(cursor)
    if len(summary) > 0:
        summary_df = pd.DataFrame(summary).set_index('seq_id')
        summary_df.index.name = 'transcript_id'
        vs_pairs = ['{}_vs_{}'.format(*compare.split('|')) for compare in cmp_list]
        columns = ['gene_id', 'gene_name', 'description'] + vs_pairs + ['sum']
        summary_df = summary_df.join(t_anno_df)
        summary_df = summary_df.reindex(columns, axis=1)
        summary_df.to_csv(os.path.join(output_dir, 'diff_summary_{}.xls'.format(diff_method)), sep='\t')
    else:
        output = os.path.join(output_dir, 'diff_summary_{}.xls'.format(diff_method))
        cmd = 'touch {}'.format(output)
        os.system(cmd)
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_lncrna_target(map_dict, task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    with open(os.path.join(output_dir, 'tar_predict_detail.xls'), 'w') as f:
        f.write(
            "lncRNA id\tlncRNA gene id\tType\tTarget gene id\tTarget gene name\tGene description\tRegulate function\tcorr\tpvalue\tpadjust\tLocation\tDistance\tlncRNA position\tTarget gene position\n")
        if os.path.exists(os.path.join(map_dict['cistrans_dir'], 'cis_annot.xls')):
            with open(os.path.join(map_dict['cistrans_dir'], 'cis_annot.xls'), 'r') as f_cis:
                f_cis.readline()
                for line in f_cis:
                    cols = line.strip("\n").split("\t")
                    f.write("\t".join(cols[0:2] + [cols[4], cols[2], cols[5], cols[-1], cols[3]] + cols[-8:-1]) + "\n")
        if os.path.exists(os.path.join(map_dict['cistrans_dir'], 'trans_annot.xls')):
            with open(os.path.join(map_dict['cistrans_dir'], 'trans_annot.xls'), 'r') as f_trans:
                f_trans.readline()
                for line in f_trans:
                    cols = line.strip("\n").split("\t")
                    f.write("\t".join(cols[0:2] + [cols[4], cols[2], cols[5], cols[-1], cols[3]] + cols[-4:-1]) + "\n")
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_lncrna_structure(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_lncrna_famliy(task_id, os.path.join(output_dir, '01_LncRNA_Family'))
    set_mirna_pre(task_id, os.path.join(output_dir, '02_miRNA_Pre'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_lncrna_famliy(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    lnc_fam_id = database['lncrna_family'].find_one({'task_id': task_id})['main_id']
    cursor = database['lncrna_family_detail'].find({'lnc_fam_id': lnc_fam_id})
    df = pd.DataFrame(list(cursor))
    columns = ['lncrna_id', 'lncrna_gene_id', 'family_name', 'family_id', 'lncrna_start', 'lncrna_end', 'e_value',
               'score']
    df = df.reindex(columns, axis=1)
    df.to_csv(os.path.join(output_dir, 'lncRNA_family.xls'), sep='\t', index=False)


def set_mirna_pre(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    mir_pre_id = database['mirna_precursor'].find_one({'task_id': task_id})['main_id']
    cursor = database['mirna_precursor_detail'].find({'mir_pre_id': mir_pre_id})
    df = pd.DataFrame(list(cursor))
    columns = ['lncrna_id', 'lncrna_gene_id', 'pre_mirna_name', 'identity', 'pre_mirna_length', 'alignment_length',
               'e_value', 'score', 'lncrna_start', 'lncrna_end', 'pre_mirna_start', 'pre_mirna_end', 'mismatch']
    df = df.reindex(columns, axis=1)
    df.to_csv(os.path.join(output_dir, 'miRNA_precursor.xls'), sep='\t', index=False)


if __name__ == '__main__':
    import unittest


    class Empty:
        pass


    class TestFunction(unittest.TestCase):
        '''
        This is test for the package. Just run this script to do test.
        '''

        def test(self):
            from mbio.packages.whole_transcriptome.catalogue import lncrna
            self.task_id = 'tsg_36238'
            self.work_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191122/WholeTranscriptome_tsg_36238'
            self.output_dir = os.path.join(self.work_dir, 'output')
            self.tools = {'transfer_l': Empty()}
            self.tools['transfer_l'].output_dir = os.path.join(self.work_dir, 'Transfer/output')
            self.modules = {'target_cistrans': Empty()}
            self.modules['target_cistrans'].output_dir = os.path.join(self.work_dir, 'TargetCistrans/output')
            self.long_task_info = database['task'].find_one({'task_id': 'tsg_36117'})
            lncrna_dir = os.path.join(self.output_dir, 'lncrna')
            if os.path.isdir(lncrna_dir):
                shutil.rmtree(lncrna_dir)
            os.mkdir(lncrna_dir)

            map_dict = {'large_gush_dir': os.path.join(self.tools['transfer_l'].output_dir, 'large_gush')}
            lncrna.set_lncrna_analysis(map_dict, self.task_id, os.path.join(lncrna_dir, '01_lncRNA_Analysis'))
            map_dict = {'t_count': os.path.join(self.tools['transfer_l'].output_dir, 'exp_make/count/T.reads.txt')}
            lncrna.set_express(map_dict, self.task_id, os.path.join(lncrna_dir, '02_Express'))
            way = self.long_task_info['options']['exp_way']
            map_dict = {'t_count': os.path.join(lncrna_dir, '02_Express/01_Exp_Annalysis/lncRNA_count.xls'),
                        't_exp': os.path.join(lncrna_dir, '02_Express/01_Exp_Annalysis/lncRNA_{}_anno.xls'.format(way))}
            lncrna.set_diff_express(map_dict, self.task_id, os.path.join(lncrna_dir, '03_Diff_Express'))
            map_dict = {'cistrans_dir': self.modules['target_cistrans'].output_dir}
            lncrna.set_lncrna_target(map_dict, self.task_id, os.path.join(lncrna_dir, '04_LncRNA_Target'))
            lncrna.set_lncrna_structure(self.task_id, os.path.join(lncrna_dir, '05_LncRNA_Structure'))


    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
