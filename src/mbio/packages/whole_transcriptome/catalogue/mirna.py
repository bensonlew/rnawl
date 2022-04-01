# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import csv
import glob
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
logger = logging.getLogger('miRNA')
logger.setLevel(logging.DEBUG)

database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]


def set_background(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    species_id = database['species_information'].find_one({'task_id': task_id})['main_id']
    cursor = database['species_information_detail'].find({'species_id': species_id})
    genome_df = pd.DataFrame(list(cursor))
    genome_df = genome_df.reindex(['chromosome', 'length', 'gene', 'proteincoding', 'pseudogene', 'other_rna'], axis=1)
    genome_df.to_csv(os.path.join(output_dir, 'genome_info.xls'), sep='\t', index=False)
    cursor = database['specimen'].find({'task_id': task_id, 'library': 'small'})
    sample_df = pd.DataFrame(list(cursor))
    if 'productive_name' in sample_df.columns:
        sample_df = sample_df.reindex(['productive_name', 'old_name', 'new_name', 'group', 'library'], axis=1)
        sample_df = sample_df.rename(
            {'productive_name': 'Sample productive name', 'old_name': 'Sample initial name', 'new_name': 'Sample analysis name', 'group': 'Group name',
             'library': 'Library type'}, axis=1)
    else:
        sample_df = sample_df.reindex(['old_name', 'new_name', 'group', 'library'], axis=1)
        sample_df = sample_df.rename(
            {'old_name': 'Sample initial name', 'new_name': 'Sample analysis name', 'group': 'Group name',
             'library': 'Library type'}, axis=1)
    sample_df['Sample description'] = str()
    sample_df.to_csv(os.path.join(output_dir, 'sample_info.xls'), sep='\t', index=False)
    cursor = database['software_database'].find({})
    software_df = pd.DataFrame(list(cursor))
    software_df = software_df.reindex(['software_database', 'version', 'usage', 'source'], axis=1)
    software_df.to_csv(os.path.join(output_dir, 'software_info.xls'), sep='\t', index=False)
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_basic_analysis(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_qc(task_id, os.path.join(output_dir, '01_QC'))
    set_align(task_id, os.path.join(output_dir, '02_Align'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_qc(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    qc_id = database['qc'].find_one({'task_id': task_id, 'library': 'small'})['main_id']
    cursor = database['qc_detail'].find({'qc_id': qc_id})
    df = pd.DataFrame(list(cursor))
    columns = ['sample', 'raw_reads', 'raw_bases', 'clean_reads', 'clean_bases', 'error_rate', 'q20_rate', 'q30_rate',
               'gc_rate', 'useful_reads']
    df = df.reindex(columns, axis=1)
    df = df.rename(
        {'sample': 'Sample', 'raw_reads': 'Raw reads', 'raw_bases': 'Raw bases', 'clean_reads': 'Clean reads',
         'clean_bases': 'Clean bases', 'error_rate': 'Error rate(%)', 'q20_rate': 'Q20(%)', 'q30_rate': 'Q30(%)',
         'gc_rate': 'GC content(%)', 'useful_reads': 'Useful reads(18-32nt)'}, axis=1)
    df.to_csv(os.path.join(output_dir, 'QC_stat.xls'), sep='\t', index=False)


def set_align(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    mapping_id = database['mapping'].find_one({'task_id': task_id, 'library': 'small'})['main_id']
    cursor = database['mapping_stat'].find({'mapping_id': mapping_id})
    stat_df = pd.DataFrame(list(cursor))
    columns = ['specimen_name', 'total_reads', 'mapped_reads', 'mapped_reads_for', 'mapped_reads_rec']
    stat_df = stat_df.reindex(columns, axis=1)
    stat_df.to_csv(os.path.join(output_dir, 'align_stat.xls'), sep='\t', index=False)
    samples = [document['old_name'] for document in database['specimen'].find({'task_id': task_id, 'library': 'small'})]
    cursor = database['mapping_detail'].find({'mapping_id': mapping_id})
    dist_df = pd.DataFrame(list(cursor))
    columns = ['ref'] + ['{}_total'.format(sample) for sample in samples]
    dist_df = dist_df.reindex(columns, axis=1)
    dist_df.rename({'ref': 'Chromosome'}, axis=1, inplace=True)
    dist_df.to_csv(os.path.join(output_dir, 'chr_distribution.xls'), sep='\t', index=False)


def set_srna_analysis(map_dict, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_known_mirna(map_dict, os.path.join(output_dir, '01_Known_miRNA'))
    set_novel_mirna(map_dict, os.path.join(output_dir, '02_Novel_miRNA'))
    set_srna_stat(map_dict, os.path.join(output_dir, '03_sRNA_stat'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_known_mirna(map_dict, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    out_structure_dir = os.path.join(output_dir, 'known_pre_structure')
    os.mkdir(out_structure_dir)
    for pdf in os.listdir(map_dict['ref_structure_dir']):
        os.link(os.path.join(map_dict['ref_structure_dir'], pdf), os.path.join(out_structure_dir, pdf))
    os.link(map_dict['ref_detail_table'], os.path.join(output_dir, 'known_miRNA_detail.xls'))
    os.link(map_dict['ref_hairpin_fasta'], os.path.join(output_dir, 'known_hairpin.fa'))
    os.link(map_dict['ref_mature_fasta'], os.path.join(output_dir, 'known_mature.fa'))


def set_novel_mirna(map_dict, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    out_structure_dir = os.path.join(output_dir, 'novel_pre_structure')
    os.mkdir(out_structure_dir)
    for pdf in os.listdir(map_dict['new_structure_dir']):
        os.link(os.path.join(map_dict['new_structure_dir'], pdf), os.path.join(out_structure_dir, pdf))
    os.link(map_dict['new_detail_table'], os.path.join(output_dir, 'novel_miRNA_detail.xls'))
    os.link(map_dict['new_hairpin_fasta'], os.path.join(output_dir, 'novel_hairpin.fa'))
    os.link(map_dict['new_mature_fasta'], os.path.join(output_dir, 'novel_mature.fa'))


def set_srna_stat(map_dict, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    os.link(map_dict['mirna_stat_table'], os.path.join(output_dir, 'miRNA_stat.xls'))
    os.link(map_dict['srna_stat_table'], os.path.join(output_dir, 'sRNA_stat.xls'))


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
    specimen_group_dict = database['specimen_group'].find_one({'task_id': task_id, 'library': 'small'})
    main_dict = database['exp'].find_one({'task_id': task_id, 'level': 'T'})
    exp_id, way = main_dict['main_id'], 'tpm'
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'miRNA'})
    t_detail_df = pd.DataFrame(list(cursor)).set_index('transcript_id')
    t_detail_df = t_detail_df.reindex(samples, axis=1)
    for group, specimens in zip(specimen_group_dict['category_names'], specimen_group_dict['specimen_names']):
        t_detail_df[group] = t_detail_df.reindex(specimens, axis=1).mean(axis=1)
    s_exp_df = t_detail_df
    s_exp_df.index.name = 'miRNA_name'
    s_exp_df.to_csv(os.path.join(output_dir, 'miRNA_{}.xls'.format(way)), sep='\t')
    row_index = s_exp_df.index
    s_count_df = pd.read_table(map_dict['s_count'], index_col=0)
    s_count_df.index.name = 'miRNA_name'
    s_count_df = s_count_df.reindex(row_index)
    s_count_df.to_csv(os.path.join(output_dir, 'miRNA_count.xls'), sep='\t')


def set_exp_corr(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    main_dict = database['exp_corr'].find_one({'task_id': task_id, 'library': 'small'})
    corr_id, samples = main_dict['main_id'], main_dict['samples']
    cursor = database['exp_corr_detail'].find({'corr_id': corr_id})
    df = pd.DataFrame(list(cursor)).set_index('sample')
    df = df.reindex(samples, axis=1)
    df.to_csv(os.path.join(output_dir, 'sample_correlation.xls'), sep='\t')


def set_exp_pca(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    ratio_dict = database['exp_pca'].find_one({'task_id': task_id, 'library': 'small'})['ratio_dict']
    names = ['PC{}'.format(i) for i in range(1, len(ratio_dict) + 1)]
    df = pd.DataFrame([(name, ratio_dict[name]) for name in names], columns=['Name', 'Proportion of Variance'])
    df.to_csv(os.path.join(output_dir, 'explained_variance_ratio.xls'), sep='\t', index=False)


def set_diff_express(map_dict, task_id, output_dir, main_id=None):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    samples = [document['old_name'] for document in database['specimen'].find({'task_id': task_id, 'library': 'small'})]
    way = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['way']
    s_count_df = pd.read_table(map_dict['s_count'], index_col=0).reindex(samples, axis=1)
    s_count_df.index.name = 'miRNA_id'
    s_count_df = s_count_df.rename({sample: '{}_count'.format(sample) for sample in samples}, axis=1)
    s_exp_df = pd.read_table(map_dict['s_exp'], index_col=0).reindex(samples, axis=1)
    s_exp_df.index.name = 'miRNA_id'
    s_exp_df = s_exp_df.rename({sample: '{}_{}'.format(sample, way) for sample in samples}, axis=1)
    exp_id = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['main_id']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'miRNA'})
    total_df = pd.DataFrame(list(cursor)).set_index('transcript_id').reindex(list(), axis=1)
    total_columns = list()
    if main_id:
        main_dict = database['diff'].find_one({'main_id': ObjectId(main_id)})
    else:
        main_dict = database['diff'].find_one({'task_id': task_id, 'category': 'miRNA'})
    diff_id, params, cmp_list = main_dict['main_id'], main_dict['params'], main_dict['cmp_list']
    diff_method = json.loads(params)['diff_method'].lower()
    df_dict = dict()
    for compare in cmp_list:
        ctrl, case = compare.split('|')
        vs_pair = '{}_vs_{}'.format(ctrl, case)
        div_pair = '{}/{}'.format(case, ctrl)
        cursor = database['diff_detail'].find({'diff_id': diff_id, 'compare': compare})
        diff_df = pd.DataFrame(list(cursor)).set_index('seq_id')
        diff_df.index.name = 'miRNA_id'
        diff_df = diff_df.join(s_count_df)
        diff_df = diff_df.join(s_exp_df)
        columns = list(s_count_df.columns) + list(s_exp_df.columns) + ['group1', 'group2', 'fc', 'log2fc', 'pvalue',
                                                                       'padjust', 'significant', 'regulate']
        diff_df = diff_df.reindex(columns, axis=1)
        diff_df.rename({'group1': '{}_mean_{}'.format(ctrl, way),
                        'group2': '{}_mean_{}'.format(case, way)}, axis=1, inplace=True)
        diff_df.to_csv(os.path.join(output_dir, '{}_vs_{}_{}.xls'.format(ctrl, case, diff_method)), sep='\t')
        columns = ['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
        diff_df = diff_df.reindex(columns, axis=1)
        mapper = {'fc': '{}_fc({})'.format(vs_pair, div_pair), 'log2fc': '{}_log2fc({})'.format(vs_pair, div_pair),
                  'pvalue': '{}_pvalue'.format(vs_pair), 'padjust': '{}_padjust'.format(vs_pair),
                  'significant': '{}_significant'.format(vs_pair), 'regulate': '{}_regulate'.format(vs_pair)}
        diff_df = diff_df.rename(mapper, axis=1)
        total_df = total_df.join(diff_df)
        total_columns.extend(diff_df.columns)
    total_df = total_df.reindex(total_columns, axis=1)
    total_df.to_csv(os.path.join(output_dir, 'total_diff_stat_{}.xls'.format(diff_method)), sep='\t')
    cursor = database['diff_summary'].find({'diff_id': diff_id})
    summary_df = pd.DataFrame(list(cursor))
    if len(summary_df) and 'seq_id' in summary_df:
        summary_df.set_index('seq_id', inplace=True)
        summary_df.index.name = 'miRNA_id'
        vs_pairs = ['{}_vs_{}'.format(*compare.split('|')) for compare in cmp_list]
        columns = vs_pairs + ['sum']
        summary_df = summary_df.reindex(columns, axis=1)
        summary_df.to_csv(os.path.join(output_dir, 'diff_summary_{}.xls'.format(diff_method)), sep='\t')
    else:
        open(os.path.join(output_dir, 'diff_summary_{}.xls'.format(diff_method)), 'w').close()
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_mirna_target(map_dict, arg_dict, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    def check_soft(params, soft):
        if 'm_' + soft in params and params['m_' + soft] == 'yes':
            return True
        elif 'c_' + soft in params and params['c_' + soft] == 'yes':
            return True
        elif 'l_' + soft in params and params['l_' + soft] == 'yes':
            return True
        else:
            return False

    columns = ['query', 'target', 'category', 'gene_id', 'gene_name', 'gene_description']
    if check_soft(arg_dict, 'miranda'):
        columns.extend(['start_miranda', 'end_miranda', 'score_miranda', 'energy_miranda'])
    if check_soft(arg_dict, 'targetscan'):
        columns.extend(['utr_start_targetscan', 'utr_end_targetscan', 'msa_start_targetscan', 'msa_end_targetscan',
                        'sead_match_targetscan'])
        columns.extend(['pct', 'context_score'])
    if check_soft(arg_dict, 'psrobot'):
        columns.extend(['start_psrobot', 'end_psrobot', 'score_psrobot'])
    if check_soft(arg_dict, 'targetfinder'):
        columns.extend(['start_targetfinder', 'end_targetfinder', 'score_targetfinder'])
    if check_soft(arg_dict, 'rnahybrid'):
        columns.extend(['start_rnahybrid', 'end_rnahybrid', 'energy_rnahybrid', 'pvalue_rnahybrid'])
    columns.extend(['mirtarbase'])

    columns_rename_dict = {'query': 'miRNA name',
                           'target': 'Target RNA id',
                           'category': 'Target RNA type',
                           'gene_id': 'target gene id',
                           'gene_name': 'target gene name',
                           'gene_description': 'target gene description'}
    with open(os.path.join(output_dir, 'tar_predict_detail.xls'), 'w') as fo:
        columns_show = [columns_rename_dict.get(x, x) for x in columns]
        fo.write('\t'.join(columns_show) + '\n')
        for sub_dir in ['m_known', 'm_novel', 'c_known', 'c_novel', 'l_known', 'l_novel']:
            for sub_file in ['known_target.xls', 'novol_target.xls']:
                category = {'m': 'mRNA', 'c': 'circRNA', 'l': 'lncRNA'}.get(sub_dir.split('_')[0])
                # target_kind = {'known': 'ref', 'novel': 'new'}.get(sub_dir.split('_')[1])
                target_kind = sub_dir.split('_')[1]
                # mi_kind = {'known': 'ref', 'novol': 'new'}.get(sub_file.split('_')[0])
                mi_kind = sub_file.split('_')[0]
                file_path = os.path.join(map_dict['target_dir'], sub_dir, sub_file)
                if os.path.exists(file_path):
                    with open(file_path) as in_handler:
                        for line_dic in csv.DictReader(in_handler, delimiter='\t'):
                            line = list()
                            for col in columns:
                                if col == 'category':
                                    line.append(category)
                                else:
                                    line.append(line_dic.get(col, ''))
                            fo.write('\t'.join(line) + '\n')
                else:
                    logger.info('fail to find {} '.format(os.path.join(sub_dir, file_path)))
    if not os.path.isdir(os.path.join(output_dir, 'target_seqs')):
        os.mkdir(os.path.join(output_dir, 'target_seqs'))
    if not os.path.isdir(os.path.join(output_dir, 'target_align_detail')):
        os.mkdir(os.path.join(output_dir, 'target_align_detail'))
    cp_flag_dict = dict()
    for sub_dir in ['c_known', 'c_novel', 'l_known', 'l_novel', 'm_known', 'm_novel']:
        fa_type = {'m': 'mRNA', 'c': 'circRNA', 'l': 'lncRNA'}.get(sub_dir.split('_')[0]) + '_' + \
                  {'known': 'known', 'novel': 'novel'}.get(sub_dir.split('_')[1])
        CopyFile().linkfile(os.path.join(map_dict['target_dir'], sub_dir, 'target.fa'),
                            os.path.join(output_dir, 'target_seqs', '{}_target.fa'.format(fa_type)),
                            mode="cp")
        if not os.path.exists(os.path.join(output_dir, 'target_seqs', '{}_target.fa'.format(fa_type))):
            cp_flag_dict[fa_type] = False
        elif os.path.getsize(os.path.join(output_dir, 'target_seqs', '{}_target.fa'.format(fa_type))) == 0:
            os.remove(os.path.join(output_dir, 'target_seqs', '{}_target.fa'.format(fa_type)))
            cp_flag_dict[fa_type] = False
        else:
            cp_flag_dict[fa_type] = True
        for align_file in glob.glob(os.path.join(map_dict['target_dir'], sub_dir, '*.gz')):
            base_name = os.path.basename(align_file)
            if cp_flag_dict[fa_type]:
                CopyFile().linkfile(align_file,
                                    os.path.join(output_dir, 'target_align_detail', fa_type + "_" + base_name),
                                    mode="cp")
                CopyFile().linkfile(align_file,
                                    os.path.join(output_dir, 'target_align_detail', fa_type + "_" + base_name),
                                    mode="cp")
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_mirna_structure(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_mirna_bias(task_id, os.path.join(output_dir, '01_miRNA_bias'))
    set_mirna_edit(task_id, os.path.join(output_dir, '02_miRNA_edit'))
    set_mirna_family(task_id, os.path.join(output_dir, '03_miRNA_family'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_mirna_bias(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    bias_id = database['bias'].find_one({'task_id': task_id})['main_id']
    cursor = database['bias_first'].find({'bias_id': bias_id, 'type': 'all'})
    first_df = pd.DataFrame(list(cursor))
    columns = ['length', 'a', 'g', 'c', 'u']
    first_df = first_df.reindex(columns, axis=1)
    first_df.to_csv(os.path.join(output_dir, 'first_bias_per.xls'), sep='\t', index=False)
    cursor = database['bias_location'].find({'bias_id': bias_id, 'type': 'all'})
    locat_df = pd.DataFrame(list(cursor))
    columns = ['location', 'a', 'g', 'c', 'u']
    locat_df = locat_df.reindex(columns, axis=1)
    locat_df.to_csv(os.path.join(output_dir, 'loc_bias.xls'), sep='\t', index=False)


def set_mirna_edit(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    samples = [document['old_name'] for document in database['specimen'].find({'task_id': task_id, 'library': 'small'})]
    edit_id = database['mirna_edit'].find_one({'task_id': task_id})['main_id']
    cursor = database['mirna_edit_detail'].find({'edit_id': edit_id})
    df = pd.DataFrame(list(cursor))
    columns = ['miRNA_name', 'pre_miRNA_name', 'W', 'E', 'PP', 'MP']
    for sample in samples:
        columns.extend([column.format(sample) for column in ('{}_PV', '{}_PA', '{}_MN', '{}_TN')])
    df = df.reindex(columns, axis=1)
    mapper = dict()
    for i in column:
        if i.endswith('_MN'):
            s = i[:-3]
            mapper[i] = 'mismatch_number({})'.format(s)
        if i.endswith('_TN'):
            s = i[:-3]
            mapper[i] = 'total_number({})'.format(s)
    df.rename(mapper, axis=1, inplace=True)
    df.to_csv(os.path.join(output_dir, 'edit_detail.xls'), sep='\t', index=False)


def set_mirna_family(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    family_id = database['family'].find_one({'task_id': task_id})['main_id']
    cursor = database['family_detail'].find({'family_id': family_id})
    lines = list()
    for document in cursor:
        lines.append('{}\t{}\t{}\n'.format(document['mirna'], ';'.join(document['mir_pre_list']), document['family']))
    open(os.path.join(output_dir, 'family_species.xls'), 'w').writelines(lines)


if __name__ == '__main__':
    import unittest


    class Empty:
        pass


    class TestFunction(unittest.TestCase):
        '''
        This is test for the package. Just run this script to do test.
        '''

        def test(self):
            from mbio.packages.whole_transcriptome.catalogue import mirna
            self.task_id = 'tsg_36238'
            self.work_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191122/WholeTranscriptome_tsg_36238'
            self.output_dir = os.path.join(self.work_dir, 'output')
            self.tools = {'transfer_l': Empty(), 'transfer_s': Empty()}
            self.tools['transfer_l'].output_dir = os.path.join(self.work_dir, 'Transfer/output')
            self.tools['transfer_s'].output_dir = os.path.join(self.work_dir, 'Transfer1/output')
            self.modules = {'target_mirna': Empty()}
            self.modules['target_mirna'].output_dir = os.path.join(self.work_dir, 'TargetMirna/output')
            self.target_mirna_arg_dict = {'c_min_support': 1, 'c_miranda': 'yes', 'c_miranda_energy': -20.0,
                                          'c_miranda_score': 160.0, 'c_miranda_strict': 'on', 'c_ps_robot_score': 2.5,
                                          'c_psrobot': 'no', 'c_rnahybird_energy': -20.0, 'c_rnahybird_num': 100,
                                          'c_rnahybird_pvalue': 0.01, 'c_rnahybrid': 'no', 'c_targetfinder': 'no',
                                          'c_targetfinder_score': 4.0, 'c_targetscan': 'no', 'l_min_support': 1,
                                          'l_miranda': 'yes', 'l_miranda_energy': -20.0, 'l_miranda_score': 160.0,
                                          'l_miranda_strict': 'on', 'l_ps_robot_score': 2.5, 'l_psrobot': 'no',
                                          'l_rnahybird_energy': -20.0, 'l_rnahybird_num': 100,
                                          'l_rnahybird_pvalue': 0.01, 'l_rnahybrid': 'no', 'l_targetfinder': 'no',
                                          'l_targetfinder_score': 4.0, 'l_targetscan': 'no', 'm_min_support': 1,
                                          'm_miranda': 'yes', 'm_miranda_energy': -20.0, 'm_miranda_score': 160.0,
                                          'm_miranda_strict': 'on', 'm_ps_robot_score': 2.5, 'm_psrobot': 'no',
                                          'm_rnahybird_energy': -20.0, 'm_rnahybird_num': 100,
                                          'm_rnahybird_pvalue': 0.01, 'm_rnahybrid': 'no', 'm_targetfinder': 'no',
                                          'm_targetfinder_score': 4.0, 'm_targetscan': 'no',
                                          'submit_location': 'mi_target', 'task_id': 'tsg_36238', 'task_type': 2}
            self.long_task_info = database['task'].find_one({'task_id': 'tsg_36023'})
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
            way = self.long_task_info['options']['exp_way']
            map_dict = {'s_count': os.path.join(mirna_dir, '04_Express/01_Exp_Annalysis/miRNA_count.xls'),
                        's_exp': os.path.join(mirna_dir, '04_Express/01_Exp_Annalysis/miRNA_{}.xls'.format(way))}
            mirna.set_diff_express(map_dict, self.task_id, os.path.join(mirna_dir, '05_Diff_Express'))
            map_dict = {'target_dir': self.modules['target_mirna'].output_dir}
            arg_dict = self.target_mirna_arg_dict
            mirna.set_mirna_target(map_dict, arg_dict, os.path.join(mirna_dir, '06_miRNA_Target'))
            mirna.set_mirna_structure(self.task_id, os.path.join(mirna_dir, '07_miRNA_Structure'))


    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
