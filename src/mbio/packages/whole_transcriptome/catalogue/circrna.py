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

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('circRNA')
logger.setLevel(logging.DEBUG)

database = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]

def set_background(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    species_id = database['species_information'].find_one({'task_id': task_id})['main_id']
    cursor = database['species_information_detail'].find({'species_id': species_id})
    genome_df = pd.DataFrame(list(cursor))
    genome_df = genome_df.reindex(
        ['chromosome', 'length', 'gene', 'qc_percent', 'proteincoding', 'pseudogene', 'other_rna'], axis=1)
    genome_df = genome_df.rename(
        {'chromosome': 'Chromosome', 'length': 'Length(Mb)', 'gene': 'Gene', 'qc_percent': 'GC',
         'proteincoding': 'Protein coding', 'pseudogene': 'Pseudogene', 'other_rna': 'Other RNA'}, axis=1)
    genome_df.to_csv(os.path.join(output_dir, 'genome_info.xls'), sep='\t', index=False)
    cursor = database['specimen'].find({'task_id': task_id, 'library': 'circle'})
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
    qc_id = database['qc'].find_one({'task_id': task_id, 'library': 'circle'})['main_id']
    cursor = database['qc_detail'].find({'qc_id': qc_id})
    df = pd.DataFrame(list(cursor))
    columns = ['sample', 'raw_reads', 'raw_bases', 'clean_reads', 'clean_bases', 'error_rate', 'q20_rate', 'q30_rate',
               'gc_rate']
    df = df.reindex(columns, axis=1)
    df = df.rename(
        {'sample': 'Sample', 'raw_reads': 'Raw reads', 'raw_bases': 'Raw bases', 'clean_reads': 'Clean reads',
         'clean_bases': 'Clean bases', 'error_rate': 'Error rate(%)', 'q20_rate': 'Q20(%)', 'q30_rate': 'Q30(%)',
         'gc_rate': 'GC content(%)'}, axis=1)
    df.to_csv(os.path.join(output_dir, 'QC_stat.xls'), sep='\t', index=False)


def set_align(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    mapping_id = database['mapping'].find_one({'task_id': task_id, 'library': 'circle'})['main_id']
    cursor = database['mapping_detail'].find({'mapping_id': mapping_id})
    df = pd.DataFrame(list(cursor))
    columns = ['specimen_name', 'total_reads', 'mapping_reads', 'multiple_mapped', 'uniq_mapped']
    df = df.reindex(columns, axis=1)
    df.to_csv(os.path.join(output_dir, 'align_stat.xls'), sep='\t', index=False)
    set_quality_assessment(task_id, os.path.join(output_dir, 'Quality_Assessment'))


def set_quality_assessment(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    distribution_id = database['assessment_distribution'].find_one({'task_id': task_id, 'library': 'circle'})['main_id']
    cursor = database['assessment_distribution_detail'].find({'distribution_id': distribution_id})
    rd_df = pd.DataFrame(list(cursor))
    columns = ['specimen_name', 'cds', 'utr5', 'utr3', 'introns', 'intergenic']
    rd_df = rd_df.reindex(columns, axis=1)
    rd_df.to_csv(os.path.join(output_dir, 'region_distribution.xls'), sep='\t', index=False)
    main_dict = database['assessment_chrom_distribution'].find_one({'task_id': task_id, 'library': 'circle'})
    row_index = main_dict['distribution']
    col_index = main_dict['specimen_name']
    chrom_distribution_id = main_dict['main_id']
    cursor = database['assessment_chrom_distribution_detail'].find({'chrom_distribution_id': chrom_distribution_id})
    data = dict()
    for document in cursor:
        sample = document['specimen_name']
        data[sample] = {chr_value['chr_name']: chr_value['value'] for chr_value in document['chr_values']}
    cd_df = pd.DataFrame(data)
    cd_df = cd_df.reindex(row_index, axis=0)
    cd_df = cd_df.reindex(col_index, axis=1)
    cd_df.index.name = 'chromosome'
    raw_indices = cd_df.index.tolist()
    new_indices = list()
    for i in raw_indices:
        try:
            j = int(i)
        except:
            j = i
        new_indices.append(j)
    new_indices.sort()
    cd_df = cd_df.reindex(map(str, new_indices))
    cd_df = cd_df.fillna(0)
    cd_df.to_csv(os.path.join(output_dir, 'chr_distribution.xls'), sep='\t')


def set_circrna_predict(map_dict, task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    os.link(map_dict['fasta'], os.path.join(output_dir, 'circRNA.fa'))
    main_dict = database['circrna_identify'].find_one({'task_id': task_id})
    detail_id, circbase = main_dict['main_id'], main_dict['circbase']
    cursor = database['circrna_identify_detail'].find({'detail_id': detail_id})
    deta_df = pd.DataFrame(list(cursor))
    columns = ['circrna_id', 'host_gene_id', 'chr', 'strand', 'circrna_start', 'circrna_end', 'signal', 'circrna_type']
    mapper = {'circrna_id': 'circRNA id', 'host_gene_id': 'Host gene id', 'chr': 'Chromosome', 'strand': 'Strand',
              'circrna_start': 'circRNA start', 'circrna_end': 'circRNA end', 'signal': 'Signal',
              'circrna_type': 'circRNA type'}
    if circbase:
        columns.extend(['circbase', 'database_id'])
        mapper.update({'circbase': 'circBase', 'database_id': 'Database id'})
    deta_df = deta_df.reindex(columns, axis=1)
    deta_df.rename(mapper, axis=1, inplace=True)
    deta_df.to_csv(os.path.join(output_dir, 'circRNA_predict_detail.xls'), sep='\t', index=False)
    stat_df = pd.DataFrame(list(deta_df['circRNA type'].value_counts().items()), columns=['Type', 'Number'])
    stat_df = stat_df.append({'Type': 'total', 'Number': stat_df['Number'].sum()}, ignore_index=True)
    stat_df.to_csv(os.path.join(output_dir, 'circRNA_stat.xls'), sep='\t', index=False)
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
    library = database['task'].find_one({'task_id': task_id})['rna']['circRNA']
    samples = [document['old_name'] for document in database['specimen'].find({'task_id': task_id, 'library': library})]
    main_dict = database['exp'].find_one({'task_id': task_id, 'level': 'T'})
    exp_id, way = main_dict['main_id'], main_dict['way']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'circRNA'})
    t_detail_df = pd.DataFrame(list(cursor)).set_index('transcript_id')
    t_detail_df = t_detail_df.reindex(['gene_id'] + samples, axis=1)
    t_detail_df['gene_id'] = t_detail_df['gene_id'].fillna(str())
    c_exp_df = t_detail_df
    c_exp_df.index.name = 'circRNA_id'
    c_exp_df.rename({'gene_id': 'host_gene_id'}, axis=1, inplace=True)
    c_exp_df.to_csv(os.path.join(output_dir, 'circRNA_rpm.xls'), sep='\t')
    row_index = c_exp_df.index
    col_index = c_exp_df.columns
    gene_id_df = c_exp_df.reindex(['gene_id'], axis=1)
    c_count_df = pd.read_table(map_dict['c_count'], index_col=0)
    c_count_df.index.name = 'circRNA_id'
    c_count_df = c_count_df.reindex(row_index)
    c_count_df = c_count_df.join(gene_id_df).reindex(col_index, axis=1)
    c_count_df.rename({'gene_id': 'host_gene_id'}, axis=1, inplace=True)
    c_count_df.to_csv(os.path.join(output_dir, 'circRNA_count.xls'), sep='\t')


def set_exp_corr(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    library = database['task'].find_one({'task_id': task_id})['rna']['circRNA']
    main_dict = database['exp_corr'].find_one({'task_id': task_id, 'library': library})
    corr_id, samples = main_dict['main_id'], main_dict['samples']
    cursor = database['exp_corr_detail'].find({'corr_id': corr_id})
    df = pd.DataFrame(list(cursor)).set_index('sample')
    df = df.reindex(samples, axis=1)
    df.to_csv(os.path.join(output_dir, 'sample_correlation.xls'), sep='\t')


def set_exp_pca(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    library = database['task'].find_one({'task_id': task_id})['rna']['circRNA']
    ratio_dict = database['exp_pca'].find_one({'task_id': task_id, 'library': library})['ratio_dict']
    names = ['PC{}'.format(i) for i in range(1, len(ratio_dict) + 1)]
    df = pd.DataFrame([(name, ratio_dict[name]) for name in names], columns=['Name', 'Proportion of Variance'])
    df.to_csv(os.path.join(output_dir, 'explained_variance_ratio.xls'), sep='\t', index=False)


def set_diff_express(map_dict, task_id, output_dir, main_id=None):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    library = database['task'].find_one({'task_id': task_id})['rna']['circRNA']
    samples = [document['old_name'] for document in database['specimen'].find({'task_id': task_id, 'library': library})]
    way = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['way']
    c_count_df = pd.read_table(map_dict['c_count'], index_col=0).reindex(samples, axis=1)
    c_count_df.index.name = 'circRNA_id'
    c_count_df = c_count_df.rename({sample: '{}_count'.format(sample) for sample in samples}, axis=1)
    c_exp_df = pd.read_table(map_dict['c_exp'], index_col=0).reindex(samples, axis=1)
    c_exp_df.index.name = 'circRNA_id'
    c_exp_df = c_exp_df.rename({sample: '{}_rpm'.format(sample) for sample in samples}, axis=1)
    exp_id = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['main_id']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'circRNA'})
    anno_df = pd.DataFrame(list(cursor))
    anno_df = anno_df.reindex(['transcript_id', 'gene_id'], axis=1).drop_duplicates()
    anno_df.set_index('transcript_id', inplace=True)
    anno_df.index.name = 'circRNA_id'
    total_df = anno_df.copy()
    total_columns = ['gene_id']
    if main_id:
        main_dict = database['diff'].find_one({'main_id': ObjectId(main_id)})
    else:
        main_dict = database['diff'].find_one({'task_id': task_id, 'category': 'circRNA'})
    diff_id, params, cmp_list = main_dict['main_id'], main_dict['params'], main_dict['cmp_list']
    diff_method = json.loads(params)['diff_method'].lower()
    df_dict = dict()
    for compare in cmp_list:
        ctrl, case = compare.split('|')
        vs_pair = '{}_vs_{}'.format(ctrl, case)
        div_pair = '{}/{}'.format(case, ctrl)
        cursor = database['diff_detail'].find({'diff_id': diff_id, 'compare': compare})
        diff_df = pd.DataFrame(list(cursor)).set_index('seq_id')
        diff_df.index.name = 'circRNA_id'
        diff_df = diff_df.join(anno_df)
        diff_df = diff_df.join(c_count_df)
        diff_df = diff_df.join(c_exp_df)
        columns = ['gene_id'] + list(c_count_df.columns) + list(c_exp_df.columns) + [
            'group1', 'group2', 'fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate'
        ]
        diff_df = diff_df.reindex(columns, axis=1)
        diff_df = diff_df.rename({'gene_id': 'host_gene_id', 'group1': '{}_mean_{}'.format(ctrl, way),
                                  'group2': '{}_mean_{}'.format(case, way)}, axis=1)
        diff_df.to_csv(os.path.join(output_dir, '{}_vs_{}_{}.xls'.format(ctrl, case, diff_method)), sep='\t')
        columns = ['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
        diff_df = diff_df.reindex(columns, axis=1)
        mapper = {'fc': '{}_fc({})'.format(vs_pair, div_pair), 'log2fc': '{}_log2fc({})'.format(vs_pair, div_pair),
                  'pvalue': '{}_pvalue'.format(vs_pair), 'padjust': '{}_padjust'.format(vs_pair),
                  'significant': '{}_significant'.format(vs_pair), 'regulate': '{}_regulate'.format(vs_pair)}
        diff_df = diff_df.rename(mapper, axis=1)
        diff_df.reset_index(inplace=True)
        diff_df = diff_df.drop_duplicates()
        diff_df.set_index('circRNA_id', inplace=True)
        total_df = total_df.join(diff_df)
        total_columns.extend(diff_df.columns)
    total_df = total_df.reindex(total_columns, axis=1)
    total_df.rename({'gene_id': 'Host_gene_id'}, axis=1, inplace=True)
    total_df.to_csv(os.path.join(output_dir, 'total_diff_stat_{}.xls'.format(diff_method)), sep='\t')
    cursor = database['diff_summary'].find({'diff_id': diff_id})
    summary_df = pd.DataFrame(list(cursor))
    if len(summary_df):
        summary_df = summary_df.set_index('seq_id')
        summary_df.index.name = 'circRNA_id'
        vs_pairs = ['{}_vs_{}'.format(*compare.split('|')) for compare in cmp_list]
        columns = vs_pairs + ['sum']
        summary_df = summary_df.reindex(columns, axis=1)
        summary_df.to_csv(os.path.join(output_dir, 'diff_summary_{}.xls'.format(diff_method)), sep='\t')
    else:
        summary_df.to_csv(os.path.join(output_dir, 'diff_summary_{}.xls'.format(diff_method)), sep='\t', index=False)
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


if __name__ == '__main__':
    import unittest


    class Empty:
        pass


    class TestFunction(unittest.TestCase):
        '''
        This is test for the package. Just run this script to do test.
        '''

        def test(self):
            from mbio.packages.whole_transcriptome.catalogue import circrna
            self.task_id = 'tsg_36238'
            self.work_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191122/WholeTranscriptome_tsg_36238'
            self.output_dir = os.path.join(self.work_dir, 'output')
            self.tools = {'transfer_l': Empty()}
            self.tools['transfer_l'].output_dir = os.path.join(self.work_dir, 'Transfer/output')
            self.long_task_info = database['task'].find_one({'task_id': 'tsg_36117'})
            self.lib_dict = {'long': True, 'small': False, 'circle': False}

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


    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
