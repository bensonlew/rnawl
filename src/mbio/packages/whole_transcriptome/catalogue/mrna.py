# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import os
import re
import shutil
import sys

import pandas as pd
from biocluster.config import Config
from mbio.packages.project_demo.interaction_rerun.interaction_delete import linkfile,linkdir
from bson.objectid import ObjectId
import glob
logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('mRNA')
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
    cursor = database['specimen'].find({'task_id': task_id, 'library': 'long'})
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


def set_basic_analysis(map_dict, task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_qc(task_id, os.path.join(output_dir, '01_QC'))
    set_align(task_id, os.path.join(output_dir, '02_Align'))
    set_assemble(map_dict, task_id, os.path.join(output_dir, '03_Assemble'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_qc(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    qc_id = database['qc'].find_one({'task_id': task_id, 'library': 'long'})['main_id']
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
    mapping_id = database['mapping'].find_one({'task_id': task_id, 'library': 'long'})['main_id']
    cursor = database['mapping_detail'].find({'mapping_id': mapping_id})
    df = pd.DataFrame(list(cursor))
    columns = ['specimen_name', 'total_reads', 'mapping_reads', 'multiple_mapped', 'uniq_mapped']
    df = df.reindex(columns, axis=1)
    df.to_csv(os.path.join(output_dir, 'align_stat.xls'), sep='\t', index=False)
    set_quality_assessment(task_id, os.path.join(output_dir, 'Quality_Assessment'))


def set_quality_assessment(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    distribution_id = database['assessment_distribution'].find_one({'task_id': task_id, 'library': 'long'})['main_id']
    cursor = database['assessment_distribution_detail'].find({'distribution_id': distribution_id})
    rd_df = pd.DataFrame(list(cursor))
    columns = ['specimen_name', 'cds', 'utr5', 'utr3', 'introns', 'intergenic']
    rd_df = rd_df.reindex(columns, axis=1)
    rd_df.to_csv(os.path.join(output_dir, 'region_distribution.xls'), sep='\t', index=False)
    main_dict = database['assessment_chrom_distribution'].find_one({'task_id': task_id, 'library': 'long'})
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


def set_assemble(map_dict, task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_sequence(map_dict, os.path.join(output_dir, 'Sequence'))
    assembly_id = database['assembly'].find_one({'task_id': task_id})['main_id']
    cursor = database['assembly_code'].find({'assembly_id': assembly_id})
    desc_df = pd.DataFrame([{'class_code': '=', 'description': 'Complete match of intron chain'},
                            {'class_code': 'c', 'description': 'Contained'},
                            {'class_code': 'j',
                             'description': 'Potentially novel isoform (fragment):at least one splice junction is shared with a reference transcript'},
                            {'class_code': 'e',
                             'description': 'Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron,indiction a possible pre-mRNA fragment'},
                            {'class_code': 'i', 'description': 'Transfrag falling entirely within a reference intron'},
                            {'class_code': 'o', 'description': 'Generic exonic overlap with a referfence transcript'},
                            {'class_code': 'p',
                             'description': 'Possible polymerase run-on fragment(within 2Kbases of a reference transcript'},
                            {'class_code': 'r',
                             'description': 'Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case'},
                            {'class_code': 'u', 'description': 'Unknown,intergenic transcript'},
                            {'class_code': 'x', 'description': 'Exonic overlap with reference on the opposite strand'},
                            {'class_code': 's',
                             'description': 'An intron of the transfrag overlaps a reference intron on the opposite strand(likely due to read mapping errors'},
                            {'class_code': '.',
                             'description': 'Tracking file only,indicates multiple classifications'}])
    code_df = pd.DataFrame(list(cursor))
    code_df = pd.merge(code_df, desc_df)
    code_df = code_df.reindex(['class_code', 'description', 'num'], axis=1)
    code_df.to_csv(os.path.join(output_dir, 'classcode_stat.xls'), sep='\t', index=False)
    cursor = database['assembly_step'].find({'assembly_id': assembly_id})

    def get_row_index(step):
        row_index = list()
        for start, end in zip(range(1, step * 9 + 1, step), range(step, step * 10, step)):
            row_index.append('{}~{}'.format(start, end))
        row_index.extend(['>{}'.format(step * 9), 'total'])
        return row_index

    for document in cursor:
        step = document['step']
        data = {step: dict()}
        for step_dict in document['step_data']:
            data[step].update(step_dict)
        step_df = pd.DataFrame(data)
        row_index = get_row_index(step)
        step_df = step_df.reindex(row_index)
        step_df = step_df.rename({step: 'Number'}, axis=1)
        step_df.index.name = 'Length'
        step_df.to_csv(os.path.join(output_dir, 'length_distribution_{}.xls'.format(step)), sep='\t')


def set_sequence(map_dict, output_dir):
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)
    os.link(map_dict['new_gtf'], os.path.join(output_dir, 'new_transcript.gtf'))
    os.link(map_dict['all_gtf'], os.path.join(output_dir, 'all_transcript.gtf'))
    os.link(map_dict['new_fasta'], os.path.join(output_dir, 'new_transcript.fa'))
    os.link(map_dict['all_fasta'], os.path.join(output_dir, 'all_transcript.fa'))
    # 序列文件
    os.system("cat {} {} > {}".format(
        map_dict['ref_pep'],
        map_dict['new_pep'],
        os.path.join(output_dir, 'all_pep.fa')
    ))
    os.system("cat {} {} > {}".format(
        map_dict['ref_cds'],
        map_dict['new_cds'],
        os.path.join(output_dir, 'all_cds.fa')
    ))
    with open(os.path.join(output_dir, 'all_id.xls'), 'w') as fo, \
         open(map_dict['all_id'], 'r') as fin:
        fo.write("gene_id\ttranscript_id\tprotein_id\n")
        for line in fin:
            cols = line.strip("\n").split("\t")
            fo.write("{}\t{}\t{}\n".format(cols[1], cols[0], cols[4]))



def set_annotation(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    set_anno_detail(task_id, os.path.join(output_dir, '01_Anno_Detail'))
    set_anno_statistics(task_id, os.path.join(output_dir, '02_Anno_Statistics'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_anno_detail(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    def parse_nr(nr_series):
        p = re.compile(r'(\S+)\((.*?)\)')
        nr_hit_names, nr_descriptions = list(), list()
        for nr in nr_series:
            m = re.match(p, nr)
            if m:
                nr_hit_name, nr_description = m.groups()
            else:
                nr_hit_name, nr_description = str(), str()
            nr_hit_names.append(nr_hit_name)
            nr_descriptions.append(nr_description)
        return nr_hit_names, nr_descriptions

    def parse_go(go_series):
        p = re.compile(r'(\S+)\((\S+):(.*?)\)')
        go_ids, go_types, go_descriptions = list(), list(), list()
        for go in go_series:
            _go_ids = list()
            _go_types = list()
            _go_descriptions = list()
            for s in go.split(';'):
                m = re.search(p, s)
                if m:
                    go_id, go_type, go_description = m.groups()
                else:
                    go_id, go_type, go_description = str(), str(), str()
                _go_ids.append(go_id)
                _go_types.append(go_type)
                _go_descriptions.append(go_description)
            go_ids.append(';'.join(_go_ids))
            go_types.append(';'.join(_go_types))
            go_descriptions.append(';'.join(_go_descriptions))
        return go_ids, go_types, go_descriptions

    def parse_paths(paths_series):
        p = re.compile(r'(\S+)\((.*?)\)')
        pathway_ids, pathway_descriptions = list(), list()
        for paths in paths_series:
            m = re.match(p, paths)
            if m:
                pathway_id, pathway_description = m.groups()
            else:
                pathway_id, pathway_description = str(), str()
            pathway_ids.append(pathway_id)
            pathway_descriptions.append(pathway_description)
        return pathway_ids, pathway_descriptions

    def parse_cog(cog_series):
        p = re.compile(r'(\S+)\((.*?)\)')
        cog_ids, cog_functional_categories = list(), list()
        for cog in cog_series:
            m = re.match(p, cog)
            if m:
                cog_id, cog_functional_category = m.groups()
            else:
                cog_id, cog_functional_category = str(), str()
            cog_ids.append(cog_id)
            cog_functional_categories.append(cog_functional_category)
        return cog_ids, cog_functional_categories

    def parse_swissprot(swissprot_series):
        p = re.compile(r'(\S+?)\((.*?)\)')
        swissprot_hit_names, swissprot_descriptions = list(), list()
        for swissprot in swissprot_series:
            m = re.match(p, swissprot)
            if m:
                swissprot_hit_name, swissprot_description = m.groups()
            else:
                swissprot_hit_name, swissprot_description = str(), str()
            swissprot_hit_names.append(swissprot_hit_name)
            swissprot_descriptions.append(swissprot_description)
        return swissprot_hit_names, swissprot_descriptions

    def parse_pfam(pfam_series):
        p = re.compile(r'(\S+)\((\S+):(.*?)\)')
        pfam_ids, pfam_domains, pfam_descriptions = list(), list(), list()
        for pfam in pfam_series:
            _pfam_ids = list()
            _pfam_domains = list()
            _pfam_descriptions = list()
            for s in pfam.split(';'):
                m = re.search(p, s)
                if m:
                    pfam_id, pfam_domain, pfam_description = m.groups()
                else:
                    pfam_id, pfam_domain, pfam_description = str(), str(), str()
                _pfam_ids.append(pfam_id)
                _pfam_domains.append(pfam_domain)
                _pfam_descriptions.append(pfam_description)
            pfam_ids.append(';'.join(_pfam_ids))
            pfam_domains.append(';'.join(_pfam_domains))
            pfam_descriptions.append(';'.join(_pfam_descriptions))
        return pfam_ids, pfam_domains, pfam_descriptions

    exp_id = database['exp'].find_one({'task_id': task_id, 'level': 'G'})['main_id']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'mRNA'})
    g_df = pd.DataFrame(list(cursor))
    columns = ['gene_id', 'gene_name', 'description', 'nr', 'go', 'ko_id', 'ko_name', 'paths', 'cog', 'cog_description',
               'swissprot', 'pfam', 'kind']
    g_df = g_df.reindex(columns, axis=1)
    nr_hit_names, nr_descriptions = parse_nr(g_df['nr'])
    go_ids, go_types, go_descriptions = parse_go(g_df['go'])
    pathway_ids, pathway_descriptions = parse_paths(g_df['paths'])
    cog_ids, cog_functional_categories = parse_cog(g_df['cog'])
    swissprot_hit_names, swissprot_descriptions = parse_swissprot(g_df['swissprot'])
    pfam_ids, pfam_domains, pfam_descriptions = parse_pfam(g_df['pfam'])
    g_df['nr_hit_name'] = nr_hit_names
    g_df['nr_description'] = nr_descriptions
    g_df['go_id'] = go_ids
    g_df['go_type'] = go_types
    g_df['go_term'] = go_descriptions
    g_df['pathway_id'] = pathway_ids
    g_df['pathway_description'] = pathway_descriptions
    g_df['cog_id'] = cog_ids
    g_df['cog_functional_category'] = cog_functional_categories
    g_df['swissprot_hit_name'] = swissprot_hit_names
    g_df['swissprot_description'] = swissprot_descriptions
    g_df['pfam_id'] = pfam_ids
    g_df['pfam_domain'] = pfam_domains
    g_df['pfam_description'] = pfam_descriptions
    columns = ['gene_id', 'gene_name', 'description', 'nr_hit_name', 'nr_description', 'go_id', 'go_type', 'go_term',
               'ko_id', 'ko_name', 'pathway_id', 'pathway_description', 'cog_id', 'cog_functional_category',
               'swissprot_hit_name', 'swissprot_description', 'pfam_id', 'pfam_domain', 'pfam_description', 'kind']
    g_df = g_df.reindex(columns, axis=1)
    ref_g_df = g_df[g_df['kind'] == 'ref'].drop('kind', axis=1)
    new_g_df = g_df[g_df['kind'] == 'new'].drop('kind', axis=1)
    all_g_df = g_df.drop('kind', axis=1)
    ref_g_df.to_csv(os.path.join(output_dir, 'ref_gene_anno_detail.xls'), sep='\t', index=False)
    new_g_df.to_csv(os.path.join(output_dir, 'new_gene_anno_detail.xls'), sep='\t', index=False)
    all_g_df.to_csv(os.path.join(output_dir, 'all_gene_anno_detail.xls'), sep='\t', index=False)

    exp_id = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['main_id']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'mRNA'})
    t_df = pd.DataFrame(list(cursor))
    columns = ['transcript_id', 'gene_id', 'gene_name', 'description', 'nr', 'go', 'ko_id', 'ko_name', 'paths', 'cog',
               'cog_description', 'swissprot', 'pfam', 'kind']
    t_df = t_df.reindex(columns, axis=1)
    nr_hit_names, nr_descriptions = parse_nr(t_df['nr'])
    go_ids, go_types, go_descriptions = parse_go(t_df['go'])
    pathway_ids, pathway_descriptions = parse_paths(t_df['paths'])
    cog_ids, cog_functional_categories = parse_cog(t_df['cog'])
    swissprot_hit_names, swissprot_descriptions = parse_swissprot(t_df['swissprot'])
    pfam_ids, pfam_domains, pfam_descriptions = parse_pfam(t_df['pfam'])
    t_df['nr_hit_name'] = nr_hit_names
    t_df['nr_description'] = nr_descriptions
    t_df['go_id'] = go_ids
    t_df['go_type'] = go_types
    t_df['go_term'] = go_descriptions
    t_df['pathway_id'] = pathway_ids
    t_df['pathway_description'] = pathway_descriptions
    t_df['cog_id'] = cog_ids
    t_df['cog_functional_category'] = cog_functional_categories
    t_df['swissprot_hit_name'] = swissprot_hit_names
    t_df['swissprot_description'] = swissprot_descriptions
    t_df['pfam_id'] = pfam_ids
    t_df['pfam_domain'] = pfam_domains
    t_df['pfam_description'] = pfam_descriptions
    columns = ['transcript_id', 'gene_id', 'gene_name', 'description', 'nr_hit_name', 'nr_description', 'go_id',
               'go_type', 'go_term', 'ko_id', 'ko_name', 'pathway_id', 'pathway_description', 'cog_id',
               'cog_functional_category',
               'swissprot_hit_name', 'swissprot_description', 'pfam_id', 'pfam_domain', 'pfam_description', 'kind']
    t_df = t_df.reindex(columns, axis=1)
    ref_t_df = t_df[t_df['kind'] == 'ref'].drop('kind', axis=1)
    new_t_df = t_df[t_df['kind'] == 'new'].drop('kind', axis=1)
    all_t_df = t_df.drop('kind', axis=1)
    ref_t_df.to_csv(os.path.join(output_dir, 'ref_transcript_anno_detail.xls'), sep='\t', index=False)
    new_t_df.to_csv(os.path.join(output_dir, 'new_transcript_anno_detail.xls'), sep='\t', index=False)
    all_t_df.to_csv(os.path.join(output_dir, 'all_transcript_anno_detail.xls'), sep='\t', index=False)


def set_anno_statistics(task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    stat_id = database['annotation_stat'].find_one({'task_id': task_id})['main_id']
    for kind in ('ref', 'new', 'all'):
        data = list()
        cursor = database['annotation_stat_detail'].find({'stat_id': stat_id, 'kind': kind})
        for i, row in pd.DataFrame(list(cursor)).iterrows():
            data.append({
                'type': row['type'],
                'e_g': '{}({}%)'.format(row['exp_g_num'], round(row['exp_g_per'], 2)),
                'e_t': '{}({}%)'.format(row['exp_t_num'], round(row['exp_t_per'], 2)),
                'a_g': '{}({}%)'.format(row['g_num'], round(row['g_per'], 2)),
                'a_t': '{}({}%)'.format(row['t_num'], round(row['t_per'], 2))})
        df = pd.DataFrame(data).set_index('type')
        df = df.reindex(['e_g', 'e_t', 'a_g', 'a_t'], axis=1)
        df = df.rename({'e_g': 'Expre_Gene number (percent)', 'e_t': 'Expre_Transcript number (percent)',
                        'a_g': 'All_Gene number (percent)', 'a_t': 'All_Transcript number (percent)'}, axis=1)
        df.to_csv(os.path.join(output_dir, '{}_anno_stat.xls'.format(kind)), sep='\t')


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
    specimen_group_dict = database['specimen_group'].find_one({'task_id': task_id, 'library': 'long'})
    g_anno_df = pd.read_table(map_dict['g_anno'], index_col='gene_id')
    main_dict = database['exp'].find_one({'task_id': task_id, 'level': 'G'})
    exp_id, way = main_dict['main_id'], main_dict['way']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'mRNA'})
    g_detail_df = pd.DataFrame(list(cursor)).set_index('gene_id')
    g_detail_df = g_detail_df.reindex(samples, axis=1)
    for group, specimens in zip(specimen_group_dict['category_names'], specimen_group_dict['specimen_names']):
        g_detail_df[group] = g_detail_df.reindex(specimens, axis=1).mean(axis=1)
    g_exp_df = g_detail_df.join(g_anno_df)
    g_exp_df.to_csv(os.path.join(output_dir, 'gene_{}_anno.xls'.format(way)), sep='\t')
    row_index = g_exp_df.index
    g_count_df = pd.read_table(map_dict['g_count'], index_col=0)
    g_count_df.index.name = 'gene_id'
    g_count_df = g_count_df.reindex(row_index)
    g_count_df = g_count_df.join(g_anno_df)
    g_count_df.to_csv(os.path.join(output_dir, 'gene_count_anno.xls'), sep='\t')
    t_anno_df = pd.read_table(map_dict['t_anno'], index_col='transcript_id')
    main_dict = database['exp'].find_one({'task_id': task_id, 'level': 'T'})
    exp_id, way = main_dict['main_id'], main_dict['way']
    cursor = database['exp_detail'].find({'exp_id': exp_id, 'category': 'mRNA'})
    t_detail_df = pd.DataFrame(list(cursor)).set_index('transcript_id')
    t_detail_df = t_detail_df.reindex(samples, axis=1)
    for group, specimens in zip(specimen_group_dict['category_names'], specimen_group_dict['specimen_names']):
        t_detail_df[group] = t_detail_df.reindex(specimens, axis=1).mean(axis=1)
    t_exp_df = t_detail_df.join(t_anno_df)
    t_exp_df.to_csv(os.path.join(output_dir, 'transcript_{}_anno.xls'.format(way)), sep='\t')
    row_index = t_exp_df.index
    t_count_df = pd.read_table(map_dict['t_count'], index_col=0)
    t_count_df.index.name = 'transcript_id'
    t_count_df = t_count_df.reindex(row_index)
    t_count_df = t_count_df.join(t_anno_df)
    t_count_df.to_csv(os.path.join(output_dir, 'transcript_count_anno.xls'), sep='\t')


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
    if 't_count' in map_dict:
        way = database['exp'].find_one({'task_id': task_id, 'level': 'T'})['way']
        t_count_df = pd.read_table(map_dict['t_count'], index_col='transcript_id').reindex(samples, axis=1)
        t_count_df = t_count_df.rename({sample: '{}_count'.format(sample) for sample in samples}, axis=1)
        t_exp_df = pd.read_table(map_dict['t_exp'], index_col='transcript_id').reindex(samples, axis=1)
        t_exp_df = t_exp_df.rename({sample: '{}_{}'.format(sample, way) for sample in samples}, axis=1)
        t_anno_df = pd.read_table(map_dict['t_anno'], index_col='transcript_id')
        if main_id:
            main_dict = database['diff'].find_one({'main_id': ObjectId(main_id)})
        else:
            main_dict = database['diff'].find_one({'task_id': task_id, 'category': 'mRNA', 'level': 'T'})
        diff_id, params, cmp_list = main_dict['main_id'], main_dict['params'], main_dict['cmp_list']
        diff_method = json.loads(params)['diff_method'].lower()
        total_df = t_anno_df.copy()
        total_columns = list()
        df_dict = dict()
        for compare in cmp_list:
            ctrl, case = compare.split('|')
            vs_pair = '{}_vs_{}'.format(ctrl, case)
            div_pair = '{}/{}'.format(case, ctrl)
            logger.info('diff_id: {}, compare: {}'.format(diff_id, compare))
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
            diff_df.to_csv(os.path.join(output_dir, '{}_vs_{}_{}_T_anno.xls'.format(ctrl, case, diff_method)), sep='\t')
            columns = ['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
            diff_df = diff_df.reindex(columns, axis=1)
            mapper = {'fc': '{}_fc({})'.format(vs_pair, div_pair), 'log2fc': '{}_log2fc({})'.format(vs_pair, div_pair),
                      'pvalue': '{}_pvalue'.format(vs_pair), 'padjust': '{}_padjust'.format(vs_pair),
                      'significant': '{}_significant'.format(vs_pair), 'regulate': '{}_regulate'.format(vs_pair)}
            diff_df = diff_df.rename(mapper, axis=1)
            total_df = total_df.join(diff_df)
            total_columns.extend(diff_df.columns)
        anno_columns = list(t_anno_df.columns)
        anno_columns.remove('gene_id')
        anno_columns.remove('gene_name')
        anno_columns.remove('description')
        total_columns = ['gene_id', 'gene_name', 'description'] + total_columns + anno_columns
        total_df = total_df.reindex(total_columns, axis=1)
        total_df.to_csv(os.path.join(output_dir, 'total_diff_stat.T.{}_anno.xls'.format(diff_method)), sep='\t')
        cursor = database['diff_summary'].find({'diff_id': diff_id})
        summary = list(cursor)
        if len(summary) > 0:
            summary_df = pd.DataFrame(summary).set_index('seq_id')
            summary_df.index.name = 'transcript_id'
            vs_pairs = ['{}_vs_{}'.format(*compare.split('|')) for compare in cmp_list]
            columns = ['gene_id', 'gene_name', 'description'] + vs_pairs + ['sum'] + anno_columns
            summary_df = summary_df.join(t_anno_df)
            summary_df = summary_df.reindex(columns, axis=1)
            summary_df.to_csv(os.path.join(output_dir, 'diff_summary_T_{}_anno.xls'.format(diff_method)), sep='\t')
        else:
            output = os.path.join(output_dir, 'diff_summary_T_{}_anno.xls'.format(diff_method))
            cmd = 'touch {}'.format(output)
            os.system(cmd)
    if 'g_count' in map_dict:
        way = database['exp'].find_one({'task_id': task_id, 'level': 'G'})['way']
        g_count_df = pd.read_table(map_dict['g_count'], index_col='gene_id').reindex(samples, axis=1)
        g_count_df = g_count_df.rename({sample: '{}_count'.format(sample) for sample in samples}, axis=1)
        g_exp_df = pd.read_table(map_dict['g_exp'], index_col='gene_id').reindex(samples, axis=1)
        g_exp_df = g_exp_df.rename({sample: '{}_{}'.format(sample, way) for sample in samples}, axis=1)
        g_anno_df = pd.read_table(map_dict['g_anno'], index_col='gene_id')
        if main_id:
            main_dict = database['diff'].find_one({'main_id': ObjectId(main_id)})
        else:
            main_dict = database['diff'].find_one({'task_id': task_id, 'category': 'mRNA', 'level': 'G'})
        diff_id, params, cmp_list = main_dict['main_id'], main_dict['params'], main_dict['cmp_list']
        diff_method = json.loads(params)['diff_method'].lower()
        total_df = g_anno_df.copy()
        total_columns = list()
        df_dict = dict()
        for compare in cmp_list:
            ctrl, case = compare.split('|')
            vs_pair = '{}_vs_{}'.format(ctrl, case)
            div_pair = '{}/{}'.format(case, ctrl)
            cursor = database['diff_detail'].find({'diff_id': diff_id, 'compare': compare})
            diff_df = pd.DataFrame(list(cursor)).set_index('seq_id')
            diff_df.index.name = 'gene_id'
            diff_df = diff_df.join(g_count_df)
            diff_df = diff_df.join(g_exp_df)
            columns = list(g_count_df.columns) + list(g_exp_df.columns) + ['group1', 'group2', 'fc', 'log2fc', 'pvalue',
                                                                           'padjust', 'significant', 'regulate']
            diff_df = diff_df.reindex(columns, axis=1)
            diff_df.rename({'group1': '{}_mean_{}'.format(ctrl, way),
                            'group2': '{}_mean_{}'.format(case, way)}, axis=1, inplace=True)
            diff_df = diff_df.join(g_anno_df)
            diff_df.to_csv(os.path.join(output_dir, '{}_vs_{}_{}_G_anno.xls'.format(ctrl, case, diff_method)), sep='\t')
            columns = ['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
            diff_df = diff_df.reindex(columns, axis=1)
            mapper = {'fc': '{}_fc({})'.format(vs_pair, div_pair), 'log2fc': '{}_log2fc({})'.format(vs_pair, div_pair),
                      'pvalue': '{}_pvalue'.format(vs_pair), 'padjust': '{}_padjust'.format(vs_pair),
                      'significant': '{}_significant'.format(vs_pair), 'regulate': '{}_regulate'.format(vs_pair)}
            diff_df = diff_df.rename(mapper, axis=1)
            total_df = total_df.join(diff_df)
            total_columns.extend(diff_df.columns)
        anno_columns = list(g_anno_df.columns)
        anno_columns.remove('gene_name')
        anno_columns.remove('description')
        total_columns = ['gene_name', 'description'] + total_columns + anno_columns
        total_df = total_df.reindex(total_columns, axis=1)
        total_df.to_csv(os.path.join(output_dir, 'total_diff_stat.G.{}_anno.xls'.format(diff_method)), sep='\t')
        cursor = database['diff_summary'].find({'diff_id': diff_id})
        summary = list(cursor)
        if len(summary) > 0:
            summary_df = pd.DataFrame(summary).set_index('seq_id')
            summary_df.index.name = 'gene_id'
            vs_pairs = ['{}_vs_{}'.format(*compare.split('|')) for compare in cmp_list]
            columns = ['gene_name', 'description'] + vs_pairs + ['sum'] + anno_columns
            summary_df = summary_df.join(g_anno_df)
            summary_df = summary_df.reindex(columns, axis=1)
            summary_df.to_csv(os.path.join(output_dir, 'diff_summary_G_{}_anno.xls'.format(diff_method)), sep='\t')
        else:
            output = os.path.join(output_dir, 'diff_summary_G_{}_anno.xls'.format(diff_method))
            cmd = 'touch {}'.format(output)
            os.system(cmd)
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def set_as(map_dict, task_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    input_dir = map_dict['input_dir']
    # alternative splicing event statistics table
    os.link(os.path.join(input_dir, 'sample.event.count.JC.txt'), os.path.join(output_dir, 'event_type.JC.xls'))
    os.link(os.path.join(input_dir, 'sample.event.count.JCEC.txt'), os.path.join(output_dir, 'event_type.JCEC.xls'))
    # alternative splicing event detail table
    g_anno_df = pd.read_table(map_dict['g_anno'], usecols=['gene_id', 'description'])
    desc_dict = dict()
    for i, row in g_anno_df.iterrows():
        gene_id = row['gene_id']
        description = str() if pd.isna(row['description']) else row['description']
        desc_dict[gene_id] = description
    for document in database['splicing_rmats'].find({'task_id': task_id}):
        splicing_id = document['main_id']
        s1, s2 = document['compare_plan'].split('|')
        vs_pair_dir = os.path.join(output_dir, '{}_vs_{}'.format(s1, s2))
        export_rmats_detail(splicing_id, desc_dict, s2=s1, s1=s2, output_dir=vs_pair_dir)
    # alternative splicing intra-group event statistics table and pattern statistics table
    for document in database['splicing_rmats_stats'].find({'task_id': task_id}):
        stat_id = document['main_id']
        logger.debug('stat_id -> ({})'.format(stat_id))
        logger.debug('group -> ({})'.format(document['group']))
        for k, v in document['group'].items():
            if v == 's1':
                s1 = k
            if v == 's2':
                s2 = k
        else:
            logger.debug('s1 -> ({})'.format(s1))
            logger.debug('s2 -> ({})'.format(s2))
        vs_pair_dir = os.path.join(output_dir, '{}_vs_{}'.format(s1, s2))
        export_rmats_diff_stats(stat_id, vs_pair_dir)
        export_rmats_psi(stat_id, vs_pair_dir)
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))


def export_rmats_detail(splicing_id, gid2des_dct, s2, s1, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    dct = {k: {'JC': dict(), 'JCEC': dict()} for k in ('SE', 'MXE', 'A3SS', 'A5SS', 'RI')}
    for d in database['splicing_rmats_detail'].find({'splicing_id': ObjectId(splicing_id)}):
        eid = d['event_id']
        gname = d['gname'] if d['gname'] != 'nan' else '-'
        for jt in ('JC', 'JCEC'):
            if d['gid'] in gid2des_dct:
                desc = gid2des_dct[d['gid']]
            else:
                desc = ''
            dct[d['type']][jt][eid] = {'AS ID': d['event_id'], 'Gene ID': d['gid'], 'Gene name': gname,
                                       'Gene description': desc,
                                       'Novel AS': d['novel_as'], 'Chr': d['chr'], 'Strand': d['strand']}
            if d['type'] == 'SE':
                dct['SE'][jt][eid].update({
                    'InclusionTranscripts': d['InclusionTranscripts'],
                    'SkippingTranscripts': d['SkippingTranscripts'],
                    'ExonStart': d['es'],
                    'ExonEnd': d['ee'],
                    'UpstreamES': d['up_es'],
                    'UpstreamEE': d['up_ee'],
                    'DownstreamES': d['down_es'],
                    'DownstreamEE': d['down_ee']
                })
            elif d['type'] == 'MXE':
                dct['MXE'][jt][eid].update({
                    '1stExonTranscripts': d['1stExonTranscripts'],
                    '2ndExonTranscripts': d['2ndExonTranscripts'],
                    '1stExonStart': d['firstes'],
                    '1stExonEnd': d['firstee'],
                    '2ndExonStart': d['secondes'],
                    '2ndExonEnd': d['secondee']
                })
            elif d['type'] == 'A3SS':
                dct['A3SS'][jt][eid].update({
                    'LongExonTranscripts': d['LongExonTranscripts'],
                    'ShortExonTranscripts': d['ShortExonTranscripts'],
                    'LongExonStart': d['les'],
                    'LongExonEnd': d['lee'],
                    'ShortES': d['ses'],
                    'ShortEE': d['see'],
                    'FlankingES': d['fes'],
                    'FlankingEE': d['fee']
                })
            elif d['type'] == 'A5SS':
                dct['A5SS'][jt][eid].update({
                    'LongExonTranscripts': d['LongExonTranscripts'],
                    'ShortExonTranscripts': d['ShortExonTranscripts'],
                    'LongExonStart': d['les'],
                    'LongExonEnd': d['lee'],
                    'ShortES': d['ses'],
                    'ShortEE': d['see'],
                    'FlankingES': d['fes'],
                    'FlankingEE': d['fee']
                })
            elif d['type'] == 'RI':
                dct['RI'][jt][eid].update({
                    'RetainTranscripts': d['RetainTranscripts'],
                    'AbandonTranscripts': d['AbandonTranscripts'],
                    'RiExonStart': d['ries'],
                    'RiExonEnd': d['riee']
                })
        dct[d['type']]['JC'][eid].update({
            'Diff significant': d['diff_jc'],
            'IncLevelDiff ({}/{},ΔPSI)'.format(s2, s1): d['inc_diff_jc'],
            'P Value JunctionCountOnly': d['pvalue_jc'],
            'FDR JunctionCountOnly': d['fdr_jc'],
            'IJC {}'.format(s2): d['ijc_s1'],
            'SJC {}'.format(s2): d['sjc_s1'],
            'IJC {}'.format(s1): d['ijc_s2'],
            'SJC {}'.format(s1): d['sjc_s2'],
            'IncFormLen': d['inclen_jc'],
            'SkipFormLen': d['skiplen_jc'],
            'IncLevel1 (PSI {})'.format(s2): d['inc1_jc'],
            'IncLevel2 (PSI {})'.format(s1): d['inc2_jc'],
            'Average IncLevel1 (Average PSI {})'.format(s2): d['aver_inc1_jc'],
            'Average IncLevel2 (Average PSI {})'.format(s1): d['aver_inc2_jc'],
            'Increase Inclusion {}'.format(s2): d['upinc_s1_jc'],
            'Increase Exclusion {}'.format(s2): d['upexc_s1_jc'],
            'Increase Inclusion {}'.format(s1): d['upinc_s2_jc'],
            'Increase Exclusion {}'.format(s1): d['upexc_s2_jc']
        })
        dct[d['type']]['JCEC'][eid].update({
            'Diff significant': d['diff_all'],
            'IncLevelDiff ({}/{},ΔPSI)'.format(s2, s1): d['inc_diff_all'],
            'P Value ReadsOnTargetAndJunctionCounts': d['pvalue_all'],
            'FDR ReadsOnTargetAndJunctionCounts': d['fdr_all'],
            'IC {}'.format(s2): d['ic_s1'],
            'SC {}'.format(s2): d['sc_s1'],
            'IC {}'.format(s1): d['ic_s2'],
            'SC {}'.format(s1): d['sc_s2'],
            'IncFormLen': d['inclen_all'],
            'SkipFormLen': d['skiplen_all'],
            'IncLevel1 (PSI {})'.format(s2): d['inc1_all'],
            'IncLevel2 (PSI {})'.format(s1): d['inc2_all'],
            'Average IncLevel1 (Average PSI {})'.format(s2): d['aver_inc1_all'],
            'Average IncLevel2 (Average PSI {})'.format(s1): d['aver_inc2_all'],
            'Increase Inclusion {}'.format(s2): d['upinc_s1_all'],
            'Increase Exclusion {}'.format(s2): d['upexc_s1_all'],
            'Increase Inclusion {}'.format(s1): d['upinc_s2_all'],
            'Increase Exclusion {}'.format(s1): d['upexc_s2_all']
        })

    def get_columns(event_type, junction_type):
        lst = ['AS ID', 'Gene ID', 'Gene name', 'Gene description', 'Novel AS', 'Chr', 'Strand',
               'Diff significant', 'IncLevelDiff ({}/{},ΔPSI)'.format(s1, s2)]
        lst.extend({
                       'JC': ['P Value JunctionCountOnly', 'FDR JunctionCountOnly'],
                       'JCEC': ['P Value ReadsOnTargetAndJunctionCounts', 'FDR ReadsOnTargetAndJunctionCounts']
                   }[junction_type])
        lst.extend({
                       'SE': ['InclusionTranscripts', 'SkippingTranscripts', 'ExonStart', 'ExonEnd', 'UpstreamES',
                              'UpstreamEE', 'DownstreamES', 'DownstreamEE'],
                       'MXE': ['1stExonTranscripts', '2ndExonTranscripts', '1stExonStart', '1stExonEnd', '2ndExonStart',
                               '2ndExonEnd'],
                       'A3SS': ['LongExonTranscripts', 'ShortExonTranscripts', 'LongExonStart', 'LongExonEnd',
                                'ShortES', 'ShortEE', 'FlankingES', 'FlankingEE'],
                       'A5SS': ['LongExonTranscripts', 'ShortExonTranscripts', 'LongExonStart', 'LongExonEnd',
                                'ShortES', 'ShortEE', 'FlankingES', 'FlankingEE'],
                       'RI': ['RetainTranscripts', 'AbandonTranscripts', 'RiExonStart', 'RiExonEnd']
                   }[event_type])
        lst.extend({
                       'JC': ['IJC {}'.format(s2), 'SJC {}'.format(s2), 'IJC {}'.format(s1), 'SJC {}'.format(s1)],
                       'JCEC': ['IC {}'.format(s2), 'SC {}'.format(s2), 'IC {}'.format(s1), 'SC {}'.format(s1)]
                   }[junction_type])
        lst.extend([
            'IncFormLen', 'SkipFormLen', 'IncLevel1 (PSI {})'.format(s2), 'IncLevel2 (PSI {})'.format(s1),
            'Average IncLevel1 (Average PSI {})'.format(s2), 'Average IncLevel2 (Average PSI {})'.format(s1),
            'Increase Inclusion {}'.format(s2), 'Increase Exclusion {}'.format(s2),
            'Increase Inclusion {}'.format(s1), 'Increase Exclusion {}'.format(s1)
        ])
        return lst

    if not os.path.isdir(os.path.join(output_dir, 'JC')):
        os.mkdir(os.path.join(output_dir, 'JC'))
    if not os.path.isdir(os.path.join(output_dir, 'JCEC')):
        os.mkdir(os.path.join(output_dir, 'JCEC'))
    for et, jt2dcts in dct.items():
        logger.debug('event_type -> ({})'.format(et))
        for jt, dcts in jt2dcts.items():
            logger.debug('junction_type -> ({})'.format(jt))
            df = pd.DataFrame(dcts.values())
            df = df.reindex(get_columns(et, jt), axis=1)
            df.to_csv(os.path.join(output_dir, '{}/{}_detail.xls'.format(jt, et)), sep='\t', index=False)
    else:
        logger.debug('succeed in exporting files to {}'.format(output_dir))


def export_rmats_diff_stats(stat_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    document = database['splicing_rmats_diff_stats'].find_one({'stat_id': ObjectId(stat_id)})
    index = ['JunctionCountOnly(JC)', 'ReadsOnTargetAndJunctionCounts(JCEC)', 'JC&JCEC', 'JC|JCEC']
    df = pd.DataFrame(document['diff_stats'], index=index).T
    df = df.rename({i: i.upper() for i in df.index})
    df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
    df.index.name = 'AS type'
    df.to_csv(os.path.join(output_dir, 'diff_event_stat.xls'), sep='\t')
    logger.debug('succeed in exporting {}'.format(os.path.join(output_dir, 'diff_event_stat.txt')))


def export_rmats_psi(stat_id, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    document = database['splicing_rmats_psi'].find_one({'stat_id': ObjectId(stat_id)})
    index = ['Exclusion (Increase exclusion in case，ΔPSI<0)',
             'Inclusion (Increase inclusion in case, ΔPSI>0)',
             'Total events']
    df = pd.DataFrame(document['s1_jc'], index=index).T
    df = df.rename({i: i.upper() for i in df.index})
    df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
    df.index.name = 'AS type'
    jc_ofile = os.path.join(output_dir, 'diff_pattern_stat.JC.xls')
    df.to_csv(jc_ofile, sep='\t')
    logger.debug('succeed in exporting {}'.format(jc_ofile))
    df = pd.DataFrame(document['s1_all'], index=index).T
    df = df.rename({i: i.upper() for i in df.index})
    df = df.reindex(['SE', 'MXE', 'A3SS', 'A5SS', 'RI', 'TOTAL'])
    df.index.name = 'AS type'
    jcec_ofile = os.path.join(output_dir, 'diff_pattern_stat.JCEC.xls')
    df.to_csv(jcec_ofile, sep='\t')
    logger.debug('succeed in exporting {}'.format(jcec_ofile))


def set_snp_indel(map_dict, output_dir):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    upload_dir = map_dict['upload']
    if os.path.exists(os.path.join(upload_dir, 'snp_anno.xls')):
        snp_detail = pd.read_table(os.path.join(upload_dir, 'snp_anno.xls'), low_memory=False)
        snp_detail.drop(['Total depth', 'type'], axis=1, inplace=True)
        snp_detail.to_csv(os.path.join(output_dir, 'snp_anno.xls'), index=False, sep='\t')
    if os.path.exists(os.path.join(upload_dir, 'indel_anno.xls')):
        indel_detail = pd.read_table(os.path.join(upload_dir, 'indel_anno.xls'), low_memory=False)
        indel_detail.drop(['Total depth', 'type'], axis=1, inplace=True)
        indel_detail.to_csv(os.path.join(output_dir, 'indel_anno.xls'), index=False, sep='\t')
    os.link(os.path.join(upload_dir, 'snp_position_distribution.xls'),
            os.path.join(output_dir, 'snp_position_distribution.xls'))
    os.link(os.path.join(upload_dir, 'snp_transition_tranversion_statistics.xls'),
            os.path.join(output_dir, 'snp_transition_tranversion_statistics.xls'))
    os.link(os.path.join(upload_dir, 'snp_depth_statistics.xls'), os.path.join(output_dir, 'snp_depth_statistics.xls'))
    os.link(os.path.join(upload_dir, 'indel_position_distribution.xls'),
            os.path.join(output_dir, 'indel_position_distribution.xls'))
    logger.info('succeed in calling {}'.format(sys._getframe().f_code.co_name))

def set_genesets_analysis(idir,odir):
    if os.path.isdir(odir):
        shutil.rmtree(odir)
    else:
        os.mkdir(odir)
    diff_genesets = []
    for gt in os.listdir(os.path.join(idir)):
        if gt != "cluster":
            diff_genesets.append(gt)
    target_geneset = diff_genesets[0]
    #聚类分析
    os.makedirs(os.path.join(odir, '01_{}_Cluster_Analysis'.format(target_geneset)))
    cluster_dir = os.path.join(odir, '01_{}_Cluster_Analysis'.format(target_geneset))
    raw_exp_file = os.path.join(idir,"cluster",target_geneset,"expression_matrix.xls")
    linkfile(raw_exp_file,os.path.join(cluster_dir,"expression_matrix.xls"))
    subcluster = glob.glob(os.path.join(idir,"cluster",target_geneset,'seq.subcluster_*.xls'))
    for file in subcluster:
        sub = os.path.basename(file).strip().split("_")[1]
        linkfile(file, os.path.join(cluster_dir, "subcluster_" + sub + ".xls"))

    #COG注释
    os.makedirs(os.path.join(odir, '02_{}_COG_Annotation'.format(target_geneset)))
    cog_class_dir = os.path.join(odir, '02_{}_COG_Annotation'.format(target_geneset))
    raw_cog_class_file = os.path.join(idir,target_geneset,"diff_cog_class","cog_class_table.xls")
    linkfile(raw_cog_class_file,os.path.join(cog_class_dir,"cog_class_table.xls"))

    #GO注释
    os.makedirs(os.path.join(odir, '03_{}_GO_Annotation'.format(target_geneset)))
    go_class_dir = os.path.join(odir, '03_{}_GO_Annotation'.format(target_geneset))
    raw_go_class_file = os.path.join(idir, target_geneset, "diff_go_class", "go_class_table.xls")
    linkfile(raw_go_class_file, os.path.join(go_class_dir, "go_class_table.xls"))

    #KEGG注释
    os.makedirs(os.path.join(odir, '04_{}_KEGG_Annotation'.format(target_geneset)))
    kegg_class_dir = os.path.join(odir, '04_{}_KEGG_Annotation'.format(target_geneset))
    raw_kegg_class_file = os.path.join(idir, target_geneset, "diff_kegg_class", "kegg_stat.xls")
    linkfile(raw_kegg_class_file, os.path.join(kegg_class_dir, "kegg_stat.xls"))
    raw_kegg_pathways_file = os.path.join(idir, target_geneset, "diff_kegg_class", "pathways.tar.gz")
    linkfile(raw_kegg_pathways_file, os.path.join(kegg_class_dir, "pathways.tar.gz"))

    #GO富集
    os.makedirs(os.path.join(odir, '05_{}_GO_Enrich'.format(target_geneset)))
    go_enrich_dir = os.path.join(odir, '05_{}_GO_Enrich'.format(target_geneset))
    raw_go_enrich_file = os.path.join(idir, target_geneset, "diff_go_enrich", "go_enrich_geneset_list_gene.xls")
    linkfile(raw_go_enrich_file, os.path.join(go_enrich_dir, "go_enrich_stat.xls"))

    #KEGG富集
    os.makedirs(os.path.join(odir, '06_{}_KEGG_Enrich'.format(target_geneset)))
    kegg_enrich_dir = os.path.join(odir, '06_{}_KEGG_Enrich'.format(target_geneset))
    raw_kegg_enrich_file = os.path.join(idir, target_geneset, "diff_kegg_enrich","enrich", "{}_gene.list.DE.list.check.kegg_enrichment.xls".format(target_geneset))
    linkfile(raw_kegg_enrich_file, os.path.join(kegg_enrich_dir, "kegg_erich_stat.xls"))
    raw_kegg_enrich_pathways_file = os.path.join(idir, target_geneset, "diff_kegg_enrich","class", "pathways.tar.gz")
    linkfile(raw_kegg_enrich_pathways_file, os.path.join(kegg_enrich_dir, "pathways.tar.gz"))

    #基因集富集
    os.makedirs(os.path.join(odir,"07_GenesetVenn"))





if __name__ == '__main__':
    import unittest


    class Empty(object):
        def __init__(self):
            self.work_dir = None
            self.output_dir = None

        pass


    class TestFunction(unittest.TestCase):
        """
        This is test for the package. Just run this script to do test.
        """

        def test(self):
            from mbio.packages.whole_transcriptome.catalogue import mrna
            self.task_id = 'tsg_36088'
            self.work_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191106/WholeTranscriptome_tsg_36088'
            self.output_dir = os.path.join(self.work_dir, 'output')
            self.tools = {'transfer_l': Empty()}
            self.tools['transfer_l'].output_dir = os.path.join(self.work_dir, 'Transfer/output')
            self.modules = {'rmats': Empty(), 'whole_snp': Empty()}
            self.modules['rmats'].output_dir = os.path.join(self.work_dir, 'Rmats/output')
            self.modules['whole_snp'].work_dir = os.path.join(self.work_dir, 'WholeSnp')
            self.long_task_info = database['task'].find_one({'task_id': 'tsg_36023'})
            mrna_dir = os.path.join(self.output_dir, 'mrna')
            if os.path.isdir(mrna_dir):
                shutil.rmtree(mrna_dir)
            os.mkdir(mrna_dir)
            mrna.set_background(self.task_id, os.path.join(mrna_dir, '01_Background'))
            map_dict = {
                'new_gtf': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/new.gtf'),
                'all_gtf': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/all.gtf'),
                'new_fasta': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/new.fasta'),
                'all_fasta': os.path.join(self.tools['transfer_l'].output_dir, 'assembly/all.fasta')}
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
                't_anno': os.path.join(mrna_dir, '03_Annotation/01_Anno_Detail/all_transcript_anno_detail.xls')}
            mrna.set_diff_express(map_dict, self.task_id, os.path.join(mrna_dir, '05_Diff_Express'))
            map_dict = {'input_dir': self.modules['rmats'].output_dir,
                        'g_anno': os.path.join(mrna_dir, '03_Annotation/01_Anno_Detail/all_gene_anno_detail.xls')}
            mrna.set_as(map_dict, self.task_id, os.path.join(mrna_dir, '06_AS'))
            map_dict = {'upload': os.path.join(self.modules['whole_snp'].work_dir, 'upload')}
            mrna.set_snp_indel(map_dict, os.path.join(mrna_dir, '07_SNP_InDel'))


    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
