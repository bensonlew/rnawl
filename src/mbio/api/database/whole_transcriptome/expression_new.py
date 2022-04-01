# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import pickle
import re
import unittest
import os
import pandas as pd
from bson.objectid import ObjectId

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class ExpressionNew(ApiBase):
    def __init__(self, bind_object):
        super(ExpressionNew, self).__init__(bind_object)

    def add_exp_g(self, map_dict, categories, way, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Expression_G_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'expdetail', 'task_type': 2}, sort_keys=True,
                            separators=(',', ':'))
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'level': 'G',
            'way': way,
            'status': 'start',
            'version': 'v1',
            'batch_version': 0,
            'is_rmbe': 'false'
        }
        main_id = self.create_db_table('exp', [main_dict])
        self.add_exp_detail(map_dict['mRNA_{}'.format(way)], main_id, 'G', 'mRNA', 'all', way)
        if 'lncRNA' in categories:
            self.add_exp_detail(map_dict['lncRNA_{}'.format(way)], main_id, 'G', 'lncRNA', 'all', way)
        self.update_db_record('exp', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_exp_t(self, map_dict, categories, way, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Expression_T_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'expdetail', 'task_type': 2}, sort_keys=True,
                            separators=(',', ':'))
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'level': 'T',
            'way': way,
            'status': 'start',
            'version': 'v1',
            'batch_version': 0,
            'is_rmbe': 'false'
        }
        main_id = self.create_db_table('exp', [main_dict])
        self.add_exp_detail(map_dict['mRNA_{}'.format(way)], main_id, 'T', 'mRNA', 'all', way)
        if 'lncRNA' in categories:
            self.add_exp_detail(map_dict['lncRNA_{}'.format(way)], main_id, 'T', 'lncRNA', 'all', way)
        if 'circRNA' in categories:
            self.add_exp_detail(map_dict['circRNA_rpm'], main_id, 'T', 'circRNA', 'all', 'rpm')
        if 'miRNA' in categories:
            self.add_exp_detail(map_dict['miRNA_tpm'], main_id, 'T', 'miRNA', 'all', 'tpm')
        self.update_db_record('exp', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_exp_mrna(self, map_dict, exp_method, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Expression_mRNA_{}_{}'.format(exp_method, time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'expdetail', 'task_type': 2}, sort_keys=True,
                            separators=(',', ':'))
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'category': 'mRNA',
            'status': 'start'
        }
        main_id = self.create_db_table('sg_exp', [main_dict])
        self.add_exp_detail(map_dict['G_fpkm'], main_id, 'G', 'mRNA', 'all', 'fpkm')
        self.add_exp_detail(map_dict['G_tpm'], main_id, 'G', 'mRNA', 'all', 'tpm')
        self.add_exp_detail(map_dict['T_fpkm'], main_id, 'T', 'mRNA', 'all', 'fpkm')
        self.add_exp_detail(map_dict['T_tpm'], main_id, 'T', 'mRNA', 'all', 'tpm')
        self.update_db_record('sg_exp', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_exp_lncrna(self, map_dict, exp_method, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Expression_lncRNA_{}_{}'.format(exp_method, time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'expdetail', 'task_type': 2}, sort_keys=True,
                            separators=(',', ':'))
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'category': 'lncRNA',
            'status': 'start'
        }
        main_id = self.create_db_table('sg_exp', [main_dict])
        self.add_exp_detail(map_dict['G_fpkm'], main_id, 'G', 'lncRNA', 'all', 'fpkm')
        self.add_exp_detail(map_dict['G_tpm'], main_id, 'G', 'lncRNA', 'all', 'tpm')
        self.add_exp_detail(map_dict['T_fpkm'], main_id, 'T', 'lncRNA', 'all', 'fpkm')
        self.add_exp_detail(map_dict['T_tpm'], main_id, 'T', 'lncRNA', 'all', 'tpm')
        self.update_db_record('sg_exp', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_exp_detail(self, exp_matrix, exp_id, level, category, kind, way):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind, 'way': way}
        ))
        df = pd.read_table(exp_matrix)
        df = df.fillna(str())
        df['batch_version'] = 0
        docs = df.to_dict('r')
        self.create_db_table('exp_detail', docs, {'exp_id': exp_id})

    def add_exp_graph(self, map_dict, category, exp_id, level, kind, group_id, group_dict, task_id, project_sn,
                      graph_id=None):
        if graph_id:
            graph_id = ObjectId(graph_id)
        else:
            time_now = datetime.datetime.now()
            name = 'Graph_{}_{}_{}'.format(level, category, time_now.strftime('%Y%m%d_%H%M%S'))
            params = json.dumps({
                'task_id': task_id,
                'submit_location': 'expgraph_{}'.format(category.lower()),
                'task_type': 2,
                'category': category,
                'level': level,
                'kind': kind,
                'group_id': str(group_id),
                'group_dict': group_dict,
                'exp_id': str(exp_id)
            }, sort_keys=True, separators=(',', ':'))
            main_dict = {
                'task_id': task_id,
                'project_sn': project_sn,
                'name': name,
                'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                'desc': 'Workflow result',
                'params': params,
                'category': category,
                'status': 'start',
                'version': 'v1'
            }
            graph_id = self.create_db_table('exp_graph', [main_dict])
        self.add_exp_graph_box(map_dict['sample_box'], graph_id, level, category, kind, 'sample')
        self.add_exp_graph_box(map_dict['group_box'], graph_id, level, category, kind, 'group')
        self.add_exp_graph_density(map_dict['sample_density'], graph_id, level, category, kind, 'sample')
        self.add_exp_graph_density(map_dict['group_density'], graph_id, level, category, kind, 'group')
        self.add_exp_graph_volin(map_dict['sample_volin'], graph_id, level, category, kind, 'sample')
        self.add_exp_graph_volin(map_dict['group_volin'], graph_id, level, category, kind, 'group')
        self.update_db_record('exp_graph', graph_id, insert_dict={'main_id': graph_id, 'status': 'end'})
        self.use_group_compare(group_id)

    def add_exp_graph_box(self, data_pk, graph_id, level, category, kind, field):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind, 'field': field}
        ))
        docs = pickle.load(open(data_pk))
        if docs:
            self.create_db_table('exp_graph_box', docs, {'graph_id': graph_id, 'type': field})

    def add_exp_graph_density(self, data_pk, graph_id, level, category, kind, field):
        if os.path.isfile(data_pk):
            self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                {'level': level, 'category': category, 'kind': kind, 'field': field}
            ))
            docs = pickle.load(open(data_pk))
            if docs:
                self.create_db_table('exp_graph_density', docs, {'graph_id': graph_id, 'type': field})

    def add_exp_graph_volin(self, data_pk, graph_id, level, category, kind, field):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind, 'field': field}
        ))
        docs = pickle.load(open(data_pk))
        if docs:
            self.create_db_table('exp_graph_volin', docs, {'graph_id': graph_id, 'type': field})

    def add_exp_venn(self, map_dict, category, exp_id, level, kind, group_id, group_dict, threshold, task_id,
                     project_sn, venn_id=None):
        if venn_id:
            venn_id = ObjectId(venn_id)
        else:
            time_now = datetime.datetime.now()
            name = 'Venn_{}_{}_{}'.format(level, category, time_now.strftime('%Y%m%d_%H%M%S'))
            params = json.dumps({
                'task_id': task_id,
                'submit_location': 'expvenn_{}'.format(category.lower()),
                'task_type': 2,
                'category': category,
                'level': level,
                'kind': kind,
                'group_id': str(group_id),
                'group_dict': group_dict,
                'threshold': str(threshold),
                'exp_id': str(exp_id)
            }, sort_keys=True, separators=(',', ':'))
            main_dict = {
                'task_id': task_id,
                'project_sn': project_sn,
                'name': name,
                'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                'desc': 'Workflow result',
                'params': params,
                'category': category,
                'status': 'start',
                'version': 'v1'
            }
            venn_id = self.create_db_table('exp_venn', [main_dict])
        self.add_exp_venn_detail(map_dict['venn'], venn_id, level, category, kind)
        self.update_db_record('exp_venn', venn_id, insert_dict={'main_id': venn_id, 'status': 'end'})
        self.use_group_compare(group_id)

    def add_exp_venn_detail(self, venn_table, venn_id, level, category, kind):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind}
        ))
        df = pd.read_table(venn_table)
        docs = df.to_dict('r')
        self.create_db_table('exp_venn_detail', docs, {'venn_id': venn_id})

    def add_exp_corr(self, map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, corr_method, take_log,
                     dist_method, clus_method, task_id, project_sn, corr_id=None):
        if corr_id:
            corr_id = ObjectId(corr_id)
        else:
            time_now = datetime.datetime.now()
            name = 'Corr_{}_{}'.format(level, time_now.strftime('%Y%m%d_%H%M%S'))
            params = json.dumps({
                'task_id': task_id,
                'submit_location': 'expcorr',
                'task_type': 2,
                'library': library,
                'level': level,
                'kind': kind,
                'group_id': str(group_id),
                'group_dict': group_dict,
                'take_mean': str(take_mean),
                'corr_method': corr_method,
                'take_log': str(take_log),
                'dist_method': dist_method,
                'clus_method': clus_method,
                'exp_id': str(exp_id)
            }, sort_keys=True, separators=(',', ':'))
            main_dict = {
                'task_id': task_id,
                'project_sn': project_sn,
                'name': name,
                'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                'desc': 'Workflow result',
                'params': params,
                'library': library,
                'status': 'start',
                'version': 'v1'
            }
            corr_id = self.create_db_table('exp_corr', [main_dict])
        clust_tree = open(map_dict['tree']).read().strip()
        samples = re.findall('[(,]([^(]*?):', clust_tree)
        insert_dict = {'main_id': corr_id, 'status': 'end', 'clust_tree': clust_tree, 'samples': samples}
        self.add_exp_corr_detail(map_dict['corr'], corr_id, level, library, kind, way)
        self.update_db_record('exp_corr', corr_id, insert_dict=insert_dict)
        self.use_group_compare(group_id)

    def add_exp_corr_detail(self, corr_table, corr_id, level, library, kind, way):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'library': library, 'kind': kind, 'way': way}
        ))
        df = pd.read_table(corr_table)
        docs = df.to_dict('r')
        self.create_db_table('exp_corr_detail', docs, {'corr_id': corr_id})

    def add_exp_pca(self, map_dict, library, exp_id, level, kind, group_id, group_dict, way, take_mean, task_id, project_sn,
                    pca_id=None):
        if pca_id:
            pca_id = ObjectId(pca_id)
        else:
            time_now = datetime.datetime.now()
            name = 'PCA_{}_{}'.format(level, time_now.strftime('%Y%m%d_%H%M%S'))
            params = json.dumps({
                'task_id': task_id,
                'submit_location': 'exppca',
                'task_type': 2,
                'library': library,
                'level': level,
                'kind': kind,
                'group_id': str(group_id),
                'group_dict': group_dict,
                'take_mean': str(take_mean),
                'exp_id': str(exp_id)
            }, sort_keys=True, separators=(',', ':'))
            main_dict = {
                'task_id': task_id,
                'project_sn': project_sn,
                'name': name,
                'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                'desc': 'Workflow result',
                'params': params,
                'library': library,
                'status': 'start',
                'version': 'v1'
            }
            pca_id = self.create_db_table('exp_pca', [main_dict])
        ratio_df = pd.read_table(map_dict['evr'], header=None, index_col=0)
        ratio_dict = dict(ratio_df.to_records())
        insert_dict = {'main_id': pca_id, 'ratio_dict': ratio_dict, 'status': 'end'}
        self.add_exp_pca_detail(map_dict['pca'], pca_id, level, library, kind, way)
        if 'ellipse' in map_dict:
            self.add_exp_pca_circ_detail(map_dict['ellipse'], pca_id, level, library, kind, way)
        else:
            self.update_db_record('exp_pca', pca_id, insert_dict={'ellipse': False})
        self.update_db_record('exp_pca', pca_id, insert_dict=insert_dict)
        self.use_group_compare(group_id)

    def add_exp_pca_detail(self, pca_table, pca_id, level, library, kind, way):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'library': library, 'kind': kind, 'way': way}
        ))
        df = pd.read_table(pca_table)
        docs = df.to_dict('r')
        self.create_db_table('exp_pca_detail', docs, {'pca_id': pca_id})

    def add_exp_pca_circ_detail(self, ellipse_table, pca_id, level, library, kind, way):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'library': library, 'kind': kind, 'way': way}
        ))
        df = pd.read_table(ellipse_table)
        doc_dict = dict()
        for idx, row in df.iterrows():
            if row['PC'] in doc_dict:
                doc_dict[row['PC']].update({row['group']: json.dumps(row[2:].to_dict())})
            else:
                doc_dict[row['PC']] = {'name': row['PC'], row['group']: json.dumps(row[2:].to_dict())}
        else:
            docs = doc_dict.values()
            self.create_db_table('exp_pca_circ_detail', docs, {'type': 'circle', 'exp_pca_ellipse_id': pca_id})
            self.update_db_record('exp_pca', pca_id, insert_dict={'ellipse': True})

    def add_diff(self, map_dict, category, level, kind, background, filter, threshold, group_id, group_dict, control_id,
                 stat_type, stat_cutoff, fc, diff_method, correct_method, way, task_id, project_sn, exp_id, diff_id=None):
        if diff_id:
            diff_id = ObjectId(diff_id)
        else:
            time_now = datetime.datetime.now()
            name = 'DE_{}_{}_{}'.format(level, category, time_now.strftime('%Y%m%d_%H%M%S'))
            if diff_method.lower() in ['degseq', 'deseq2', 'edger', 'limma']:
                param_dict = {
                    'task_id': task_id,
                    'submit_location': 'diff_{}'.format(category).lower(),
                    'task_type': 2,
                    'category': category,
                    'level': level,
                    'kind': kind,
                    'filter': filter,
                    'threshold': str(threshold),
                    'group_id': str(group_id),
                    'group_dict': group_dict,
                    'control_id': str(control_id),
                    'stat_type': stat_type,
                    'stat_cutoff': str(stat_cutoff),
                    'fc': str(fc),
                    'diff_method': diff_method,
                    'correct_method': correct_method,
                    'exp_id': exp_id,
                    'is_batch': 'False'
                }
            else:
                param_dict = {
                    'task_id': task_id,
                    'submit_location': 'diff_{}'.format(category).lower(),
                    'task_type': 2,
                    'category': category,
                    'level': level,
                    'kind': kind,
                    'filter': filter,
                    'threshold': str(threshold),
                    'group_id': str(group_id),
                    'group_dict': group_dict,
                    'control_id': str(control_id),
                    'prob': str(stat_cutoff),
                    'fc': str(fc),
                    'diff_method': diff_method,
                    'correct_method': correct_method,
                    'exp_id': exp_id,
                    'is_batch': 'False'
                }
            if background:
                param_dict['background'] = background
            params = json.dumps(param_dict, sort_keys=True, separators=(',', ':'))
            main_dict = {
                'task_id': task_id,
                'project_sn': project_sn,
                'name': name,
                'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
                'desc': 'Workflow result',
                'params': params,
                'level': level,
                'category': category,
                'status': 'start',
                'version': 'v1.1'
            }
            diff_id = self.create_db_table('diff', [main_dict])
        # exp_id = self.db['exp'].find_one({'task_id': task_id, 'level': level})['main_id']
        match_dict = {'exp_id': ObjectId(exp_id), 'category': category}
        if kind != 'all':
            match_dict['kind'] = kind
        # exp_count = len(list(self.db['exp_detail'].find(match_dict)))
        exp_count = self.db['exp_detail'].find(match_dict).count()
        compare_names = self.db['specimen_group_compare'].find_one(
            {'task_id': task_id, 'main_id': ObjectId(control_id)})['compare_names']
        cmp_list = json.loads(compare_names)
        sig_status = dict()
        cmp_detail = dict()
        volcano_df = pd.read_table(map_dict['volcano'])
        relation_df = self.get_relation_df(task_id)
        self.bind_object.logger.info('succeed in preparing data for dumping into mongo')
        for compare in cmp_list:
            issig_volcanno_df = volcano_df[(volcano_df['significant'] == 'yes') & (volcano_df['compare'] == compare)]
            up_count = issig_volcanno_df[issig_volcanno_df['regulate'] == 'up'].shape[0]
            down_count = issig_volcanno_df[issig_volcanno_df['regulate'] == 'down'].shape[0]
            nosig_count = exp_count - up_count - down_count
            sig_status[compare] = ['nosig_{}'.format(nosig_count), 'down_{}'.format(down_count),
                                   'up_{}'.format(up_count)]
            ctrl, case = compare.split('|')
            cmp_detail[compare] = group_dict[ctrl] + group_dict[case]
            self.add_diff_detail(map_dict['{}_vs_{}'.format(ctrl, case)], diff_id, level, category, kind)
        self.add_diff_summary(map_dict['summary'], diff_id, level, category, kind)
        self.add_diff_volcano(map_dict['volcano'], diff_id, level, category, kind, relation_df)
        self.add_diff_scatter(map_dict['scatter'], diff_id, level, category, kind, relation_df)
        insert_dict = {'main_id': diff_id, 'cmp_list': cmp_list, 'sig_status': sig_status, 'cmp_detail': cmp_detail,
                       'status': 'end'}
        self.update_db_record('diff', diff_id, insert_dict=insert_dict)
        self.use_group_compare(group_id, control_id)

    def add_diff_detail(self, detail_table, diff_id, level, category, kind):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind}
        ))
        df = pd.read_table(detail_table)
        docs = df.to_dict('r')
        self.create_db_table('diff_detail', docs, {'diff_id': diff_id})

    def add_diff_summary(self, summary_table, diff_id, level, category, kind):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind}
        ))
        df = pd.read_table(summary_table)
        docs = df.to_dict('r')
        self.create_db_table('diff_summary', docs, {'diff_id': diff_id})

    def add_diff_volcano(self, volcano_table, diff_id, level, category, kind, relation_df):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind}
        ))
        seq_id = 'transcript_id' if level == 'T' else 'gene_id'
        relation_df = relation_df.rename({seq_id: 'seq_id'}, axis=1).reindex(['seq_id', 'gene_name'], axis=1)
        volcano_df = pd.read_table(volcano_table)
        df = pd.merge(volcano_df, relation_df, how='left', on='seq_id')
        sig_df = df[df['significant'] == 'yes']
        nosig_df = df[df['significant'] != 'yes']
        flag = nosig_df.shape[0] > 50000
        while flag:
            nosig_df = nosig_df.sample(frac=0.9)
            flag = nosig_df.shape[0] > 50000
        df = pd.concat([sig_df, nosig_df])
        docs = df.to_dict('r')
        self.create_db_table('diff_volcano', docs, {'diff_id': diff_id})

    def add_diff_scatter(self, scatter_table, diff_id, level, category, kind, relation_df):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'level': level, 'category': category, 'kind': kind}
        ))
        seq_id = 'transcript_id' if level == 'T' else 'gene_id'
        relation_df = relation_df.rename({seq_id: 'seq_id'}, axis=1).reindex(['seq_id', 'gene_name'], axis=1)
        scatter_df = pd.read_table(scatter_table)
        df = pd.merge(scatter_df, relation_df, how='left', on='seq_id')
        sig_df = df[df['significant'] == 'yes']
        nosig_df = df[df['significant'] != 'yes']
        flag = nosig_df.shape[0] > 50000
        while flag:
            nosig_df = nosig_df.sample(frac=0.9)
            flag = nosig_df.shape[0] > 50000
        df = pd.concat([sig_df, nosig_df])
        docs = df.to_dict('r')
        self.create_db_table('diff_scatter', docs, {'diff_id': diff_id})

    def get_relation_df(self, task_id):
        relation_data = list()
        exp_id = self.db['exp'].find_one({'task_id': task_id, 'level': 'T'})['main_id']
        for document in self.db['exp_detail'].find({'exp_id': exp_id, 'level': 'T'}):
            relation_dict = {'transcript_id': document['transcript_id']}
            relation_dict['gene_id'] = document['gene_id'] if 'gene_id' in document else str()
            relation_dict['gene_name'] = document['gene_name'] if 'gene_name' in document else str()
            relation_data.append(relation_dict)
        else:
            relation_df = pd.DataFrame(relation_data)
            return relation_df

    def get_group_dict(self, group_id):
        group_dict = dict()
        document = self.db['specimen_group'].find_one({'main_id': ObjectId(group_id)})
        for group, samples in zip(document['category_names'], document['specimen_names']):
            group_dict[group] = sorted(list(set(samples)))
        return group_dict

    def use_group_compare(self, group_id=None, control_id=None):
        if group_id and group_id != 'all':
            self.update_db_record('specimen_group', ObjectId(group_id), insert_dict={'is_use': 1})
        if control_id:
            self.update_db_record('specimen_group_compare', ObjectId(control_id), insert_dict={'is_use': 1})


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test_exp_by_level(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression')
        map_dict = {
            'mRNA_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/G.fpkm.txt',
            'mRNA_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/G.tpm.txt',
            'lncRNA_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/G.fpkm.txt',
            'lncRNA_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/G.tpm.txt',
        }
        wf.test_api.add_exp_g(
            map_dict=map_dict,
            categories='mRNA,lncRNA',
            way='tpm',
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )
        map_dict = {
            'mRNA_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/T.fpkm.txt',
            'mRNA_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/T.tpm.txt',
            'lncRNA_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/T.fpkm.txt',
            'lncRNA_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/T.tpm.txt',
            'circRNA_rpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/circrna/T.rpm.txt',
        }
        wf.test_api.add_exp_t(
            map_dict=map_dict,
            categories='mRNA,lncRNA,circRNA',
            way='tpm',
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_exp_by_category(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression')
        map_dict = {
            'G_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/G.fpkm.txt',
            'G_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/G.tpm.txt',
            'T_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/T.fpkm.txt',
            'T_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/mrna/T.tpm.txt',
        }
        wf.test_api.add_exp_mrna(
            map_dict=map_dict,
            exp_method='Salmon',
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )
        map_dict = {
            'G_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/G.fpkm.txt',
            'G_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/G.tpm.txt',
            'T_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/T.fpkm.txt',
            'T_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191022/Single_exp_make_7164_9741/ExpMake/output/lncrna/T.tpm.txt',
        }
        wf.test_api.add_exp_lncrna(
            map_dict=map_dict,
            exp_method='Salmon',
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_exp_graph(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression')
        map_dict = {
            'sample_box': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_9598_4720/Formation/output/exp_graph/sample_box_data.pk',
            'group_box': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_9598_4720/Formation/output/exp_graph/group_box_data.pk',
            'sample_density': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_9598_4720/Formation/output/exp_graph/sample_density_data.pk',
            'group_density': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_9598_4720/Formation/output/exp_graph/group_density_data.pk',
            'sample_volin': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_9598_4720/Formation/output/exp_graph/sample_volin_data.pk',
            'group_volin': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_9598_4720/Formation/output/exp_graph/group_volin_data.pk',
        }
        wf.test_api.add_exp_graph(
            map_dict=map_dict,
            category='mRNA',
            exp_id='5d8af5d917b2bf2b172ef70f',
            level='T',
            kind='all',
            group_id='5d8b0645e096b94bf4754c44',
            group_dict={'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']},
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )
        map_dict = {
            'sample_box': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_4733_9632/Formation/output/exp_graph/sample_box_data.pk',
            'group_box': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_4733_9632/Formation/output/exp_graph/group_box_data.pk',
            'sample_density': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_4733_9632/Formation/output/exp_graph/sample_density_data.pk',
            'group_density': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_4733_9632/Formation/output/exp_graph/group_density_data.pk',
            'sample_volin': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_4733_9632/Formation/output/exp_graph/sample_volin_data.pk',
            'group_volin': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_4733_9632/Formation/output/exp_graph/group_volin_data.pk',
        }
        wf.test_api.add_exp_graph(
            map_dict=map_dict,
            category='lncRNA',
            exp_id='5d8af5d917b2bf2b172ef70f',
            level='T',
            kind='all',
            group_id='5d8b0645e096b94bf4754c44',
            group_dict={'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']},
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_exp_venn(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression')
        map_dict = {
            'venn': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_9598_4720/Formation/output/exp_venn/venn.txt',
        }
        wf.test_api.add_exp_venn(
            map_dict=map_dict,
            category='mRNA',
            exp_id='5d8af5d917b2bf2b172ef70f',
            level='T',
            kind='all',
            group_id='5d8b0645e096b94bf4754c44',
            group_dict={'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']},
            threshold=1.0,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )
        map_dict = {
            'venn': '/mnt/ilustre/users/sanger-dev/workspace/20190925/Single_formation_4733_9632/Formation/output/exp_venn/venn.txt',
        }
        wf.test_api.add_exp_venn(
            map_dict=map_dict,
            category='lncRNA',
            exp_id='5d8af5d917b2bf2b172ef70f',
            level='T',
            kind='all',
            group_id='5d8b0645e096b94bf4754c44',
            group_dict={'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']},
            threshold=1.0,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_exp_corr(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression')
        map_dict = {
            'corr': '/mnt/ilustre/users/sanger-dev/workspace/20190926/Single_exp_corr_1773_4137/ExpCorr/output/corr.txt',
            'tree': '/mnt/ilustre/users/sanger-dev/workspace/20190926/Single_exp_corr_1773_4137/ExpCorr/output/tree.txt'
        }
        wf.test_api.add_exp_corr(
            map_dict=map_dict,
            library='long',
            level='T',
            kind='all',
            group_id='5d8b0645e096b94bf4754c44',
            group_dict={'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']},
            way='tpm',
            take_mean='False',
            corr_method='pearson',
            take_log='False',
            dist_method='euclidean',
            clus_method='complete',
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_exp_pca(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression')
        map_dict = {
            'evr': '/mnt/ilustre/users/sanger-dev/workspace/20190926/Single_exp_pca_8343_2545/ExpPca/output/explained_variance_ratio.txt',
            'pca': '/mnt/ilustre/users/sanger-dev/workspace/20190926/Single_exp_pca_8343_2545/ExpPca/output/pca.txt',
            'ellipse': '/mnt/ilustre/users/sanger-dev/workspace/20190926/Single_exp_pca_8343_2545/ExpPca/output/ellipse.txt'
        }
        wf.test_api.add_exp_pca(
            map_dict=map_dict,
            library='long',
            level='T',
            kind='all',
            group_id='5d8b0645e096b94bf4754c44',
            group_dict={'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']},
            way='tpm',
            take_mean='False',
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_diff(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'expression_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.expression')
        map_dict = {
            'control_vs_Acute': '/mnt/ilustre/users/sanger-dev/workspace/20190930/Single_diff_exp_1011_6574/DiffExp/output/control_vs_Acute.detail.txt',
            'control_vs_Stable': '/mnt/ilustre/users/sanger-dev/workspace/20190930/Single_diff_exp_1011_6574/DiffExp/output/control_vs_Stable.detail.txt',
            'summary': '/mnt/ilustre/users/sanger-dev/workspace/20190930/Single_diff_exp_1011_6574/DiffExp/output/summary.txt',
            'volcano': '/mnt/ilustre/users/sanger-dev/workspace/20190930/Single_diff_exp_1011_6574/DiffExp/output/volcano.txt',
            'scatter': '/mnt/ilustre/users/sanger-dev/workspace/20190930/Single_diff_exp_1011_6574/DiffExp/output/scatter.txt'
        }
        wf.test_api.add_diff(
            map_dict=map_dict,
            category='mRNA',
            level='T',
            kind='all',
            background='mRNA+lncRNA',
            filter='none',
            threshold='NA',
            group_id='5d8b0645e096b94bf4754c44',
            group_dict={'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']},
            control_id='5d8c8d1a17b2bf43445da3cd',
            stat_type='padjust',
            stat_cutoff=0.05,
            fc=2.0,
            diff_method='DESeq2',
            correct_method='BH',
            way='tpm',
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    # suite.addTests([TestFunction('test_exp_by_level'),
    #                 TestFunction('test_exp_graph'),
    #                 TestFunction('test_exp_venn'),
    #                 TestFunction('test_exp_corr'),
    #                 TestFunction('test_exp_pca'),
    #                 TestFunction('test_diff')])
    suite.addTests([TestFunction('test_exp_graph'), TestFunction('test_exp_venn')])
    unittest.TextTestRunner(verbosity=2).run(suite)
