# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import os
import unittest

import pandas as pd
from bson.son import SON

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Annotation(ApiBase):
    def __init__(self, bind_object):
        super(Annotation, self).__init__(bind_object)

    def add_annotation_go(self, map_dict, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Annotation_go_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'annotationgo', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.create_db_table('annotation_go', [main_dict])
        for anno_type, seq_type in [('T', 'new'), ('G', 'new'), ('T', 'ref'), ('G', 'ref')]:
            self.add_annotation_go_level(map_dict['{}_{}_2'.format(anno_type, seq_type)], anno_type, seq_type, 2,
                                         main_id)
            self.add_annotation_go_detail(map_dict['{}_{}_2'.format(anno_type, seq_type)], anno_type, seq_type, 2,
                                          main_id)
            self.add_annotation_go_detail(map_dict['{}_{}_3'.format(anno_type, seq_type)], anno_type, seq_type, 3,
                                          main_id)
            self.add_annotation_go_detail(map_dict['{}_{}_4'.format(anno_type, seq_type)], anno_type, seq_type, 4,
                                          main_id)
            self.add_annotation_go_graph(map_dict['{}_{}_2'.format(anno_type, seq_type)], anno_type, seq_type, 2,
                                         main_id)
            self.add_annotation_go_graph(map_dict['{}_{}_3'.format(anno_type, seq_type)], anno_type, seq_type, 3,
                                         main_id)
            self.add_annotation_go_graph(map_dict['{}_{}_4'.format(anno_type, seq_type)], anno_type, seq_type, 4,
                                         main_id)
            self.add_annotation_go_list(map_dict['{}_{}_gos'.format(anno_type, seq_type)], anno_type, seq_type, main_id)
        else:
            self.add_annotation_go_all(map_dict['T_all_2'], 'T', 'all', 2, main_id)
            self.add_annotation_go_all(map_dict['T_all_3'], 'T', 'all', 3, main_id)
            self.add_annotation_go_all(map_dict['T_all_4'], 'T', 'all', 4, main_id)
            self.add_annotation_go_all(map_dict['G_all_2'], 'G', 'all', 2, main_id)
            self.add_annotation_go_all(map_dict['G_all_3'], 'G', 'all', 3, main_id)
            self.add_annotation_go_all(map_dict['G_all_4'], 'G', 'all', 4, main_id)
        self.update_db_record('annotation_go', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_annotation_go_level(self, level_path, anno_type, seq_type, level, go_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        # level in [2, 3, 4]
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'anno_type': anno_type, 'seq_type': seq_type, 'level': level}
        ))
        data = list()
        for line in open(level_path).readlines()[1:]:
            eles = line.strip().split('\t')
            docx = {
                'go_id': go_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'level': level,
                'parent_name': eles[0],
                'term_type': eles[1],
                'go': eles[2],
                'num': int(eles[3]),
                'percent': round(float(eles[4]), 4)
            }
            data.append(docx)
        else:
            self.create_db_table('annotation_go_level', data)

    def add_annotation_go_detail(self, level_path, anno_type, seq_type, level, go_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        # level in [2, 3, 4]
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'anno_type': anno_type, 'seq_type': seq_type, 'level': level}
        ))
        data = list()
        for line in open(level_path).readlines()[1:]:
            eles = line.strip().split('\t')
            docx = {
                'go_id': go_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'level': level,
                'goterm': eles[0],
                'goterm_2': eles[1],
                'goid_2': eles[2],
                'seq_number': int(eles[-3]),
                'percent': round(float(eles[-2]), 4)
            }
            if level == 2:
                docx['seq_list'] = eles[-1]
            if level >= 3:
                docx['goterm_3'] = eles[3]
                docx['goid_3'] = eles[4]
            if level == 4:
                docx['goterm_4'] = eles[5]
                docx['goid_4'] = eles[6]
            data.append(docx)
        else:
            self.create_db_table('annotation_go_detail', data)

    def add_annotation_go_graph(self, level_path, anno_type, seq_type, level, go_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        # level in [2, 3, 4]
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'anno_type': anno_type, 'seq_type': seq_type, 'level': level}
        ))
        data = list()
        term_dict = dict()
        lines = open(level_path).readlines()
        for i in range(1, len(lines)):
            eles = lines[i].strip().split('\t')
            if level == 2:
                term_type = eles[0]
                go_term = eles[1]
                if go_term not in term_dict:
                    term_dict[go_term] = list()
                    term_dict[go_term].append(i)
                else:
                    term_dict[go_term].append(i)
            if level == 3:
                term_type = eles[0]
                go_term = eles[3]
                if go_term not in term_dict:
                    term_dict[go_term] = list()
                    term_dict[go_term].append(i)
                else:
                    term_dict[go_term].append(i)
            if level == 4:
                term_type = eles[0]
                go_term = eles[5]
                if go_term not in term_dict:
                    term_dict[go_term] = list()
                    term_dict[go_term].append(i)
                else:
                    term_dict[go_term].append(i)
        for go_term in term_dict:
            for j in term_dict[go_term]:
                eles = lines[j].strip().split('\t')
                term_type = eles[0]
            docx = {
                'go_id': go_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'level': level,
                'term_type': term_type,
                'go_term': go_term,
                'seq_number': eles[-3],
                'percent': eles[-2]
            }
            data.append(docx)
        else:
            self.create_db_table('annotation_go_graph', data)

    def add_annotation_go_list(self, list_path, anno_type, seq_type, go_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'anno_type': anno_type, 'seq_type': seq_type}
        ))
        data = list()
        for line in open(list_path):
            eles = line.strip().split('\t')
            docx = {
                'go_id': go_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'gene_id': eles[0],
                'gos_list': eles[1],
            }
            data.append(docx)
        else:
            self.create_db_table('annotation_go_list', data)

    def add_annotation_go_all(self, level_path, anno_type, seq_type, level, go_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        # level in [2, 3, 4]
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
            {'anno_type': anno_type, 'seq_type': seq_type, 'level': level}
        ))
        graph_data, detail_data, query_ids = list(), list(), list()
        func_dict, term_dict = dict(), dict()
        lines = open(level_path).readlines()
        for line in lines[1:]:
            eles = line.strip().split('\t')
            func = '{}|||{}|||{}'.format(eles[0], eles[1], eles[2])
            term = '{}|||{}'.format(eles[0], eles[1])
            if level == 3:
                func += '|||{}|||{}'.format(eles[3], eles[4])
                term = '{}|||{}'.format(eles[0], eles[3])
            if level == 4:
                func += '|||{}|||{}|||{}|||{}'.format(eles[3], eles[4], eles[5], eles[6])
                term = '{}|||{}'.format(eles[0], eles[5])
            func_dict[func] = eles[-1].split(';')
            if term not in term_dict:
                term_dict[term] = set(eles[-1].split(';'))
            else:
                for seq_id in eles[-1].split(';'):
                    if seq_id not in term_dict[term]:
                        term_dict[term].add(seq_id)
            query_ids.extend(eles[-1].split(';'))
        else:
            query_ids = list(set(query_ids))
        for term in sorted(term_dict):
            terms = term.split('|||')
            seq_list = term_dict[term]
            percent = float(len(seq_list)) / len(query_ids)
            data = dict([
                ('go_id', go_id),
                ('anno_type', anno_type),
                ('seq_type', seq_type),
                ('level', level),
                ('term_type', terms[0]),
                ('go_term', terms[1]),
                ('seq_number', len(seq_list)),
                ('percent', round(percent, 4))
            ])
            data = SON(data)
            graph_data.append(data)
        else:
            self.create_db_table('annotation_go_graph', graph_data)
        for func in func_dict:
            funcs = func.split('|||')
            data = [
                ('go_id', go_id),
                ('anno_type', anno_type),
                ('seq_type', seq_type),
                ('level', level),
                ('goterm', funcs[0]),
                ('goterm_2', funcs[1]),
                ('goid_2', funcs[2])
            ]
            if level >= 3:
                data.append(('goterm_3', funcs[3]))
                data.append(('goid_3', funcs[4]))
            if level == 4:
                data.append(('goterm_4', funcs[5]))
                data.append(('goid_4', funcs[6]))
            seq_list = func_dict[func]
            percent = float(len(func_dict[func])) / len(query_ids)
            data.append(('seq_number', len(seq_list)))
            data.append(('percent', round(percent, 4)))
            if level == 2:
                data.append(('seq_list', ';'.join(seq_list)))
            data = SON(data)
            detail_data.append(data)
        else:
            self.create_db_table('annotation_go_detail', detail_data)

    def add_annotation_kegg(self, map_dict, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Annotation_kegg_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'annotationkegg', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.create_db_table('annotation_kegg', [main_dict])
        for anno_type, seq_type in [('T', 'new'), ('G', 'new'), ('T', 'ref'), ('G', 'ref')]:
            self.add_annotation_kegg_categories(map_dict['{}_{}_c'.format(anno_type, seq_type)], anno_type, seq_type,
                                                main_id)
            self.add_annotation_kegg_level(map_dict['{}_{}_l'.format(anno_type, seq_type)], anno_type, seq_type,
                                           main_id)
            self.add_annotation_kegg_pic(map_dict['{}_{}_l'.format(anno_type, seq_type)],
                                         map_dict['{}_{}_p'.format(anno_type, seq_type)], anno_type, seq_type, main_id)
            self.add_annotation_kegg_table(map_dict['{}_{}_t'.format(anno_type, seq_type)], anno_type, seq_type,
                                           main_id)
        else:
            self.add_annotation_kegg_categories_all(map_dict['T_ref_c'], map_dict['T_new_c'], 'T', 'all', main_id)
            self.add_annotation_kegg_categories_all(map_dict['G_ref_c'], map_dict['G_new_c'], 'G', 'all', main_id)
            self.add_annotation_kegg_level(map_dict['T_all_l'], 'T', 'all', main_id)
            self.add_annotation_kegg_level(map_dict['G_all_l'], 'G', 'all', main_id)
            self.add_annotation_kegg_pic(map_dict['T_all_l'], map_dict['T_all_p'], 'T', 'all', main_id)
            self.add_annotation_kegg_pic(map_dict['G_all_l'], map_dict['G_all_p'], 'G', 'all', main_id)
        self.update_db_record('annotation_kegg', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_annotation_kegg_categories(self, categories_path, anno_type, seq_type, kegg_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        data = list()
        first_type = list()
        for line in open(categories_path):
            eles = line.strip('\n').split('\t')
            type_abr = ''.join([x[0] for x in eles[0].split(' ')])
            first_type.append(type_abr)
            docx = {
                'kegg_id': kegg_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'first_category': eles[0],
                'second_category': eles[1],
                'num': int(eles[2]),
                'seq_list': eles[3]
            }
            data.append(docx)
        else:
            if data:
                self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                    {'anno_type': anno_type, 'seq_type': seq_type}
                ))
                self.create_db_table('annotation_kegg_categories', data)
                self.update_db_record('annotation_kegg', kegg_id, categories=list(set(first_type)))

    def add_annotation_kegg_level(self, level_path, anno_type, seq_type, kegg_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        data = list()
        lines = open(level_path).readlines()
        for line in lines[1:]:
            eles = line.strip('\n').split('\t')
            docx = {
                'kegg_id': kegg_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'pathway_id': eles[0],
                'first_category': eles[1],
                'second_category': eles[2],
                'pathway_definition': eles[3],
                'number_of_seqs': int(eles[4]),
                'seq_list': eles[5],
                'hyperlink': eles[-1]
            }
            data.append(docx)
        else:
            if data:
                self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                    {'anno_type': anno_type, 'seq_type': seq_type}
                ))
                self.create_db_table('annotation_kegg_level', data)

    def add_annotation_kegg_table(self, table_path, anno_type, seq_type, kegg_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        data = list()
        lines = open(table_path).readlines()
        for line in lines[1:]:
            eles = line.strip('\n').split('\t')
            docx = {
                'kegg_id': kegg_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'transcript_id': eles[0],
                'ko_id': eles[1],
                'ko_name': eles[2],
                'hyperlink': eles[3],
                'paths': eles[4]
            }
            data.append(docx)
        else:
            if data:
                self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                    {'anno_type': anno_type, 'seq_type': seq_type}
                ))
                self.create_db_table('annotation_kegg_table', data)

    def add_annotation_kegg_pic(self, level_path, pic_path, anno_type, seq_type, kegg_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        data = list()
        lines = open(level_path).readlines()
        for line in lines[1:]:
            eles = line.strip('\n').split('\t')
            pathway_id = eles[0]
            mark_path = os.path.join(pic_path, '{}.html.mark'.format(pathway_id))
            if os.path.exists(mark_path):
                for mark_line in open(mark_path):
                    mark_eles = mark_line.strip().split('\t')
                    if len(mark_eles) == 8:
                        [png, shape, bg_colors, fg_colors, coords, title, kos, href] = mark_eles
                        title = title.replace('\\n', '\n')
                        docx = {
                            'kegg_id': kegg_id,
                            'anno_type': anno_type,
                            'seq_type': seq_type,
                            'pathway_id': pathway_id,
                            'shape': shape,
                            'bg_colors': bg_colors,
                            'fg_colors': fg_colors,
                            'coords': coords,
                            'href': href,
                            'kos': kos,
                            'title': title
                        }
                        if bg_colors != '' and len(bg_colors.split(',')) > 0:
                            docx.update({'bg_type': len(bg_colors.split(','))})
                        if fg_colors != '' and len(fg_colors.split(',')) > 0:
                            docx.update({'fg_type': len(fg_colors.split(','))})
                        data.append(docx)
        else:
            if data:
                self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                    {'anno_type': anno_type, 'seq_type': seq_type}
                ))
                self.create_db_table('annotation_kegg_pic', data)

    def add_annotation_kegg_categories_all(self, ref_path, new_path, anno_type, seq_type, kegg_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        data = list()
        first_type = list()
        cate_dict = dict()
        for line in open(ref_path):
            eles = line.strip('\n').split('\t')
            type_abr = ''.join([x[0] for x in eles[0].split(' ')])
            first_type.append(type_abr)
            if eles[0] == 'Metabolism' and eles[1] == 'Global and overview maps':
                pass
            else:
                cate_key = '{}|||{}'.format(eles[0], eles[1])
                cate_dict[cate_key] = eles[3].split(';')
        for line in open(new_path):
            eles = line.strip('\n').split('\t')
            type_abr = ''.join([x[0] for x in eles[0].split(' ')])
            first_type.append(type_abr)
            if eles[0] == 'Metabolism' and eles[1] == 'Global and overview maps':
                pass
            else:
                cate_key = '{}|||{}'.format(eles[0], eles[1])
                if cate_key in cate_dict:
                    for seq_id in eles[3].split(';'):
                        if seq_id not in cate_dict[cate_key]:
                            cate_dict[cate_key].append(seq_id)
                else:
                    cate_dict[cate_key] = eles[3].split(';')
        for category in ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
                         'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']:
            for cate_key in cate_dict:
                first_category, second_category = cate_key.split('|||')
                if category == first_category:
                    num = len(cate_dict[cate_key])
                    seq_list = ';'.join(cate_dict[cate_key])
                    docx = {
                        'kegg_id': kegg_id,
                        'anno_type': anno_type,
                        'seq_type': seq_type,
                        'first_category': first_category,
                        'second_category': second_category,
                        'num': num,
                        'seq_list': seq_list
                    }
                    data.append(docx)
        else:
            if data:
                self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                    {'anno_type': anno_type, 'seq_type': seq_type}
                ))
                self.create_db_table('annotation_kegg_categories', data)
                self.update_db_record('annotation_kegg', kegg_id, categories=list(set(first_type)))

    def add_annotation_cog(self, map_dict, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Annotation_cog_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'annotationcog', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.create_db_table('annotation_cog', [main_dict])
        self.add_annotation_cog_detail(map_dict['T_new'], 'T', 'new', main_id)
        self.add_annotation_cog_detail(map_dict['G_new'], 'G', 'new', main_id)
        self.add_annotation_cog_detail(map_dict['T_ref'], 'T', 'ref', main_id)
        self.add_annotation_cog_detail(map_dict['G_ref'], 'G', 'ref', main_id)
        self.add_annotation_cog_detail_all(map_dict['T_ref'], map_dict['T_new'], 'T', 'all', main_id)
        self.add_annotation_cog_detail_all(map_dict['G_new'], map_dict['G_ref'], 'G', 'all', main_id)
        self.update_db_record('annotation_cog', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_annotation_cog_detail(self, sum_path, anno_type, seq_type, cog_id):
        # anno_type in ['G', 'T']
        # seq_type in ['ref', 'new']
        data = list()
        lines = open(sum_path).readlines()
        for line in lines[1:]:
            eles = line.strip().split('\t')
            docx = {
                'cog_id': cog_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'type': eles[0],
                'function_categories': '[{}] {}'.format(eles[2], eles[1]),
                'cog': int(eles[3]),
                'cog_list': eles[4] if len(eles) > 4 else None
            }
            data.append(docx)
        else:
            if data:
                self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                    {'anno_type': anno_type, 'seq_type': seq_type}
                ))
                self.create_db_table('annotation_cog_detail', data)

    def add_annotation_cog_detail_all(self, ref_path, new_path, anno_type, seq_type, cog_id):
        # anno_type in ['G', 'T']
        # seq_type in ['all']
        first = [
            'INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING',
            'METABOLISM', 'POORLY CHARACTERIZED'
        ]
        func_type = {
            'INFORMATION STORAGE AND PROCESSING': sorted(['J', 'A', 'K', 'L', 'B']),
            'CELLULAR PROCESSES AND SIGNALING': sorted(['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O']),
            'METABOLISM': sorted(['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q']),
            'POORLY CHARACTERIZED': sorted(['R', 'S']),
        }
        func_decs = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'A': 'RNA processing and modification', 'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'B': 'Chromatin structure and dynamics',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular structures',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'C': 'Energy production and conversion', 'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only', 'S': 'Function unknown'
        }
        data = list()
        func_dict = dict()
        for line in open(ref_path).readlines()[1:]:
            eles = line.strip().split('\t')
            func = eles[2]
            func_dict[func] = eles[4].split(';') if len(eles) > 4 else list()
        for line in open(new_path).readlines()[1:]:
            eles = line.strip().split('\t')
            func = eles[2]
            if func in func_dict:
                func_dict[func].extend(eles[4].split(';') if len(eles) > 4 else list())
            else:
                func_dict[func] = eles[4].split(';') if len(eles) > 4 else list()
        for function in func_dict:
            categories = func_decs[function]
            function_categories = '[{}] {}'.format(function, categories)
            seq_ids = [seq_id for seq_id in list(set(func_dict[function])) if seq_id]
            for key, value in func_type.items():
                if function in value:
                    cog_type = key
            docx = {
                'cog_id': cog_id,
                'anno_type': anno_type,
                'seq_type': seq_type,
                'type': cog_type,
                'function_categories': function_categories,
                'cog': len(seq_ids),
                'cog_list': ';'.join(seq_ids),
            }
            data.append(docx)
        else:
            if data:
                self.bind_object.logger.info('start inserting documents into mongo with ({})'.format(
                    {'anno_type': anno_type, 'seq_type': seq_type}
                ))
                self.create_db_table('annotation_cog_detail', data)

    def add_annotation_stat(self, map_dict, task_id, project_sn):
        time_now = datetime.datetime.now()
        name = 'Annotation_stat_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        params = json.dumps({'task_id': task_id, 'submit_location': 'annotationstat', 'task_type': 2}, sort_keys=True)
        main_dict = {
            'task_id': task_id,
            'project_sn': project_sn,
            'name': name,
            'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'desc': 'Workflow result',
            'params': params,
            'status': 'start',
            'version': 'v1'
        }
        main_id = self.create_db_table('annotation_stat', [main_dict])
        self.add_annotation_stat_detail(map_dict, main_id, task_id)
        self.update_db_record('annotation_stat', main_id, insert_dict={'main_id': main_id, 'status': 'end'})

    def add_annotation_stat_detail(self, map_dict, stat_id, task_id):
        self.bind_object.logger.info('start inserting documents into mongo with ({})'.format({'kind': 'all'}))
        all_t_id_set = set()
        all_g_id_set = set()
        for line in open(map_dict['all_t2g']):
            transcript_id, gene_id = line.strip().split('\t')[:2]
            all_t_id_set.add(transcript_id)
            all_g_id_set.add(gene_id)
        ref_t_id_set = set()
        ref_g_id_set = set()
        for line in open(map_dict['ref_t2g']):
            transcript_id, gene_id = line.strip().split('\t')[:2]
            ref_t_id_set.add(transcript_id)
            ref_g_id_set.add(gene_id)
        new_t_id_set = set()
        new_g_id_set = set()
        for line in open(map_dict['new_t2g']):
            transcript_id, gene_id = line.strip().split('\t')[:2]
            new_t_id_set.add(transcript_id)
            new_g_id_set.add(gene_id)

        t_dict = dict()
        g_dict = dict()
        for level, x_dict in (('T', t_dict), ('G', g_dict)):
            for db in ('go', 'kegg', 'cog', 'nr', 'swissprot', 'pfam'):
                x_dict[db] = set(line.strip() for line in open(map_dict['{}_{}'.format(level, db)]))

        # deprecated
        # t_count_df = pd.read_table(map_dict['T_count'], index_col=0)
        # exp_t_id_set = set(t_count_df[t_count_df.sum(axis=1) > 0].index)
        # g_count_df = pd.read_table(map_dict['G_count'], index_col=0)
        # exp_g_id_set = set(g_count_df[g_count_df.sum(axis=1) > 0].index)

        samples = [document['old_name'] for document in
                   self.db['specimen'].find({'task_id': task_id, 'library': 'long'})]
        t_exp_df = pd.read_table(map_dict['T_exp'])
        t_exp_df = t_exp_df.reindex(['transcript_id', 'category'] + samples, axis=1)
        t_exp_df = t_exp_df[t_exp_df['category'] == 'mRNA'].set_index('transcript_id').drop('category', axis=1)
        exp_t_id_set = set(t_exp_df[t_exp_df.sum(axis=1) > 0].index)
        g_exp_df = pd.read_table(map_dict['G_exp'])
        g_exp_df = g_exp_df.reindex(['gene_id', 'category'] + samples, axis=1)
        g_exp_df = g_exp_df[g_exp_df['category'] == 'mRNA'].set_index('gene_id').drop('category', axis=1)
        exp_g_id_set = set(g_exp_df[g_exp_df.sum(axis=1) > 0].index)

        docs = list()
        for kind in ('all', 'ref', 'new'):
            for idx in ('go', 'kegg', 'cog', 'nr', 'swissprot', 'pfam', 'total_anno', 'total'):
                docs.append({'type': idx, 'exp_g_num': 0, 'exp_t_num': 0, 'g_num': 0, 't_num': 0, 'kind': kind})
        else:
            num_df = pd.DataFrame(docs).set_index(['type', 'kind'])

        zero_indexes = list()
        for kind, kind_t_id_set, kind_g_id_set in (('all', all_t_id_set, all_g_id_set),
                                                   ('ref', ref_t_id_set, ref_g_id_set),
                                                   ('new', new_t_id_set, new_g_id_set)):
            t_total_anno_set, g_total_anno_set = set(), set()
            for db in ('go', 'kegg', 'cog', 'nr', 'swissprot', 'pfam'):
                db_t_id_set = t_dict[db] & kind_t_id_set
                db_g_id_set = g_dict[db] & kind_g_id_set
                exp_db_t_id_set = db_t_id_set & exp_t_id_set
                exp_db_g_id_set = db_g_id_set & exp_g_id_set
                num_df.loc[(db, kind), 't_num'] = len(db_t_id_set)
                num_df.loc[(db, kind), 'g_num'] = len(db_g_id_set)
                num_df.loc[(db, kind), 'exp_t_num'] = len(exp_db_t_id_set)
                num_df.loc[(db, kind), 'exp_g_num'] = len(exp_db_g_id_set)
                t_total_anno_set.update(db_t_id_set)
                g_total_anno_set.update(db_g_id_set)
            else:
                num_df.loc[('total_anno', kind), 't_num'] = len(t_total_anno_set)
                num_df.loc[('total_anno', kind), 'g_num'] = len(g_total_anno_set)
                num_df.loc[('total_anno', kind), 'exp_t_num'] = len(t_total_anno_set & exp_t_id_set)
                num_df.loc[('total_anno', kind), 'exp_g_num'] = len(g_total_anno_set & exp_g_id_set)
                num_df.loc[('total', kind), 't_num'] = len(kind_t_id_set)
                num_df.loc[('total', kind), 'g_num'] = len(kind_g_id_set)
                num_df.loc[('total', kind), 'exp_t_num'] = len(kind_t_id_set & exp_t_id_set)
                num_df.loc[('total', kind), 'exp_g_num'] = len(kind_g_id_set & exp_g_id_set)
                for col in ('t_num', 'g_num', 'exp_t_num', 'exp_g_num'):
                    if not num_df.loc[('total', kind), col]:
                        num_df.loc[('total', kind), col] += 1
                        zero_indexes.append((kind, col))
        else:
            all_per_df = (num_df.loc[(slice(None), 'all'),] / num_df.loc[('total', 'all')] * 100).rename(
                lambda x: x.replace('num', 'per'), axis=1)
            ref_per_df = (num_df.loc[(slice(None), 'ref'),] / num_df.loc[('total', 'ref')] * 100).rename(
                lambda x: x.replace('num', 'per'), axis=1)
            new_per_df = (num_df.loc[(slice(None), 'new'),] / num_df.loc[('total', 'new')] * 100).rename(
                lambda x: x.replace('num', 'per'), axis=1)
            per_df = pd.concat([all_per_df, ref_per_df, new_per_df])
            for kind, col in zero_indexes:
                num_df.loc[('total', kind), col] = 0
            df = pd.concat([num_df, per_df], axis=1)

        # for t_id in all_t_id_set:
        #     is_anno = False
        #     columns = ['t_num', 'exp_t_num'] if t_id in exp_t_id_set else ['t_num']
        #     for col in columns:
        #         for db in ('go', 'kegg', 'cog', 'nr', 'swissprot', 'pfam'):
        #             if t_id in t_dict[db]:
        #                 is_anno = True
        #                 df.loc[(db, 'all'), col] += 1
        #                 if t_id in ref_t_id_set:
        #                     df.loc[(db, 'ref'), col] += 1
        #                 if t_id in new_t_id_set:
        #                     df.loc[(db, 'new'), col] += 1
        #         else:
        #             if is_anno:
        #                 df.loc[('total_anno', 'all'), col] += 1
        #                 if t_id in ref_t_id_set:
        #                     df.loc[('total_anno', 'ref'), col] += 1
        #                 if t_id in new_t_id_set:
        #                     df.loc[('total_anno', 'new'), col] += 1
        #             df.loc[('total', 'all'), col] += 1
        #             if t_id in ref_t_id_set:
        #                 df.loc[('total', 'ref'), col] += 1
        #             if t_id in new_t_id_set:
        #                 df.loc[('total', 'new'), col] += 1
        #
        # for g_id in all_g_id_set:
        #     is_anno = False
        #     columns = ['g_num', 'exp_g_num'] if g_id in exp_g_id_set else ['g_num']
        #     for col in columns:
        #         for db in ('go', 'kegg', 'cog', 'nr', 'swissprot', 'pfam'):
        #             if g_id in g_dict[db]:
        #                 is_anno = True
        #                 df.loc[(db, 'all'), col] += 1
        #                 if g_id in ref_g_id_set:
        #                     df.loc[(db, 'ref'), col] += 1
        #                 if g_id in new_g_id_set:
        #                     df.loc[(db, 'new'), col] += 1
        #         else:
        #             if is_anno:
        #                 df.loc[('total_anno', 'all'), col] += 1
        #                 if g_id in ref_g_id_set:
        #                     df.loc[('total_anno', 'ref'), col] += 1
        #                 if g_id in new_g_id_set:
        #                     df.loc[('total_anno', 'new'), col] += 1
        #             df.loc[('total', 'all'), col] += 1
        #             if g_id in ref_g_id_set:
        #                 df.loc[('total', 'ref'), col] += 1
        #             if g_id in new_g_id_set:
        #                 df.loc[('total', 'new'), col] += 1

        type_dict = {'go': 'GO', 'kegg': 'KEGG', 'cog': 'COG', 'nr': 'NR', 'swissprot': 'Swiss-Prot', 'pfam': 'Pfam',
                     'total_anno': 'Total_anno', 'total': 'Total'}

        data = df.reset_index().to_dict('r')
        for i, document in enumerate(data):
            document['type'] = type_dict[document['type']]
            data[i] = document

        self.create_db_table('annotation_stat_detail', data, {'stat_id': stat_id})


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test_go(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.annotation')
        map_dict = {
            'T_new_2': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_lev2_tran.stat.xls',
            'T_new_3': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_lev3_tran.stat.xls',
            'T_new_4': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_lev4_tran.stat.xls',
            'T_new_gos': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_list_tran.xls',
            'G_new_2': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_lev2_gene.stat.xls',
            'G_new_3': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_lev3_gene.stat.xls',
            'G_new_4': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_lev4_gene.stat.xls',
            'G_new_gos': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/go/go_list_gene.xls',

            'T_ref_2': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_lev2_tran.stat.xls',
            'T_ref_3': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_lev3_tran.stat.xls',
            'T_ref_4': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_lev4_tran.stat.xls',
            'T_ref_gos': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_list_tran.xls',
            'G_ref_2': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_lev2_gene.stat.xls',
            'G_ref_3': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_lev3_gene.stat.xls',
            'G_ref_4': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_lev4_gene.stat.xls',
            'G_ref_gos': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/go/go_list_gene.xls',

            'T_all_2': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/go/go_lev2_tran.stat.xls',
            'T_all_3': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/go/go_lev3_tran.stat.xls',
            'T_all_4': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/go/go_lev4_tran.stat.xls',
            'G_all_2': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/go/go_lev2_gene.stat.xls',
            'G_all_3': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/go/go_lev3_gene.stat.xls',
            'G_all_4': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/go/go_lev4_gene.stat.xls',
        }
        wf.test_api.add_annotation_go(
            map_dict=map_dict,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_kegg(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.annotation')
        map_dict = {
            'T_new_c': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_layer_tran.xls',
            'T_new_l': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_pathway_tran.xls',
            'T_new_p': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_pathway_tran_dir',
            'T_new_t': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_gene_tran.xls',
            'G_new_c': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_layer_gene.xls',
            'G_new_l': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_pathway_gene.xls',
            'G_new_p': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_pathway_gene_dir',
            'G_new_t': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/kegg/kegg_gene_gene.xls',

            'T_ref_c': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_layer_tran.xls',
            'T_ref_l': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_pathway_tran.xls',
            'T_ref_p': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_pathway_tran_dir',
            'T_ref_t': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_gene_tran.xls',
            'G_ref_c': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_layer_gene.xls',
            'G_ref_l': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_pathway_gene.xls',
            'G_ref_p': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_pathway_gene_dir',
            'G_ref_t': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/kegg/kegg_gene_gene.xls',

            'T_all_l': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/kegg/kegg_pathway_tran.xls',
            'T_all_p': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/kegg/kegg_pathway_tran_dir',
            'G_all_l': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/kegg/kegg_pathway_gene.xls',
            'G_all_p': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/allannot_class/kegg/kegg_pathway_gene_dir'
        }
        wf.test_api.add_annotation_kegg(
            map_dict=map_dict,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_cog(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.annotation')
        map_dict = {
            'T_new': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/cog/summary.T.tsv',
            'G_new': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/newannot_class/cog/summary.G.tsv',
            'T_ref': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/cog/summary.T.tsv',
            'G_ref': '/mnt/ilustre/users/sanger-dev/workspace/20190924/WholeTranscriptome_workflow_1623_1285/Annotation/output/refannot_class/cog/summary.G.tsv'
        }
        wf.test_api.add_annotation_cog(
            map_dict=map_dict,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )

    def test_stat(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('whole_transcriptome.annotation')
        map_dict = {
            'all_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/all_tran2gene.txt',
            'ref_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/refannot_class/all_tran2gene.txt',
            'new_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/newannot_class/all_tran2gene.txt',

            'T_go': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/go/go_venn_tran.txt',
            'T_kegg': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/kegg/kegg_venn_tran.txt',
            'T_cog': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/cog/cog_venn_tran.txt',
            'T_nr': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/nr/nr_venn_tran.txt',
            'T_swissprot': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/swissprot/swissprot_venn_tran.txt',
            'T_pfam': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/pfam/pfam_venn_tran.txt',
            'G_go': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/go/go_venn_gene.txt',
            'G_kegg': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/kegg/kegg_venn_gene.txt',
            'G_cog': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/cog/cog_venn_gene.txt',
            'G_nr': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/nr/nr_venn_gene.txt',
            'G_swissprot': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/swissprot/swissprot_venn_gene.txt',
            'G_pfam': '/mnt/ilustre/users/sanger-dev/workspace/20191014/Longrna_workflow_6700_7886/Annotation/output/allannot_class/pfam/pfam_venn_gene.txt',

            'T_count': '/mnt/ilustre/users/sanger-dev/workspace/20191015/Single_exp_make_6370_6765/ExpMake/output/count/T.reads.txt',
            'G_count': '/mnt/ilustre/users/sanger-dev/workspace/20191015/Single_exp_make_6370_6765/ExpMake/output/count/G.reads.txt'
        }
        wf.test_api.add_annotation_stat(
            map_dict=map_dict,
            task_id='whole_transcriptome',
            project_sn='whole_transcriptome'
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_go')])
    unittest.TextTestRunner(verbosity=2).run(suite)
