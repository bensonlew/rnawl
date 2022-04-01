# -*- coding: utf-8 -*-
# __author__ = 'liubinxu,qinjincheng'

from mbio.api.database.ref_rna_v2.api_base import ApiBase
from biocluster.api.database.base import report_check
import datetime
import os
import json
import re
from bson.objectid import ObjectId
from bson.son import SON
import pandas as pd
import types
import gridfs
from biocluster.config import Config
import unittest

class Annotation(ApiBase):
    def __init__(self, bind_object):
        super(Annotation, self).__init__(bind_object)
        self.task_id = self.bind_object.sheet.id
        self.project_sn = self.bind_object.sheet.project_sn
        self.t2g_dict = dict()
        self.t2r_dict = dict()
        self.has_new = True
        self.annot_type = 'origin'
        self.time_now = datetime.datetime.now()
        # the following attributes were deprecated
        self.result_dir = str()
        self.result_file = dict()

    @report_check
    def run(self, annot_merge_output_dir, params_dict, taxonomy, exp_level):
        self.set_results(annot_merge_output_dir)
        if self.has_new:
            self.set_relation(self.results['ref']['t2g2r2l2p'], self.results['new']['t2g2r2l2p'])
        else:
            self.set_relation(self.results['ref']['t2g2r2l2p'])
        self.params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        self.taxonomy = taxonomy
        self.exp_level = exp_level

        # stat
        self.remove_table_by_main_record(
            main_table='sg_annotation_stat',
            detail_table='sg_annotation_stat_detail',
            detail_table_key='stat_id',
            task_id=self.task_id,
            type=self.annot_type
        )
        stat_id = self.add_annotation_stat()
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.results['ref']['stat'], root_dir=self.ref_dir, seq_type ='ref', exp_level=exp_level)
        if self.has_new:
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.results['new']['stat'], root_dir=self.new_dir, seq_type ='new', exp_level=exp_level)
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.results['all']['stat'], root_dir=self.all_dir, seq_type ='all', exp_level=exp_level)
        self.update_db_record('sg_annotation_stat', stat_id, status='end', main_id=stat_id)

        # query
        self.remove_table_by_main_record(
            main_table='sg_annotation_query',
            detail_table='sg_annotation_query_detail',
            detail_table_key='query_id',
            task_id=self.task_id,
            type=self.annot_type
        )
        query_id = self.add_annotation_query(stat_id)
        self.add_annotation_query_detail(query_id=query_id, query_path=self.results['ref']['query'], seq_type ='ref', exp_level=exp_level)
        if self.has_new:
            self.add_annotation_query_detail(query_id=query_id, query_path=self.results['new']['query'], seq_type ='new', exp_level=exp_level)
        self.update_db_record('sg_annotation_query',query_id, status='end', main_id=query_id)

        # cog
        self.remove_table_by_main_record(
            main_table='sg_annotation_cog',
            detail_table='sg_annotation_cog_detail',
            detail_table_key='cog_id',
            task_id=self.task_id,
            type=self.annot_type
        )
        self.cog_params = json.dumps(
            dict([(k, params_dict.get(k, None)) for k in ('cog_evalue', 'cog_similarity', 'cog_identity')]),
            sort_keys=True, separators=(',', ':')
        )
        cog_id = self.add_annotation_cog()
        if self.has_new:
            if exp_level.lower() == 'transcript':
                self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.results['new']['cog']['txpt'], seq_type='new', anno_type='T')
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.results['new']['cog']['gene'], seq_type='new', anno_type='G')
        if exp_level.lower() == 'transcript':
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.results['ref']['cog']['txpt'], seq_type='ref', anno_type='T')
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.results['ref']['cog']['gene'], seq_type='ref', anno_type='G')
        if self.has_new:
            if exp_level.lower() == 'transcript':
                self.add_annotation_cog_detail_all(cog_id=cog_id, ref_cog_path=self.results['ref']['cog']['txpt'], new_cog_path=self.results['new']['cog']['txpt'], seq_type='all', anno_type='T')
            self.add_annotation_cog_detail_all(cog_id=cog_id, ref_cog_path=self.results['ref']['cog']['gene'], new_cog_path=self.results['new']['cog']['gene'], seq_type='all', anno_type='G')
        self.update_db_record('sg_annotation_cog', cog_id, status='end', main_id=cog_id)

        # go
        go_detail_tables = ['sg_annotation_go_detail', 'sg_annotation_go_graph', 'sg_annotation_go_level', 'sg_annotation_go_list']
        self.remove_table_by_main_record(
            main_table='sg_annotation_go',
            detail_table=go_detail_tables,
            detail_table_key='go_id',
            task_id=self.task_id,
            type=self.annot_type
        )
        self.go_params = json.dumps(
            dict([(k, params_dict.get(k, None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')]),
            sort_keys=True, separators=(',', ':')
        )
        go_id = self.add_annotation_go()
        if self.has_new:
            seq_type = 'new'
            if exp_level.lower() == 'transcript':
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type='T', level=2, level_path=self.results['new']['go']['txpt']['level2'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='T', level=2, level_path=self.results['new']['go']['txpt']['level2'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='T', level=3, level_path=self.results['new']['go']['txpt']['level3'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='T', level=4, level_path=self.results['new']['go']['txpt']['level4'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='T', level=2, level_path=self.results['new']['go']['txpt']['level2'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='T', level=3, level_path=self.results['new']['go']['txpt']['level3'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='T', level=4, level_path=self.results['new']['go']['txpt']['level4'])
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type='T', list_path=self.results['new']['go']['txpt']['list'])
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type='G', level=2, level_path=self.results['new']['go']['gene']['level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='G', level=2, level_path=self.results['new']['go']['gene']['level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='G', level=3, level_path=self.results['new']['go']['gene']['level3'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='G', level=4, level_path=self.results['new']['go']['gene']['level4'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='G', level=2, level_path=self.results['new']['go']['gene']['level2'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='G', level=3, level_path=self.results['new']['go']['gene']['level3'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='G', level=4, level_path=self.results['new']['go']['gene']['level4'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type='G', list_path=self.results['new']['go']['gene']['list'])
        seq_type = 'ref'
        if exp_level.lower() == 'transcript':
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type='T', level=2, level_path=self.results['ref']['go']['txpt']['level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='T', level=2, level_path=self.results['ref']['go']['txpt']['level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='T', level=3, level_path=self.results['ref']['go']['txpt']['level3'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='T', level=4, level_path=self.results['ref']['go']['txpt']['level4'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='T', level=2, level_path=self.results['ref']['go']['txpt']['level2'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='T', level=3, level_path=self.results['ref']['go']['txpt']['level3'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='T', level=4, level_path=self.results['ref']['go']['txpt']['level4'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type='T', list_path=self.results['ref']['go']['txpt']['list'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type='G', level=2, level_path=self.results['ref']['go']['gene']['level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='G', level=2, level_path=self.results['ref']['go']['gene']['level2'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='G', level=3, level_path=self.results['ref']['go']['gene']['level3'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type='G', level=4, level_path=self.results['ref']['go']['gene']['level4'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='G', level=2, level_path=self.results['ref']['go']['gene']['level2'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='G', level=3, level_path=self.results['ref']['go']['gene']['level3'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type='G', level=4, level_path=self.results['ref']['go']['gene']['level4'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type='G', list_path=self.results['ref']['go']['gene']['list'])
        if self.has_new:
            seq_type = 'all'
            if exp_level.lower() == 'transcript':
                self.add_annotation_go_all(go_id=go_id, seq_type=seq_type, anno_type='T', level=2, level_path=self.results['all']['go']['txpt']['level2'])
                self.add_annotation_go_all(go_id=go_id, seq_type=seq_type, anno_type='T', level=3, level_path=self.results['all']['go']['txpt']['level3'])
                self.add_annotation_go_all(go_id=go_id, seq_type=seq_type, anno_type='T', level=4, level_path=self.results['all']['go']['txpt']['level4'])
            self.add_annotation_go_all(go_id=go_id, seq_type=seq_type, anno_type='G', level=2, level_path=self.results['all']['go']['gene']['level2'])
            self.add_annotation_go_all(go_id=go_id, seq_type=seq_type, anno_type='G', level=3, level_path=self.results['all']['go']['gene']['level3'])
            self.add_annotation_go_all(go_id=go_id, seq_type=seq_type, anno_type='G', level=4, level_path=self.results['all']['go']['gene']['level4'])
        self.update_db_record('sg_annotation_go', go_id, status='end', main_id=go_id)

        # kegg
        kegg_detail_tables = ['sg_annotation_kegg_categories', 'sg_annotation_kegg_level', 'sg_annotation_kegg_table', 'sg_annotation_kegg_pic']
        self.remove_table_by_main_record(
            main_table='sg_annotation_kegg',
            detail_table=kegg_detail_tables,
            detail_table_key='kegg_id',
            task_id=self.task_id,
            type=self.annot_type,
        )
        self.kegg_params = json.dumps(
            dict([(k, params_dict.get(k, None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')]),
            sort_keys=True, separators=(',', ':')
        )
        kegg_id = self.add_annotation_kegg()
        self.categories = set()
        if self.has_new:
            seq_type = 'new'
            if exp_level.lower() == 'transcript':
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', layer_path=self.results['new']['kegg']['txpt']['layer'])
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', pathway_file=self.results['new']['kegg']['txpt']['pathway']['file'])
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', pathway_file=self.results['new']['kegg']['txpt']['pathway']['file'], pathway_dir=self.results['new']['kegg']['txpt']['pathway']['dir'])
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', table_path=self.results['new']['kegg']['txpt']['table'])
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', layer_path=self.results['new']['kegg']['gene']['layer'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', pathway_file=self.results['new']['kegg']['gene']['pathway']['file'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', pathway_file=self.results['new']['kegg']['gene']['pathway']['file'], pathway_dir=self.results['new']['kegg']['gene']['pathway']['dir'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', table_path=self.results['new']['kegg']['gene']['table'])
        seq_type = 'ref'
        if exp_level.lower() == 'transcript':
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', layer_path=self.results['ref']['kegg']['txpt']['layer'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', pathway_file=self.results['ref']['kegg']['txpt']['pathway']['file'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', pathway_file=self.results['ref']['kegg']['txpt']['pathway']['file'], pathway_dir=self.results['ref']['kegg']['txpt']['pathway']['dir'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type='T', table_path=self.results['ref']['kegg']['txpt']['table'])
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', layer_path=self.results['ref']['kegg']['gene']['layer'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', pathway_file=self.results['ref']['kegg']['gene']['pathway']['file'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', pathway_file=self.results['ref']['kegg']['gene']['pathway']['file'], pathway_dir=self.results['ref']['kegg']['gene']['pathway']['dir'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type='G', table_path=self.results['ref']['kegg']['gene']['table'])
        if self.has_new:
            if exp_level.lower() == 'transcript':
                self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type='all', anno_type='T', ref_layer_path=self.results['ref']['kegg']['txpt']['layer'], new_layer_path=self.results['new']['kegg']['txpt']['layer'])
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type='all', anno_type='T', pathway_file=self.results['all']['kegg']['txpt']['pathway']['file'])
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type='all', anno_type='T', pathway_file=self.results['all']['kegg']['txpt']['pathway']['file'], pathway_dir=self.results['all']['kegg']['txpt']['pathway']['dir'])
            self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type='all', anno_type='G', ref_layer_path=self.results['ref']['kegg']['gene']['layer'], new_layer_path=self.results['new']['kegg']['gene']['layer'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type='all', anno_type='G', pathway_file=self.results['all']['kegg']['gene']['pathway']['file'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type='all', anno_type='G', pathway_file=self.results['all']['kegg']['gene']['pathway']['file'], pathway_dir=self.results['all']['kegg']['gene']['pathway']['dir'])
        self.update_db_record('sg_annotation_kegg', kegg_id, status='end', categories=list(self.categories), main_id=kegg_id)

    def set_results(self, annot_merge_output_dir):
        self.ref_dir = os.path.join(annot_merge_output_dir, 'refannot_class')
        self.new_dir = os.path.join(annot_merge_output_dir, 'newannot_class')
        self.all_dir = os.path.join(annot_merge_output_dir, 'allannot_class')
        self.ref_results = {
            't2g2r2l2p': os.path.join(self.ref_dir, 'all_tran2gene.txt'),
            'stat': os.path.join(self.ref_dir, 'all_stat.xls'),
            'query': os.path.join(self.ref_dir, 'all_annot.xls'),
            'cog': {
                'txpt': os.path.join(self.ref_dir, 'cog/summary.T.tsv'),
                'gene': os.path.join(self.ref_dir, 'cog/summary.G.tsv')
            },
            'go': {
                'txpt': {
                    'level2': os.path.join(self.ref_dir, 'go/go_lev2_tran.stat.xls'),
                    'level3': os.path.join(self.ref_dir, 'go/go_lev3_tran.stat.xls'),
                    'level4': os.path.join(self.ref_dir, 'go/go_lev4_tran.stat.xls'),
                    'list': os.path.join(self.ref_dir, 'go/go_list_tran.xls')
                },
                'gene': {
                    'level2': os.path.join(self.ref_dir, 'go/go_lev2_gene.stat.xls'),
                    'level3': os.path.join(self.ref_dir, 'go/go_lev3_gene.stat.xls'),
                    'level4': os.path.join(self.ref_dir, 'go/go_lev4_gene.stat.xls'),
                    'list': os.path.join(self.ref_dir, 'go/go_list_gene.xls')
                }
            },
            'kegg': {
                'txpt': {
                    'layer': os.path.join(self.ref_dir, 'kegg/kegg_layer_tran.xls'),
                    'pathway': {
                        'file': os.path.join(self.ref_dir, 'kegg/kegg_pathway_tran.xls'),
                        'dir': os.path.join(self.ref_dir, 'kegg/kegg_pathway_tran_dir')
                    },
                    'table': os.path.join(self.ref_dir, 'kegg/kegg_gene_tran.xls')
                },
                'gene': {
                    'layer': os.path.join(self.ref_dir, 'kegg/kegg_layer_gene.xls'),
                    'pathway': {
                        'file': os.path.join(self.ref_dir, 'kegg/kegg_pathway_gene.xls'),
                        'dir': os.path.join(self.ref_dir, 'kegg/kegg_pathway_gene_dir')
                    },
                    'table': os.path.join(self.ref_dir, 'kegg/kegg_gene_gene.xls')
                }
            }
        }
        self.new_results = {
            't2g2r2l2p': os.path.join(self.new_dir, 'all_tran2gene.txt'),
            'stat': os.path.join(self.new_dir, 'all_stat.xls'),
            'query': os.path.join(self.new_dir, 'all_annot.xls'),
            'cog': {
                'txpt': os.path.join(self.new_dir, 'cog/summary.T.tsv'),
                'gene': os.path.join(self.new_dir, 'cog/summary.G.tsv')
            },
            'go': {
                'txpt': {
                    'level2': os.path.join(self.new_dir, 'go/go_lev2_tran.stat.xls'),
                    'level3': os.path.join(self.new_dir, 'go/go_lev3_tran.stat.xls'),
                    'level4': os.path.join(self.new_dir, 'go/go_lev4_tran.stat.xls'),
                    'list': os.path.join(self.new_dir, 'go/go_list_tran.xls')
                },
                'gene': {
                    'level2': os.path.join(self.new_dir, 'go/go_lev2_gene.stat.xls'),
                    'level3': os.path.join(self.new_dir, 'go/go_lev3_gene.stat.xls'),
                    'level4': os.path.join(self.new_dir, 'go/go_lev4_gene.stat.xls'),
                    'list': os.path.join(self.new_dir, 'go/go_list_gene.xls')
                }
            },
            'kegg': {
                'txpt': {
                    'layer': os.path.join(self.new_dir, 'kegg/kegg_layer_tran.xls'),
                    'pathway': {
                        'file': os.path.join(self.new_dir, 'kegg/kegg_pathway_tran.xls'),
                        'dir': os.path.join(self.new_dir, 'kegg/kegg_pathway_tran_dir')
                    },
                    'table': os.path.join(self.new_dir, 'kegg/kegg_gene_tran.xls')
                },
                'gene': {
                    'layer': os.path.join(self.new_dir, 'kegg/kegg_layer_gene.xls'),
                    'pathway': {
                        'file': os.path.join(self.new_dir, 'kegg/kegg_pathway_gene.xls'),
                        'dir': os.path.join(self.new_dir, 'kegg/kegg_pathway_gene_dir')
                    },
                    'table': os.path.join(self.new_dir, 'kegg/kegg_gene_gene.xls')
                }
            }
        }
        self.all_results = {
            'stat': os.path.join(self.all_dir, 'all_stat.xls'),
            'go': {
                'txpt': {
                    'level2': os.path.join(self.all_dir, 'go/go_lev2_tran.stat.xls'),
                    'level3': os.path.join(self.all_dir, 'go/go_lev3_tran.stat.xls'),
                    'level4': os.path.join(self.all_dir, 'go/go_lev4_tran.stat.xls')
                },
                'gene': {
                    'level2': os.path.join(self.all_dir, 'go/go_lev2_gene.stat.xls'),
                    'level3': os.path.join(self.all_dir, 'go/go_lev3_gene.stat.xls'),
                    'level4': os.path.join(self.all_dir, 'go/go_lev4_gene.stat.xls')
                }
            },
            'kegg': {
                'txpt': {
                    'pathway': {
                        'file': os.path.join(self.all_dir, 'kegg/kegg_pathway_tran.xls'),
                        'dir': os.path.join(self.all_dir, 'kegg/kegg_pathway_tran_dir')
                    }
                },
                'gene': {
                    'pathway': {
                        'file': os.path.join(self.all_dir, 'kegg/kegg_pathway_gene.xls'),
                        'dir': os.path.join(self.all_dir, 'kegg/kegg_pathway_gene_dir')
                    }
                }
            }
        }
        self.results = {'ref': self.ref_results, 'all': self.all_results}
        if self.has_new:
            self.results.update({'new': self.new_results})
        def check(dct):
            for k, v in dct.items():
                if isinstance(v, dict):
                    check(v)
                else:
                    if os.path.exists(v):
                        self.bind_object.logger.info('succeed in finding {}'.format(v))
                    else:
                        self.bind_object.set_error('fail to find %s', variables=(v), code="53703756")
        check(self.results)

    def set_relation(self, ref_tran2gene, new_tran2gene=None):
        self.bind_object.logger.info('start reading {}'.format(ref_tran2gene))
        for line in open(ref_tran2gene):
            items = line.strip().split('\t')
            self.t2g_dict[items[0]] = items[1]
            if items[2] == 'yes':
                self.t2r_dict[items[0]] = True
            else:
                self.t2r_dict[items[0]] = False
        if self.has_new:
            self.bind_object.logger.info('start reading {}'.format(new_tran2gene))
            for line in open(new_tran2gene):
                items = line.strip().split('\t')
                self.t2g_dict[items[0]] = items[1]
                if items[2] == 'yes' and items[2] not in self.t2r_dict:
                    self.t2r_dict[items[0]] = True
                else:
                    self.t2r_dict[items[0]] = False
            self.bind_object.logger.info('succeed in building relationship dict')

    @report_check
    def add_annotation_stat(self):
        insert_data = {
            'task_id': self.task_id,
            'project_sn': self.project_sn,
            'created_ts': self.time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'name': 'Annotation_stat_{}'.format(self.time_now.strftime('%Y%m%d_%H%M%S')),
            'params': self.params,
            'taxonomy': self.taxonomy,
            'exp_level': self.exp_level,
            'type': self.annot_type,
            'has_new': self.has_new,
            'status': 'start',
            'desc': 'annotation statistics main table'
        }
        collection = self.db['sg_annotation_stat']
        stat_id = collection.insert_one(insert_data).inserted_id
        if stat_id:
            self.bind_object.logger.info('succeed in inserting record into sg_annotation_stat')
        return stat_id

    @report_check
    def add_annotation_stat_detail(self, stat_id, stat_path, root_dir, seq_type, exp_level):
        def list_str(venn_file):
            return ','.join(set(line.strip() for line in open(os.path.join(root_dir, venn_file))))
        map_dict = {
            'transcript': {
                'nr': 'nr/nr_venn_tran.txt',
                'swissprot': 'swissprot/swissprot_venn_tran.txt',
                'cog': 'cog/cog_venn_tran.txt',
                'kegg': 'kegg/kegg_venn_tran.txt',
                'pfam': 'pfam/pfam_venn_tran.txt',
                'go': 'go/go_venn_tran.txt',
            },
            'gene': {
                'nr': 'nr/nr_venn_gene.txt',
                'swissprot': 'swissprot/swissprot_venn_gene.txt',
                'cog': 'cog/cog_venn_gene.txt',
                'kegg': 'kegg/kegg_venn_gene.txt',
                'pfam': 'pfam/pfam_venn_gene.txt',
                'go': 'go/go_venn_gene.txt',
            },
        }
        insert_data = list()
        records = list(pd.read_table(stat_path).to_dict('r'))
        for n, items in enumerate(records):
            data = {
                'stat_id': stat_id,
                'type': items['type'],
                'seq_type': seq_type,
                'gene': int(items['gene']),
                'gene_percent': round(float(items['gene_percent']), 4)
            }
            if items['type'] in map_dict['gene']:
                data.update({'gene_list': list_str(map_dict['gene'][items['type']])})
            if exp_level == 'transcript':
                data.update({
                    'transcript': int(items['transcript']),
                    'transcript_percent': round(float(items['transcript_percent']), 4)
                })
                if items['type'] in map_dict['transcript']:
                    data.update({'transcript_list': list_str(map_dict['transcript'][items['type']])})
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_stat_detail', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_stat_detail')

    @report_check
    def add_annotation_query(self, stat_id):
        name = 'Annotation_query_{}'.format(self.time_now.strftime('%Y%m%d_%H%M%S'))
        insert_data = {
            'task_id': self.task_id,
            'project_sn': self.project_sn,
            'created_ts': self.time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'name': 'Annotation_query_{}'.format(self.time_now.strftime('%Y%m%d_%H%M%S')),
            'params': self.params,
            'type': self.annot_type,
            'stat_id': stat_id,
            'status': 'start',
            'desc': 'annotation query main table'
        }
        collection = self.db['sg_annotation_query']
        query_id = collection.insert_one(insert_data).inserted_id
        if query_id:
            self.bind_object.logger.info('succeed in inserting record into sg_annotation_query')
        return query_id

    @report_check
    def add_annotation_query_detail(self, query_id, query_path, seq_type, exp_level):
        insert_data = list()
        lines = open(query_path).readlines()
        for line in lines[1:]:
            items = line.strip().split('\t')
            is_gene = True
            if self.t2r_dict.has_key(items[1]) and self.t2r_dict[items[1]] == False:
                is_gene = False
            if exp_level.lower() == "gene" and is_gene == False:
                continue
            data = [
                ('query_id', query_id),
                ('seq_type', seq_type),
                ('transcript_id', items[1]),
                ('gene_id', items[0]),
                ('is_gene', is_gene),
                ('gene_name', items[3]),
            ]
            try:
                data.append(('length', items[4]))
            except:
                data.append(('length', None))
            try:
                data.append(('description', items[5]))
            except:
                data.append(('description', None))
            try:
                data.append(('cog', items[6]))
                data.append(('cog_description', items[7]))
            except:
                data.append(('cog', None))
                data.append(('cog_description', None))
            try:
                data.append(('ko_id', items[8]))
            except:
                data.append(('ko_id', None))
            try:
                data.append(('ko_name', items[9]))
            except:
                data.append(('ko_name', None))
            try:
                data.append(('pathways', items[10]))
            except:
                data.append(('pathways', None))
            try:
                data.append(('pfam', items[11]))
            except:
                data.append(('pfam', None))
            try:
                data.append(('go', items[12]))
            except:
                data.append(('go', None))
            try:
                data.append(('nr', items[13]))
            except:
                data.append(('nr', None))
            try:
                data.append(('swissprot', items[14]))
            except:
                data.append(('swissprot', None))
            try:
                data.append(('enterz', items[15]))
            except:
                data.append(('enterz', None))
            data = SON(data)
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_query_detail', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_query_detail')

    @report_check
    def add_annotation_cog(self):
        insert_data = {
            'task_id': self.task_id,
            'project_sn': self.project_sn,
            'created_ts': self.time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'name': 'Annotation_cog_{}_{}'.format(self.annot_type, self.time_now.strftime('%Y%m%d_%H%M%S')),
            'params': self.cog_params,
            'type': self.annot_type,
            'status': 'start',
            'desc': 'annotation cog main table'
        }
        collection = self.db['sg_annotation_cog']
        cog_id = collection.insert_one(insert_data).inserted_id
        if cog_id:
            self.bind_object.logger.info('succeed in inserting record into sg_annotation_cog')
        return cog_id

    @report_check
    def add_annotation_cog_detail(self, cog_id, cog_path, seq_type, anno_type):
        insert_data = list()
        lines = open(cog_path).readlines()
        for line in lines[1:]:
            items = line.strip().split('\t')
            data = [
                ('cog_id', cog_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('type', items[0]),
                ('function_categories', '[{}] {}'.format(items[2], items[1])),
                ('cog', int(items[3])),
             ]
            try:
                data.append(('cog_list', items[4]))
            except:
                data.append(('cog_list', None))
            insert_data.append(SON(data))
        else:
            self.create_db_table('sg_annotation_cog_detail', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_cog_detail')

    @report_check
    def add_annotation_cog_detail_all(self, cog_id, ref_cog_path, new_cog_path, seq_type, anno_type):
        func_type = {
            'INFORMATION STORAGE AND PROCESSING': ['A', 'B', 'J', 'K', 'L'],
            'CELLULAR PROCESSES AND SIGNALING': ['D', 'M', 'N', 'O', 'T', 'U', 'V', 'W', 'Y', 'Z'],
            'METABOLISM': ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q'],
            'POORLY CHARACTERIZED': ['R', 'S'],
        }
        func_decs = {
            'A': 'RNA processing and modification',
            'B': 'Chromatin structure and dynamics',
            'C': 'Energy production and conversion',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism',
            'G': 'Carbohydrate transport and metabolism',
            'H': 'Coenzyme transport and metabolism',
            'I': 'Lipid transport and metabolism',
            'J': 'Translation, ribosomal structure and biogenesis',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only',
            'S': 'Function unknown',
            'T': 'Signal transduction mechanisms',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'V': 'Defense mechanisms',
            'W': 'Extracellular structures',
            'Y': 'Nuclear structure',
            'Z': 'Cytoskeleton'
        }
        insert_data = list()
        func_dict = {'COG': dict(), 'NOG': dict()}
        ref_lines = open(ref_cog_path).readlines()
        new_lines = open(new_cog_path).readlines()
        for ref_line in ref_lines[1:]:
            items = ref_line.strip().split('\t')
            items[1] = '[{}] {}'.format(items[2], items[1])
            m = re.match(r'\[(.+)\].+$', items[1])
            if m:
                func = m.group(1)
                try:
                    func_dict['COG'][func] = items[4].split(';')
                except:
                    func_dict['COG'][func] = list()
                try:
                    func_dict['NOG'][func] = items[5].split(';')
                except:
                    func_dict['NOG'][func] = list()
        for new_line in new_lines[1:]:
            items = new_line.strip().split('\t')
            items[1] = '[{}] {}'.format(items[2], items[1])
            m = re.match(r'\[(.+)\].+$', items[1])
            if m:
                func = m.group(1)
                if func in func_dict['COG']:
                    try:
                        tmp_set = set(func_dict['COG'][func])
                        tmp_set.update(set(items[4].split(';')))
                        func_dict['COG'][func] = list(tmp_set)
                    except:
                        pass
                else:
                    try:
                        func_dict['COG'][func] = items[4].split(';')
                    except:
                        func_dict['COG'][func] = list()
                if func in func_dict['NOG']:
                    try:
                        tmp_set = set(func_dict['NOG'][func])
                        tmp_set.update(set(items[5].split(';')))
                        func_dict['NOG'][func] = list(tmp_set)
                    except:
                        pass
                else:
                    try:
                        func_dict['NOG'][func] = items[5].split(';')
                    except:
                        func_dict['NOG'][func] = list()
        for func in func_dict['COG'].keys():
            function_categories = '[{}] {}'.format(func, func_decs[func])
            try:
                cog_list = list(set(func_dict['COG'][func]))
            except:
                cog_list = list()
            try:
                nog_list = list(set(func_dict['NOG'][func]))
            except:
                nog_list = list()
            for cog_class in func_type.keys():
                if func in func_type[cog_class]:
                    cog_type = cog_class
            data = [
                ('cog_id', cog_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('type', cog_type),
                ('function_categories', function_categories),
                ('cog', len([x for x in cog_list if x])),
                ('cog_list', ';'.join(cog_list)),
            ]
            data = SON(data)
            if len(cog_list) + len(nog_list) != 0:
                insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_cog_detail', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_cog_detail')

    @report_check
    def add_annotation_go(self):
        insert_data = {
            'task_id': self.task_id,
            'project_sn': self.project_sn,
            'created_ts': self.time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'name': 'Annotation_go_{}_{}'.format(self.annot_type, self.time_now.strftime('%Y%m%d_%H%M%S')),
            'type': self.annot_type,
            'params': self.go_params,
            'status': 'start',
            'desc': 'annotation go main table',
        }
        collection = self.db['sg_annotation_go']
        go_id = collection.insert_one(insert_data).inserted_id
        if go_id:
            self.bind_object.logger.info('succeed in inserting record into sg_annotation_go')
        return go_id

    @report_check
    def add_annotation_go_level(self, go_id, seq_type, anno_type, level, level_path):
        insert_data = list()
        lines = open(level_path).readlines()
        for line in lines[1:]:
            items = line.strip().split('\t')
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('parent_name', items[0]),
                ('term_type', items[1]),
                ('go', items[2]),
                ('num', int(items[3])),
                ('percent', round(float(items[4]), 4))
            ]
            data = SON(data)
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_go_level', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_go_level')

    @report_check
    def add_annotation_go_detail(self, go_id, seq_type, anno_type, level, level_path):
        insert_data = list()
        lines = open(level_path).readlines()
        for line in lines[1:]:
            items = line.strip().split('\t')
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('goterm', items[0]),
                ('goterm_2', items[1]),
                ('goid_2', items[2]),
                ('seq_number', int(items[-3])),
                ('percent', round(float(items[-2]), 4)),
            ]
            if level == 2:
                data.append(('seq_list', items[-1]))
            if level >= 3:
                data.append(('goterm_3', items[3]))
                data.append(('goid_3', items[4]))
            if level == 4:
                data.append(('goterm_4', items[5]))
                data.append(('goid_4', items[6]))
            data = SON(data)
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_go_detail', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_go_detail')

    @report_check
    def add_annotation_go_graph(self, go_id, seq_type, anno_type, level, level_path):
        insert_data = list()
        term_dict = dict()
        lines = open(level_path).readlines()
        for i in range(1, len(lines)):
            items = lines[i].strip().split('\t')
            if level == 2:
                term_type = items[0]
                go_term = items[1]
                if go_term not in term_dict:
                    term_dict[go_term] = list()
                    term_dict[go_term].append(i)
                else:
                    term_dict[go_term].append(i)
            if level == 3:
                term_type = items[0]
                go_term = items[3]
                if go_term not in term_dict:
                    term_dict[go_term] = list()
                    term_dict[go_term].append(i)
                else:
                    term_dict[go_term].append(i)
            if level == 4:
                term_type = items[0]
                go_term = items[5]
                if go_term not in term_dict:
                    term_dict[go_term] = list()
                    term_dict[go_term].append(i)
                else:
                    term_dict[go_term].append(i)
        for go_term in term_dict:
            for j in term_dict[go_term]:
                items = lines[j].strip().split('\t')
                term_type = items[0]
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('term_type', term_type),
                ('go_term', go_term),
                ('seq_number', items[-3]),
                ('percent', items[-2])
            ]
            data = SON(data)
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_go_graph', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_go_graph')

    @report_check
    def add_annotation_go_list(self, go_id, seq_type, anno_type, list_path):
        insert_data = list()
        for line in open(list_path):
            items = line.strip().split('\t')
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('gene_id', items[0]),
                ('gos_list', items[1]),
            ]
            data = SON(data)
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_go_list', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_go_list')

    @report_check
    def add_annotation_go_all(self, go_id, seq_type, anno_type, level, level_path):
        graph_data, detail_data, query_ids = list(), list(), list()
        func_dict, term_dict = dict(), dict()
        lines = open(level_path).readlines()
        for line in lines[1:]:
            items = line.strip().split('\t')
            func = '{}|||{}|||{}'.format(items[0], items[1], items[2])
            term = '{}|||{}'.format(items[0], items[1])
            if level == 3:
                func += '|||{}|||{}'.format(items[3], items[4])
                term = '{}|||{}'.format(items[0], items[3])
            if level == 4:
                func += '|||{}|||{}|||{}|||{}'.format(items[3], items[4], items[5], items[6])
                term = '{}|||{}'.format(items[0], items[5])
            func_dict[func] = items[-1].split(';')
            if term not in term_dict:
                term_dict[term] = set(items[-1].split(';'))
            else:
                for seq_id in items[-1].split(';'):
                    if seq_id not in term_dict[term]:
                        term_dict[term].add(seq_id)
            query_ids.extend(items[-1].split(';'))
        else:
            query_ids = list(set(query_ids))
        for term in sorted(term_dict):
            terms = term.split('|||')
            seq_list = term_dict[term]
            percent = float(len(seq_list)) / len(query_ids)
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('term_type', terms[0]),
                ('go_term', terms[1]),
                ('seq_number', len(seq_list)),
                ('percent', round(percent, 4))
            ]
            data = SON(data)
            graph_data.append(data)
        else:
            self.create_db_table('sg_annotation_go_graph', graph_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_go_graph')
        for func in func_dict:
            funcs = func.split('|||')
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
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
            self.create_db_table('sg_annotation_go_detail', graph_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_go_detail')

    @report_check
    def add_annotation_kegg(self):
        insert_data = {
            'task_id': self.task_id,
            'project_sn': self.project_sn,
            'created_ts': self.time_now.strftime('%Y-%m-%d %H:%M:%S'),
            'name': 'Annotation_kegg_{}_{}'.format(self.annot_type, self.time_now.strftime('%Y%m%d_%H%M%S')),
            'type': self.annot_type,
            'params': self.kegg_params,
            'status': 'start',
            'desc': 'annotation kegg main table',
            'categories': ['M','CP','EIP','GIP','OS','HD','DD']
        }
        collection = self.db['sg_annotation_kegg']
        kegg_id = collection.insert_one(insert_data).inserted_id
        if kegg_id:
            self.bind_object.logger.info('succeed in inserting record into sg_annotation_kegg')
        return kegg_id

    @report_check
    def add_annotation_kegg_categories(self, kegg_id, seq_type, anno_type, layer_path):
        insert_data = list()
        for line in open(layer_path):
            items = line.strip().split('\t')
            self.categories.add(''.join([x[0] for x in items[0].split(' ')]))
            data = [
                ('kegg_id', kegg_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('first_category', items[0]),
                ('second_category', items[1]),
                ('num', int(items[2])),
                ('seq_list', items[3]),
            ]
            data = SON(data)
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_kegg_categories', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_kegg_categories')

    @report_check
    def add_annotation_kegg_level(self, kegg_id, seq_type, anno_type, pathway_file):
        documents = list()
        lines = open(pathway_file).readlines()
        for line in lines[1:]:
            items = line.strip('\n').split('\t')
            insert_data = {
                'kegg_id': kegg_id,
                'seq_type': seq_type,
                'anno_type': anno_type,
                'pathway_id': items[0],
                'first_category': items[1],
                'second_category': items[2],
                'pathway_definition': items[3],
                'number_of_seqs': int(items[4]),
                'seq_list': items[5],
                'hyperlink': items[-1]
            }
            documents.append(insert_data)
        else:
            self.create_db_table('sg_annotation_kegg_level', documents)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_kegg_level')

    @report_check
    def add_annotation_kegg_pic(self, kegg_id, seq_type, anno_type, pathway_file, pathway_dir):
        insert_data = list()
        lines = open(pathway_file).readlines()
        for line in lines[1:]:
            items = line.strip('\n').split('\t')
            mark_file = os.path.join(pathway_dir, '{}.html.mark'.format(items[0]))
            if os.path.exists(mark_file):
                with open(mark_file, 'r') as mark_handle:
                    for mark_line in mark_handle.readlines():
                        mark_items = mark_line.strip('\n').split('\t')
                        if len(mark_items) == 8:
                            png, shape, bg_color, fg_color, coords, title, kos, href = mark_items
                            title = title.strip()
                        else:
                            continue
                        data = {
                            'kegg_id': kegg_id,
                            'seq_type': seq_type,
                            'anno_type': anno_type,
                            'pathway_id': items[0],
                            'shape': shape,
                            'bg_colors': bg_color,
                            'fg_colors': fg_color,
                            'coords': coords,
                            'href': href,
                            'kos': kos,
                            'title': title
                        }
                        if bg_color != str() and len(bg_color.split(',')) > 0:
                            data.update({'bg_type': len(bg_color.split(','))})
                        if fg_color != str() and len(fg_color.split(',')) > 0:
                            data.update({'fg_type': len(fg_color.split(','))})
                        insert_data.append(data)
            else:
                self.bind_object.logger.debug('fail to find {}'.format(mark_file))
        else:
            self.create_db_table('sg_annotation_kegg_pic', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_kegg_pic')

    @report_check
    def add_annotation_kegg_table(self, kegg_id, seq_type, anno_type, table_path):
        insert_data = list()
        lines = open(table_path).readlines()
        for line in lines[1:]:
            items = line.strip('\n').split('\t')
            data = {
                'kegg_id': kegg_id,
                'seq_type': seq_type,
                'anno_type': anno_type,
                'transcript_id': items[0],
                'ko_id': items[1],
                'ko_name': items[2],
                'hyperlink': items[3],
                'paths': items[4],
            }
            insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_kegg_table', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_kegg_table')

    @report_check
    def add_annotation_kegg_categories_all(self, kegg_id, seq_type, anno_type, ref_layer_path, new_layer_path):
        insert_data = list()
        c2s_dict = dict()
        ref_lines = open(ref_layer_path).readlines()
        new_lines = open(new_layer_path).readlines()
        for line in ref_lines:
            items = line.strip().split('\t')
            if items[0] == 'Metabolism' and items[1] == 'Global and overview maps':
                pass
            else:
                c_key = '{}|||{}'.format(items[0], items[1])
                c2s_dict[c_key] = items[3].split(';')
        for line in new_lines:
            items = line.strip().split('\t')
            if items[0] == 'Metabolism' and items[1] == 'Global and overview maps':
                pass
            else:
                c_key = '{}|||{}'.format(items[0], items[1])
                if c_key not in c2s_dict:
                    c2s_dict[c_key] = items[3].split(';')
                else:
                    tmp_set = set(c2s_dict[c_key])
                    tmp_set.update(items[3].split(';'))
                    c2s_dict[c_key] = list(tmp_set)
        for fc_bak in ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing', 'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']:
            for c_key in c2s_dict:
                categories = c_key.split("|||")
                first_category = categories[0]
                second_category = categories[1]
                if fc_bak == first_category:
                    num = len(c2s_dict[c_key])
                    seq_list = ';'.join(c2s_dict[c_key])
                    data = [
                        ('kegg_id', kegg_id),
                        ('seq_type', seq_type),
                        ('anno_type', anno_type),
                        ('first_category', first_category),
                        ('second_category', second_category),
                        ('num', num),
                        ('seq_list', seq_list)
                    ]
                    data = SON(data)
                    insert_data.append(data)
        else:
            self.create_db_table('sg_annotation_kegg_categories', insert_data)
            self.bind_object.logger.info('succeed in inserting records in sg_annotation_kegg_categories')

    # the following functions were deprecated
    def set_result_dir(self, merge_annot_path):
        merge_annot_path =  merge_annot_path + "/"
        self.result_dir = merge_annot_path
        self.result_file['new_stat_path'] = os.path.join(merge_annot_path, "newannot_class/all_stat.xls")
        self.result_file['new_venn_path'] = os.path.join(merge_annot_path, "newannot_class")
        self.result_file['ref_stat_path'] = os.path.join(merge_annot_path, "refannot_class/all_stat.xls")
        self.result_file['ref_venn_path'] = os.path.join(merge_annot_path, "refannot_class/")
        self.result_file['all_stat_path'] = os.path.join(merge_annot_path, "allannot_class/all_stat.xls")
        self.result_file['all_venn_path'] = os.path.join(merge_annot_path, "allannot_class/")

        for db in ["nr", "swissprot"]:
            db_gene_name = db + "genes_blast_path"
            db_trans_name = db + "trans_blast_path"
            self.result_file[db_gene_name] = merge_annot_path + "newannot_class/" + db + "/{}_blast_gene.xls".format(db)
            self.result_file[db_trans_name] = merge_annot_path + "newannot_class/" + db + "/{}_blast_tran.xls".format(db)
            self.result_file[db_gene_name + "ref"] = merge_annot_path + "refannot_class/" + db + "/{}_blast_gene.xls".format(db)
            self.result_file[db_trans_name + "ref"] = merge_annot_path + "refannot_class/" + db + "/{}_blast_tran.xls".format(db)

        self.result_file['gene_pfam_path'] = merge_annot_path + "newannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path'] = merge_annot_path + "newannot_class/pfam/pfam_domain_tran.xls"
        self.result_file['query_path'] = merge_annot_path + "newannot_class/all_annot.xls"

        self.result_file['gene_pfam_path' + "ref"] = merge_annot_path + "refannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path' + "ref"] = merge_annot_path + "refannot_class/pfam/pfam_domain_tran.xls"
        self.result_file['query_path' + "ref"] = merge_annot_path + "refannot_class/all_annot.xls"

        self.result_file['gene_pfam_path' + "all"] = merge_annot_path + "allannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path' + "all"] = merge_annot_path + "allannot_class/pfam/pfam_domain_tran.xls"

        self.result_file['n_sum_path'] = merge_annot_path + "newannot_class/cog/cog_summary_tran.xls"
        self.result_file['n_gene_sum_path'] = merge_annot_path + "newannot_class/cog/cog_summary_gene.xls"

        self.result_file['n_sum_path' + 'ref'] = merge_annot_path + "refannot_class/cog/cog_summary_tran.xls"
        self.result_file['n_gene_sum_path' + 'ref'] = merge_annot_path + "refannot_class/cog/cog_summary_gene.xls"

        self.result_file['stat_level2'] = merge_annot_path + "newannot_class/go/go_lev2_tran.stat.xls"
        self.result_file['stat_level3'] = merge_annot_path + "newannot_class/go/go_lev3_tran.stat.xls"
        self.result_file['stat_level4'] = merge_annot_path + "newannot_class/go/go_lev4_tran.stat.xls"
        self.result_file['gene_stat_level2'] = merge_annot_path + "newannot_class/go/go_lev2_gene.stat.xls"
        self.result_file['gene_stat_level3'] = merge_annot_path + "newannot_class/go/go_lev3_gene.stat.xls"
        self.result_file['gene_stat_level4'] = merge_annot_path + "newannot_class/go/go_lev4_gene.stat.xls"
        self.result_file['gos_path'] = merge_annot_path + "newannot_class/go/go_list_tran.xls"
        self.result_file['gene_gos_path'] = merge_annot_path + "newannot_class/go/go_list_gene.xls"

        self.result_file['stat_level2' + '_ref'] = merge_annot_path + "refannot_class/go/go_lev2_tran.stat.xls"
        self.result_file['stat_level3' + '_ref'] = merge_annot_path + "refannot_class/go/go_lev3_tran.stat.xls"
        self.result_file['stat_level4' + '_ref'] = merge_annot_path + "refannot_class/go/go_lev4_tran.stat.xls"
        self.result_file['gene_stat_level2' + '_ref'] = merge_annot_path + "refannot_class/go/go_lev2_gene.stat.xls"
        self.result_file['gene_stat_level3' + '_ref'] = merge_annot_path + "refannot_class/go/go_lev3_gene.stat.xls"
        self.result_file['gene_stat_level4' + '_ref'] = merge_annot_path + "refannot_class/go/go_lev4_gene.stat.xls"
        self.result_file['gos_path' + '_ref'] = merge_annot_path + "refannot_class/go/go_list_tran.xls"
        self.result_file['gene_gos_path' + '_ref'] = merge_annot_path + "refannot_class/go/go_list_gene.xls"

        self.result_file['stat_level2' + '_all'] = merge_annot_path + "allannot_class/go/go_lev2_tran.stat.xls"
        self.result_file['stat_level3' + '_all'] = merge_annot_path + "allannot_class/go/go_lev3_tran.stat.xls"
        self.result_file['stat_level4' + '_all'] = merge_annot_path + "allannot_class/go/go_lev4_tran.stat.xls"
        self.result_file['gene_stat_level2' + '_all'] = merge_annot_path + "allannot_class/go/go_lev2_gene.stat.xls"
        self.result_file['gene_stat_level3' + '_all'] = merge_annot_path + "allannot_class/go/go_lev3_gene.stat.xls"
        self.result_file['gene_stat_level4' + '_all'] = merge_annot_path + "allannot_class/go/go_lev4_gene.stat.xls"
        self.result_file['gos_path' + '_all'] = merge_annot_path + "allannot_class/go/go_list_tran.xls"
        self.result_file['gene_gos_path' + '_all'] = merge_annot_path + "allannot_class/go/go_list_gene.xls"

        self.result_file['layer_path'] = merge_annot_path + "newannot_class/kegg/kegg_layer_tran.xls"
        self.result_file['pathway_path'] = merge_annot_path + "newannot_class/kegg/kegg_pathway_tran.xls"
        self.result_file['png_path'] = merge_annot_path + "newannot_class/kegg/kegg_pathway_tran_dir/"
        self.result_file['table_path'] = merge_annot_path + "newannot_class/kegg/kegg_gene_tran.xls"
        self.result_file['gene_layer_path'] = merge_annot_path + "newannot_class/kegg/kegg_layer_gene.xls"
        self.result_file['gene_pathway_path'] = merge_annot_path + "newannot_class/kegg/kegg_pathway_gene.xls"
        self.result_file['gene_png_path'] = merge_annot_path + "newannot_class/kegg/kegg_pathway_gene_dir/"
        self.result_file['gene_table_path'] = merge_annot_path + "newannot_class/kegg/kegg_gene_gene.xls"

        self.result_file['layer_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_layer_tran.xls"
        self.result_file['pathway_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_pathway_tran.xls"
        self.result_file['png_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_pathway_tran_dir/"
        self.result_file['table_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_gene_tran.xls"
        self.result_file['gene_layer_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_layer_gene.xls"
        self.result_file['gene_pathway_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_pathway_gene.xls"
        self.result_file['gene_png_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_pathway_gene_dir/"
        self.result_file['gene_table_path' + '_ref'] = merge_annot_path + "refannot_class/kegg/kegg_gene_gene.xls"

        self.result_file['pathway_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_tran.xls"
        self.result_file['png_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_tran_dir/"
        self.result_file['table_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_gene_tran.xls"

        self.result_file['gene_pathway_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_gene.xls"
        self.result_file['gene_png_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_gene_dir/"
        self.result_file['gene_table_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_gene_gene.xls"

        for key, value in self.result_file.items():
            if self.has_new == False:
                if key.endswith("_ref") or key.startswith("ref_"):
                    if os.path.exists(value):
                        pass
                    else:
                        self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="53703757")
            else:
                if os.path.exists(value):
                    pass
                else:
                    self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="53703758")
        self.bind_object.logger.info("数据路径正确，文件完整 {}")

    def get_pic(self, path, kos_path, png_path):
        """
        画通路图
        """
        fs = gridfs.GridFS(self.mongodb)
        pid = re.sub("map", "ko", path)
        with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
            result = self.png_coll.find_one({"pathway_id": pid})
            if result:
                kgml_id = result['pathway_ko_kgml']
                png_id = result['pathway_map_png']
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
        cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path, kos_path, png_path, "pathway.kgml", "pathway.png")
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "{}画图出错".format(path)
            os.system("cp {} {}".format("pathway.png", png_path))

    @report_check
    def add_annotation_kegg_level_all(self, kegg_id, seq_type, anno_type, r_level_path, n_level_path):
        """
        r_level_path: pathway_table.xls(ref)
        n_level_path: pathway_table.xls(new)
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="53703759")
        if not os.path.exists(r_level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_level_path), code="53703760")
        if not os.path.exists(n_level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_level_path), code="53703761")
        data_list = []
        path_def = {}
        r_path_list = {}
        n_path_list = {}
        fs = gridfs.GridFS(self.mongodb)
        with open(r_level_path, "rb") as r, open(n_level_path, "rb") as n:
            lines = r.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[1] == "Metabolism" and line[2] == "Global and overview maps":
                    pass
                else:
                    path = line[0] + "|||" + line[1] + "|||" + line[2] + "|||" + line[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + line[0]
                    k_cols = line[7].split(sp)
                    k_ids = k_cols[1].split("%09yellow")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
                            k_list.append(k_id)
                    seqlist = line[5].split(";")
                    path_def[line[0]] = path
                    r_path_list[line[0]] = []
                    r_path_list[line[0]].append(seqlist)
                    r_path_list[line[0]].append(k_list)
            for item in items[1:]:
                item = item.strip().split("\t")
                if item[1] == "Metabolism" and item[2] == "Global and overview maps":
                    pass
                else:
                    path = item[0] + "|||" + item[1] + "|||" + item[2] + "|||" + item[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + item[0]
                    k_cols = item[7].split(sp)
                    k_ids = k_cols[1].split("%09green")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
                            k_list.append(k_id)
                    seqlist = item[5].split(";")
                    n_path_list[item[0]] = []
                    n_path_list[item[0]].append(seqlist)
                    n_path_list[item[0]].append(k_list)
                    if item[0] not in path_def:
                        path_def[item[0]] = path
        for map_id in path_def:
            link = []
            r_kos, n_kos, b_kos = [], [], []
            ref, new, both = [], [], []
            paths = path_def[map_id].split("|||")
            first_category = paths[1]
            second_category = paths[2]
            pathway_definition = paths[3]
            try:
                seq_list = r_path_list[map_id][0]
                r_ko = r_path_list[map_id][1]
            except:
                seq_list = []
                r_ko = []
            try:
                for q in n_path_list[map_id][0]:
                    seq_list.append(q)
                n_ko = n_path_list[map_id][1]
            except:
                n_ko = []
            seq_list = list(set(seq_list))
            for r_c in r_ko:
                r = r_c.split("%09tomato")[0]
                if r_c not in n_ko:
                    r_id = r + '%09' + 'yellow'
                    link.append(r_id)
                    r_kos.append(r)
                else:
                    b_id = r + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(r)
            for n_c in n_ko:
                n = n_c.split("%09tomato")[0]
                if n_c not in r_ko:
                    n_id = n + '%09' + 'green'
                    link.append(n_id)
                    n_kos.append(n)
                else:
                    b_id = n + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(n)
            link = list(set(link))
            b_kos = list(set(b_kos))
            link = 'http://www.genome.jp/kegg-bin/show_pathway?' + map_id + '/' + '/'.join(link)
            png_path = os.getcwd() + '/' + map_id + ".png"
            pdf_path = os.getcwd() + '/' + map_id + ".pdf"
            kos_path = os.path.join(os.getcwd(), "KOs.txt")
            with open(kos_path, "w") as w:
                w.write("#KO\tbg\tfg\n")
                for k in n_kos:
                    w.write(k + "\t" + "#00CD00" + "\t" + "NA" + "\n")
                for k in r_kos:
                    w.write(k + "\t" + "#FFFF00" + "\t" + "NA" + "\n")
                for k in b_kos:
                    w.write(k + "\t" + "#FFFF00,#00CD00" + "\t" + "NA" + "\n")
            self.get_pic(map_id, kos_path, png_path)
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print '图片格式pdf转png出错'
            # pdfid = fs.put(open(pdf_path, 'rb'))
            # graph_png_id = fs.put(open(png_path, 'rb'))
            insert_data = {
                'kegg_id': kegg_id,
                'seq_type': seq_type,
                'anno_type': anno_type,
                'pathway_id': map_id,
                'first_category': first_category,
                'second_category': second_category,
                'pathway_definition': pathway_definition,
                'number_of_seqs': len(seq_list),
                'seq_list': ";".join(seq_list),
                # 'graph_id': pdfid,
                "hyperlink": link,
                # 'graph_png_id': graph_png_id
            }
            data_list.append(insert_data)
            # os.remove(pdf)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_level']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释层级all信息出错：%s、%s" , variables=(r_level_path, n_level_path), code="53703762")
            else:
                self.bind_object.logger.info("导入kegg注释层级all信息成功：%s、%s" % (level_path, png_dir))

    @report_check
    def add_annotation_cog_table(self, cog_id, table_path, seq_type, anno_type):
        '''
        table_path:cog_table.xls
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                self.bind_object.set_error('cog_id必须为ObjectId对象或其对应的字符串！', code="53703763")
        if not os.path.exists(table_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(table_path), code="53703764")
        data_list = list()
        with open(table_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    # ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('query_name', line[0]),
                    ('query_length', line[1]),
                    ('hsp_start', line[2]),
                    ('hsp_end', line[3]),
                    ('hsp_strand', line[4]),
                    ('hit_name', line[5]),
                    ('hit_description', line[6]),
                    ('hit_length', line[7]),
                    ('hit_start', line[8]),
                    ('hit_end', line[9]),
                    ('group', line[10]),
                    ('group_description', line[11]),
                    ('group_categories', line[12]),
                    ('region_start', line[13]),
                    ('region_end', line[14]),
                    ('region_coverage', line[15]),
                    ('region_identities', line[16]),
                    ('region_positives', line[17])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_cog_table']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入cog注释table信息：%s出错!" , variables=(table_path), code="53703765")
        else:
            self.bind_object.logger.info("导入cog注释table信息：%s成功!" % (table_path))

    def run_webroot(self, result_dir, trans2gene, trans2gene_ref, params_dict, task_id, stat_id, last_id, taxonomy, exp_level):
        """
        用于注释重运行导表
        result_dir: 新序列注释的结果文件夹
        trans2gene: 转录本基因对应关系文件
        params_dict: 参数
        stat_id: 统计结果主表ID
        last_id: 上次重运行ID
        taxonomy: 物种分类
        """
        self.bind_object.logger.info("开始导表webroot数据路径为 {}".format(result_dir))
        self.set_result_dir(result_dir)
        self.set_relation(trans2gene, trans2gene_ref)
        # self.get_trans2gene(trans2gene)
        self.bind_object.logger.info("开始导表task_id为 {}".format(self.task_id))

        stat_id = ObjectId(stat_id)
        if last_id:
            last_id = ObjectId(last_id)
        else:
            pass
        self.task_id = task_id
        #task_id = self.task_id
        self.bind_object.logger.info("开始导表task_id为 {}".format(task_id))
        # task_id = "denovo_rna_v2"
        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))

        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['ref_stat_path'], root_dir=self.result_file['ref_venn_path'], seq_type ="ref", exp_level=exp_level)
        if self.has_new:
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['new_stat_path'], root_dir=self.result_file['new_venn_path'], seq_type ="new", exp_level=exp_level)
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['all_stat_path'], root_dir=self.result_file['all_venn_path'], seq_type ="all", exp_level=exp_level)

        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)
        if last_id:
            self.bind_object.logger.info("删除表格为 {}".format(last_id))
            self.remove_table_by_main_record(main_table='sg_annotation_stat', _id=last_id, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
            self.bind_object.logger.info("删除表格成功 {}".format(last_id))

        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))

        nr_old_id = self.get_table_by_main_record(main_table='sg_annotation_nr', task_id=task_id, type=self.annot_type)
        self.bind_object.logger.info("查找表格{}".format(nr_old_id))
        nr_id = self.add_annotation_nr(name=None, params=params_select_nr, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="new", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path"], exp_level=exp_level)
        self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="ref", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path" + "ref"], exp_level=exp_level)
        self.bind_object.logger.info("插入新表 {}".format(nr_id))
        self.update_db_record('sg_annotation_nr', nr_id, status="end", main_id=nr_id)


        if nr_old_id:
            self.bind_object.logger.info( "删除旧表 {}".format(nr_old_id['_id']))
            self.remove_table_by_main_record(main_table='sg_annotation_nr', _id=nr_old_id['_id'], detail_table=['sg_annotation_nr_detail', 'sg_annotation_nr_pie'], detail_table_key='nr_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        swissprot_old_id = self.get_table_by_main_record(main_table='sg_annotation_swissprot', task_id=task_id, type=self.annot_type)
        self.bind_object.logger.info("查找表格成{}".format(swissprot_old_id))
        params_select_swissprot = dict([(k,params_dict.get(k,None)) for k in ('swissprot_evalue', 'swissprot_similarity', 'swissprot_identity')])
        params_select_swissprot = json.dumps(params_select_swissprot, sort_keys=True, separators=(',', ':'))
        swissprot_id = self.add_annotation_swissprot(name=None, params=params_select_swissprot, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="new", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path"], exp_level=exp_level)
        self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="ref", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path"  + "ref"], exp_level=exp_level)
        self.update_db_record('sg_annotation_swissprot', swissprot_id, status="end", main_id=swissprot_id)
        if swissprot_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_swissprot', _id=swissprot_old_id['_id'],  detail_table=['sg_annotation_swissprot_detail', 'sg_annotation_swissprot_pie'], detail_table_key='swissprot_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        pfam_old_id = self.get_table_by_main_record(main_table='sg_annotation_pfam', task_id=task_id, type=self.annot_type)
        params_select_pfam = dict([('pfam_evalue',params_dict['pfam_evalue'])])
        params_select_pfam = json.dumps(params_select_pfam, sort_keys=True, separators=(',', ':'))
        pfam_id = self.add_annotation_pfam(name=None, params=params_select_pfam, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T", exp_level=exp_level)
        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'  + "ref"], seq_type="ref", anno_type="T", exp_level=exp_level)
        #self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path'], seq_type="new", anno_type="G")
        if exp_level.lower() == "transcript":
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path' + "ref"], seq_type="ref", anno_type="T")
        self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path' + "ref"], seq_type="ref", anno_type="G")
        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path' + "all"], seq_type="all", anno_type="T")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=self.result_file['gene_pfam_path' + "all"], seq_type="all", anno_type="G")

        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_pfam', pfam_id, status="end", main_id=pfam_id)

        if pfam_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_pfam', _id=pfam_old_id['_id'], detail_table=['sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], detail_table_key='pfam_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        query_old_id = self.get_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.annot_type)
        self.remove_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.annot_type, detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_query_detail(query_id=query_id, query_path=self.result_file['query_path'], seq_type ="new", anno_type="T", exp_level=exp_level)
        self.add_annotation_query_detail(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type ="ref", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)

        if query_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_query', _id=query_old_id['_id'], detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        cog_old_id = self.get_table_by_main_record(main_table='sg_annotation_cog', task_id=task_id, type=self.annot_type)
        params_select_cog = dict([(k,params_dict.get(k,None)) for k in ('cog_evalue', 'cog_similarity', 'cog_identity')])
        params_select_cog = json.dumps(params_select_cog, sort_keys=True, separators=(',', ':'))
        cog_id = self.add_annotation_cog(name=None, params=params_select_cog, result_dir=self.result_dir)
        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_sum_path'], seq_type="new", anno_type="T")
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path'], seq_type="new", anno_type="G")
        #r_sum_path = ref_anno_path + "/cog/cog_summary.xls"
        if exp_level.lower() == "transcript":
            self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_sum_path' + 'ref'], seq_type="ref", anno_type="T")
        #r_gene_sum_path = ref_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=self.result_file['n_gene_sum_path' + 'ref'], seq_type="ref", anno_type="G")
        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_cog_detail_all(cog_id=cog_id, ref_cog_path=self.result_file['n_sum_path' + 'ref'], new_cog_path=self.result_file['n_sum_path'], seq_type="all", anno_type="T")
            self.add_annotation_cog_detail_all(cog_id=cog_id, ref_cog_path=self.result_file['n_gene_sum_path' + 'ref'], new_cog_path=self.result_file['n_gene_sum_path'], seq_type="all", anno_type="G")
        self.update_db_record('sg_annotation_cog', cog_id, status="end", main_id=cog_id)
        if cog_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_cog', _id=cog_old_id['_id'], detail_table=['sg_annotation_cog_detail'], detail_table_key='cog_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        go_old_id = self.get_table_by_main_record(main_table='sg_annotation_go', task_id=task_id, type=self.annot_type)

        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, level_path=self.result_file['stat_level3'])
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, level_path=self.result_file['stat_level4'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, level_path=self.result_file['stat_level3'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, level_path=self.result_file['stat_level4'])
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", list_path=self.result_file['gos_path'])
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, level_path=self.result_file['gene_stat_level3'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, level_path=self.result_file['gene_stat_level4'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, level_path=self.result_file['gene_stat_level3'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, level_path=self.result_file['gene_stat_level4'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", list_path=self.result_file['gene_gos_path'])
        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, level_path=self.result_file['stat_level3_ref'])
            self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, level_path=self.result_file['stat_level4_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, level_path=self.result_file['stat_level3_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, level_path=self.result_file['stat_level4_ref'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", list_path=self.result_file['gos_path_ref'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, level_path=self.result_file['gene_stat_level3_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, level_path=self.result_file['gene_stat_level4_ref'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, level_path=self.result_file['gene_stat_level3_ref'])
        self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, level_path=self.result_file['gene_stat_level4_ref'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", list_path=self.result_file['gene_gos_path_ref'])

        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=2, level_path=self.result_file['stat_level2_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=3, level_path=self.result_file['stat_level3_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=4, level_path=self.result_file['stat_level4_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=3, level_path=self.result_file['gene_stat_level3_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=4, level_path=self.result_file['gene_stat_level4_all'])
        self.update_db_record('sg_annotation_go', go_id, status="end", main_id=go_id)
        if go_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_go', _id=go_old_id['_id'], detail_table=['sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], detail_table_key='go_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        kegg_old_id = self.get_table_by_main_record(main_table='sg_annotation_kegg', task_id=task_id, type=self.annot_type)
        params_select_kegg = dict([(k,params_dict.get(k,None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=params_select_kegg, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", layer_path=self.result_file['layer_path'])
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", pathway_file=self.result_file['pathway_path'], figure_dir=self.result_file['png_path'])
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path'])
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", pathway_file=self.result_file['pathway_path'], pathway_dir=self.result_file['png_path'])
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", layer_path=self.result_file['gene_layer_path'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", pathway_file=self.result_file['gene_pathway_path'], figure_dir=self.result_file['gene_png_path'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", pathway_file=self.result_file['gene_pathway_path'], pathway_dir=self.result_file['gene_png_path'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path'])

        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", layer_path=self.result_file['layer_path_ref'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", pathway_file=self.result_file['pathway_path_ref'], figure_dir=self.result_file['png_path_ref'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", pathway_file=self.result_file['pathway_path_ref'], pathway_dir=self.result_file['png_path_ref'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path_ref'])
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", layer_path=self.result_file['gene_layer_path_ref'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", pathway_file=self.result_file['gene_pathway_path_ref'], pathway_dir=self.result_file['gene_png_path_ref'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", pathway_file=self.result_file['gene_pathway_path_ref'], figure_dir=self.result_file['gene_png_path_ref'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path_ref'])

        if self.has_new:
            r_cate_path = self.result_file['layer_path_ref']
            n_cate_path = self.result_file['layer_path']
            r_gene_cate_path = self.result_file['gene_layer_path_ref']
            n_gene_cate_path = self.result_file['gene_layer_path']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="T", ref_layer_path=r_cate_path, new_layer_path=n_cate_path)
            self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="G", ref_layer_path=r_gene_cate_path, new_layer_path=n_gene_cate_path)
            pathway_path = self.result_file['pathway_path' + '_all']
            png_path = self.result_file['png_path' + '_all']
            gene_pathway_path = self.result_file['gene_pathway_path' + '_all']
            gene_png_path =  self.result_file['gene_png_path' + '_all']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="T", pathway_file=pathway_path, figure_dir=png_path)
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="T", pathway_file=pathway_path, pathway_dir=png_path)
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="G", pathway_file=gene_pathway_path, figure_dir=gene_png_path)
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="G", pathway_file=gene_pathway_path, pathway_dir=gene_png_path)

        self.update_db_record('sg_annotation_kegg', kegg_id, status="end", main_id=kegg_id)
        if kegg_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_kegg', _id=kegg_old_id['_id'], detail_table=['sg_annotation_kegg_categories', 'sg_annotation_kegg_level','sg_annotation_kegg_table'], detail_table_key='kegg_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

    def add_annotation(self, name=None, params=None, ref_anno_path=None, new_anno_path=None, pfam_path=None, merge_tran_output=None, merge_gene_output=None):
        """
        ref_anno_path: 已知序列注释的结果文件夹
        new_anno_path: 新序列注释的结果文件夹
        pfam_path:转录本的pfam_domain
        merge_tran_output: 转录本的merge_annot tool输出结果路径
        merge_gene_output: 基因的merge_annot tool输出结果路径
        """
        new_stat_path = new_anno_path + "/anno_stat/all_annotation_statistics.xls"
        new_venn_path = new_anno_path + "/anno_stat/venn"
        stat_id = self.add_annotation_stat(name=None, params=params, seq_type="new", database="nr,swissprot,pfam,cog,go,kegg")
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=new_stat_path, root_dir=new_venn_path)
        blast_id = self.add_annotation_blast(name=None, params=params, stat_id=stat_id)
        blast_path = new_anno_path + "/anno_stat/blast"
        if os.path.exists(blast_path):
            for db in ["nr", "swissprot"]:
                trans_blast_path = blast_path + "/" + db + '.xls'
                gene_blast_path = blast_path + '/gene_' + db + '.xls'
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="transcript", database=db, blast_path=trans_blast_path)
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="gene", database=db, blast_path=gene_blast_path)
        else:
            self.bind_object.set_error("没有blast的结果文件夹", code="53703766")
        nr_id = self.add_annotation_nr(name=None, params=params, stat_id=stat_id)
        nr_path = new_anno_path + "/anno_stat/blast_nr_statistics"
        if os.path.exists(nr_path):
            evalue_path = nr_path + "/nr_evalue.xls"
            similar_path = nr_path + "/nr_similar.xls"
            gene_evalue_path = nr_path + "/gene_nr_evalue.xls"
            gene_similar_path = nr_path + "/gene_nr_similar.xls"
            self.add_annotation_nr_pie(nr_id=nr_id, evalue_path=evalue_path, similar_path=similar_path, seq_type="new", anno_type="transcript")
            self.add_annotation_nr_pie(nr_id=nr_id, evalue_path=gene_evalue_path, similar_path=gene_similar_path,  seq_type="new", anno_type="gene")
        else:
            self.bind_object.set_error("新序列NR注释结果文件不存在", code="53703767")
        swissprot_id = self.add_annotation_swissprot(name=None, params=params, stat_id=stat_id)
        swissprot_path = new_anno_path + "/anno_stat/blast_swissprot_statistics"
        if os.path.exists(swissprot_path):
            evalue_path = swissprot_path + "/swissprot_evalue.xls"
            similar_path = swissprot_path + "/swissprot_similar.xls"
            gene_evalue_path = swissprot_path + "/gene_swissprot_evalue.xls"
            gene_similar_path = swissprot_path + "/gene_swissprot_similar.xls"
            self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=evalue_path, similar_path=similar_path, seq_type="new", anno_type="transcript")
            self.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=gene_evalue_path, similar_path=gene_similar_path, seq_type="new", anno_type="gene")
        else:
            self.bind_object.set_error("新序列Swiss-Prot注释结果文件不存在", code="53703768")
        pfam_id = self.add_annotation_pfam(name=None, params=params, stat_id=stat_id)
        gene_pfam_path = new_anno_path + "/anno_stat/pfam_stat/gene_pfam_domain"
        if os.path.exists(pfam_path) and os.path.exists(gene_pfam_path):
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
        else:
            self.bind_object.set_error("pfam注释结果文件不存在", code="53703769")
        ref_stat_path = ref_anno_path + "/anno_stat/all_annotation_statistics.xls"
        ref_venn_path = ref_anno_path + "/anno_stat/venn"
        if os.path.exists(ref_stat_path) and os.path.exists(ref_venn_path):
            stat_id = self.add_annotation_stat(name=None, params=params, seq_type="ref" , database="cog,go,kegg")
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=ref_stat_path, root_dir=ref_venn_path)
        else:
            self.bind_object.set_error("已知序列注释统计文件和venn图文件夹不存在", code="53703770")
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id)
        query_path = ref_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = ref_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail_webroot(query_id=query_id, query_path=query_path, anno_type="transcript")
        self.add_annotation_gene_query_detail(query_id=query_id, query_path=gene_query_path, anno_type="gene")
        query_path = new_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = new_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail_webroot(query_id=query_id, query_path=query_path, anno_type="transcript")
        self.add_annotation_gene_query_detail(query_id=query_id, query_path=gene_query_path, anno_type="gene")
        cog_id = self.add_annotation_cog(name=name, params=params)
        r_sum_path = ref_anno_path + "/cog/cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=r_sum_path, seq_type="ref", anno_type="transcript")
        r_gene_sum_path = ref_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=r_gene_sum_path, seq_type="ref", anno_type="gene")
        n_sum_path = new_anno_path + "/cog/cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=n_sum_path, seq_type="new", anno_type="transcript")
        n_gene_sum_path = new_anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        self.add_annotation_cog_detail(cog_id=cog_id, cog_path=n_gene_sum_path, seq_type="new", anno_type="gene")
        self.add_annotation_cog_detail_all(cog_id=cog_id, ref_cog_path=r_sum_path, new_cog_path=n_sum_path, seq_type="all", anno_type="transcript")
        self.add_annotation_cog_detail_all(cog_id=cog_id, ref_cog_path=r_gene_sum_path, new_cog_path=n_gene_sum_path, seq_type="all", anno_type="gene")

        def add_go(go_id, go_path, gene_go_path, seq_type):
            if os.path.exists(go_path) and os.path.exists(gene_go_path):
                stat_level2 = go_path + "/go12level_statistics.xls"
                stat_level3 = go_path + "/go123level_statistics.xls"
                stat_level4 = go_path + "/go1234level_statistics.xls"
                gene_stat_level2 = gene_go_path + "/gene_go12level_statistics.xls"
                gene_stat_level3 = gene_go_path + "/gene_go123level_statistics.xls"
                gene_stat_level4 = gene_go_path + "/gene_go1234level_statistics.xls"
                gos_path = go_path + "/query_gos.list"
                gene_gos_path = gene_go_path + "/gene_gos.list"
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, level_path=stat_level2)
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, level_path=gene_stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, level_path=stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, level_path=stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, level_path=stat_level4)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, level_path=gene_stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, level_path=gene_stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, level_path=gene_stat_level4)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, level_path=stat_level2)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, level_path=stat_level3)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, level_path=stat_level4)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, level_path=gene_stat_level2)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, level_path=gene_stat_level3)
                self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, level_path=gene_stat_level4)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="transcript", list_path=gos_path)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="gene", list_path=gene_gos_path)
            else:
                self.bind_object.set_error("GO注释的结果文件不存在", code="53703771")
        go_id = self.add_annotation_go(name=name, params=params)
        go_path = ref_anno_path + "/go"
        gene_go_path = ref_anno_path + "/anno_stat/go_stat"
        add_go(go_id=go_id, go_path=go_path, gene_go_path=gene_go_path, seq_type="ref")
        go_path = new_anno_path + "/go"
        gene_go_path = new_anno_path + "/anno_stat/go_stat"
        add_go(go_id=go_id, go_path=go_path, gene_go_path=gene_go_path, seq_type="new")
        r_stat_level2 = ref_anno_path + "/go/go12level_statistics.xls"
        r_stat_level3 = ref_anno_path + "/go/go123level_statistics.xls"
        r_stat_level4 = ref_anno_path + "/go/go1234level_statistics.xls"
        n_stat_level2 = new_anno_path + "/go/go12level_statistics.xls"
        n_stat_level3 = new_anno_path + "/go/go123level_statistics.xls"
        n_stat_level4 = new_anno_path + "/go/go1234level_statistics.xls"
        r_gene_stat_level2 = ref_anno_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        r_gene_stat_level3 = ref_anno_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        r_gene_stat_level4 = ref_anno_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        n_gene_stat_level2 = new_anno_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        n_gene_stat_level3 = new_anno_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        n_gene_stat_level4 = new_anno_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=2, level_path=r_stat_level2, n_go_path=n_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=3, level_path=r_stat_level3, n_go_path=n_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=4, level_path=r_stat_level4, n_go_path=n_stat_level4)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=2, level_path=r_gene_stat_level2, n_go_path=n_gene_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=3, level_path=r_gene_stat_level3, n_go_path=n_gene_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=4, level_path=r_gene_stat_level4, n_go_path=n_gene_stat_level4)

        def add_kegg(kegg_id, kegg_path, gene_kegg_path, seq_type):
            if os.path.exists(kegg_path) and os.path.exists(gene_kegg_path):
                layer_path = kegg_path + "/kegg_layer.xls"
                pathway_path = kegg_path + "/pathway_table.xls"
                png_path = kegg_path + "/pathways"
                table_path = kegg_path + "/kegg_table.xls"
                gene_layer_path = gene_kegg_path + "/gene_kegg_layer.xls"
                gene_pathway_path = gene_kegg_path + "/gene_pathway_table.xls"
                gene_png_path = gene_kegg_path + "/gene_pathway"
                gene_table_path = gene_kegg_path + "/gene_kegg_table.xls"
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", layer_path=layer_path)
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", layer_path=gene_layer_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", pathway_file=pathway_path, figure_dir=png_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", pathway_file=gene_pathway_path, figure_dir=gene_png_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", table_path=table_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", table_path=gene_table_path)
            else:
                self.bind_object.set_error("KEGG注释文件不存在", code="53703772")
        kegg_id = self.add_annotation_kegg(name=None, params=params)
        kegg_path = ref_anno_path + "/kegg"
        gene_kegg_path = ref_anno_path + "/anno_stat/kegg_stat"
        add_kegg(kegg_id=kegg_id, kegg_path=kegg_path, gene_kegg_path=gene_kegg_path, seq_type="ref")
        kegg_path = new_anno_path + "/kegg"
        gene_kegg_path = new_anno_path + "/anno_stat/kegg_stat"
        add_kegg(kegg_id=kegg_id, kegg_path=kegg_path, gene_kegg_path=gene_kegg_path, seq_type="new")
        r_cate_path = ref_anno_path + "/kegg/kegg_layer.xls"
        n_cate_path = new_anno_path + "/kegg/kegg_layer.xls"
        r_gene_cate_path = ref_anno_path + "/anno_stat/kegg_stat/gene_kegg_layer.xls"
        n_gene_cate_path = new_anno_path + "/anno_stat/kegg_stat/gene_kegg_layer.xls"
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="transcript", ref_layer_path=r_cate_path, new_layer_path=n_cate_path)
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="gene", ref_layer_path=r_gene_cate_path, new_layer_path=n_gene_cate_path)
        pathway_path = merge_tran_output + "/pathway_table.xls"
        png_path = merge_tran_output + "/all_pathways"
        gene_pathway_path = merge_gene_output + "/pathway_table.xls"
        gene_png_path = merge_gene_output + "/all_pathways"
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="transcript", pathway_file=pathway_path, figure_dir=png_path)
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="gene", pathway_file=gene_pathway_path, figure_dir=gene_png_path)

    @report_check
    def add_stat_detail(self, old_stat_id, stat_id, nr_evalue, gene_nr_evalue, sw_evalue, gene_sw_evalue):
        """
        注释重运行时注释统计导表sg_annotation_stat_detail
        """
        if not isinstance(old_stat_id, ObjectId):
            if isinstance(old_stat_id, types.StringTypes):
                old_stat_id = ObjectId(old_stat_id)
            else:
                self.bind_object.set_error('old_stat_id必须为ObjectId对象或其对应的字符串！', code="53703773")
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="53703774")
        collection = self.db["sg_annotation_stat_detail"]
        results = collection.find({"stat_id": old_stat_id})
        data_list, data = [], []
        for result in results:
            db = result["type"]
            if db == "total":
                total_tran = result["transcript"]
                total_gene = result["gene"]
            # 增加数据库类型，区分有注释的基因和有分类的基因
            # if db in ["pfam", "total_anno", "total_anno_nsp","total_class", "total"]:
            if db in ["pfam", "total_anno",  "total"]:
                data = [
                    ('stat_id', stat_id),
                    ('type', result["type"]),
                    ('transcript', result["transcript"]),
                    ('gene', result["gene"]),
                    ('transcript_percent', result["transcript_percent"]),
                    ('gene_percent', result["gene_percent"]),
                    ('gene_list', result["gene_list"]),
                    ('transcript_list', result["transcript_list"])
                ]
                data = SON(data)
                data_list.append(data)
        nr_ids = self.stat(stat_path=nr_evalue)
        gene_nr_ids = self.stat(stat_path=gene_nr_evalue)
        data = [
            ('stat_id', stat_id),
            ('type', "nr"),
            ('transcript', len(nr_ids)),
            ('gene', len(gene_nr_ids)),
            ('transcript_percent', round(float(len(nr_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_nr_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_nr_ids)),
            ('transcript_list', ",".join(nr_ids))
        ]
        data = SON(data)
        data_list.append(data)
        sw_ids = self.stat(stat_path=sw_evalue)
        gene_sw_ids = self.stat(stat_path=gene_sw_evalue)
        data = [
            ('stat_id', stat_id),
            ('type', "swissprot"),
            ('transcript', len(sw_ids)),
            ('gene', len(gene_sw_ids)),
            ('transcript_percent', round(float(len(sw_ids))/total_tran, 4)),
            ('gene_percent', round(float(len(gene_sw_ids))/total_gene, 4)),
            ('gene_list', ",".join(gene_sw_ids)),
            ('transcript_list', ",".join(sw_ids))
        ]
        data = SON(data)
        data_list.append(data)
        try:
            collection = self.db['sg_annotation_stat_detail']
            collection.insert_many(data_list)
        except:
            self.bind_object.set_error("导入注释统计信息出错", code="53703775")
        else:
            self.bind_object.logger.info("导入注释统计信息成功")

    def stat(self, stat_path):
        with open(stat_path, "rb") as f:
            id_list = []
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                q_id = line[5]
                id_list.append(q_id)
        id_list = list(set(id_list))
        return id_list

    @report_check
    def add_annotation_blast(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="53703776")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationBlast_' + self.annot_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'blast最佳比对结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_blast']
        blast_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_blast!")
        return blast_id

    @report_check
    def add_annotation_blast_detail(self, blast_id, seq_type, anno_type, database, blast_path):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="53703777")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="53703778")
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'blast_id': blast_id,
                        'seq_type': seq_type,
                        'anno_type': anno_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'query_id': line[5],
                        'q_len': int(line[6]),
                        'q_begin': line[7],
                        'q_end': line[8],
                        'q_frame': line[9],
                        'hit_name': line[10],
                        'hit_len': int(line[11]),
                        'hsp_begin': line[12],
                        'hsp_end': line[13],
                        'hsp_frame': line[14],
                        'description': line[15]
                    }
                    if  anno_type == 'T':
                        try:
                            data.update({'gene_id': self.t2g_dict[line[5]]})
                            data.update({'is_gene': self.t2r_dict[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    collection = self.db['sg_annotation_blast_detail']
                    collection.insert_one(data).inserted_id
        self.bind_object.logger.info("导入blast信息：%s成功!" % (blast_path))

    @report_check
    def add_annotation_blast_nr_detail(self, blast_id, seq_type, anno_type, database, blast_path, exp_level):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="53703779")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="53703780")
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'nr_id': blast_id,
                        'seq_type': seq_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'transcript_id': line[5],
                        'q_len': int(line[6]),
                        'hit_name': line[10],
                        'description': line[15]
                    }
                    if anno_type == 'T':
                        try:
                            data.update({'gene_id': self.t2g_dict[line[5]]})
                            data.update({'is_gene': self.t2r_dict[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    else:
                        pass
                    if exp_level.lower() == "gene" and self.t2r_dict[line[5]] == False:
                        continue
                    else:
                        data_list.append(SON(data))
        collection = self.db['sg_annotation_nr_detail']
        if data_list:
            try:
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入nrblast信息：%s成功!" % (blast_path))
            except Exception, e:
                self.bind_object.set_error("导入注释统计信息失败:%s" , variables=( blast_path), code="53703781")

    @report_check
    def add_annotation_blast_swissprot_detail(self, blast_id, seq_type, anno_type, database, blast_path, exp_level):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="53703782")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="53703783")
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip().split('\t')
                query_name = line[5]
                hit_name = line[10]
                if flag == query_name:
                    pass
                else:
                    flag = query_name
                    data = {
                        'swissprot_id': blast_id,
                        'seq_type': seq_type,
                        # 'anno_type': anno_type,
                        'database': database,
                        'score': float(line[0]),
                        'e_value': float(line[1]),
                        'hsp_len': int(line[2]),
                        'identity_rate': round(float(line[3]), 4),
                        'similarity_rate': round(float(line[4]), 4),
                        'transcript_id': line[5],
                        'q_len': int(line[6]),
                        # 'q_begin': line[7],
                        # 'q_end': line[8],
                        # 'q_frame': line[9],
                        'hit_name': line[10],
                        # 'hit_len': int(line[11]),
                        # 'hsp_begin': line[12],
                        # 'hsp_end': line[13],
                        # 'hsp_frame': line[14],
                        'description': line[15]
                    }
                    if anno_type == 'T':
                        try:
                            data.update({'gene_id': self.t2g_dict[line[5]]})
                            data.update({'is_gene': self.t2r_dict[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    else:
                        pass
                    if exp_level.lower() == "gene" and self.t2r_dict[line[5]] == False:
                        continue
                    else:
                        data_list.append(SON(data))
        collection = self.db['sg_annotation_swissprot_detail']
        if data_list:
            try:
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入swissprot blast信息：%s成功!" % (blast_path))
            except Exception, e:
                self.bind_object.set_error("导入swissprot注释统计信息失败:%s" , variables=( blast_path), code="53703784")

    @report_check
    def add_annotation_nr(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="53703785")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationNr_' + self.annot_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.annot_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'nr注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_nr']
        nr_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_nr!")
        return nr_id

    @report_check
    def add_annotation_nr_pie(self, nr_id, species_path, evalue_path, similar_path, seq_type, anno_type):
        if not isinstance(nr_id, ObjectId):
            if isinstance(nr_id, types.StringTypes):
                nr_id = ObjectId(nr_id)
            else:
                self.bind_object.set_error('nr_id必须为ObjectId对象或其对应的字符串！', code="53703786")
        if not os.path.exists(evalue_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="53703787")
        if not os.path.exists(similar_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="53703788")
        evalue, evalue_list,species_list_a, species_a, similar, similar_list = [], [], [], [], [], []
        data_list = []
        with open(evalue_path, "r") as f1, open(similar_path, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            for line1 in lines1[1:]:
                line1 = line1.strip().split('\t')
                value = {"key": line1[0], "value": int(line1[1])}
                try:
                    value_list = {"key": line1[0], "value": line1[2]}
                except:
                    value_list = {"key": line1[0], "value": None}
                evalue.append(value)
                evalue_list.append(value_list)
            for line2 in lines2[1:]:
                line2 = line2.strip().split('\t')
                similarity = {"key": line2[0], "value": int(line2[1])}
                try:
                    similarity_list = {"key": line2[0], "value": line2[2]}
                except:
                    similarity_list = {"key": line2[0], "value": None}
                similar.append(similarity)
                similar_list.append(similarity_list)

        with open(species_path, "r") as f3:
            lines3 = f3.readlines()
            for line3 in lines3[1:15]:
                line3 = line3.strip().split("\t")
                if anno_type == "T":
                    species = {"key": line3[0], "value": int(line3[1])}
                    try:
                        species_list = {"key": line3[0], "value": line3[5]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                elif anno_type == "G":
                    species = {"key": line3[0], "value": int(line3[2])}
                    try:
                        species_list = {"key": line3[0], "value": line3[6]}
                    except:
                        species_list = {"key": line3[0], "value": None}
                species_a.append(species)
                species_list_a.append(species_list)
            other = 0
            for line3 in lines3[16:]:
                line3 = line3.strip().split("\t")
                if anno_type == "T":
                    other += int(line3[1])
                elif anno_type == "G":
                    other += int(line3[2])
            species = {"key": 'other', "value": other}
            species_list = {"key": 'other', "value": None}
            species_a.append(species)
            species_list_a.append(species_list)

            data = [
                ('nr_id', nr_id),
                #('seq_type', seq_type),
                ('anno_type', anno_type),
                ('e_value', evalue),
                ('similar', similar),
                ('evalue_list', evalue_list),
                ('similar_list', similar_list),
                ('species', species_a),
                ('species_list', species_list_a),
            ]

        data = SON(data)
        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_nr_pie']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入nr库注释作图信息evalue,similar：%s、%s出错!" , variables=(evalue_path, similar_path), code="53703789")
            else:
                self.bind_object.logger.info("导入nr库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))

    @report_check
    def add_annotation_swissprot(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="53703790")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationSwissprot_' + self.annot_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.annot_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'swissprot注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_swissprot']
        swissprot_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_swissprot!")
        return swissprot_id

    @report_check
    def add_annotation_swissprot_pie(self, swissprot_id, evalue_path, similar_path, seq_type, anno_type):
        """
        """
        if not isinstance(swissprot_id, ObjectId):
            if isinstance(swissprot_id, types.StringTypes):
                swissprot_id = ObjectId(swissprot_id)
            else:
                self.bind_object.set_error('swissprot_id必须为ObjectId对象或其对应的字符串！', code="53703791")
        if not os.path.exists(evalue_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="53703792")
        if not os.path.exists(similar_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="53703793")
        evalue, evalue_list, similar, similar_list = [], [], [], []
        data_list = []
        with open(evalue_path, "r") as f1, open(similar_path, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            for line1 in lines1[1:]:
                line1 = line1.strip().split('\t')
                value = {"key": line1[0], "value": int(line1[1])}
                try:
                    value_list = {"key": line1[0], "value": line1[2]}
                except:
                    value_list = {"key": line1[0], "value": None}
                evalue.append(value)
                evalue_list.append(value_list)
            for line2 in lines2[1:]:
                line2 = line2.strip().split('\t')
                similarity = {"key": line2[0], "value": int(line2[1])}
                try:
                    similarity_list = {"key": line2[0], "value": line2[2]}
                except:
                    similarity_list = {"key": line2[0], "value": None}
                similar.append(similarity)
                similar_list.append(similarity_list)
        data = [
            ('swissprot_id', swissprot_id),
            # ('seq_type', seq_type),
            ('anno_type', anno_type),
            ('e_value', evalue),
            ('similar', similar),
            ('evalue_list', evalue_list),
            ('similar_list', similar_list),
        ]
        data = SON(data)
        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_swissprot_pie']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入swissprot库注释作图信息evalue,similar：%s、%s出错!" , variables=(evalue_path, similar_path), code="53703794")
            else:
                self.bind_object.logger.info("导入swissprot库注释作图信息evalue,similar：%s、%s成功!" % (evalue_path, similar_path))

    @report_check
    def add_annotation_pfam(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="53703795")
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationPfam_' + self.annot_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.annot_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'pfam注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_pfam']
        pfam_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_annotation_pfam!")
        return pfam_id

    def add_annotation_pfam_detail(self, pfam_id, pfam_path, seq_type, anno_type, exp_level):
        """
        pfam_path: pfam_domain
        """
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="53703796")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="53703797")
        data_list = []
        with open(pfam_path, "r") as f:
            lines = f.readlines()
            last_seq_id = ''
            last_pfam_id = ''
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[2] != last_seq_id or line[3] != last_pfam_id:
                    # 过滤不在参考范围的转录本注释
                    if exp_level.lower() == "gene" or line[0] not in self.t2r_dict:
                        continue
                    data = [
                        ('pfam_id', pfam_id),
                        ('seq_type', seq_type),
                        # ('anno_type', anno_type),
                        ('transcript_id', line[0]),
                        ('pfam', line[2]),
                        ('domain', line[3]),
                        ('description', line[4]),
                        ('protein_id', line[1]),
                        ('e_value', float(line[9])),
                        ('length', int(line[6])-int(line[5])),
                        ('protein_start', int(line[5])),
                        ('protein_end', int(line[6])),
                        ('pfam_start', int(line[7])),
                        ('pfam_end', int(line[8])),
                    ]
                    if  anno_type == 'T':
                        try:
                            data.append(('gene_id', self.t2g_dict[line[0]]))
                            data.append(('is_gene', self.t2r_dict[line[0]]))
                        except:
                            data.append(('gene_id', line[0]))
                            data.append(('is_gene', True))
                        data = SON(data)
                        data_list.append(data)
                    else:
                        data = SON(data)
                        data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_pfam_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="53703798")
            else:
                self.bind_object.logger.info("导入pfam注释信息:%s成功" % pfam_path)

    @report_check
    def add_annotation_pfam_bar(self, pfam_id, pfam_path, seq_type, anno_type):
        pfam = []
        domain = {}
        with open(pfam_path, "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[3] not in pfam:
                    pfam.append(line[3])
                    domain[line[3]] = 1
                else:
                    domain[line[3]] += 1
        if not isinstance(pfam_id, ObjectId):
            if isinstance(pfam_id, types.StringTypes):
                pfam_id = ObjectId(pfam_id)
            else:
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="53703799")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="537037100")
        data_list = []
        dom = zip(domain.values(), domain.keys())
        dom_sort = sorted(dom, reverse=True)
        for num,dom in dom_sort:
            data = [
                ('pfam_id', pfam_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('domain', dom),
                ('num', num)
            ]
            data = SON(data)
            data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_pfam_bar']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="537037101")
            else:
                self.bind_object.logger.info("导入pfam注释信息:%s成功！" % pfam_path)

    @report_check
    def add_annotation_query_detail_webroot(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537037102")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537037103")
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('anno_type', anno_type),
                    ('transcript_id', line[0]),
                    ('gene_id', line[1]),
                ]
                try:
                    data.append(('gene_name', line[2]))
                except:
                    data.append(('gene_name', None))
                try:
                    data.append(('length', line[3]))
                except:
                    data.append(('length', None))
                try:
                    data.append(('cog', line[4]))
                    data.append(('cog_description', line[6]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[5]))
                    data.append(('nog_description', line[7]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[8]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[9]))
                except:
                    data.append(('ko_name', None))

                try:
                    data.append(('pathways', line[10]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[11]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[12]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[13]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[14]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入转录本注释统计信息：%s出错!" , variables=(query_path), code="537037104")
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537037105")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537037106")
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                line = lines[j].strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('gene_id', line[0]),
                    ('anno_type', anno_type),
                ]
                try:
                    data.append(('gene_name', line[1]))
                except:
                    data.append(('gene_name', None))
                try:
                    data.append(('cog', line[2]))
                    data.append(('cog_description', line[4]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[3]))
                    data.append(('nog_description', line[5]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[6]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[7]))
                except:
                    data.append(('ko_name', None))
                try:
                    data.append(('pathways', line[8]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[9]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[10]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[11]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[12]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入基因注释查询信息：%s出错!" , variables=(query_path), code="537037107")
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_denovo_detail(self, query_id, query_path, seq_type, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537037108")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537037109")
        data_list = []
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                line = lines[j].strip().split('\t')
                data = [
                    ('query_id', query_id),
                    ('gene_id', line[0]),
                    ('tran_id', line[1]),
                   # if $line[2] == "yes":
                   #     ('is_gene', True),
                    ('anno_type', anno_type),
                ]

                try:
                    data.append(('cog', line[1]))
                    data.append(('cog_description', line[3]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                try:
                    data.append(('nog', line[2]))
                    data.append(('nog_description', line[4]))
                except:
                    data.append(('nog', None))
                    data.append(('nog_description', None))
                try:
                    data.append(('ko_id', line[5]))
                except:
                    data.append(('ko_id', None))
                try:
                    data.append(('ko_name', line[6]))
                except:
                    data.append(('ko_name', None))
                try:
                    data.append(('pathways', line[7]))
                except:
                    data.append(('pathways', None))
                try:
                    data.append(('pfam', line[8]))
                except:
                    data.append(('pfam', None))
                try:
                    data.append(('go', line[9]))
                except:
                    data.append(('go', None))
                try:
                    data.append(('nr', line[10]))
                except:
                    data.append(('nr', None))
                try:
                    data.append(('swissprot', line[11]))
                except:
                    data.append(('swissprot', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入基因注释查询信息：%s出错!" , variables=(query_path), code="537037110")
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'ref_rna_v2.refrna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wheet)
        wf.sheet.id = 'tsg_34573'
        wf.sheet.project_sn = '188_5d12fb95c7db8'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('ref_rna_v2.annotation')
        annot_merge_output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20190626/Refrna_tsg_34573/AnnotMerge/output'
        params_dict = {
            'task_id': 'tsg_34573',
            'submit_location': 'annotationstat',
            'task_type': 2,
            'nr_evalue': 1e-5,
            'nr_identity': 0,
            'nr_similarity': 0,
            'swissprot_evalue': 1e-5,
            'swissprot_identity': 0,
            'swissprot_similarity': 0,
            'cog_evalue': 1e-5,
            'cog_identity': 0,
            'cog_similarity': 0,
            'kegg_evalue': 1e-5,
            'kegg_identity': 0,
            'kegg_similarity': 0,
            'pfam_evalue': 1e-5
        }
        taxonomy = 'Fungi'
        exp_level = 'transcript'
        wf.test_api.run(annot_merge_output_dir, params_dict, taxonomy, exp_level)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
