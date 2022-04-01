# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

from __future__ import print_function

import time

from biocluster.config import Config
from concurrent.futures import ThreadPoolExecutor

medical_add ={"stop_record":"sg_cancel",
              "diffgene":"sg_diff",
              "diffgene_batch_analysis":"sg_diff_geneset_pipline",
              "diffgenevenn":"sg_diff_geneset_venn",
              "diffgenecluster":"sg_diff_geneset_cluster",
              "genesetcorrelation":"sg_exp_corrsf",
              "snpsearch":"sg_snp_search",
              'splicingrmats_stat': 'sg_splicing_rmats_stats',
              'asprofile': 'sg_asprofile',
              'asprofilediff':'sg_asprofile_diff',
              'asprofilesearch':'sg_asprofile_search',
              'genefusion':"sg_gene_fusion",
              'genefusionvenn':"sg_gene_fusion_venn",
              'genesetreactome':"sg_geneset_reactome_class",
              'genesetdo':"sg_geneset_do_class",
              'genesetreactome_rich':"sg_geneset_reactome_enrich",
              'genesetdo_rich':'sg_geneset_do_enrich',
              'genesetdisgenet_rich':'sg_geneset_disgenet_enrich',
              'genesetgsea': 'sg_geneset_gsea',
              'genesetgsva': 'sg_geneset_gsva',
              'diffgenecirc':'sg_diff_geneset_circ',
              'diffgenedisgenet_rich':'sg_diff_geneset_idsgenet_enrich',
              'seq_download':'sg_seq_extract',
              "geneset_seq_download" :"sg_seq_extract"

              }
SUB2NAME = {'cerna': 'cerna',
            'cerna_sankey': "cerna_sankey",
            u'diff_circrna': u'diff',
            u'diff_lncrna': u'diff',
            u'diff_mirna': u'diff',
            u'diff_mrna': u'diff',
            u'expcorr': u'exp_corr',
            u'expgraph_circrna': u'exp_graph',
            u'expgraph_lncrna': u'exp_graph',
            u'expgraph_mirna': u'exp_graph',
            u'expgraph_mrna': u'exp_graph',
            u'exppca': u'exp_pca',
            u'expvenn_circrna': u'exp_venn',
            u'expvenn_lncrna': u'exp_venn',
            u'expvenn_mirna': u'exp_venn',
            u'expvenn_mrna': u'exp_venn',
            u'genesetcirc': u'geneset_circ',
            u'genesetcluster': u'geneset_cluster',
            u'genesetcog': u'geneset_cog_class',
            u'genesetcorrsf': u'exp_corrsf',
            u'genesetgo': u'geneset_go_class',
            u'genesetgo_acyclic': u'geneset_go_dag',
            u'genesetgo_rich': u'geneset_go_enrich',
            'genesetgsea': 'geneset_gsea',
            'genesetipath': 'geneset_ipath',
            u'genesetkegg': u'geneset_kegg_class',
            u'genesetkegg_rich': u'geneset_kegg_enrich',
            u'genesetvenn': u'geneset_venn',
            'masigpro': 'masigpro',
            'ppinetwork': 'geneset_ppi',
            u'snp': u'snp',
            'splicingrmats': 'splicing_rmats',
            'splicingrmats_diffcomp': 'splicing_rmats_diffcomp',
            'splicingrmats_model': 'splicing_rmats_model',
            'splicingrmats_stat': 'splicing_rmats_stats',
            'stem': 'stem',
            u'targe': u'target_mi',
            u'targe_cistrans': u'target_cistrans',
            u'target_mi': u'target_mi',
            'tfpredict': 'tf_predict',
            'wgcnapipeline': 'wgcna_pipeline',
            'wgcnaprepare': 'wgcna_prepare',
            'wgcnamodule': 'wgcna_module',
            'wgcnarelate': 'wgcna_relate',
            'wgcnanetwork': 'wgcna_network'}


class DeleteRecordsMongo(object):
    def __init__(self, task_id, project_type, submit_location, status='failed'):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self.db_static = Config().get_mongo_client(mtype=project_type, dydb_forbid=True)[
            Config().get_mongo_dbname(project_type, dydb_forbid=True)]
        self.task_id = task_id
        self.status = status
        self._project_type = project_type
        self.submit_location = submit_location
        self.delete_set = self.get_delete_targets()
        self.sub2name = {
            'expcorr': 'sg_exp_corr',
            'exppca': 'sg_exp_pca',
            'wgcnarelate': 'sg_wgcna_relate',
            'genesetgo_acyclic': 'sg_geneset_go_dag',
            'ppinetwork': 'sg_geneset_ppi',
            'snp_detail': 'sg_snp',
            'splicingrmats_type': 'sg_splicing_rmats_stats',
            'genesetgo_rich': 'sg_geneset_go_enrich',
            'wgcnamodule': 'sg_wgcna_module',
            'expgraph': 'sg_exp_graph',
            'splicingrmats_diffcomp': 'sg_splicing_rmats_diffcomp',
            'genesetcorrsf': 'sg_exp_corrsf',
            'wgcnanetwork': 'sg_wgcna_network',
            'splicingrmats_model': 'sg_splicing_rmats_model',
            'genesetipath': 'sg_geneset_ipath',
            'expvenn': 'sg_exp_venn',
            'genesetcog': 'sg_geneset_cog_class',
            'genesetkegg_rich': 'sg_geneset_kegg_enrich',
            'splicingrmats': 'sg_splicing_rmats',
            'wgcnaprepare': 'sg_wgcna_prepare',
            'tfpredict': 'sg_tf_predict',
            'wgcnapipeline': 'sg_wgcna_pipeline',
            'genesetgo': 'sg_geneset_go_class',
            'diff_detail': 'sg_diff',
            'exp_detail': 'sg_exp',
            'blast': 'sg_blast',
            'snp': 'sg_snp',
            'ssr': 'sg_ssr',
            'genesetvenn': 'sg_geneset_venn',
            'tfbspredict': 'sg_tfbs_predict',
            'genesetkegg': 'sg_geneset_kegg_class',
            'genesetcluster': 'sg_geneset_cluster',
            'genesetcirc': 'sg_geneset_circ',
            'annotationstat': 'sg_annotation_stat',
            'geneset_upload': 'sg_geneset',
            'splicingrmats_stat': 'sg_splicing_rmats_stats',
            'minet': 'sg_minet',
            'targe': 'sg_target',
            'specimen_mapping': 'sg_mapping',
            'mapping_circs': 'sg_circos',
            'lncrna_family': 'sg_lncrna_family',
            'loccon': 'sg_lncrna_loc_cons',
            'ortholog': 'sg_lncrna_ortholog',
            'lncrna_seq_cons': 'sg_lncrna_seq_cons',
            'mirna_precursor': 'sg_mirna_precursor',
            'target_cis': 'sg_target_cis',
            'trans': 'sg_target_trans',
            'cerna': 'sg_cerna',
            'masigpro': 'sg_masigpro',
            'stem': 'sg_stem',
            'genesetgsea': 'sg_geneset_gsea',
            'genesetgsva': 'sg_geneset_gsva',
            'asprofile_detail': 'sg_asprofile',
            'asprofilediff':'sg_asprofile_diff',
            'expbatch': 'sg_exp_batch'
        }
        if project_type == 'whole_transcriptome':
            self.sub2name.update(SUB2NAME)
        if project_type == 'medical_transcriptome':
            self.sub2name.update(medical_add)
        self.delete_id = list()
        self.delete_targets = self.get_all_delete_target()

    def run(self):
        self.deal_sg_status()
        target_collections = self.delete_targets
        if not target_collections:
            self.println('The deletion target is empty, and the operation ends')
            return
        with ThreadPoolExecutor(8) as pool:
            pool.map(self.remove_table_by_task_id, target_collections)

    def deal_sg_status(self):
        conn = self.db['sg_status']
        if self.submit_location != 'all':
            a = conn.find({'task_id': self.task_id, 'status': self.status, 'submit_location': self.submit_location},
                          {'_id': 1})
        else:
            a = conn.find({'task_id': self.task_id, 'status': self.status},
                          {'_id': 1})
        if type(a) == dict:
            a = [a]
        for each in a:
            conn.update({'_id': each['_id']}, {'$set': {'status': 'deleted'}})

    def get_some_delete_target(self, submit_location):
        delete_targets = []
        target_collection = ''
        if self.submit_location in self.sub2name:
            target_collection = self.sub2name[self.submit_location]
        for i in self.delete_set:
            if target_collection == i[0]:
                delete_targets.append(i)
                break
        if not delete_targets:
            raise Exception('There is no corresponding table for the incoming submit location')
        return delete_targets

    def get_all_delete_target(self):
        delete_targets = []
        conn = self.db['sg_status']
        if self.submit_location.lower() == 'all':
            cursor = conn.find({'task_id': self.task_id, 'status': self.status},
                               {'submit_location': 1, 'table_id': 1})
        else:
            cursor = conn.find(
                {'task_id': self.task_id, 'status': self.status, 'submit_location': self.submit_location},
                {'submit_location': 1, 'table_id': 1})
        if type(cursor) == dict:
            cursor = [cursor]
        for document in cursor:
            self.delete_id.append(document['table_id'])
            if document['submit_location'] in self.sub2name:
                target_collection = self.sub2name[document['submit_location']]
                for collection in self.delete_set:
                    if target_collection == collection[0]:
                        delete_targets.append(collection)
        self.println(delete_targets)
        self.println(self.delete_id)
        if not delete_targets:
            raise Exception('There is no corresponding table for the incoming submit location')
        return delete_targets

    def get_delete_targets(self):
        find_result = self.db_static['sg_table_relation'].find_one({})
        if find_result:
            target = find_result['target']
        else:
            raise Exception(
                'The collection <sg_table_relation> does not exist for the {} project'.format(self._project_type))
        return target

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        conn = self.db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        result = conn.find_one(kwargs)
        if result:
            conn.delete_many(kwargs)
            self.println('Success to delete records in {} by query {}'.format(table_name, kwargs))
        else:
            self.println('No record to delete from {} by query {}'.format(table_name, kwargs))

    def remove_table_by_task_id(self, tuple_arg):
        main_table, detail_table, detail_table_key = tuple_arg
        all_collections = self.db.collection_names()
        if main_table not in list(all_collections):
            self.println('Warning: {} was not found in {}'.format(main_table, self._project_type))
            return
        conn = self.db[main_table]
        cursor = conn.find({'task_id': self.task_id, 'status': self.status}, {'_id': 1})
        if type(cursor) == dict:
            cursor = [cursor]
        remove_num = 0
        for each in cursor:
            if each['_id'] in self.delete_id:
                remove_num += 1
                self.remove_db_record(main_table, _id=each['_id'])
                if detail_table:
                    if not type(detail_table) == list:
                        detail_table = [detail_table]
                        detail_table_key = [detail_table_key]
                    else:
                        if not type(detail_table_key) == list:
                            detail_table_key = [detail_table_key] * len(detail_table)
                    for table, table_key in zip(detail_table, detail_table_key):
                        if not table_key:
                            raise Exception('You should specify schedule key whose value is main table <_id>')
                        self.remove_db_record(table, query_dict={table_key: each['_id']})
        if remove_num:
            self.println('Found {} main table(s) in {} by task_id {}.'.format(remove_num, main_table, self.task_id))
            self.println('And, finished to remove records of {} and {} by task_id {}'.format(main_table, detail_table,
                                                                                             self.task_id))

    def println(self, value):
        if hasattr(self, 'logger'):
            self.logger.info(value)
        else:
            print(value)


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 5:
        exit('Usage: python delete_demo.py <task_id> <project_type> <submit_location> <status>')
    submit_location = sys.argv[3]
    project_type = sys.argv[2]
    task_id = sys.argv[1]
    status = sys.argv[4]
    start_time = time.time()
    del_demo = DeleteRecordsMongo(task_id, project_type, submit_location, status)
    del_demo.run()
    end_time = time.time()
    print('total time: {}'.format(end_time - start_time))
