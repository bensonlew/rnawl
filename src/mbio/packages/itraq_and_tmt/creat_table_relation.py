# -*- coding: utf-8 -*-
from biocluster.config import Config
from biocluster.api.database.base import Base
from mainapp.controllers.project.itraq_and_tmt_controller import ItraqTmtController
from collections import defaultdict
from mbio.api.database.itraq_and_tmt.api_base import ApiBase
from bson.objectid import ObjectId
import types


class CreatTableRelation(ItraqTmtController):
    def __init__(self, project_type='itraq_tmt', creat_targets=None):
        super(CreatTableRelation, self).__init__()
        self._project_type = project_type
        self.creat_targets = creat_targets
        # self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

    def creat_table(self):
        target=self.get_target()
        main_id=self.get_main()
        for each in target:
            main_table, detail_table, detail_table_key = each
            up_dict=defaultdict(list)
            # if detail_table:
            #     if not type(detail_table) == list:
            #         detail_table = [detail_table]
            #         detail_table_key = [detail_table_key]
            #     else:
            #         if not type(detail_table_key) == list:
            #             detail_table_key = [detail_table_key]*len(detail_table)
            #     for table, table_key in zip(detail_table, detail_table_key):
            #         if not table_key:
            #             raise Exception('you should specify detail_table_key whose value is main table "_id"')
            #         up_dict[main_table].append([table,table_key])
            # else:
            #     up_dict[main_table]=[]
            self.update_db_record('sg_table_relation_v2', main_id, target=target)

    def get_main(self):
        mongo_data = [
            ('name', 'table_relation'),
            ('task_id', 'itraq_tmt'),
        ]
        main_table_id = self.itraq_tmt.insert_main_table('sg_table_relation_v2', mongo_data)
        return main_table_id

    def get_target(self):
        target = \
            [(u'sg_specimen', None, None),  # yes
             (u'sg_peptide_error', None, None),  # yes
             (u'sg_peptide_len', None, None),  # yes
             (u'sg_peptide_num', None, None),  # yes
             (u'sg_protein_coverage', None, None),  # yes
             (u'sg_protein_info', None, None),  # yes
             (u'sg_protein_mw', None, None),  # yes
             (u'sg_protein_mw', None, None),  # yes
             (u'sg_proteinset_info', None, None),  # yes
             (u'sg_query_seq', None, None),  # yes
             (u'sg_software', None, None),  # yes
             (u'sg_annotation_stat', u'sg_annotation_stat_detail', u'stat_id'),  # yes
             (u'sg_annotation_pfam', [u'sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], u'pfam_id'),  # yes
             (u'sg_annotation_query', u'sg_annotation_query_detail', u'query_id'),  # yes
             (u'sg_annotation_go',
              [u'sg_annotation_go_detail', 'sg_annotation_go_graph', 'sg_annotation_go_level', 'sg_annotation_go_list'],
              u'go_id'),
             (u'sg_annotation_cog', u'sg_annotation_cog_detail', u'cog_id'),  # yes
             (u'sg_annotation_subloc', [u'sg_annotation__subloc_bar', u'sg_annotation__subloc_detail'], u'subloc_id'),
             # yes
             (u'sg_annotation_kegg',
              [u'sg_annotation_kegg_categories', 'sg_annotation_kegg_level', 'sg_annotation_kegg_table',
               u'sg_annotation_kegg_important'], u'kegg_id'),  # yes
             (u'sg_proteinset_venn', None, None),  # yes
             (u'sg_proteinset_kegg_class', [u'sg_proteinset_kegg_class_detail', u'sg_proteinset_kegg_class_pathway',
                                            u'sg_proteinset_kegg_class_statistic'], u'kegg_id'),  # yes
             (u'sg_specimen_group', None, None),  # yes
             (u'sg_proteinset', u'sg_proteinset_detail', u'proteinset_id'),  # yes
             (u'sg_search_db', u'sg_search_db_detail', u'searchdb_id'),  # yes
             (u'sg_proteinset_ipath', u'sg_proteinset_ipath_detail', u'ipath_id'),  # yes
             (u'sg_exp_venn', u'sg_exp_venn_detail', u'venn_id'),  # yes
             (u'sg_express_corr_', u'sg_exp_corr_detail', u'corr_id'),  # yes
             (u'sg_proteinset_cog_class', u'sg_proteinset_cog_class_detail', u'proteinset_cog_id'),  # yes
             (u'sg_status', None, None),  # yes
             (u'sg_specimen_group_compare', None, None),  # yes
             (
             u'sg_proteinset_cluster', [u'sg_proteinset_cluster_detail', u'sg_proteinset_cluster_tree'], u'cluster_id'),
             # yes
             (u'sg_proteinset_circ', [u'sg_proteinset_circ_detail', u'sg_proteinset_circ_graph'], u'circ_id'),  # yes
             (u'sg_task', None, None),  # yes
             (u'sg_proteinset_kegg_enrich', u'sg_proteinset_kegg_enrich_detail', u'kegg_enrich_id'),  # yes
             (u'sg_proteinset_go_enrich', u'sg_proteinset_go_enrich_detail', u'go_enrich_id'),  # yes
             (u'sg_diff', [u'sg_diff_detail', u'sg_diff_summary'], u'diff_id'),  # yes
             (u'sg_express', u'sg_express_detail', u'express_id'),  # yes
             (u'sg_result_table_deposit', None, None),  # yes
             (u'sg_proteinset_go_class', u'sg_proteinset_go_class_detail', u'go_regulate_id'),  # yes
             (u'sg_proteinset_go_class2', u'sg_proteinset_go_class2_detail', u'go_regulate_id'),  # yes
             (u'sg_express_pca', [u'sg_express_pca_gene_detail', u'sg_express_pca_sample_detail'], u'pca_id'), ]  # yes
        return target

    def update_db_record(self, table_name, record_id=None, query_dict=None, insert_dict=None, **kwargs):
        if record_id is not None:
            if isinstance(record_id, types.StringTypes):
                record_id = ObjectId(record_id)
            elif isinstance(record_id, ObjectId):
               record_id = record_id
            else:
                raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = self.db[table_name]
        if query_dict:
            if record_id is not None:
                query_dict.update({"_id": record_id})
        else:
            if record_id is not None:
                query_dict = {"_id": record_id}
            else:
                raise Exception("Please provide query dict")
        if insert_dict:
            kwargs.update(insert_dict)
        conn.update(query_dict, {"$set": kwargs}, upsert=True)

if __name__ == '__main__':
    CreatTableRelation().creat_table()