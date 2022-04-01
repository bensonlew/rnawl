# -*- coding: utf-8 -*-
from biocluster.config import Config
from biocluster.api.database.base import Base
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from collections import defaultdict
from mbio.api.database.prok_rna.api_base import ApiBase
from bson.objectid import ObjectId
import types


class CreatTableRelation(ProkRNAController):
    def __init__(self, project_type='prok_rna', creat_targets=None):
        super(CreatTableRelation, self).__init__()
        self._project_type = project_type
        self.creat_targets = creat_targets
        # self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

    def creat_table(self):
        target=self.get_target()
        main_id=self.get_main()
        # for each in target:
        #     main_table, detail_table, detail_table_key = each
        #     up_dict=defaultdict(list)
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
        self.update_db_record('sg_table_relation', main_id, target=target)

    def get_main(self):
        mongo_data = [
            ('name', 'table_relation'),
            ('task_id', 'prok_rna'),
        ]
        main_table_id = self.prok_rna.insert_main_table('sg_table_relation', mongo_data)
        return main_table_id

    def get_target(self):
        target = \
            [(u'sg_geneset', u'sg_geneset_detail', u'geneset_id'),
            (u'sg_specimen_rfam', None, None),
              (u'sg_geneset_kegg_enrich', [u'sg_geneset_kegg_enrich_detail', u'sg_geneset_kegg_enrich_gene'],
               u'kegg_enrich_id'),
              (u'sg_specimen_group', None, None),
              (u'sg_srna_fold', u'sg_srna_fold_detail', u'fold_id'),
              (u'sg_structure_promote', u'sg_structure_promote_detail', u'promote_id'),
              (u'sg_query_seq', None, None),
              (u'sg_exp', u'sg_exp_detail', u'exp_id'),
              (u'sg_geneset_go_class', [u'sg_geneset_go_class_gene', u'sg_geneset_go_class_detail'], u'go_regulate_id'),
              (u'sg_geneset_kegg_class',
               [u'sg_geneset_kegg_class_statistic', u'sg_geneset_kegg_class_detail', u'sg_geneset_kegg_class_pic',
                u'sg_geneset_kegg_class_gene'], u'kegg_id'),
              (u'sg_assessment_chrom_distribution', u'sg_assessment_chrom_distribution_detail',
               u'chrom_distribution_id'),
              (
              u'sg_geneset_go_enrich', [u'sg_geneset_go_enrich_detail', u'sg_geneset_go_enrich_gene'], u'go_enrich_id'),
              (u'sg_geneset_ipath', u'sg_geneset_ipath_detail', u'ipath_id'),
              (u'sg_srna_predict', [u'sg_srna_predict_detail', u'sg_srna_predict_bed'], u'predict_id'),
              (u'sg_annotation_stat', u'sg_annotation_stat_detail', u'stat_id'),
              (u'sg_geneset_circ', [u'sg_geneset_circ_graph', u'sg_geneset_circ_detail'], u'circ_id'),
              (u'sg_exp_corr', u'sg_exp_corr_detail', u'corr_id'),
              (u'sg_specimen_group_compare', None, None),
              (u'sg_result_table_deposit', None, None),
              (u'sg_status', None, None),
              (u'sg_srna_anno',
               [u'sg_srna_anno_rfam', u'sg_srna_anno_venn', u'sg_srna_anno_stat', u'sg_srna_anno_detail'], u'annot_id'),
              (u'sg_srna_target', [u'sg_srna_target_stat', u'sg_srna_target_hybrid', u'sg_srna_target_plex'],
               u'target_id'),
              (u'sg_annotation_go', [u'sg_annotation_go_list', u'sg_annotation_go_graph', u'sg_annotation_go_stat',
                                     u'sg_annotation_go_level'], u'go_id'),
              (u'sg_specimen_mapping', None, None),
              (u'sg_assessment_coverage', u'sg_assessment_coverage_detail', u'coverage_id'),
              (u'sg_geneset_cog_class', [u'sg_geneset_cog_class_detail', u'sg_geneset_cog_class_gene'],
               u'geneset_cog_id'),
              (u'sg_specimen', u'sg_specimen_graphic', u'specimen_id'),
              (u'sg_assessment_saturation', u'sg_assessment_saturation_curve', u'saturation_id'),
              (u'sg_assessment_distribution', u'sg_assessment_distribution_detail', u'distribution_id'),
              (u'sg_exp_graph', [u'sg_exp_graph_box', u'sg_exp_graph_density', u'sg_exp_graph_volin'], u'graph_id'),
              (u'sg_diff', [u'sg_diff_volcano', u'sg_diff_summary', u'sg_diff_detail', u'sg_diff_scatter'], u'diff_id'),
              (u'sg_species_information', None, None), (u'sg_task', None, None),
              (u'sg_annotation_kegg',
               [u'sg_annotation_kegg_pic', u'sg_annotation_kegg_categories', u'sg_annotation_kegg_table',
                u'sg_annotation_kegg_level'], u'kegg_id'),
              (u'sg_geneset_venn', None, None),
              (u'sg_annotation_cog',
               [u'sg_annotation_cog_stat', u'sg_annotation_cog_summary', u'sg_annotation_cog_detail'], u'cog_id'),
              (u'sg_structure_utr', [u'sg_structure_utr_detail', u'sg_structure_utr_len'], u'utr_id'),
              (u'sg_annotation_query', u'sg_annotation_query_detail', u'query_id'),
              (u'sg_structure_operon',
               [u'sg_structure_operon_detail', u'sg_structure_operon_len', u'sg_structure_operon_num'], u'operon_id'),
              (u'sg_annotation_swissprot', u'sg_annotation_swissprot_detail', u'swissprot_id'),
              (u'sg_exp_pca', u'sg_exp_pca_detail', u'pca_id'),
              (u'sg_geneset_cluster', [u'sg_geneset_cluster_tree', u'sg_geneset_cluster_detail'], u'cluster_id'),
              (u'sg_structure_tsstts', u'sg_structure_tsstts_detail', u'tsstts_id'),
              (u'sg_exp_venn', u'sg_exp_venn_detail', u'venn_id'),
              (u'sg_annotation_pfam',
               [u'sg_annotation_pfam_bar', u'sg_annotation_pfam_summary', u'sg_annotation_pfam_detail'], u'pfam_id'),
              (u'sg_snp', [u'sg_snp_stat', u'sg_snp_detail'], u'snp_id'),
              (u'sg_annotation_nr', u'sg_annotation_nr_detail', u'nr_id')]

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
