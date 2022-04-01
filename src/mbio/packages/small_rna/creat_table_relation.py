# -*- coding: utf-8 -*-
from biocluster.config import Config
from biocluster.api.database.base import Base
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from collections import defaultdict
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from bson.objectid import ObjectId
import types


class CreatTableRelation(RefRnaV2Controller):
    def __init__(self, project_type='ref_rna_v2', creat_targets=None):
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
            ('task_id', 'ref_rna_v2'),
        ]
        main_table_id = self.ref_rna_v2.insert_main_table('sg_table_relation_v2', mongo_data)
        return main_table_id

    def get_target(self):
        target = \
            [(u'sg_specimen', None, None),  # yes
             (u'sg_assessment_chrom_distribution', 'sg_assessment_chrom_distribution_detail', u'chrom_distribution_id'),
             # yes
             (u'sg_assessment_distribution', 'sg_assessment_distribution_detail', u'distribution_id'),  # yes
             (u'sg_assessment_coverage', 'sg_assessment_coverage_detail', u'coverage_id'),  # yes
             (u'sg_assessment_saturation', 'sg_assessment_saturation_curve', u'saturation_id'),  # yes
             (u'sg_annotation_stat', u'sg_annotation_stat_detail', u'stat_id'),  # yes
             (u'sg_annotation_blast', u'sg_annotation_blast_detail', u'blast_id'),  # yes
             (u'sg_annotation_nr', u'sg_annotation_nr_detail', u'nr_id'),  # yes
             (u'sg_annotation_swissprot', 'sg_annotation_swissprot_pie', u'swissprot_id'),  # yes
             (u'sg_annotation_pfam', [u'sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], u'pfam_id'),  # yes
             (u'sg_annotation_query', u'sg_annotation_query_detail', u'query_id'),  # yes
             (u'sg_geneset_cog', u'sg_annotation_cog_detail', u'cog_id'),  # yes
             (u'sg_annotation_go',
              [u'sg_annotation_go_detail', 'sg_annotation_go_graph', 'sg_annotation_go_level', 'sg_annotation_go_list'],
              u'go_id'),  # yes
             (u'sg_annotation_kegg',
              [u'sg_annotation_kegg_categories', 'sg_annotation_kegg_level', 'sg_annotation_kegg_table'], u'kegg_id'),
             (u'sg_exp', u'sg_exp_detail', u'express_id'),  # yes
             (u'sg_exp_corr', u'sg_exp_corr_detail', u'corr_id'),  # yes
             (u'sg_exp_corrsf', u'sg_exp_corrsf_detail', u'expcorr_id'),  # yes
             (u'sg_exp_graph', [u'sg_exp_graph_box', u'sg_exp_graph_density', u'sg_exp_graph_volin'], u'graph_id'),
             (u'sg_diff', [u'sg_diff_detail', u'sg_diff_summary', u'sg_diff_scatter', u'sg_diff_volcano'], u'diff_id'),
             (u'sg_exp_pca', u'sg_exp_pca_detail', u'pca_id'),  # yes
             (u'sg_exp_venn', u'sg_exp_venn_detail', u'venn_id'),  # yes
             (u'sg_geneset', u'sg_geneset_detail', u'geneset_id'),  # yes
             (u'sg_geneset_circ', [u'sg_geneset_circ_detail', u'sg_geneset_circ_graph'], u'circ_id'),  # yes
             (u'sg_geneset_cluster', [u'sg_geneset_cluster_detail', u'sg_geneset_cluster_tree'], u'cluster_id'),  # yes
             (u'sg_geneset_cog_class', u'sg_geneset_cog_class_detail', u'geneset_cog_id'),
             (u'sg_geneset_go_class', u'sg_geneset_go_class_detail', u'go_regulate_id'),  # yes
             (u'sg_geneset_go_enrich', u'sg_geneset_go_enrich_detail', u'go_enrich_id'),  # yes
             (u'sg_geneset_info', None, None),
             (u'sg_geneset_ipath', u'sg_geneset_ipath_detail', u'ipath_id'),  # yes
             (u'sg_geneset_kegg_class',
              [u'sg_geneset_kegg_class_detail', u'sg_geneset_kegg_class_pathway', u'sg_geneset_kegg_class_statistic'],
              u'kegg_id'),  # yes
             (u'sg_geneset_kegg_enrich', u'sg_geneset_kegg_enrich_detail', u'kegg_enrich_id'),  # yes
             (u'sg_geneset_ppi',
              [u'sg_geneset_ppi_edge', u'sg_geneset_ppi_node', u'sg_geneset_ppi_nodes_centrality',
               u'sg_geneset_ppi_nodes_degreestat', u'sg_geneset_ppi_stat', ],
              u'ppi_id'),  # yes
             (u'sg_geneset_venn', None, None),  # yes
             (u'sg_result_table_deposit', None, None),  # yes
             (u'sg_rmats_diffcomp', u'sg_rmats_diffcomp_detail', u'rmats_diffcomp_id'),  # yes
             (u'sg_query_seq', None, None),
             (u'sg_software_para', None, None),
             (u'sg_snp', [u'sg_snp_detail', u'sg_snp_stat'], u'snp_id'),  # yes
             (u'sg_species_information', u'sg_species_information_detail', u'species_id'),  # yes
             (u'sg_specimen', u'sg_specimen_graphic', u'specimen_id'),  # yes
             (u'sg_specimen_group', u'sg_specimen_group_compare', u'specimen_group_id'),  # yes
             (u'sg_specimen_mapping', None, None),  # yes
             (u'sg_splicing_rmats', u'sg_splicing_rmats_detail', u'splicing_id'),  # yes
             (u'sg_splicing_rmats_stats',
              [u'sg_splicing_rmats_psi', u'sg_splicing_rmats_model', u'sg_splicing_rmats_graph'], u'stat_id'),  # yes
             (u'sg_status', None, None),  # yes
             (u'sg_task', None, None),  # yes
             (u'sg_tf_predict', u'sg_tf_predict_detail', u'tf_predict_id'),  # yes
             (u'sg_tf_stat', [u'sg_tf_stat_bar_detail', u'sg_tf_stat_circos_detail'], u'tf_stat_id'),  # yes
             (u'sg_transcripts', [u'sg_transcripts_relations', u'sg_transcripts_seq_type', u'sg_transcripts_step'],
              u'transcripts_id'),  # yes
             (u'sg_wgcna_module',
              [u'sg_wgcna_module_cluster_detail', u'sg_wgcna_module_corr_detail', u'sg_wgcna_module_eigengenes_detail',
               u'sg_wgcna_module_membership_detail', u'sg_wgcna_module_stat_detail', u'sg_wgcna_module_tree_detail',
               u'sg_wgcna_module_unmerged_detail'], u'module_id'),  # yes
             (u'sg_wgcna_network', [u'sg_wgcna_network_edge_detail', u'sg_wgcna_network_node_detail'], u'network_id'),
             (u'sg_wgcna_pipeline', None, None),
             (u'sg_wgcna_prepare', u'sg_wgcna_prepare_detail', u'prepare_id'),
             (u'sg_wgcna_relate', [u'sg_wgcna_relate_gene_detail', u'sg_wgcna_relate_module_detail'], u'relate_id')
             ]  # yes
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