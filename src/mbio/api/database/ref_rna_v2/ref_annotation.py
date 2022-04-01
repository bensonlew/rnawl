# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import json
import unittest
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from mbio.api.database.ref_rna_v2.api_base import ApiBase
import random


class RefAnnotation(ApiBase):
    def __init__(self, bind_object):
        super(RefAnnotation, self).__init__(bind_object)
        self.result_dir = ''
        self.result_file = {}
        self.trans_gene = {}
        self.trans_isgene = {}
        self.task_id = self.bind_object.sheet.id
        self.has_new = True
        self.anno_type = 'origin'
        self.species_name = ""
        self.version = "v2"
        self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
        self.new_gene_set = set()
        self.new_trans_set = set()
        self.known_gene_set = set()
        self.known_trans_set = set()

        #self._db_name = Config().MONGODB + '_ref_rna'

    def set_result_dir(self, merge_annot_path):
        '''
        根据注释模块结果导入结果路径
        '''
        merge_annot_path =  merge_annot_path + "/"
        self.result_dir = merge_annot_path
        self.bind_object.logger.info("导入**** {}".format(merge_annot_path))
        self.result_file['new_stat_path'] = os.path.join(merge_annot_path, "newannot_class/all_stat.xls")
        self.result_file['new_venn_path'] = os.path.join(merge_annot_path, "newannot_class")
        self.result_file['ref_stat_path'] = os.path.join(merge_annot_path, "refannot_class/all_stat.xls")
        self.result_file['ref_venn_path'] = os.path.join(merge_annot_path, "refannot_class/")
        self.result_file['all_stat_path'] = os.path.join(merge_annot_path, "allannot_class/all_stat.xls")
        self.result_file['all_venn_path'] = os.path.join(merge_annot_path, "allannot_class/")
        # blast_path = annotation_mudule_dir + "/anno_stat/blast"
        # blast_ref_path = annotation_mudule_dir + "/anno_stat/blast"
        # for db in ["nr", "swissprot", "string", "kegg"]:
        for db in ["nr", "swissprot"]:
            db_gene_name = db + "genes_blast_path"
            db_trans_name = db + "trans_blast_path"
            self.result_file[db_gene_name] = merge_annot_path + "newannot_class/" + db + "/{}_blast_gene.xls".format(db)
            self.result_file[db_trans_name] = merge_annot_path + "newannot_class/" + db + "/{}_blast_tran.xls".format(db)
            self.result_file[db_gene_name + "ref"] = merge_annot_path + "refannot_class/" + db + "/{}_blast_gene.xls".format(db)
            self.result_file[db_trans_name + "ref"] = merge_annot_path + "refannot_class/" + db + "/{}_blast_tran.xls".format(db)

        # nr_path = annotation_mudule_dir + "/anno_stat/blast_nr_statistics"

        self.result_file['gene_pfam_path'] = merge_annot_path + "newannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path'] = merge_annot_path + "newannot_class/pfam/pfam_domain_tran.xls"
        self.result_file['query_path'] = merge_annot_path + "newannot_class/all_annot.xls"
        # self.result_file['gene_query_path'] = merge_annot_path + "/anno_stat/gene_anno_detail.xls"

        self.result_file['gene_pfam_path' + "ref"] = merge_annot_path + "refannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path' + "ref"] = merge_annot_path + "refannot_class/pfam/pfam_domain_tran.xls"
        self.result_file['query_path' + "ref"] = merge_annot_path + "refannot_class/all_annot.xls"

        self.result_file['gene_pfam_path' + "all"] = merge_annot_path + "allannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path' + "all"] = merge_annot_path + "allannot_class/pfam/pfam_domain_tran.xls"
        # self.result_file['query_path' + "ref"] = merge_annot_path + "refannot_class/all_annot.xls"


        self.result_file['n_sum_path'] = merge_annot_path + "newannot_class/cog/cog_summary_tran.xls"
        self.result_file['n_gene_sum_path'] = merge_annot_path + "newannot_class/cog/cog_summary_gene.xls"

        self.result_file['n_sum_path' + 'ref'] = merge_annot_path + "refannot_class/cog/cog_summary_tran.xls"
        self.result_file['n_gene_sum_path' + 'ref'] = merge_annot_path + "refannot_class/cog/cog_summary_gene.xls"

        # go_path = merge_annot_path + "/go"
        # gene_go_path = merge_annot_path + "/anno_stat/go_stat"

        self.result_file['stat_level2'] = merge_annot_path + "newannot_class/go/go_lev2_tran.stat.xls"
        self.result_file['stat_level3'] = merge_annot_path + "newannot_class/go/go_lev3_tran.stat.xls"
        self.result_file['stat_level4'] = merge_annot_path + "newannot_class/go/go_lev4_tran.stat.xls"
        self.result_file['gene_stat_level2'] = merge_annot_path + "newannot_class/go/go_lev2_gene.stat.xls"
        self.result_file['gene_stat_level3'] = merge_annot_path + "newannot_class/go/go_lev3_gene.stat.xls"
        self.result_file['gene_stat_level4'] = merge_annot_path + "newannot_class/go/go_lev4_gene.stat.xls"
        self.result_file['gos_path'] = merge_annot_path + "newannot_class/go/go_list_tran.xls"
        self.result_file['gene_gos_path'] = merge_annot_path + "newannot_class/go/go_list_gene.xls"
        '''
        self.result_file['n_stat_level2'] = merge_annot_path + "/go/go12level_statistics.xls"
        self.result_file['n_stat_level3'] = merge_annot_path + "/go/go123level_statistics.xls"
        self.result_file['n_stat_level4'] = merge_annot_path + "/go/go1234level_statistics.xls"

        self.result_file['n_gene_stat_level2'] = merge_annot_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        self.result_file['n_gene_stat_level3'] = merge_annot_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        self.result_file['n_gene_stat_level4'] = merge_annot_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        '''
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

        # self.result_file['layer_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_layer_tran.xls"
        self.result_file['pathway_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_tran.xls"
        self.result_file['png_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_tran_dir/"
        self.result_file['table_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_gene_tran.xls"
        # self.result_file['gene_layer_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_layer_gene.xls"
        self.result_file['gene_pathway_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_gene.xls"
        self.result_file['gene_png_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_gene_dir/"
        self.result_file['gene_table_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_gene_gene.xls"

        for key, value in self.result_file.items():
            if self.has_new == False:
                if key.endswith("_ref") or key.startswith("ref_"):
                    if os.path.exists(value):
                        pass
                    else:
                        self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="537018108")
            else:
                if os.path.exists(value):
                    pass
                else:
                    self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="537018109")
        self.bind_object.logger.info("数据路径正确，文件完整 {}")

    def set_result_dir_v3(self, merge_annot_path):
        '''
        根据注释模块结果导入结果路径
        '''
        merge_annot_path =  merge_annot_path + "/"
        self.result_dir = merge_annot_path
        self.bind_object.logger.info("导入**** {}".format(merge_annot_path))
        self.result_file['new_stat_path'] = os.path.join(merge_annot_path, "newannot_class/all_stat.xls")
        self.result_file['new_venn_path'] = os.path.join(merge_annot_path, "newannot_class")
        self.result_file['ref_stat_path'] = os.path.join(merge_annot_path, "refannot_class/all_stat.xls")
        self.result_file['ref_venn_path'] = os.path.join(merge_annot_path, "refannot_class/")
        self.result_file['all_stat_path'] = os.path.join(merge_annot_path, "allannot_class/all_stat.xls")
        self.result_file['all_venn_path'] = os.path.join(merge_annot_path, "allannot_class/")
        # blast_path = annotation_mudule_dir + "/anno_stat/blast"
        # blast_ref_path = annotation_mudule_dir + "/anno_stat/blast"
        # for db in ["nr", "swissprot", "string", "kegg"]:
        for db in ["nr", "swissprot"]:
            db_gene_name = db + "genes_blast_path"
            db_trans_name = db + "trans_blast_path"
            self.result_file[db_gene_name] = merge_annot_path + "newannot_class/" + db + "/{}_blast_gene.xls".format(db)
            self.result_file[db_trans_name] = merge_annot_path + "newannot_class/" + db + "/{}_blast_tran.xls".format(db)
            self.result_file[db_gene_name + "ref"] = merge_annot_path + "refannot_class/" + db + "/{}_blast_gene.xls".format(db)
            self.result_file[db_trans_name + "ref"] = merge_annot_path + "refannot_class/" + db + "/{}_blast_tran.xls".format(db)

        # nr_path = annotation_mudule_dir + "/anno_stat/blast_nr_statistics"

        self.result_file['gene_pfam_path'] = merge_annot_path + "newannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path'] = merge_annot_path + "newannot_class/pfam/pfam_domain_tran.xls"
        self.result_file['query_path'] = merge_annot_path + "newannot_class/all_annot.xls"
        # self.result_file['gene_query_path'] = merge_annot_path + "/anno_stat/gene_anno_detail.xls"

        self.result_file['gene_pfam_path' + "ref"] = merge_annot_path + "refannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path' + "ref"] = merge_annot_path + "refannot_class/pfam/pfam_domain_tran.xls"
        self.result_file['query_path' + "ref"] = merge_annot_path + "refannot_class/all_annot.xls"

        self.result_file['gene_pfam_path' + "all"] = merge_annot_path + "allannot_class/pfam/pfam_domain_gene.xls"
        self.result_file['pfam_path' + "all"] = merge_annot_path + "allannot_class/pfam/pfam_domain_tran.xls"
        # self.result_file['query_path' + "ref"] = merge_annot_path + "refannot_class/all_annot.xls"


        self.result_file['n_sum_path'] = merge_annot_path + "newannot_class/cog/summary.T.tsv"
        self.result_file['n_gene_sum_path'] = merge_annot_path + "newannot_class/cog/summary.G.tsv"

        self.result_file['n_sum_path' + 'ref'] = merge_annot_path + "refannot_class/cog/summary.T.tsv"
        self.result_file['n_gene_sum_path' + 'ref'] = merge_annot_path + "refannot_class/cog/summary.G.tsv"

        # go_path = merge_annot_path + "/go"
        # gene_go_path = merge_annot_path + "/anno_stat/go_stat"

        self.result_file['stat_level2'] = merge_annot_path + "newannot_class/go/go_lev2_tran.stat.xls"
        self.result_file['stat_level3'] = merge_annot_path + "newannot_class/go/go_lev3_tran.stat.xls"
        self.result_file['stat_level4'] = merge_annot_path + "newannot_class/go/go_lev4_tran.stat.xls"
        self.result_file['gene_stat_level2'] = merge_annot_path + "newannot_class/go/go_lev2_gene.stat.xls"
        self.result_file['gene_stat_level3'] = merge_annot_path + "newannot_class/go/go_lev3_gene.stat.xls"
        self.result_file['gene_stat_level4'] = merge_annot_path + "newannot_class/go/go_lev4_gene.stat.xls"
        self.result_file['gos_path'] = merge_annot_path + "newannot_class/go/go_list_tran.xls"
        self.result_file['gene_gos_path'] = merge_annot_path + "newannot_class/go/go_list_gene.xls"
        '''
        self.result_file['n_stat_level2'] = merge_annot_path + "/go/go12level_statistics.xls"
        self.result_file['n_stat_level3'] = merge_annot_path + "/go/go123level_statistics.xls"
        self.result_file['n_stat_level4'] = merge_annot_path + "/go/go1234level_statistics.xls"

        self.result_file['n_gene_stat_level2'] = merge_annot_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        self.result_file['n_gene_stat_level3'] = merge_annot_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        self.result_file['n_gene_stat_level4'] = merge_annot_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        '''
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

        # self.result_file['layer_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_layer_tran.xls"
        self.result_file['pathway_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_tran.xls"
        self.result_file['png_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_tran_dir/"
        self.result_file['table_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_gene_tran.xls"
        # self.result_file['gene_layer_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_layer_gene.xls"
        self.result_file['gene_pathway_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_gene.xls"
        self.result_file['gene_png_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_pathway_gene_dir/"
        self.result_file['gene_table_path' + '_all'] = merge_annot_path + "allannot_class/kegg/kegg_gene_gene.xls"

        for key, value in self.result_file.items():
            if self.has_new == False:
                if key.endswith("_ref") or key.startswith("ref_"):
                    if os.path.exists(value):
                        pass
                    else:
                        self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="537018110")
            else:
                if os.path.exists(value):
                    pass
                else:
                    self.bind_object.set_error('%s对应的结果文件%s 不存在，请检查', variables=(key, value), code="537018111")
        self.bind_object.logger.info("数据路径正确，文件完整 {}")



    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('assemble_id必须为ObjectId对象或其对应的字符串！', code="537018112")
        return object_id

    def get_trans2gene(self, trans2gene, trans2gene_ref=None):
        if trans2gene and os.path.exists(trans2gene):
            pass
        elif trans2gene_ref and os.path.exists(trans2gene_ref):
            pass
        else:
            self.bind_object.set_error('转录本基因对应的结果文件%s不存在，请检查', variables=(trans2gene), code="537018113")
        self.bind_object.logger.info("读入基因转录本对应关系文件 {}".format(trans2gene))
        if trans2gene:
            with open(trans2gene, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.trans_gene[line[0]] = line[1]
                    self.new_trans_set.add(line[0])
                    if line[2] == "yes":
                        self.trans_isgene[line[0]] = True
                        self.new_gene_set.add(line[1])
                    else:
                        self.trans_isgene[line[0]] = False
        if trans2gene_ref:
            with open(trans2gene_ref, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    self.trans_gene[line[0]] = line[1]
                    self.known_trans_set.add(line[0])
                    if line[2] == "yes":
                        self.known_gene_set.add(line[1])
                        self.trans_isgene[line[0]] = True
                    else:
                        self.trans_isgene[line[0]] = False

        self.bind_object.logger.info("读入基因转录本对应关系文件结束")


    def run(self, annotation_mudule_dir, trans2gene, trans2gene_ref, params_dict, taxonomy='Animals', exp_level='transcript', version = "v2", gene_exp = None, trans_exp = None):
        """
        annotation_mudule_dir
        annotation_ref_dir
        merge_gen
        merge_tran
        new_anno_path: 新序列注释的结果文件夹
        pfam_path:转录本的pfam_domain
        merge_tran_output: 转录本的merge_annot tool输出结果路径
        merge_gene_output: 基因的merge_annot tool"
        """
        self.bind_object.logger.info("开始到表情数据路径为 {}".format(annotation_mudule_dir))
        if version == "v3":
            self.set_result_dir_v3(annotation_mudule_dir)
        else:
            self.set_result_dir(annotation_mudule_dir)
        self.get_trans2gene(trans2gene, trans2gene_ref)

        task_id = self.task_id
        self.version = version

        params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
        self.remove_table_by_main_record(main_table='sg_annotation_stat', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
        stat_id = self.add_annotation_stat(name=None, params=params, seq_type="new", database="nr,swissprot,pfam,cog,go,kegg", taxon=taxonomy, result_dir=self.result_dir, exp_level=exp_level)

        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['ref_stat_path'], venn_path=self.result_file['ref_venn_path'], seq_type = "ref", exp_level=exp_level, gene_exp = gene_exp, trans_exp = trans_exp, task_id=task_id)
        if self.has_new:
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['new_stat_path'], venn_path=self.result_file['new_venn_path'], seq_type = "new", exp_level=exp_level, gene_exp = gene_exp, trans_exp = trans_exp, task_id=task_id)
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['all_stat_path'], venn_path=self.result_file['all_venn_path'], seq_type = "all", exp_level=exp_level, gene_exp = gene_exp, trans_exp = trans_exp, task_id=task_id)
        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)

        # self.update_db_record('sg_assembly', obj_id, status="end", main_id=obj_id)
        # blast_id = self.add_annotation_blast(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        #for db in ["nr", "swissprot"]:
            #self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="T", database=db, blast_path=self.result_file[db + "trans_blast_path"])
            #self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="G", database=db, blast_path=self.result_file[db + "gene_blast_path"])

        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))
        '''
        self.remove_table_by_main_record(main_table='sg_annotation_nr', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_nr_detail', 'sg_annotation_nr_pie'], detail_table_key='nr_id')
        nr_id = self.add_annotation_nr(name=None, params=params_select_nr, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="new", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path"],exp_level=exp_level)
        self.add_annotation_blast_nr_detail(blast_id=nr_id, seq_type="ref", anno_type="T", database='nr', blast_path=self.result_file["nr" + "trans_blast_path" + "ref"], exp_level=exp_level)

        self.update_db_record('sg_annotation_nr', nr_id, status="end",  main_id=nr_id)


        self.remove_table_by_main_record(main_table='sg_annotation_swissprot', task_id=task_id, type=self.anno_type,  detail_table=['sg_annotation_swissprot_detail', 'sg_annotation_swissprot_pie'], detail_table_key='swissprot_id')
        params_select_swissprot = dict([(k,params_dict.get(k,None)) for k in ('swissprot_evalue', 'swissprot_similarity', 'swissprot_identity')])
        params_select_swissprot = json.dumps(params_select_swissprot, sort_keys=True, separators=(',', ':'))
        swissprot_id = self.add_annotation_swissprot(name=None, params=params_select_swissprot, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="new", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path"], exp_level=exp_level)
        self.add_annotation_blast_swissprot_detail(blast_id=swissprot_id, seq_type="ref", anno_type="T", database='swissprot', blast_path=self.result_file["swissprot" + "trans_blast_path" + "ref"], exp_level=exp_level)
        self.update_db_record('sg_annotation_swissprot', swissprot_id, status="end", main_id=swissprot_id)

        self.remove_table_by_main_record(main_table='sg_annotation_pfam', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], detail_table_key='pfam_id')
        params_select_pfam = dict([('pfam_evalue',params_dict['pfam_evalue'])])
        params_select_pfam = json.dumps(params_select_pfam, sort_keys=True, separators=(',', ':'))
        pfam_id = self.add_annotation_pfam(name=None, params=params_select_pfam, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path'], seq_type="new", anno_type="T", exp_level=exp_level)
        self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=self.result_file['pfam_path' + "ref"], seq_type="ref", anno_type="T", exp_level=exp_level)
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

        self.update_db_record('sg_annotation_pfam', pfam_id, status="end", main_id=pfam_id)
        '''

        self.remove_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        if self.has_new:
            self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path'], seq_type = "new", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)
        #self.add_annotation_gene_query_denovo_detail(query_id=query_id, query_path=self.result_file['gene_query_path'], anno_type="G")
        self.remove_table_by_main_record(main_table='sg_annotation_cog', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_cog_detail'], detail_table_key='cog_id')
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
                self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=self.result_file['n_sum_path' + 'ref'], n_cog_path=self.result_file['n_sum_path'], seq_type="all", anno_type="T")
            self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=self.result_file['n_gene_sum_path' + 'ref'], n_cog_path=self.result_file['n_gene_sum_path'], seq_type="all", anno_type="G")
        self.update_db_record('sg_annotation_cog', cog_id, status="end", main_id=cog_id)


        self.remove_table_by_main_record(main_table='sg_annotation_go', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], detail_table_key='go_id')

        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path'])
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path'])
        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4_ref'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2_ref'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3_ref'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4_ref'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path_ref'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4_ref'])
        # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2_ref'])
        # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3_ref'])
        # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4_ref'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path_ref'])

        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=2, r_go_path=self.result_file['stat_level2_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=3, r_go_path=self.result_file['stat_level3_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=4, r_go_path=self.result_file['stat_level4_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=2, r_go_path=self.result_file['gene_stat_level2_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=3, r_go_path=self.result_file['gene_stat_level3_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=4, r_go_path=self.result_file['gene_stat_level4_all'])

        self.update_db_record('sg_annotation_go', go_id, status="end", main_id=go_id)

        self.remove_table_by_main_record(main_table='sg_annotation_kegg', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_kegg_categories', 'sg_annotation_kegg_level','sg_annotation_kegg_table', 'sg_annotation_kegg_pic'], detail_table_key='kegg_id')
        params_select_kegg = dict([(k,params_dict.get(k,None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=params_select_kegg, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path'])
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path'])
                # self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
            # self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path'])

        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path_ref'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path_ref'], png_dir=self.result_file['png_path_ref'])
            # self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path_ref'], png_dir=self.result_file['png_path_ref'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path_ref'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path_ref'], png_dir=self.result_file['gene_png_path_ref'])
        # self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path_ref'], png_dir=self.result_file['gene_png_path_ref'])
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path_ref'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path_ref'])

        if self.has_new:
            r_cate_path = self.result_file['layer_path_ref']
            n_cate_path = self.result_file['layer_path']
            r_gene_cate_path = self.result_file['gene_layer_path_ref']
            n_gene_cate_path = self.result_file['gene_layer_path']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="T", r_cate_path=r_cate_path, n_cate_path=n_cate_path)
            self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="G", r_cate_path=r_gene_cate_path, n_cate_path=n_gene_cate_path)

            #self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type='ref', anno_type="T", categories_path=self.result_file['layer_path_all'])
            #self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type='ref', anno_type="G", categories_path=self.result_file['gene_layer_path_all'])

            pathway_path = self.result_file['pathway_path' + '_all']
            png_path = self.result_file['png_path' + '_all']
            gene_pathway_path = self.result_file['gene_pathway_path' + '_all']
            gene_png_path =  self.result_file['gene_png_path' + '_all']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="T", level_path=pathway_path, png_dir=png_path)
                # self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="T", level_path=pathway_path, png_dir=png_path)
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="G", level_path=gene_pathway_path, png_dir=gene_png_path)
            # self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="G", level_path=gene_pathway_path, png_dir=gene_png_path)
        self.update_db_record('sg_annotation_kegg', kegg_id, status="end", main_id=kegg_id)

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
        self.get_trans2gene(trans2gene, trans2gene_ref)
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

        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['ref_stat_path'], venn_path=self.result_file['ref_venn_path'], seq_type = "ref", exp_level=exp_level)
        if self.has_new:
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['new_stat_path'], venn_path=self.result_file['new_venn_path'], seq_type = "new", exp_level=exp_level)
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=self.result_file['all_stat_path'], venn_path=self.result_file['all_venn_path'], seq_type = "all", exp_level=exp_level)

        self.update_db_record('sg_annotation_stat', stat_id, status="end", main_id=stat_id)
        if last_id:
            self.bind_object.logger.info("删除表格为 {}".format(last_id))
            self.remove_table_by_main_record(main_table='sg_annotation_stat', _id=last_id, detail_table=['sg_annotation_stat_detail'], detail_table_key='stat_id')
            self.bind_object.logger.info("删除表格成功 {}".format(last_id))

        params_select_nr = dict([(k,params_dict.get(k,None)) for k in ('nr_evalue', 'nr_similarity', 'nr_identity')])
        params_select_nr = json.dumps(params_select_nr, sort_keys=True, separators=(',', ':'))

        nr_old_id = self.get_table_by_main_record(main_table='sg_annotation_nr', task_id=task_id, type=self.anno_type)
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

        swissprot_old_id = self.get_table_by_main_record(main_table='sg_annotation_swissprot', task_id=task_id, type=self.anno_type)
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

        pfam_old_id = self.get_table_by_main_record(main_table='sg_annotation_pfam', task_id=task_id, type=self.anno_type)
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

        query_old_id = self.get_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type)
        self.remove_table_by_main_record(main_table='sg_annotation_query', task_id=task_id, type=self.anno_type, detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id, result_dir=self.result_dir)
        if self.has_new:
            self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path'], seq_type = "new", anno_type="T", exp_level=exp_level)
        self.add_annotation_query_denovo_detail(query_id=query_id, query_path=self.result_file['query_path' + "ref"], seq_type = "ref", anno_type="T", exp_level=exp_level)
        self.update_db_record('sg_annotation_query',query_id, status="end", main_id=query_id)

        if query_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_query', _id=query_old_id['_id'], detail_table=['sg_annotation_query_detail'], detail_table_key='query_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        cog_old_id = self.get_table_by_main_record(main_table='sg_annotation_cog', task_id=task_id, type=self.anno_type)
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
                self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=self.result_file['n_sum_path' + 'ref'], n_cog_path=self.result_file['n_sum_path'], seq_type="all", anno_type="T")
            self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=self.result_file['n_gene_sum_path' + 'ref'], n_cog_path=self.result_file['n_gene_sum_path'], seq_type="all", anno_type="G")
        self.update_db_record('sg_annotation_cog', cog_id, status="end", main_id=cog_id)
        if cog_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_cog', _id=cog_old_id['_id'], detail_table=['sg_annotation_cog_detail'], detail_table_key='cog_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        go_old_id = self.get_table_by_main_record(main_table='sg_annotation_go', task_id=task_id, type=self.anno_type)

        go_id = self.add_annotation_go(name=None, params=params_select_nr, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2'])
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3'])
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4'])
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path'])
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path'])
        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, level_path=self.result_file['stat_level2_ref'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2_ref'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3_ref'])
            # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=2, go_path=self.result_file['stat_level2_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=3, go_path=self.result_file['stat_level3_ref'])
            self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="T", level=4, go_path=self.result_file['stat_level4_ref'])
            self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="T", gos_path=self.result_file['gos_path_ref'])
        self.add_annotation_go_level(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, level_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3_ref'])
        self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4_ref'])
        # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=2, go_path=self.result_file['gene_stat_level2_ref'])
        # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=3, go_path=self.result_file['gene_stat_level3_ref'])
        # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="G", level=4, go_path=self.result_file['gene_stat_level4_ref'])
        self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="G", gos_path=self.result_file['gene_gos_path_ref'])

        if self.has_new:
            if exp_level.lower() == "transcript":
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=2, r_go_path=self.result_file['stat_level2_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=3, r_go_path=self.result_file['stat_level3_all'])
                self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="T", level=4, r_go_path=self.result_file['stat_level4_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=2, r_go_path=self.result_file['gene_stat_level2_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=3, r_go_path=self.result_file['gene_stat_level3_all'])
            self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="G", level=4, r_go_path=self.result_file['gene_stat_level4_all'])
        self.update_db_record('sg_annotation_go', go_id, status="end", main_id=go_id)
        if go_old_id:
            self.remove_table_by_main_record(main_table='sg_annotation_go', _id=go_old_id['_id'], detail_table=['sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], detail_table_key='go_id')
        else:
            self.bind_object.logger.info("未找到旧表不做删除")

        kegg_old_id = self.get_table_by_main_record(main_table='sg_annotation_kegg', task_id=task_id, type=self.anno_type)
        params_select_kegg = dict([(k,params_dict.get(k,None)) for k in ('kegg_evalue', 'kegg_similarity', 'kegg_identity')])
        params_select_kegg = json.dumps(params_select_kegg, sort_keys=True, separators=(',', ':'))
        kegg_id = self.add_annotation_kegg(name=None, params=params_select_kegg, result_dir=self.result_dir)
        if self.has_new:
            seq_type = "new"
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path'])
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path'])
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path'], png_dir=self.result_file['png_path'])
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path'], png_dir=self.result_file['gene_png_path'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path'])

        seq_type = "ref"
        if exp_level.lower() == "transcript":
            self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", categories_path=self.result_file['layer_path_ref'])
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path_ref'], png_dir=self.result_file['png_path_ref'])
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", level_path=self.result_file['pathway_path_ref'], png_dir=self.result_file['png_path_ref'])
            self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="T", table_path=self.result_file['table_path_ref'])
        self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", categories_path=self.result_file['gene_layer_path_ref'])
        self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path_ref'], png_dir=self.result_file['gene_png_path_ref'])
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", level_path=self.result_file['gene_pathway_path_ref'], png_dir=self.result_file['gene_png_path_ref'])
        self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="G", table_path=self.result_file['gene_table_path_ref'])

        if self.has_new:
            r_cate_path = self.result_file['layer_path_ref']
            n_cate_path = self.result_file['layer_path']
            r_gene_cate_path = self.result_file['gene_layer_path_ref']
            n_gene_cate_path = self.result_file['gene_layer_path']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="T", r_cate_path=r_cate_path, n_cate_path=n_cate_path)
            self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="G", r_cate_path=r_gene_cate_path, n_cate_path=n_gene_cate_path)
            pathway_path = self.result_file['pathway_path' + '_all']
            png_path = self.result_file['png_path' + '_all']
            gene_pathway_path = self.result_file['gene_pathway_path' + '_all']
            gene_png_path =  self.result_file['gene_png_path' + '_all']
            if exp_level.lower() == "transcript":
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="T", level_path=pathway_path, png_dir=png_path)
                self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="T", level_path=pathway_path, png_dir=png_path)
            self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="G", level_path=gene_pathway_path, png_dir=gene_png_path)
            self.add_annotation_kegg_pic(kegg_id=kegg_id, seq_type="all", anno_type="G", level_path=gene_pathway_path, png_dir=gene_png_path)

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
        self.add_annotation_stat_detail(stat_id=stat_id, stat_path=new_stat_path, venn_path=new_venn_path, task_id=self.task_id)
        blast_id = self.add_annotation_blast(name=None, params=params, stat_id=stat_id)
        blast_path = new_anno_path + "/anno_stat/blast"
        if os.path.exists(blast_path):
            for db in ["nr", "swissprot"]:
                trans_blast_path = blast_path + "/" + db + '.xls'
                gene_blast_path = blast_path + '/gene_' + db + '.xls'
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="transcript", database=db, blast_path=trans_blast_path)
                self.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="gene", database=db, blast_path=gene_blast_path)
        else:
            self.bind_object.set_error("没有blast的结果文件夹", code="537018114")
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
            self.bind_object.set_error("新序列NR注释结果文件不存在", code="537018115")
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
            self.bind_object.set_error("新序列Swiss-Prot注释结果文件不存在", code="537018116")
        pfam_id = self.add_annotation_pfam(name=None, params=params, stat_id=stat_id)
        gene_pfam_path = new_anno_path + "/anno_stat/pfam_stat/gene_pfam_domain"
        if os.path.exists(pfam_path) and os.path.exists(gene_pfam_path):
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_detail(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=pfam_path, seq_type="new", anno_type="transcript")
            self.add_annotation_pfam_bar(pfam_id=pfam_id, pfam_path=gene_pfam_path, seq_type="new", anno_type="gene")
        else:
            self.bind_object.set_error("pfam注释结果文件不存在", code="537018117")
        ref_stat_path = ref_anno_path + "/anno_stat/all_annotation_statistics.xls"
        ref_venn_path = ref_anno_path + "/anno_stat/venn"
        if os.path.exists(ref_stat_path) and os.path.exists(ref_venn_path):
            stat_id = self.add_annotation_stat(name=None, params=params, seq_type="ref" , database="cog,go,kegg")
            self.add_annotation_stat_detail(stat_id=stat_id, stat_path=ref_stat_path, venn_path=ref_venn_path, task_id=self.task_id)
        else:
            self.bind_object.set_error("已知序列注释统计文件和venn图文件夹不存在", code="537018118")
        query_id = self.add_annotation_query(name=None, params=params, stat_id=stat_id)
        query_path = ref_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = ref_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail(query_id=query_id, query_path=query_path, anno_type="transcript")
        self.add_annotation_gene_query_detail(query_id=query_id, query_path=gene_query_path, anno_type="gene")
        query_path = new_anno_path + "/anno_stat/trans_anno_detail.xls"
        gene_query_path = new_anno_path + "/anno_stat/gene_anno_detail.xls"
        self.add_annotation_query_detail(query_id=query_id, query_path=query_path, anno_type="transcript")
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
        self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=r_sum_path, n_cog_path=n_sum_path, seq_type="all", anno_type="transcript")
        self.add_annotation_cog_detail_all(cog_id=cog_id, r_cog_path=r_gene_sum_path, n_cog_path=n_gene_sum_path, seq_type="all", anno_type="gene")

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
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, go_path=stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, go_path=stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, go_path=stat_level4)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, go_path=gene_stat_level2)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, go_path=gene_stat_level3)
                self.add_annotation_go_detail(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, go_path=gene_stat_level4)
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=2, go_path=stat_level2)
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=3, go_path=stat_level3)
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="transcript", level=4, go_path=stat_level4)
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=2, go_path=gene_stat_level2)
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=3, go_path=gene_stat_level3)
                # self.add_annotation_go_graph(go_id=go_id, seq_type=seq_type, anno_type="gene", level=4, go_path=gene_stat_level4)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="transcript", gos_path=gos_path)
                self.add_annotation_go_list(go_id=go_id, seq_type=seq_type, anno_type="gene", gos_path=gene_gos_path)
            else:
                self.bind_object.set_error("GO注释的结果文件不存在", code="537018119")
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
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=2, r_go_path=r_stat_level2, n_go_path=n_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=3, r_go_path=r_stat_level3, n_go_path=n_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="transcript", level=4, r_go_path=r_stat_level4, n_go_path=n_stat_level4)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=2, r_go_path=r_gene_stat_level2, n_go_path=n_gene_stat_level2)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=3, r_go_path=r_gene_stat_level3, n_go_path=n_gene_stat_level3)
        self.add_annotation_go_all(go_id=go_id, seq_type="all", anno_type="gene", level=4, r_go_path=r_gene_stat_level4, n_go_path=n_gene_stat_level4)

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
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", categories_path=layer_path)
                self.add_annotation_kegg_categories(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", categories_path=gene_layer_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", level_path=pathway_path, png_dir=png_path)
                self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", level_path=gene_pathway_path, png_dir=gene_png_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="transcript", table_path=table_path)
                self.add_annotation_kegg_table(kegg_id=kegg_id, seq_type=seq_type, anno_type="gene", table_path=gene_table_path)
            else:
                self.bind_object.set_error("KEGG注释文件不存在", code="537018120")
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
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="transcript", r_cate_path=r_cate_path, n_cate_path=n_cate_path)
        self.add_annotation_kegg_categories_all(kegg_id=kegg_id, seq_type="all", anno_type="gene", r_cate_path=r_gene_cate_path, n_cate_path=n_gene_cate_path)
        pathway_path = merge_tran_output + "/pathway_table.xls"
        png_path = merge_tran_output + "/all_pathways"
        gene_pathway_path = merge_gene_output + "/pathway_table.xls"
        gene_png_path = merge_gene_output + "/all_pathways"
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="transcript", level_path=pathway_path, png_dir=png_path)
        self.add_annotation_kegg_level(kegg_id=kegg_id, seq_type="all", anno_type="gene", level_path=gene_pathway_path, png_dir=gene_png_path)


    @report_check
    def add_annotation_stat(self, name=None, params=None, seq_type=None, database=None, result_dir=None, taxon='Animals', exp_level="transcript"):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if self.anno_type == "latest":
            anno_str = "new"
        else:
            anno_str = self.anno_type

        insert_data = {
            'version': self.version,
            'exp_levl':exp_level.lower(),
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationStat_' + anno_str + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            "taxonomy": taxon,
            'has_new': self.has_new,
            'species_name': self.species_name,
            'status': 'start',
            'desc': '注释统计主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            # 'seq_type': seq_type,
            'database': database
        }
        collection = self.db['sg_annotation_stat']
        stat_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_stat!")
        return stat_id

    @report_check
    def add_annotation_stat_detail(self, stat_id, stat_path, venn_path, seq_type, exp_level, gene_exp=None, trans_exp=None, task_id=None):
        """
        database: 进行统计的数据库
        stat_path: all_annotation_statistics.xls
        venn_path: venn图目录
        """
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="537018121")
        if not os.path.exists(stat_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(stat_path), code="537018122")
        if not os.path.exists(venn_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(venn_path), code="537018123")

        if self.version == "v3":
            exp_g_set = set()
            if gene_exp is not None:
                with open(gene_exp, 'r') as gene_exp_f:
                    gene_exp_f.readline()
                    for line in gene_exp_f:
                        exps = line.strip().split("\t")[1:]
                        if sum(map(float, exps)) != 0:
                            exp_g_set.add(line.strip().split("\t")[0])
                        else:
                            pass
            exp_t_set = set()
            if trans_exp is not None:
                with open(trans_exp, 'r') as trans_exp_f:
                    trans_exp_f.readline()
                    for line in trans_exp_f:
                        exps = line.strip().split("\t")[1:]
                        if sum(map(float, exps)) != 0:
                            exp_t_set.add(line.strip().split("\t")[0])
                        else:
                            pass
            if seq_type == "ref":
                exp_t_set = exp_t_set & self.known_trans_set
                exp_g_set = exp_g_set & self.known_gene_set
            elif seq_type == "new":
                exp_t_set = exp_t_set & self.new_trans_set
                exp_g_set = exp_g_set & self.new_gene_set


        gene_dict = dict()
        trans_dict = dict()

        data_list = []
        tail = "gene"

        database_venn = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }
        tail = "tran"
        database_venn_tran = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }
        v3_db2db = {
            "pfam": 'Pfam',
            "kegg": "KEGG",
            "swissprot": "Swiss-Prot",
            "cog": "COG",
            "go": "GO",
            "nr": "NR",
            "annotation": "Total_anno",
            "total": "Total",
            "Pfam": 'Pfam',
            "KEGG": "KEGG",
            "Swiss-Prot": "Swiss-Prot",
            "COG": "COG",
            "GO": "GO",
            "NR": "NR",
            "Total_anno": "Total_anno",
            "Total": "Total"
        }
        with open(stat_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if database_venn.has_key(v3_db2db[line[0]]):
                    db = v3_db2db[line[0]]
                    venn = venn_path + "/" + database_venn_tran[v3_db2db[line[0]]]
                    gene_venn = venn_path + "/" + database_venn[v3_db2db[line[0]]]
                    if os.path.exists(venn) and os.path.exists(gene_venn):
                        with open(venn, "rb") as f:
                            venn_list = f.readline().strip('\n')
                            for line in f:
                                venn_list += ',{}'.format(line.strip('\n'))
                            trans_dict[db] = set(venn_list.split(","))
                        with open(gene_venn, "rb") as f:
                            gene_venn_list = f.readline().strip('\n')
                            for line in f:
                                gene_venn_list += ',{}'.format(line.strip('\n'))
                            gene_dict[db] = set(gene_venn_list.split(","))
        gene_intersection = reduce(lambda x, y: x & y, [gene_dict[i] for i in gene_dict])
        trans_intersection = reduce(lambda x, y: x & y, [trans_dict[i] for i in trans_dict])
        trans2enum = {t: str(n) for n, t in enumerate(trans_intersection)}
        genes2enum = {g: str(n) for n, g in enumerate(gene_intersection)}


        with open(stat_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if exp_level.lower() == "gene":
                    data = [
                        ('stat_id', stat_id),
                        ('type', v3_db2db[line[0]]),
                        ('gene', int(line[2])),
                        ('gene_percent', round(float(line[4]), 4)),
                        ('seq_type', seq_type),
                        ('anno_type', 'all')
                    ]
                else:
                    data = [
                        ('stat_id', stat_id),
                        ('type', v3_db2db[line[0]]),
                        ('transcript', int(line[1])),
                        ('gene', int(line[2])),
                        ('transcript_percent', round(float(line[3]), 4)),
                        ('gene_percent', round(float(line[4]), 4)),
                        ('seq_type', seq_type),
                        ('anno_type', 'all')
                    ]
                venn_list, gene_venn_list = None, None
                database = ["nr", "swissprot", "pfam", "kegg", "go", "string", "cog"]
                if database_venn.has_key(v3_db2db[line[0]]):
                    db = v3_db2db[line[0]]
                    if len(gene_dict[db]) > 150000 or len(trans_dict[db]) > 150000:
                        self.update_db_record('sg_annotation_stat', stat_id, created_geneset='no')
                        gene_list_all = list()
                        for i in gene_dict[db]:
                            if i in genes2enum:
                                gene_list_all.append(genes2enum[i])
                            else:
                                gene_list_all.append(i)

                        data.append(("gene_list", ','.join(gene_list_all)))
                        if exp_level.lower() == "transcript":
                            trans_list_all = list()
                            for i in trans_dict[db]:
                                if i in trans2enum:
                                    trans_list_all.append(trans2enum[i])
                                else:
                                    trans_list_all.append(i)
                            data.append(("transcript_list", ','.join(trans_list_all)))
                    else:
                        data.append(("gene_list", ','.join(list(gene_dict[db]))))
                        if exp_level.lower() == "transcript":
                            data.append(("transcript_list", ','.join(list(trans_dict[db]))))
                    # venn = venn_path + "/" + database_venn_tran[v3_db2db[line[0]]]
                    # gene_venn = venn_path + "/" + database_venn[v3_db2db[line[0]]]
                    # if os.path.exists(venn) and os.path.exists(gene_venn):
                    #     with open(venn, "rb") as f:
                    #         venn_list = f.readline().strip('\n')
                    #         for line in f:
                    #             venn_list += ',{}'.format(line.strip('\n'))
                    #         trans_dict[db] = set(venn_list.split(","))
                    #     with open(gene_venn, "rb") as f:
                    #         gene_venn_list = f.readline().strip('\n')
                    #         for line in f:
                    #             gene_venn_list += ',{}'.format(line.strip('\n'))
                    #         gene_dict[db] = set(gene_venn_list.split(","))
                    #     if len(gene_dict[db]) > 100000:
                    #         gene_venn_list = ','.join(random.sample(gene_dict[db], 100000))
                    #     data.append(("gene_list", gene_venn_list))
                    #     if exp_level.lower() == "transcript":
                    #         if len(trans_dict[db]) > 100000:
                    #             venn_list = ','.join(random.sample(trans_dict[db], 100000))
                    #         data.append(("transcript_list", venn_list))
                    # else:
                    #     self.bind_object.set_error("%s对应的venn.txt文件%s %s不存在", variables=(line[0], venn, gene_venn), code="537018124")
                data = SON(data)
                data_list.append(data)

        if self.version == "v3":
            anno_gene_list = set()
            anno_trans_list = set()
            go_trans_set = exp_t_set & set(trans_dict['GO'])
            go_gene_set = exp_g_set & set(gene_dict["GO"])
            KEGG_trans_set = exp_t_set & set(trans_dict['KEGG'])
            KEGG_gene_set = exp_g_set & set(gene_dict["KEGG"])
            COG_trans_set = exp_t_set & set(trans_dict['COG'])
            COG_gene_set = exp_g_set & set(gene_dict["COG"])
            NR_trans_set = exp_t_set & set(trans_dict['NR'])
            NR_gene_set = exp_g_set & set(gene_dict["NR"])
            Swiss_Prot_trans_set = exp_t_set & set(trans_dict["Swiss-Prot"])
            Swiss_Prot_gene_set = exp_g_set & set(gene_dict["Swiss-Prot"])
            Pfam_trans_set = exp_t_set & set(trans_dict['Pfam'])
            Pfam_gene_set = exp_g_set & set(gene_dict["Pfam"])
            intersection_all_gene = go_gene_set & KEGG_gene_set & COG_gene_set & NR_gene_set & Swiss_Prot_gene_set & Pfam_gene_set
            intersection_all_trans = go_trans_set & KEGG_trans_set & COG_trans_set & NR_trans_set & Swiss_Prot_trans_set & Pfam_trans_set
            trans2enum = {t: str(n) for n, t in enumerate(intersection_all_trans)}
            genes2enum = {g: str(n) for n, g in enumerate(intersection_all_gene)}

            for db in ['GO', 'KEGG', 'COG', 'NR', 'Swiss-Prot', 'Pfam']:
                trans_set = exp_t_set & set(trans_dict[db])
                gene_set = exp_g_set & set(gene_dict[db])
                # if len(trans_set) < 200000 and len(gene_set) < 200000:
                if len(trans_set) > 150000 or len(gene_set) > 150000:
                    self.update_db_record('sg_annotation_stat', stat_id, created_geneset='no')
                    trans_list = list()
                    for i in trans_set:
                        if i in trans2enum:
                            trans_list.append(trans2enum[i])
                        else:
                            trans_list.append(i)
                    trans_set_new = set(trans_list)
                    gene_list = list()
                    for j in gene_set:
                        if j in genes2enum:
                            gene_list.append(genes2enum[j])
                        else:
                            gene_list.append(j)
                    gene_set_new = set(gene_list)
                    data = [
                        ('stat_id', stat_id),
                        ('type', db),
                        ('transcript', len(trans_set)),
                        ('gene', len(gene_set)),
                        ('transcript_percent', round(float(len(trans_set))/float(len(exp_t_set) if len(exp_t_set) else 1), 4)),
                        ('gene_percent', round(float(len(gene_set))/float(len(exp_g_set) if len(exp_g_set) else 1), 4)),
                        ('seq_type', seq_type),
                        ('anno_type', 'exp'),
                        ('gene_list', ",".join(list(gene_set_new))),
                        ('transcript_list', ",".join(list(trans_set_new)))
                    ]
                else:
                    data = [
                        ('stat_id', stat_id),
                        ('type', db),
                        ('transcript', len(trans_set)),
                        ('gene', len(gene_set)),
                        ('transcript_percent',
                         round(float(len(trans_set)) / float(len(exp_t_set) if len(exp_t_set) else 1), 4)),
                        ('gene_percent',
                         round(float(len(gene_set)) / float(len(exp_g_set) if len(exp_g_set) else 1), 4)),
                        ('seq_type', seq_type),
                        ('anno_type', 'exp'),
                        ('gene_list', ",".join(list(gene_set))),
                        ('transcript_list', ",".join(list(trans_set)))
                    ]
                anno_gene_list |= gene_set
                anno_trans_list |= trans_set
                data = SON(data)
                data_list.append(data)
            data = [
                ('stat_id', stat_id),
                ('type', 'Total_anno'),
                ('transcript', len(anno_trans_list)),
                ('gene', len(anno_gene_list)),
                ('transcript_percent', round(float(len(anno_trans_list))/float(len(exp_t_set) if len(exp_t_set) else 1), 4)),
                ('gene_percent', round(float(len(anno_gene_list))/float(len(exp_g_set) if len(exp_g_set) else 1), 4)),
                ('seq_type', seq_type),
                ('anno_type', 'exp')
            ]
            data = SON(data)
            data_list.append(data)
            data = [
                ('stat_id', stat_id),
                ('type', 'Total'),
                ('transcript', len(exp_t_set)),
                ('gene', len(exp_g_set)),
                ('transcript_percent', "1.0"),
                ('gene_percent', "1.0"),
                ('seq_type', seq_type),
                ('anno_type', 'exp')
            ]
            data = SON(data)
            data_list.append(data)

        collection = self.db['sg_annotation_stat_detail']
        collection.insert_many(data_list)

        # try:
        #     collection = self.db['sg_annotation_stat_detail']
        #     collection.insert_many(data_list)
        # except Exception, e:
        #     self.bind_object.set_error("导入注释统计信息失败:%s" , variables=( stat_path), code="537018125")
        # else:
        #     self.bind_object.logger.info("导入注释统计信息成功：%s" % (stat_path))

    @report_check
    def add_stat_detail(self, old_stat_id, stat_id, nr_evalue, gene_nr_evalue, sw_evalue, gene_sw_evalue):
        """
        注释重运行时注释统计导表sg_annotation_stat_detail
        """
        if not isinstance(old_stat_id, ObjectId):
            if isinstance(old_stat_id, types.StringTypes):
                old_stat_id = ObjectId(old_stat_id)
            else:
                self.bind_object.set_error('old_stat_id必须为ObjectId对象或其对应的字符串！', code="537018126")
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="537018127")
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
            self.bind_object.set_error("导入注释统计信息出错", code="537018128")
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
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="537018129")
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationBlast_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
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
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="537018130")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="537018131")
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
                            data.update({'gene_id': self.trans_gene[line[5]]})
                            data.update({'is_gene': self.trans_isgene[line[5]]})
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
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="537018132")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="537018133")
        data_list = []
        with open(blast_path, 'r') as f:
            lines = f.readlines()
            flag = None
            for line in lines[1:]:
                line = line.strip("\n").split('\t')
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
                            data.update({'gene_id': self.trans_gene[line[5]]})
                            data.update({'is_gene': self.trans_isgene[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    else:
                        pass
                    if exp_level.lower() == "gene" and self.trans_isgene[line[5]] == False:
                        continue
                    else:
                        data_list.append(SON(data))
        collection = self.db['sg_annotation_nr_detail']
        if data_list:
            try:
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入nrblast信息：%s成功!" % (blast_path))
            except Exception, e:
                self.bind_object.set_error("导入注释统计信息失败:%s" , variables=( blast_path), code="537018134")


    @report_check
    def add_annotation_blast_swissprot_detail(self, blast_id, seq_type, anno_type, database, blast_path, exp_level):
        if not isinstance(blast_id, ObjectId):
            if isinstance(blast_id, types.StringTypes):
                blast_id = ObjectId(blast_id)
            else:
                self.bind_object.set_error('blast_id必须为ObjectId对象或其对应的字符串！', code="537018135")
        if not os.path.exists(blast_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(blast_path), code="537018136")
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
                            data.update({'gene_id': self.trans_gene[line[5]]})
                            data.update({'is_gene': self.trans_isgene[line[5]]})
                        except:
                            data.update({'gene_id': line[5]})
                            data.update({'is_gene': True})
                    else:
                        pass
                    if exp_level.lower() == "gene" and self.trans_isgene[line[5]] == False:
                        continue
                    else:
                        data_list.append(SON(data))
        collection = self.db['sg_annotation_swissprot_detail']
        if data_list:
            try:
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入swissprot blast信息：%s成功!" % (blast_path))
            except Exception, e:
                self.bind_object.set_error("导入swissprot注释统计信息失败:%s" , variables=( blast_path), code="537018137")

    @report_check
    def add_annotation_nr(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="537018138")
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationNr_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
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
                self.bind_object.set_error('nr_id必须为ObjectId对象或其对应的字符串！', code="537018139")
        if not os.path.exists(evalue_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="537018140")
        if not os.path.exists(similar_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="537018141")
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
                self.bind_object.set_error("导入nr库注释作图信息evalue,similar：%s、%s出错!" , variables=(evalue_path, similar_path), code="537018142")
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
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="537018143")
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationSwissprot_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
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
                self.bind_object.set_error('swissprot_id必须为ObjectId对象或其对应的字符串！', code="537018144")
        if not os.path.exists(evalue_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(evalue_path), code="537018145")
        if not os.path.exists(similar_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(similar_path), code="537018146")
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
                self.bind_object.set_error("导入swissprot库注释作图信息evalue,similar：%s、%s出错!" , variables=(evalue_path, similar_path), code="537018147")
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
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="537018148")
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationPfam_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
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
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="537018149")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="537018150")
        data_list = []
        with open(pfam_path, "r") as f:
            lines = f.readlines()
            last_seq_id = ''
            last_pfam_id = ''
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[2] != last_seq_id or line[3] != last_pfam_id:
                    # 过滤不在参考范围的转录本注释
                    if exp_level.lower() == "gene" or line[0] not in self.trans_isgene:
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
                            data.append(('gene_id', self.trans_gene[line[0]]))
                            data.append(('is_gene', self.trans_isgene[line[0]]))
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
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="537018151")
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
                self.bind_object.set_error('pfam_id必须为ObjectId对象或其对应的字符串！', code="537018152")
        if not os.path.exists(pfam_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(pfam_path), code="537018153")
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
                self.bind_object.set_error("导入pfam注释信息:%s失败！" , variables=( pfam_path), code="537018154")
            else:
                self.bind_object.logger.info("导入pfam注释信息:%s成功！" % pfam_path)

    @report_check
    def add_annotation_cog(self, name=None, params=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationCog_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'cog注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        collection = self.db['sg_annotation_cog']
        cog_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_cog!")
        return cog_id

    @report_check
    def add_annotation_cog_detail(self, cog_id, cog_path, seq_type, anno_type):
        '''
        cog_path: cog_summary.xls
        seq_type: ref/new
        anno_type: transcript/gene
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                self.bind_object.set_error('cog_id必须为ObjectId对象或其对应的字符串！', code="537018155")
        if not os.path.exists(cog_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(cog_path), code="537018156")
        data_list = list()
        with open(cog_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('cog_id', cog_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('type', line[0]),
                    ('function_categories', "[" + line[2] + "]" + " " + line[1]),
                    ('cog', int(line[3])),
                 ]
                try:
                    data.append(('cog_list', line[4]))
                except:
                    data.append(('cog_list', None))
                # try:
                #     data.append(('nog_list', line[5]))
                # except:
                #     data.append(('nog_list', None))
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入cog注释信息：%s出错!" , variables=(cog_path), code="537018157")
            else:
                self.bind_object.logger.info("导入cog注释信息：%s成功!" % (cog_path))

    @report_check
    def add_annotation_cog_detail_all(self, cog_id, r_cog_path, n_cog_path, seq_type, anno_type):
        '''
        r_cog_path: cog_summary.xls(ref)
        n_cog_path: cog_summary.xls(new)
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                self.bind_object.set_error('cog_id必须为ObjectId对象或其对应的字符串！', code="537018158")
        if not os.path.exists(r_cog_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_cog_path), code="537018159")
        if not os.path.exists(n_cog_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_cog_path), code="537018160")
        first = ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING', 'METABOLISM', 'POORLY CHARACTERIZED']
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
        data_list = list()
        cate = {}
        funlist = {'COG': {}, 'NOG': {}}
        cog_fun, nog_fun = {}, {}
        with open(r_cog_path, 'r') as f, open(n_cog_path, 'r') as n:
            lines = f.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                line[1] = "[" + line[2] + "]" + " " + line[1]
                m = re.match(r"\[(.+)\].+$", line[1])
                if m:
                    fun1 = m.group(1)
                    try:
                        funlist['COG'][fun1] = line[4].split(";")
                    except:
                        funlist['COG'][fun1] = []
                    try:
                        funlist['NOG'][fun1] = line[5].split(";")
                    except:
                        funlist['NOG'][fun1] = []
            for item in items[1:]:
                item = item.strip().split('\t')
                item[1] = "[" + item[2] + "]" + " " + item[1]
                m = re.match(r"\[(.+)\].+$", item[1])
                if m:
                    fun1 = m.group(1)
                    if fun1 in funlist['COG']:
                        try:
                            cog_ids = item[4].split(";")
                            if cog_ids:
                                for cog in cog_ids:
                                    if cog not in funlist['COG'][fun1]:
                                        funlist['COG'][fun1].append(cog)
                        except:
                            pass
                    else:
                        try:
                            funlist['COG'][fun1] = item[4].split(";")
                        except:
                            funlist['COG'][fun1] = []
                    if fun1 in funlist['NOG']:
                        try:
                            nog_ids = item[5].split(";")
                            if nog_ids:
                                for nog in nog_ids:
                                    if nog not in funlist['NOG'][fun1]:
                                        funlist['NOG'][fun1].append(nog)
                        except:
                            pass
                    else:
                        try:
                            funlist['NOG'][fun1] = item[5].split(";")
                        except:
                            funlist['NOG'][fun1] = []
        for g in funlist['COG'].keys():
            detail = func_decs[g]
            category = '[' + g + ']' + ' ' + detail
            try:
                cog_list = list(set(funlist['COG'][g]))
            except:
                cog_list = []
            try:
                nog_list = list(set(funlist['NOG'][g]))
            except:
                nog_list = []
            for cog_class in func_type.keys():
                if g in func_type[cog_class]:
                    cog_type = cog_class
            data = [
                ('cog_id', cog_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('type', cog_type),
                ('function_categories', category),
                ('cog', len([x for x in cog_list if x ])),
                # ('nog', len([x for x in nog_list if x])),
                ('cog_list', ';'.join(cog_list)),
                # ('nog_list', ';'.join(nog_list)),
            ]
            data = SON(data)
            if len(cog_list) + len(nog_list) != 0:
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_cog_detail']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入cog注释all出错：%s， %s" , variables=(r_cog_path, n_cog_path), code="537018161")
            else:
                self.bind_object.logger.info("导入cog注释all成功：%s, %s" % (r_cog_path, n_cog_path))

    @report_check
    def add_annotation_cog_table(self, cog_id, table_path, seq_type, anno_type):
        '''
        table_path:cog_table.xls
        '''
        if not isinstance(cog_id, ObjectId):
            if isinstance(cog_id, types.StringTypes):
                cog_id = ObjectId(cog_id)
            else:
                self.bind_object.set_error('cog_id必须为ObjectId对象或其对应的字符串！', code="537018162")
        if not os.path.exists(table_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(table_path), code="537018163")
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
            self.bind_object.set_error("导入cog注释table信息：%s出错!" , variables=(table_path), code="537018164")
        else:
            self.bind_object.logger.info("导入cog注释table信息：%s成功!" % (table_path))

    @report_check
    def add_annotation_go(self, name=None, params=None, result_dir=None):
        """
        go注释导表函数
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationGo_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'go注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        collection = self.db['sg_annotation_go']
        go_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_go!")
        return go_id

    def add_annotation_go_detail(self, go_id, seq_type, anno_type, level, go_path):
        """
        go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="537018165")
        if not os.path.exists(go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(go_path), code="537018166")
        data_list = list()
        with open(go_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('goterm', line[0]),
                    ('goterm_2', line[1]),
                    ('goid_2', line[2]),
                    ('seq_number', int(line[-3])),
                    ('percent', round(float(line[-2]), 4)),
                    #('seq_list', line[-1])
                ]
                if level == 2:
                    data.append(('seq_list', line[-1]))
                if level >= 3:
                    data.append(('goterm_3', line[3]))
                    data.append(('goid_3', line[4]))
                if level == 4:
                    data.append(('goterm_4', line[5]))
                    data.append(('goid_4', line[6]))
                data = SON(data)
                data_list.append(data)
            if data_list:
                try:
                    collection = self.db['sg_annotation_go_detail']
                    collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.set_error("导入go注释信息：%s出错!" , variables=(go_path), code="537018167")
                else:
                    self.bind_object.logger.info("导入go注释信息：%s成功!" % (go_path))

    @report_check
    def add_annotation_go_graph(self, go_id, seq_type, anno_type, level, go_path):
        """
        go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="537018168")
        if not os.path.exists(go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(go_path), code="537018169")
        data_list = list()
        self.bind_object.logger.info("开始导入go注释画图信息：%s" % (go_path))
        with open(go_path, 'r') as f:
            lines = f.readlines()
            term = {}
            term_list = []
            for i in range(1, len(lines)):
                line = lines[i].strip().split('\t')
                if level == 2:
                    term_type = line[0]
                    go_term = line[1]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
                if level == 3:
                    term_type = line[0]
                    go_term = line[3]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
                if level == 4:
                    term_type = line[0]
                    go_term = line[5]
                    if go_term not in term:
                        term[go_term] = []
                        term[go_term].append(i)
                    else:
                        term[go_term].append(i)
        with open(go_path, 'r') as f:
            lines = f.readlines()
            for item in term:
                seq_list = []
                for j in term[item]:
                    line = lines[j].strip().split("\t")
                    term_type = line[0]
                    '''
                    for seq in line[-1].split(";"):
                        if seq not in seq_list:
                            seq_list.append(seq)
                    '''
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('term_type', term_type),
                    ('go_term', item),
                    ('seq_number', line[-3]),
                    ('percent', line[-2]),
                    #('seq_list', seq_list)
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_go_graph']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入go注释画图信息出错：%s" , variables=(go_path), code="537018170")
            else:
                self.bind_object.logger.info("导入go注释画图信息成功：%s" % (go_path))

    @report_check
    def add_annotation_go_level(self, go_id, seq_type, anno_type, level, level_path):
        """
        level_path: go2level.xls
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="537018171")
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="537018172")
        data_list = list()
        with open(level_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('level', level),
                    ('term_type', line[1]),
                    ('parent_name', line[0]),
                    ('num', int(line[3])),
                    ('percent', round(float(line[4]), 4)),
                    ('go', line[2]),
                    #('seq_list', line[5]),
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_go_level']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入go注释第二层级信息：%s出错!" , variables=(level_path), code="537018173")
            else:
                self.bind_object.logger.info("导入go注释第二层级信息：%s成功!" % (level_path))

    @report_check
    def add_annotation_go_list(self, go_id, seq_type, anno_type, gos_path):
        """
        gos_path: go.list
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id须为ObjectId对象或其他对应的字符串！', code="537018174")
        if not os.path.exists(gos_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(gos_path), code="537018175")
        data_list = []
        with open(gos_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                data = [
                    ('go_id', go_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('gene_id', line[0]),
                    ('gos_list', line[1]),
                ]
                data = SON(data)
                data_list.append(data)
            if data_list:
                try:
                    collection = self.db['sg_annotation_go_list']
                    collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.set_error("导入gos_list注释信息：%s出错:%s" , variables=(gos_path), code="537018176")
                else:
                    self.bind_object.logger.info("导入gos_list注释信息：%s成功!" % (gos_path))

    def add_annotation_go_all(self, go_id, seq_type, anno_type, level, r_go_path):
        """
        r_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(ref)
        n_go_path: go1234level_statistics.xls/go123level_statistics.xls/go12level_statistics.xls(new)
        """
        if not isinstance(go_id, ObjectId):
            if isinstance(go_id, types.StringTypes):
                go_id = ObjectId(go_id)
            else:
                self.bind_object.set_error('go_id必须为ObjectId对象或其对应的字符串！', code="537018177")
        if not os.path.exists(r_go_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_go_path), code="537018178")
        data_list1, data_list2, query_ids = list(), list(), list()
        funlist, termlist = {}, {}
        with open(r_go_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                fun = line[0] + "|||" + line[1] + "|||" + line[2]
                term = line[0] + "|||" + line[1]
                if level == 3:
                    fun += "|||" + line[3] + "|||" + line[4]
                    term = line[0] + "|||" + line[3]
                if level == 4:
                    fun += "|||" + line[3] + "|||" + line[4]
                    fun += "|||" + line[5] + "|||" + line[6]
                    term = line[0] + "|||" + line[5]
                funlist[fun] = line[-1].split(";")
                if term not in termlist:
                    termlist[term] = set(line[-1].split(";"))
                else:
                    for q in line[-1].split(";"):
                        if q not in termlist[term]:
                            termlist[term].add(q)
                query_ids.extend(line[-1].split(";"))
        query_ids = list(set(query_ids))
        for term in sorted(termlist):
            terms = term.split("|||")
            seq_list = termlist[term]
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
                #('seq_list', ";".join(seq_list))
            ]
            data = SON(data)
            data_list1.append(data)
        if data_list1:
            try:
                collection = self.db['sg_annotation_go_graph']
                # collection.insert_many(data_list1)
            except Exception, e:
                self.bind_object.set_error("导入go注释画图all信息出错：%s" % (r_go_path))
                print "导入go注释画图all出错：%s" % (r_go_path)
            else:
                self.bind_object.logger.info("导入go注释画图all信息成功：%s" % (r_go_path))
        for fun in funlist:
            terms = fun.split("|||")
            data = [
                ('go_id', go_id),
                ('seq_type', seq_type),
                ('anno_type', anno_type),
                ('level', level),
                ('goterm', terms[0]),
                ('goterm_2', terms[1]),
                ('goid_2', terms[2])
            ]
            if level >= 3:
                data.append(('goterm_3', terms[3]))
                data.append(('goid_3', terms[4]))
            if level == 4:
                data.append(('goterm_4', terms[5]))
                data.append(('goid_4', terms[6]))
            seq_list = funlist[fun]
            percent = float(len(funlist[fun])) / len(query_ids)
            data.append(('seq_number', len(seq_list)))
            data.append(('percent', round(percent, 4)))
            if level == 2:
                data.append(('seq_list', ";".join(seq_list)))
            data = SON(data)
            data_list2.append(data)
        if data_list2:
            try:
                collection = self.db['sg_annotation_go_detail']
                collection.insert_many(data_list2)
            except Exception, e:
                self.bind_object.set_error("导入go注释all出错：%s" , variables=(r_go_path), code="537018179")
            else:
                self.bind_object.logger.info("导入go注释all信息成功：%s" % (r_go_path))

    @report_check
    def add_annotation_kegg(self, name=None, params=None, result_dir=None):
        """
        kegg注释导表函数
        """
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationKegg_' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': 'kegg注释结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'categories': ["M","CP","EIP","GIP","OS","HD","DD"]
        }
        collection = self.db['sg_annotation_kegg']
        kegg_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_kegg!")
        return kegg_id

    @report_check
    def add_annotation_kegg_categories(self, kegg_id, seq_type, anno_type, categories_path):
        """
        categories_path:kegg_layer.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="537018180")
        if not os.path.exists(categories_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(categories_path), code="537018181")
        data_list = list()
        first_type = []
        with open(categories_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                type_abr = ''.join([x[0] for x in line[0].split(' ')])
                first_type.append(type_abr)
                data = [
                    ('kegg_id', kegg_id),
                    ('seq_type', seq_type),
                    ('anno_type', anno_type),
                    ('first_category', line[0]),
                    ('second_category', line[1]),
                    ('num', int(line[2])),
                    ('seq_list', line[3]),
                ]
                data = SON(data)
                data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_categories']
                collection.insert_many(data_list)
                self.update_db_record('sg_annotation_kegg', kegg_id, categories=list(set(first_type)))
            except Exception, e:
                self.bind_object.set_error("导入kegg注释分类信息：%s出错!" , variables=(categories_path), code="537018182")
            else:
                self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % categories_path)

    @report_check
    def add_annotation_kegg_level(self, kegg_id, seq_type, anno_type, level_path, png_dir):
        """
        level_path: pathway_table.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="537018183")
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="537018184")
        if not os.path.exists(png_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(png_dir), code="537018185")
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                # fs = gridfs.GridFS(self.db)
                pid = re.sub('path:', '', line[0])
                # pdfid = fs.put(open(png_dir + '/' + pid + '.pdf', 'rb'))
                # graph_png_id = fs.put(open(png_dir + '/' + pid + '.png', 'rb'))
                insert_data = {
                    'kegg_id': kegg_id,
                    'seq_type': seq_type,
                    'anno_type': anno_type,
                    'pathway_id': line[0],
                    'first_category': line[1],
                    'second_category': line[2],
                    'pathway_definition': line[3],
                    'number_of_seqs': int(line[4]),
                    'seq_list': line[5],
                    # 'graph_id': pdfid,
                    # 'graph_png_id': graph_png_id,
                    'hyperlink': line[-1]
                }
                data_list.append(insert_data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_level']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释层级信息：%s、%s出错!" , variables=(level_path, png_dir), code="537018186")
            else:
                self.bind_object.logger.info("导入kegg注释层级信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def add_annotation_kegg_pic(self, kegg_id, seq_type, anno_type, level_path, png_dir):
        """
        level_path: pathway_table.xls
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="537018187")
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="537018188")
        if not os.path.exists(png_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(png_dir), code="537018189")
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 8:
                                [png, shape, bg_color, fg_color, coords, title, kos, href] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            insert_data = {
                                'kegg_id': kegg_id,
                                'seq_type': seq_type,
                                'anno_type': anno_type,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[0]))

        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_pic']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!" , variables=(level_path, png_dir), code="537018190")
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))


    @report_check
    def add_annotation_kegg_table(self, kegg_id, seq_type, anno_type, table_path):
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="537018191")
        if not os.path.exists(table_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(table_path), code="537018192")
        with open(table_path, 'rb') as r:
            data_list = []
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                insert_data = {
                    'kegg_id': kegg_id,
                    'seq_type': seq_type,
                    'anno_type': anno_type,
                    'transcript_id': line[0],
                    'ko_id': line[1],
                    'ko_name': line[2],
                    'hyperlink': line[3],
                    'paths': line[4],
                }
                data_list.append(insert_data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_table']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释table信息：%s出错!" , variables=(table_path), code="537018193")
            else:
                self.bind_object.logger.info("导入kegg注释table信息：%s成功!" % (table_path))

    @report_check
    def add_annotation_kegg_categories_all(self, kegg_id, seq_type, anno_type, r_cate_path, n_cate_path):
        """
        r_cate_path:kegg_layer.xls(ref)
        n_cate_path:kegg_layer.xls(new)
        """
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="537018194")
        if not os.path.exists(r_cate_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_cate_path), code="537018195")
        if not os.path.exists(n_cate_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_cate_path), code="537018196")
        data_list = list()
        cate = {}
        with open(r_cate_path, 'r') as f, open(n_cate_path, "r") as n:
            lines = f.readlines()
            items = n.readlines()
            for line in lines:
                line = line.strip().split('\t')
                if line[0] == "Metabolism" and line[1] == "Global and overview maps":
                    pass
                else:
                    cate_d = line[0] + "|||" + line[1]
                    cate[cate_d] = line[3].split(";")
            for item in items:
                item = item.strip().split('\t')
                if item[0] == "Metabolism" and item[1] == "Global and overview maps":
                    pass
                else:
                    cate_d = item[0] + "|||" + item[1]
                    if cate_d not in cate:
                        cate[cate_d] = item[3].split(";")
                    else:
                        ids = item[3].split(";")
                        for q in ids:
                            if q not in cate[cate_d]:
                                cate[cate_d].append(q)
        for f in ["Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes", "Organismal Systems", "Human Diseases", "Drug Development"]:
            for c in cate:
                ca = c.split("|||")
                first_category = ca[0]
                second_category = ca[1]
                if f == first_category:
                    num = len(cate[c])
                    seq_list = ";".join(cate[c])
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
                    data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_annotation_kegg_categories']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释分类alls出错：%s, %s" , variables=(r_cate_path, n_cate_path), code="537018197")
            else:
                self.bind_object.logger.info("导入kegg注释分类all成功：%s, %s" % (r_cate_path, n_cate_path))

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
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="537018198")
        if not os.path.exists(r_level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(r_level_path), code="537018199")
        if not os.path.exists(n_level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(n_level_path), code="537018200")
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
                self.bind_object.set_error("导入kegg注释层级all信息出错：%s、%s" , variables=(r_level_path, n_level_path), code="537018201")
            else:
                self.bind_object.logger.info("导入kegg注释层级all信息成功：%s、%s" % (level_path, png_dir))

    @report_check
    def add_annotation_query(self, name=None, params=None, stat_id=None, result_dir=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn
        if not isinstance(stat_id, ObjectId):
            print stat_id
            if isinstance(stat_id, types.StringTypes):
                stat_id = ObjectId(stat_id)
            else:
                self.bind_object.set_error('stat_id必须为ObjectId对象或其对应的字符串！', code="537018202")
        insert_data = {
            'version': self.version,
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'AnnotationQuery' + self.anno_type + '_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'type': self.anno_type,
            'params': params,
            'result_dir': result_dir,
            'status': 'start',
            'desc': '注释查询主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'stat_id': stat_id
        }
        collection = self.db['sg_annotation_query']
        query_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add ref_annotation_query!")
        return query_id

    def get_kegg_map2class(self):
        '''
        返回kegg mapid 与 classI classII 对应关系字典
        key: map00010
        value: (classi_name, classii_name)
        '''
        map2class = dict()
        with open(self.kegg_json, "rb") as f:
            root = json.load(f)
        classI = root['children']

        for class1 in classI:
            class1_name = class1['name']
            for class2 in class1['children']:
                class2_name = class2['name']
                paths = class2['children']
                for path in paths:
                    path_name = "map" + str(path['name']).split(" ")[0]
                    map2class[path_name] = (class1_name, class2_name)
        return map2class

    @report_check
    def add_annotation_query_denovo_detail(self, query_id, query_path, seq_type, anno_type, exp_level):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537018203")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537018204")
        data_list = []
        map2class = self.get_kegg_map2class()
        with open(query_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                is_gene = True
                if self.trans_isgene.has_key(line[1]) and self.trans_isgene[line[1]] == False:
                    is_gene = False
                if exp_level.lower() == "gene" and is_gene == False:
                    continue
                data = [
                    ('query_id', query_id),
                    #('anno_type', anno_type),
                    ('seq_type', seq_type),
                    ('transcript_id', line[1]),
                    ('gene_id', line[0]),
                    ('is_gene', is_gene),
                    ('gene_name', line[3]),
                ]
                try:
                    data.append(('length', line[4]))
                except:
                    data.append(('length', None))
                try:
                    data.append(('description', line[5]))
                except:
                    data.append(('description', None))
                try:
                    data.append(('cog', line[6]))
                    data.append(('cog_description', line[7]))
                except:
                    data.append(('cog', None))
                    data.append(('cog_description', None))
                # try:
                #     data.append(('nog', line[4]))
                #     data.append(('nog_description', line[6]))
                # except:
                #     data.append(('nog', None))
                #     data.append(('nog_description', None))
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
                    paths = [pathway.split("(")[0] for pathway in line[10].split("; ")]
                    if line[10] != "":
                        pathway_class1 = [map2class[path][0] for path in paths if path in map2class]
                        pathway_class2 = [map2class[path][1] for path in paths if path in map2class]
                    data.append(('pathways_class1', ";".join(pathway_class1)))
                    data.append(('pathways_class2', ";".join(pathway_class2)))

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
                try:
                    data.append(('enterz', line[15]))
                except:
                    data.append(('enterz', None))
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db['sg_annotation_query_detail']
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入转录本注释统计信息：%s出错!" , variables=(query_path), code="537018205")
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))


    @report_check
    def add_annotation_query_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537018206")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537018207")
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
            self.bind_object.set_error("导入转录本注释统计信息：%s出错!" , variables=(query_path), code="537018208")
        else:
            self.bind_object.logger.info("导入转录本注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_detail(self, query_id, query_path, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537018209")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537018210")
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
            self.bind_object.set_error("导入基因注释查询信息：%s出错!" , variables=(query_path), code="537018211")
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

    @report_check
    def add_annotation_gene_query_denovo_detail(self, query_id, query_path, seq_type, anno_type):
        if not isinstance(query_id, ObjectId):
            if isinstance(query_id, types.StringTypes):
                query_id = ObjectId(query_id)
            else:
                self.bind_object.set_error('query_id必须为ObjectId对象或其对应的字符串！', code="537018212")
        if not os.path.exists(query_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(query_path), code="537018213")
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
            self.bind_object.set_error("导入基因注释查询信息：%s出错!" , variables=(query_path), code="537018214")
        else:
            self.bind_object.logger.info("导入基因注释统计信息：%s成功!" % (query_path))

class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "RefrnaV2_7320",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "RefrnaV2_7320",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "ref_rna_v2.ref_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wsheet)

        result_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_rna_v2/annotation_result1/"

        new_annotation = result_dir + "new_annot"
        ref_annotation = result_dir + "ref_annot"
        merge_gene = result_dir + "gene_merge"
        merge_tran = result_dir + "tran_merge"
        gene2trans = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/AnnotMerge/output/newannot_class/all_tran2gene.txt"
        gene2trans_ref = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/AnnotMerge/output/refannot_class/all_tran2gene.txt"
        merge_annot = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_rna_v2/annotation_result2/merge_annot"

        merge_annot_v3 = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/AnnotMerge/output"
        gene_exp = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/Quant/gene.tpm.matrix"
        trans_exp = "/mnt/ilustre/users/sanger-dev/workspace/20190702/Refrna_tsg_34677/Quant/transcript.tpm.matrix"

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("ref_rna_v2.ref_annotation")
        params = {
            "nr_evalue": 1e-3,
            "nr_similarity": 0,
            "nr_identity": 0,
            "swissprot_evalue":1e-3,
            "swissprot_similarity": 0,
            "swissprot_identity": 0,
            "cog_evalue": 1e-3,
            "cog_similarity": 0,
            "cog_identity": 0,
            "kegg_evalue": 1e-3,
            "kegg_similarity": 0,
            "kegg_identity": 0,
            "pfam_evalue": 1e-3,
        }
        wf.test_api.species_name = "Mus_musculus"
        wf.test_api.run(merge_annot_v3, gene2trans, gene2trans_ref, params, version="v3", gene_exp = gene_exp, trans_exp=trans_exp)
        # wf.test_api.anno_type = 'latest'
        # wf.test_api.run(test_dir, trans2gene,  params)

if __name__ == '__main__':
    unittest.main()
