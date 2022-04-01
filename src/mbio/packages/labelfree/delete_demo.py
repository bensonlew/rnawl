# -*- coding: utf-8 -*-
# __author__ = 'gdq'
import time
from biocluster.api.database.base import Base
import gevent.pool
import gevent.monkey

gevent.monkey.patch_all()


class DeleteDemoMongo(Base):
    def __init__(self, task_id, project_type='labelfree', delete_targets=None):
        super(DeleteDemoMongo, self).__init__()
        self._project_type = project_type
        self.task_id = task_id
        self.delete_targets = delete_targets

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        """
        根据kwargs删除table_name中查询到的相关记录
        :param table_name: 要查询的表名，如sg_exp
        :param kwargs: 查询条件， 如 diff_id = ObjectId('xxx'), gene_id='ABC'
        :param quey_dict: dict info for querying
        :return:
        """
        conn = self.db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        result = conn.find_one(kwargs)
        if result:
            conn.delete_many(kwargs)
        else:
            print('No record to delete from {} by {}'.format(table_name, kwargs))

    def remove_table_by_task_id(self, tuple_arg):
        """
        根据task_id删除指定的表，并且删除带有该记录_id的详情表
        :param main_table: 要删除的表的名称，如sg_exp
        :param task_id: task_id对应的值，是删除记录的条件
        :param detail_table: 主表关联的详情表名称如sg_exp_detail，如果有多个详情表，可以是列表,如[sg_xx_1, sg_xx_2]
        :param detail_table_key: 指示是那个字段(如exp_id)对应的值为主表记录的_id, 是删除记录的重要依据。可以是列表，与detail_table一一对应。
        :return:
        """
        main_table, detail_table, detail_table_key = tuple_arg
        conn = self.db[main_table]
        a = conn.find({"task_id": self.task_id}, {'_id': 1})
        if type(a) == dict:
            a = [a]
        # delete main table record
        for each in a:
            self.remove_db_record(main_table, _id=each['_id'])
            if detail_table:
                if not type(detail_table) == list:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if not type(detail_table_key) == list:
                        detail_table_key = [detail_table_key]*len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    if not table_key:
                        raise Exception('you should specify detail_table_key whose value is main table "_id"')
                    self.remove_db_record(table, query_dict={table_key: each['_id']})
        print('Have removed records of {} and {} by {}'.format(main_table, detail_table, self.task_id))

    def get_target_collection_info(self):
        delete_target = \
        [(u'sg_specimen', None, None),  # yes
         (u'sg_peptide_error', None, None),  # yes
         (u'sg_peptide_len', None, None),  # yes
         (u'sg_peptide_num', None, None),  # yes
         (u'sg_protein_coverage', None, None),  # yes
         (u'sg_protein_info', None, None),  # yes
         (u'sg_protein_mw', None, None),        # yes
         (u'sg_protein_mw', None, None),        # yes
         (u'sg_proteinset_info', None, None),        # yes
         (u'sg_query_seq', None, None),        # yes
         (u'sg_software', None, None),        # yes
         (u'sg_annotation_stat', u'sg_annotation_stat_detail', u'stat_id'),  # yes
         (u'sg_annotation_pfam', [u'sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], u'pfam_id'),  # yes
         (u'sg_annotation_query', u'sg_annotation_query_detail', u'query_id'),  # yes
         (u'sg_annotation_go', [u'sg_annotation_go_detail', 'sg_annotation_go_graph', 'sg_annotation_go_level', 'sg_annotation_go_list'], u'go_id'),
         (u'sg_annotation_cog', u'sg_annotation_cog_detail', u'cog_id'),  # yes
         (u'sg_annotation_subloc', [u'sg_annotation__subloc_bar', u'sg_annotation__subloc_detail'], u'subloc_id'),  # yes
         (u'sg_annotation_kegg', [u'sg_annotation_kegg_categories', 'sg_annotation_kegg_level', 'sg_annotation_kegg_table', u'sg_annotation_kegg_important'], u'kegg_id'),  # yes
         (u'sg_proteinset_venn', None, None),  # yes
         (u'sg_proteinset_kegg_class', [u'sg_proteinset_kegg_class_detail', u'sg_proteinset_kegg_class_pathway', u'sg_proteinset_kegg_class_statistic'], u'kegg_id'),  # yes
         (u'sg_specimen_group', None, None), # yes
         (u'sg_proteinset', u'sg_proteinset_detail', u'proteinset_id'),  # yes
         (u'sg_search_db', u'sg_search_db_detail', u'searchdb_id'),  # yes
         (u'sg_proteinset_ipath', u'sg_proteinset_ipath_detail', u'ipath_id'),  # yes
         (u'sg_exp_venn', u'sg_exp_venn_detail', u'venn_id'),  # yes
         (u'sg_express_corr_', u'sg_exp_corr_detail', u'corr_id'),  # yes
         (u'sg_proteinset_cog_class', u'sg_proteinset_cog_class_detail', u'proteinset_cog_id'),  # yes
         (u'sg_status', None, None),  # yes
         (u'sg_specimen_group_compare', None, None), # yes
         (u'sg_proteinset_cluster', [u'sg_proteinset_cluster_detail', u'sg_proteinset_cluster_tree'], u'cluster_id'),  # yes
         (u'sg_proteinset_circ', [u'sg_proteinset_circ_detail', u'sg_proteinset_circ_graph'], u'circ_id'),  # yes
         (u'sg_task', None, None),  # yes
         (u'sg_proteinset_kegg_enrich', u'sg_proteinset_kegg_enrich_detail', u'kegg_enrich_id'),  # yes
         (u'sg_proteinset_go_enrich', u'sg_proteinset_go_enrich_detail', u'go_enrich_id'),  # yes
         (u'sg_diff', [u'sg_diff_detail', u'sg_diff_summary'], u'diff_id'),  # yes
         (u'sg_express', u'sg_express_detail', u'express_id'),  # yes
         (u'sg_result_table_deposit', None, None),  # yes
         (u'sg_proteinset_go_class', u'sg_proteinset_go_class_detail', u'go_regulate_id'),  # yes
         (u'sg_proteinset_go_class2', u'sg_proteinset_go_class2_detail', u'go_regulate_id'),  # yes
         (u'sg_express_pca', [u'sg_express_pca_gene_detail', u'sg_express_pca_sample_detail'], u'pca_id'),# yes
         (u'sg_proteinset_string_picture', [u"sg_proteinset_string_picture_annotation", u"sg_proteinset_string_picture_bitscore",u"sg_proteinset_string_picture_interaction",], u'string_id'), #yes
         (u"sg_proteinset_pfam", u"sg_proteinset_pfam_stat", u"pfam_id"), # yes
         (u"sg_proteinset_subloc", u"sg_proteinset_subloc_stat", u"subloc_id"), # yes
         ]
        return delete_target

    def run(self):
        if not self.delete_targets:
            target_collections = self.get_target_collection_info()
        else:
            target_collections = self.delete_targets
        if not target_collections:
            return
        print('we will delete the following tables:', target_collections)
        pool = gevent.pool.Pool(100)
        pool.map(self.remove_table_by_task_id, target_collections)


if __name__ == '__main__':
    import sys
    task_id = sys.argv[1]
    start_time = time.time()
    # task_id = 'copy_demo_test'
    del_demo = DeleteDemoMongo(task_id)
    del_demo.run()
    end_time = time.time()
    print("total time: {}".format(end_time - start_time))

