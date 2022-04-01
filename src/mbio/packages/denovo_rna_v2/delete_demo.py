# -*- coding: utf-8 -*-
# __author__ = 'gdq'
import time
from biocluster.api.database.base import Base
import gevent.pool
import gevent.monkey

#gevent.monkey.patch_all()


class DeleteDemoMongo(Base):
    def __init__(self, task_id, project_type='denovo_rna_v2', delete_targets=None):
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
        if type(a) == dict():
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
        [(u'sg_specimen', None, None), # yes
         (u'sg_assembly', [u'sg_assembly_len', 'sg_assembly_stat'], u'assembly_id'),  # yes
         (u'sg_align', u'sg_align_detail', u'align_id'),  # yes
         (u'sg_annotation_stat', u'sg_annotation_stat_detail', u'stat_id'),  # yes
         (u'sg_annotation_nr', [u'sg_annotation_nr_detail', 'sg_annotation_nr_pie'], u'nr_id'),  # yes
         (u'sg_annotation_swissprot', [u'sg_annotation_swissprot_detail', 'sg_annotation_swissprot_pie'], u'swissprot_id'),  # yes
         (u'sg_annotation_pfam', [u'sg_annotation_pfam_detail', 'sg_annotation_pfam_bar'], u'pfam_id'),  # yes
         (u'sg_annotation_query', u'sg_annotation_query_detail', u'query_id'),  # yes
         (u'sg_geneset_cog_class', u'sg_annotation_cog_detail', u'cog_id'),  # yes
         (u'sg_annotation_go', [u'sg_annotation_go_detail', 'sg_annotation_go_graph','sg_annotation_go_level', 'sg_annotation_go_list'], u'go_id'),  # yes
         (u'sg_geneset_kegg_class', [u'sg_annotation_kegg_categories', 'sg_annotation_kegg_level','sg_annotation_kegg_table'], u'kegg_id'),  # yes
         (u'sg_snp', [u'sg_snp_detail', u'sg_snp_dp', u'sg_snp_hh', u'sg_snp_tt'], u'snp_id'),  # yes
         (u'sg_geneset_venn', None, None),  # yes
         (u'sg_geneset_kegg_class', [u'sg_geneset_kegg_class_detail', u'sg_geneset_kegg_class_pathway', u'sg_geneset_kegg_class_statistic'], u'kegg_id'),  # yes
         (u'sg_specimen_group', None, None), # yes
         (u'sg_geneset', u'sg_geneset_detail', u'geneset_id'),  # yes
         (u'sg_exp_venn', u'sg_exp_venn_detail', u'venn_id'),  # yes
         (u'sg_exp_corr', u'sg_exp_corr_detail', u'corr_id'),  # yes
         (u'sg_cds', None, None),  # yes
         (u'sg_geneset_cog_class', u'sg_geneset_cog_class_detail', u'geneset_cog_id'),  # yes
         (u'sg_status', None, None),  # yes
         (u'sg_ssr', [u'sg_ssr_detail', u'sg_ssr_class'], u'ssr_id'),  # yes
         (u'sg_exp_graph', [u'sg_exp_graph_box', u'sg_exp_graph_density', u'sg_exp_graph_volin'], u'graph_id'),  # yes
         (u'sg_specimen_group_compare', None, None), # yes
         (u'sg_geneset_cluster', [u'sg_geneset_cluster_detail', u'sg_geneset_cluster_tree'], u'cluster_id'),  # yes
         (u'sg_task', None, None),  # yes
         (u'sg_tf', u'sg_tf_detail', u'tf_id'),  # yes
         (u'sg_geneset_kegg_enrich', u'sg_geneset_kegg_enrich_detail', u'kegg_enrich_id'),  # yes
         (u'sg_geneset_go_enrich', u'sg_geneset_go_enrich_detail', u'go_enrich_id'),  # yes
         (u'sg_specimen_graphic', None, None), # yes
         (u'sg_diff', [u'sg_diff_detail', u'sg_diff_summary', u'sg_diff_volcano'], u'diff_id'),  # yes
         (u'sg_exp', u'sg_exp_detail', u'exp_id'),  # yes
         (u'sg_result_table_deposit', None, None),  # yes
         (u'sg_geneset_go_class', u'sg_geneset_go_class_detail', u'go_regulate_id'),  # yes
         (u'sg_exp_pca', u'sg_exp_pca_detail', u'pca_id'),  # yes
         (u'sg_t2g', u'sg_t2g_detail', u't2g_id'),  # yes
         (u'sg_t2g_detail', None, None)]  # yes
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

