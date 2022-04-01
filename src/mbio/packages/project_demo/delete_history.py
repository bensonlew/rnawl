# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, shicaiping'
## 用于清除历史项目的部分详情表，保留主表及图形数据

import logging
import sys
import time
from biocluster.config import Config
from concurrent.futures import ThreadPoolExecutor

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


class DeleteHistoryMongo(object):
    def __init__(self, task_id, project_type, delete_targets=None):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self.db_static = Config().get_mongo_client(mtype=project_type, dydb_forbid=True)[Config().get_mongo_dbname(project_type, dydb_forbid=True)]
        self.task_id = task_id
        self._project_type = project_type
        self.project2table = {'ref_rna_v2': ['sg_snp_detail', 'sg_annotation_cog_detail', 'sg_annotation_go_detail',
                                             'sg_annotation_go_graph', 'sg_annotation_go_level',
                                             'sg_annotation_go_list', 'sg_annotation_kegg_categories',
                                             'sg_annotation_kegg_level', 'sg_annotation_kegg_pic',
                                             'sg_annotation_kegg_table', 'sg_annotation_nr_detail',
                                             'sg_annotation_pfam_bar', 'sg_annotation_pfam_detail',
                                             'sg_annotation_query_detail', 'sg_annotation_swissprot_detail',
                                             'sg_diff_detail', 'sg_diff_summary', 'sg_exp_corrsf_detail',
                                             'sg_exp_detail', 'sg_genes_detail', 'sg_species_information_detail',
                                             'sg_splicing_rmats_detail', 'sg_splicing_rmats_diffcomp_detail',
                                             'sg_tf_predict_detail', 'sg_wgcna_module_eigengenes_detail',
                                             'sg_wgcna_module_membership_detail', 'sg_wgcna_module_tree_detail'
                                             ]}
        if not delete_targets:
            self.delete_targets = self.get_delete_targets()
        else:
            self.delete_targets = delete_targets

    def get_delete_targets(self):
        if self._project_type == 'whole_transcriptome':
            document = self.db_static['table_relation'].find_one({})
            if document and 'target' in document:
                return self.db_static['table_relation'].find_one({})['target']
        find_result = self.db_static['sg_table_relation'].find_one({})
        if find_result:
            target = find_result['target']
        else:
            raise Exception('can not find sg_table_relation with project type ({})'.format(self._project_type))
        return target

    def run(self):
        target_collections = self.delete_targets
        if not target_collections:
            logging.info('targets for deleting is empty, abord')
            return
        with ThreadPoolExecutor(8) as pool:
            pool.map(self.remove_table_by_task_id, target_collections)
        #self.deal_sg_status()

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
        all_collections = self.db.collection_names()
        if main_table not in list(all_collections):
            logging.warn('{} was not found with project type ({})'.format(main_table, self._project_type))
            return
        collection = self.db[main_table]
        cursor = collection.find({'task_id': self.task_id}, {'_id': 1})
        if type(cursor) == dict:
            cursor = [cursor]
        remove_num = 0
        delete_table = list()
        for main_document in cursor:
            if detail_table:
                if not type(detail_table) == list:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if not type(detail_table_key) == list:
                        detail_table_key = [detail_table_key] * len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    if table in self.project2table[project_type]:
                        delete_table.append(table)
                        if not table_key:
                            raise Exception('you should specify detail_table_key whose value is main table "_id"')
                        self.remove_db_record(table, query_dict={table_key: main_document['_id']})
                        remove_num += 1
        if remove_num:
            logging.info(
                'found {} main table(s) in {} by task_id ({}), removed records of {}'.format(remove_num, main_table,
                                                                                             self.task_id,
                                                                                             delete_table))

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        """
        根据kwargs删除table_name中查询到的相关记录
        :param table_name: 要查询的表名，如sg_exp
        :param kwargs: 查询条件， 如 diff_id = ObjectId('xxx'), gene_id='ABC'
        :param quey_dict: dict info for querying
        :return:
        """
        collection = self.db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        document = collection.find_one(kwargs)
        if document:
            collection.delete_many(kwargs)
            logging.info('succeed in deleting records in {} by query {}'.format(table_name, kwargs))
        else:
            logging.info('find no record to delete from {} by query {}'.format(table_name, kwargs))

    def deal_sg_status(self):
        collection = self.db['sg_status']
        cursor = collection.find({"task_id": self.task_id}, {'_id': 1})
        if type(cursor) == dict:
            cursor = [cursor]
        for document in cursor:
            collection.update({'_id': document['_id']}, {'$set': {'status': 'deleted'}})


if __name__ == '__main__':
    logging.info('Usage: python delete_history.py <task_id> <project_type>')
    task_id = sys.argv[1]
    project_type = sys.argv[2]
    start_time = time.time()
    inst = DeleteHistoryMongo(task_id, project_type, delete_targets=None)
    inst.run()
    end_time = time.time()
    logging.info('Elapsed time: {}'.format(end_time - start_time))
