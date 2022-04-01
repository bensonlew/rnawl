# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from __future__ import print_function
import time
import os
import json
import shutil
from bson import ObjectId
from biocluster.config import Config


def linkfile(from_path,to_path):
    if os.path.exists(to_path):
        os.remove(to_path)
    os.link(from_path,to_path)

def linkdir(olddir, newdir, mode='link'):
    """
    移动一个目录到另一个目录
    ; olddir 需要移动的目录参数
    ；newdir 需要移动的目的位置
    """
    if not os.path.isdir(olddir):
        raise Exception('需要移动到output目录的文件夹不存在{}'.format(olddir))
    if os.path.exists(newdir):
        shutil.rmtree(newdir)
    os.mkdir(newdir)
    allfiles = os.listdir(olddir)
    oldfiles = [os.path.join(olddir, i) for i in allfiles]
    newfiles = [os.path.join(newdir, i) for i in allfiles]
    for i in range(len(allfiles)):
        if os.path.isfile(oldfiles[i]):
            os.link(oldfiles[i], newfiles[i])
        else:
            linkdir(oldfiles[i], newfiles[i])


class InteractionDelete(object):
    def __init__(self,bind_object, project_type, table=None, main_id=None):
        self.bind_object = bind_object
        self.project_type = project_type
        print("项目是{}".format(self.project_type))
        self.db_static = Config().get_mongo_client(mtype=self.project_type, dydb_forbid=True)[Config().get_mongo_dbname(self.project_type, dydb_forbid=True)]
        self.db =  Config().get_mongo_client(mtype=self.project_type)[Config().get_mongo_dbname(self.project_type)]
        self.main_id = main_id
        self.status_table = self.get_interaction_status_info()
        self.delete_targets = self.get_delete_targets()



    def get_interaction_status_info(self):
        sg_status_col = self.db['sg_status']
        target_record = sg_status_col.find_one({"table_id":ObjectId(self.main_id) })
        return target_record

    def delete_interactions_records(self):
        collection_table = self.status_table["type_name"]
        self.remove_table_by_main_id(collection_table,self.main_id)



    def get_delete_targets(self):
        if self.project_type == 'whole_transcriptome':
            document = self.db_static['table_relation'].find_one({})
            if document and 'target' in document:
                target = self.db_static['table_relation'].find_one({})['target']
        else:
            find_result = self.db_static['sg_table_relation'].find_one({})
            if find_result:
                target = find_result['target']
            else:
                self.bind_object.set_error('can not find sg_table_relation with project type ({})'.format(self.project_type))
        table_infos = {}
        for i in target:
            main_table_name, releted_tables, related_id = i[0], i[1], i[2]
            if main_table_name in table_infos:
                table_infos[main_table_name].append([releted_tables, related_id])
            else:
                table_infos[main_table_name] = [[releted_tables, related_id]]
        return table_infos

    def remove_table_by_main_id(self,collection_table,main_id):
        # main_table, detail_table, detail_table_key = self.delete_targets
        all_collections = self.db.collection_names()
        if collection_table not in all_collections:
            self.bind_object.set_error('{} was not found with project type ({})'.format(collection_table, self.project_type))
        deleted_infos = self.delete_targets[collection_table]
        for deleted_info in deleted_infos:
            main_table = collection_table
            detail_table, detail_table_key = deleted_info
            # self.remove_db_record(main_table, _id=ObjectId(main_id))
            if detail_table:
                if not type(detail_table) == list:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if not type(detail_table_key) == list:
                        detail_table_key = [detail_table_key] * len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    if not table_key:
                        raise Exception('you should specify detail_table_key whose value is main table "_id"')
                    self.remove_db_record(table, query_dict={table_key: ObjectId(main_id)})


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
            self.bind_object.logger.info('succeed in deleting records in {} by query {}'.format(table_name, kwargs))
        else:
            self.bind_object.logger.info('find no record to delete from {} by query {}'.format(table_name, kwargs))






