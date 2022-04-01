# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

from mbio.workflows.itraq_and_tmt.itraq_test_api import ItraqApiWorkflow
from biocluster.wsheet import Sheet
import sys
from biocluster.config import Config

task_id = 'tsg_32701'#sys.argv[1]
project_sn = '188_5bc7fe410c7ff'#sys.argv[2]
project_type = 'prok_rna'#sys.argv[3]
table = ['sg_structure_utr']#sys.argv[4]
result_dir = '/mnt/ilustre/users/sanger-dev/workspace/20181022/Prokrna_tsg_32701/Rockhopper/output'

data = {
    "id": task_id,
    #+ str(random.randint(1,10000)),
    #"id": "denovo_rna_v2",
    "project_sn": project_sn,
    #+ str(random.randint(1,10000)),
    "type": "workflow",
    "name": "prok_rna.prokrna_test_api",
    "options": {},
    }
wsheet = Sheet(data=data)
wf = ItraqApiWorkflow(wsheet)
wf.IMPORT_REPORT_DATA = True
wf.IMPORT_REPORT_AFTER_END = False

db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


def get_delete_targets():
    find_result = db['sg_table_relation'].find_one({})
    if find_result:
        target = find_result['target']
    else:
        raise Exception('{}没有"sg_table_relation"这张表'.format(self._project_type))
    return target

def remove_db_record(table_name, query_dict=None, **kwargs):
    """
    根据kwargs删除table_name中查询到的相关记录
    :param table_name: 要查询的表名，如sg_exp
    :param kwargs: 查询条件， 如 diff_id = ObjectId('xxx'), gene_id='ABC'
    :param quey_dict: dict info for querying
    :return:
    """
    conn = db[table_name]
    if query_dict:
        kwargs.update(query_dict)
    result = conn.find_one(kwargs)
    if result:
        conn.delete_many(kwargs)
        print('Success to delete records in {} by query {}'.format(table_name, kwargs))
    else:
        print('No record to delete from {} by query {}'.format(table_name, kwargs))

def remove_table_by_task_id(tuple_arg):
    """
    根据task_id删除指定的表，并且删除带有该记录_id的详情表
    :param main_table: 要删除的表的名称，如sg_exp
    :param task_id: task_id对应的值，是删除记录的条件
    :param detail_table: 主表关联的详情表名称如sg_exp_detail，如果有多个详情表，可以是列表,如[sg_xx_1, sg_xx_2]
    :param detail_table_key: 指示是那个字段(如exp_id)对应的值为主表记录的_id, 是删除记录的重要依据。可以是列表，与detail_table一一对应。
    :return:
    """
    main_table, detail_table, detail_table_key = tuple_arg
    all_collections = db.collection_names()
    if main_table not in list(all_collections):
        print("Warning: {} was not found in {}".format(main_table, self._project_type))
        return
    conn = db[main_table]
    a = conn.find({"task_id": task_id}, {'_id': 1})
    if type(a) == dict:
        a = [a]
    # delete main table record
    remove_num = 0
    for each in a:
        remove_num += 1
        remove_db_record(main_table, _id=each['_id'])
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
                remove_db_record(table, query_dict={table_key: each['_id']})
    if remove_num:
        print('Found {} main table(s) in {} by task_id {}.'.format(remove_num, main_table,task_id))
        print('And, finished to remove records of {} and {} by task_id {}'.format(main_table, detail_table, task_id))

target = get_delete_targets()
for turple in target:
    if turple[0] in table:
        remove_table_by_task_id(turple)

api = wf.api.api("prok_rna.gene_structure")
api.add_rock_structure(result_dir)