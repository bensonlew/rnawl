# -*- coding: utf-8 -*-
# __author__ = 'gdq'
import time
from biocluster.config import Config
from concurrent.futures import ThreadPoolExecutor


class MongoIndex(object):
    def __init__(self, project_type):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self._project_type = project_type
        self.targets = self.get_targets()

    def get_targets(self):
        find_result = self.db['sg_table_relation'].find_one({})
        if find_result:
            target = find_result['target']
        else:
            raise Exception('{}没有"sg_table_relation"这张表'.format(self._project_type))
        return target

    def deal_target(self, tuple_arg):
        main_table, detail_table, detail_table_key = tuple_arg
        index_info = ''
        if not detail_table:
            index_info += main_table + '\ttask_id' + '\t主表不进行片键\n'
        else:
            index_info += main_table + '\ttask_id' + '\t主表不进行片键\n'
            if not type(detail_table) == list:
                detail_table = [detail_table]
            for de in detail_table:
                index_info += de + '\t' + detail_table_key + '\t{' + detail_table_key + ':1}\n'
        return index_info


if __name__ == '__main__':
    import sys

    print("Usage: python get_project_mongo_index.py <project_type>")
    project_type = sys.argv[1]
    start_time = time.time()
    index = MongoIndex(project_type)
    with open(project_type + '_mongo_index.txt', 'w') as index_w:
        index_w.write('#collection_name\tindex\tpianjian\n')
        for i in index.targets:
            index_w.write(index.deal_target(i))
    end_time = time.time()
    print("total time: {}".format(end_time - start_time))


