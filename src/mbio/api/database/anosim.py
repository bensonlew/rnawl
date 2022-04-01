# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# last modified guhaidong 20171115
# import os
import json
import datetime
import re
from biocluster.api.database.base import Base, report_check
from types import StringTypes
from bson.objectid import ObjectId
# from biocluster.config import Config
# import re
# import datetime
# from bson.son import SON
# from types import StringTypes


class Anosim(Base):
    def __init__(self, bind_object):
        super(Anosim, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def add_beta_anosim_result(self, dir_path, main_id=None, main=False, group_id=None,
                               task_id=None, otu_id=None, name=None, params=None, level=9):

        if task_id == None:
            task_id = self.bind_object.sheet.id
            task_id = "_".join(task_id.split('_')[0:2])  # get task_id by guhaidong 20171115
        def insert_table_detail(file_path, table_type, update_id,
                                coll_name='sg_beta_multi_anosim_detail',
                                main_coll='sg_beta_multi_anosim',
                                update_column=True, db=self.db, comment='', stats=False,
                                columns=None):
            collection = db[coll_name]
            data_list = []
            with open(file_path, 'rb') as f:
                all_lines = f.readlines()
                if comment:
                    flag = 0
                    for line in all_lines:
                        if len(line) > 0:
                            if line[0] == '#':
                                flag += 1
                            else:
                                break
                        else:
                            self.bind_object.set_error('表中存在空的行', code="51000101")
                    all_lines = all_lines[flag:]
                if isinstance(columns, list):
                    pass
                else:
                    columns = all_lines[0].rstrip().split('\t')[1:]
                for line in all_lines[1:]:
                    values = line.rstrip().split('\t')

                    insert_data = {
                        'anosim_id': update_id,
                        'type': table_type
                        # 'name': values[0]
                        }
                    if stats:
                        insert_data['group1'] = values[0]
                        for m, n in enumerate(values):
                            if re.match(r'^\[[ \.0-9]+\]$', n):
                                values[m] = n.strip('[ ]')
                    else:
                        insert_data['name'] = values[0]
                    values_dict = dict(zip(columns, values[1:]))
                    insert_data = dict(insert_data, **values_dict)
                    data_list.append(insert_data)
                    # collection.insert_one(insert_data)
                try:
                    collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.logger.error("导入%s数据库出错:%s" % (coll_name, e))
                else:
                    self.bind_object.logger.info("导入%s数据库成功!" % coll_name)
                if update_column:
                    main_collection = db[main_coll]
                    default_column = {'specimen': 'detail_column', 'factor': 'factor_column', 'vector': 'vector_column',
                                      'species': 'species_column', 'rotation': 'rotation_column', 'box': 'box_column',
                                      'anosim': 'anosim_column', 'stats': 'stats_column'}
                    if table_type in default_column:
                        main_collection.update_one({'_id': update_id, 'task_id': task_id},
                                                   {'$set': {default_column[table_type]: ','.join(columns)}},
                                                   upsert=False)  # add 'task_id' by guhaidong 20171115

        def insert_text_detail(file_path, data_type, main_id,
                               coll_name='sg_beta_multi_anosim_json_detail', db=self.db):
            collection = db[coll_name]
            with open(file_path, 'rb') as f:
                data = f.read()
                insert_data = {
                    'multi_analysis_id': main_id,
                    'type': data_type,
                    'json_value': data
                    }
                collection.insert_one(insert_data)
        _main_collection = self.db['sg_beta_multi_anosim']

        if main:
            if level and level not in range(1, 10):
                self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
                self.bind_object.set_error("level水平错误", code="51000102")
            #　if task_id is None:
            #　    task_id = self.bind_object.sheet.id  # remove by guhaidong 20171115
            if isinstance(group_id, ObjectId):
                pass
            elif group_id is not None:
                group_id = ObjectId(group_id)
            else:
                group_id = ObjectId(self.bind_object.option('group_id'))
            if not isinstance(otu_id, ObjectId) and otu_id is not None:
                otu_id = ObjectId(otu_id)
            if 'params' in self.bind_object.sheet._data:
                params = self.bind_object.sheet.option('params')
            insert_mongo_json = {
                'project_sn': self.bind_object.sheet.project_sn,
                'task_id': task_id,
                'otu_id': otu_id,
                # 'level_id': int(level),
                'name': 'Anosim&Adonis_Origin',
                'group_id': group_id,
                'params': (json.dumps(params, sort_keys=True, separators=(',', ':'))
                           if isinstance(params, dict) else params),
                'status': 'end',
                'desc': '',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                }
            anosim_id = _main_collection.insert_one(insert_mongo_json).inserted_id
            main_id = anosim_id
            #_main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id}})
        else:
            if not main_id:
                self.bind_object.set_error('不写入主表时，需要提供主表ID', code="51000103")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        result = _main_collection.find_one({'_id': main_id, 'task_id': task_id})  # add task_id by guhaidong   20171115
        if result:
            anosim_path = dir_path.rstrip('/') + '/Anosim/format_results.xls'
            adonis_path = dir_path.rstrip('/') + '/Anosim/adonis_results.txt' # add by qingchen.zhang 20190528
            box_path = dir_path.rstrip('/') + '/AnosimBox/box_data.xls'
            # stats_path = dir_path.rstrip('/') + '/Box/Stats.xls'
            #_main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id}})
            insert_table_detail(anosim_path, 'anosim', update_id=main_id, update_column=False,
                                columns=['statistic', 'pvalue', 'permutation_number'])
            insert_table_detail(adonis_path, 'adonis', update_id=main_id, update_column=False,
                                columns=['Df', 'Sums_Of_Sqs', 'Mean_Sqs',  'F_Model', 'R2', 'Pr'])
            insert_table_detail(box_path, 'box', update_id=main_id, update_column=False)
            # insert_table_detail(stats_path, 'stats', update_id=main_id, comment='#', stats=True,
            #                     update_column=False, columns=['group2', 't_statistic',
            #                                                   'param_pvalue', 'param_correct_pvalue',
            #                                                   'nonparam_pvalue', 'nonparam_correct_pvalue'])
        else:
            self.bind_object.logger.error('提供的_id：%s在sg_beta_multi_anosim中无法找到表, taskid: %s' % (str(main_id), task_id))
            self.bind_object.set_error("sg_beta_multi_anosim表找不到数据", code="51000104")
        if main:
            return main_id

    @report_check
    def add_beta_anosim_main(self, name, otu_id, params=None):
        if otu_id != 0 and not isinstance(otu_id, ObjectId):
            if isinstance(otu_id, StringTypes):
                otu_id = ObjectId(otu_id)
            else:
                self.bind_object.set_error("otu_id必须为ObjectId对象或其对应的字符串!", code="51000712")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": otu_id})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "anosim分析"
        insert_data = {
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "desc": desc,
            "name": name if name else "anosim",
            "otu_id": otu_id,
            "params": params,
            "project_sn": project_sn,
            "status": "end",
            "task_id": task_id
        }
        collection = self.db["sg_beta_multi_anosim"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        # collection.update({"_id": inserted_id}, {"$set": {"main_id": inserted_id}})
        return inserted_id
