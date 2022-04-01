# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modified = guhaidong 20171122

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from bson import ObjectId
# from biocluster.config import Config


class EnvMetagenomic(Base):
    def __init__(self, bind_object=None):
        super(EnvMetagenomic, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.env_lab = ""
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_env_table(self, file_path, name_to_id, task_id=None, name=None,
                      is_used="0", env_true={}):  # add is_used 20171122
        if isinstance(name_to_id, types.DictType):
            pass
        else:
            self.bind_object.set_error('错误的name_to_id类型，或者没有提供', code="52800901")
        if not task_id:
            task_id = self.bind_object.sheet.id
        names, details = self._get_table_info(file_path, name_to_id)
        self.env_lab = ','.join(names)
        main_insert_data = {
            'project_sn': self.bind_object.sheet.project_sn,
            'task_id': task_id,
            'name': name if name else 'origin_env_table',
            'env_names': ','.join(names),
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'is_used': is_used
        }  # add is_used 20171122
        if env_true:
            main_insert_data['env_true_names'] = [env_true[n] for n in names]
        main_collection = self.db['env']
        collections = self.db['env_detail']
        data_list = []
        main_env_id = main_collection.insert_one(main_insert_data).inserted_id
        for sp_id, sp_dict in details.items():  # "details" is a dict, and "details"'s values is dict too
            sp_dict['env_id'] = main_env_id
            if not isinstance(sp_id, ObjectId):
                sp_id = ObjectId(sp_id)
            sp_dict['specimen_id'] = sp_id
            data_list.append(sp_dict)
            # collections.insert_one(sp_dict).inserted_id
        collections.insert_many(data_list)
        return main_env_id

    def _get_table_info(self, file_path, name_to_id):
        with open(file_path, 'rb') as f:
            lines = f.readlines()
            env_names = lines[0].rstrip().split('\t')[1:]
            info_dic = {}
            for line in lines[1:]:
                values = line.rstrip().split('\t')
                if values[0] not in name_to_id:
                    self.bind_object.set_error('分组文件或环境因子表中的样本：%s没有在fastq文件中找到',
                                               variables=(values[0]), code="52800902")
                info_dic[name_to_id[values[0]]] = dict(zip(env_names, values[1:]))
            return env_names, info_dic

    def get_env_lab(self):  #  add env_labs by guhaidong 20171114
        return self.env_lab
