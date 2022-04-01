# -*- coding: utf-8 -*-
# __author__ = 'zzg'

from bson.objectid import ObjectId
import datetime
from mbio.api.database.whole_transcriptome.api_base import ApiBase
import types

class Common(ApiBase):
    def __init__(self, bind_object):
        super(Common, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_main_table(self, name, main_id = None,project_sn = None,task_id = None):
        dbname = name
        if main_id is None:
            dbname = name
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=  project_sn if project_sn else "123",
                task_id=task_id if project_sn else "1234" ,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc=name,
                params=[],
                status="start",
            )
            main_id = self.create_db_table(dbname, [main_info])
        else:
            dbname = name
            main_id = ObjectId(main_id)
            self.update_db_record(dbname, main_id, status="end", main_id=main_id)

    def update_main(self, collection, main_id, update):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！')

        self.db[collection].update({"_id": main_id}, {"$set":update})

    def add_synteny(self, name, dual_synteny_path,circle_path,dot_path,main_id = None,project_sn = None,task_id = None):
        dbname = name
        if main_id is None:
            dbname = name
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=  project_sn if project_sn else "123",
                task_id=task_id if project_sn else "1234" ,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc=name,
                params=[],
                status="start",
            )
            main_id = self.create_db_table(dbname, [main_info])
        else:
            dbname = name
            main_id = ObjectId(main_id)

        self.update_db_record(dbname, main_id, status="end", main_id=main_id,dual_synteny_plot=dual_synteny_path,circle_plot=circle_path,dot_plot=dot_path)