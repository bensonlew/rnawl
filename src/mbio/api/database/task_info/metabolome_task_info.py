# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import datetime,json
from biocluster.api.database.base import Base, report_check
# from biocluster.config import Config
from bson import SON


class MetabolomeTaskInfo(Base):
    def __init__(self, bind_object):
        super(MetabolomeTaskInfo, self).__init__(bind_object)
        self._project_type = 'metabolome'

    #@report_check
    def add_task_info(self, db_name=None,diff_params=None, two_sample_diff_params=None):
        if self.db["sg_task"].find_one({"task_id":self.bind_object.sheet.id}):
            return
        if db_name:
            self._db_name = db_name
        json_data = [
            ('task_id', self.bind_object.sheet.id),
            ('member_id', self.bind_object.sheet.member_id),
            ('project_sn', self.bind_object.sheet.project_sn),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('member_type', self.bind_object.sheet.member_type),
            ('cmd_id', int(self.bind_object.sheet.cmd_id)),
            ('is_demo', 0),
            ('demo_id', self.bind_object.sheet.id),
            ('type',"")
        ]
        if diff_params:
            json_data.append(('diff_params',json.dumps(diff_params, sort_keys=True, separators=(',', ':'))))
        if two_sample_diff_params:
            json_data.append(('diff_params2',json.dumps(two_sample_diff_params, sort_keys=True, separators=(',', ':'))))
        self.db['sg_task'].insert_one(SON(json_data))
        self.bind_object.logger.info('任务信息导入sg_task成功。')




    def add_project_type(self, task_id, project_sn, project_type,hmdb="F",mix_table = None):
        project_type = project_type
        if mix_table:
            self.db['sg_task'].update({'task_id': task_id,'project_sn': project_sn,}, {'$set': {'type': project_type, 'hmdb':hmdb,
                                                                               "soft_series":"3.0" ,"version":"3.0","mix":mix_table}},multi=True)
        else:
            self.db['sg_task'].update({'task_id': task_id,'project_sn': project_sn,}, {'$set': {'type': project_type,'hmdb':hmdb,"soft_series":"3.0","version":"3.0"}},multi=True)

    def workflow_metabset_add_sg_status(self, collection_names,task_id,is_workflow_use=False):
        if is_workflow_use:
            data = []
            for collection_name in collection_names:
                temp_finds = self.db[collection_name].find({"task_id": task_id})
                for temp_find in temp_finds:
                    params = ""
                    submit_location = ""
                    temp_json = ''
                    if temp_find and "params" in temp_find.keys():
                        params = temp_find['params']
                        try:
                            temp_json = json.loads(temp_find['params'])
                        except:
                            pass
                        else:
                            if isinstance(temp_json, dict) and "submit_location" in temp_json.keys():
                                submit_location = temp_json['submit_location']

                    insert_data = {
                        "table_id": temp_find["_id"],
                        "table_name": temp_find["name"] if temp_find and "name" in temp_find.keys() else "",
                        "task_id": temp_find["task_id"],
                        "type_name": collection_name,
                        "params": params,
                        "submit_location": submit_location,
                        "status": "end",
                        "is_new": "new",
                        "desc": temp_find["desc"],
                        "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    }

                    #zouguanqing 20200304
                    if temp_json:
                        add_item = ['metab_set','metabset']
                        for each_i in add_item:
                            if each_i in temp_json.keys():
                                insert_data['metabset_id'] = temp_json[each_i]
                    insert_data['is_workflow'] = 'T'
                    data.append(insert_data)
       
            self.db["sg_status"].insert_many(data)

    def update_mongo(self,db_name,search, change):
        db = self.db[db_name]
        ret = db.find_one(search)
        if ret:
            db.update({"_id":ret["_id"]},{"$set":change})
            return str(ret["_id"])
        else:
            return 'not find'
    
