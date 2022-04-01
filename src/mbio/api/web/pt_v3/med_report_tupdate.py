# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# lasted modified by hongdong @ 20180205
import datetime
import traceback
import sys
from bson.objectid import ObjectId
from biocluster.core.function import filter_error_info
from ..meta.update_status import UpdateStatus


class MedReportTupdate(UpdateStatus):

    def __init__(self, data):
        super(MedReportTupdate, self).__init__(data)
        self._client = "client03"
        # self._key = "hM4uZcGs9d"
        # self._url = "http://api.tsg.com/task/add_file"
        self._project_type = 'pt_v3'
        self._binds_id = "5f4324a19b79000093008059"
        self._interface_id = 66
        self._env_name = "offline"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._url = "http://apicenter.nsg.com/index/in"
        
    def update(self):
        pass

    def _get_params(self, collection_name, main_id):
        try:
            return self.db[collection_name].find_one({"_id": main_id})
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()

    def update_status(self):
        if "status" in self.data["sync_task_log"]["task"].keys():
            status = self.data["sync_task_log"]["task"]["status"]
        else:
            return
        desc = ''
        for i in self.data['sync_task_log']['log']:
            if 'step_name' not in i:
                desc = i['desc']
        desc = filter_error_info(desc)
        # create_time = str(self.data["sync_task_log"]["task"]["created_ts"])

        if not self.update_info:
            return
        batch_id = self.update_info.pop("batch_id") if 'batch_id' in self.update_info else None
        for obj_id, collection_name in self.update_info.items():
            obj_id = ObjectId(obj_id)
            collection = self.db[collection_name]
            if str(collection_name) == 'sg_file_check':
                temp_find_ = self._get_params(collection_name, obj_id)
                if status in ["finish", "end"]:
                    data = {
                        "status": "end",
                        "desc": temp_find_['log'] if temp_find_ and "log" in temp_find_.keys() else desc,
                        "end_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "log": ""
                    }
                    collection.update({"_id": obj_id}, {'$set': data}, upsert=True, multi=True)
                elif status == 'failed':
                    # self.logger.info(desc)
                    data = {
                        "status": "failed",
                        "desc": temp_find_['log'] if temp_find_ and desc == 'FileCheck.FileCheckV4 error' else desc,
                        "end_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        'log': ''
                    }
                    collection.update({"_id": obj_id}, {'$set': data}, upsert=True, multi=True)
                else:
                    return
            else:
                if status != "start":
                    data = {
                        "status": "end" if status == 'finish' else status,
                        "desc": desc,
                        "end_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    }
                    collection.update({"_id": obj_id}, {'$set': data}, upsert=True, multi=True)
                    self.family_update(batch_id)
                if str(collection_name) == 'sg_datasplit':
                    sg_status_col = self.db['sg_datasplit_records']
                    temp_find = self._get_params(collection_name, obj_id)
                    if status in ["finish", "end"]:  # 只能有一次finish状态
                        insert_data = {
                            "status": "end",
                            "desc": "运行结束",
                            "batch_id": obj_id,
                            "board_batch": temp_find['board_batch'] if temp_find and "board_batch" in temp_find.keys()
                            else "--",
                            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            "member_id": temp_find['member_id'] if temp_find and "member_id" in temp_find.keys() else "--"
                        }
                        sg_status_col.insert_one(insert_data)
                    elif status == "failed":
                        insert_data = {
                            "status": "failed",
                            "desc": desc,
                            "batch_id": obj_id,
                            "board_batch": temp_find['board_batch'] if temp_find and "board_batch" in temp_find.keys()
                            else "--",
                            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            "member_id": temp_find['member_id'] if temp_find and "member_id" in temp_find.keys() else "--"
                        }
                        sg_status_col.insert_one(insert_data)
                else:
                    return

    def family_update(self, batch_id):
        """
        用于更新家系分析的状态表-- add by hd 20180303 修复没有FM call snp的时候 默认插入进度0/0
        :param batch_id:
        :return:
        """
        if not batch_id:
            return
        batch_id = ObjectId(batch_id)
        collection = self.db['sg_analysis_status']
        try:
            result_ = collection.find({'batch_id': batch_id, 'type': "pt", "is_show": "1"}).sort([("created_ts", 1)])
            if result_.count() == 0:
                collection.insert_one({
                    "member_id": "656",
                    "batch_id": batch_id,
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "desc": "任务结束",
                    "status": "end",
                    "type": "pt",
                    "snp_all_counts": 0,
                    "snp_end_counts": 0,
                    "is_show": "1",
                    "family_end_counts": 0,
                    "family_all_counts": 0,
                    "end_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                })
            family_all_counts = self.db['sg_father'].find({"batch_id": batch_id}).count()
            family_end_counts = self.db['sg_father'].find({"batch_id": batch_id,
                                                           'status': {'$in': ['end', 'failed']}}).count()
            collection.update({'_id': result_[0]['_id']}, {'$set': {'family_end_counts': family_end_counts,
                                                                    "family_all_counts": family_all_counts}},
                              upsert=True, multi=True)
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            print e
            sys.stdout.flush()



