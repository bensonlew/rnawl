# -*- coding: utf-8 -*-
# __author__ = 'hongdong.xuan'
from bson.objectid import ObjectId
from biocluster.core.function import filter_error_info
from ..meta.update_status import UpdateStatus


class MedReportTupdate(UpdateStatus):

    def __init__(self, data):
        super(MedReportTupdate, self).__init__(data)
        self._client = "client03"
        # self._key = "hM4uZcGs9d"
        # self._url = "http://api.tsanger.com/task/add_file"
        # self.database = self._mongo_client[self.config.MONGODB + '_nipt']
        self._env_name = "offline"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._url = "http://apicenter.nsg.com/index/in"
        self._binds_id = "5f4324a19b79000093008059"
        self._interface_id = 66

    def update(self):
        pass

    def update_status(self):
        status = self.data["sync_task_log"]["task"]["status"]
        self.logger.info("status %s" % (status))
        desc = ''
        for i in self.data['sync_task_log']['log']:
            if 'name' not in i:
                desc = i['desc']
        desc = filter_error_info(desc)
        create_time = str(self.data["sync_task_log"]["task"]["created_ts"])

        if not self.update_info:
            return

        for obj_id, collection_name in self.update_info.items():
            obj_id = ObjectId(obj_id)
            collection = self.db[collection_name]
            if status != "start":
                data = {
                    "status": "end" if status == 'finish' else status,
                    "desc": desc,
                    "end_time": create_time
                }
                collection.update({"_id": obj_id}, {'$set': data}, upsert=True)
            # sg_status_col = self.database[collection_name]
            # if status == "start":
            #     tmp_col = self.database[collection_name]
            #     try:
            #         temp_find = tmp_col.find_one({"_id": obj_id})
            #         tb_name = temp_find["name"]
            #         temp_params = temp_find['params']
            #         submit_location = json.loads(temp_params)['submit_location']
            #     except:
            #         tb_name = ""
            #         temp_params = ''
            #         submit_location = ''
            #     insert_data = {
            #         "status": "start",
            #         "desc": desc,
            #         "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            #     }
            #     sg_status_col.update({"_id": obj_id}, {'$set': insert_data}, upsert=True)
            # elif status == "finish":  # 只能有一次finish状态
            #     insert_data = {
            #         "status": 'end',
            #         "desc": desc,
            #         "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            #     }
            #     sg_status_col.update({"_id": obj_id}, {'$set': insert_data}, upsert=True)
            # else:
            #     insert_data = {
            #         "status": status,
            #         "desc": desc,
            #         "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            #     }
            #     sg_status_col.update({"_id": obj_id}, {'$set': insert_data}, upsert=True)
            # self._mongo_client.close()
