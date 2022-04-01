# -*- coding: utf-8 -*-
# __author__ = 'shijin'
from ..meta.update_status import UpdateStatus as US


class UpdateStatus(US):

    def __init__(self, data):
        super(UpdateStatus, self).__init__(data)
        self.mongodb = self._mongo_client["samplebase"]

    # 格式变更，解析错误
    # def update_status(self):
    #     status = self.data["content"]["stage"]["status"]
    #     desc = filter_error_info(self.data["content"]["stage"]["error"])
    #     if not self.update_info:
    #         return
    #     for obj_id, collection_name in json.loads(self.update_info).items():
    #         obj_id = ObjectId(obj_id)
    #         collection = self.mongodb[collection_name]
    #         if status == "start":
    #             insert_data = {
    #                 "status": 'start',
    #                 "desc": desc,
    #                 "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #             }
    #             collection.find_one_and_update({"table_id": obj_id}, {'$set': insert_data}, upsert=True)
    #         elif status == "finish":  # 只能有一次finish状态
    #             insert_data = {
    #                 "status": 'end',
    #                 "desc": desc,
    #                 "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #             }
    #             collection.find_one_and_update({"table_id": obj_id}, {'$set': insert_data}, upsert=True)
    #         else:
    #             insert_data = {
    #                 "status": status,
    #                 "desc": desc,
    #                 "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #             }
    #             collection.find_one_and_update({"table_id": obj_id}, {'$set': insert_data}, upsert=True)
