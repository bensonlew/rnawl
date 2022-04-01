# -*- coding: utf-8 -*-
import web
import json
import os
import datetime
from bson import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.core.base import Base
from mainapp.controllers.core.basic import Basic
from biocluster.core.function import filter_error_info


class DeleteTask(Base):
    def __init__(self):
        """
        1.用于页面删除同时删除mongo数据；
        2.用于页面删除task或者批量删除，同时删除mongo数据；
        :return:
        """
        super(DeleteTask, self).__init__()
        self.input_data = web.input()
        self._project_type = self.input_data.type

    def POST(self):
        requires = ['type', 'task_id']#分别是指项目类型，任务id
        for i in requires:
            if not hasattr(self.input_data, i):
                return json.dumps({"success":False, "info": "Lack params: %s" %(i)})
        workflow_id = "Delete_" + "{}_{}".format(self.input_data.task_id, datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        task_id = ",".join(self.input_data.task_id)

        data_json = {
            'id': workflow_id,
            'stage_id':0,
            'name': "project_demo.delete_task",
            'type': "workflow",
            'client': self.input_data.client,
            "IMPORT_REPORT_DATA": False,
            "IMPORT_REPORT_AFTER_END": False,
            'options': {
                "task_id": task_id,
                "project_type": self.input_data.type
            }
        }
        if hasattr(self.input_data, "location"):#需要删除交互分析还是工作流
            data_json["location"] = self.input_data.location
        if hasattr(self.input_data, "main_id"):#
            data_json["location"] = self.input_data.location

        workflow_client =Basic(data=data_json, instant=False)
        try:
            run_info = workflow_client.run()
            run_info['info'] = filter_error_info(run_info['info'])
            return json.dumps(run_info)
        except Exception as e:
            return json.dumps({"success": False, "info": "Running error: %s" %(filter_error_info(str(e)))})

