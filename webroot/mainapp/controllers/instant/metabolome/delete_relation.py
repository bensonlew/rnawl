# -*- coding: utf-8 -*-
import web
import json
import os
import datetime
from mainapp.libs.signature import check_sig
from biocluster.core.function import filter_error_info
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome


class DeleteRelationAction(MetabolomeController):
    def __init__(self):
        """
        1.用于页面删除同时删除mongo数据；
        2.用于页面删除task或者批量删除，同时删除mongo数据；
        :return:
        """
        super(DeleteRelationAction, self).__init__(instant=True)

    def POST(self):
        data = web.input()
        requires = ['task_id']#分别是指项目类型，任务id
        for i in requires:
            if not hasattr(data, i):
                return json.dumps({"success":False, "info": "Lack params: %s" %(i)})
        name = "metabolome.report.delete_relation"
        table_name = "DeleteRelation_" + "{}_{}".format(data.task_id, datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        metabolome = Metabolome()
        r = metabolome.conmon_find_one('sg_task',{'task_id': data.task_id})
        project_sn = r['project_sn']
        project_type = r['type']

        self.set_sheet_data(name=name, options={"task_id": data.task_id},
                            main_table_name=table_name,
                            task_id=data.task_id, project_sn=project_sn,
                            module_type="workflow")
        post_info = super(DeleteRelationAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'name': table_name}}
        return json.dumps(post_info)
