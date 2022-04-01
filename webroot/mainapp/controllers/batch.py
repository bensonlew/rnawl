# -*- coding: utf-8 -*-
# __author__ = 'quan.guo'
import web
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
import json


class Batch(object):

    @check_sig
    def GET(self):
        data = web.input()
        if hasattr(data, "batch_id"):
            wid = data.batch_id
        else:
            raise web.badrequest
        model = Workflow()
        data_list = model.batch_list(wid)
        output = []
        for d in data_list:
            output.append({
                "workflow_id": d["workflow_id"],
                "is_end": d["is_end"],
                "server": d["server"],
                "is_error": d["is_error"],
                "has_run": d["has_run"],
                "pid": d["pid"],
                "cpu_used": d["cpu_used"],
                "memory_used": d["memory_used"],
                "work_dir": d["work_dir"],
                "error": d["error"],
                "error_data": d["error_data"]
                })
        return json.dumps(output)
