# -*- coding: utf-8 -*-
# __author__ = 'shenghe'

from biocluster.wpm.client import wait
import json


class TaskWait(object):
    def GET(self, task_id):
        print "WAIT TASK:", task_id
        try:
            result = wait(task_id, timeout=7200)
            if result:
                return json.dumps({"info": "等待结束", "success": True})
            else:
                return json.dumps({"info": "等待超时", "success": False})
        except Exception as e:
            return json.dumps({"info": str(e), "success": False})
