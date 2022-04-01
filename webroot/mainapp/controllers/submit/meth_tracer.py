# -*- coding: utf-8 -*-
# __author__ = 'luotong'

from mainapp.controllers.core.basic import Basic
from mainapp.libs.signature import check_sig
from biocluster.core.function import filter_error_info
from mainapp.models.workflow import Workflow
from biocluster.config import Config
from biocluster.wpm.client import *
import web
import json
import datetime
import random
from pymongo import MongoClient
from bson.objectid import ObjectId


PACKAGE_URL = "meth_tracer"
class MethTracerAction(object):
    """
    MethTracer设置
    """
    def __init__(self):
        super(MethTracerAction, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()

        data['task_id'] = "meth_tracer" + str(random.randint(1000, 9999))
        requires = ['path', 'model', 'email']
        for i in requires:
            if not (hasattr(data, i)):
                return json.dumps({"success": False, "info": "缺少%s参数!" % i})

        workflow_id = "MethTracer_" + "{}_{}".format(data['task_id'], datetime.datetime.now().strftime("%H%M%S%f")[:-3])
        data['path'] = data['path'].replace("/web/html/webroot/static/upload/data/", "/mnt/ilustre/institute/ctmeth/")
        _client = MongoClient('10.100.200.131',27017)
        #_client = MongoClient('192.168.10.186',27017)

        _db = _client.cmeth
        _job = _db.job
        _posts = {
                  'step':0,
                  'status':1,
                  'note': '',
                  'default': 0,
                  'model': int(data.model),
                  'params':{
                        #'input':'/mnt/ilustre/users/sanger-dev/sg-users/luotong/ctMethTracer/10.test/knn/test.bed',
                        'input':data['path'],
                        'model':int(data['model']),
                  },
                  'email':data['email'],
                  'has_sent_email': 1,
                  'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
              }

        print _posts['params']
        _jobID = _job.insert_one(_posts).inserted_id
        #_job.update({'_id':_jobID},{"$set":{'jobID':str(_jobID)}})

        json_data = {
          'id': workflow_id,
          'stat_id': 0,
          # 'type': 'workflow',
          'type': "tool",
          # 'name': "copy_demo.demo_init",  # 需要配置
          'name': "meth_tracer",
          'client': data.client,
          "IMPORT_REPORT_DATA": False,
          "IMPORT_REPORT_AFTER_END": False,
          'options': {
              'path': data['path'],
              'model': data['model'],
              'jobID': str(_jobID),
              #'email': data['email'],
          },
        }

        try:
            workflow_client = Basic(data=json_data, instant=False)
            info = workflow_client.run()
            info['info'] = filter_error_info(info['info'])

            if "success" in info.keys() and info["success"]:
                res = {"success": True, "info": "MethTracer success!", "jobid": str(_jobID)}
                return json.dumps(res)
            else:
                print 'error'
                return json.dumps({"success": False, "info": "MethTracer fail."})
        except:
            print 'error2'
            return json.dumps({"success": False, "info": "MethTracer fail."})
