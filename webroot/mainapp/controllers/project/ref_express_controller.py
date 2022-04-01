# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
import web
import random
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from meta_controller import MetaController
from biocluster.config import Config
from bson.objectid import ObjectId
import json

class RefExpressController(MetaController):
    def __init__(self, instant=False):
        super(RefExpressController, self).__init__(instant)
        #self.mongodb = Config().MONGODB + '_ref_rna'
        self.mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]

    def _update_status_api(self):
        """
        根据client决定接口api为ref.update_status/ref.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'ref.ref_update_status'
        else:
            return 'ref.ref_status'

    def set_sheet_data(self, *arg, **kwarg):
        print arg, kwarg
        super(RefExpressController, self).set_sheet_data(*arg, **kwarg)

    def get_main_info(self, express_id=None, task_id=None):  # add function by khl 20170414
        collection=self.mongodb["sg_express"]
        if express_id:
             if isinstance(express_id, types.StringTypes):
                 express_id = ObjectId(express_id)
             elif isinstance(express_id, ObjectId):
                 express_id = express_id
             else:
                 raise Exception("输入express_id参数必须为字符串或者ObjectId类型!")
             result = collection.find_one({"_id":express_id})
        elif task_id:
            result = collection.find_one({"task_id":task_id})
        #result = collection.find_one({"_id":express_id})
        else:
            raise Exception("必须输入{}或{}参数！".format("express_id", "task_id"))
        return result
    
    @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return run_info
        except Exception, e:
            return {"success": False, "info": "运行出错: %s" % e }

    def get_express_id(self, task_id, _type, express_method): #add by khl 20170426
        """
        暂时实现的功能是根据表达量的软件如RSEM和表达量水平FPKM 获取表达量的主表(tsanger_ref_rna["sg_express"]) id
        :params _type: 表达量水平fpkm/tpm
        :params query_type: gene or transcript
        :params express_method: 表达量方法 featurecounts/rsem
        """
        collection=self.mongodb["sg_express"]
        db=collection.find({"task_id":task_id})
        try:
            for d in db:
                _id = d["_id"]
                params=json.loads(d['params'])
                print params
                print params['type']
                print params['express_method']
                if _type == params['type'] and express_method == params['express_method'].lower():
                      return _id
                else:
                      pass
        except Exception:
                print "没有找到{}和{}对应的express_id".format(_type,express_method)

