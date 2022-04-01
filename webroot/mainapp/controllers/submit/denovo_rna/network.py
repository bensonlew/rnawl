# -*- coding: utf-8 -*-
# __author__ = 'konghualei'
import web
import json
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from mainapp.libs.param_pack import *
from biocluster.config import Config
from mbio.api.database.denovo_network import *
from mainapp.models.mongo.submit.denovo_rna.denovo_network import DenovoNetwork
import types
#from mainapp.models.mongo.denovo import Denovo
from mainapp.models.mongo.meta import Meta
from mainapp.controllers.project.denovo_controller import DenovoController
from mainapp.models.workflow import Workflow
import random
import datetime

class Network(DenovoController):
    def __init__(self):
        super(Network, self).__init__(instant=False)

    @check_sig
    def POST(self):
        """web用法"""
        data = web.input()
        client = data.client if hasattr(data, 'client') else web.ctx.env.get('HTTP_CLIENT')
        #print data
        return_result = self.check_options(data)
        if return_result:
            info = {'success': False, 'info': "+".join(return_result)}
            return json.dumps(info)
        my_param = dict()
        my_param['express_id'] = data.express_id
        my_param['softpower'] = data.softpower
        my_param['module'] = data.module
        #my_param['network'] = data.network
        my_param['submit_location'] = data.submit_location
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        express_info = Meta(db=self.mongodb).get_main_info(data.express_id, 'sg_denovo_express')
        if express_info:
            main_table_name = "Network_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            task_id = express_info['task_id']
            project_sn = express_info['project_sn']
            network_express_id = DenovoNetwork().add_network(params = params, softpower_png = None, module_png = None, module = None, name = None, project_sn = project_sn, task_id = task_id)
            update_info ={str(network_express_id):"sg_denovo_network"}
            update_info = json.dumps(update_info)
            options = {
                    'express_file': data.express_id,
                    'update_info': update_info,
                    'softpower': data.softpower,
                    'module': data.module,
                    'network_express_id': str(network_express_id)
            }
            to_file = "denovo.export_express_matrix(express_file)",
            self.set_sheet_data(name='denovo_rna.report.network', options=options,main_table_name=main_table_name,module_type='workflow',to_file=to_file, main_id =network_express_id, collection_name = 'sg_denovo_network')
            #self.set_sheet_data(name='denovo_rna.report.network', options=options, main_table_name=main_table_name, module_type='workflow', to_file=to_file, params=my_param)
            task_info = super(Network, self).POST()
            task_info['content'] = {'ids': {'id': str(network_express_id), 'name':main_table_name}}
            #print task_info
            return json.dumps(task_info)
        else:
            info = {'success': False, "info": "express_id不存在，请确认参数是否正确！"}
            return json.dumps(info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ['express_id', 'softpower', 'module', 'submit_location']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数！")
        if float(data.softpower) > 20 or float(data.softpower) < 1:
            success.append("softpower不在范围内")
        if float(data.module) > 1 or float(data.module) < 0:
            success.append("module权重值不在范围内")
        ids=str(data.express_id)
        #print type(ids)
        if not isinstance(ids, ObjectId) and not isinstance(ids, types.StringTypes):
            success.append("传入的id: {}不是一个ObjectId对象或字符串类型".format(ids))
        return success
