# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

import web
import json
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from biocluster.config import Config
from mainapp.models.mongo.submit.denovo_rna.denovo_cluster import DenovoExpress
import types
from mainapp.models.mongo.meta import Meta
from mainapp.models.workflow import Workflow
from mainapp.controllers.project.denovo_controller import DenovoController
import datetime


class Cluster(DenovoController):
    def __init__(self):
        super(Cluster, self).__init__()

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        #print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        my_param = dict()
        my_param["submit_location"] = data.submit_location
        my_param["express_id"] = data.express_id
        my_param["distance_method"] = data.distance_method
        my_param["cluster_method"] = data.cluster_method
        my_param["log"] = data.log
        my_param["sub_num"] = data.sub_num
        my_param['gene_list'] = data.gene_list
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        express_info = Meta(db=self.mongodb).get_main_info(data.express_id, 'sg_denovo_express')
        if express_info:
            main_table_name = "Cluster_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            task_id = express_info["task_id"]
            project_sn = express_info["project_sn"]
            cluster_id = DenovoExpress().add_cluster(params=params, project_sn=project_sn, task_id=task_id)
            update_info = {str(cluster_id): 'sg_denovo_cluster'}
            update_info = json.dumps(update_info)
            options = {
                "express_file": data.express_id,
                "update_info": update_info,
                "cluster_id": str(cluster_id),
                "distance_method": data.distance_method,
                "sub_num": data.sub_num,
                "cluster_method": data.cluster_method,
                "log": data.log,
                "gene_list": data.gene_list
            }
            to_file = "denovo.export_express_matrix(express_file)"
            self.set_sheet_data(name='denovo_rna.report.cluster', options=options, main_table_name=main_table_name, module_type='workflow', to_file=to_file, main_id=cluster_id, collection_name="sg_denovo_cluster")
            task_info = super(Cluster, self).POST()
            task_info['content'] = {'ids': {'id': str(cluster_id), 'name': main_table_name}}
            #print task_info
            return json.dumps(task_info)
        else:
            info = {"success": False, "info": "express_id不存在，请检查参数是否正确！"}
            return json.dumps(info)

    def check_options(self, data):
        """
        检查网页端传来的参数是否正确
        """
        params_name = ["submit_location", "distance_method", "cluster_method", "log", "sub_num", "express_id"]
        success = []
        for names in params_name:
            if not hasattr(data, names):
                success.append("缺少参数！")
        express_id = str(data.express_id)
        if not isinstance(express_id, ObjectId) and not isinstance(express_id, types.StringType):
            success.append("传入的express_id {}不是一个ObjectId对象或字符串类型".format(express_id))
        return success
