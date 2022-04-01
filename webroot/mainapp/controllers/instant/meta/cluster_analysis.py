# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import param_pack, group_detail_sort
from mainapp.libs.signature import check_sig
from bson import SON
from bson import ObjectId


class ClusterAnalysis(MetaController):
    def __init__(self):  # 20170106 2 lines
        super(ClusterAnalysis, self).__init__(instant=True)

    @check_sig
    def POST(self):
        # return_info = super(ClusterAnalysis, self).POST() # 20170106 3 lines
        # if return_info:
        #     return return_info
        data = web.input()
        postArgs = ["group_detail", "otu_id", "level_id", "task_type", "group_id"] # guanqing.zou 20180514
        #postArgs = ["group_detail", "otu_id", "level_id", "task_type", "group_id", "group_method", "combine_value"]  # modify by guanqign.zou 2018.04.04
        group_method = ""
        combine_value = "0.00"

        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '{}parameters missing!'.format(arg)}
                return json.dumps(info)

        ###guanqing.zou 20180514 
        if  hasattr(data,"group_method"):
            group_method = data.group_method
        if  hasattr(data,"combine_value"):
            combine_value = data.combine_value
                
        task_name = 'meta.report.cluster_analysis'
        otu_info = self.meta.get_otu_table_info(data.otu_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！", 'code':'C2200401'}
            return json.dumps(info)
        task_info = self.meta.get_task_info(otu_info['task_id'])
        params_json = {
            'otu_id': data.otu_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'submit_location': data.submit_location,
            'task_type': data.task_type,
            'group_method': group_method,   # by guanqing.zou 2018.5.14
            'combine_value': combine_value  #by guanqing.zou 2018.5.14
        }
        main_table_name = 'CommunityBarPie_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]  #modified by hongdongxuan 201703221 将ClusteringAnalysis_改为CommunityBarPie
        newick_id = None
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('from_id', data.otu_id),  # maybe ObjectId(data.otu_id)
            ('name', main_table_name),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('newick_id', newick_id),
            ('status', 'start'),
            ('desc', 'otu table after Cluster Analysis'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("show", 0),
            ("type", "otu_group_analyse")
        ]
        main_table_id = self.meta.insert_none_table('sg_otu')
        update_info = {str(main_table_id): 'sg_otu'}
        options = {
            "input_otu_id": data.otu_id,
            "in_otu_table": data.otu_id,
            "group_detail": data.group_detail,
            "level": str(data.level_id),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
            'main_table_data': SON(mongo_data),
            'method' : group_method,            #guanqing.zou 2018.5.14
            'combine_value' : combine_value  #guanqing.zou 2018.5.14
        }
        to_file = "meta.export_otu_table_by_level(in_otu_table)"
        print "set_sheet_data"
        self.set_sheet_data(name=task_name, options=options, main_table_name="CommunityAnalysis/" + main_table_name,
                            module_type='workflow', to_file=to_file)
        task_info = super(ClusterAnalysis, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
