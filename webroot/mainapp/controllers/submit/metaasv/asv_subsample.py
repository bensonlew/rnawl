# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from bson.objectid import ObjectId
from collections import OrderedDict
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import param_pack, group_detail_sort, filter_json_sort
from mainapp.libs.signature import check_sig


class AsvSubsampleAction(MetaasvController):
    """
    ASV 分类学分析
    """
    def __init__(self):
        super(AsvSubsampleAction, self).__init__(instant=False)
    
    @check_sig
    def POST(self):
        data = web.input()
        print("+++++++++++++++++++++++%s++++++++++++++++++++++"%data)
        postArgs = ['size', 'submit_location', "asv_id", "group_detail", "group_id"]
        for arg in postArgs:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'parameters missing:%s' % arg}
                return json.dumps(info)

        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        # task_info = self.metaasv.get_task_info(otu_info['task_id'])
        project_sn = otu_info['project_sn']
        task_id = otu_info['task_id']
        task_type = 'workflow'
        main_table_name = 'ASV_' + \
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        my_param = dict()
        my_param["group_id"] = data.group_id
        my_param['asv_id'] = data.asv_id
        my_param["submit_location"] = data.submit_location
        my_param["size"] = data.size
        # my_param["filter_json"] = filter_json_sort(data.filter_json)

        if hasattr(data, "filter_json"):
            if data.filter_json == "":
                filter_json = data.filter_json
                set_id = ""
            else:
                filter_json = json.loads(data.filter_json)
                filter_json = OrderedDict(sorted(filter_json.items()))
                if "asv_filter" in filter_json: ##20200709 增加此功能用于前端判断asv集参与哪些分析
                    set_list = []
                    new_asv_filter = []
                    asv_filter = filter_json["asv_filter"]
                    for asv_dict in asv_filter:
                        new_asv_dict = OrderedDict(sorted(asv_dict.items()))
                        if "set_id" in asv_dict:
                            new_set = asv_dict["set_id"]
                            if new_set not in set_list:
                                set_list.append(new_set)
                        if new_asv_dict not in new_asv_filter:
                            new_asv_filter.append(new_asv_dict)
                    new_asv_filter.sort()
                    filter_json["asv_filter"] = new_asv_filter
                    set_id = ",".join(set_list)
                else:
                    set_id = ""
            # my_param["filter_json"] = json.dumps(filter_json, sort_keys=True, separators=(',', ':'))
            my_param["filter_json"] = filter_json
        else:
            set_id = ""
        my_param["group_detail"] = group_detail_sort(data.group_detail)
        my_param["task_type"] = str(data.task_type)
        params = param_pack(my_param)
        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('asv_id', ObjectId(data.asv_id)),
            ('name', main_table_name),
            ("params", params),
            ('status', 'start'),
            ("level_id", int(9)),
            ('desc', '正在计算'),
            ('set_id', set_id),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]

        main_table_id = self.metaasv.insert_main_table('asv', mongo_data)
        update_info = {str(main_table_id): 'asv'}

        # convert_filter_json = self.metaasv.convert_json(filter_json)
        options = {
            "in_otu_table": data.asv_id,
            "input_otu_id": data.asv_id,
            "group_detail": data.group_detail,
            "group_table": data.group_id,
            "level": "9",
            "size": data.size,
            'update_info': json.dumps(update_info),
            # "params": params,
            'main_id': str(main_table_id),
        }
        if hasattr(data, "filter_json"):
            options["filter_json"] = data.filter_json
        to_file = ["metaasv.export_otu_table_by_level_remove_blank(in_otu_table)", "metaasv.export_group_table_by_detail(group_table)"]
        task_name = 'metaasv.report.asv_subsample'
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name="ASV/" + main_table_name, module_type=task_type,
                            to_file=to_file)
        task_info = super(AsvSubsampleAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
