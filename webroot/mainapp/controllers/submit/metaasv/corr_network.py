# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from bson import SON
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson import ObjectId


class CorrNetworkAction(MetaasvController):
    """
    metaasv 单因素相关性网络
    """
    def __init__(self):
        super(CorrNetworkAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['asv_id', 'level_id', 'submit_location', 'group_detail', 'group_id', 'ratio_method']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' , "variables":[ argu]}
                return json.dumps(info)
        if hasattr(data, "significance"):
            significance = float(data.significance)
        if int(data.level_id) not in range(1, 10):
            variables = []
            variables.append(data.level_id)
            info = {'success': False, 'info': 'level{}不在规定范围内{}'.format(data.level_id), 'variables':variables}
            return json.dumps(info)
        if int(data.level_id) in [1, 2]:
            info = {'success': False, 'info': 'OTU表中的物种数目太少，不能进行该分析，请选择Phylum以下的分类水平！'}
            return json.dumps(info)
        group_detail = json.loads(data.group_detail)
        if not isinstance(group_detail, dict):
            info = {'success': False, 'info': '传入的group_detail不是一个字典！'}
            return json.dumps(info)
        sample_num = 0
        for key in group_detail:
            sample_num += len(group_detail[key])
        if sample_num < 3:
            info = {'success': False, 'info': '样本的总个数少于3个，不能进行物种相关性网络分析！'}
            return json.dumps(info)
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "OTU不存在，请确认参数是否正确！"}
            return json.dumps(info)
        task_name = 'metaasv.report.corr_network'
        module_type = 'workflow'
        task_info = self.metaasv.get_task_info(otu_info['task_id'])
        level_name = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"]
        main_table_name = 'CorrNetwork' + data.ratio_method.capitalize() + level_name[int(data.level_id) - 1] + '_' + \
                          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        params_json = {
            'asv_id': data.asv_id,
            'level_id': int(data.level_id),
            'group_id': data.group_id,
            'group_detail': group_detail_sort(data.group_detail),
            'ratio_method': data.ratio_method,
            'submit_location': data.submit_location,
            'task_type': str('2'),
        }
        if hasattr(data, 'abundance'):
            params_json['abundance'] = str(data.abundance)
        if hasattr(data, 'lable'):
            params_json['lable'] = data.lable
        if hasattr(data, 'color_level'):
            params_json['color_level'] = int(data.color_level)
        if hasattr(data, "coefficient"):
            params_json["coefficient"] = str(data.coefficient)
        if hasattr(data, "significance"):
            params_json["significance"] = str(data.significance)
        if data.group_id == 'all':
            group__id = data.group_id
        else:
            group__id = ObjectId(data.group_id)

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('asv_id', ObjectId(data.asv_id)),
            ('group_id', group__id),
            ('status', 'start'),
            ('desc', 'corr_network分析'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('level_id', int(data.level_id)),
            ('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('name', main_table_name),
        ]

        main_table_id = self.metaasv.insert_none_table('corr_network')
        update_info = {str(main_table_id): 'corr_network'}
        options = {
            'otutable': data.asv_id,
            'grouptable': data.group_id,
            'group_detail': data.group_detail,
            'method': data.ratio_method,
            'level': int(data.level_id),
            'main_table_data': SON(mongo_data),
            'update_info': json.dumps(update_info),
            'main_id': str(main_table_id),
        }
        if hasattr(data,'lable'):
            options['lable'] = data.lable
        if hasattr(data,'color_level'):
            options['color_level'] = data.color_level
        if hasattr(data, "coefficient"):
            options["coefficient"] = float(data.coefficient)
        if hasattr(data, "significance"):
            options["significance"] = float(data.significance)
        if hasattr(data, "abundance"):
            options["abundance"] = int(data.abundance)

        to_file = ["metaasv.export_otu_table_by_detail(otutable)", "metaasv.export_group_table_by_detail(grouptable)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="CorrNetwork/" + main_table_name,
                            module_type=module_type, to_file=to_file)
        task_info = super(CorrNetworkAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        print "task_info", task_info
        return json.dumps(task_info)
