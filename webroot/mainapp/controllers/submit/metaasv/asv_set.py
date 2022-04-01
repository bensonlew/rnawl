# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.signature import check_sig
from bson import SON
import datetime


class AsvSetAction(MetaasvController):
    """
    Metaasv 创建基因集
    """
    def __init__(self):
        super(AsvSetAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['table_name', 'table_id', 'submit_location', 'asv_name', 'desc', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing: %s' % argu}
                return json.dumps(info)

        task_name = 'metaasv.report.asv_set'
        module_type = 'workflow'
        params_json = {
            'table_name': data.table_name,
            'table_id': data.table_id,
            'submit_location': data.submit_location,
            'asv_name': data.asv_name,
            'desc': data.desc,
            'task_type': data.task_type
        }
        if hasattr(data, "pvalue"):
            params_json['pvalue']=  float(data.pvalue)
        if hasattr(data, "qvalue"):
            params_json['qvalue']=  float(data.qvalue)
        if hasattr(data, "lda"):
            params_json['lda']=  float(data.lda)
        if hasattr(data, "top"):
            params_json['top']=  int(data.top)
        if hasattr(data, "species_name"):
            params_json['species_name']= str(data.species_name)
        if hasattr(data, "label"):
            params_json['label']= str(data.label)
        main_table_name = 'ASV_set' + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        asv_id = self.metaasv.get_asv_id_from_table(data.table_name, data.table_id)
        print(asv_id)
        otu_info = self.metaasv.get_otu_table_info(asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        asv_name_list = self.metaasv.find_asv_name(otu_info['task_id'])
        if data.asv_name in asv_name_list:
            info = {"success": False, "info": "ASV集名称已经存在，请重新命名！"}
            return json.dumps(info)
        mongo_data = [
            ('status', 'start'),
            ('desc', 'computing'),
            ('name', data.asv_name),
            ('project_sn', otu_info['project_sn']),
            ('task_id', otu_info['task_id']),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.metaasv.insert_main_table('asv_set', mongo_data)
        update_info = {str(main_table_id): 'asv_set'}
        options = {
            'table_name': data.table_name,
            'table_id': data.table_id,
            'main_id': str(main_table_id),
            'update_info': json.dumps(update_info),
            'main_table_data': SON(mongo_data)
        }
        if hasattr(data, "pvalue"):
            options['pvalue']=  float(data.pvalue)
        if hasattr(data, "qvalue"):
            options['qvalue']=  float(data.qvalue)
        if hasattr(data, "lda"):
            options['lda']=  float(data.lda)
        if hasattr(data, "top"):
            options['top']=  int(data.top)
        if hasattr(data, "species_name"):
            options['species_name']= str(data.species_name)
        if hasattr(data, "label"):
            options['label']= str(data.label)
        to_file = []
        self.set_sheet_data(name=task_name, options=options, main_table_name="ASV_set/" + main_table_name,
                            module_type=module_type, to_file=to_file, main_id=str(main_table_id),collection_name="asv_set")
        task_info = super(AsvSetAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)
