# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import web
import json
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.libs.param_pack import group_detail_sort
import types
from mainapp.models.mongo.metagenomic import Metagenomic
import datetime
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from bson import SON
from mainapp.libs.signature import check_sig
import time
from .comm_creat import CommCreatAction

class HclusterAction(MetagenomicController):
    def __init__(self):
        super(HclusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # default_argu = ['analysis_type', 'level_id', 'submit_location', 'group_id', 'group_detail','anno_id','geneset_id','second_level','anno_type']
        default_argu = ['submit_location', 'group_id', 'group_detail', 'geneset_id', 'anno_type',
                        'hcluster_method', 'distance_method', 'method', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'PARAMETERS MISSING: %s' % argu}
                return json.dumps(info)
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        task_name = 'metagenomic.report.hcluster'
        module_type = 'workflow'
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        # return data.group_detail
        group_detail = group_detail_sort(data.group_detail)
        params_json = {
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_id': data.group_id,
            'group_detail': group_detail,
            'geneset_id': data.geneset_id,
            'anno_type': data.anno_type,
            'method': data.method,
            'hcluster_method': data.hcluster_method,
            'distance_method': data.distance_method
        }
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        if data.anno_type != 'gene':
            params_json['level_id'] = int(data.level_id)
            params_json['anno_id'] = data.anno_id
            # params_abu_json['level_id'] = data.level_id
            # params_abu_json['anno_id'] = data.anno_id
            level = '_' + self.level_id(data.level_id)
            if hasattr(data, 'second_level'):
                params_json['second_level'] = data.second_level
                # params_abu_json['second_level'] = data.second_level
        else:
            level = ''
        # params_abu = json.dumps(params_abu_json, sort_keys=True, separators=(',', ':'))
        # abu_mongo.append(('params', params_abu))
        [geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        to_file = []
        options = {
            # 'second_level': data.second_level,
            'distance_method': data.distance_method,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            'geneset_table': geneset_table,
            'hcluster_method': data.hcluster_method,
            'group_detail': data.group_detail,
            'task_id': task_info['task_id'],
            'method': data.method,
            'anno_type': data.anno_type
        }
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        if data.anno_type != 'gene':
            options['level_id'] = self.level_id(data.level_id)
            anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            options['anno_id'] = data.anno_id
            options['anno_table'] = self.use_s3(anno_info['anno_file'])
            if data.anno_type != 'nr':
                options['lowest_level'] = self.use_s3(anno_info['lowest_level'])
            if hasattr(data, 'second_level'):
                if data.second_level != "":
                    options['second_level'] = self.level_convert(data.second_level, data.level_id)
        else:
            options['gene_list'] = gene_list
        options['group_table'] = data.geneset_id
        to_file.append('metagenomic.export_group_table_by_detail(group_table)')
        main_table_name = HclusterAction.get_main_table_name(data.anno_type, level).replace(" ","_")
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        # main_table_id = metagenomic.insert_none_table('hcluster_tree')
        main_table_id = self.metagenomic.insert_main_table('hcluster_tree', mongo_data)
        update_info = {str(main_table_id): 'hcluster_tree'}
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=module_type, project_sn=task_info['project_sn'], to_file=to_file,
                            params=params_json, task_id=task_info['task_id'])
        task_info = super(HclusterAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def _update_status_api(self):
        """
        根据client决定接口api为metagenomic.update_status/metagenomic.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'metagenomic.update_status'
        else:
            return 'metagenomic.tupdate_status'
    
    def check_objectid(self, in_id):
        """检查一个id是否可以被ObjectId"""
        if isinstance(in_id, types.StringTypes):
            try:
                in_id = ObjectId(in_id)
            except InvalidId:
                return False
        elif isinstance(in_id, ObjectId):
            pass
        else:
            return False
        return in_id

    @staticmethod
    def get_main_table_name(anno_type, level):
        time_now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        database_list = ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy", "gene",'go','phi','qs','mvirdb','tcdb','pfam','cyps','probio']
        if anno_type in database_list:
            database = anno_type.upper()
        else:
            raise Exception('PARAMETERS ERROR: wrong value of anno_type')
        return 'Hcluster_%s%s_%s' % (database, level, time_now)
