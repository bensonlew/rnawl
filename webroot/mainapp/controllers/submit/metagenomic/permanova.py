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

class PermanovaAction(MetagenomicController):
    def __init__(self):
        super(PermanovaAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        # default_argu = ['analysis_type', 'level_id', 'submit_location', 'group_id', 'group_detail','anno_id','geneset_id','second_level','anno_type']
        default_argu = ['submit_location', 'group_detail', 'geneset_id', 'anno_type',
                        'permutation', 'distance_method', 'method', 'task_type', 'group_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'PARAMETERS MISSING: %s' % argu}
                return json.dumps(info)
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        if (not hasattr(data, 'env_id')) and (not hasattr(data, 'group_detail_list')):
            info = {'success': False, 'info': '参数缺失：必须提供分组或环境因子中的至少一个', 'code':'C2402201', 'variables':''}
            return json.dumps(info)
        distance_name = ['euclidean', 'binary_euclidean', 'manhattan', 'binary_manhattan', 'gowerM', 'binary_gowerM',
                              'altGower', 'binary_altGower', 'canberraNZ', 'binary_canberraNZ', 'bray_curtis', 'binary_bray_curtis',
                              'kulczynski', 'binary_kulczynski', 'morisita_horn', 'binary_morisita_horn', 'morisita',
                              'binomial', 'binary_binomial', 'cao', 'binary_cao', 'chao', 'jaccard', 'binary_jaccard',
                              'raup_crick', 'mountford', 'mahalanobis']
        if data.distance_method not in distance_name:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of distance_method (%s)' % data.distance_method}
            return json.dumps(info)
        task_name = 'metagenomic.report.permanova'
        module_type = 'workflow'
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        # return data.group_detail
        group_detail = group_detail_sort(data.group_detail)
        params_json = {
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_detail': group_detail,
            'geneset_id': data.geneset_id,
            'anno_type': data.anno_type,
            'method': data.method,
            'group_id': data.group_id,
            'permutation': int(data.permutation),
            'distance_method': data.distance_method
        }
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
        if hasattr(data, 'group_detail_list'):
            params_json['group_detail_list'] = eval(data.group_detail_list)
        if hasattr(data, 'env_id'):
            params_json['env_id'] = data.env_id
            env_id = self.check_objectid(data.env_id)
            if not env_id:
                info = {'success': False, 'info': 'PARAMETERS ERROR: wrong type of env_id (%s), can not transformed to ObjectId' % data.env_id}
                return json.dumps(info)
            if hasattr(data, 'env_labs'):
                params_json['env_labs'] = data.env_labs
            else:
                info = {'success': False, 'info': '没有选择任何环境因子列', 'code':'C2402202', 'variables':''}
                return json.dumps(info)
        [geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        to_file = []
        options = {
            # 'second_level': data.second_level,
            'distance_method': data.distance_method,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'geneset_table': geneset_table,
            'permutation': int(data.permutation),
            'group_detail': data.group_detail,
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
        if hasattr(data, 'group_detail_list'):
            # options['group_list'] = data.group_list
            options['group_detail_list'] = data.group_detail_list
        if hasattr(data, 'env_id'):
            options['env_id'] = data.env_id
            if hasattr(data, 'env_labs'):
                options['env_labs'] = data.env_labs
        options['group_table'] = data.geneset_id
        options['envtable'] = data.geneset_id
        to_file.append('metagenomic.export_group_table_by_detail(group_table)')
        to_file.append('metagenomic.export_env_group_table(envtable)')
        main_table_name = PermanovaAction.get_main_table_name(data.anno_type, level)
        main_table_name = main_table_name.replace(" ","_")
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        # main_table_id = metagenomic.insert_none_table('hcluster_tree')
        main_table_id = self.metagenomic.insert_main_table('permanova', mongo_data)
        update_info = {str(main_table_id): 'permanova'}
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=module_type, project_sn=task_info['project_sn'], to_file=to_file,
                            params=params_json, task_id=task_info['task_id'])
        task_info = super(PermanovaAction, self).POST()
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
        return 'PERMANOVA_%s%s_%s' % (database, level, time_now)
