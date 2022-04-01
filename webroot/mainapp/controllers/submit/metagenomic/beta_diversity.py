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
from .comm_creat import CommCreatAction


class BetaDiversityAction(MetagenomicController):
    def __init__(self):
        super(BetaDiversityAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        # default_argu = ['analysis_type', 'level_id', 'submit_location', 'group_id', 'group_detail','anno_id','geneset_id','second_level','anno_type']
        default_argu = ['analysis_type', 'submit_location', 'group_id', 'group_detail', 'geneset_id', 'anno_type',
                        'method', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'PARAMETERS MISSING: %s' % argu}
                return json.dumps(info)
        add_personal = CommCreatAction()
        info = add_personal.judge_database(data)
        if info != "":
            return json.dumps(info)
        if hasattr(data, 'scale'):
            if data.scale not in ["T","F"]:
                info = {'success': False, 'info': 'PARAMETERS ERROR: value of scale must in [T, F] !'}
                return json.dumps(info)
        task_name = 'metagenomic.report.beta_diversity'
        module_type = 'workflow'
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        group_detail = group_detail_sort(data.group_detail)
        params_json = {
            'analysis_type': data.analysis_type,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'group_id': data.group_id,
            'group_detail': group_detail,
            'geneset_id': data.geneset_id,
            'anno_type': data.anno_type,
            'method': data.method
        }
        env_id = None
        env_labs = ''
        distance_method = ''
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('table_type', data.analysis_type),
            ('status', 'start'),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        if data.anno_type != 'gene':
            params_json['level_id'] = int(data.level_id)
            params_json['anno_id'] = data.anno_id
            level = '_' + self.level_id(data.level_id)
            if hasattr(data, 'second_level'):
                params_json['second_level'] = data.second_level
        else:
            level = ''
        sample_len = sum([len(i) for i in group_detail.values()])
        if data.analysis_type == 'pca':
            if hasattr(data, 'scale'):
                params_json['scale'] = data.scale
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of env_id, can not transformed to ObjectId (%s)' % data.env_id}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列','code':'C2401001', 'variables':''}
                    return json.dumps(info)
        elif data.analysis_type == 'pcoa' or data.analysis_type == 'nmds':
            if sample_len < 3:
                info = {'success': False, 'info': '样本数量少于3，不可进行此分析！','code':'C2401002', 'variables':''}
                return json.dumps(info)
            if not hasattr(data, 'distance_method'):
                info = {'success': False, 'info': 'PARAMETERS MISSING: distance_method'}
                return json.dumps(info)
            params_json['distance_method'] = data.distance_method
            distance_method = data.distance_method
        elif data.analysis_type == 'dbrda':
            if not hasattr(data, 'distance_method'):
                info = {'success': False, 'info': 'PARAMETERS MISSING: distance_method'}
                return json.dumps(info)
            params_json['distance_method'] = data.distance_method
            distance_method = data.distance_method
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of env_id, can not transformed to ObjectId (%s)' % data.env_id}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列', 'code':'C2401001', 'variables':''}
                    return json.dumps(info)
            else:
                # info = {'success': False, 'info': 'dbrda分析缺少参数:env_id!'}
                info = {'success': False, 'info': 'PARAMETERS MISSING: env_id, when analysis_type==dbrda'}  # modified by hongdongxuan 20170310
                return json.dumps(info)
        elif data.analysis_type == 'rda_cca':
            if hasattr(data, 'env_id'):
                params_json['env_id'] = data.env_id
                env_id = self.check_objectid(data.env_id)
                if not env_id:
                    info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of env_id, can not transformed to ObjectId (%s)' % data.env_id}
                    return json.dumps(info)
                if hasattr(data, 'env_labs'):
                    params_json['env_labs'] = data.env_labs
                    env_labs = data.env_labs
                else:
                    info = {'success': False, 'info': '没有选择任何环境因子列','code':'C2401001', 'variables':''}
                    return json.dumps(info)
            else:
                # info = {'success': False, 'info': 'rda_cca分析缺少参数:env_id!'}
                info = {'success': False, 'info': 'PARAMETERS MISSING: env_id, when analysis_type==rda_cca'}  # modified by hongdongxuan 20170310
                return json.dumps(info)
        else:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of analysis_type (%s)' % data.analysis_type}
            return json.dumps(info)
        [geneset_table,gene_list] = metagenomic.export_geneset_table(data.geneset_id,data.method)
        options = {
            'analysis_type': data.analysis_type,
            'group_detail': data.group_detail,
            'distance_method': distance_method,
            'env_labs': env_labs,
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            'group_id': data.group_id,
            'geneset_id': data.geneset_id,
            'geneset_table': geneset_table,
            'method': data.method,
            'anno_type': data.anno_type
        }
        to_file = []
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        mongo_data.append(('env_id', env_id))
        if hasattr(data, 'scale'):
            options['scale'] = data.scale
        if env_id:
            to_file.append('metagenomic.export_env_table(env_file)')
            options['env_file'] = data.env_id
            options['env_id'] = data.env_id
        if data.anno_type != 'gene':
            options['level_id'] = self.level_id(data.level_id)
            anno_info = metagenomic.get_anno_info(data.anno_id,data.anno_type)
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
        main_table_name = BetaDiversityAction.get_main_table_name(data.anno_type, level, data.analysis_type)
        main_table_name = main_table_name.replace(" ","_")
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        # main_table_id = metagenomic.insert_main_table('sg_beta_multi_analysis', mongo_data)
        # main_table_id = metagenomic.insert_none_table('beta_diversity')
        main_table_id = self.metagenomic.insert_main_table('beta_diversity', mongo_data)
        update_info = {str(main_table_id): 'beta_diversity'}
        options['update_info'] = json.dumps(update_info)
        options['main_id'] = str(main_table_id)
        options['main_table_data'] = SON(mongo_data)
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name,
                            module_type=module_type, to_file=to_file, project_sn=task_info['project_sn'],
                            task_id=task_info['task_id'], params=params_json)
        task_info = super(BetaDiversityAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

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
    def get_main_table_name(anno_type, level, analysis_type):
        time_now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        database_list = ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy", "gene",'go','phi','qs','mvirdb','tcdb','pfam','cyps','probio']
        if anno_type in database_list:
            database = anno_type.upper()
        else:
            raise Exception('错误的注释类型')
        if analysis_type == 'pca':
            return 'PCA_%s%s_%s' % (database, level, time_now)
        elif analysis_type == 'pcoa':
            return 'PCoA_%s%s_%s' % (database, level, time_now)
        elif analysis_type == 'nmds':
            return 'NMDS_%s%s_%s' % (database, level, time_now)
        # elif analysis_type == 'plsda':
        #     return 'PLS-DA_%s%s_%s' % (database, level, time_now)
        elif analysis_type == 'dbrda':
            return 'db-RDA_%s%s_%s' % (database, level, time_now)
        elif analysis_type == 'rda_cca':
            return 'RDACCA_%s%s_%s' % (database, level, time_now)
        else:
            raise Exception('PARAMETERS ERROR: wrong value of analysis_type')
