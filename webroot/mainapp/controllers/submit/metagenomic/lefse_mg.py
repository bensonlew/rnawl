# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import web
import json
import datetime
from bson.objectid import ObjectId
from bson.errors import InvalidId
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.libs.param_pack import group_detail_sort
from bson import SON
from mainapp.libs.signature import check_sig
from .comm_creat import CommCreatAction

class LefseMgAction(MetagenomicController):
    def __init__(self):
        super(LefseMgAction, self).__init__(instant = False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu=['submit_location', 'group_id', 'group_detail','second_group_detail', 'geneset_id', 'anno_type','method','strict']
        for argu in default_argu:
            if not hasattr(data,argu):
                info = {"success": False, "info": 'PARAMETERS MISSING: %s' % argu}
                return json.dumps(info)
        database_type = ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy", "gene", "annopersonal"]
        if not data.anno_type in database_type:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of anno_type (%s)' % data.anno_type}
            return json.dumps(info)
        group_detail = json.loads(data.group_detail)
        if not isinstance(group_detail, dict):
            info = {"success": False, "info": 'PARAMETERS ERROR: wrong type of group_detail, dict expected!'}
            return json.dumps(info)
        if data.group_id == 'all':
            info = {"success": False, "info": '分组方案至少选择两个分组！', 'code':'C2401701', 'variables':''}
            return json.dumps(info)
        elif len(group_detail) < 2:
            info = {"success": False, "info": '请选择至少两组以上的分组方案!', 'code':'C2401702', 'variables':''}
            return json.dumps(info)
        key1 = list(group_detail.values())
        for i in range(len(key1)):
            if (len(key1[i]) < 3):
                info = {"success": False, "info": '组内样本不能少于3个，请检查！', 'code':'C2401703', 'variables':''}
                return json.dumps(info)
        if data.second_group_detail != '"null"':
            second_group_detail = json.loads(data.second_group_detail)
            first = 0
            second = 0
            for i in group_detail.values():
                first += len(i)
            for n in second_group_detail.values():
                second += len(n)
            if not isinstance(second_group_detail, dict):
                info = {"success": False, "info": 'PARAMETERS ERROR: wrong type of second_group_detail, dict expected!'}
                return json.dumps(info)
            if first != second:
                info = {"success": False, "info": '二级分组与一级分组的样本数不相同，请检查！', 'code':'C2401704', 'variables':''}
                return json.dumps(info)
            else:
                fist_list = []
                second_list = []
                key1 = list(group_detail.values())
                key2 = list(second_group_detail.values())
                for i in range(len(key1)):
                    for j in range(len(key1[i])):
                        fist_list.append(key1[i][j])
                fist_list.sort()
                for i in range(len(key2)):
                    if (len(key2[i]) < 3):
                        info = {"success": False, "info": '二级组内样本不能少于3个，请检查！', 'code':'C2401705', 'variables':''}
                        return json.dumps(info)
                    else:
                        for j in range(len(key2[i])):
                            second_list.append(key2[i][j])
                second_list.sort()
                for a, b in zip(fist_list, second_list):
                    if a == b:
                        continue
                    else:
                        info = {"success": False, "info": '二级分组与一级分组的样本名称不一致，请检查！', 'code':'C2401706', 'variables':''}
                        return json.dumps(info)
        if int(data.strict) not in [1, 0]:
            info = {"success": False, "info": "PARAMETERS ERROR: wrong value of strict (%s), 0 or 1 expected!" % data.strict}
            return json.dumps(info)
        if data.anno_type != 'gene' and data.anno_type != 'nr':
            if not hasattr(data,'level_id'):
                info = {"success": False, "info": 'PARAMETERS MISSING: level_id'}
                return json.dumps(info)
        if data.anno_type == 'nr':
            if not hasattr(data,'start_level'):
                info = {"success": False, "info": 'PARAMETERS MISSING: start_level'}
                return json.dumps(info)
            else:
                if int(data.start_level) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
                    info = {"success": False, "info": 'PARAMETERS ERROR: wrong value of start_level'}
                    return json.dumps(info)
            if not hasattr(data,'end_level'):
                info = {"success": False, "info": 'PARAMETERS MISSING: start_level'}
                return json.dumps(info)
            else:
                if int(data.end_level) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
                    info = {"success": False, "info": 'PARAMETERS ERROR: wrong value of end_level'}
                    return json.dumps(info)
            if not data.domain_type in ["Viruses","Viroids","Eukaryota","Archaea","Bacteria"]:
                info = {"success": False, "info": 'PARAMETERS ERROR: wrong value of domain_type (%s)' % data.domain_type}
                return json.dumps(info)
        task_name ='metagenomic.report.lefse'
        task_type= 2
        module_type = 'workflow'
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_info = metagenomic.get_task_info(geneset_info['task_id'])
        group_detail = group_detail_sort(data.group_detail)
        params_json = {
            'submit_location': data.submit_location,
            'task_type': task_type,
            'group_id': data.group_id,
            'group_detail': group_detail,
            'geneset_id': data.geneset_id,
            'anno_type': data.anno_type,
            'method': data.method,
            'strict':int(data.strict)
        }
        group_id = data.group_id if data.group_id in ['all', 'All', 'ALL'] else ObjectId(data.group_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('geneset_id', ObjectId(data.geneset_id)),
            ('status', 'start'),
            ('group_id', group_id),
            ('desc', 'processing'),
            ('anno_type', data.anno_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_name = ''
        params_json['second_group_id'] = data.second_group_id
        if data.second_group_detail != '"null"':
            second_group_detail = json.loads(data.second_group_detail)
            params_json['second_group_detail'] =group_detail_sort(second_group_detail)
        else:
            params_json['second_group_detail'] = json.loads(data.second_group_detail)
        if data.anno_type not in ['gene'] and data.anno_type not in  ['nr']:
            params_json['anno_id'] = data.anno_id
            params_json['level_id'] = int(data.level_id)
            level_name = self.level_id(data.level_id)# 此处是level缩写，需要再转义成真实level name
            mongo_data.append(('anno_id', ObjectId(data.anno_id)))
            #mongo_data.append(('level_id', data.level_id))
            if hasattr(data, 'second_level'):
                params_json['second_level'] = data.second_level
               # mongo_data.append(('second_level', data.second_level))
            #main_table_name = 'lefse' +  '_' + data.anno_type.upper() + '_' + level_name + '_' + \
            #          datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            mongo_data.append(('name', main_table_name))
        if data.anno_type in ['nr']:
            params_json['anno_id'] = data.anno_id
            params_json['start_level'] = int(data.start_level)
            params_json['end_level'] = int(data.end_level)
            params_json['domain_type'] = data.domain_type
            mongo_data.append(('anno_id', ObjectId(data.anno_id)))
            #mongo_data.append(('start_level', data.start_level))
            #mongo_data.append(('end_level', data.end_level))
            main_table_name = 'lefse' + '_' + data.anno_type.upper() + '_' +datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            mongo_data.append(('name', main_table_name))
        if data.anno_type in ['gene']:
            main_table_name = 'lefse' + '_' + data.anno_type.upper() + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            mongo_data.append(('name', main_table_name))
        [geneset_table, gene_list] = metagenomic.export_geneset_table(data.geneset_id, data.method)
        options = {
            #'update_info':json.dumps(update_info),
            'params': json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            "group_detail": data.group_detail,
            "second_group_detail": data.second_group_detail,
            "lefse_gname": metagenomic.get_group_name(data.group_id, second_group=data.second_group_detail),
            'group_id': data.group_id,
            'geneset_id': data.geneset_id,
            'group_detail': data.group_detail,
            'geneset_table': geneset_table,
            'strict':int(data.strict),
            'lda_filter':0,
            'group':data.geneset_id
        }
        to_file = []
        #个性化数据库调用
        add_personal = CommCreatAction()
        data,options, params_json,info, to_file = add_personal.add_personal(data,options,params_json,geneset_table,to_file)
        if info!= "":
            return json.dumps(info)
        #个性化数据库用调用
        to_file.append('metagenomic.export_cascading_table_by_detail(group)')
        if data.anno_type not in ['gene'] and data.anno_type not in ['nr']:
            options['lefse_type'] = 'metagenome_func'
            level_name = self.level_id(data.level_id)
            options['level_type'] = level_name
            anno_info = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            options['anno_id'] = data.anno_id
            options['anno_table'] = self.use_s3(anno_info['anno_file'])
            options['lowest_level'] = self.use_s3(anno_info['lowest_level'])
            if hasattr(data, 'second_level'):
                if data.second_level != "":
                    options['level_type_name'] = self.level_convert(data.second_level, data.level_id)
            main_table_name = 'lefse' +  '_' + data.anno_type.upper() + '_' + level_name + '_' + \
                      datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            main_table_name = main_table_name.replace(" ","_")
            mongo_data.append(('name', main_table_name))
        elif data.anno_type in ['nr']:
            options['lefse_type'] = "metagenome_taxon"
            options['start_level'] = int(data.start_level)
            options['end_level'] = int(data.end_level)
            options['domain_type'] = data.domain_type
            anno_info2 = metagenomic.get_anno_info(data.anno_id, data.anno_type)
            options['anno_id'] = data.anno_id
            options['anno_table'] = self.use_s3(anno_info2['anno_file'])
        elif data.anno_type in ['gene']:
            options['lefse_type'] = 'metagenome_func'
            options['gene_list'] = gene_list
        mongo_data.append(('params',json.dumps(params_json,sort_keys=True,separators=(',', ':'))))
        main_table_id =self.metagenomic.insert_main_table("lefse",mongo_data)
        options['main_id'] = str(main_table_id)
        options['main_table_data'] = SON(mongo_data)
        update_info = {str(main_table_id):'lefse'}
        options['update_info'] = json.dumps(update_info)
        self.set_sheet_data(name=task_name, options=options,
                                main_table_name=main_table_name.strip().split("_")[0] + '/' + main_table_name, task_id=task_info['task_id'],
                                project_sn=task_info['project_sn'], module_type=module_type, params=params_json,
                                to_file=to_file)
        task_info = super(LefseMgAction, self).POST()
        if task_info['success']:
           task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

