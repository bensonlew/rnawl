# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

import os
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bacgenome_controller import BacgenomeController
from mainapp.libs.signature import check_sig
from bson import SON


class RefalignDnabacAction(BacgenomeController):
    def __init__(self):
        super(RefalignDnabacAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['task_id', 'specimen_list', 'submit_location', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': 'parameters missing:%s' % argu}
                return json.dumps(info)
        project_sn = self.bacgenome.get_projectsn(data.task_id)
        task_name = 'bacgenome.report.refalign_dnabac'
        module_type = 'workflow'
        params_json = {
            'task_id': data.task_id,
            'specimen_id': data.specimen_list,
            'project_sn': project_sn,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        if hasattr(data,'evalue'):  # zouguanqing  20190326
            params_json['evalue'] = data.evalue
        if not hasattr(data, 'ref') and not hasattr(data, 'ref_path'):
            info = {'success': False, 'info': '参考基因组参数缺少!', 'code':'C1100501'}
            return json.dumps(info)
        elif hasattr(data, 'ref') and hasattr(data, 'ref_path'):
            info = {'success': False, 'info': '不允许同时提供两种形式参考基因组参数缺少!', 'code':'C1100502'}
            return json.dumps(info)
        spe_list = sorted(data.specimen_list.split(","))
        update_info = {}
        main_id = {}
        main_name = {}
        ref_accession = ""
        if hasattr(data, 'ref'):
            params_json['ref'] = data.ref  # 此处ref的应该包含ancession对应的数据库名称
            params_json['ref_name'] = data.ref_name
            ref_accession = data.ref  # 页面参数名称
        elif hasattr(data, 'ref_path'):
            params_json['ref_path'] = data.ref_path
            params_json['file_dir_id'] = data.file_dir_id
            ref_accession = (os.path.basename(data.ref_path).split("."))[0].rsplit("_",1)[0]
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        for i in spe_list:
            mongo_params = {
                'task_id': data.task_id,
                'submit_location': data.submit_location,
                'task_type': int(data.task_type),
                'specimen_id': i,
            }
            mongo_data = [
                ('project_sn', project_sn),
                ('status', 'start'),
                ('task_id', data.task_id),
                ('desc', '正在计算'),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", mongo_params)
            ]
            if hasattr(data, 'ref'):
                mongo_params['ref'] = data.ref  # 数据库名称
                mongo_params['ref_name'] = data.ref_name
                mongo_data.append(('ref', data.ref))
            elif hasattr(data, 'ref_path'):
                mongo_params['ref_path'] = data.ref_path
                mongo_params['file_dir_id'] = data.file_dir_id
                mongo_data.append(('ref_path', data.ref_path))
            mongo_data.append(('params',mongo_params))
            main_table_name = 'AnnoRef_' + i + '_' + ref_accession + \
                              '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            mongo_params = json.dumps(mongo_params, sort_keys=True, separators=(',', ':'))
            mongo_data.append(("specimen_id",i))
            mongo_data.append(("name",main_table_name))
            mongo_data.append(("created_ts",datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            mongo_data.append(("params",mongo_params))
            main_name[i] = str(main_table_name)
            main_table_id = self.bacgenome.insert_main_table('anno_ref', mongo_data)
            update_info[str(main_table_id)] = 'anno_ref'
            main_id[i] = str(main_table_id)
        sample_seqpath = self.bacgenome.get_genefile_bysample(data.task_id, data.specimen_list, type="faa", predict = "gene")
        options = {'sample_seqpath': str(sample_seqpath),
                   'task_id': data.task_id,
                   'update_info': json.dumps(update_info),
                   'main_name': str(main_name),
                   'params': params,
                   'main_id': str(main_id)
                   }
        if hasattr(data,'evalue'):  #zouguanqing 20190326
            options['evalue'] = float(data.evalue)
        if hasattr(data, 'ref'):
            if hasattr(data, 'ref_name') and (data.ref_name != ""):
                options['ref'] = data.ref_name
            else:
                if hasattr(data, 'ref_name') and (data.ref_name == ""):
                    info = {'success': False, 'info': '传入的ref_name为空，请检查ref_namecam参数!', 'code':'C1100503'}
                    return json.dumps(info)
        elif hasattr(data, 'ref_path'):
            options['ref_path'] =  data.ref_path
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Anno_ref/",
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=data.task_id,
                            params=params)
        task_info = super(RefalignDnabacAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': main_id, 'name': main_table_name}}
        return json.dumps(task_info)
