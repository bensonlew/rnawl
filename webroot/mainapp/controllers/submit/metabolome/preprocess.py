# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.controllers.project.metabolome_controller import MetabolomeController


class PreprocessAction(MetabolomeController):
    def __init__(self):
        super(PreprocessAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        # default_argu = ['metab_table', 'group', 'fillna', 'rsd', 'norm', 'scale', 'submit_location']
        default_argu = ['group', 'fillna', 'rsd', 'norm', 'scale', 'submit_location','fill_percent','task_id']
        task_info = self.Metabolome.conmon_find_one('sg_task',{'task_id':data.task_id})
        if 'version' in task_info:
            task_version = task_info['version']
        else:
            task_version = '1.0'

        if task_version != '1.0':
            default_argu.append('fill_type')

        for argu in default_argu:
            if not hasattr(data, argu):
                info = {"success": False, "info": "parameters missing:%s" % argu}
                return json.dumps(info)

        metab_table_path = self.Metabolome.get_metab_path(data.task_id)
        if len(metab_table_path.split(',')) == 1:
            project_type = "GC"
        else:
            project_type = "LC"
        if data.fillna not in ["min", "median", "mean", "rf", "none"]:
            variables = []
            variables.append(data.fillna)
            info = {"success": False, "info": "缺失值方法错误：%s" % data.fillna, 'code':'C2301701', 'variables':variables}
            return json.dumps(info)
        if data.rsd not in  ["none",""] :
            try:
                rsd_value = int(data.rsd.rstrip("%"))
            except:
                variables = []
                variables.append(data.rsd)
                info = {"success": False, "info": "质控参数填写错误：%s" % data.rsd, 'code':'C2301702', 'variables':variables}
                return json.dumps(info)
            if rsd_value < 25 or rsd_value > 30:
                variables = []
                variables.append(rsd_value)
                info = {"success": False, "info": "质控参数必须在25-30之间：%s" % rsd_value, 'code':'C2301703', 'variables':variables}
                return json.dumps(info)
        else:
            rsd_value = data.rsd
        if data.norm not in ["median", "mean", "sum", "sample", "inner", "none"]:
            info = {"success": False, "info": "所选归一化方法不在范围内", 'code':'C2301704'}
            return json.dumps(info)
        if data.norm == "sample" and not hasattr(data, "sample_name"):
            info = {"success": False, "info": "归一化方法为指定样本时必须选择样本", 'code':'C2301705'}
            return json.dumps(info)
        if data.norm == "inner" and not hasattr(data, "inner_ref"):
            info = {"success": False, "info": "归一化方法为内参时必须选择内参代谢物", 'code':'C2301706'}
            return json.dumps(info)
        if data.scale not in ["log2", "log10", "none", "defined"]:
            variables = []
            variables.append(data.scale)
            info = {"success": False, "info": "Log取值方法不正确：%s" % data.scale, 'code':'C2301707', 'variables':variables}
            return json.dumps(info)
        if data.scale == "defined" and not hasattr(data, "log"):
            info = {"success": False, "info": "自定义Log取值时必须填入数值", 'code':'C2301708'}
            return json.dumps(info)
        # elif data.scale == "undefined":
        #     scale_value = "log" + data.log
        # else:
        #     scale_value = data.scale
        task_name = 'metabolome.report.preprocess'
        module_type = 'workflow'

        params_json = {
            'fillna': data.fillna,
            # 'rsd': data.rsd,
            'rsd': str(rsd_value) if rsd_value !="none" else "",  # change @ 20180710
            'norm': data.norm,
            'scale': data.scale,
            'submit_location': data.submit_location,
            'task_id': data.task_id,
            "task_type": int(data.task_type),
            "fill_percent" : data.fill_percent
        }
        if task_version != "1.0":
            fill_type = data.fill_type
            params_json['fill_type'] = fill_type
        else:
            fill_type = 'all'


        task_info = self.Metabolome.get_task_info(data.task_id)
        mongo_data = [
            # ('project_sn', data.project_sn),
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('type', project_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('is_raw', 0),
            ('setted', 1),  # 目前所有表都可以用于代谢集分析

        ]
        if task_version != "1.0":
            mongo_data.append(('version', '3.0'))
        else:
            mongo_data.append(('version', '1.0'))

        if data.norm == 'sample':
            params_json['sample_name'] = data.sample_name
        if data.norm == 'inner':
            params_json['inner_ref'] = data.inner_ref
        if data.scale == 'defined':
            params_json['log'] = data.log
        main_table_name = "Metab_table_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data.append(('name', main_table_name))
        metabolome = Metabolome()
        ### rere path modify
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            rere_table_path = self._get_file_path('/pos', data.task_id, '1.Preprocess/' + main_table_name)
            if project_type == 'LC':
                rere_table_path += ',' + self._get_file_path('/neg', data.task_id, '1.Preprocess/' + main_table_name)
                rere_table_path += ',' + self._get_file_path('/mix', data.task_id, '1.Preprocess/' + main_table_name)
        else:
            rere_table_path = self._get_file_path('/pos', data.task_id, 'Preprocess/' + main_table_name)
            if project_type == 'LC':
                rere_table_path += ',' + self._get_file_path('/neg', data.task_id, 'Preprocess/' + main_table_name)
                rere_table_path += ',' + self._get_file_path('/mix', data.task_id, 'Preprocess/' + main_table_name)
        print rere_table_path
        mongo_data.append(('table_path', rere_table_path))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = metabolome.insert_main_table('exp', mongo_data)
        metabolome.insert_main_id('exp', main_table_id)
        update_info = {str(main_table_id): 'exp'}
        if rsd_value not in ["none",'']:                     # change @ 2018080920 by shaohua.yuan
            use_rsd = str(rsd_value) + "%"
        else:
            use_rsd = "none"
        options = {
            'update_info': json.dumps(update_info),
            'group_table': data.group,
            'ana_method': project_type,
            'fillna': data.fillna,
            # 'rsd': data.rsd,
            'rsd': use_rsd,  # change @ 2018080920 by shaohua.yuan
            'norm': data.norm,
            'scale': data.scale,
            'main_table_id': str(main_table_id),
            'name': main_table_name,
            'rm_nan' : float(data.fill_percent) , #20190605
            'fill_type' : fill_type,  # v3 20200305
            'task_version' : task_version
        }
        if project_type == "GC":
            options['pos_rout'] = metab_table_path
        elif project_type == "LC":
            table_dir = metab_table_path.split(",")
            options['pos_rout'] = table_dir[0]
            options['neg_rout'] = table_dir[1]
        if hasattr(data, "sample_name"):
            options['sample_name'] = data.sample_name
        elif hasattr(data, "inner_ref"):
            options['inner_ref'] = data.inner_ref
        if hasattr(data, "log"):
            options['log'] = data.log

        if task_version != '1.0':  #1.0时 ，不画v3版新增的cv图和统计
            exp_find = self.Metabolome.conmon_find_one('exp',{'task_id':data.task_id,'name':'raw'})
            if exp_find:
                raw_exp_id = exp_find["main_id"]
                exp_cv_check = self.Metabolome.conmon_find_one('exp_cv', {'exp_id': raw_exp_id})
                if not exp_cv_check:
                    options['run_cv_raw'] = 'True'
            else:
                info = {'success':False, 'info': '该项目没有raw表'}
                return json.dumps(info)

            options['raw_exp_id'] = str(raw_exp_id)

        to_file = ['metabolome.export_group_table(group_table)']
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = '1.Preprocess/' + main_table_name
        else:
            m_table_name = 'Preprocess/' + main_table_name
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,
                            # project_sn=data.project_sn,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        post_info = super(PreprocessAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(post_info)
