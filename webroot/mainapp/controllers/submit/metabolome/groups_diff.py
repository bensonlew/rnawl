# -*- coding: utf-8 -*-


import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from bson.objectid import ObjectId


class GroupsDiffAction(MetabolomeController):
    def __init__(self):
        super(GroupsDiffAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'group_id', 'group_detail', 'metab_table',
                        'test_method', 'post_hoc', 'coverage']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.groups_diff'
        module_type = 'workflow'
        task_id = data.task_id

        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        if project_type == 'LC':
            task = metabolome.conmon_find_one('sg_task',{'task_id':task_id, "type" : "LC","mix" : "T"})
            if task:
                mix = 'T'  #T or F
            else:
                mix = 'F'
        metab_table_main_id = metabolome.metab_table_main_id(data.metab_table, task_id)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        specimen_group = metabolome.conmon_find_one("specimen_group",{"_id":ObjectId(data.group_id)})
        if specimen_group:
            group_name = specimen_group['group_name']
        else:
            info = {'success':False, "info": "没有找到specimen_group数据库数据"}
            return json.dumps(info)

        group_detail = json.loads(data.group_detail)
        group_num = len(group_detail.keys())
        if group_num < 3:
            info = {'success': False, 'info': '分组数小于3组!'}
            return json.dumps(info)


        if not data.test_method in ["ow", 'kw']:  #ow : 单因素方差分析  kw： kw秩合检验
            info = {'success': False, 'info': '差异检验方法错误!'}
            return json.dumps(info)

        if not data.post_hoc in ["scheffe",'tukeykramer','gameshowell','welchuncorrected']:
          info = {'success': False, 'info': 'Post hoc 检验方法错误'}
          return json.dumps(info)

        if project_type == "LC":
            if mix == 'T':
                metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
                metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
                name = "GroupsDiff_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            else:
                metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'pos')
                metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'pos')
                name = "GroupsDiff_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                metab_neg_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'neg')
                metab_neg_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'neg')

        else:
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
            name = "GroupsDiff_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "test_method": data.test_method,
            "post_hoc" : data.post_hoc,
            "coverage" :data.coverage,
            "task_id" : data.task_id
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("metab_table_id", ObjectId(data.metab_table)),
            ("metab_table_main_id", ObjectId(metab_table_main_id))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))


        options = {
            "metab_table": metab_table_path,
            "metab_desc": metab_desc_path,
            "group": data.group_id,
            "group_detail": data.group_detail,
            "test_method": data.test_method,
            "post_hoc" : data.post_hoc,
            "coverage" : data.coverage,
            "group_name" : group_name,
            "name": name
        }
        if project_type == "LC":
            if mix == 'T':
                options["table_type"] = 'mix'
            else:
                options["table_type"] = 'pos,neg'
                options['metab_table_neg'] = metab_neg_table_path
                options['metab_desc_neg'] = metab_neg_desc_path
        else:
            options["table_type"] = 'pos'

        main_table_id = metabolome.insert_main_table("groups_diff", mongo_data)
        update_info = {str(main_table_id): "groups_diff"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_id"] = str(main_table_id)
        to_file = []
        to_file.append('metabolome.export_group_by_detail(group)')
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = '4.ExpDiff/03.GroupsDiff/' + name
        else:
            m_table_name = name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(GroupsDiffAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

