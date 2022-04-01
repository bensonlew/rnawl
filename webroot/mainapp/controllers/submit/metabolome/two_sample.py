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


class TwoSampleAction(MetabolomeController):
    def __init__(self):
        super(TwoSampleAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'group_id', 'group_detail', 'metab_table',
                        'diff_group_id', 'diff_group_detail', 'test_method', 'tail', 'correct']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.two_sample'
        module_type = 'workflow'
        task_id = data.task_id

        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        if project_type == 'LC':
            task = metabolome.conmon_find_one('sg_task',{'task_id':task_id, "type" : "LC","mix" : "T"})
            if task:
                mix = 'T'  #T or F
            else:
                mix = 'F'
        ##metab_table_main_id = metabolome.metab_table_main_id(data.metab_table, task_id)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        group_detail = json.loads(data.group_detail)
        diff_group_detail = json.loads(data.diff_group_detail)
        for eachvs in diff_group_detail:
            vsgroups= eachvs.split("_vs_")
            if len(vsgroups) > 2:
                info = {'success': False, 'info': '分组名不能带_vs_,与差异对比组冲突', 'code':'C2300501'}
                return json.dumps(info)
            for each in vsgroups:
                if not each in group_detail.keys():
                    variables = []
                    variables.append(each)
                    info = {'success': False, 'info': '差异组group-{}不在分组方案中，请重新选择!'.format(each), 'code':'C2300502', 'variables':variables}
                    return json.dumps(info)
        group_len1 = 0
        grouo_len3 = 0
        group_len4 = 0
        for eachgroup in group_detail.keys():
            sam_len = len(group_detail[eachgroup])
            if sam_len == 1:
                group_len1 += 1
            if sam_len >= 3:
                grouo_len3 += 1
            if sam_len >= 4:
                group_len4 += 1

        if data.test_method == ["chiq", "fisher"]:
            if not len(group_detail) == group_len1:
                info = {'success': False, 'info': 'chiq，fisher差异检验方法时每组样本数应为1!', 'code':'C2300505'}
                return json.dumps(info)
        if not data.test_method in ["chiq", 'fisher']:
            info = {'success': False, 'info': '差异检验方法错误!', 'code':'C2300509'}
            return json.dumps(info)
        if not data.tail in ["two-tailed", "left-tailed", "right-tailed"]:    #["two.side","less","greater"]: ##
            info = {'success': False, 'info': '单双尾检验方法错误!', 'code':'C2300510'}
            return json.dumps(info)
        #if not data.correct in [""]:
        #   info = {'success': False, 'info': '多重检验校正方法错误'}
        #   return json.dumps(info)

        if project_type == "LC":
            # if not hasattr(data, "table_type"):
            #     info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code':'C2300511'}
            #     return json.dumps(info)
            # else:
            if mix == 'T':
                metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
                metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
                name = "TwoSample_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            else:
                metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'pos')
                metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'pos')
                name = "TwoSample_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                metab_neg_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'neg')
                metab_neg_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'neg')

        else:
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
            name = "TwoSample_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            #data.table_type = "pos"

        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "diff_group_id": data.diff_group_id,
            "diff_group_detail": json.loads(data.diff_group_detail),
            "test_method": data.test_method,
            "correct" :data.correct,
            "task_id": task_id,
            "tail": data.tail
        }
        if data.test_method == 'fisher':
            params_json['tail'] = data.tail
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("metab_table_id", ObjectId(data.metab_table)),
            #("metab_table_main_id", ObjectId(metab_table_main_id)),
            ("version","3.0")
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        # prepare option for workflow
        tmp_diff = json.loads(data.diff_group_detail)
        diff_group_name = ";".join(tmp_diff)
        options = {
            "metab_table": metab_table_path,
            "metab_desc": metab_desc_path,
            "group": data.group_id,
            "diff_group_name": diff_group_name,
            "group_detail": data.group_detail,
            "test_method": data.test_method,
            "tail": data.tail,
            "correct" : data.correct,  #######
            # "name": name
        }
        if data.test_method == 'fisher':
            options['tail'] = data.tail
        if project_type == "LC":
            if mix == 'T':
                options["table_type"] = 'mix'
            else:
                options["table_type"] = 'pos,neg'
                options['metab_table_neg'] = metab_neg_table_path
                options['metab_desc_neg'] = metab_neg_desc_path
        else:
            options["table_type"] = 'pos'

        main_table_id = metabolome.insert_main_table("diff_sample", mongo_data)
        update_info = {str(main_table_id): "diff_sample"}
        options["update_info"] = json.dumps(update_info)
        ###options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        to_file = []
        to_file.append('metabolome.export_group_by_detail(group)')
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            # options["save_pdf"] = 1
            m_table_name = '4.ExpDiff/02.TwoSamExpDiff/' + name
        else:
            m_table_name = name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(TwoSampleAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

    def ext_samples(self, group_detail):
        samples = []
        for each in group_detail.values():
            for i in each:
                if i not in samples:
                    samples.append(i)
        return samples
