# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0525

import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from bson.objectid import ObjectId


class ExpDiffAction(MetabolomeController):
    def __init__(self):
        super(ExpDiffAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'group_id', 'group_detail', 'metab_table',
                        'diff_group_id', 'diff_group_detail', 'test_method', 'tail', 'pca_trans', 'pca_confidence',
                        'plsda_trans', 'plsda_confidence', 'plsda_replace', 'oplsda_trans', 'oplsda_confidence',
                        'oplsda_replace']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.exp_diff'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        metab_table_main_id = metabolome.metab_table_main_id(data.metab_table, task_id)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        group_detail = json.loads(data.group_detail)
        select_samples = self.ext_samples(group_detail)
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
        # if data.test_method in ["t-test", "welch"]:
        #     if not len(group_detail) == grouo_len3:
        #         info = {'success': False, 'info': 't-test检验和welch差异检验方法时每组样本数需不少于3!', 'code':'C2300503'}
        #         return json.dumps(info)
        # elif data.test_method == "wilcox":
        #     if not len(group_detail) == group_len4:
        #         info = {'success': False, 'info': 'wilcox差异检验方法时每组样本数需不少于4!', 'code':'C2300504'}
        #         return json.dumps(info)
        if data.test_method == ["chi_sq", "fisher"]:
            if not len(group_detail) == group_len1:
                info = {'success': False, 'info': 'chi_sq，fisher差异检验方法时每组样本数应为1!', 'code':'C2300505'}
                return json.dumps(info)
        for each in [data.pca_trans, data.plsda_trans, data.oplsda_trans]:
            if not each in ["UV", "Par", "Ctr", ""]:
                info = {'success': False, 'info': '数据转换方法错误!', 'code':'C2300506'}
                return json.dumps(info)
        for each in [data.pca_confidence, data.plsda_confidence, data.oplsda_confidence]:
            if not 0 < float(each) < 1:
                info = {'success': False, 'info': '置信度范围为0-1!', 'code':'C2300507'}
                return json.dumps(info)
        for each in [data.plsda_replace, data.oplsda_replace]:
            if not 19 < int(each) < 1001:
                info = {'success': False, 'info': '置换次数范围为20-1000!', 'code':'C2300508'}
                return json.dumps(info)
        if not data.test_method in ["t-test", 'welch', 'wilcox']:
            info = {'success': False, 'info': '差异检验方法错误!', 'code':'C2300509'}
            return json.dumps(info)
        if not data.tail in ["two-tailed", "left-tailed", "right-tailed"]:
            info = {'success': False, 'info': '单双尾检验方法错误!', 'code':'C2300510'}
            return json.dumps(info)
        if project_type == "LC":
            if not hasattr(data, "table_type"):
                info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code':'C2300511'}
                return json.dumps(info)
            else:
                metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, data.table_type)
                metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, data.table_type)
                name = "ExpDiff_" + data.table_type + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
            name = "ExpDiff_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            data.table_type = "pos"
        trans_data = ";".join([data.pca_trans, data.plsda_trans, data.oplsda_trans])
        confidence = ";".join([data.pca_confidence, data.plsda_confidence, data.oplsda_confidence])
        replacement = ";".join(["0", data.plsda_replace, data.oplsda_replace])
        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "diff_group_id": data.diff_group_id,
            "diff_group_detail": json.loads(data.diff_group_detail),
            "test_method": data.test_method,
            "tail": data.tail,
            "pca_trans": data.pca_trans,
            "pca_confidence": data.pca_confidence,
            "plsda_trans": data.plsda_trans,
            "plsda_confidence": data.plsda_confidence,
            "plsda_replace": data.plsda_replace,
            "oplsda_trans": data.oplsda_trans,
            "oplsda_confidence": data.oplsda_confidence,
            "oplsda_replace": data.oplsda_replace,
            "task_id": task_id,
            "table_type": data.table_type
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("metab_table_id", ObjectId(data.metab_table)),
            ("metab_table_main_id", ObjectId(metab_table_main_id)),
            ("table_type", data.table_type),
            ("diff_dir", "")
        ]
        # prepare option for workflow
        tmp_diff = json.loads(data.diff_group_detail)
        diff_group_name = ";".join(tmp_diff)
        options = {
            "metab_table": metab_table_path,
            "metab_desc": metab_desc_path,
            "group": data.group_id,
            "group_name": diff_group_name,
            "group_detail": data.group_detail,
            "mul_type": "pca;plsda;oplsda",
            "confidence": confidence,
            "perm": replacement,
            "data_trans": trans_data,
            "test_method": data.test_method,
            "side_type": data.tail,
            "name": name
        }
        if project_type == "LC":
            options["table_type"] = data.table_type
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("exp_diff", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "exp_diff"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        to_file = []
        to_file.append('metabolome.export_group_by_detail(group)')
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = '4.ExpDiff/01.TwoGroupExpDiff/' + name
        else:
            m_table_name = name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(ExpDiffAction, self).POST()
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
