# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0708

import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
#from mbio.api.to_file.metabolome import *
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from bson.objectid import ObjectId


class MetabsetVipAction(MetabolomeController):
    def __init__(self):
        super(MetabsetVipAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'diff_table', 'group_id', 'group_detail',
                        'metab_table', 'metab_set', 'group_method', 'diff_group_id', 'diff_group_detail',
                        'metab_cluster_method','vip', 'vip_from','top_vip', 'scale']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing:%s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.metabset_vip'
        module_type = 'workflow'
        task_id = data.task_id
        task_info = metabolome.conmon_find_one('sg_task',{"task_id": task_id})
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        diff_dir_ori = metabolome.get_diff_dir(data.diff_table, task_id)
        sp_diff_dir = diff_dir_ori.split(',')
        has_neg_dir = False
        if len(sp_diff_dir) == 2:
            diff_dir = sp_diff_dir[0]
            diff_dir2 = sp_diff_dir[1]
            has_neg_dir = True
        else:
            diff_dir = diff_dir_ori


        # if hasattr(data,'table_type'):
        #     table_type = data.table_type
        if task_info['type'] == 'GC':
            table_type = 'pos'
        else:
            if  'mix' in task_info and task_info['mix'] == 'T':
                table_type = 'mix'
            else:
                if hasattr(data,'table_type'):
                    table_type = data.table_type


        if table_name == "raw":
            if table_type == "mix":
                info = {"success": False, "info": "原始表没有合并表,请选择正确参数!", 'code':'C2301601'}
                return json.dumps(info)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        set_name_r = metabolome.conmon_find_one('metab_set',{"_id":ObjectId(data.metab_set)})
        if set_name_r:
            set_name_pls = set_name_r['name']+'_'
        else:
            set_name_pls = ''
        name = "MetabsetVip_" + set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if project_type == "LC":
            # if not hasattr(data, "table_type"):
            #     info = {"success": False, "info": "原始表没有合并表,请选择正确参数!", 'code':'C2301601'}
            #     return json.dumps(info)
            # else:
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
        if not data.metab_cluster_method in ["hierarchy", "kmeans", "none"]:
            info = {"success": False, "info": "原始表没有合并表,请选择正确参数!", 'code':'C2301601'}
            return json.dumps(info)
        if data.metab_cluster_method == "hierarchy":
            if not data.metab_dist or not data.metab_cluster:
                info = {"success": False, "info": "hierarchy聚类方法时必须输入距离算法和层级聚类方式!", 'code':'C2301604'}
                return json.dumps(info)
        if data.metab_cluster_method == "kmeans":
            if not data.metab_dist or not data.n_cluster:
                info = {"success": False, "info": "kmeans聚类方法时必须输入距离算法和代谢物子聚类数目!", 'code':'C2301605'}
                return json.dumps(info)
            if not int(data.n_cluster) > 1:
                info = {"success": False, "info": "代谢物子聚类数目必须大于等于2!", 'code':'C2301606'}
                return json.dumps(info)
        '''
        if data.metab_cluster_method == "none":
            if data.metab_dist or data.n_cluster or data.metab_cluster:
                raise OptionError("无聚类方法时不需要距离算法和代谢物子聚类数目！")
        '''
        if data.vip_from not in ["plsda", "oplsda"]:
            info = {'success': False, 'info': 'VIP值来源错误!', 'code':'C2301607'}
            return json.dumps(info)
        if not float(data.vip) > 0:
            info = {'success': False, 'info': 'VIP值为正整数!', 'code':'C2301608'}
            return json.dumps(info)
        if not 0 < int(data.top_vip) <= 300:
            info = {'success': False, 'info': 'top vip范围为1-300!', 'code':'C2301609'}
            return json.dumps(info)
        params_json = {
            "metab_table": data.metab_table,
            "group_id": data.group_id,
            "group_detail": json.loads(data.group_detail),
            "diff_table": data.diff_table,
            "diff_group_id": data.diff_group_id,
            "diff_group_detail": json.loads(data.diff_group_detail),
            "metab_set": data.metab_set,
            "group_method": data.group_method,
            "metab_cluster_method": data.metab_cluster_method,
            "task_id": task_id,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "vip_from":data.vip_from,
            "vip": data.vip,
            "top_vip": data.top_vip
        }
        if hasattr(data,'table_type'):
            params_json['table_type'] = data.table_type

        if data.group_method != "none":
            diff_detail = {}
            tmp_detail = json.loads(data.group_detail)
            for each in tmp_detail.keys():
                diff_detail[each] = [each]
        else:
            diff_detail = json.loads(data.group_detail)
        diff_detail = json.dumps(diff_detail)
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("main_id", ""),
            ("diff_detail",diff_detail),
            ("diff_id",ObjectId(data.diff_table)),
            #("table_type",table_type),
            ("metab_set",ObjectId(data.metab_set))
        ]
        if hasattr(data,'table_type'):
            mongo_data.append(("table_type", data.table_type))

        # prepare option for workflow
        if hasattr(data, "metab_dist") and data.metab_dist == "manhattan":
            metab_dist = "cityblock"
        elif hasattr(data, "metab_dist"):
            metab_dist = data.metab_dist
        options = {
            "diff_dir": diff_dir,
            "vip_type": data.vip_from,
            "vip_cut": data.vip,
            "vip_top": data.top_vip,
            "mct": data.metab_cluster_method,
            "group_method":data.group_method,
            "group": data.group_id,
            "group_detail": data.group_detail,
            #"group_name": data.diff_group_detail,
            "metab_set_table": data.metab_set,
            "metab_set_id": data.metab_set,
            "metab_desc": metab_desc_path,
            "metab_table": metab_table_path,
            "name": name
            #"table_type": table_type
        }
        if has_neg_dir:
            options['diff_dir2'] = diff_dir2

        if data.metab_cluster_method == "kmeans":
            options["n_cluster"]  = data.n_cluster
            options["mcd"]  = data.metab_dist
            params_json["n_cluster"] = int(data.n_cluster)
            params_json["metab_dist"] = data.metab_dist
        if data.metab_cluster_method == "hierarchy":
            options["mcm"]  = data.metab_cluster
            options["mcd"]  = metab_dist
            params_json["metab_cluster"] = data.metab_cluster
            params_json["metab_dist"] = data.metab_dist
        if hasattr(data, "scale"):
            params_json["scale"] = data.scale
            options["scale"]  = data.scale
        to_file = []
        to_file.append('metabolome.export_group_by_detail(group)')
        to_file.append('metabolome.export_metab_set(metab_set_table)')
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("metabset_vip", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "metabset_vip"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = "5.Metabset/03.MetabsetVip/" + name
        else:
            m_table_name = "Metabset/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn,to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(MetabsetVipAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        self.Metabolome.insert_set_info(data.metab_set, "metabset_vip", main_table_id)
        return json.dumps(task_info)

