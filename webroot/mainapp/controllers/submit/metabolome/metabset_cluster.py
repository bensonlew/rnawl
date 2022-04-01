# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0525

import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
#from mbio.api.to_file.metabolome import *
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from bson import ObjectId



class MetabsetClusterAction(MetabolomeController):
    def __init__(self):
        super(MetabsetClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'metab_set', 'group_id', 'group_detail',
                        'metab_table', 'group_method', 'metab_cluster_method', 'sam_cluster_method','top_meta','scale']
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)

        # if data.metab_cluster_method =='no' and data.sam_cluster_method == 'no':
        #     info = {"success": False, "info": "metab_cluster_method is no and sam_cluster_method is no." }
        #     return info

        metabolome = Metabolome()
        task_name = 'metabolome.report.metabset_cluster'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        # if table_name == "raw":
        #     if data.table_type == "mix":
        #         info = {"success": False, "info": "原始表没有合并表,请选择正确参数!", 'code':'C2300701'}
        #         return json.dumps(info)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        group_detail = json.loads(data.group_detail)
        set_name_r = metabolome.conmon_find_one('metab_set',{"_id":ObjectId(data.metab_set)})
        if set_name_r:
            set_name_pls = set_name_r['name']+'_'
        else:
            set_name_pls = ''
        name = "MetabsetCluster_"+ set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "metab_table": data.metab_table,
            "group_id": data.group_id,
            "group_detail": group_detail,
            "metab_cluster_method": data.metab_cluster_method,
            "sam_cluster_method": data.sam_cluster_method,
            "group_method": data.group_method,
            "task_id": task_id,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "top_meta" : int(data.top_meta),
            "metab_set": data.metab_set,
        }
        if project_type == "LC":
            # if not hasattr(data, "table_type"):
            #     info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code':'C2300702'}
            #     return json.dumps(info)
            # else:
            table_type = 'mix'
            ###params_json["table_type"] = 'mix'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id,'mix')
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
            table_type = 'pos'
            ###params_json["table_type"] = 'pos'
        if hasattr(data,'table_type'):
            params_json['table_type'] = data.table_type

        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("specimen", ""),
            ("metab_tree", ""),
            ("sample_tree", ""),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("table_type", table_type),
            ("metab_set",ObjectId(data.metab_set))

        ]
        # prepare option for workflow
        '''
        if len(group_detail.values()) == 1:
            if data.sam_cluster_method != "no":
                info = {"success": False, "info": "样本个数为1聚类方式应选无！", 'code':'C2300703'}
                return json.dumps(info)
        '''
        if not data.group_method in ["sum", "average", "median","no"]:
            info = {"success": False, "info": "分组计算方法只能为sum, average, median,no!", 'code':'C2300703'}
            return json.dumps(info)
        if not data.metab_cluster_method in ["hierarchy", "kmeans", "no"]:
            info = {"success": False, "info": "代谢物聚类算法只能为hierarchy,kmeans,无！", 'code':'C2300704'}
            return json.dumps(info)
        if not data.sam_cluster_method in ["hierarchy", "no"]:
            info = {"success": False, "info": "样本聚类算法只能为hierarchy,无！", 'code':'C2300705'}
            return json.dumps(info)
        if data.metab_cluster_method == "no":
            data.metab_cluster_method = ""
        if data.sam_cluster_method  == "no":
            data.sam_cluster_method  = ""
        if hasattr(data, "sam_dist") and data.sam_dist == "manhattan":
            sam_dist = "cityblock"
        elif hasattr(data, "sam_dist"):
            sam_dist = data.sam_dist
        if hasattr(data, "metab_dist") and data.metab_dist == "manhattan":
            metab_dist = "cityblock"
        elif hasattr(data, "metab_dist"):
            metab_dist = data.metab_dist
        options = {
            "metab_table": metab_table_path,
            "metab_desc": metab_desc_path,
            "metab_cluster_method": data.metab_cluster_method,
            "sam_cluster_method": data.sam_cluster_method,
            "group_method": data.group_method,
            "metab_set_table": data.metab_set,
            "metab_set_id": data.metab_set,
            "group_table": data.group_id,
            "top_meta": data.top_meta,
            "group_detail": data.group_detail,
            "name": name
        }
        if data.metab_cluster_method == "hierarchy":
            if not hasattr(data, "metab_dist") or not hasattr(data, "metab_cluster"):
                info = {"success": False, "info": "代谢物hierarchy聚类方法时必须输入距离算法和层级聚类方式！", 'code':'C2300706'}
                return json.dumps(info)
            options["metab_dist"] = metab_dist
            options["metab_cluster"] = data.metab_cluster
            params_json["metab_dist"] = data.metab_dist
            params_json["metab_cluster"] = data.metab_cluster
        if data.sam_cluster_method == "hierarchy":
            if not hasattr(data, "sam_dist") or not hasattr(data, "sam_cluster"):
                info = {"success": False, "info": "样本hierarchy聚类方法时必须输入距离算法和层级聚类方式！", 'code':'C2300707'}
                return json.dumps(info)
            options["sam_dist"] = sam_dist
            options["sam_cluster"] = data.sam_cluster
            params_json["sam_dist"] = data.sam_dist
            params_json["sam_cluster"] = data.sam_cluster
        if data.metab_cluster_method == "kmeans":
            if not hasattr(data, "metab_dist") or not hasattr(data, "n_cluster"):
                info = {"success": False, "info": "kmeans聚类方法时必须输入距离算法和代谢物子聚类数目！", 'code':'C2300708'}
                return json.dumps(info)
            if not int(data.n_cluster) > 1:
                info = {"success": False, "info": "代谢物子聚类数目必须大于等于2！", 'code':'C2300709'}
                return json.dumps(info)
            options["metab_dist"] = data.metab_dist
            params_json["metab_dist"] = data.metab_dist
        if hasattr(data,'n_cluster'):
            options["n_cluster"] = data.n_cluster
            params_json["n_cluster"] = int(data.n_cluster)
        #if not  1 < int(data.top_meta) <= 300:
        #    info = {"success": False, "info": "Top代谢物范围2-300！！", 'code':'C2300710'}
        #    return json.dumps(info)
        if hasattr(data, "scale"):
            params_json["scale"] = data.scale
            options["scale"]  = data.scale
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("metabset_cluster", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "metabset_cluster"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        to_file = []
        to_file.append('metabolome.export_group_by_detail(group_table)')
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = "5.Metabset/02.MetabsetCluster/" + name
        else:
            m_table_name = "Metabset/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn,to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(MetabsetClusterAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        self.Metabolome.insert_set_info(data.metab_set, "metabset_cluster", main_table_id)
        return json.dumps(task_info)

