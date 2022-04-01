# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0529

import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from bson import ObjectId


class MetabsetCorrAction(MetabolomeController):
    def __init__(self):
        super(MetabsetCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ["submit_location", 'task_type', 'task_id', 'group_id', 'group_detail', 'metab_table',
                        'metab_set', 'coefficient', 'cluster_type', 'top']
        # check arg
        for arg in default_argu:
            print arg
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.metabset_corr'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        if table_name == "raw":
            if data.table_type == "mix":
                info = {"success": False, "info": "原始表没有合并表,请选择正确参数!", 'code':'C2300801'}
                return json.dumps(info)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        group_detail = json.loads(data.group_detail)
        select_samples = self.ext_samples(group_detail)
        if len(select_samples) < 2:
            info = {"success": False, "info": "至少两个样本!"}
            return json.dumps(info)
        set_name_r = metabolome.conmon_find_one('metab_set',{"_id":ObjectId(data.metab_set)})
        if set_name_r:
            set_name_pls = set_name_r['name']+'_'
        else:
            set_name_pls = ''
        name = "MetabsetCorr_"+set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "coefficient": data.coefficient,
            "cluster_type": data.cluster_type,
            "top": int(data.top),
            "metab_set": data.metab_set,
            "task_id": task_id
        }
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("metab", ""),
            ("tree", ""),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        if project_type == "LC":
            # if not hasattr(data, "table_type"):
            #     info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code':'C2300803'}
            #     return json.dumps(info)

            if hasattr(data, "table_type"):
                params_json["table_type"] = data.table_type
                metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, data.table_type)
                metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, data.table_type)
                table_type = data.table_type
            else:
                metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
                metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
                table_type = 'mix'

        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
            params_json["table_type"] = 'pos'
            table_type = 'pos'

        print metab_table_path
        mongo_data.append(("table_type",table_type))
        mongo_data.append(("metab_set",ObjectId(data.metab_set)))
        # prepare option for workflow
        if hasattr(data, "dist") and data.dist == "manhattan":
            dist = "cityblock"
        elif hasattr(data, "dist"):
            dist = data.dist
        options = {
            "metab_table": metab_table_path,
            "metab_desc": metab_desc_path,
            "samples": ",".join(select_samples),
            "cluster_type": data.cluster_type,
            "coefficient": data.coefficient,
            "top": int(data.top),
            "metab_set_table": data.metab_set,
            "metab_set_id": data.metab_set,
            "name": name
            #"group_id": data.group_id,
            #"group_table": data.group_id

        }
        if not data.coefficient in ["spearman", "pearson", "kendall"]:
            info = {'success': False, 'info': '相关性算法类型错误!', 'code':'C2300804'}
            return json.dumps(info)
        if int(data.top) < 2 or int(data.top) > 300:
            info = {'success': False, 'info': 'top代谢物个数只能在2-300!', 'code':'C2300805'}
            return json.dumps(info)
        if data.cluster_type != "no":
            if not hasattr(data, "dist"):
                info = {'success': False, 'info': '进行聚类时必须输入距离算法!', 'code':'C2300806'}
                return json.dumps(info)
            else:
                options["dist"] = dist
                params_json["dist"] = data.dist
        if data.cluster_type == "kmeans":
            if not hasattr(data, "n_cluster"):
                info = {'success': False, 'info': '聚类算法为kmean时必须输入代谢物子聚类数目!', 'code':'C2300807'}
                return json.dumps(info)
            else:
                options["n_cluster"] = int(data.n_cluster)
                params_json["n_cluster"] = int(data.n_cluster)
            if int(data.n_cluster) > int(data.top):
                info = {'success': False, 'info': '代谢物子聚类数目必须小于代谢物个数!', 'code':'C2300808'}
                return json.dumps(info)
        elif data.cluster_type == "hierarchy":
            if not hasattr(data, "cluster"):
                info = {'success': False, 'info': '聚类算法为hierarchy时必须输入层级聚类方式!', 'code':'C2300809'}
                return json.dumps(info)
            else:
                options["cluster"] = data.cluster
                params_json["cluster"] = data.cluster
            if not data.cluster in ["complete", "average", "single"]:
                info = {'success': False, 'info': '层次聚类算法错误!', 'code':'C2300810'}
                return json.dumps(info)
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("metabset_corr", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "metabset_corr"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        to_file = []
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        #to_file.append('metabolome.export_group_table(group_table)')
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = "5.Metabset/04.MetabsetCorr/" + name
        else:
            m_table_name = "Metabset/"+name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(MetabsetCorrAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        self.Metabolome.insert_set_info(data.metab_set, "metabset_corr", main_table_id)
        return json.dumps(task_info)

    def ext_samples(self, group_detail):
        samples = []
        for each in group_detail.values():
            for i in each:
                if i not in samples:
                    samples.append(i)
        return samples
