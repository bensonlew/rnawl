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


class ExpCorrAction(MetabolomeController):
    def __init__(self):
        super(ExpCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ["submit_location", 'task_type', 'task_id', 'group_id', 'group_detail', 'metab_table', 'dist',
                        'coefficient', 'cluster',"transform"]
        # check arg
        print(data)
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)
        if data.transform not in ["None","none","UV","Ctr","Par"]:
            info = {"success":False, "info":"transform Vaule must be in [UV, Ctr, Par, None, none"}
            return info
        metabolome = Metabolome()
        task_name = 'metabolome.report.exp_corr'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)
        group_detail = json.loads(data.group_detail)
        print group_detail
        select_samples = self.ext_samples(group_detail)
        if len(select_samples) < 2:
            info = {"success": False, "info": "样本数必须大于1!", 'code':'C2300401'}
            return json.dumps(info)
        name = "ExpCorr_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "group_id": data.group_id,
            "group_detail": group_detail,
            "dist": data.dist,
            "coefficient": data.coefficient,
            "cluster": data.cluster,
            "task_id": task_id,
            "transform":data.transform
        }
        if project_type == "LC":
            if hasattr(data, "table_type"):
                table_type = data.table_type
                params_json["table_type"] = table_type
            else:
                table_type = 'mix'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, table_type)
        elif project_type == "GC":
            table_type = 'pos'
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
        if not data.coefficient in ["spearman", "pearson", "kendall"]:
            info = {'success': False, 'info': '相关性算法类型错误!', 'code':'C2300403'}
            return json.dumps(info)
        if not data.cluster in ["complete", "average", "single", "none"]:
            info = {'success': False, 'info': '层次聚类算法错误!', 'code':'C2300404'}
            return json.dumps(info)
        if not data.dist in ["euclidean", "braycurtis", "manhattan", "jaccard", "canberra", "hamming", "chebyshev",
                             "minkowski", "cosine"]:
            info = {'success': False, 'info': '距离算法错误!', 'code':'C2300405'}
            return json.dumps(info)
        #metab_table_path = metab_table_path.replace("/mnt/ilustre/tsanger-data/","s3://")
        print metab_table_path
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("specimen", select_samples),
            ("tree", ""),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("table_type",table_type)
        ]
        # prepare option for workflow
        if data.dist == "manhattan":
            dist = "cityblock"
        else:
            dist = data.dist
        options = {
            "metab_table": metab_table_path,
            "samples": ",".join(select_samples),
            "dist": dist,
            "coefficient": data.coefficient,
            "cluster": data.cluster,
            "transform":data.transform,  #zgq 20190605
            "name": name
        }
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("exp_corr", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "exp_corr"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
            m_table_name = '2.SampleComp/01.ExpCorr/' + name
        else:
            m_table_name = name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn,
                            task_id=task_id, params=params_json)
        task_info = super(ExpCorrAction, self).POST()
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
