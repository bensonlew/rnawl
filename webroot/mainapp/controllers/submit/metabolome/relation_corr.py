# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhou'
# last_modifiy = modified 2021.11.30

import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig


class RelationCorrAction(MetabolomeController):
    def __init__(self):
        super(RelationCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'metab_table', 'metab_set_table', 'group_detail',
                        'group_id', 'trans_exp_main_id', 'trans_geneset_main_id', "group_method", 'coefficient',
                        'metab_cluster_is', 'trans_cluster_is', 'top_value']
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.relation_corr'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        name = "RelationCorr_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if not int(data.group_method) in [0, 2, 3]:
            info = {'success': False, 'info': '分组计算类型选择错误', 'code': 'C2300308'}
            return json.dumps(info)
        if not data.coefficient in ["spearman", "pearson", "kendall"]:
            info = {'success': False, 'info': '相关性算法类型错误!', 'code': 'C2300308'}
            return json.dumps(info)
        dist_list = ["euclidean", "braycurtis", "manhattan", "jaccard", "canberra", "hamming", "chebyshev", "minkowski", "cosine"]
        if hasattr(data, "dist_method") and data.dist_method not in dist_list:
            info = {'success': False, 'info': '代谢物距离算法错误!', 'code': 'C2300304'}
            return json.dumps(info)
        if hasattr(data, "dist_method") and data.dist_method == "manhattan":
            dist_method = "cityblock"
        elif hasattr(data, "dist_method"):
            dist_method = data.dist_method
        if hasattr(data, "cluster_method") and data.cluster_method not in ["hierarchy", "kmeans"]:
            info = {"success": False, "info": "代谢物聚类算法只能为hierarchy,kmeans！", 'code': 'C2300309'}
            return json.dumps(info)
        params_json = {
            'metab_table': data.metab_table,
            'metab_set_table': data.metab_set_table,
            "group_detail": json.loads(data.group_detail),
            'group_id': data.group_id,
            'trans_exp_main_id': data.trans_exp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            'coefficient': data.coefficient,
            'group_method': int(data.group_method),
            "metab_cluster_is": data.metab_cluster_is,
            "trans_cluster_is": data.trans_cluster_is,
            "top_value": int(data.top_value),
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id
        }
        if data.metab_cluster_is == "yes" or data.trans_cluster_is == "yes":
            params_json["cluster_method"] = data.cluster_method
            params_json["subcluster_num"] = int(data.subcluster_num)
            if data.cluster_method == "hierarchy":
                params_json["dist_method"] = data.dist_method
                params_json["cluster_way"] = data.cluster_way
            elif data.cluster_method == "kmeans":
                params_json["dist_method"] = data.dist_method
        else:
            pass

        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('status', 'start'),
            ("name", name),
            ("desc", "Running"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = metabolome.insert_main_table('relation_corr', mongo_data)
        metabolome.insert_main_id('relation_corr', main_table_id)

        if project_type == "LC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
        options = {
            'metab_table': metab_table_path,
            "metab_desc": metab_desc_path,
            'metab_set_table': data.metab_set_table,
            'metab_anno': data.task_id,
            'trans_exp_main_id': data.trans_exp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            "group_detail": data.group_detail,
            'group_method': int(data.group_method),
            'coefficient': data.coefficient,
            "metab_top": int(data.top_value),
            "trans_top": int(data.top_value),
            'main_table_id': str(main_table_id),
            "name": name
        }
        if table_name == "MetabTable_Origin":
            scale = metabolome.get_scale_type(data.metab_table, task_id)
            if scale and scale == "none":
                options["log10"] = True
            else:
                options["log10"] = False
        elif table_name == "raw":
            options["log10"] = True
        to_file = []
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        to_file.append('metabolome.export_group_by_detail(group_detail)')
        to_file.append('metabolome.export_overview(metab_anno)')
        info2 = metabolome.get_hmdb_anno(data.task_id)
        if info2 == "":
            pass
        else:
            options["metab_hmdb_anno"] = data.task_id
            to_file.append('metabolome.export_hmdb_level(metab_hmdb_anno)')
        if data.metab_cluster_is == "yes":
            if not hasattr(data, "subcluster_num"):
                info = {"success": False, "info": "选择聚类时必须输入子聚类数目！", 'code': 'C2300313'}
                return json.dumps(info)
            if not int(data.subcluster_num) > 1:
                info = {"success": False, "info": "代谢物子聚类数目必须大于等于2！", 'code': 'C2300314'}
                return json.dumps(info)
            options["metab_cluster_method"] = data.cluster_method
            options["metab_n_cluster"] = int(data.subcluster_num)
            if data.cluster_method == "hierarchy":
                if not hasattr(data, "dist_method") or not hasattr(data, "cluster_way"):
                    info = {"success": False, "info": "代谢物hierarchy聚类方法时必须输入距离算法和层级聚类方式！", 'code': 'C2300311'}
                    return json.dumps(info)
                options["metab_dist"] = dist_method
                options["metab_cluster"] = data.cluster_way
            if data.cluster_method == "kmeans":
                if not hasattr(data, "dist_method"):
                    info = {"success": False, "info": "kmeans聚类方法时必须输入距离算法！", 'code': 'C2300313'}
                    return json.dumps(info)
                options["metab_dist"] = dist_method
        else:
            options["metab_cluster_method"] = "no"
        if data.trans_cluster_is == "yes":
            if not hasattr(data, "subcluster_num"):
                info = {"success": False, "info": "选择聚类时必须输入子聚类数目！", 'code': 'C2300313'}
                return json.dumps(info)
            if not int(data.subcluster_num) > 1:
                info = {"success": False, "info": "代谢物子聚类数目必须大于等于2！", 'code': 'C2300314'}
                return json.dumps(info)
            options["trans_cluster_method"] = data.cluster_method
            options["trans_n_cluster"] = int(data.subcluster_num)
            if data.cluster_method == "hierarchy":
                if not hasattr(data, "dist_method") or not hasattr(data, "cluster_way"):
                    info = {"success": False, "info": "代谢物hierarchy聚类方法时必须输入距离算法和层级聚类方式！", 'code': 'C2300311'}
                    return json.dumps(info)
                options["trans_dist"] = dist_method
                options["trans_cluster"] = data.cluster_way
            if data.cluster_method == "kmeans":
                if not hasattr(data, "dist_method"):
                    info = {"success": False, "info": "kmeans聚类方法时必须输入距离算法！", 'code': 'C2300313'}
                    return json.dumps(info)
                options["trans_dist"] = dist_method
        else:
            options["trans_cluster_method"] = "no"
        update_info = {str(main_table_id): 'relation_corr'}
        options["update_info"] = json.dumps(update_info)
        task_info = metabolome.get_task_info(data.task_id)
        m_table_name = "Relation/" + name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(RelationCorrAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)