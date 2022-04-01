# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhou'
# last_modifiy = modified 2021.11.30

import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from biocluster.config import Config
import os, re
from biocluster.file import download


class RelationKeggClusterAction(MetabolomeController):
    def __init__(self):
        super(RelationKeggClusterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'metab_keggp_table', 'metab_set_table',
                        'trans_keggp_main_id', 'trans_geneset_main_id', 'padjust_method', "select", "pathway_method",
                        'set_method']
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        metabolome = Metabolome()
        task_name = 'metabolome.report.relation_kegg_cluster'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, member_id = metabolome.get_project_info(task_id)
        name = "RelationKeggCluster_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if hasattr(data, "pathway_method") and data.pathway_method not in ["hierarchy", "kmeans", "none"]:
            info = {"success": False, "info": "通路聚类算法只能为hierarchy,kmeans,无！", 'code': 'C2300309'}
            return json.dumps(info)
        if hasattr(data, "set_method") and data.set_method not in ["hierarchy", "kmeans", "none"]:
            info = {"success": False, "info": "集合聚类算法只能为hierarchy,kmeans,无！", 'code': 'C2300309'}
            return json.dumps(info)
        dist_list = ["euclidean", "braycurtis", "manhattan"]
        if hasattr(data, "pathway_dist") and data.pathway_dist not in dist_list:
            info = {'success': False, 'info': '通路距离算法选择错误!', 'code': 'C2300304'}
            return json.dumps(info)
        if hasattr(data, "set_dist") and data.set_dist not in dist_list:
            info = {'success': False, 'info': '集合距离算法选择错误!', 'code': 'C2300304'}
            return json.dumps(info)
        params_json = {
            'metab_keggp_table': data.metab_keggp_table,
            'metab_set_table': data.metab_set_table,
            'trans_keggp_main_id': data.trans_keggp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            'padjust_method': data.padjust_method,
            "select": data.select,
            "pathway_method": data.pathway_method,
            'set_method': data.set_method,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'task_id': task_id
        }
        if data.pathway_method == "hierarchy":
            params_json["pathway_dist"] = data.pathway_dist
            params_json["pathway_ctype"] = data.pathway_ctype
            params_json["pathway_n_cluster"] = data.pathway_n_cluster
        if data.pathway_method == "kmeans":
            params_json["pathway_dist"] = data.pathway_dist
            params_json["pathway_n_cluster"] = data.pathway_n_cluster
        if data.set_method == "hierarchy":
            params_json["set_dist"] = data.set_dist
            params_json["set_ctype"] = data.set_ctype
        if data.set_method == "kmeans":
            params_json["set_dist"] = data.set_dist
        mongo_data = [
            ('project_sn', project_sn),
            ('task_id', task_id),
            ('status', 'start'),
            ("name", name),
            ("desc", "Running"),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = metabolome.insert_main_table('relation_kegg_heatmap', mongo_data)
        metabolome.insert_main_id('relation_kegg_cluster', main_table_id)
        options = {
            'metab_set_table': data.metab_set_table,
            'anno_overview': data.task_id,
            'ko_overview': data.task_id,
            'trans_keggp_main_id': data.trans_keggp_main_id,
            'trans_geneset_main_id': data.trans_geneset_main_id,
            'padjust_method': data.padjust_method,
            "select": data.select,
            "pathway_method": data.pathway_method,
            'set_method': data.set_method,
            'main_table_id': str(main_table_id),
            "name": name
        }
        species = metabolome.get_kegg_species(data.task_id)
        if species:
            options["species"] = species
        else:
            info = {'success': False, 'info': 'not result'}
            return json.dumps(info)
        update_info = {str(main_table_id): 'relation_kegg_heatmap'}
        options["update_info"] = json.dumps(update_info)
        if hasattr(data, "pathway_dist") and data.pathway_dist == "manhattan":
            pathway_dist = "cityblock"
        elif hasattr(data, "pathway_dist"):
            pathway_dist = data.pathway_dist
        if hasattr(data, "set_dist") and data.set_dist == "manhattan":
            set_dist = "cityblock"
        elif hasattr(data, "set_dist"):
            set_dist = data.set_dist
        if data.pathway_method == "hierarchy":
            if not int(data.pathway_n_cluster) > 2:
                info = {"success": False, "info": "通路子聚类数目必须大于等于2！", 'code': 'C2300314'}
                return json.dumps(info)
            options["pathway_dist"] = pathway_dist
            options["pathway_ctype"] = data.pathway_ctype
            options["pathway_n_cluster"] = int(data.pathway_n_cluster)
        if data.pathway_method == "kmeans":
            if not int(data.pathway_n_cluster) > 2:
                info = {"success": False, "info": "通路子聚类数目必须大于等于2！", 'code': 'C2300314'}
                return json.dumps(info)
            options["pathway_dist"] = pathway_dist
            options["pathway_n_cluster"] = int(data.pathway_n_cluster)
        if data.set_method == "hierarchy":
            options["set_dist"] = set_dist
            options["set_ctype"] = data.set_ctype.lower()
        if data.set_method == "kmeans":
            options["set_dist"] = set_dist
        to_file = []
        to_file.append('metabolome.export_mul_metab_set(metab_set_table)')
        to_file.append('metabolome.export_overview(anno_overview)')
        to_file.append('metabolome.export_overview_ko(ko_overview)')
        m_table_name = "Relation/" + name.strip().split("_")[0] + '/' + name
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(RelationKeggClusterAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)
