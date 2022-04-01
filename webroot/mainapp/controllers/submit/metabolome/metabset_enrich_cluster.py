# -*- coding: utf-8 -*-
import os
import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
import sqlite3
import unittest

class MetabsetEnrichClusterAction(MetabolomeController):
    def __init__(self):
        super(MetabsetEnrichClusterAction, self).__init__(instant=False)
    @check_sig
    def POST(self):
        data = web.input()
        print data
        for arg in ['metabset', 'correct', 'task_id', 'submit_location', 'task_type', 'submit_location','bg','correct','pathway_cluster_method',
                    'pathway_dist','pathway_ctype', 'metab_cluster_method', 'metabset_dist', 'metabset_ctype', 'n_cluster', 'select','scale']:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        task_name = 'metabolome.report.enrich_cluster'
        module_type = 'workflow'
        params_json = {
            'metabset': data.metabset,
            'correct': data.correct,
            'task_type': int(data.task_type),
            'submit_location': data.submit_location,
            'bg' : data.bg,
            'pathway_cluster_method' : data.pathway_cluster_method,
            'pathway_dist' : data.pathway_dist,
            'pathway_ctype' : data.pathway_ctype,
            'metab_cluster_method' : data.metab_cluster_method,
            'metabset_dist' : data.metabset_dist,
            'metabset_ctype' : data.metabset_ctype,
            'n_cluster' : data.n_cluster,
            'select' : data.select,
            'scale' : data.scale,
            'task_id' : data.task_id
        }
        task_info = self.Metabolome.get_task_info(data.task_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("metab_set",data.metabset)
        ]
        metabolome = Metabolome()

        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            main_table_name = "MetabsetEnrichHeatmap_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            m_table_name = "5.Metabset/11.MetabsetEnrichHeatmap/" + main_table_name
        else:
            main_table_name = "MetabsetEnrichCluster_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            m_table_name = "Metabset/" + main_table_name.strip().split("_")[0] + '/' + main_table_name
        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        main_table_id = metabolome.insert_main_table('metabset_kegg_heatmap', mongo_data)
        metabolome.insert_main_id('metabset_kegg_heatmap', main_table_id)
        update_info = {str(main_table_id): 'metabset_kegg_heatmap'}

        options = {
            'metabset': data.metabset,
            'correct': data.correct,
            'anno_overview': data.task_id,
            'ko_overview': data.task_id,  # add by ghd @20191015
            'main_table_id': str(main_table_id),
            'update_info': json.dumps(update_info),
            'bg': data.bg,
            'pathway_method' : data.pathway_cluster_method,
            'pathway_dist' : data.pathway_dist,
            'pathway_ctype' : data.pathway_ctype,
            'metabset_method' : data.metab_cluster_method,
            'metabset_dist' : data.metabset_dist,
            'metabset_ctype' : data.metabset_ctype,
            'n_cluster' : data.n_cluster,
            'select' : data.select,
            'scale' : data.scale,
            "name": main_table_name
        }
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
        if data.bg == "species":
            species = metabolome.get_kegg_species(data.task_id)
            if species:
                options["species"] = species
            else:
                info = {'success': False, 'info': 'not result'}
                return json.dumps(info)
        to_file = ['metabolome.export_mul_metab_set(metabset)', 'metabolome.export_overview(anno_overview)', 'metabolome.export_overview_ko(ko_overview)']  # add ko_overview by ghd @20191015
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        post_info = super(MetabsetEnrichClusterAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        self.Metabolome.insert_set_info(data.metabset, "metabset_kegg_heatmap", main_table_id)
        return json.dumps(post_info)
