# -*- coding: utf-8 -*-
import os
import web
import json
import datetime
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig


class RelationKeggpViewAction(MetabolomeController):
    def __init__(self):
        super(RelationKeggpViewAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        args = ['metab_expdiff', 'metab_set', 'gene_expdiff','gene_set', 'metab_annot', 'gene_annot', 'task_id', 'metab_group', 'gene_group', 'submit_location', 'task_type']
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': 'Argument:{} required'.format(arg)}
                return json.dumps(info)
        task_name = 'metabolome.report.relation_keggp_view'
        module_type = 'workflow'
        data.task_type = int(data.task_type)
        params_json = {k: v for k, v in data.items() if k in args}
        task_info = self.Metabolome.get_task_info(data.task_id)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', data.task_id),
            ('status', 'start'),
            ('desc', '正在计算'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ]
        main_table_name = "RelationKeggpView_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        m_table_name = "Relation/" + main_table_name.strip().split("_")[0] + '/' + main_table_name

        mongo_data.append(('name', main_table_name))
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(',', ':'))))
        rere_table_path = self._get_file_path('/pathway_img', data.task_id, m_table_name)
        mongo_data.append(('result_dir', rere_table_path))
        metabolome = Metabolome()
        main_table_id = metabolome.insert_main_table('relation_keggpview', mongo_data)
        relation_info = metabolome.get_mongo_common("sg_relation_analysis", {"task_id": data.task_id, "delete_by": ''})
        anno_overview = metabolome.get_mongo_common("anno_overview", {"task_id": data.task_id})
        metabolome.insert_main_id('relation_keggpview', main_table_id)
        update_info = {str(main_table_id): 'relation_keggpview'}
        options = {
            'update_info': json.dumps(update_info),
            'main_table_id': str(main_table_id),
            'metab_expdiff': data.metab_expdiff,
            'metab_set': data.metab_set,
            'gene_expdiff': data.gene_expdiff,
            'gene_set': data.gene_set,
            'metab_annot': str(anno_overview['_id']),
            'gene_annot': data.gene_annot,
            'relation_proj_type': relation_info["relate_project_type"],
            'relation_task_id': relation_info["relate_task_id"],
            'gene_group': data.gene_group,
            'metab_group': data.metab_group
        }
        to_file = None
        self.set_sheet_data(name=task_name, options=options,
                            main_table_name=m_table_name,
                            task_id=data.task_id,
                            project_sn=task_info['project_sn'],
                            module_type=module_type, params=params_json, to_file=to_file)
        post_info = super(RelationKeggpViewAction, self).POST()
        if post_info['success']:
            post_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        self.Metabolome.insert_set_info(data.metab_set, "relation_keggpview", main_table_id)
        return json.dumps(post_info)
