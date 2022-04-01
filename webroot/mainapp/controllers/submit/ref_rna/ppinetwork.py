# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig


class PpinetworkAction(RefRnaController):
    """
    蛋白质互作网络的接口
    """

    species_list = [30611, 9598, 61853, 9593, 9606, 9544, 9483, 30608, 9601, 9478, 10141, 10020, 10090, 9986, 10116,
                    43179, 37347, 9685, 9913, 9739, 9669, 9796, 132908, 59463, 9646, 9823, 9785, 9813, 9371, 9361,
                    28377, 9031, 13735, 9103, 59729, 8049, 31033, 8090, 8083, 69293, 99883, 8128, 7955, 13616, 9258,
                    9305, 9315, 7897, 7757, 7719, 51511, 6239, 7227, 4932, 15368, 4513, 4641, 4533, 4538, 4555,
                    4558, 4577, 59689, 3702, 3711, 3847, 3694, 4081, 4113, 29760, 88036, 3218, 3055, 45157]

    def __init__(self):
        super(PpinetworkAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'gene_type', 'species', 'submit_location', 'combine_score', 'task_type']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if str(data.species) == '':
            info = {"success": False, "info": "网页端输入的species参数为空，不合法！".format(data.species)}
            return json.dumps(info)
        if int(data.species) not in self.species_list:
            info = {"success": False, "info": "不能进行蛋白质互作分析，因为数据库中不存在该物种的蛋白质互作组数据！"}
            return json.dumps(info)
        task_name = 'ref_rna.report.ppinetwork'
        task_type = 'workflow'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            "geneset_id": data.geneset_id,
            "species": data.species,
            "combine_score": data.combine_score,
            "gene_type": data.gene_type,
            "task_id": data.task_id
        }
        geneset_info = self.ref_rna.get_main_info(data.geneset_id, 'sg_geneset')
        if not geneset_info:
            info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
            return json.dumps(info)
        task_info = self.ref_rna.get_task_info(geneset_info['task_id'])
        main_table_name = 'PPINetwork_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('geneset_id', ObjectId(data.geneset_id)),
            ('desc', 'ppi_network分析中...'),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ("gene_type", data.gene_type)
        ]
        main_table_id = self.ref_rna.insert_main_table('sg_ppinetwork', mongo_data)
        update_info = {str(main_table_id): "sg_ppinetwork"}

        options = {
            'update_info': json.dumps(update_info),
            "ppi_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "diff_exp_gene": data.geneset_id,
            "species": data.species,
            "combine_score": data.combine_score
        }
        to_file = "ref_rna.export_gene_list_ppi(diff_exp_gene)"
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type=task_type,
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(PpinetworkAction, self).POST()

        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        geneset_info = self.ref_rna.insert_geneset_info(data.geneset_id, "sg_ppinetwork", str(main_table_id))
        #if geneset_info:
            #print "geneset_info插入成功"
        return json.dumps(task_info)
