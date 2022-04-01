# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import web
import json
import datetime
import unittest
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig


class GenesetIpathAction(ProkRNAController):
    def __init__(self):
        super(GenesetIpathAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'submit_location', 'type', 'task_type','task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        task_name = 'prok_rna.report.geneset_ipath'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "task_id": data.task_id,
            "type": data.type,
        }
        # 判断传入的基因集id是否存在
        geneset_info = {}
        for geneset in data.geneset_id.split(","):
            geneset_info = self.prok_rna.get_main_info(geneset, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!"}
                return json.dumps(info)
        if len(data.geneset_id.split(",")) > 2:
            info = {"success": False, "info": "最多选择两个蛋白集!"}
            return json.dumps(info)
        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.prok_rna.get_task_info(geneset_info['task_id'])

        table_name = "Kegg"
        collection_name = "sg_geneset_ipath"
        to_file = ['prok_rna.export_multi_gene_list_new(geneset_kegg)',
                   "prok_rna.export_kegg_table(kegg_table)"]

        # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        option = {"geneset_kegg": data.geneset_id,
                  "kegg_table": data.geneset_id.split(",")[0]}

        main_table_name = 'Geneset' + table_name + "Ipath_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('submit_type', submit_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ("pathways", ['Metabolic_pathways.svg','Regulatory_pathways.svg', 'Biosynthesis_of_secondary_metabolities.svg'])
        ]
        main_table_id = self.prok_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "type":data.type,
            }
        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetIpathAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        geneset_info = self.prok_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/prok_rna/geneset_ipath "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna",
            task_type="2",
            submit_location="genesetipath",
            geneset_id="5b70dcb5a4e1af118c63916c,5b70dd58a4e1af13f89b4362",
            type="origin"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
