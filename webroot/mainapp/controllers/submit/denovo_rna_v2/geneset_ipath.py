# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import web
import json
import datetime
import unittest
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import os


class GenesetIpathAction(DenovoRnaV2Controller):
    def __init__(self):
        super(GenesetIpathAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'submit_location', 'type', 'geneset_type', 'task_type','task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                var = []
                var.append(argu)
                info = {'success': False, 'info': "Lack argument: %s" % (argu), "code": 'C1601301', "variables": var}
                return json.dumps(info)

        task_name = 'denovo_rna_v2.report.geneset_ipath'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "task_id": data.task_id,
            "geneset_type": data.geneset_type,
            "type": data.type,
        }
        # 判断传入的基因集id是否存在
        geneset_info = {}
        for geneset in data.geneset_id.split(","):
            geneset_info = self.denovo_rna_v2.get_main_info(geneset, 'sg_geneset', data.task_id)
            if geneset_info['type'] != data.geneset_type:
                var = []
                var.append(geneset_info['type'])
                var.append(data.geneset_type)
                info = {"success": False, "info": "基因集属性为%s 但筛选基因集属性为%s"%(geneset_info['type'], data.geneset_type), "code": 'C1601302', "variables": var}
                return json.dumps(info)
            if not geneset_info:
                var = []
                var.append(geneset)
                info = {"success": False, "info": "geneset %s不存在，请确认参数是否正确！!"%(geneset), "code": 'C1601303', "variables": var}
                return json.dumps(info)
        if len(data.geneset_id.split(",")) > 2:
            info = {"success": False, "info": "最多选择两个蛋白集!", "code": 'C1601304', "variables": ''}
            return json.dumps(info)
        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.denovo_rna_v2.get_task_info(geneset_info['task_id'])

        table_name = "Kegg"
        collection_name = "sg_geneset_ipath"
        to_file = ['denovo_rna_v2.export_multi_gene_list(geneset_kegg)',
                   "denovo_rna_v2.export_kegg_table(kegg_table)"]

        # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        option = {"geneset_kegg": data.geneset_id,
                  "kegg_table": data.geneset_id.split(",")[0]}

        main_table_name = 'geneset' + table_name + "Ipath_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ("version", "v2"),
            ('geneset_type', data.geneset_type),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ("pathways", ['Metabolism.svg', 'Secondary_metabolites.svg', 'Antibiotics.svg', 'Microbial_metabolism.svg'])
        ]
        main_table_id = self.denovo_rna_v2.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type,
            "t2g": task_info['assemble_t2g'],
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
        geneset_info = self.denovo_rna_v2.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/denovo_rna_v2/geneset_ipath "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="denovo_rna_v2_upgrade",
            task_type="2",
            submit_location="genesetipath",
            geneset_type="G",
            geneset_id="5d53a00217b2bf317a20efb0",
            type="origin"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
