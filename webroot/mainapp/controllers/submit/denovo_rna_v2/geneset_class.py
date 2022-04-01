# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os

class GenesetClassAction(DenovoRnaV2Controller):
    def __init__(self):
        super(GenesetClassAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        default_argu = ['geneset_id', 'geneset_type', 'submit_location', 'anno_type', 'type', 'task_type']
        for arg in default_argu:
            if not hasattr(data, arg):
                var = []
                var.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % (arg), "code": 'C1600901', "variables": var}
                return json.dumps(info)

        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "anno_type": data.anno_type,
            "geneset_type": data.geneset_type,
            "task_id": data.task_id,
            "type": data.type,
        }
        # 判断传入的基因集id是否存在
        #print(data.geneset_id)
        geneset_info = {}
        for gd in data.geneset_id.split(","):
            geneset_info = self.denovo_rna_v2.get_main_info(gd, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset not found", "code": 'C1600902', "variables": ''}
                return json.dumps(info)
        task_info = self.denovo_rna_v2.get_task_info(geneset_info['task_id'])
        if "version" in task_info:
           version = task_info["version"]
        else:
             version = "v1"
        # version="v2"


        if data.anno_type == "go":
            table_name = "Go"
            collection_name = "sg_geneset_go_class"
            to_file = 'denovo_rna_v2.export_go_class(geneset_go)'
            option = {"geneset_go": data.geneset_id}
            task_name = 'denovo_rna_v2.report.geneset_class'
        elif data.anno_type == "cog":
            table_name = "Cog"
            collection_name = "sg_geneset_cog_class"
            if version == "v2":
                to_file = 'denovo_rna_v3.export_cog_class(geneset_cog)'
                task_name = 'denovo_rna_v3.report.geneset_class'
            else:
                to_file = 'denovo_rna_v2.export_cog_class(geneset_cog)'
                task_name = 'denovo_rna_v2.report.geneset_class'
            option = {"geneset_cog": data.geneset_id}
        elif data.anno_type == "kegg":
            table_name = "Kegg"
            collection_name = "sg_geneset_kegg_class"
            if version == "v2":
                task_name = 'denovo_rna_v3.report.geneset_class'
            else:
                task_name = 'denovo_rna_v2.report.geneset_class'
            to_file = ['denovo_rna_v2.export_multi_gene_list(geneset_kegg)', "denovo_rna_v2.export_kegg_table(kegg_table)'"]
            option = {"geneset_kegg": data.geneset_id, "kegg_table": data.geneset_id.split(",")[0]}
        else:
            info = {'success': False, 'info': '不支持的功能分类!', "code": 'C1600903', "variables": ''}
            return json.dumps(info)

        main_table_name = 'Geneset' + table_name + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        #print(main_table_name)
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        if version == "v2":
            mongo_data.append(('version', 'v2'))
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
            "anno_type": data.anno_type,
            "type": data.type,
            }
        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=main_table_name,
                            module_type="workflow",
                            to_file=to_file,
                            project_sn=task_info['project_sn'],
                            task_id=task_info['task_id'])

        task_info = super(GenesetClassAction, self).POST()

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


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/denovo_rna_v2/geneset_class "
            cmd += "-b http://bcl.tsg.com "
            args = dict(
                submit_location="geneset_class",
                #db_type="ref_rna_v2",
                #geneset_id="5d7992ee17b2bf1d4000b32a",
                geneset_id="5d76432e17b2bf0ffc38262f",
                #geneset_type="G",
                #anno_type="kegg",
                # kegg_table="5b1924afa4e1af33064178b3",
                anno_type="go",
                type="origin",
                geneset_type='T',
                task_id="tsg_35517",
                task_type="2",
                # group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
                # draw_in_groups="yes"
                #type="origin",
                #method="bh"
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()