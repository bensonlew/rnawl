# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
# from mbio.api.to_file.lnc_rna import *


class GenesetClassAction(LncRnaController):
    def __init__(self):
        super(GenesetClassAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()

        # 参数检测
        default_argu = ('geneset_id', 'geneset_type', 'submit_location', 'anno_type', 'task_type')
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "lack argument: %s" % arg, 'code': 'C2901004', 'variables': [arg]}
                return json.dumps(info)

        # 判断传入的基因集id是否存在
        geneset_info = {}
        geneset_name = list()
        for gd in data.geneset_id.split(","):
            geneset_info = self.lnc_rna.get_main_info(gd, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset not found", 'code': 'C2901005', 'variables': ''}
                return json.dumps(info)
            geneset_name.append(geneset_info['name'])

        # 插入主表部分
        if data.anno_type == "go":
            table_name = "Go"
            collection_name = "sg_geneset_go_class"
            to_file = 'lnc_rna.export_go_class(geneset_go)'
            option = {"geneset_go": data.geneset_id}
        elif data.anno_type == "cog":
            table_name = "Cog"
            collection_name = "sg_geneset_cog_class"
            to_file = 'lnc_rna.export_cog_class(geneset_cog)'
            option = {"geneset_cog": data.geneset_id}
        elif data.anno_type == "kegg":
            table_name = "Kegg"
            collection_name = "sg_geneset_kegg_class"
            to_file = ['lnc_rna.export_multi_gene_list(geneset_kegg)',
                       "lnc_rna.export_kegg_table(kegg_table)'"]
            option = {"geneset_kegg": data.geneset_id, "kegg_table": data.geneset_id.split(",")[0]}
        else:
            info = {'success': False, 'info': '不支持的功能分类!', 'code': 'C2901006', 'variables': ''}
            return json.dumps(info)

        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "anno_type": data.anno_type,
            "geneset_type": data.geneset_type,
            "task_id": data.task_id,
            # "type": data.type,
        }

        #   add geneset names
        if len(geneset_name) > 1:
            main_table_name = "|".join(geneset_name) + "_" + table_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        else:
            main_table_name = geneset_name[0] + "_" + table_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        task_name = 'lnc_rna.report.geneset_class'
        task_info = self.lnc_rna.get_task_info(geneset_info['task_id'])
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.lnc_rna.insert_main_table(collection_name, mongo_data)
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
            "type": "origin", # 因lncrna分析中没有区分所以写成固定值
        }
        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetClassAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        geneset_info = self.lnc_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))

        # if 'group_id' in data and str(data.group_id).lower() != 'all':
        #     _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        # if 'control_id' in data:
        #     _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_cog(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/lnc_rna/geneset_class "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                submit_location="genesetcog",
                # db_type="lnc_rna",
                geneset_id="5c90a82017b2bf147206fe98,5c90894a17b2bf407167bc7f",
                geneset_type="T",
                anno_type="cog",
                # geneset_go="5b07cc2517b2bf2d03d010e1",
                task_id="lnc_rna",
                task_type="2",
                # type="origin"
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)

        def test_go(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/lnc_rna/geneset_class "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                submit_location="genesetgo",
                # db_type="lnc_rna",
                geneset_id="5c90a82017b2bf147206fe98,5c90894a17b2bf407167bc7f",
                geneset_type="T",
                anno_type="go",
                # geneset_go="5b07cc2517b2bf2d03d010e1",
                task_id="lnc_rna",
                task_type="2",
                # type="origin"
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()
