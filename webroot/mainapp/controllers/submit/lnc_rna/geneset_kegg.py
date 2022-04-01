# -*- coding: utf-8 -*-
# 2019-03-21
import web
import json
import datetime
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import os
import unittest


class GenesetKeggAction(LncRnaController):
    def __init__(self):
        super(GenesetKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()

        # 参数检测
        default_argv = ('geneset_id', 'geneset_type', 'submit_location', 'task_type', 'task_id')
        for argv in default_argv:
            if not hasattr(data, argv):
                info = {'success': False, 'info': '%s参数缺少!' % argv, 'code': 'C2901601', 'variables': [argv]}
                return json.dumps(info)

        # 判断传入的基因集id是否存在
        geneset_info = {}
        geneset_documents = []
        geneset_name = list()
        for geneset in data.geneset_id.split(","):
            geneset_info = self.lnc_rna.get_main_info(geneset, 'sg_geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!", 'code': 'C2901602', 'variables': ''}
                return json.dumps(info)
            geneset_documents.append(geneset_info)
            geneset_name.append(geneset_info['name'])
        source = 'non_diff_exp'
        if len(geneset_documents) == 1:
            my_result = geneset_documents[0]
            if 'source' in my_result and my_result['source'] == 'diff_exp':
                source = 'diff_exp'
        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.lnc_rna.get_task_info(geneset_info['task_id'])

        collection_name = "sg_geneset_kegg_class"
        to_file = [# 'lnc_rna.export_multi_gene_list(geneset_kegg)',
                   'lnc_geneset.export_multi_gene_list(geneset_kegg)',
                   "lnc_rna.export_kegg_table(kegg_table)",
                   "lnc_rna.export_kegg_level_table(kegg_table_2)",
                   "lnc_rna.export_add_info(add_info)"]
        if len(geneset_name) > 1:
            main_table_name = '|'.join(geneset_name) + '_Kegg_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        else:
            main_table_name = geneset_name[0] + '_Kegg_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "geneset_type": data.geneset_type,
            "task_id": data.task_id,
            "anno_type": data.anno_type
            # "type": 'origin'
        }
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
        ]
        main_table_id = self.lnc_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}
        # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_kegg": json.dumps({'geneset_id': data.geneset_id, 'source': source}),
            "kegg_table": data.geneset_id.split(",")[0],
            "kegg_table_2": data.geneset_id.split(",")[0],
            "geneset_type": data.geneset_type,
            "type": "origin",
            "add_info": geneset_info['task_id'] + "\t" + data.geneset_type,
            "source": source
        }

        if "database_version" in  task_info:
            kegg_version = task_info["database_version"].get("kegg", "2017")
            options.update({"kegg_version": kegg_version})
        self.set_sheet_data(
            name='lnc_rna.report.geneset_kegg',
            module_type="workflow",
            project_sn=task_info['project_sn'],
            task_id=task_info['task_id'],
            main_table_name=main_table_name,
            to_file=to_file,
            options=options)

        task_info = super(GenesetKeggAction, self).POST()

        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        geneset_info = self.lnc_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))

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
            cmd += "s/lnc_rna/geneset_kegg "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                submit_location="genesetkegg",
                # db_type="lnc_rna",
                # geneset_id="5c90a82017b2bf147206fe98,5c90894a17b2bf407167bc7f",
                geneset_id="5c90894a17b2bf407167bc7f",
                geneset_type="T",
                anno_type="kegg",
                # geneset_kegg="5b0bfa56a4e1af2d43a6c12a",
                task_id="lnc_rna",
                task_type="2",
                # type="origin"
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()
