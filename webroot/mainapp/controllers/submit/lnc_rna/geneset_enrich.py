# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.core.basic import Basic
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import os
import unittest


class GenesetEnrichAction(LncRnaController):
    def __init__(self):
        super(GenesetEnrichAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        web.input()
        # 参数检测
        default_argu = ['geneset_id', 'geneset_type', 'submit_location', 'task_type', 'anno_type', 'method', "type"]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2901301', 'variables': [argu]}
                return json.dumps(info)

        # 判断传入的基因集id是否存在
        geneset_info = self.lnc_rna.get_main_info(data.geneset_id, 'sg_geneset', data.task_id)
        if not geneset_info:
            info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!", 'code': 'C2901302', 'variables': ''}
            return json.dumps(info)
        task_info = self.lnc_rna.get_task_info(geneset_info['task_id'])
        geneset_name = geneset_info['name']

        if data.anno_type == "go":
            table_name = "Go"
            collection_name = "sg_geneset_go_enrich"
            to_file = ['lnc_rna.export_gene_list(geneset_list)',
                       'lnc_rna.export_all_list(all_list)',
                       'lnc_rna.export_go_list(go_list)']
            infile = {"go_list": data.geneset_id}
        elif data.anno_type == "kegg":
            table_name = "Kegg"
            collection_name = "sg_geneset_kegg_enrich"
            to_file = ['lnc_rna.export_gene_list(geneset_list)',
                       'lnc_rna.export_multi_gene_list(geneset_kegg)',
                       'lnc_rna.export_all_list(all_list)',
                       'lnc_rna.export_kegg_table(kegg_table)',
                       "lnc_rna.export_add_info(add_info)"]
            infile = {"kegg_table": data.geneset_id, "add_info": geneset_info['task_id'] + "\t" + data.geneset_type}
        else:
            info = {'success': False, 'info': '不支持该富集分子!', 'code': 'C2901303', 'variables': ''}
            return json.dumps(info)

        main_table_name = geneset_name + "_" + table_name + "_Enrich_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "anno_type": data.anno_type,
            "method": data.method,
            "geneset_type": data.geneset_type,
            "type": data.type,
            "task_id": data.task_id,
        }
        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ("geneset_type", data.geneset_type),
            ("geneset_name", geneset_name),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.lnc_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "geneset_kegg": data.geneset_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_type": data.geneset_type,
            "anno_type": data.anno_type,
            "method": data.method,
            "geneset_list": data.geneset_id,
            "type": data.type,
            "all_list": data.geneset_id,
            "geneset_id": data.geneset_id,
        }
        options.update(infile)

        if "database_version" in task_info and data.anno_type == "kegg":
            kegg_version = task_info["database_version"].get("kegg", "2017")
            options.update({"kegg_version": kegg_version})
        if "database_version" in task_info and data.anno_type == "go":
            go_version = task_info["database_version"].get("go", "20200628").split("_")[0]
            options.update({"go_version": go_version})

        task_name = 'lnc_rna.report.geneset_enrich'
        task_type = 'workflow'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=main_table_name,
                            module_type=task_type,
                            to_file=to_file,
                            task_id=task_info["task_id"],
                            project_sn=task_info["project_sn"])

        task_info = super(GenesetEnrichAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        self.lnc_rna.insert_geneset_info(data.geneset_id, collection_name, str(main_table_id))
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

        def test_go(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/lnc_rna/geneset_enrich "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                submit_location="genesetgo_rich",
                geneset_id="5c90894a17b2bf407167bc7f",
                geneset_type="T",
                anno_type="go",
                task_id="lnc_rna",
                task_type="2",
                type="origin",
                method="bh"
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)

        def test_kegg(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/lnc_rna/geneset_enrich "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                submit_location="genesetkegg_rich",
                geneset_id="5c90894a17b2bf407167bc7f",
                geneset_type="T",
                anno_type="kegg",
                task_id="lnc_rna",
                task_type="2",
                type="origin",
                # kegg： 'by', 'bh', 'bonferroni', 'holm'
                method="bh"
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)

    unittest.main()

"""Sheet_Data: {'UPDATE_STATUS_API': 'ref_rna.tupdate_status', 'output': u'tsanger:rerewrweset/files/m_188/
188_59f029eed6609/tsg_29390/interaction_results/GenesetGoEnrich20180528_104127467', 'id': 'tsg_29390_75287_472188', 
'stage_id': 0, 'interaction': True, 'name': 'ref_rna.report.geneset_enrich', 'db_type': 'ref_rna', 
'client': u'client03', 'project_sn': u'188_59f029eed6609', 'IMPORT_REPORT_DATA': True, 'to_file': 
['ref_rna.export_gene_list(genset_list)', 'ref_rna.export_all_list(all_list)', 'ref_rna.export_go_list(go_list)'], 
'type': 'workflow', 'options': {'geneset_type': u'gene', 'method': u'fdr', 'anno_type': u'go', 'task_type': u'', 
'genset_list': u'5ae93968a4e1af7f1cd07e16', 'main_table_id': '5b0b6c57a4e1af239b6ad5ff', 'all_list': 
u'5ae93968a4e1af7f1cd07e16', 'submit_location': u'genesetgo_rich', 'go_list': u'5ae93968a4e1af7f1cd07e16', 
'update_info': '{"5b0b6c57a4e1af239b6ad5ff": "sg_geneset_go_enrich"}'}}"""
