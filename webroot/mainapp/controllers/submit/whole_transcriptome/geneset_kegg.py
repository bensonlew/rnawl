# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

import datetime
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
from mbio.api.to_file.whole_transcriptome.geneset import *


class GenesetKeggAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetKeggAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'level', 'submit_location', 'task_type', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2901601', 'variables': variables}
                return json.dumps(info)

        task_name = 'whole_transcriptome.report.geneset_kegg'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "level": data.level,
            "task_id": data.task_id,
        }
        # 判断传入的基因集id是否存在
        geneset_info = {}
        geneset_documents = []
        if len(data.geneset_id.split(",")) >= 3:
            info = {"success": False, "info": "只能选择一个或两个基因集", "code": "C2901604"}
            return json.dumps(info)
        for geneset in data.geneset_id.split(","):
            geneset_info = self.whole_transcriptome.get_main_info(geneset, 'geneset', data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset不存在，请确认参数是否正确！!", 'code': 'C2901602', 'variables': ''}
                return json.dumps(info)
            geneset_documents.append(geneset_info)
        source = 'non_diff_exp'
        if len(geneset_documents) == 1:
            my_result = geneset_documents[0]
            if 'source' in my_result:
                if my_result['source'] == 'DE_mR_detail' or my_result['source'] == 'DE_mR_G_detail' or my_result['source'] == 'DE_mR_T_detail':
                    # 目的是过滤来自单组分析统计表额基因集
                    if self.whole_transcriptome.check_geneset_regulate(my_result["main_id"]):
                        source = 'diff'
                    else:
                        print("没有regulate_list, 核查基因集来源")

        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.whole_transcriptome.get_task_info(geneset_info['task_id'])

        table_name = "KEGG"
        collection_name = "geneset_kegg_class"
        # ref_rna_v2_new.geneset
        to_file = ['whole_transcriptome.geneset.export_multi_gene_list(geneset_kegg)',
                   "whole_transcriptome.geneset.export_kegg_table(kegg_table)",
                   "whole_transcriptome.geneset.export_kegg_level_table(kegg_table_2)",
                   "whole_transcriptome.geneset.export_add_info(add_info)"]

        # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        option = {
            "task_id": data.task_id,
            "geneset_kegg": json.dumps({'geneset_id': data.geneset_id, 'source': source}),
            "kegg_table": data.geneset_id.split(",")[0],
            "kegg_table_2": data.geneset_id.split(",")[0],
            "add_info": geneset_info['task_id'] + "\t" + data.level}

        main_table_name = table_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ('version', 'v1')
        ]
        main_table_id = self.whole_transcriptome.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
            "geneset_type": data.level,
            # "type": data.type,
            "source": source
        }

        if "database_version" in task_info["long_task"]:
            kegg_version = task_info["long_task"]["database_version"].get("kegg", "2018")
            options.update({"kegg_version": kegg_version})

        options.update(option)
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetKeggAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }}
        geneset_info = self.whole_transcriptome.insert_geneset_info(data.geneset_id, collection_name,
                                                                    str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
            cmd += "s/whole_transcriptome/geneset_kegg "
            cmd += "-b http://bcl.tsg.com "
            args = dict(
                submit_location="genesetkegg",
                db_type="whole_transcriptome",
                geneset_id="5db032b317b2bf3cdfb10638",
                level="T",
                anno_type="kegg",
                # geneset_kegg="5b0bfa56a4e1af2d43a6c12a",
                task_id="whole_transcriptome",
                task_type="2",
                # type="origin"

            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)


    unittest.main()
