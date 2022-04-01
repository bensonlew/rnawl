# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class GenesetIpathAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetIpathAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['geneset_id', 'submit_location', 'level', 'task_type', 'task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2901501', 'variables': variables}
                return json.dumps(info)

        task_name = 'whole_transcriptome.report.geneset_ipath'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "geneset_id": data.geneset_id,
            "task_id": data.task_id,
            "level": data.level,
        }
        # 判断传入的基因集id是否存在
        geneset_info = {}
        for geneset in data.geneset_id.split(","):
            geneset_info = self.whole_transcriptome.get_main_info(geneset, 'geneset', data.task_id)
            if geneset_info['level'] != data.level:
                variables = []
                variables.append(geneset_info['level'])
                variables.append(data.level)
                info = {"success": False, "info": "基因集属性为%s 但筛选基因集属性为%s" % (geneset_info['level'], data.level),
                        'code': 'C2901502', 'variables': variables}
                return json.dumps(info)
            if not geneset_info:
                variables = []
                variables.append(geneset)
                info = {"success": False, "info": "geneset %s不存在，请确认参数是否正确！!" % geneset, 'code': 'C2901503',
                        'variables': variables}
                return json.dumps(info)
        if len(data.geneset_id.split(",")) > 2:
            info = {"success": False, "info": "最多选择两个蛋白集!", 'code': 'C2901504', 'variables': ''}
            return json.dumps(info)
        # 两个geneset_id到task返回的project_sn等信息是一样的
        task_info = self.whole_transcriptome.get_task_info(geneset_info['task_id'])

        table_name = "Kegg"
        collection_name = "geneset_ipath"
        to_file = ['whole_transcriptome.advance.export_multi_gene_list(geneset_kegg)',
                   "whole_transcriptome.advance.export_kegg_table(kegg_table)"]

        # level只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        option = {"geneset_kegg": data.geneset_id,
                  "kegg_table": data.geneset_id.split(",")[0]}

        # main_table_name = table_name + "Ipath_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        main_table_name = "iPath_Analysis_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('name', main_table_name),
            ('version', "3.0"),
            ('level', data.level),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
            ("pathways", ['Metabolism.svg', 'Secondary_metabolites.svg', 'Antibiotics.svg', 'Microbial_metabolism.svg'])
            # ("pathways", ['Metabolic_pathways.svg','Regulatory_pathways.svg', 'Biosynthesis_of_secondary_metabolities.svg'])
        ]
        main_table_id = self.whole_transcriptome.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}

        options = {
            "task_id": data.task_id,
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "geneset_id": data.geneset_id,
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
        geneset_info = self.whole_transcriptome.insert_geneset_info(data.geneset_id, collection_name,
                                                                    str(main_table_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/whole_transcriptome/geneset_ipath "
        cmd += "-b http://bcl.tsg.com  "
        args = dict(
            task_id='tsg_36088',
            task_type="2",
            submit_location="genesetipath",
            level="G",
            geneset_id="5ddc80be17b2bf305f486f62,5ddc812017b2bf3060486f62"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
