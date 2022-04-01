# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import datetime
import json
import os
import unittest

import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig


class GenesetCircAction(WholeTranscriptomeController):
    def __init__(self):
        super(GenesetCircAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['enrich_id', 'enrich_type', 'diff_id', 'compare_group', 'submit_location', 'task_type','task_id']
        print("input par is {}".format(data))
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2900901', 'variables': variables}
                print("return info is {}".format(json.dumps(info)) )
                return json.dumps(info)

        task_name = 'whole_transcriptome.report.geneset_circ'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "enrich_id": data.enrich_id,
            "diff_id": data.diff_id,
            "compare_group": data.compare_group,
            "task_id" : data.task_id,
            "enrich_type": data.enrich_type,
            "level": data.level
        }

        if data.enrich_type == "GO":
            if data.go_type not in ['All', 'ALL', 'BP', 'CC', 'MF']:
                info = {"success": False, "info": u"不存在%s这种GO富集结果！!", "variables":[data.go_type], "code" : "C2900907"}
                return info
            else:
                params_json.update(dict(go_type=data.go_type.upper()))

        # 弦图筛选参数
        if hasattr(data, "anno_list"):
            if data.enrich_type == 'GO':
                params_json.update(dict(anno_list=data.anno_list.upper()))
            elif data.enrich_type == 'KEGG':
                params_json.update(dict(anno_list=data.anno_list))
        elif hasattr(data, "p_thre") and hasattr(data, "anno_num_thre"):
            params_json.update(dict(p_thre=float(data.p_thre), p_type=data.p_type, anno_num_thre=int(data.anno_num_thre)))
        else:
            info = {"success": False, "info": u"目标列表和筛选阈值至少填入一个", "code" : "C2900908"}
            return info

        # 判断传入的富集id是否存在
        if data.enrich_type == "GO":
            enrich_info = self.whole_transcriptome.get_main_info(data.enrich_id, 'geneset_go_enrich', data.task_id)
        elif data.enrich_type == "KEGG":
            enrich_info = self.whole_transcriptome.get_main_info(data.enrich_id, 'geneset_kegg_enrich', data.task_id)
        else:
            info = {"success": False, "info": "不存在%s这种富集结果！!", "variables":[data.enrich_type], "code" : "C2900909"}
            return info

        if not enrich_info:
            info = {"success": False, "info": "富集结果不存在，请确认参数是否正确！!", 'code': 'C2900902', 'variables': ''}
            return json.dumps(info)

        #判断是否enrich表与diff表基因集属性是否一致，新基因相关
        if not self.whole_transcriptome.consistent_in_enrich_and_diff(data.enrich_id, data.diff_id, data.enrich_type):
            info = {"success": False, "info": "富集结果里存在新基因，但差异结果中没有", "code" : "C2900910"}
            return json.dumps(info)

        task_info = self.whole_transcriptome.get_task_info(enrich_info['task_id'])

        # 判断传入的差异分析是否存在 ，暂时没写
        # table_name = "Kegg"
        collection_name = "geneset_circ"
        to_file = []
        if data.enrich_type == "GO":
            to_file = ['whole_transcriptome.geneset.export_go_enrich_matrix(enrich_table)',
                       "whole_transcriptome.geneset.export_compare_exp_fc(diff_fc)",
                       "whole_transcriptome.geneset.get_gene_detail(gene_detail)"
            ]
        elif data.enrich_type == "KEGG":
            to_file = ['whole_transcriptome.geneset.export_kegg_enrich_matrix(enrich_table)',
                       "whole_transcriptome.geneset.export_compare_exp_fc(diff_fc)",
                       "whole_transcriptome.geneset.get_gene_detail(gene_detail)"
            ]

        # proteinset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        main_table_name = "EnrichCirc_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

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
            "enrich_id": data.enrich_id,
            "diff_id": data.diff_id,
            "compare_group": data.compare_group,
            "task_id" : data.task_id,
            "enrich_type": data.enrich_type
        }

        option = {
            "diff_fc": "diff_fc",
            "gene_detail": task_info['task_id'],
            "enrich_table": "enrich_table"
        }
        options.update(option)

        if hasattr(data, "anno_list"):
            if data.enrich_type == 'GO':
                options.update(dict(anno_list=data.anno_list.upper()))
            elif data.enrich_type == 'KEGG':
                options.update(dict(anno_list=data.anno_list))
        elif hasattr(data, "p_type") and hasattr(data, "anno_num_thre") and data.p_type == "padjust":
            options.update(dict(padj_thre=float(data.p_thre), anno_num_thre=int(data.anno_num_thre)))
        elif hasattr(data, "p_type") and hasattr(data, "anno_num_thre") and data.p_type == "pvalue":
            options.update(dict(p_thre=float(data.p_thre), anno_num_thre=int(data.anno_num_thre)))
        else:
            pass

        if data.enrich_type == "GO":
            options.update(dict(go_type=data.go_type))
        else:
            pass
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'])

        task_info = super(GenesetCircAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
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
        cmd += "s/whole_transcriptome/geneset_circ "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_36088",
            task_type="2",
            submit_location="genesetcirc",
            enrich_id="5ddb609017b2bf40bab6d8c9",
            enrich_type="GO",
            compare_group="control|Stable",
            diff_id="5dd5173a17b2bf7328fe8237",
            p_thre="1",
            p_type="padjust",
            anno_num_thre="15",
            level='T',
            go_type='ALL'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()


