# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import web
import json
import datetime
import unittest
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
# from mbio.api.to_file.medical_transcriptome_new.medical_transcriptome import *
from mbio.api.to_file.medical_transcriptome_new.geneset_circ import *
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
import os
import re

class GenesetCircAction(MedicalTranscriptomeController):
    def __init__(self):
        super(GenesetCircAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        default_argu = ['enrich_id', 'enrich_type', 'diff_id', 'compare_group',
                        'submit_location', 'task_type', 'task_id', 'level']
        print("input par is {}".format(data))
        for argu in default_argu:
            if not hasattr(data, argu):
                variables = []
                variables.append(argu)
                info = {'success': False, 'info': '%s参数缺少!' % argu, 'code': 'C2900901', 'variables': variables}
                print("return info is {}".format(json.dumps(info)) )
                return json.dumps(info)

        task_name = 'medical_transcriptome.report.geneset_circ'
        params_json = {
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "enrich_id": data.enrich_id,
            "diff_id": data.diff_id,
            "compare_group": data.compare_group,
            "task_id" : data.task_id,
            "enrich_type": data.enrich_type,
            "level": data.level,
        }

        if data.enrich_type == "GO":
            if data.go_type not in ['All', 'ALL', 'BP', 'CC', 'MF']:
                info = {"success": False, "info": u"不存在%s这种GO富集结果！!", "variables":[data.go_type], "code" : "C2900907"}
                return info
            else:
                params_json.update(dict(go_type=data.go_type.upper()))

        if data.level == 'T':
            if data.enrich_type.upper() not in ['GO', 'KEGG', 'REACTOME']:
                info = {"success": False, "info": u"转录本数据不支持%s富集分析！!", "variables": [data.enrich_type]}
                return info

        # 弦图筛选参数
        if hasattr(data, "anno_list"):
            if data.enrich_type.upper() in ['DO', 'REACTOME', 'DISGENET', 'GO']:
                params_json.update(dict(anno_list=data.anno_list.upper()))
            elif data.enrich_type == 'KEGG':
                params_json.update(dict(anno_list=data.anno_list))

        elif hasattr(data, "p_thre") and hasattr(data, "anno_num_thre"):
            params_json.update(dict(p_thre=float(data.p_thre), p_type=data.p_type, anno_num_thre=int(data.anno_num_thre)))
        else:
            info = {"success": False, "info": u"目标列表和筛选阈值至少填入一个", "code" : "C2900908"}
            return info

        if hasattr(data, "gene_list_path"):
            params_json.update({
                "gene_list_path": data.gene_list_path,
                "gene_list_id": data.gene_list_id
            })

        # 当选择根据指定id绘图时，对输入信息进行判定  add by fwy 20200506
        if hasattr(data, "anno_list"):
            if data.anno_list == "":
                info = {"success": False, "info": u"请输入目标GO Term/KEGG pathway/Reactome pathway/DO Term/DisGeNET IDs 进行作图"}
                return info
            else:
                if data.enrich_type == 'GO':
                    if not "go" in data.anno_list.lower():
                        info = {"success": False, "info": "Please enter the GO id"}
                        return info
                elif data.enrich_type == 'KEGG':
                    if not re.match("mmu|map|rna|hsa", data.anno_list.lower()):
                    # if not "map" in data.anno_list.lower():
                        info = {"success": False, "info": "Please enter the kegg id"}
                        return info
                elif data.enrich_type.upper() == 'DISGENET':
                    if not re.match('c', data.anno_list.lower()):
                        info = {"success": False, "info": "Please enter the DisGeNET id"}
                        return info
                elif data.enrich_type.upper() == 'DO':
                    if not "doid" in data.anno_list.lower():
                        info = {"success": False, "info": "Please enter the DO id"}
                        return info
                elif data.enrich_type.upper() == 'REACTOME':
                    if not re.match('r-', data.anno_list.lower()):
                        info = {"success": False, "info": "Please enter the Reactome id"}
                        return info

        # 判断传入的富集id是否存在
        if data.enrich_type == "GO":
            enrich_info = self.medical_transcriptome.get_main_info(data.enrich_id, 'sg_geneset_go_enrich', data.task_id)
        elif data.enrich_type.upper() == "DISGENET":
            enrich_info = self.medical_transcriptome.get_main_info(data.enrich_id, 'sg_geneset_disgenet_enrich', data.task_id)
        elif data.enrich_type.upper() == "REACTOME":
            enrich_info = self.medical_transcriptome.get_main_info(data.enrich_id, 'sg_geneset_reactome_enrich', data.task_id)
        elif data.enrich_type.upper() == "DO":
            enrich_info = self.medical_transcriptome.get_main_info(data.enrich_id, 'sg_geneset_do_enrich', data.task_id)
        elif data.enrich_type == "KEGG":
            enrich_info = self.medical_transcriptome.get_main_info(data.enrich_id, 'sg_geneset_kegg_enrich', data.task_id)
        else:
            info = {"success": False, "info": "不存在%s这种富集结果！!", "variables":[data.enrich_type], "code" : "C2900909"}
            return info

        if not enrich_info:
            info = {"success": False, "info": "富集结果不存在，请确认参数是否正确！!", 'code': 'C2900902', 'variables': ''}
            return json.dumps(info)

        #判断是否enrich表与diff表基因集属性是否一致，新基因相关
        if not self.medical_transcriptome.consistent_in_enrich_and_diff(data.enrich_id, data.diff_id, data.enrich_type.upper()):
            info = {"success": False, "info": "富集结果里存在新基因，但差异结果中没有", "code" : "C2900910"}
            return json.dumps(info)

        task_info = self.medical_transcriptome.get_task_info(enrich_info['task_id'])

        geneset_id = json.loads(enrich_info['params'])['geneset_id']
        geneset_info = self.medical_transcriptome.get_main_info(ObjectId(geneset_id), 'sg_geneset', data.task_id)
        if not geneset_info:
            info = {"success": False, "info": "基因集不存在，请确认参数是否正确！!", 'variables': ''}
            return json.dumps(info)
        geneset_name = geneset_info['name']

        # 判断传入的差异分析是否存在 ，暂时没写
        # table_name = "Kegg"
        collection_name = "sg_geneset_circ"
        to_file = []
        if data.enrich_type == "GO":
            to_file = ['medical_transcriptome_new.geneset_circ.export_go_enrich_matrix(enrich_table)',
                       "medical_transcriptome_new.geneset_circ.export_compare_exp_fc(diff_fc)",
                       "medical_transcriptome_new.geneset_circ.get_gene_detail_new(gene_detail)"
            ]
        elif data.enrich_type == "KEGG":
            to_file = ['medical_transcriptome_new.geneset_circ.export_kegg_enrich_matrix(enrich_table)',
                       "medical_transcriptome_new.geneset_circ.export_compare_exp_fc(diff_fc)",
                       "medical_transcriptome_new.geneset_circ.get_gene_detail_new(gene_detail)"
            ]
        elif data.enrich_type.upper() == "DISGENET":
            to_file = ['medical_transcriptome_new.geneset_circ.export_disgenet_enrich_matrix(enrich_table)',
                       "medical_transcriptome_new.geneset_circ.export_compare_exp_fc(diff_fc)",
                       "medical_transcriptome_new.geneset_circ.get_gene_detail_new(gene_detail)"
            ]
        elif data.enrich_type.upper() == "DO":
            to_file = ['medical_transcriptome_new.geneset_circ.export_do_enrich_matrix(enrich_table)',
                       "medical_transcriptome_new.geneset_circ.export_compare_exp_fc(diff_fc)",
                       "medical_transcriptome_new.geneset_circ.get_gene_detail_new(gene_detail)"
            ]
        elif data.enrich_type.upper() == "REACTOME":
            to_file = ['medical_transcriptome_new.geneset_circ.export_reactome_enrich_matrix(enrich_table)',
                       "medical_transcriptome_new.geneset_circ.export_compare_exp_fc(diff_fc)",
                       "medical_transcriptome_new.geneset_circ.get_gene_detail_new(gene_detail)"
            ]
        # proteinset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        main_table_name = geneset_name + "_EnrichCirc_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('version', 'v1'),
            ('status', 'start'),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':')))
        ]
        main_table_id = self.medical_transcriptome.insert_main_table(collection_name, mongo_data)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}
        update_info = {str(main_table_id): collection_name}

        options = {
            "submit_location": data.submit_location,
            "task_type": data.task_type,
            'main_table_data': main_table_data,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "enrich_id": data.enrich_id,
            "diff_id": data.diff_id,
            "compare_group": data.compare_group,
            "task_id": data.task_id,
            "enrich_type": data.enrich_type,
        }
        task_id = task_info['task_id']
        option = {
            "diff_fc": "diff_fc",
            "gene_detail": task_id+"|"+data.level,
            "enrich_table": "enrich_table"
        }
        options.update(option)

        if hasattr(data, "anno_list"):
            if data.enrich_type.upper() in ['DO', 'REACTOME', 'DISGENET', 'GO']:
                options.update(dict(anno_list=data.anno_list.upper()))
            elif data.enrich_type == 'KEGG':
                options.update(dict(anno_list=data.anno_list))
        elif hasattr(data, "p_type") and hasattr(data, "anno_num_thre") and data.p_type == "padjust":
            options.update(dict(padj_thre=float(data.p_thre), anno_num_thre=int(data.anno_num_thre)))
        elif hasattr(data, "p_type") and hasattr(data, "anno_num_thre") and data.p_type == "pvalue":
            options.update(dict(p_thre=float(data.p_thre), anno_num_thre=int(data.anno_num_thre)))
        else:
            pass

        if hasattr(data, "gene_list_path"):
            options.update({
                "gene_list": data.gene_list_path
            })

        if data.enrich_type == "GO":
            options.update(dict(go_type=data.go_type))
        else:
            pass
        # print(options)
        self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, module_type="workflow",
                            to_file=to_file, project_sn=task_info['project_sn'], task_id=task_info['task_id'],
                            new_task_id=new_task_id)

        task_info = super(GenesetCircAction, self).POST()

        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
                }}
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/medical_transcriptome/geneset_circ "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="genesetcirc",
            enrich_id="5f87a40b17b2bf0ba6415c81",
            enrich_type="do",
            compare_group="S1|S3",
            diff_id="5f45c6d117b2bf78d9c9c16d",
            p_thre="1",
            p_type="padjust",
            anno_num_thre="15",
            level='G',
            go_type='ALL',
            # gene_list_path='/mnt/ilustre/users/sanger-dev/workspace/20200508/GenesetCirc_tsg_37104_6636_7030/xx',
            # gene_list_id='aaaaa'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
