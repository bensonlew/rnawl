# -*- coding: utf-8 -*-
from biocluster.config import Config
import datetime
import unittest
import json
import web
import os
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.ref_rna_v2 import *
from mbio.api.to_file.ref_rna_v2 import *
from bson.objectid import ObjectId

class ToolCircularbarAction(RefRnaV2Controller):
    def __init__(self):
        super(ToolCircularbarAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['anno_id', 'geneset', 'web_page', 'tool_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        if hasattr(data, 'submit_type'):
            submit_type = int(data.submit_type)
        else:
            submit_type = 0
        sg_task = self.ref_rna_v2.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        task_id = data.task_id
        # create main table record
        # anno_id_list = str(data.anno_id).split(',')
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            anno_id=data.anno_id,
            geneset=data.geneset,
            web_page=data.web_page,
            tool_type=data.tool_type,

        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        if data.web_page == "geneset_go":
            name = "AnnotationGo_Circularbar" + '_'
        elif data.web_page == "geneset_kegg":
            name = "AnnotationKegg_Circularbar" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='AnnotationGo Circularbar main table',
            params=params,
            submit_type=submit_type,
            status="start"
        )
        main_id = self.ref_rna_v2.insert_main_table('sg_tool_lab_circularbar', main_info)
        # prepare option for workflow
        options = {
            "task_id": data.task_id,
            "main_id": str(main_id),
            "relate_name": name,
            "anno_id_geneset": str(data.anno_id)+"____"+str(data.geneset)+"____"+str(data.web_page),
            'params': params,
            'tool_type': data.tool_type,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_circularbar"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        to_files = ["ref_rna_tool.export_circularbar(anno_id_geneset)"]
        task_name = 'ref_rna_v3.report.tool_circularbar'
        self.set_sheet_data(name=task_name,
                            options=options,
                            to_file=to_files,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolCircularbarAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # 更新基因集的使用信息
        return json.dumps(task_info)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/lustre/users/sanger-dev/wpm2/sanger_bioinfo/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += '-dbversion {} '.format(1)
        cmd += "s/ref_rna_v3/tool_circularbar "
        cmd += "-b http://wpm2.sanger.com "
        # args = dict(
        #     task_id="g4f8_fhf43ch9vfoet93n4n8gi7",
        #     task_type="2",
        #     submit_location="AnnotationGo",
        #     anno_id="6115c6912f51e2fd6f32b686",
        #     geneset="Con_vs_CM_G",
        #     web_page="geneset_go",
        #     tool_type='circularbar'
        # )
        args = dict(
            task_id="g4f8_fhf43ch9vfoet93n4n8gi7",
            task_type="2",
            submit_location="AnnotationKEGG",
            anno_id="6115dfe02f51e2fd6f32b692",
            geneset="Con_vs_CM_G",
            web_page="geneset_kegg",
            tool_type='circularbar'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
