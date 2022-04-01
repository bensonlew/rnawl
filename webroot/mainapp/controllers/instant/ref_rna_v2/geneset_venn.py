# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class GenesetVennAction(RefRnaV2Controller):
    def __init__(self):
        super(GenesetVennAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['geneset_ids', 'exp_level']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg, 'code': 'C2901901', 'variables': variables}
                return json.dumps(info)
        sg_task = self.ref_rna_v2.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        task_id = data.task_id
        # create main table record
        geneset_id_list = str(data.geneset_ids).split(',')
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_level=data.exp_level,
            geneset_ids=geneset_id_list,
            geneset_id=','.join(geneset_id_list),
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "GenesetVenn" + '_' + data.exp_level + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v3",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Geneset venn analysis main table',
            params=params,
            status="end"
        )
        main_id = self.ref_rna_v2.insert_main_table('sg_geneset_venn', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "task_id": data.task_id,
            "main_id": str(main_id),
            "geneset_id": data.geneset_ids,
            "submit_location": data.submit_location,
            "main_table_data": main_table_data,
            "update_info": json.dumps({str(main_id): "sg_geneset_venn"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'ref_rna_v2.report.geneset_venn'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetVennAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                },
            'info': "成功存储参数，可进行venn分析"
        }
        # 更新基因集的使用信息
        self.ref_rna_v2.insert_geneset_info(data.geneset_ids, "sg_geneset_venn", str(main_id))
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
        cmd += "i/ref_rna_v2/geneset_venn "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="RefrnaV2_7320",
            task_type="1",
            submit_location="GenesetVenn",
            exp_level="T",
            geneset_ids="5b0686eca4e1af2636e97f7e,5b06899da4e1af3a72687201",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
