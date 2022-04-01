# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig
import unittest
import os


class GenesetVennAction(MedicalTranscriptomeController):
    def __init__(self):
        super(GenesetVennAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['geneset_id', 'level']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:  %s" % arg, 'code': 'C2901901', 'variables': variables}
                return json.dumps(info)
        sg_task = self.medical_transcriptome.get_task_info(data.task_id)
        project_sn = sg_task["project_sn"]
        task_id = data.task_id
        # create main table record
        # geneset_id_list = str(data.geneset_ids).split(',')
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            level=data.level,
            # geneset_ids=geneset_id_list,
            geneset_id=data.geneset_id,
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "GenesetVenn" + '_' + data.level + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v1",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Geneset venn analysis main table',
            params=params,
            status="end"
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_geneset_venn', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        options = {
            "task_id": data.task_id,
            "main_id": str(main_id),
            "geneset_id": data.geneset_id,
            "submit_location": data.submit_location,
            "main_table_data": main_table_data,
            "update_info": json.dumps({str(main_id): "sg_geneset_venn"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.geneset_venn'
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
        self.medical_transcriptome.insert_geneset_info(data.geneset_id, "sg_geneset_venn", str(main_id))
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
        cmd += "i/medical_transcriptome/geneset_venn "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="1",
            submit_location="GenesetVenn",
            level="G",
            geneset_ids="5f642b3117b2bf6711b00244,5f642b1f17b2bf6711afc4e1,5f642b1a17b2bf6711afc0ef",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
