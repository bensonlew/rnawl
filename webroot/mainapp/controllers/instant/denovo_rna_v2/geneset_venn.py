# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class GenesetVennAction(DenovoRnaV2Controller):
    def __init__(self):
        super(GenesetVennAction, self).__init__(instant=True)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['geneset_ids', 'exp_level']
        # check arg
        for argu in basic_args:
            if not hasattr(data, argu):
                var = []
                var.append(argu)
                info = {'success': False, 'info': "Lack argument: %s" % (argu), "code": 'C1601601', "variables": var}
                return json.dumps(info)
        sg_task = self.denovo_rna_v2.get_task_info(data.task_id)
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
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "GenesetVenn" + '_' + data.exp_level + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            version="v2",
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Geneset venn analysis main table',
            params=params,
            status="end"
        )
        main_id = self.denovo_rna_v2.insert_main_table('sg_geneset_venn', main_info)

        # prepare option for workflow
        options = {
            "task_id": data.task_id,
            "main_id": str(main_id),
            "geneset_id": data.geneset_ids,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_geneset_venn"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'denovo_rna_v2.report.geneset_venn'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
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
        self.denovo_rna_v2.insert_geneset_info(data.geneset_ids, "sg_geneset_venn", str(main_id))
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
        cmd += "i/denovo_rna_v2/geneset_venn "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="denovo_rna_v2_upgrade",
            task_type="1",
            submit_location="GenesetVenn",
            exp_level="T",
            geneset_ids="5cb64c0d17b2bf6187b11fbd,5cb64c1c17b2bf6187b1d8a4",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
