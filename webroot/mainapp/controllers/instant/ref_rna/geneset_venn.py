# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os


class GenesetVennAction(RefRnaController):
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
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
        sg_task = self.ref_rna.get_task_info(data.task_id)
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
            task_id=task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='Geneset venn analysis main table',
            params=params,
            status="end"
        )
        main_id = self.ref_rna.insert_main_table('sg_geneset_venn', main_info)

        # 运行workflow 并传回参数
        task_info = dict()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                },
            'info': "成功存储参数，可进行venn分析"
        }
        # 更新基因集的使用信息
        self.ref_rna.insert_geneset_info(data.geneset_ids, name, str(main_id))
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
        cmd += "i/ref_rna/geneset_venn "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="test_refrna_genesetvenn",
            task_type="1",
            submit_location="GenesetVenn",
            exp_level="G",
            geneset_ids="5a35da8ba4e1af38e6fb83d7,5a35da8ba4e1af38e6fb83d9,5a35da8ba4e1af38e6fb83db",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
