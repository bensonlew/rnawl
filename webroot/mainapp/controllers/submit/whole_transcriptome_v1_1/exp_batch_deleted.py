# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from biocluster.config import Config
import datetime
import unittest
import json
import web
import os
from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.whole_transcriptome import *
from mbio.api.to_file.whole_transcriptome import *
from bson.objectid import ObjectId


class ExpBatchDeletedAction(WholeTranscriptomeController):

    def __init__(self):
        super(ExpBatchDeletedAction, self).__init__(instant=False)
    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['library', 'level']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                print info
                return json.dumps(info)
        task_info = self.whole_transcriptome.get_task_info(task_id=data.task_id)
        project_sn = task_info['project_sn']
        params_json = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'library': data.library,
            'level': data.level
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        name = "Deleted_ExpBatch" + '_' + data.library + '_' + data.level + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            version="v1.1",
            desc='Deleted ExpBatch Table',
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            status="start",
        )

        main_id = self.whole_transcriptome.insert_main_table('sg_deleted_expbatch', main_info)
        # new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        # main_table_data = {'run_id': new_task_id}

        options = {
            'task_id': data.task_id,
            'update_info': json.dumps({str(main_id): "sg_deleted_expbatch"}),
            'main_id': str(main_id),
            'library': data.library,
            'level': data.level,
            # 'main_table_data': main_table_data,
        }

        # prepare to file
        task_name = 'whole_transcriptome_v1_1.report.deleted_expbatch'

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            # new_task_id=new_task_id,
                            task_id=data.task_id)
        task_info = super(ExpBatchDeletedAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        return json.dumps(task_info)

        # 更新基因集的使用信息
        # self.whole_transcriptome.insert_geneset_info(data.geneset_id, 'geneset_cluster', str(main_id))
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
        cmd += "s/whole_transcriptome_v1_1/exp_batch_deleted "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_248932",
            task_type="2",
            submit_location="expatch_deleted",
            library='long',
            level='T'




        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
