# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
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
from collections import OrderedDict


class DiffAsprofileAction(RefRnaV2Controller):

    def __init__(self):
        super(DiffAsprofileAction, self).__init__(instant=False)
    @check_sig
    def POST(self):
        project_type = 'ref_rna_v2'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ["diff_type", 'group_id', 'group_dict']


        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                print info
                return json.dumps(info)

        task_id = data.task_id
        task_info = self.ref_rna_v2.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        connect = db['sg_asprofile']
        result = connect.find_one({'task_id':task_id,'status':'end'})
        asprofile_path = result["asprofile_path"]





        params_json = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'diff_type': data.diff_type,
            'group_id': data.group_id,
            'group_dict': json.loads(data.group_dict)

        }
        # if hasattr(data,'sample'):
        #     params_json.update({
        #         "sample": data.sample
        #     })
        # if hasattr(data, 'group_dict'):
        #     group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        #     params_json.update({
        #         'group_dict': group_dict
        #     })
        if hasattr(data, 'filter'):
            params_json.update({
                'filter': data.filter
            })
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))


        name = "Diff_ASprofile" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if data.diff_type == 'sample':
            group_dict = json.loads(data.group_dict)
            sample_list = list()
            for s in group_dict:
                sample_list.extend(group_dict[s])
        if data.diff_type == 'group':
            group_dict = json.loads(data.group_dict)
            sample_list = [x for x in group_dict]
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            version="v3.1",
            desc='Diff_ASprofile main table',
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            sample_group = sample_list,
            status="start",
        )

        main_id = self.ref_rna_v2.insert_main_table('sg_asprofile_diff', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        options = {
            'AS_result_merge': asprofile_path,
            'update_info': json.dumps({str(main_id): "sg_asprofile_diff"}),
            'main_id': str(main_id),
            'main_table_data': main_table_data,
        }
        if data.diff_type == 'sample':
            group_dict = json.loads(data.group_dict)
            sample_list = list()
            for s in group_dict:
                sample_list.extend(group_dict[s])
            sample = ','.join(sample_list)

            options.update({
                "sample": sample
            })
        if data.diff_type == 'group':
            options.update({
                "group_dict": data.group_dict
            })
        # if hasattr(data, 'group_dict'):
        #     options.update({
        #         "group_dict": data.group_dict
        #     })
        if hasattr(data, 'filter'):
            options.update({
                "filter": data.filter
            })
        # prepare to file
        task_name = 'ref_rna_v3.report.diff_asprofile'

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)
        task_info = super(DiffAsprofileAction, self).POST()
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
            _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/ref_rna_v3/diff_asprofile "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_37259",

            task_type="2",
            submit_location="diff_asprofile",
            diff_type='group',
            group_id='12345678',
            filter="50",
            # sample='BY4741_1,Ni_BY_1,H4K5R_1',
            group_dict=json.dumps({"H4K5R": ["H4K5R_1", "H4K5R_2", "H4K5R_3"], "WT": ["BY4741_1","BY4741_2", "BY4741_3"]}).replace('"', '\\"')




        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
