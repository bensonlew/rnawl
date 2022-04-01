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


class AsprofileAction(RefRnaV2Controller):

    def __init__(self):
        super(AsprofileAction, self).__init__(instant=False)
    @check_sig
    def POST(self):
        project_type = 'ref_rna_v2'
        db = Config().get_mongo_client(mtype=project_type, dydb_forbid=True)[Config().get_mongo_dbname(project_type, dydb_forbid=True)]


        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']


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
        sample_gtf = task_info['sample_gtf']
        genome_doc = db['sg_genome_db'].find_one({'genome_id': task_info['genome_id']})
        ref_gtf = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish', genome_doc['gtf'])
        ref_fa = os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish', genome_doc['dna_fa'])
        # ref_fa = task_info['ref_genome']
        # ref_gtf = task_info['ref_gtf']
        connect = self.db['sg_specimen_group']
        result = connect.find_one({'task_id':task_id})
        category_names = result["category_names"]
        specimen_names = result['specimen_names']
        group_dict = json.dumps(dict(zip(category_names, specimen_names)))

        def flat(nums):
            res = []
            for i in nums:
                if isinstance(i, list):
                    res.extend(flat(i))
                else:
                    res.append(i)
            return res

        specimen_names_list = flat(specimen_names)


        params_json = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            # 'group_dict': group_dict,
            # 'sample': specimen_names_list,
            'software': 'Asprofile'

        }

        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))


        name = "ASprofile" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        event_type = ['AE', 'IR_OFF', 'IR_ON', 'MIR_OFF', 'MIR_ON', 'MSKIP_OFF', 'MSKIP_ON', 'SKIP_OFF', 'SKIP_ON', 'TSS',
                      'TTS', 'XAE', 'XIR_OFF', 'XIR_ON', 'XMIR_OFF', 'XMIR_ON', 'XMSKIP_OFF', 'XMSKIP_ON', 'XSKIP_OFF', 'XSKIP_ON']
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            version="v3.1",
            desc='ASprofile main table',
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            sample=specimen_names_list,
            status="start",
            # event_type=event_type,
            group_dict=group_dict
        )

        main_id = self.ref_rna_v2.insert_main_table('sg_asprofile', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        options = {
            'ref_gtf': ref_gtf,
            "ref_fa": ref_fa,
            'gtf_dir': sample_gtf,
            'group_table': data.task_id,
            'sample_list':','.join(specimen_names_list),
            'update_info': json.dumps({str(main_id): "sg_asprofile"}),
            'main_id': str(main_id),
            'main_table_data': main_table_data,
        }

        # prepare to file
        task_name = 'ref_rna_v3.report.asprofile'
        to_files = ["ref_rna_v2.export_group_detail(group_table)",]

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            to_file=to_files,
                            new_task_id=new_task_id,
                            task_id=data.task_id)
        task_info = super(AsprofileAction, self).POST()
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
        cmd += "s/ref_rna_v3/asprofile "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_37259",

            task_type="2",
            submit_location="asprofile",




        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
