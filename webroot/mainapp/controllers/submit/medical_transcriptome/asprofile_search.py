# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from biocluster.config import Config
import datetime
import unittest
import json
import web
import os
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.medical_transcriptome import *
from mbio.api.to_file.medical_transcriptome import *
from bson.objectid import ObjectId
from collections import OrderedDict


class AsprofileSearchAction(MedicalTranscriptomeController):

    def __init__(self):
        super(AsprofileSearchAction, self).__init__(instant=False)
    @check_sig
    def POST(self):
        # project_type = 'medical_transcriptome'
        # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


        data = web.input()
        print data
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ["target_genes", 'geneset_asprofile', 'type',"id_type","sample"]


        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2901101', 'variables': variables}
                print info
                return json.dumps(info)

        record_counts = self.db["sg_asprofile_search"].find(
            {'task_id': data.task_id, "status": {"$in": ["end", "start"]}}).count()
        if record_counts >= 4:
            info = {"success": False, "info": u"已有四条运行记录，不可再投递"}
            return json.dumps(info)

        task_id = data.task_id
        task_info = self.medical_transcriptome.get_task_info(task_id=task_id)
        project_sn = task_info['project_sn']
        connect = self.db['sg_asprofile']
        result = connect.find_one({'task_id':task_id,'status':'end'})
        asprofile_path = result["asprofile_path"]

        params_json = {
            'task_id': data.task_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type),
            'geneset_asprofile': data.geneset_asprofile,
            'target_genes': data.target_genes,
            'type':  data.type,
            'id_type': data.id_type,
            'sample': data.sample
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
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))

        # self.medical_transcriptome.delete_search_records("sg_asprofile_search",data.task_id,3)#这个3的意思是最多保留三条
        name = "Asprofile_search" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            version="v1",
            desc='ASprofile Search main table',
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            status="start",
        )

        main_id = self.medical_transcriptome.insert_main_table('sg_asprofile_search', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        options = {
            'as_result_merge': asprofile_path,
            'target_genes':data.target_genes,
            "geneset_asprofile":data.geneset_asprofile,
            "type": data.type,
            "id_type": data.id_type,
            "sample": data.sample,
            "anno": task_id,
            'main_table_data': main_table_data,
            'update_info': json.dumps({str(main_id): "sg_asprofile_search"}),
            'main_id': str(main_id),
        }

        to_files = ["medical_transcriptome.get_gene_detail_whole(anno)",
                    "medical_transcriptome.export_multi_gene_list_search(geneset_asprofile)"]

        # prepare to file
        task_name = 'medical_transcriptome.report.asprofile_search'

        self.set_sheet_data(name=task_name,
                            options=options,
                            to_file=to_files,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)
        task_info = super(AsprofileSearchAction, self).POST()
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
        cmd += "s/medical_transcriptome/diff_asprofile "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="diff_asprofile",
            diff_type='group',
            group_id='12345678',
            filter="50",
            # sample='BY4741_1,Ni_BY_1,H4K5R_1',
            group_dict=json.dumps({"H1": ["H1581_1", "H1581_2", "H1581_3"], "H2": ["H1581_4","H1581_5", "H1581_6"]}).replace('"', '\\"')
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
