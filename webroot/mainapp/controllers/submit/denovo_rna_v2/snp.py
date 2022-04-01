# -*- coding: utf-8 -*-
from __future__ import print_function
import web
import json
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from bson import SON
import datetime
import unittest
from mbio.api.to_file.denovo_rna_v3 import *
import os
from biocluster.config import Config

class SnpAction(DenovoRnaV2Controller):
    def __init__(self):
        super(SnpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['submit_location', 'task_id', 'task_type']
        print(data)
        print("where is method")
        for param in params_name:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!!" % param, "code": 'C1601901', "variables": var}
                return json.dumps(info)


        # 测试接口的时候如果没有写对应的workflow文件，就会status为failed,然后更新desc的描述
        task_name = 'denovo_rna_v2.report.snp'
        task_type = 'workflow'
        # 'no'的目的是为了获取基础workflow跑出来的结果，而不是交互的,get_snp_info到sg_snp主表去获取基础流程的信息
        task_info_snp = self.denovo_rna_v2.get_snp_info(data.task_id, 'no')
        if not task_info_snp:
            info = {"success": False, "info": "因为找不到snp的call_vcf_path", "code": 'C1601904', "variables": ''}
            return json.dumps(info)
        main_table_name = 'snp_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        # get_task_info根据task_id到sg_task获得相关记录信息
        sg_task_info = self.denovo_rna_v2.get_task_info(data.task_id)
        if "version" in sg_task_info:
             version=sg_task_info["version"]
        else:
            version="v1"
        if version == "v1":
            if float(data.qual) <= 0:
                var = []
                var.append(data.qual)
                info = {"success": False, "info": "qual的值是%s，不在规定范围内!"%(data.qual), "code": 'C1601902', "variables": var}
                return json.dumps(info)
            if int(data.dp) < 1:
                var = []
                var.append(data.dp)
                info = {"success": False, "info": "dp的值是%s，不在规定范围内!"%(data.dp), "code": 'C1601903', "variables": var}
                return json.dumps(info)
            '''params_json里面的内容首先需要主表里面的运行的params，然后还有一些前端关注的参数，最后是比较排序后的字符串'''
            if float(data.qual) == int(data.qual):
                params_json = {
                    'qual': str(int(data.qual)),
                    'dp': int(data.dp),
                    'submit_location': data.submit_location,
                    'task_type': int(data.task_type),
                    'task_id' : data.task_id
                }
            else:
                params_json = {
                    'qual': str(float(data.qual)),
                    'dp': int(data.dp),
                    'submit_location': data.submit_location,
                    'task_type': int(data.task_type),
                    'task_id' : data.task_id
                }
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))

            '''mongo_data里面的内容就是和主表里面的内容保持一致'''
            call_vcf_path, bamlist_num = self.denovo_rna_v2.get_callvcf_path(data.task_id, 'no')
            mongo_data = [
                ("call_vcf_path", call_vcf_path),
                ('project_sn', sg_task_info['project_sn']),
                ('task_id', data.task_id),
                ('status', 'start'),
                ('desc', 'snp结果表'),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ('params', params),
                ('name', main_table_name),
                ('sample_names', task_info_snp['sample_names']),
            ]

            collection_name = 'sg_snp'
            main_table_id = self.denovo_rna_v2.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): 'sg_snp'}
            # 从sg_task的mongo表获取bamlist的路径
            #bamlist_path = self.denovo_rna_v2.get_bamlist_path("denovo_rna_v2")

            '''options里面的参数一部分是用于tool，一部分是用于to_file'''
            options = {
                'bamlist': bamlist_num,
                'call_vcf': call_vcf_path,
                'qual': float(data.qual),
                'dp': int(data.dp),
                "update_info": json.dumps(update_info),
                "snp_id": str(main_table_id),
            }

            self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                                task_id=task_info_snp['task_id'], project_sn=sg_task_info['project_sn'],
                                module_type=task_type, params=params, to_file=None)

            task_info = super(SnpAction, self).POST()

        else :
            task_name = 'denovo_rna_v3.report.snp'
            '''params_json里面的内容首先需要主表里面的运行的params，然后还有一些前端关注的参数，最后是比较排序后的字符串'''
            params_json = {
                    'method': data.method,
                    'submit_location': data.submit_location,
                    'task_type': int(data.task_type),
                    'task_id': data.task_id
                }
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))

            '''mongo_data里面的内容就是和主表里面的内容保持一致'''
#            call_vcf_path, bamlist_num = self.denovo_rna_v2.get_callvcf_path(data.task_id, 'no')
            mongo_data = [
#                ("call_vcf_path", call_vcf_path),
                ('project_sn', sg_task_info['project_sn']),
                ('version','v2'),
                ('task_id', data.task_id),
                ('status', 'start'),
                ('desc', 'snp结果表'),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ('params', params),
                ('name', main_table_name),
                ('sample_names', task_info_snp['sample_names']),
            ]

            collection_name = 'sg_snp'
            main_table_id = self.denovo_rna_v2.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): 'sg_snp'}
            # 从sg_task的mongo表获取bamlist的路径
            # bamlist_path = self.denovo_rna_v2.get_bamlist_path("denovo_rna_v2")

            to_files = ["denovo_rna_v3.get_after_qc_bam_path(bamlist)"]
            '''options里面的参数一部分是用于tool，一部分是用于to_file'''
            project_type = 'denovo_rna_v2'
            db_ = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
            connect = db_['sg_specimen_group']
            record = connect.find_one({'task_id': data.task_id})
            samples = record['specimen_names']
            sample_list = sum(samples, [])
            sample_list_str = ';'.join(sample_list)
            options = {
                "sample_list_str": sample_list_str,
                "bamlist":data.task_id,
                "update_info": json.dumps(update_info),
                "snp_id": str(main_table_id),
                "call_type":data.method,
                "cds_bed":sg_task_info["bedpath"],
                "allt2g":sg_task_info["assemble_t2g"],
                'ref_fasta':sg_task_info["unigene_fa"],
                'anno':os.path.join(self.denovo_rna_v2.get_annotation_stat_info(data.task_id)["result_dir"],'all_annot.xls')
            }

            self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                                task_id=task_info_snp['task_id'], project_sn=sg_task_info['project_sn'],
                                module_type=task_type,to_file=to_files,params=params)

            task_info = super(SnpAction, self).POST()

        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        else:
            pass
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/denovo_rna_v2/snp "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_35544",
            task_type="2",
            submit_location="snp",
            method="gatk"
            # exp_id="5d41458917b2bf10c8b3ffb0",
            # group_id="5d47f30b17b2bf08c619bf78",
            # exp_level="trans",
            # group_dict=r'{"A":["A1", "A2","A3"],"B": [ "B1", "B2","B3"],"C":["C1","C2","C3"]}'.replace('"', '\\"'),
            # control_id="5d47f30b17b2bf08c619bf79",
            # diff_method='DESeq2',
            # stat_type='padjust',
            # stat_cutoff='0.05',
            # fc='1.5',
            # filter_method="min",
            # tpm_filter_threshold="5"
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)
