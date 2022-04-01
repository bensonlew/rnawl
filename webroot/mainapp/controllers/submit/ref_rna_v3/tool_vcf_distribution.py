# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.models.mongo.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ToolVcfDistributionAction(RefRnaV2Controller):
    def __init__(self):
        super(ToolVcfDistributionAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['snp_id', 'tool_type']
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
        if 'ref_gtf' in sg_task:
            ref_gtf = sg_task['ref_gtf']
            ref_path = os.path.join(os.path.dirname(ref_gtf), 'assembly_level.txt')
        else:
            ref_path = ''
        s3_path = sg_task['assemble_fa']
        s3_dir = s3_path.split('/workflow_results')[0]
        # create main table record
        params = dict(
            task_id=data.task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            snp_id=data.snp_id,
            tool_type=data.tool_type
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Vcf_Distribution" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=data.task_id,
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='',
            params=params,
            submit_type=submit_type,
            status="start"
        )

        snp_main = self.ref_rna_v2.get_main_info(data.snp_id, 'sg_snp', data.task_id)
        if 'version' in snp_main:
            version = snp_main['version'].split('v')[1]
        else:
            version = 2.9
        if 'run_id' in snp_main or float(version) < 3.0:
            info = {'success': False, 'info': "暂不支持该SNP分析结果进行变异位点染色体分布图分析。"}
            return json.dumps(info)
        else:
            vcf_path = os.path.join(s3_dir, 'workflow_results', '08SNP', 'SNP_vcf', 'final.vcf')

        main_id = self.ref_rna_v2.insert_main_table('sg_tool_lab_vcf_distribution', main_info)
        # prepare option for workflow
        options = {
            "vcf_file": vcf_path,
            'tool_type': data.tool_type,
            'ref_path': ref_path,
            'params': params,
            "main_id": str(main_id),
            "relate_name": name,
            'task_id': data.task_id,
            "submit_location": data.submit_location,
            "update_info": json.dumps({str(main_id): "sg_tool_lab_vcf_distribution"})  # to update sg_status
        }

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        to_files = []
        task_name = 'ref_rna_v3.report.tool_vcf_distribution'
        self.set_sheet_data(name=task_name,
                            to_file=to_files,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ToolVcfDistributionAction, self).POST()
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
        cmd += "s/ref_rna_v3/tool_vcf_distribution "
        cmd += "-b http://wpm2.sanger.com "
        args = dict(
            task_id="sg_249157",
            task_type="2",
            submit_location="VcfDistribution",
            snp_id="5fffd24117b2bf0ec5e56c2c",
            tool_type='Vcf_Distribution'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()