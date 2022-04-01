# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from mbio.api.to_file.prok_rna import *
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId


class SnpAction(ProkRNAController):
    def __init__(self):
        super(SnpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['method_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': "Lack argument: {}".format(arg)}
                return json.dumps(info)
            if arg.lower() == "null":
                info = {'success': False, 'info': "{} : is null or NULL".format(arg)}
                return json.dumps(info)

        exp_info = self.prok_rna.get_task_info(task_id=data.task_id)
        # id2name = self.prok_rna.get_annotation_stat_info(task_id=data.task_id)["result_dir"] + '/anno_stat/all_anno_detail.xls'
        id2name = exp_info['id2name']
        # 我们到时取前三列，大工作流也需要传入这个文件
        project_sn = exp_info["project_sn"]
        ref_gtf = exp_info["ref_gtf"]
        ref_genome_custom = exp_info["ref_genome"]
        fastq = exp_info["fastq"]
        task_id = data.task_id
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            method_type=data.method_type
        )
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "SNP" + '_' + data.method_type + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        # prepare main table info
        main_info = dict(
            name=name,
            task_id=task_id,
            project_sn=project_sn,
            desc='snp_indel analysis main table',
            status="start",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
        )
        main_id = self.prok_rna.insert_main_table('sg_snp', main_info)

        # prepare option for workflow
        options = {
            "bamlist": data.task_id,
            "fq_list": fastq,
            "method_type": data.method_type,
            'ref_gtf': self.use_s3(ref_gtf),
            'ref_genome_custom': self.use_s3(ref_genome_custom),
            "snp_main_id": str(main_id),
            "id2name": self.use_s3(id2name),
            "update_info": json.dumps({str(main_id): "sg_snp"})  # to update sg_status
        }
        # prepare to file
        to_files = ["prok_rna.get_after_qc_bam_path(bamlist)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'prok_rna.report.snp'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(SnpAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        # task_info['group_dict'] = group_dict
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.prok_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.prok_rna.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/prok_rna/snp "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="prok_rna_srna",
            task_type="2",
            submit_location="snp",
            method_type="gatk",
            fq_list="/mnt/ilustre/users/sanger-dev/workspace/20180806/Single_HiseqQc_3546/HiseqQc/output/sickle_dir/fq_list.txt",
            ref_genome="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_prok_rna/test_data1/ref/GCF_000009345.1_ASM934v1_genomic.fna",
            ref_gtf="/mnt/ilustre/users/sanger-dev/workspace/20180920/Prokrna_tsg_32069/output/rock_index/reshape.gtf",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
