# -*- coding: utf-8 -*-
import web
import json
import datetime
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os,re
from bson.objectid import ObjectId
from biocluster.config import Config


class SnpAction(RefRnaV2Controller):
    def __init__(self):
        super(SnpAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        project_type = 'ref_rna_v2'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['method_type']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C3400203', 'variables': variables}
                return json.dumps(info)
            if arg.lower() == "null":
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "%s : is null or NULL" % arg, 'code': 'C3400204', 'variables': variables}
                return json.dumps(info)

        exp_info = self.ref_rna_v2.get_task_info(task_id=data.task_id)

        des, des_type = self.ref_rna_v2.get_des_type(task_id=data.task_id)
        # desc_type_info = self.ref_rna_v2.get_annotation_stat_info(task_id=data.task_id)
        # species_name = desc_type_info['species_name']
        project_sn = exp_info["project_sn"]
        ref_gtf = exp_info["ref_gtf"]
        if not os.path.exists(ref_gtf) and re.match(r'/mnt/ilustre/', ref_gtf):
            ref_gtf = exp_info["ref_gtf"].replace('ilustre','lustre')
        ref_genome_custom = exp_info["ref_genome"]
        if not os.path.exists(ref_genome_custom) and re.match(r'/mnt/ilustre/', ref_genome_custom):
            ref_genome_custom = exp_info["ref_genome"].replace('ilustre','lustre')
        if not os.path.exists(des) and re.match(r'/mnt/ilustre/', des):
            des = des.replace('ilustre','lustre')
        fastq = exp_info["fastq"]
        task_id = data.task_id
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            #project_sn=project_sn,
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
            version="v3",
            desc='snp_indel analysis main table',
            status="start",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
        )
        main_id = self.ref_rna_v2.insert_main_table('sg_snp', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(task_id)
        main_table_data = {'run_id': new_task_id}

        connect = db['sg_specimen_group']
        record = connect.find_one({'task_id': data.task_id})
        samples = record['specimen_names']
        sample_list = sum(samples, [])
        sample_list_str = ';'.join(sample_list)
        # prepare option for workflow
        options = {
            "sample_list_str": sample_list_str,
            "bamlist": data.task_id,
            "fq_list": self.use_s3(fastq),
            "method_type": data.method_type,
            'ref_gtf': self.use_s3(ref_gtf),
            'ref_genome_custom': self.use_s3(ref_genome_custom),
            "snp_main_id": str(main_id),
            # "species_name": species_name,
            "des": self.use_s3(des),
            "des_type": des_type,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_snp"})  # to update sg_status
        }
        # prepare to file
        to_files = ["ref_rna_v2.get_after_qc_bam_path(bamlist)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'ref_rna_v3.report.snp'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
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
            _ = self.ref_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.ref_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py  '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna_v3/snp "
        cmd += "-b http://192.168.12.101:9090 "
        args = dict(
            task_id="tsg_33912",
            task_type="2",
            submit_location="snp",
            method_type="samtools",
            fq_list="/mnt/ilustre/users/sanger-dev/workspace/20190419/Refrna_tsg_33912/HiseqQc/output/sickle_dir/fq_list.txt",
            ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
            ref_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
