# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import datetime
import os
import re
import unittest
import json
import web
from biocluster.config import Config
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mainapp.libs.signature import check_sig
from mbio.api.to_file.ref_rna_v2 import *


class GeneFusionVennAction(MedicalTranscriptomeController):
    def __init__(self):
        super(GeneFusionVennAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['use_group','group_dict','group_id',"fusion_id"]
        if hasattr(data, 'filter_threshold'):
            basic_args += ["filter_threshold"]
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2902601', 'variables': variables}
                return json.dumps(info)
            if arg.lower() == "null":
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "%s : is null or NULL" % arg, 'code': 'C2902602',
                        'variables': variables}
                return json.dumps(info)

        # project_sn = "ref_rna_v2"
        task_info = self.medical_transcriptome.get_task_info(task_id=data.task_id)
        project_sn = task_info["project_sn"]
        task_id = data.task_id
        # create main table record
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            use_group = data.use_group,
            group_dict = group_dict,
            fusion_id = data.fusion_id,
            group_id = data.group_id
        )
        if hasattr(data, 'filter_threshold'):
            params.update({
                'filter_threshold': data.filter_threshold
            })
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Fusion_result" + '_' + 'Venn' + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        # prepare main table info

        main_info = dict(
            name=name,
            task_id=task_id,
            project_sn=project_sn,
            desc='fusion_venn analysis main table',
            status="start",
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            params=params,
            version = "v1"
        )
        main_id = self.medical_transcriptome.insert_main_table('sg_gene_fusion_venn', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}

        # prepare option for workflow
        # group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        if str(data.group_id).lower() == 'all':
            samples = group_dict['all']
            group_dict = OrderedDict([(x, [x]) for x in samples])
        print(group_dict)

        options = {
            "group": json.dumps(group_dict),
            "group_dict": json.dumps(group_dict),
            "fusion_id" :data.fusion_id,
            "fusion_venn_main_id": str(main_id),
            "use_group":data.use_group,
            "fusion_matrix":data.fusion_id,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_gene_fusion_venn"})  # to update sg_status
        }
        if hasattr(data, 'filter_threshold'):
            options.update({"filter_threshold":data.filter_threshold})

        # prepare to file
        to_files = ["medical_transcriptome.export_fusion_matrix(fusion_matrix)",
                    'medical_transcriptome.export_group(group)']


        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.fusion_venn'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GeneFusionVennAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
            }
        }
        task_info['group_dict'] = group_dict
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
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py  '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/ref_rna_v3/gene_fusion_venn "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_37575",
            task_type="2",
            submit_location="gene_fusion",
            method_type="sentieon",
            fusion_id = "5ee9a2c417b2bf1648706db9",
            group_id = "suibian",
            # group_dict = json.dumps({"A": ["A1", "A2", "A3"], "B": ["B1","B2", "B3"]}).replace('"', '\\"'),
            group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"]}'.replace('"', '\\"'),
            use_group = "yes",
            filter_threshold ="50"
            # fq_list="/mnt/ilustre/users/sanger-dev/workspace/20190419/Refrna_tsg_33912/HiseqQc/output/sickle_dir/fq_list.txt",
            # ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
            # ref_genome="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
