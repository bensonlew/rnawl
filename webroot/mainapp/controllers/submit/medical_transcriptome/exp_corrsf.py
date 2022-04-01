# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.medical_transcriptome_controller import MedicalTranscriptomeController
from mbio.api.to_file.medical_transcriptome import *
from mainapp.libs.signature import check_sig
import unittest
import os
from bson.objectid import ObjectId
from biocluster.config import Config


class ExpCorrsfAction(MedicalTranscriptomeController):
    def __init__(self):
        super(ExpCorrsfAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['geneset_id', 'exp_id',
                       'group_dict', 'group_id',
                       'cor_cutoff', 'corr_way', 'level',
                        'padjust_way']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2900501', 'variables': variables}
                return json.dumps(info)
        exp_info = self.medical_transcriptome.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            geneset_id=data.geneset_id,
            cor_cutoff=data.cor_cutoff,
            level=data.level,
            padjust_way=data.padjust_way,
            corr_way=data.corr_way,
            group_dict=group_dict,
            sig_type=data.sig_type
        )
        if "qvalue" in data:
            params.update({'qvalue': data.qvalue})
        if "pvalue" in data:
            params.update({'pvalue': data.pvalue})
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "GeneSet_ExpCorr" + '_' + data.corr_way + '_' + data.padjust_way + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v1",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='express correlation analysis main table',
            params=params,
            status="start"
        )
        # 主表名字为了区分样本相关性分析,这个是做6-15(six-fifteen), 所以叫做sg_exp_corrsf
        main_id = self.medical_transcriptome.insert_main_table('sg_exp_corrsf', main_info)
        new_task_id = self.medical_transcriptome.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}
        result_dir = self.medical_transcriptome.get_annotation_stat_info(data.task_id)['result_dir']
        if data.level == "G":
            anno_name = 'all_annot_gene.xls'
        elif data.level == "T":
            anno_name = 'all_annot_tran.xls'
        anno_file = os.path.join(result_dir, "allannot_class", anno_name)
        if os.path.exists(anno_file):
            result_dir = anno_file
        else:
            result_dir = os.path.join(result_dir, "refannot_class", anno_name)

        project_type = 'medical_transcriptome'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        connect_exp = db['sg_exp']
        record_exp = connect_exp.find_one({'task_id': data.task_id, 'main_id': ObjectId(data.exp_id)})
        if not record_exp:
            self.bind_object.set_error("意外错误，main_id:%s的背景基因在sg_exp中未找到！", variables=(data.exp_id,), )
        is_rmbe = str(record_exp['is_rmbe']).lower()
        if is_rmbe == 'false':
            exp_id = data.exp_id
        if is_rmbe == 'true':
            exp_id = str(record_exp['batch_main_id'])
        # prepare option for workflow
        options = {
            "exp_matrix": exp_id+";"+data.geneset_id + ";" + data.level+";"+is_rmbe,
            "corr_main_id": str(main_id),
            'cor_cutoff': data.cor_cutoff,
            'gt': data.level,
            'group_dict': data.group_dict,
            'group_id': data.group_id,
            'padjust_way': data.padjust_way,
            'corr_way': data.corr_way,
            'anno': result_dir,
            'task_id': data.task_id,
            # "anno":'/mnt/ilustre/users/sanger-dev/workspace/20200807/Refrna_tsg_38276/intermediate_results/Express/ExpAnnalysis/gene.count.matrix.annot.xls',
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_exp_corrsf"})  # to update sg_status
        }
        if "qvalue" in data:
            options.update({'sig_type': 1, 'qvalue_cutoff': data.qvalue})
        if "pvalue" in data:
            options.update({'sig_type': 0, 'pvalue_cutoff': data.pvalue})
        # prepare to file
        to_files = ["medical_transcriptome.export_geneset_exp_matrix_new(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'medical_transcriptome.report.exp_corrsf'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            new_task_id=new_task_id,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpCorrsfAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        geneset_info = self.medical_transcriptome.insert_geneset_info(data.geneset_id, "sg_exp_corrsf", str(main_id))
        # if 'group_id' in data and str(data.group_id).lower() != 'all':
        #     _ = self.medical_transcriptome.update_group_is_use(data.task_id, data.group_id)
        # if 'control_id' in data:
        #     _ = self.medical_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/medical_transcriptome/exp_corrsf "
        cmd += "-b http://bcl.tsg.com "
        # cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="medical_transcriptome",
            task_type="2",
            submit_location="genesetcorrsf",
            exp_id="5f6db9d017b2bf3607b2bf32",
            group_dict=r'{"H2": ["H1581_4", "H1581_5", "H1581_6"], "H3": ["H1581_7", "H1581_8", "H1581_9"], '
                       r'"S1": ["SNU16_1", "SNU16_2", "SNU16_3"], "H1": ["H1581_1", "H1581_2", "H1581_3"], '
                       r'"S3": ["SNU16_7", "SNU16_8", "SNU16_9"], "S2": ["SNU16_4", "SNU16_5", "SNU16_6"]}'.replace('"', '\\"'),
            group_id="5f703bec17b2bf74b39b635c",
            geneset_id="5f8e6d4117b2bf609322ee78",
            pvalue="0.05",
            qvalue="0.05",
            cor_cutoff="0.5",
            level="G",
            sig_type="qvalue",
            corr_way="spearman",
            padjust_way="fdr_by",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join([str(x) for x in arg_names]), ";".join([str(x) for x in arg_values]))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()


