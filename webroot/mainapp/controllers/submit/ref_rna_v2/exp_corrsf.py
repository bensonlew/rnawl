# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.ref_rna_v2_controller import RefRnaV2Controller
from mbio.api.to_file.ref_rna_v2 import *
from biocluster.file import getsize, exists
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpCorrsfAction(RefRnaV2Controller):
    def __init__(self):
        super(ExpCorrsfAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['geneset_id', 'exp_id',
                       'group_dict', 'group_id',
                       'cor_cutoff', 'corr_way', 'exp_level',
                        'padjust_way']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C2900501', 'variables': variables}
                return json.dumps(info)
        exp_info = self.ref_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
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
            exp_level=data.exp_level,
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
            version="v3",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='express correlation analysis main table',
            params=params,
            status="start"
        )
        # 主表名字为了区分样本相关性分析,这个是做6-15(six-fifteen), 所以叫做sg_exp_corrsf
        main_id = self.ref_rna_v2.insert_main_table('sg_exp_corrsf', main_info)
        new_task_id = self.ref_rna_v2.get_new_id(data.task_id)
        main_table_data = {'run_id': new_task_id}
        result_dir = self.ref_rna_v2.get_annotation_stat_info(data.task_id)['result_dir']
        # if os.path.exists(result_dir + "/allannot_class/all_annot.xls"):
        #     result_dir = result_dir + "/allannot_class/all_annot.xls"
        # else:
        #     result_dir = result_dir + "/refannot_class/all_annot.xls"
        if exists(result_dir + "/allannot_class/all_annot.xls"):
            result_dir = result_dir + "/allannot_class/all_annot.xls"
        else:
            result_dir = result_dir + "/refannot_class/all_annot.xls"
        # prepare option for workflow
        options = {
            "exp_matrix": data.exp_id+";"+data.geneset_id,
            "corr_main_id": str(main_id),
            'cor_cutoff': data.cor_cutoff,
            'gt': data.exp_level,
            'group_dict': data.group_dict,
            'group_id': data.group_id,
            'padjust_way': data.padjust_way,
            'corr_way': data.corr_way,
            'anno': result_dir,
            'main_table_data': main_table_data,
            "update_info": json.dumps({str(main_id): "sg_exp_corrsf"})  # to update sg_status
        }
        if "qvalue" in data:
            options.update({'sig_type':1, 'qvalue_cutoff': data.qvalue})
        if "pvalue" in data:
            options.update({'sig_type':0, 'pvalue_cutoff': data.pvalue})
        # prepare to file
        to_files = ["ref_rna_v2.export_geneset_exp_matrix(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'ref_rna_v2.report.exp_corrsf'
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
        geneset_info = self.ref_rna_v2.insert_geneset_info(data.geneset_id, "sg_exp_corrsf", str(main_id))
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
        cmd += "s/ref_rna_v2/exp_corrsf "
        cmd += "-b http://192.168.12.102:9090 "
        args = dict(
            task_id="tsg_30848",
            task_type="2",
            submit_location="genesetcorrsf",
            exp_id="5b36ca73a4e1af041a95ecda",
            group_dict={"A1":["A1_1","A1_2","A1_3"],"A2":["A2_1","A2_2","A2_3"],"B1":["B1_1","B1_2","B1_3"]},
            group_id="5b36a4e9a4e1af041a7a4534",
            geneset_id="5b36cafda4e1af041a9b9d75",
            pvalue="0.05",
            qvalue="0.05",
            cor_cutoff="0.5",
            exp_level="G",
            sig_type="qvalue",
            corr_way="spearman",
            padjust_way="fdr_bh",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join([str(x) for x in arg_names]), ";".join([str(x) for x in arg_values]))
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
