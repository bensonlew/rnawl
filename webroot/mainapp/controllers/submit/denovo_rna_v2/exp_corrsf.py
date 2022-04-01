# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpCorrsfAction(DenovoRnaV2Controller):
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
                info = {'success': False, 'info': "Lack argument: %s" % arg, 'code': 'C1603101', 'variables': variables}
                return json.dumps(info)
        exp_info = self.denovo_rna_v2.get_exp_params_info(data.exp_id, data.task_id)
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
        name = "ExpCorrsf" + '_' + data.corr_way + '_' + data.padjust_way + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            version="v2",
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='express correlation analysis main table',
            params=params,
            status="start"
        )
        # 主表名字为了区分样本相关性分析,这个是做6-15(six-fifteen), 所以叫做sg_exp_corrsf
        main_id = self.denovo_rna_v2.insert_main_table('sg_exp_corrsf', main_info)
        sg_task_info = self.denovo_rna_v2.get_task_info(data.task_id)
        if "version" in sg_task_info:
             version=sg_task_info["version"]
        else:
            version="v1"
        if version == "v2":
            result_dir = self.denovo_rna_v2.get_annotation_stat_info(data.task_id)['result_dir']
            print(result_dir)
            result_dir = os.path.join(result_dir, 'all_annot.xls')
        else:
            result_dir=self.denovo_rna_v2.get_annotation_stat_info(data.task_id)['result_dir']
            result_dir=os.path.join(result_dir.replace("Annotation","AnnoQuery"), 'transcript_anno_detail.xls')

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
            "update_info": json.dumps({str(main_id): "sg_exp_corrsf"})  # to update sg_status
        }
        if "qvalue" in data:
            options.update({'sig_type':1, 'qvalue_cutoff': data.qvalue})
        if "pvalue" in data:
            options.update({'sig_type':0, 'pvalue_cutoff': data.pvalue})
        # prepare to file
        to_files = ["denovo_rna_v2.export_geneset_exp_matrix(exp_matrix)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'denovo_rna_v2.report.exp_corrsf'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpCorrsfAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        self.denovo_rna_v2.insert_geneset_info(data.geneset_id, "sg_exp_corrsf", str(main_id))
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
        cmd += "s/denovo_rna_v2/exp_corrsf "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_35544",
            task_type="2",
            submit_location="genesetcorrsf",
            exp_id="5d79928017b2bf1d40fae25a",
            group_dict=r'{"A":["A1", "A2", "A3"],"B": [ "B1", "B2", "B3"], "C": [ "C1", "C2", "C3"]}'.replace('"', '\\"'),
            group_id="5d798fb617b2bf1d40ee5486",
            geneset_id="5d7992da17b2bf1d40ffba33",
            #pvalue="0.05",
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
