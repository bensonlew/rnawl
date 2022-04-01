# -*- coding: utf-8 -*-

import datetime
import unittest
from bson.objectid import ObjectId
import web

from mainapp.controllers.project.whole_transcriptome_controller import WholeTranscriptomeController
from mainapp.libs.signature import check_sig
from mbio.api.to_file.whole_transcriptome.geneset import *


class ExpCorrsfAction(WholeTranscriptomeController):
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
        # exp = self.whole_transcriptome.get_exp_id(data.level, data.task_id)
        # exp_id = exp['main_id'].__str__()
        exp_dict = self.whole_transcriptome.get_main_info_by_record("exp", level=data.level,
                                                                      task_id=data.task_id, main_id=ObjectId(data.exp_id))
        is_rmbe = str(exp_dict['is_rmbe'])
        if 'is_rmbe' not in exp_dict or is_rmbe == 'false':
            exp_id = str(exp_dict['main_id'])
        if is_rmbe == 'true':
            exp_id = str(exp_dict['batch_main_id'])
        exp_info = self.whole_transcriptome.get_exp_params_info(data.exp_id, data.task_id)
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
            name=name,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            desc='express correlation analysis main table',
            params=params,
            status="start",
            version='v1.1'
        )
        # 主表名字为了区分样本相关性分析,这个是做6-15(six-fifteen), 所以叫做exp_corrsf
        main_id = self.whole_transcriptome.insert_main_table('exp_corrsf', main_info)
        # result_dir = self.whole_transcriptome.get_annotation_stat_info(data.task_id)['result_dir']
        # if os.path.exists(result_dir + "/allannot_class/all_annot.xls"):
        #     result_dir = result_dir + "/allannot_class/all_annot.xls"
        # else:
        #     result_dir = result_dir + "/refannot_class/all_annot.xls"
        # prepare option for workflow
        options = {
            "exp_matrix": exp_id + ";" + data.geneset_id + ';' + is_rmbe,
            "corr_main_id": str(main_id),
            'cor_cutoff': data.cor_cutoff,
            'gt': data.level,
            'group_dict': data.group_dict,
            'group_id': data.group_id,
            'padjust_way': data.padjust_way,
            'corr_way': data.corr_way,
            'anno': data.task_id,
            'category': data.task_id + "," + data.level,
            "update_info": json.dumps({str(main_id): "exp_corrsf"})  # to update status
        }
        if "qvalue" in data:
            options.update({'sig_type': 1, 'qvalue_cutoff': data.qvalue})
        if "pvalue" in data:
            options.update({'sig_type': 0, 'pvalue_cutoff': data.pvalue})
        # prepare to file
        to_files = ["whole_transcriptome_v1_1.geneset.export_geneset_exp_matrix(exp_matrix)",
                    "whole_transcriptome.geneset.get_gene_detail_whole(anno)",
                    "whole_transcriptome.geneset.get_gene_type(category)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'whole_transcriptome.report.exp_corrsf'
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
        geneset_info = self.whole_transcriptome.insert_geneset_info(data.geneset_id, "exp_corrsf",
                                                                    str(main_id))
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.whole_transcriptome.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.whole_transcriptome.update_group_compare_is_use(data.task_id, data.control_id)
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
        cmd += "s/whole_transcriptome/exp_corrsf "
        cmd += "-b http://bcl.tsg.com "
        args = dict(
            task_id="tsg_38079",
            task_type="2",
            submit_location="genesetcorrsf",
            # exp_id="5db7c6ca17b2bf4a2ff053a3",
            group_dict=json.dumps({'NFD': ['NFD1', 'NFD2', 'NFD3', 'NFD4'], 'HFD': ['HFD1', 'HFD2', 'HFD3', 'HFD4'], 'NAC_HFD': ['NAC_HFD1', 'NAC_HFD2', 'NAC_HFD3', 'NAC_HFD4']}).replace(
                '"', '\\"'),
            group_id="5f276ea517b2bf2294742998",
            geneset_id="5f2770ff17b2bf2294881920",
            # pvalue="0.05",
            qvalue="0.05",
            cor_cutoff="0.5",
            level="G",
            sig_type="qvalue",
            corr_way="pearson",
            padjust_way="fdr_bh",
            exp_id='5fa3992e17b2bf0df27f4aa0'
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join([str(x) for x in arg_names]), ";".join([str(x) for x in arg_values]))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()
