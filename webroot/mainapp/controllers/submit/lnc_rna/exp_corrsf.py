# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import unittest
import os


class ExpCorrsfAction(LncRnaController):
    def __init__(self):
        super(ExpCorrsfAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type', 'exp_id', 'group_dict', 'group_id',
                      'cor_cutoff', 'corr_way', 'exp_level', 'padjust_way']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument: %s" % arg,
                        'code': 'C2900501', 'variables': variables}
                return json.dumps(info)
        geneset_ids = []
        geneset_ids_dic = {}
        geneset_name = list()
        for arg in ('l_geneset_id', 'm_geneset_id'):
            if hasattr(data, arg):
                obj_id = getattr(data, arg)
                geneset_ids.append((arg, obj_id))
                geneset_ids_dic[arg] = obj_id
                geneset_info = self.lnc_rna.get_main_info(obj_id, 'sg_geneset', data.task_id)
                if not geneset_info:
                    info = {"success": False, "info": "基因集不存在，请确认参数是否正确！!", 'variables': ''}
                    return json.dumps(info)
                geneset_name.append(geneset_info['name'])
        if not geneset_ids:
            info = {'success': False, 'info': "请选择分析的基因集 : %s" % arg,
                    'code': 'C2900501', 'variables': ['l_geneset_id', 'm_geneset_id']}
            return json.dumps(info)

        exp_info = self.lnc_rna.get_exp_params_info(data.exp_id, data.task_id)
        project_sn = exp_info["project_sn"]
        task_id = data.task_id
        print(type(data.group_dict), " data.group_dict")
        if isinstance(data.group_dict, (str, unicode)):
            try:
                group_dict = json.loads(data.group_dict, object_pairs_hook=OrderedDict)
            except ValueError as e:
                return {'success': False, 'info': str(e)}
        elif isinstance(data.group_dict, dict):
            group_dict = data.group_dict
        else:
            return {'success': False,
                    'info': 'group_dict type: dict or string of json ===== ' + str(type(data.group_dict)) + str(
                        data.group_dict)}
        # create main table record
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            exp_id=data.exp_id,
            group_id=data.group_id,
            cor_cutoff=data.cor_cutoff,
            exp_level=data.exp_level,
            padjust_way=data.padjust_way,
            corr_way=data.corr_way,
            group_dict=group_dict,
            sig_type=data.sig_type
        )
        params.update(geneset_ids_dic)
        if "qvalue" in data:
            params.update({'qvalue': data.qvalue})
        if "pvalue" in data:
            params.update({'pvalue': data.pvalue})
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        if len(geneset_name) > 1:
            name = "|".join(geneset_name) + "_ExpCorrsf" + '_' + data.corr_way + '_' + data.padjust_way + '_'
        else:
            name = geneset_name[0] + "_ExpCorrsf" + '_' + data.corr_way + '_' + data.padjust_way + '_'

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
            exp_level=data.exp_level
        )
        # 主表名字为了区分样本相关性分析,这个是做6-15(six-fifteen), 所以叫做sg_exp_corrsf
        main_id = self.lnc_rna.insert_main_table('sg_exp_corrsf', main_info)
        result_dir = self.lnc_rna.get_annotation_stat_info(data.task_id)['result_dir']
        if os.path.exists(result_dir + "/allannot_class/all_annot.xls"):
            result_dir = result_dir + "/allannot_class/all_annot.xls"
        else:
            result_dir = result_dir + "/refannot_class/all_annot.xls"
        # prepare option for workflow
        options = {
            "exp_matrix": data.exp_id + ";" + ';'.join(i for _, i in geneset_ids),
            "genes_info": ';'.join('%s,%s' % (k, v) for k, v in geneset_ids),
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
            options.update({'sig_type': 1, 'qvalue_cutoff': data.qvalue})
        if "pvalue" in data:
            options.update({'sig_type': 0, 'pvalue_cutoff': data.pvalue})
        # prepare to file
        to_files = ["lnc_geneset.export_genesets_exp_matrix(exp_matrix)",
                    "lnc_geneset.export_genes_list(genes_info)"]

        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'lnc_rna.report.exp_corrsf'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(ExpCorrsfAction, self).POST()
        task_info['content'] = {'ids': {'id': str(main_id), 'name': name}}
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.lnc_rna.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.lnc_rna.update_group_compare_is_use(data.task_id, data.control_id)
        return json.dumps(task_info)


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_this(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "s/lnc_rna/exp_corrsf "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                task_id="lnc_rna",
                # Conf配置
                task_type="2",
                # Conf配置
                submit_location="genesetcorrsf",
                exp_id="5cb40efa17b2bf2f4b9c7cb1",
                # json string of dict
                group_dict=json.dumps({"Con": ["Con1", "Con2", "Con3", "Con4", "Con5"],
                                       "Vit": ["Vit1", "Vit2", "Vit3", "Vit4", "Vit5"]}).replace('"', '\\"'),
                group_id="5caae8e917b2bf6ddc645eaf",
                m_geneset_id="5c90894a17b2bf407167bc7f",
                l_geneset_id="5c90890017b2bf407167bc7d",
                pvalue="0.05",
                qvalue="0.05",
                cor_cutoff="0.5",
                exp_level="T",
                sig_type="qvalue",
                corr_way="spearman",
                padjust_way="fdr_bh"
            )
            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join([str(x) for x in arg_names]),
                                             ";".join([str(x) for x in arg_values]))

            print(cmd)
            os.system(cmd)


    unittest.main()
