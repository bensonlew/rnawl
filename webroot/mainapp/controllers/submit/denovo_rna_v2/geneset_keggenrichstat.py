# -*- coding: utf-8 -*-
import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.denovo_rna_v2_controller import DenovoRnaV2Controller
from mbio.api.to_file.denovo_rna_v2 import *
from mainapp.libs.signature import check_sig
import os
import unittest


class GenesetKeggenrichstatAction(DenovoRnaV2Controller):
    def __init__(self):
        super(GenesetKeggenrichstatAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        basic_args = ["task_id", "submit_location", 'task_type']
        basic_args += ['kegg_enrich_id']
        # check arg
        for arg in basic_args:
            if not hasattr(data, arg):
                variables = []
                variables.append(arg)
                info = {'success': False, 'info': "Lack argument:   %s" % arg, 'code': 'C1602701', 'variables': variables}
                return json.dumps(info)
        if data.draw_for_ids == "yes":
            if not data.pathway_ids:
                variables = []
                lag_info="pathway_ids"
                variables.append(lag_info)
                info = {'success': False, 'info': "Lack argument:   %s" % lag_info, 'code': 'C1602702','variables': variables}
                return json.dumps(info)
        else:
            other_args = ["stat_level","stat_threshold_value","stat_numbers_value"]
            for arg in other_args:
                if not hasattr(data,arg):
                    variables = []
                    variables.append(arg)
                    info = {'success': False, 'info': "Lack argument:   %s" % arg, 'code': 'C1602703','variables': variables}
                    return json.dumps(info)

        enrich_info = self.denovo_rna_v2.get_genesetkeggenrich_params_info(data.kegg_enrich_id, data.task_id)
        project_sn = enrich_info["project_sn"]
        task_id = data.task_id
        kegg_enrich_name=enrich_info["name"]
        kegg_enrich_params=json.loads(enrich_info["params"])
        geneset_id=kegg_enrich_params["geneset_id"]
        geneset_info = self.denovo_rna_v2.get_main_info(geneset_id, 'sg_geneset', data.task_id)
        #geneset_name=geneset_info["name"]

        # create main table record
        enrich_genesettype = geneset_info['type']
        version="v3"
        params = dict(
            task_id=task_id,
            submit_location=data.submit_location,
            task_type=int(data.task_type),
            kegg_enrich_id=data.kegg_enrich_id,
            geneset_type=enrich_genesettype,
            kegg_enrich_name=kegg_enrich_name,
            #geneset_name=geneset_name
            draw_for_ids=data.draw_for_ids,
            pathway_ids=data.pathway_ids,
            stat_level=data.stat_level,
            stat_threshold_value=data.stat_threshold_value,
            stat_numbers_value=data.stat_numbers_value
        )
        # for special args
        # if data.draw_for_ids == "yes":
        #     #draw_type = "draw_for_ids"
        #     params.update(dict(draw_for_ids=data.draw_for_ids))
        #     params.update(dict(pathways_ids=data.pathway_ids))
        # else:
        #     #draw_type = "draw_for_filtrate"
        #     #params.update(dict(draw_type=draw_type))
        #     params.update(dict(pathways_ids=data.pathway_ids))
        #     params.update(dict(stat_level=data.stat_level))
        #     params.update(dict(stat_threshold_value=data.stat_threshold_value))
        #     params.update(dict(stat_numbers_value=data.stat_numbers_value))

        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        name = "Genesetkeggstat" + '_'
        time_now = datetime.datetime.now()
        name += time_now.strftime("%Y%m%d_%H%M%S")
        if type(params) == dict:
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        main_info = dict(
            project_sn=project_sn,
            task_id=task_id,
            name=name,
            version=version,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            kegg_enrich_id=data.kegg_enrich_id,
            desc='geneset kegg enrich stat main table',
            params=params,
            status="start",
        )
        if data.draw_for_ids == "yes":
            main_info.update(dict(draw_for_ids=data.draw_for_ids))
        main_id = self.denovo_rna_v2.insert_main_table('sg_geneset_kegg_enrich_stat', main_info)

        # prepare option for workflow
        if data.draw_for_ids == "yes":
            options = {
                "geneset_kegg_enrich_info": data.kegg_enrich_id,
                "pathway_ids":data.pathway_ids,
                "geneset_kegg_enrich_stat_id": str(main_id),
                "update_info": json.dumps({str(main_id): "sg_geneset_kegg_enrich_stat"}),
                "task_id": data.task_id,
                #"geneset_name":geneset_name,
                "geneset_list":geneset_id
        }
            to_files = ["denovo_rna_v2.export_kegg_enrich_info_ids(geneset_kegg_enrich_info)"]
        else:
            options = {
                "geneset_kegg_enrich_info": data.kegg_enrich_id,
                "geneset_kegg_enrich_stat_id": str(main_id),
                "update_info": json.dumps({str(main_id): "sg_geneset_kegg_enrich_stat"}),
                "task_id": data.task_id,
                "stat_level":data.stat_level,
                "stat_threshold_value":data.stat_threshold_value,
                "stat_numbers_value":data.stat_numbers_value,
                #"geneset_name": geneset_name,
                "geneset_list":geneset_id
            }
            to_files = ["denovo_rna_v2.export_kegg_enrich_info_filter(geneset_kegg_enrich_info)",
                        'denovo_rna_v2.export_gene_list(geneset_list)']
        # prepare to file


        # 把参数交给workflow运行相应的tool， 其中to_file用于准备tool的输入文件
        task_name = 'denovo_rna_v2.report.geneset_keggenrichstat'
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=name,  # 设置交互分析结果目录名
                            module_type="workflow",
                            to_file=to_files,
                            project_sn=project_sn,
                            task_id=data.task_id)

        # 运行workflow 并传回参数
        task_info = super(GenesetKeggenrichstatAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_id),
                'name': name
                }
        }
        if 'group_id' in data and str(data.group_id).lower() != 'all':
            _ = self.denovo_rna_v2.update_group_is_use(data.task_id, data.group_id)
        if 'control_id' in data:
            _ = self.denovo_rna_v2.update_group_compare_is_use(data.task_id, data.control_id)
        # task_info['group_dict'] = group_dict
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
        cmd += "s/denovo_rna_v2/geneset_keggenrichstat "
        cmd += "-b http://bcl.tsg.com  "
        args = dict(
            task_id="tsg_33857",
            task_type="2",
            submit_location="genesetkegg_rich_stat",
            kegg_enrich_id="5d37bf7317b2bf3dab249669",
            stat_threshold_value="0.05",
            stat_numbers_value="5",
            stat_level="pvalue",
            draw_for_ids="",
            pathway_ids=""
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestFunction('test_this'))
    unittest.TextTestRunner(verbosity=2).run(suite)
