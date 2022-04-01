# -*- coding: utf-8 -*-
# 2019-03-21
import web
import json
import datetime
from mainapp.controllers.project.lnc_rna_controller import LncRnaController
from mainapp.libs.signature import check_sig
import os
import unittest


class EnrichHeatAction(LncRnaController):
    def __init__(self):
        super(EnrichHeatAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()

        # 参数检测
        wanted_args = ['task_id', 'submit_location', 'geneset_type', 'enrich_type',
                       'geneset_enrich_ids', 'task_type', 'pvalue_type']

        sheet_dict = {}
        params_json = {}

        for arg in wanted_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '%s参数缺少!' % arg, 'code': 'C2901601', 'variables': [arg]}
                return json.dumps(info)
            value = getattr(data, arg)
            if arg in ('task_type', 'top_num'):
                value = int(value)
            sheet_dict[arg] = value
            params_json[arg] = value

        if data.enrich_type.lower() == 'kegg':
            data_table_name = 'sg_geneset_kegg_enrich'
            collection_name = 'sg_geneset_kegg_enrich_heatmap'
        else:
            collection_name = 'sg_geneset_go_enrich_heatmap'
            if not hasattr(data, 'go_level'):
                info = {'success': False, 'info': 'go_level 参数缺少!', 'code': 'C2901601', 'variables': ['go_level']}
                return json.dumps(info)
            go_level = getattr(data, 'go_level')
            if go_level == 'ALL':
                go_level = go_level.lower()
            params_json['go_level'] = go_level
            sheet_dict['go_level'] = go_level
            data_table_name = 'sg_geneset_go_enrich'

        top_num = ''
        if hasattr(data, 'ids_list'):
            ids_list = getattr(data, 'ids_list')
            params_json['ids_list'] = ids_list
            sheet_dict['ids_list'] = ids_list.strip().split(',')
        else:
            for name in ('p_value', 'top_num'):
                if not hasattr(data, name):
                    info = {'success': False, 'info': name + ' 参数缺少!', 'code': 'C2901601', 'variables': ['p_value']}
                    return json.dumps(info)
            params_json['p_value'] = data.p_value
            sheet_dict['p_value'] = float(data.p_value)
            top_num = int(getattr(data, 'top_num'))
            params_json['top_num'] = top_num
            sheet_dict['top_num'] = top_num

        # 判断传入的基因集id是否存在
        for geneset in data.geneset_enrich_ids.split(","):
            geneset_info = self.lnc_rna.get_main_info(geneset, data_table_name, data.task_id)
            if not geneset_info:
                info = {"success": False, "info": "geneset enrich id %s 不存在，请确认参数是否正确！!" % geneset,
                        'code': 'C2901602', 'variables': ''}
                return json.dumps(info)

        # 两个geneset_id到sg_task返回的project_sn等信息是一样的
        task_info = self.lnc_rna.get_task_info(data.task_id)

        datetime_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        main_table_name = data.enrich_type.lower() + '_heat_map_' + datetime_str[:-3]

        mongo_data = [
            ('project_sn', task_info['project_sn']),
            ('task_id', task_info['task_id']),
            ('status', 'start'),
            ('type', data.enrich_type),
            ('pvalue_type', data.pvalue_type),
            ('name', main_table_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("params", json.dumps(params_json, sort_keys=True, separators=(',', ':'))),
        ]
        if top_num:
            mongo_data.append(('top_num', top_num))
        if 'go_level' in data:
            mongo_data.append(('go_level', data.go_level))

        main_table_id = self.lnc_rna.insert_main_table(collection_name, mongo_data)
        update_info = {str(main_table_id): collection_name}
        # geneset_type只有2种，T或者G,s而页面限制只能同时取G或者T，所以这里用一个geneset_id作为示例即可
        options = {
            "task_id": data.task_id,
            'update_info': json.dumps(update_info),
            "main_table_id": str(main_table_id),
            "data": sheet_dict,
            "geneset_type": data.geneset_type,
        }

        self.set_sheet_data(
            name='lnc_rna.report.enrich_heat',
            module_type="workflow",
            project_sn=task_info['project_sn'],
            task_id=task_info['task_id'],
            main_table_name=main_table_name,
            to_file='lnc_geneset.export_enrich_heatmap(data)',
            options=options)

        task_info = super(EnrichHeatAction, self).POST()

        task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test_go(self):
            cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
            cmd += 'post '
            cmd += "-fr no "
            cmd += '-c {} '.format("client03")
            cmd += "i/lnc_rna/enrich_heat "
            cmd += "-b http://192.168.12.101:9090 "
            args = dict(
                task_id="lnc_rna",
                task_type="1",
                submit_location="enrichgoheat",

                geneset_type="T",
                enrich_type="go",
                geneset_enrich_ids="5caaff6017b2bf2cef6be03b,5caaffdc17b2bf5046947a1c",
                go_level='all',  # "BP", "CC", "MF"
                # ids_list='GO:2000288,GO:1902894,GO:0006287,GO:0000183,GO:0009299,GO:0060149,GO:0060967,'
                #          'GO:0048762,GO:0048663,GO:0006270,GO:0006284,GO:0010092,GO:0000722,GO:0060964,'
                #          'GO:0006383,GO:0000730,GO:2000779,GO:0060850,GO:0045005,GO:0045599,GO:0042276,'
                #          'GO:0031297,GO:0051573,GO:0031057,GO:0000729,GO:0045814,GO:0036297,GO:0008406,'
                #          'GO:0048863,GO:0006304,GO:0019827,GO:0060147,GO:0060966,GO:0000723,GO:0032508',
                pvalue_type='pvalue',
                p_value='0.05',
                top_num='20'
            )

            arg_names, arg_values = args.keys(), args.values()
            cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
            print(cmd)
            os.system(cmd)
        #
        # def test_kegg(self):
        #     cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        #     cmd += 'post '
        #     cmd += "-fr no "
        #     cmd += '-c {} '.format("client03")
        #     cmd += "i/lnc_rna/enrich_heat "
        #     cmd += "-b http://192.168.12.101:9090 "
        #     args = dict(
        #         task_id="lnc_rna",
        #         task_type="1",
        #         submit_location="enrichkeggheat",
        #
        #         geneset_type="T",
        #         enrich_type="kegg",
        #         geneset_enrich_ids="5caaff6017b2bf2cef6be03b,5caaffdc17b2bf5046947a1c",
        #         # ids_list='',
        #         pvalue_type='pvalue',
        #         p_value='0.05',
        #         top_num='20'
        #     )
        #
        #     arg_names, arg_values = args.keys(), args.values()
        #     cmd += '-n "{}" -d "{}" '.format(";".join(arg_names), ";".join(arg_values))
        #     print(cmd)
        #     os.system(cmd)


    unittest.main()
