# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import web
import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.signature import check_sig
from bson import SON
import datetime


class OtusetAction(MetaController):
    """
    otuset 创建基因集
    """
    def __init__(self):
        super(OtusetAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        default_argu = [
            'table_name', 'table_id', 'submit_location', 'otuset_name', 'desc',
            'task_type'
        ]
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {
                    'success': False,
                    'info': 'parameters missing: %s' % argu
                }
                return json.dumps(info)
        task_name = 'meta.report.otuset'
        module_type = 'workflow'
        params_json = {
            'table_name': data.table_name,
            'table_id': data.table_id,
            'submit_location': data.submit_location,
            'task_type': int(data.task_type)
        }
        if hasattr(data, "pvalue"):
            params_json['pvalue'] = float(data.pvalue)
        if hasattr(data, "qvalue"):
            params_json['qvalue'] = float(data.qvalue)
        if hasattr(data, "lda"):
            params_json['lda'] = float(data.lda)
        if hasattr(data, "top"):
            params_json['top'] = int(data.top)
        if hasattr(data, "species_name"):
            params_json['species_name'] = str(data.species_name)
        if hasattr(data, "label"):
            params_json['label'] = str(data.label.replace("&amp;", "&"))
        main_table_name = 'Otuset' + '_' + datetime.datetime.now().strftime(
            "%Y%m%d_%H%M%S%f")[:-3]
        table_info = self.meta.get_main_info(data.table_id, data.table_name)
        task_info = self.meta.get_task_info(table_info['task_id'])
        mongo_data = [('status', 'start'),
                      ('desc', data.desc),
                      ('description', data.desc),
                      ('name', data.otuset_name),
                      ('project_sn', task_info['project_sn']),
                      ('task_id', task_info['task_id']),
                      ('created_ts',
                       datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                      ("params",
                       json.dumps(params_json,
                                  sort_keys=True,
                                  separators=(',', ':')))]
        main_table_id = self.meta.insert_main_table('sg_otuset', mongo_data)
        update_info = {str(main_table_id): 'sg_otuset'}
        options = {
            'table_name': data.table_name,
            'table_id': data.table_id,
            'main_id': str(main_table_id),
            'update_info': json.dumps(update_info),
        }
        if hasattr(data, "pvalue"):
            options['pvalue'] = float(data.pvalue)
        if hasattr(data, "qvalue"):
            options['qvalue'] = float(data.qvalue)
        if hasattr(data, "lda"):
            options['lda'] = float(data.lda)
        if hasattr(data, "top"):
            options['top'] = int(data.top)
        if hasattr(data, "species_name"):
            options['species_name'] = str(data.species_name)
        if hasattr(data, "label"):
            options['label'] = str(data.label.replace("&amp;", "&"))
        to_file = []
        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name="Otuset/" + main_table_name,
                            module_type=module_type, to_file=to_file,
                            main_id=main_table_id,
                            collection_name='sg_otuset')
        task_info = super(OtusetAction, self).POST()
        task_info['content'] = {
            'ids': {
                'id': str(main_table_id),
                'name': main_table_name
            }
        }
        return json.dumps(task_info)


if __name__ == "__main__":
    import os
    test_dict = {
        'otuset_test_venn': {
            'table_name': 'sg_otu_venn',
            'table_id': '5f9632bf17b2bf257aabc492',
            'label': 'C & E',
        },
        'otuset_test_diff': {
            'table_name': 'sg_species_difference_check',
            'table_id': '5f521b7617b2bf6ec16fbaa5',
            'pvalue': 1,
            'qvalue': 1,
        },
        'otuset_test_lefse': {
            'table_name': 'sg_species_difference_lefse',
            'table_id': '5f521b7717b2bf6ec16fbaaa',
            'lda': 1,
            'pvalue': 1,
        },
        'otuset_test_rf': {
            'table_name': 'sg_randomforest',
            'table_id': '5f97897d17b2bf1abe694530',
            'top': 20,
        }
    }
    for test, opts in test_dict.items():
        opt_dict = dict(
            otuset_name=test,
            desc=test,
            submit_location='otuset',
            task_type='2'
        )
        opt_dict.update(opts)
        n = []
        d = []
        [[n.append(str(a)), d.append(str(b))] for a, b in opt_dict.items()]
        n = ';'.join(n)
        d = ';'.join(d)

        api = 's/m/otuset'
        cmd = "python ~/biocluster/bin/webapitest.py post" +\
            " '{}' -n '{}' -d '{}' ".format(api, n, d) +\
            "-c client03 -b http://bcl.tsg.com -fr n"
        os.system(cmd)
