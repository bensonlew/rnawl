# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import web
import json
import datetime
# from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig
# from bson import SON
# from biocluster.file import download, exists


class CorePanCompAction(BacComparativeController):
    def __init__(self, bind_object=None, instant=False):
        super(CorePanCompAction, self).__init__(bind_object, instant)
        self._data = None
        self.wf = 'bac_comp_genome.report.comp'
        self.main_table = 'cor_pan_comp'

    def check_opts(self):
        data = self._data
        opts = ['annotable', 'functiontype', 'level', 'group_id', 'group_detail',
                'corepan', 'pancat']
        for opt in opts:
            if not(hasattr(data, opt)):
                return {'success': False, 'info': '参数缺失: ' + opt}
        group_detail = json.loads(data.group_detail)
        samples = []
        [samples.extend(v) for v in group_detail.values()]
        pan_samples = self.bac_comparative.get_value_of_key('pan', data.corepan, 'spe_list').split(',')
        self.sample_num = len(samples)
        if set(samples) - set(pan_samples):
            return {'success': False, 'info': '选择了泛基因组分析未使用到的样本: {}'.format(set(samples) - set(pan_samples))}

    def set_opts(self):
        data = self._data
        pancat = self.bac_comparative.get_pan_category2(data.pancat, self.sample_num)
        opts = {
            'annotable': data.annotable,  # 1. 根据前段命名修改
            'functiontype': data.functiontype,
            'level': data.level,
            'gfile': data.group_detail,
            'corepan': self.bac_comparative.get_value_of_key('pan', data.corepan, 'cluster_path'),
            'pancat': pancat,
            'result_type': data.result_type,
        }
        return opts

    def set_params(self):
        data = self._data
        params = {
            'functiontype': data.functiontype,
            'group_id': data.group_id,
            'group_detail': json.loads(data.group_detail),
            'submit_location': data.submit_location,
            'task_id': data.task_id,
            'task_type': int(data.task_type),
            #'annotable': data.annotable,
            'corepan': data.corepan,
            'pancat': data.pancat,
            'level': data.level,
            'level_low': data.level_low,
            'result_type': data.result_type,
        }
        return params

    def to_files(self):
        to_file = []
        #to_file.append('bac_comparative.export_group_file(gfile)')
        return to_file

    def info_to_status(self):
        info = {
            'group_id': self._data.group_id,
            'group_lab': self._data.group_detail,
            'pan_group_id': self._data.pancat,
            'pan_category_names': self.bac_comparative.get_value_of_key('pan_group',
                                                                        self._data.pancat,
                                                                        'category_names')
        }
        return info

    @check_sig
    def POST(self):
        self._data = web.input()

        check_info = self.check_opts()
        if check_info:
            return check_info
        options = self.set_opts()
        params = self.set_params()
        params = json.dumps(params, sort_keys=True, separators=(',' ':'))

        task_name = self.wf
        module_type = 'workflow'
        project_sn = self.bac_comparative.get_projectsn(self._data.task_id)

        # 需要根据情况修改
        main_id = self.bac_comparative.insert_none_table(self.main_table)
        update_info = {str(main_id): self.main_table}
        status_info = self.info_to_status()
        options['main_id'] = str(main_id)
        options['main_name'] = self.main_table
        options['update_info'] = json.dumps(update_info)
        options['status_info'] = json.dumps(status_info)
        created_ts = datetime.datetime.now()
        main_table = ''.join(
            map(lambda x: x.capitalize(), self.main_table.split('_'))
        )
        main_table_name = main_table + '_' + created_ts.strftime("%Y%m%d_%H%M%S")
        mongo_data = {
            'project_sn': project_sn,
            'task_id': self._data.task_id,
            'name': main_table_name,
            'params': params,
            'desc': main_table_name + '统计分析',
            'status': 'start',
            'created_ts': created_ts.strftime("%Y-%m-%d %H:%M:%S"),
            'main_id': main_id
        }
        self.bac_comparative.insert_main_table_new(self.main_table,
                                                   main_id,
                                                   mongo_data)

        self.set_sheet_data(name=task_name,
                            options=options,
                            main_table_name=main_table + '/' + main_table_name,
                            module_type=module_type,
                            project_sn=project_sn,
                            task_id=self._data.task_id,
                            params=params,
                            to_file=self.to_files()
                            )

        task_info = super(CorePanCompAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)


if __name__ == '__main__':
    '''test
    '''
    import os
    zz = {"all":["GCF_000014785.1","GCF_000015025.1","GCF_000021245.2","GCF_000162295.1","GCF_000163375.2","GCF_000163395.2","GCF_000173395.1","GCF_003399715.1","GCF_003399765.1"]}
    opt_dict = dict(
        task_id='tsg_37890',
        task_type='2',
        submit_location='core_pan_comp',
        annotable="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/cog.annot.xls",
        functiontype='cog',
        result_type='no',
        level='function',
        level_low='H',
        group_id='5f2784817f8b9aaf218b4587',
        group_detail=json.dumps(zz),
        pancat='5f21066f17b2bf368b7ba947',
        corepan="5f21062917b2bf368b7881cc",
    )
    n = []
    d = []
    [[n.append(a), d.append(b)] for a, b in opt_dict.items()]
    n = ';'.join(n)
    d = ';'.join(d)

    api = 's/bac_comparative/core_pan_comp'
    cmd = "python ~/biocluster/bin/webapitest.py post '{}' -n '{}' -d '{}' -c client03 -b http://bcl.tsg.com -fr n".format(api, n, d)
    os.system(cmd)
