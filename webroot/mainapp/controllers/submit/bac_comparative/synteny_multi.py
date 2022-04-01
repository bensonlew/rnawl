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


class SyntenyMultiAction(BacComparativeController):
    def __init__(self, bind_object=None, instant=False):
        super(SyntenyMultiAction, self).__init__(bind_object, instant)
        self._data = None
        self.wf = 'bac_comp_genome.report.mummer'
        self.main_table = 'synteny_multi'

    def check_opts(self):
        data = self._data
        opts = ['mummer', 'samples', 'seq_dir']
        for opt in opts:
            if not(hasattr(data, opt)):
                return {'success': False, 'info': '参数缺失: ' + opt}

    def set_opts(self):
        data = self._data
        opts = {
            'seq_dir': data.seq_dir,  # 1. 根据前段命名修改
            'samples': data.samples,
            'mummer': data.mummer,
        }
        if (hasattr(data, 'ref')):
            opts['ref'] = data.ref
        return opts

    def set_params(self):
        data = self._data
        sp = data.samples.split(',')
        params = {
            'submit_location': data.submit_location,
            'task_id': data.task_id,
            'task_type': int(data.task_type),
            'mummer': data.mummer,
            'samples': sorted(sp)
        }
        return params

    def to_files(self):
        to_file = []
        return to_file

    def info_to_status(self):
        info = {}
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
            'desc': main_table_name,
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

        task_info = super(SyntenyMultiAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)


if __name__ == '__main__':
    '''test
    '''
    import os
    opt_dict = {
                    "ref": "GCF_000007685.1",
                    'mummer': 'nucmer',
                    'seq_dir': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/packages/t/seq/',
                    'samples': 'GCF_000007685.1,GCF_000216235.1,GCF_000216095.1',
                    'task_id': 'tsg_36106',
                    'submit_location': 'synteny_multi'
                    }
    n = []
    d = []
    [[n.append(a), d.append(b)] for a, b in opt_dict.items()]
    n = ';'.join(n)
    d = ';'.join(d)

    api = 's/bac_comparative/synteny_multi'
    cmd = "python ~/biocluster/bin/webapitest.py post '{}' -n '{}' -d '{}' -c client03 -b http://192.168.12.101:9090 -fr n".format(api, n, d)
    os.system(cmd)
