# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import web
import json
import datetime
from bson import ObjectId
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig


class FunctionEnrichmentAction(BacComparativeController):
    def __init__(self, bind_object=None, instant=False):
        super(FunctionEnrichmentAction, self).__init__(bind_object)
        self._data = None
        self.main_table = 'kegg_enrichment'
        self.wf = 'bac_comp_genome.report.corepan_enrichment'

    def check_opts(self):
        data = self._data
        opts = ['corepan', 'pancat', 'annotable', 'corrected']
        for opt in opts:
            if not hasattr(data, opt):
                return {'success': False, 'info': '参数缺失:' + opt}


    def set_opts(self):
        data = self._data
        sample_num = len(self.bac_comparative.get_value_of_key('pan', data.corepan, 'spe_list').split(','))
        pancat = self.bac_comparative.get_pan_category2(data.pancat, sample_num)
        opts = {
            'corepan': self.bac_comparative.get_value_of_key('pan', data.corepan, 'cluster_path'),
            'pancat': pancat,
            'samples': data.samples,
            'annotable': data.annotable,
            'corrected': data.corrected,
            'main_name': self.main_table,
        }
        return opts

    def set_params(self):
        data = self._data
        params = dict(
            submit_location=data.submit_location,
            task_id=data.task_id,
            task_type=int(data.task_type),
            corepan=data.corepan,
            corrected=str(data.corrected),
            pancat=data.pancat,
            samples=data.samples,
            #annotable=data.annotable,
        )
        params = json.dumps(params, sort_keys=True, separators=(',' ':'))
        return params

    def info_to_status(self):
        info = dict(
            pan_group_id=ObjectId(self._data.pancat),
            catgory_names=self.bac_comparative.get_value_of_key('pan_group', pan_id, 'category_names'),
        )
        return info

    def to_files(self):
        to_file = []
        return to_file

    @check_sig
    def POST(self):
        self._data = web.input()
        task_name = self.wf
        module_type = 'workflow'
        project_sn = self.bac_comparative.get_projectsn(self._data.task_id)
        check_info = self.check_opts()
        if check_info:
            return check_info
        options = self.set_opts()
        params = self.set_params()

        # 导主表，返回main_id
        main_id = self.bac_comparative.insert_none_table(self.main_table)
        update_info = {str(main_id): self.main_table}
        options['main_id'] = str(main_id)
        options['update_info'] = json.dumps(update_info)
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
            'desc': main_table + '富集分析',
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

        task_info = super(FunctionEnrichmentAction, self).POST()
        #new_info = self.info_to_status()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_id), 'name': main_table_name}}
        return json.dumps(task_info)


if __name__ == '__main__':
    import os
    # sg_task_info()
    opt_dict = dict(
        task_id = 'tsg_36106',
        #pancat='{"core": [100, 100], "dis": [1, 99], "uniq":[1, 1]}',
        pancat='5dca3d6c17b2bf7352f9d088',
        corepan="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/clusters.txt",
        samples='GCF_000007685.1,GCF_000092565.1',
        annotable = "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/kegg.annot.xls",
        corrected = 'fdr',
        submit_location='function_enrichemnt'
    )
    n = []
    d = []
    [[n.append(a), d.append(b)] for a, b in opt_dict.items()]
    n = ';'.join(n)
    d = ';'.join(d)

    api = 's/bac_comparative/function_enrichment'
    cmd = "python ~/biocluster/bin/webapitest.py post '{}' -n '{}' -d '{}' -c client03 -b http://bcl.tsg.com -fr n".format(api, n, d)
    #exit(cmd)
    os.system(cmd)

    def sg_task_info():
        '''
        测试用，导入 sg_taks 信息, 第一次测试运行
        '''
        from bson import SON
        from biocluster.config import Config
        info = {
            "task_id" : "tsg_36106",
            "member_id" : "-",
            "project_sn" : "188_5dc9218341905",
            "created_ts" : datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "member_type" : 1,
            "cmd_id" : 123,
            "is_demo" : 0,
            "demo_id" : "-"
        }
        cl = Config().get_mongo_client(mtype='bac_comparative')
        db = cl[Config().get_mongo_dbname('bac_comparative')]
        db['sg_task'].insert_one(SON(info))

