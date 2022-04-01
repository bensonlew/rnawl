# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import web
import json
import datetime
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController
from mainapp.libs.signature import check_sig


class CogDiffAction(BacComparativeController):
    def __init__(self, bind_object=None, instant=False):
        super(CogDiffAction, self).__init__(bind_object, instant)
        self._data = None
        self.main_table = 'cog_diff'
        self.wf = 'bac_comp_genome.report.diff'
        self.functiontype = 'COG'

    def check_opts(self):
        data = self._data
        print data
        opts = ['type', 'annotable', 'level', 'test', 'norm']
        for opt in opts:
            if not(hasattr(data, opt)):
                return {'success': False, 'info': '参数缺失:%s' % opt}
        if data.type == 'group':
            if not(hasattr(data, 'group_detail')):
                 return {'success': False, 'info': '组间比较缺分组文件'}
        elif data.type == 'sample':
            if not(hasattr(data, 'samples')):
                return {'success': False, 'info': '样本间比较确定两个样本名'}

    def set_opts(self):
        '''接口参数
            annotable 注释表的直接给文件路径
            gfile 分组文件路径
        '''
        data = self._data    
        opts = {
            'annotable': data.annotable,  # 1. 根据前端命名修改
            'functiontype': self.functiontype,
            'level': data.level,
            'test': data.test,
            'norm': data.norm,
        }
        if data.type == 'sample':
            [opts['samp1'], opts['samp2']] = data.samples.split(',')
        elif data.type == 'group':
            opts['gfile'] = data.group_detail
        return opts

    def set_params(self):
        '''params 设置

        '''
        data = self._data
        opts = {}
        params = dict(
            submit_location=data.submit_location,
            task_id=data.task_id,
            task_type=int(data.task_type),
            #annotable=data.annotable,
            level=data.level,
            test=data.test,
            norm=data.norm,
            type=data.type,
            level_low=data.level_low
        )
        if data.type == 'group':
            params['group_detail'] = json.loads(data.group_detail)
            params['group_id'] = data.group_id
        else:
            params['samples'] = data.samples

        for k, v in params.items():
            opts[k] = str(v) if isinstance(v, (int, float)) else v
        opts['task_type'] = int(opts['task_type'])
        return opts

    def to_files(self):
        to_file = []
        #to_file.append('bac_comparative.export_group_file(gfile)')
        return to_file

    def info_to_status(self):
        if self._data.type == 'group':
            info = dict(
                group_id = self._data.group_id,
                group_lab = ','.join(sorted(json.loads(self._data.group_detail)))
            )
        else:
            info = dict(
                samples=','.join(sorted(self._data.samples.split(',')))
            )
        return info

    @check_sig
    def POST(self):
        self._data = web.input()

        check_info = self.check_opts()
        if check_info:
            return check_info
        options = self.set_opts()
        params = self.set_params()
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))

        task_name = self.wf
        module_type = 'workflow'
        project_sn = self.bac_comparative.get_projectsn(self._data.task_id)

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
            'desc': main_table + '分析',
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

        task_info = super(CogDiffAction, self).POST()
        if task_info['success']:
            task_info['content'] = {
                'ids': {'id': str(main_id), 'name': main_table_name}
                }
        return json.dumps(task_info)


if __name__ == '__main__':
    '''
    test
    '''
    import os
    def sg_task_info():
        '''
        测试用，导入 sg_taks 信息, 第一次测试需要导表
        '''
        from bson import SON
        from biocluster.config import Config
        # task_id：tsg_36106；project_sn：188_5dc9218341905
        info = {
            "task_id": "tsg_36106",
            "member_id": "-",
            "project_sn": "188_5dc9218341905",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d%H:%M:%S"),
            "member_type": 1,
            "cmd_id": 123,
            "is_demo": 0,
            "demo_id": "-"
        }
        cl = Config().get_mongo_client(mtype='bac_comparative')
        db = cl[Config().get_mongo_dbname('bac_comparative')]
        if not db['sg_task'].find_one({'task_id': 'tsg_36106'}):
            db['sg_task'].insert_one(SON(info))

    opt_dict = dict(
        task_id='tsg_36106',
        #group_file="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/group2.txt",
        gfile='5dc53932e4efe6fcf0548318',
        groups='ALL',
        annotable="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/cog.annot.xls",
        test='anova',
        testtype='two.sided',
        level='Function',
        submit_location='cog_diff',
        type='group',
        norm='T'
    )

    n = []
    d = []
    [[n.append(a), d.append(b)] for a, b in opt_dict.items()]
    n = ';'.join(n)
    d = ';'.join(d)

    api = 's/bac_comparative/cog_diff'
    cmd = "python ~/biocluster/bin/webapitest.py post '{}' -n '{}' -d '{}' -c client03 -b http://192.168.12.101:9090 -fr n".format(
        api, n, d)
    # exit(cmd)
    #sg_task_info()
    os.system(cmd)

