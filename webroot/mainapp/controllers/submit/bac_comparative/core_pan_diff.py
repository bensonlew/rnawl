# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from cog_diff import CogDiffAction
import json


class CorePanDiffAction(CogDiffAction):
    def __init__(self, bind_object=None):
        super(CorePanDiffAction, self).__init__(bind_object)
        self.main_table = 'core_pan_diff'
        self.functiontype = None

    def set_opts(self):
        '''
        接口参数设置 20200508
        '''
        data = self._data
        samples = []
        if data.type == 'pan':
            samples = data.splist.split(',')
        else:
            group_detail = json.loads(data.group_detail)
            [samples.extend(v) for v in group_detail.values()]
        pan_samples = self.bac_comparative.get_value_of_key('pan', data.corepan, 'spe_list').split(',')
        sample_num = len(samples)
        if set(samples) - set(pan_samples):
            return {'success': False, 'info': '选择了泛基因组分析未使用到的样本: {}'.format(set(samples) - set(pan_samples))}
        print data.pancat
        pancat = self.bac_comparative.get_pan_category2(data.pancat, sample_num)
        print("pancat is ")
        print pancat
        new_opts = {
            'corepan': self.bac_comparative.get_value_of_key('pan', data.corepan, 'cluster_path'),
            'pancat': self.bac_comparative.get_pan_category2(data.pancat, sample_num),
            'selectedcat': data.selectedcat,
            'functiontype': data.functiontype
        }
        if data.type == 'pan':
            self._data.samples = self._data.selectedcat
            new_opts['splist'] = data.splist
        opts = super(CorePanDiffAction, self).set_opts()
        opts.update(new_opts)
        return opts

    def set_params(self):
        data = self._data
        params = super(CorePanDiffAction, self).set_params()
        params['corepan'] = data.corepan
        params['selectedcat'] = data.selectedcat
        params['pancat'] = data.pancat
        params['functiontype'] = data.functiontype

        return params

    def info_to_status(self):
        info = super(CorePanDiffAction, self).info_to_status()
        info['pan_group_id'] = self._data.pancat
        info['pan_category_names'] = self.bac_comparative.get_value_of_key('pan_group', self._data.pancat, 'category_names')
        return info


if __name__ == '__main__':
    '''
    test
    '''
    import os

    opt_dict = dict(
        task_id='tsg_36106',
        #pancat='{"core": [100, 100], "dis": [1, 99], "uniq":[1, 1]}',
        pancat='5dca3d6b17b2bf7352f9d087',
        corepan="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/clusters.txt",
        gfile='5dc53932e4efe6fcf0548318',
        groups='ALL',
        annotable="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/cog.annot.xls",
        test='anova',
        testtype='two.sided',
        level='Function',
        #selectedcat='core',
        selectedcat='core,dispensable',
        #type='group',
        splist='ALL',
        type='pan',
        functiontype='COG',
        norm='T',
        submit_location='core_pan_diff'
    )

    n = []
    d = []
    [[n.append(a), d.append(b)] for a, b in opt_dict.items()]
    n = ';'.join(n)
    d = ';'.join(d)

    api = 's/bac_comparative/core_pan_diff'
    cmd = "python ~/biocluster/bin/webapitest.py post '{}' -n '{}' -d '{}' -c client03 -b http://bcl.tsg.com -fr n".format(
        api, n, d)
   # sg_task_info()
    print cmd
    os.system(cmd)

    def sg_task_info():
        '''
        测试用，导入 sg_taks 信息, 第一次测试需要导表
        '''
        from bson import SON
        from biocluster.config import Config
        info = {
            "task_id": "core_pantest",
            "member_id": "-",
            "project_sn": "188_5b5acb3018914",
            "created_ts": datetime.datetime.now().strftime("% Y-%m-%d %H: % M: % S"),
            "member_type": 1,
            "cmd_id": 123,
            "is_demo": 0,
            "demo_id": "-"
        }
        cl = Config().get_mongo_client(mtype='bac_comparative')
        db = cl[Config().get_mongo_dbname('bac_comparative')]
        if not d['sg_task'].find_one({'task_id': 'core_pan_diff_test'}):
            db['sg_task'].insert_one(SON(info))
