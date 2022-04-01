# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from mainapp.models.mongo.submit.med_mongo import MedMongo as PT
from mainapp.controllers.project.pt_controller import PtController
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from bson import ObjectId
from bson import SON
import datetime
import json
import web


class PtErrSiteFiterAction(PtController):
    def __init__(self):
        """
        进行错配位点过滤后重新计算的接口
        modified by hongdong 20171205
        传入的参数必须有该家系的father_id，以及错配点
        """
        super(PtErrSiteFiterAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'medical.paternity_test_v2.report.err_site_fiter'
        task_type = 'workflow'
        father_result = PT("pt_v2").get_one("sg_father_err", "_id", ObjectId(data.father_err_id))
        if not father_result:
            raise Exception("father_err_id{}在数据库中不存在！".format(data.father_err_id))
        params_json = {
            'father_err_id': data.father_err_id,
            'err_min': father_result['err_min']
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        sample_name = str(father_result['name']).split('_')
        mongo_data = [
            ('member_id', data.member_id),
            ('batch_id', father_result['batch_id']),
            ('family_id', father_result['family_id']),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('params', params),
            ('name', father_result['name']),
            ('desc', '开始进行位点过滤分析！'),
            ('status', 'start'),
            ('err_min', father_result['err_min']),
            ('father_id', father_result['father_id']),
            ('types', father_result['types']),
            ('is_filter', 'yes')
        ]
        main_table_id = PT('pt_v2').add_main_id('sg_father_err')
        update_info = {str(main_table_id): 'sg_father_err'}
        options = {
            'mom_id': sample_name[1],
            'son_id': sample_name[2],
            'err_min': father_result['err_min'],
            "update_info": json.dumps(update_info),
            'dad_id': sample_name[0],
            'main_table_data': SON(mongo_data),
            'father_id': str(father_result['father_id']),
            'member_id': data.member_id,
            'father_err_id_old': str(data.father_err_id)
        }
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params,
                             db_type='pt_v2')
        task_info = super(PtErrSiteFiterAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': father_result['name']}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传入的参数
        :param data:网页端传入的数据(是否全库查重的参数必须传递)
        :return: 检查结果
        """
        success = []
        params_name = ['father_err_id', 'member_id', 'client']
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数{}".format(names))
        if data.client not in ['client01', 'client03']:
            success.append('{}不在client合法范围内！'.format(data.client))
        return success
