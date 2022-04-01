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
import re


class ErrSiteFiterHomohybridAction(PtController):
    def __init__(self):
        """
        进行错配位点过滤后重新计算的接口-仅仅用于胎儿是纯合与杂合的情况下进行过滤
        modified by hongdong 20191111
        传入的参数必须有该家系的father_id，以及错配点
        """
        super(ErrSiteFiterHomohybridAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'medical.paternity_test_v3.report.err_site_fiter_homohybrid'
        task_type = 'workflow'
        father_result = PT("pt_v3").get_one("sg_father_err_homohybrid", "_id", ObjectId(data.father_err_id))
        if not father_result:
            raise Exception("father_err_id{}在数据库中不存在！".format(data.father_err_id))
        params_json = {
            'father_err_id': data.father_err_id,
            'err_min': father_result['err_min']
        }
        sample_name = []
        old_sample_name = []
        dp_correction = None
        if father_result['types'] == "2":  # 自由交互的结果需要特殊处理下@by hongdong20181027
            samples = PT("pt_v3").get_one("sg_free_combine", "_id", father_result['father_id'])
            if samples:
                params_ = json.loads(samples['params'])
                dp_correction = params_["dp_correction"]
                old_sample_name.extend([params_["case_id"] + "-" + params_['dad_id'], params_['mom_id'],
                                        params_['preg_id']])
                if params_["is_report"] == "true":
                    sample_name.extend([params_['new_dad_id'], params_['new_mom_id'], params_['new_preg_id']])
                else:
                    sample_name.extend(old_sample_name)
            else:
                raise Exception("sg_free_combine中找不到{}".format(father_result['father_id']))
        else:
            temp = str(father_result['name']).split('_')
            file_info = PT("pt_v3").get_one("sg_file_info", "sample_id", temp[0])
            analysis_type = file_info["analysis_type"] if file_info and "analysis_type" in file_info else ""
            if analysis_type == "wqcf" and self.is_ms(father_result['name']):
                sample_name.extend([temp[0], temp[1], "_".join(temp[2:4])])
            else:
                sample_name.extend([temp[0], temp[1], temp[2]])
        print "***********:{}:*********".format(sample_name)
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
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
        main_table_id = PT('pt_v3').add_main_id('sg_father_err_homohybrid')
        update_info = {str(main_table_id): 'sg_father_err_homohybrid'}
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
        if father_result['types'] == "2":
            options.update({
                "old_dad": old_sample_name[0],
                'old_mom': old_sample_name[1],
                'old_son': old_sample_name[2],
                'dp_correction': dp_correction
            })
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params,
                             db_type='pt_v3')
        task_info = super(ErrSiteFiterHomohybridAction, self).POST()
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

    def is_ms(self, sample_id):
        """
        判断是否是ms降噪样本
        :param sample_id:
        :return:
        """
        if re.match('(WQ.*)-S(.*)_M(.*)_error.*', sample_id):
            return True
        else:
            return False
