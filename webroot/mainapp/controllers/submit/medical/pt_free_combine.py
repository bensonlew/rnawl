# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
import web
import json
import datetime
from bson import SON
from mainapp.libs.param_pack import *
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.submit.med_mongo import MedMongo as PT
from mainapp.controllers.project.pt_controller import PtController


class PtFreeCombineAction(PtController):
    """
    该接口用于亲子鉴定的家系自由交互的分析
    laste modified by hongdong@20171208
    """
    def __init__(self):
        super(PtFreeCombineAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'medical.paternity_test_v2.report.pt_free_combine'
        task_type = 'workflow'
        params_json = {
            'case_id': data.case_id,
            'mom_id': data.mom_id,
            'preg_id': data.preg_id,
            'is_report': data.is_report,
            'err_min': data.err_min,
        }
        if hasattr(data, 'new_mom_id'):
            params_json['new_mom_id'] = data.new_mom_id
            params_json['new_dad_id'] = data.new_dad_id
            params_json['new_preg_id'] = data.new_preg_id
        if hasattr(data, "dad_id"):
            params_json['dad_id'] = data.dad_id
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        if str(data.is_report) == "true":
            task_id = "Report_" + data.new_dad_id + '_' + data.new_mom_id + '_' + data.new_preg_id \
                      + '_err_' + data.err_min + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            remark = "{}改为{};{}改为{};{}改为" \
                     "{}".format(data.case_id + '-' + data.dad_id, data.new_dad_id, data.mom_id,
                                 data.new_mom_id, data.preg_id, data.new_preg_id)
            dad_id = "{}-{}".format(data.case_id, data.dad_id)
            mom_id = data.mom_id
            preg_id = data.preg_id
            report_name = '_'.join([data.new_dad_id, data.new_mom_id])
        else:
            if not data.dad_id:
                father_name = "{}--".format(data.case_id)
            else:
                father_name = "{}-{}".format(data.case_id, data.dad_id)
            task_id = "Dedup_" + father_name + '_' + data.mom_id + '_' + data.preg_id + '_err_' + data.err_min + "_" +\
                      datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            remark = ""
            dad_id = father_name
            mom_id = data.mom_id
            preg_id = data.preg_id
            report_name = '--'
        mongo_data = [
            ('params', params),
            ('name', task_id),
            ('desc', '亲子鉴定家系自由组合正在分析！'),
            ('member_id', data.member_id),
            ('status', 'start'),
            ('dad_id', dad_id),
            ('mom_id', mom_id),
            ('preg_id', preg_id),
            ('report_name', report_name),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('remark', remark)
        ]
        main_table_id = PT('pt_v2').add_main_id('sg_free_combine')
        update_info = {str(main_table_id): 'sg_free_combine'}
        options = {
            'case_id': data.case_id,
            'mom_id': data.mom_id,
            'preg_id': data.preg_id,
            'err_min': data.err_min,
            'is_report': data.is_report,
            "update_info": json.dumps(update_info),
            'father_id': str(main_table_id),
            'main_table_data': SON(mongo_data),
            'member_id': data.member_id
        }
        if hasattr(data, 'dad_id'):
            options['dad_id'] = data.dad_id
        if hasattr(data, 'new_mom_id'):
            options['new_mom_id'] = data.new_mom_id
            options['new_dad_id'] = data.new_dad_id
            options['new_preg_id'] = data.new_preg_id
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params, db_type='pt_v2')
        task_info = super(PtFreeCombineAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': task_id}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传入的参数,必须有母本id\胎儿id\父本组id\允许错配数\member_id\是否出报告
        :param data:网页端传入的数据(是否全库查重的参数必须传递)
        :return: 检查结果
        """
        params_name = ['mom_id', 'case_id', 'preg_id', 'err_min', 'dad_id', 'member_id', 'is_report']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数{}".format(names))
        if str(data.is_report) == "true":
            if not data.new_dad_id or not data.new_mom_id or not data.new_preg_id:
                success.append("必须同时输入新的父本，母本，胎儿的编号，否则后面无法进行家系分析！")
            if not data.dad_id:
                success.append("缺少父本过滤器参数，必须输入完整的父本组，过滤器，母本，胎儿数据，否则后面无法进行家系分析！")
        return success
