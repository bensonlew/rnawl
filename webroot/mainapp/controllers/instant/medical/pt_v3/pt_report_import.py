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


class PtReportImportAction(PtController):
    def __init__(self):
        """
        进行报告存入的时候用于整合好数据汇总表
        modified by hongdong 20171207
        """
        super(PtReportImportAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'medical.paternity_test_v3.report.report_import'
        task_type = 'workflow'
        postype = 'hybrid'
        if hasattr(data, 'datatype'):
            if data.datatype not in ['hybrid', 'homozygousandhybrid']:
                info = {"success": False, "info": 'datatype值不正确！'}
                return json.dumps(info)
            else:
                if data.datatype == 'hybrid':
                    postype = 'hybrid'
                elif data.datatype == 'homozygousandhybrid':
                    postype = 'homohybrid'
        if postype == 'hybrid':
            father_result = PT("pt_v3").get_one("sg_father_err", "_id", ObjectId(data.father_err_id))
        else:
            father_result = PT("pt_v3").get_one("sg_father_err_homohybrid", "_id", ObjectId(data.father_err_id))
        if not father_result:
            raise Exception("father_err_id{}在数据库中不存在！-{}".format(data.father_err_id, postype))
        params_json = {
            'father_err_id': data.father_err_id,
            'report_time': data.report_time
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        sample_name = str(father_result['name']).split('_')
        name = "Report_import_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = [
            ('member_id', data.member_id),
            ('family_id', father_result['family_id']),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('params', params),
            ('name', name),
            ('desc', '开始进行报告数据导入！'),
            ('status', 'start'),
            ('father_err_id', data.father_err_id)
        ]
        main_table_id = PT('pt_v3').add_main_id('sg_report')
        update_info = {str(main_table_id): 'sg_report'}
        file_info = PT("pt_v3").get_one("sg_file_info", "sample_id", sample_name[0])
        analysis_type = file_info["analysis_type"] if file_info and "analysis_type" in file_info else ""
        options = {
            'mom_id': sample_name[1],
            'son_id': "_".join(sample_name[2:4]) if analysis_type == "wqcf" else sample_name[2],
            'err_min': father_result['err_min'],
            "update_info": json.dumps(update_info),
            'dad_id': sample_name[0],
            'main_table_data': SON(mongo_data),
            'report_time': data.report_time,
            'member_id': data.member_id,
            'father_err_id': str(data.father_err_id),
            'report_id': str(main_table_id),
            'postype': postype
        }
        print options
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params,
                             db_type='pt_v3')
        task_info = super(PtReportImportAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传入的参数
        :param data:网页端传入的数据(是否全库查重的参数必须传递)
        :return: 检查结果
        """
        success = []
        params_name = ['father_err_id', 'member_id', 'report_time']
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数{}".format(names))
        return success
