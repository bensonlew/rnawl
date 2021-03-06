# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# lasted modified by hongdongxuan at 20170914
import web
import json
import datetime
from mainapp.models.mongo.submit.paternity_test_mongo import PaternityTest as PT
from mainapp.controllers.project.pt_controller import PtController
from mainapp.libs.param_pack import *
from mainapp.libs.signature import check_sig
import re


class PtFamilyCombine(PtController):
    def __init__(self):
        super(PtFamilyCombine, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'paternity_test.report.pt_family_combine'
        task_type = 'workflow'
        params_json = {
            'dad_group': data.dad_group,
            'mom_id': data.mom_id,
            'preg_id': data.preg_id,
            'member_id': data.member_id,
            'dedup_all': data.dedup_all,
            'err_min': data.err_min,
        }
        if hasattr(data, 'new_mom_id'):
            params_json['new_mom_id'] = data.new_mom_id
            params_json['new_dad_id'] = data.new_dad_id
            params_json['new_preg_id'] = data.new_preg_id
        if hasattr(data, 'dad_id'):
            params_json['dad_id'] = data.dad_id
        if hasattr(data, 'dedup_start'):
            params_json['dedup_start'] = data.dedup_start
            params_json['dedup_end'] = data.dedup_end
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        if data.new_dad_id or data.new_mom_id or data.new_preg_id:
            task_id = "Report_" + data.dad_group + '_' + data.dad_id + '_' + data.mom_id + '_' + data.preg_id + '_' \
                      + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        else:
            if not data.dad_id:
                filter_name = "_"
            else:
                filter_name = data.dad_id
            task_id = "Dedup_" + data.dad_group + '_' + filter_name + '_' + data.mom_id + '_' + data.preg_id + '_' + \
                      datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = [
            ('params', params),
            ('name', task_id),
            ('desc', '??????????????????????????????'),
            ('member_id', data.member_id),
            ('status', 'start'),
            ('time', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_table_id = PT().insert_main_table('sg_pt_family_combine', mongo_data)
        update_info = {str(main_table_id): 'sg_pt_family_combine'}
        options = {
            'dad_group': data.dad_group,
            'mom_id': data.mom_id,
            'preg_id': data.preg_id,
            'err_min': data.err_min,
            'dedup_all': data.dedup_all,
            "update_info": json.dumps(update_info),
            'main_id': str(main_table_id),
            'member_id': data.member_id,
        }
        if hasattr(data, 'dad_id'):
            options['dad_id'] = data.dad_id
        if hasattr(data, 'new_mom_id'):
            options['new_mom_id'] = data.new_mom_id
            options['new_dad_id'] = data.new_dad_id
            options['new_preg_id'] = data.new_preg_id
        if hasattr(data, 'dedup_start'):
            options['dedup_start'] = data.dedup_start
            options['dedup_end'] = data.dedup_end
        sheet_data = self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params)
        # print "*********"
        # print sheet_data
        task_info = super(PtFamilyCombine, self).POST()
        return json.dumps(task_info)

    def check_options(self, data):
        """
        ??????????????????????????????,???????????????id\??????id\?????????id\???????????????\member_id\??????????????????
        :param data:????????????????????????(???????????????????????????????????????)
        :return: ????????????
        """
        params_name = ['mom_id', 'dad_group', 'preg_id', 'dad_id', 'err_min', 'member_id', 'dedup_all']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("????????????{}".format(names))
        if data.new_dad_id or data.new_mom_id or data.new_preg_id:
            if not data.new_dad_id or not data.new_mom_id or not data.new_preg_id:
                success.append("???????????????????????????????????????????????????????????????????????????????????????????????????")
            if not data.dad_id:
                success.append("??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????")
        # ???????????????????????????????????????????????????????????????????????????????????????
        if data.dedup_all == 'false':
            if data.dedup_start and data.dedup_end:
                if len(str(data.dedup_start)) != 4 or len(str(data.dedup_end)) != 4:
                    success.append("?????????????????????????????????????????????????????????????????????+??????????????????example???1708???17???08??????")
                try:
                    int(data.dedup_start) or int(data.dedup_end)
                except:
                    success.append("?????????????????????????????????????????????????????????example???1708???")
                if int(data.dedup_start) > int(data.dedup_end):
                    success.append("?????????????????????????????????????????????????????????????????????")
            # else:
            #     success.append("??????????????????????????????????????????????????????????????????")
        else:
            if data.dedup_start or data.dedup_end:
                success.append("????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????")
            pass
        return success



