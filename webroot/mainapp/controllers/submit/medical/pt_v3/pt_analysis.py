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


class PtAnalysisAction(PtController):
    def __init__(self):
        """
        亲子鉴定的家系分析的接口
         modified by hongdong 20180821
        当data中有samples这个字段的时候，该接口就是进行call snp工作，如果不存在该字段就进行亲子的家系分析
        """
        super(PtAnalysisAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'medical.paternity_test_v3.report.pt_analysis'
        task_type = 'workflow'
        options = {
            'fastq_path': data.fastq_path,
            'batch_id': data.batch_id,
            'member_id': data.member_id
        }
        if hasattr(data, 'samples'):  # call snp
            mongo_data = [
                ('member_id', data.member_id),
                ('batch_id', ObjectId(data.batch_id)),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ('desc', '亲子鉴定call snp正在分析'),
                ('status', 'start'),
                ('type', 'pt'),
                ('snp_all_counts', len(json.loads(data.samples))),
                ('snp_end_counts', 0),
                # ('family_end_counts', 0),
                # ('family_all_counts', 0),
                ('is_show', "1")  # 如果运行多次的时候，会有重复的数据在页面显示，这里加个字段进行过滤，1的时候，在
                # 页面展示，2的时候 不展示
            ]
            if str(data.ana_type) == 'F-M':
                mongo_data.append(('family_end_counts', 0))
                mongo_data.append(('family_all_counts', 0))
            main_table_id = PT('pt_v3').add_main_id('sg_analysis_status')
            update_info = {str(main_table_id): 'sg_analysis_status'}
            collection = 'sg_analysis_status'
            params = ''
            analysis_name = 'snp'  # 该信息会在workspace中进行体现用于区分分析的类型，好处是可以快速定位到结果目录
            options.update({
                "update_info": json.dumps(update_info),
                'main_table_data': SON(mongo_data),
                'samples': data.samples
            })
        else:  # 家系分析
            if hasattr(data, 'old_mom_id'):
                son_id = data.son_id
            else:
                son_id = "{}_{}".format(data.son_id, self.find_sample_id(data.mom_id))
            params_json = {
                'dad_id': data.dad_id,
                'mom_id': data.mom_id,
                'son_id': son_id
            }
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
            family_result = PT("pt_v3").get_one("sg_family", "case_name", str(data.dad_id).strip().split('-')[0])
            if not family_result:
                info = {"success": False, "info": "家系编号{}在数据库中不存在！".format(str(data.dad_id).strip().split('_')[0])}
                return json.dumps(info)
            else:
                family_id = family_result['_id']
            main_table_name = data.dad_id + '_' + data.mom_id + '_' + son_id
            mongo_data = [
                ('member_id', data.member_id),
                ('batch_id', ObjectId(data.batch_id)),
                ('family_id', family_id),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ('params', params),
                ('name', main_table_name),
                ('dad_id', data.dad_id),
                ('mom_id', str(data.mom_id).strip().split('_')[0]),
                ('preg_id', son_id),
                ('desc', '亲子鉴定家系分析'),
                ('status', 'start'),
                ('report_name', "_".join([data.dad_id, str(data.mom_id).strip().split('_')[0]])),
                ('remarks', ""),
                ('result', ''),
                ('problem_sample_num', ''),
                ('batch_dedup', ''),
                ('ms_match', "YES"),
                ('problem_sample_name', ''),
                ('is_report', "2"),
                ('is_show', "1")  # 用于后面的多个end状态的时候取那条记录展示，1的时候展示，2的时候不展示
            ]
            main_table_id = PT('pt_v3').add_main_id('sg_father')
            update_info = {str(main_table_id): 'sg_father', 'batch_id': str(data.batch_id)}
            collection = 'sg_father'
            analysis_name = 'pt_{}_{}_{}'.format(data.dad_id, "-".join(str(data.mom_id).split('-')[1:]),
                                                 "-".join(str(data.son_id).split('-')[1:]))
            options.update({
                'mom_id': data.mom_id,
                'son_id': data.son_id,
                'err_min_num': data.err_min_num,
                "update_info": json.dumps(update_info),
                'dad_id': data.dad_id,
                'main_table_data': SON(mongo_data),
                'family_id': str(family_id)
            })
            if hasattr(data, 'old_mom_id'):
                options.update({"old_mom_id": data.old_mom_id})
        self.set_sheet_data_(name=task_name, options=options, module_type=task_type, params=params,
                             db_type='pt_v3', analysis_name=analysis_name)
        task_info = super(PtAnalysisAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': collection}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传入的参数,必须有母本id\胎儿id\父本组id\member_id\batch_id,检查父本与母本的tab文件是否存在,两者之一不存在就报错
        :param data:网页端传入的数据(是否全库查重的参数必须传递)
        :return: 检查结果
        """
        success = []
        if not hasattr(data, 'samples'):
            params_name = ['mom_id', 'dad_id', 'son_id', 'member_id', 'batch_id', 'fastq_path', 'err_min_num', 'client']
        else:
            params_name = ['member_id', 'batch_id', 'fastq_path', 'client']
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数{}".format(names))
        if data.client not in ['client01', 'client03']:
            success.append('{}不在client合法范围内！'.format(data.client))
        if hasattr(data, "samples"):
            if len(json.loads(data.samples)) == 0:
                success.append("流程进行call snp分析，但是samples列表为空，没有样本进行call snp")
        return success

    def find_sample_id(self, sample_id):
        """
        根据正则表达式去匹配查找对应的样本id，这里要注意的是'WQ181495-2-M-1' 与'WQ181495-M-1' 其实是同一个家系
        re.match("(WQ.*)-(\d)-(M|F|S)(.*)", a)
        :return:
        """
        m = re.match("(WQ.*)-(\d)-(M|F|S)(.*)", sample_id)
        if m:
            new_id = m.group(3) + m.group(4)
        else:
            new_id = "-".join(sample_id.split('-')[1:])
        return new_id
