# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from mainapp.models.mongo.submit.med_mongo import MedMongo as PT
from mainapp.controllers.project.nipt_controller import NiptController
from mainapp.libs.signature import check_sig
from mainapp.libs.param_pack import *
from bson import ObjectId
from bson import SON
import datetime
import json
import web


class NiptAction(NiptController):
    def __init__(self):
        """
        无创产筛分析的接口
        modified by hongdong 20171201
        用于激发工作流
        """
        super(NiptAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        task_name = 'medical.nipt_v2.nipt'
        task_type = 'workflow'
        mongo_data = [
            ('batch_id', ObjectId(data.batch_id)),
            ("type", "nipt"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ("status", "start"),
            ("urgent", data.urgent),  # nipt-emergency, nipt-hurry, nipt-normal
            ("member_id", data.member_id),
            ("desc", "产前筛查正在分析中"),
            ("end_counts", 0),
            ("all_counts", len(str(data.sample_list).strip().split(','))),
            ("is_show", "1")  # 为1页面展示，为2页面不展示
        ]
        main_table_id = PT('pt_v2').add_main_id('sg_analysis_status')
        update_info = {str(main_table_id): 'sg_analysis_status'}
        collection = 'sg_analysis_status'
        options = {
            'fastq_path': data.fastq_path,
            # 'batch_id': data.batch_id,
            'member_id': data.member_id,
            'update_info': json.dumps(update_info),
            'ref_group': 1,
            # 'urgent': data.urgent,
            'main_table_data': SON(mongo_data),
            'single': 'false' if str(data.single) == 'PE' else "true",
            'sample_list': data.sample_list,
            'sanger_type': "sanger" if str(data.client) == 'client01' else 'tsanger'
        }
        self.set_sheet_data(name=task_name, options=options, module_type=task_type, db_type='pt_v2')
        task_info = super(NiptAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': collection}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传入的参数fastq_path:产筛的fastq文件所在目录，urgent：用于区分是加急还是非加急，sample_list：是所有要进行产筛
        分析的样本
        :param data:网页端传入的数据
        :return: 检查结果
        """
        success = []
        params_name = ['member_id', 'batch_id', 'fastq_path', 'client', 'urgent', 'sample_list', 'single']
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数{}".format(names))
        if data.client not in ['client01', 'client03']:
            success.append('{}不在client合法范围内！'.format(data.client))
        if data.single not in ['PE', 'SE']:
            success.append('{}不在[PE,SE]合法范围内！'.format(data.single))
        return success

