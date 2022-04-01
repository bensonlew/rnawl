# -*- coding: utf-8 -*-
# __author__ = 'khl'  
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
#from mainapp.controllers.project.ref_express_controller import RefExpressController
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from bson import ObjectId
from mainapp.libs.signature import check_sig


class ExpressCorrAction(RefRnaController):
    def __init__(self):
        super(ExpressCorrAction, self).__init__(instant=True)
    def GET(self):
        return 'khl'
    @check_sig
    def POST(self):
        data = web.input()
        postArgs = ['group_id', 'group_detail', 'submit_location', 'express_id', 'method']
        # postArgs删除了 hclust_method 参数
        for args in postArgs:
            if not hasattr(data, args):
                info = {'success': False, 'info': '%s参数缺少!' % args}
                return json.dumps(info)
        express_info = self.ref_rna.get_main_info(data.express_id, 'sg_express')
        task_name = 'ref_rna.report.express_corr'
        task_type = 'workflow'
        my_param = dict()
        my_param['express_id'] = data.express_id
        my_param['group_id'] = data.group_id
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['method'] = data.method
        # my_param['hclust_method'] = data.hclust_method
        my_param['submit_location'] = data.submit_location
        my_param["task_type"] = task_type
        group_detail = my_param["group_detail"]
        if data.group_id in ["all", "ALL", "All"]:
            sample_num = len(group_detail.values()[0])
            if sample_num<2:
                info = {'success': False, 'info': '表达量相关性分析至少选择两个样本!'}
                return json.dumps(info)
        else:
            sample_total = []
            for keys,values in group_detail.items():
                for sam in values:
                    if sam not in sample_total:
                        sample_total.append(sam)
                    else:
                        info = {'success': False, 'info': '不同分组中有重复的样本名字!'}
                        return json.dumps(info)
            if len(sample_total)<2:
                info = {'success': False, 'info': '表达量相关性分析至少选择两个样本!'}
                return json.dumps(info)

        if express_info:
            # task_info = self.ref_rna.get_task_info(express_info['task_id'])
            # print task_info
            # main_table_name = "ExpressCorr_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            params = json.loads(express_info['params'])
            express_method = params["express_method"].lower()
            express_level = params["type"].lower()
            re_name_info = {"featurecounts": "FeaCount", "rsem": "RSEM"}
            re_express_method = re_name_info[express_method]
            main_table_name = "ExpCor_{}_{}_".format(re_express_method,express_level) +  str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            mongo_data = [
                ('project_sn', express_info['project_sn']),
                ('task_id', express_info['task_id']),
                ('status', 'start'),
                ('desc',"样本间相关性分析"),
                ('name', main_table_name),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
            ]
            collection_name = "sg_express_correlation"
            main_table_id = self.ref_rna.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): collection_name}
            
            params=json.loads(express_info['params'])
            express_level = params['type']
            
            options = {
                "express_file": data.express_id,
                "method": data.method,
                "correlation_id": str(main_table_id),
                "group_id": data.group_id,
                "group_detail": data.group_detail,
                "corr_pca": "corr",
                "type": "gene",  # 给workflow传参使用
                "update_info": json.dumps(update_info),
                "express_level":express_level
            }
            to_file = ["ref_rna.export_express_matrix_level(express_file)", "ref_rna.export_group_table_by_detail(group_id)"]
            self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                                task_id=express_info['task_id'], project_sn=express_info['project_sn'],
                                params=my_param, to_file=to_file)
            task_info = super(ExpressCorrAction, self).POST()
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }}
            #print "Shenghe: " + json.dumps(task_info)
            return json.dumps(task_info)

        else:
            info = {"success": False, "info": "表达量表不存在，请确认参数是否正确！!"}
            return json.dumps(info)
