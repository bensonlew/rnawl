# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170425'
import web
import json
import datetime
from mainapp.controllers.project.meta_controller import MetaController
#from mainapp.controllers.project.ref_express_controller import RefExpressController
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mainapp.libs.param_pack import group_detail_sort
from mbio.api.to_file.ref_rna import *
from bson import ObjectId
from mainapp.libs.signature import check_sig

class ExpressVennAction(RefRnaController):
    def __init__(self):
        super(ExpressVennAction, self).__init__(instant=True)
    
    def GET(self):
        return 'khl'

    @check_sig
    def POST(self):
        data = web.input()
        # args = ["express_id", "group_id", "group_detail", "type", "submit_location"]
        args = ["express_id", "group_id", "group_detail", "type", "threshold", "submit_location"]
        for arg in args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '%s参数缺少!' % arg}
                return json.dumps(info)
        group_detal_dict = json.loads(data.group_detail)
        if data.group_id not in ['all', 'All', 'ALL']:
            if len(group_detal_dict) < 2:
                info = {"success": False, "info": "进行Venn分析，分组方案的分组类别必须大于等于2且小于等于6！"}
                return json.dumps(info)
            if len(group_detal_dict) > 6:
                info = {"success": False, "info": "进行Venn分析，分组方案的分组类别必须大于等于2且小于等于6！"}
                return json.dumps(info)
        else:
            if len(group_detal_dict.values()[0])<2 or len(group_detal_dict.values()[0])>=7:
                info = {"success": False, "info": "进行Venn分析，分组方案的样本数量必须大于等于2且小于等于6！"}
                return json.dumps(info)
        #print data.group_detail
        task_name = "ref_rna.report.express_venn"
        task_type = "workflow"
        my_param = dict()
        my_param["express_id"] = data.express_id
        my_param["group_detail"] = group_detail_sort(data.group_detail)
        my_param["group_id"] = data.group_id
        my_param["submit_location"] = data.submit_location
        my_param["task_type"] = task_type
        my_param["type"] = data.type
        my_param['threshold'] = data.threshold #根据fpkm/tpm过滤

        express_info = self.ref_rna.get_main_info(data.express_id, 'sg_express')
        task_info = self.ref_rna.get_task_info(express_info['task_id'])
        # main_table_name = 'ExpressVenn_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        # 'ExpVenn_基因\转录本_软件_tpm/fpkm_筛选阈值_日期_时间'

        group_detail_dict = json.loads(data.group_detail)

        if express_info:
            re_name_info = {"gene": "G", "transcript": "T", "featurecounts": "FeaCount", "rsem": "RSEM", "edger": "ER",
                            "deseq2": "DS", "degseq": "DG"}
            re_query_type = re_name_info[data.type.lower()]
            params = json.loads(express_info['params'])
            re_express_method = re_name_info[params['express_method'].lower()]
            express_level = params['type']
            threshold =  str(data.threshold)
            main_table_name = "ExpVenn_{}_{}_{}_{}_".format(re_query_type,re_express_method,express_level,threshold) + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3])
            # 基因\转录本_软件_tpm / fpkm_筛选阈值_日期_时间
            mongo_data = [
                ('project_sn', task_info['project_sn']),
                ('task_id', task_info['task_id']),
                ('status', 'start'),
                ('desc',"样本间特异性venn图分析"),
                ('name', main_table_name),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
            ]
            collection_name = "sg_express_venn"
            main_table_id = self.ref_rna.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): collection_name}
            
            params = json.loads(express_info['params'])
            express_level = params['type']
            
            options = {
                "express_file": data.express_id,
                "group_id": data.group_id,
                "group_detail": data.group_detail,
                "update_info": json.dumps(update_info),
                "venn_id": str(main_table_id),
                "express_level":express_level,  #传fpkm还是tpm值给workflow
                "type": data.type,  #给workflow传参用的
                "threshold": data.threshold  #过滤表达量
            }
            to_file = ["ref_rna.export_express_matrix_level(express_file)", "ref_rna.export_group_table_by_detail(group_id)"]
            self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name,
                                task_id=express_info['task_id'], project_sn=express_info['project_sn'],
                                params=my_param, to_file=to_file)
            task_info = super(ExpressVennAction, self).POST()
            task_info['content'] = {
                'ids': {
                    'id': str(main_table_id),
                    'name': main_table_name
                }
            }
            return json.dumps(task_info)
        else:
            info = {"success": False, "info": "表达量表不存在，请确认参数是否正确！!"}
            return json.dumps(info)
