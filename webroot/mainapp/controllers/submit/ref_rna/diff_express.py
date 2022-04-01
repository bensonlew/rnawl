# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng, modify by khl 20170428'
import web
import json
import random
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from mainapp.libs.param_pack import *
from biocluster.config import Config
from mainapp.models.mongo.submit.ref_rna.ref_diff import RefDiff
import types
from mainapp.models.mongo.meta import Meta
from mainapp.models.workflow import Workflow
from mainapp.controllers.project.ref_express_controller import RefExpressController
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from mbio.api.to_file.ref_rna import *
from mainapp.models.mongo.submit.ref_rna import *
import datetime


class DiffExpressAction(RefRnaController):
    def __init__(self):
        super(DiffExpressAction, self).__init__()
    
    def GET(self):
        return 'khl'

    @check_sig
    def POST(self):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        #print data.control_id

        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)
        my_param = dict()

        if str(data.group_id) != "all":
            return_control_id_group_detail = self.check_group_id_control_id(data)
            if return_control_id_group_detail:
                info = {"success": False, "info": '+'.join(return_control_id_group_detail)}
                return json.dumps(info)

        task_type = ''
        task_name = 'ref_rna.report.diff_express'

        my_param['express_id'] = data.express_id
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        my_param['control_id'] = data.control_id
        my_param['fc'] = data.fc
        my_param['pvalue_padjust'] = data.pvalue_padjust
        my_param['pvalue'] = data.pvalue
        my_param['diff_method'] = data.diff_method
        my_param['type'] = data.type  #基因还是转录本
        my_param['task_type'] = task_type
        my_param['submit_location'] = data.submit_location
        my_param['task_id'] = data.task_id

        if data.pvalue ==0 and data.fc == 1:
            info = {"success":False,"info":'{}值为0和fc为1不能同时存在!'.format(data.pvalue_padjust)}
            return json.dumps(info)
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))
        express_info = self.ref_rna.get_main_info(data.express_id, 'sg_express')
        task_info = self.ref_rna.get_task_info(express_info['task_id'])
        express_params=json.loads(express_info["params"])
        express_method = express_params["express_method"]
        value_type = data.type  #gene or transcript
        diff_method = data.diff_method.lower()

        if express_info:
            express_params = json.loads(express_info['params'])
            express_level = express_params['type']  #fpkm or tpm
            # DiffExp_基因\转录本_软件_tpm\fpkm_差异分析软件_日期_时间
            re_name_info = {"gene":"G","transcript":"T","featurecounts":"FeaCount","rsem":"RSEM","edger":"ER","deseq2":"DS","degseq":"DG"}
            re_value_type = re_name_info[value_type.lower()]
            re_express_method = re_name_info[express_method.lower()]
            re_express_level = express_level.lower()
            re_diff_method = re_name_info[diff_method]
            main_table_name = "DiffExp_{}_{}_{}_{}_".format(re_value_type,re_express_method,re_express_level,re_diff_method) + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            task_id = express_info["task_id"]
            project_sn = express_info["project_sn"]

            mongo_data = [
                ('project_sn', task_info['project_sn']),
                ('task_id', task_info['task_id']),
                ('status', 'start'),
                ('desc',"表达量差异主表"),
                ('name', main_table_name),
                ("value_type",express_level),
                ("express_id",ObjectId(data.express_id)),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
            ]
            #if express_info["is_duplicate"] and express_info['trans'] and express_info['genes']:
            # if "is_duplicate" in express_info.keys():
            #     mongo_data.append(
            #         ("is_duplicate", express_info["is_duplicate"]))
            if str(data.group_id) == 'all':  #判断is_duplicate参数
                mongo_data.append(('is_duplicate',False))
            else:
                mongo_data.append(('is_duplicate',True))

            if data.type == 'gene':
                mongo_data.extend([
                    ("trans", False),
                    ('genes', True)
                ])   #参数值方便前端取数据
            if data.type == 'transcript':
                mongo_data.extend([
                    ('trans',True),
                    ('genes',False)
                ])
            collection_name = "sg_express_diff"
            main_table_id = self.ref_rna.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): collection_name}
            update_info = json.dumps(update_info)
            try:
                class_code_id = self.ref_rna.get_class_code_id(task_id)
            except Exception:
                pass
            options = {
                "express_file": data.express_id,
                "update_info": update_info,
                "type":data.type,  #gene/ transcript
                "control_file": data.control_id,
                'fc': data.fc,
                "express_method": express_method,
                "group_id_id":data.group_id,
                'class_code': class_code_id,
                "diff_method":data.diff_method,
                "diff_express_id": str(main_table_id),
                "log":"None",
                "express_level":express_level,
                "pvalue_padjust":data.pvalue_padjust,
                "pvalue":data.pvalue,
                "class_code_type":"express_diff"
                # "group_id": data.group_id,
                # "group_detail":data.group_detail,
            }
            # options 中的type参数在差异分析的流程中已写死，传给export_class_code函数
            to_file = ["ref_rna.export_express_matrix_level(express_file)",  "ref_rna.export_control_file(control_file)", "ref_rna.export_class_code(class_code)"]
            if data.group_id != 'all':
                options.update({
                    "group_id": data.group_id,
                    "group_detail": data.group_detail,
                })
                to_file.append("ref_rna.export_group_table_by_detail(group_id)")
            self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, 
                            to_file=to_file,params=my_param, project_sn=task_info['project_sn'], task_id=task_info['task_id'])
            task_info = super(DiffExpressAction, self).POST()
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
            return json.dumps(task_info)
        
        else:
            info = {"success": False, "info": "express_id不存在，请确认参数是否正确！!"}
            return json.dumps(info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ['express_id', 'fc','group_detail', 'group_id', 'control_id', \
                        'submit_location','pvalue_padjust','pvalue','diff_method','type']
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数{}!".format(names))
        for ids in [data.express_id, data.group_id, data.control_id]:
            ids = str(ids)
            if not isinstance(ids, ObjectId) and not isinstance(ids, types.StringTypes):
                success.append("传入的id：{}不是一个ObjectId对象或字符串类型".format(ids))
        return success

    def check_group_id_control_id(self, data):
        """检测control_id的样本分组信息是否和group_detail表一一对应"""
        control_id, group_detail = data.control_id,json.loads(data.group_detail)
        compare_names = self.ref_rna.get_control_id(control_id)
        group_names = group_detail.keys()
        #print group_names
        success = []
        unique_compare_names=[]
        for i in compare_names:
            compare_id = i.split("|")
            for j in compare_id:
                if j not in group_names:
                    if j not in unique_compare_names:
                        unique_compare_names.append(j)
                        success.append("分组方案和对照组方案不一致，分组方案没有选择{}".format(str(j),str(j)))
                    else:
                        pass
        group_size = list()
        for group in group_detail:
            group_size.append(len(group_detail[group]))
        group_size.sort()
        if group_size[0] == group_size[-1] == 1:
            if data.diff_method == "DESeq2":
                success.append('DESeq2不适合处理单样本和单样本比较， 请选择DEGseq或edgeR')
        elif group_size[0] == 1 and group_size[1] >= 2:
            if data.diff_method == "DESeq2" or data.diff_method == "DEGseq":
                success.append('你的分组方案表明可能要进行单样本和多样本的比较，'
                               '此时请选择edgeR做差异分析或者重新设计分组方案')
        elif group_size[0] >= 2:
            if data.diff_method == "DEGseq":
                success.append('只涉及包含多样本的组与组间比较时，我们只推荐DESeq2或者edgeR')
        return success
