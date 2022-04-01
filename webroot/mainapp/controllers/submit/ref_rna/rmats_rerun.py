# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/12 09:59

import re, os, Bio, argparse, sys, fileinput, urllib2
import web
import json
import random
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from mainapp.libs.param_pack import *
from biocluster.config import Config

import types
from mainapp.models.mongo.meta import Meta
from mainapp.models.workflow import Workflow

from mbio.api.to_file.ref_rna import *
from mainapp.models.mongo.ref_rna import RefRna
from mainapp.controllers.project.ref_rna_controller import RefRnaController
from pymongo import MongoClient
import subprocess


def isFloat(string):
    try:
        temp = float(string)
        return True
    except Exception as e:
        return False


class RmatsRerunAction(RefRnaController):
    '''
    rmats接口
    
    '''
    
    def __init__(self):
        super(RmatsRerunAction, self).__init__(instant=False)
    
    def GET(self):
        return 'jlf'
    
    @check_sig
    def POST(self):
        data = web.input()
        post_args = ['group_id', 'group_detail', 'cut_off', 'submit_location', 'splicing_id', "compare_plan", "task_type", "control_id"]
        # print data.control_id
        for arg in post_args:
            if not hasattr(data, arg):
                info = {'success': False, 'info': '%s参数缺少!' % arg}
                return json.dumps(info)
        task_name = 'ref_rna.report.rmats'
        # task_type = 'workflow'
        # control_doc = RefRna().get_main_info(ObjectId(data.control_id), 'sg_specimen_group_compare')
        # print control_doc
        group_doc = RefRna().get_main_info(ObjectId(data.group_id), 'sg_specimen_group')
        #print group_doc
        #print(data.cut_off)

        # compare_plans = re.sub(r'[\"\[\]\s+\']', '', control_doc['compare_names']).split(',')
        # print compare_plans
        # if len(compare_plans) > 1:
        #     info = {"success": False, "info": "只能选一组对照组方案"}
        #     return json.dumps(info)
        #print(data.compare_plan)
        case_group_name = data.compare_plan.split('|')[0]
        #print case_group_name
        control_group_name = data.compare_plan.split('|')[1]
        #print control_group_name

        # case_group_members = group_doc["specimen_names"][
        #     group_doc["category_names"].index(case_group_name)].keys()
        group_dict = json.loads(data.group_detail)
        case_group_members = group_dict[case_group_name]
        #print case_group_members
        # control_group_members = group_doc["specimen_names"][
        #     group_doc["category_names"].index(control_group_name)].keys()
        control_group_members = group_dict[control_group_name]
        #print control_group_members

        case_group_bam_lst = []
        control_group_bam_lst = []
        for case_sp in case_group_members:
            doc = RefRna().get_main_info(ObjectId(case_sp), 'sg_specimen')
            #print('现在装载%s样本的bam信息: %s' % (case_sp, doc))
            case_group_bam_lst.append(
                doc['bam_path'])
        for control_sp in control_group_members:
            doc = RefRna().get_main_info(ObjectId(control_sp), 'sg_specimen')
            #print('现在装载%s样本的bam信息: %s' % (control_sp, doc))
            control_group_bam_lst.append(
                doc['bam_path'])
        case_group_bam_str = ','.join([p.strip() for p in case_group_bam_lst])
        control_group_bam_str = ','.join([p.strip() for p in control_group_bam_lst])
        # if len(case_group_bam_lst) <= 1 or len(control_group_bam_lst) <= 1:
        #     info = {"success": False, "info": "每组必须有2个以上的样本"}
        #     return json.dumps(info)
        # if not isFloat(data.cut_off):
        #     info = {"success": False, "info": "cut_off必须为浮点数"}
        #     return json.dumps(info)
        if float(data.cut_off) >= 1 or float(data.cut_off) < 0:
            info = {"success": False, "info": "cut_off必须在[0,1)之间"}
            return json.dumps(info)
        my_param = {}
        #print data.group_detail
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        my_param['compare_plan'] = data.compare_plan
        my_param['splicing_id'] = data.splicing_id
        my_param['cut_off'] = data.cut_off
        my_param['control_id'] = data.control_id
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = data.task_type
        splicing_info = RefRna().get_main_info(ObjectId(data.splicing_id), 'sg_splicing_rmats')
        if splicing_info:
            # old_params = json.loads(splicing_info['params'])
            task_info = self.ref_rna.get_task_info(splicing_info['task_id'])
            main_table_name = "SplicingRmats_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
            task_id = splicing_info["task_id"]
            project_sn = splicing_info["project_sn"]
            # if not splicing_info["ref_gtf"]:
            #     raise Exception('关联主表没有设置ref gtf路径')
            # if not 'ref_gtf' in old_params.keys():
            #     raise Exception('关联主表没有设置ref gtf路径')
            # if not 'seq_type' in old_params.keys():
            #     raise Exception('关联主表没有设置seq_type')
            # if not 'read_len' in old_params.keys():
            #     raise Exception('关联主表没有设置read_len')
            # if not 'novel' in old_params.keys():
            #     raise Exception('关联主表没有设置novel')
            # if not 'analysis_mode' in old_params.keys():
            #     raise Exception('关联主表没有设置analysis_mode')
            # if not 'lib_type' in old_params.keys():
            #     raise Exception('关联主表没有设置lib_type')
            # my_param["ref_gtf"] = old_params['ref_gtf']
            # my_param["seq_type"] = old_params['seq_type']
            # my_param["read_len"] = old_params['read_len']
            # my_param["novel"] = old_params['novel']
            # my_param["analysis_mode"] = old_params['analysis_mode']
            # my_param["lib_type"] = old_params['lib_type']
            # if my_param['analysis_mode'] == 'P' and len(case_group_bam_lst) != len(control_group_bam_lst):
            #     info = {"success": False, "info": "分析模式为paired的时候，实验组和对照组的样本数必须相等"}
            #     return json.dumps(info)
            
            # chr_set = [e.strip() for e in subprocess.check_output(
            #     'awk -F \'\\t\'  \'$0!~/^#/{print $1}\' %s  | uniq | sort |uniq ' % old_params['ref_gtf'],
            #     shell=True).strip().split('\n')]
            ref_gtf = splicing_info["ref_gtf"]
            chr_set = splicing_info['chr_set']
            group = {case_group_name: 's1', control_group_name: 's2'}
            mongo_data = [
                ('project_sn', task_info['project_sn']),
                ('task_id', task_info['task_id']),
                ('status', 'start'),
                ('desc', "rmats主表"),
                ('name', main_table_name),
                ('chr_set', chr_set),
                ('ref_gtf', ref_gtf),
                ('group', group),
                ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("params", json.dumps(my_param, sort_keys=True, separators=(',', ':')))
            ]
            collection_name = "sg_splicing_rmats"
            main_table_id = self.ref_rna.insert_main_table(collection_name, mongo_data)
            update_info = {str(main_table_id): collection_name}
            update_info = json.dumps(update_info)
            options = {
                "update_info": update_info,
                "splicing_id": str(main_table_id),
                # "ref_gtf": old_params['ref_gtf'],
                "ref_gtf": ref_gtf,
                # "chr_set": splicing_info['chr_set'].__str__,
                "cut_off": float(data.cut_off),
                'control_file': data.control_id,
                "group_id": data.group_id,
                "group_detail": data.group_detail
            }
            to_file = ["ref_rna.export_control_file(control_file)"]
            if data.group_id != 'all':
                to_file.append("ref_rna.export_group_table_by_detail(group_id)")
                options.update({

                    "case_group_bam_str": case_group_bam_str,
                    "control_group_bam_str": control_group_bam_str,
                    "case_group_name": case_group_name,
                    "control_group_name": control_group_name
                })
            else:
                info = {"success": False, "info": "不能选单个样本作为一个组 组必须包含多个样本"}
                return json.dumps(info)
            #print 'workflow设置的option为: %s' % options
            self.set_sheet_data(name=task_name, options=options, main_table_name=main_table_name, task_id=task_id,
                                project_sn=project_sn, module_type='workflow', to_file=to_file, params=my_param)
            task_info = super(RmatsRerunAction, self).POST()
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
            #print task_info
            return json.dumps(task_info)
        else:
            info = {"success": False, "info": "splicing_id不存在，请确认参数是否正确！!"}
            return json.dumps(info)
    
    # def check_options(self, data):
    #     params_name = ['group_id', 'group_detail', 'cut_off', 'analysis_mode', 'submit_location', 'splicing_id',
    #                    'control_id']
    #     success = []
    #     for names in params_name:
    #         if not (hasattr(data, names)):
    #             success.append("缺少参数!")
    #
    #     return success
