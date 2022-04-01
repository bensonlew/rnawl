# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# created at 20171109

import web
import json
import datetime
import os
from mainapp.libs.signature import check_sig
from mainapp.controllers.project.pt_controller import PtController as MedController
from mainapp.models.mongo.submit.med_mongo import MedMongo
from bson import ObjectId
import xml.etree.ElementTree as ET


class DataSplitAction(MedController):
    """
    医学数据拆分程序接口
    该接口关联可以同时进行多个分析流程的启动，也可以根据project_types进行设置，针对具体的某一个流程的分析，project_types按照项目
    类型来设定:亲子鉴定：pt；无创产筛：nipt；当project_types=all的时候所有的流程都会启动
    表的更新逻辑：
    1）默认的参数是board_batch， member_id， project_types， start_type（可选）
    2）由于自动拆分与页面手动拆分的时候，传入的板号信息是不同的，自动拆分的传入的是具体的板子的路径，而页面手动点击的时候传入的就是
    板号名字（前端不想查找，所以就我们这边进行个判断）
       2.1当start_type存在且等于auto的时候，也就是第一次分析，并且是自动拆分，首先检查该板子的信息在sg_datasplit中有没有，
       如果没有，则报错终止，提示要先进行文件的上传，并将该板子的一些信息更新到sg_datasplit中用于后面的手动拆分；如果存在，
       将板子的具体的路径更新到board_path这个字段（用于后面分析的时候能够找到板子所在的路径），下面激活拆分
       2.2 当start_type不存在，也就是页面进行手动拆分，首先检查该板子的信息在sg_datasplit中有没有，没有的话一样报错，这个情况基
       本上不存在，我们只讨论存在，手动拆分肯定是进行重运行，首先要将旧表的desc与status个更新下，以保证页面上显示的拆分状态在拆分，
       然后启动后面的工作流
    last modified by HONGDONG 20171211
    """
    def __init__(self):
        super(DataSplitAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        params_name = ['board_batch', 'member_id', "project_types"]
        for param in params_name:
            if not hasattr(data, param):
                info = {'success': False, 'info': '缺少{}参数'.format(param)}
                return json.dumps(info)
        if str(data.project_types) not in ['pt', 'nipt', 'all']:
            info = {'success': False, 'info': '项目类型{}不合法！'.format(data.analysis_types)}
            return json.dumps(info)
        if hasattr(data, "start_type"):  # 启动类型分为自动与手动
            if str(data.start_type) not in ['auto']:
                info = {'success': False, 'info': '启动类型{}不合法！'.format(data.start_type)}
                return json.dumps(info)
            else:
                board_batch_name = os.path.basename(str(data.board_batch).rstrip())
        else:
            board_batch_name = data.board_batch
        self.check_can_split(board_batch_name)
        board = MedMongo("pt_v2").get_one('sg_datasplit', 'board_batch', board_batch_name)
        if not board:
            info = {'success': False, 'info': '样本信息未查到，请先上传样本信息表'}
            # 在sg_board_info表中插入一条板记录标记该板已下机,提示异常情况
            split_data = {
                "board_batch": board_batch_name,   # 板号
                "board_path": data.board_batch,    # 板子具体的路径
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "wq_sample": '-',
                "ws_sample": '-',
                "wq_num": '-',
                "ws_num": '-',
                "status": 'produced',  # 状态包括：unproduced（未下机）、produced（已下机）、start（开始拆分）、
                # end（拆分完成）、failed（失败）。
                "start_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "desc": '数据下机时自动拆分失败：样本信息未查到，请检查是否上传样本相关信息表',
                "end_time": '-'
                }
            datasplit_batch_id = MedMongo("pt_v2").add_one('sg_datasplit', split_data)
            record_data = {
                "batch_id": datasplit_batch_id,
                "board_batch": data.board_batch,
                "status": "failed",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "member_id": data.member_id,
                "desc": '数据下机时自动拆分失败：样本信息未查到，请检查是否上传样本相关信息表'
            }
            MedMongo("pt_v2").add_one('sg_datasplit_records', record_data)
            return json.dumps(info)
        else:
            if hasattr(data, "start_type") and str(data.start_type) == 'auto':
                MedMongo("pt_v2").update_table("sg_datasplit", "_id", board['_id'], {'board_path': data.board_batch})
            init_data = {'status': "start", "desc": "", "created_ts": datetime.datetime.now().strftime("%Y-%m-%"
                                                                                                       "d %H:%M:%S")}
            MedMongo("pt_v2").update_table("sg_datasplit", "_id", board['_id'], init_data)

        if data.project_types in ['all', 'pt']:
            ana = MedMongo("pt_v2").get_one_bymore('sg_analysis_status',
                                                   {'batch_id': ObjectId(board['_id']), 'type': 'pt'})
            if ana:
                analysis_data = {
                    "member_id": data.member_id,
                    "desc": '重新进行分析',
                    "is_show": '2'   # 更新为2的时候 页面就不展示
                    # "snp_end_counts": 0,
                    # "family_end_counts": 0
                }
                MedMongo("pt_v2").update_table("sg_analysis_status", "batch_id", ObjectId(board['_id']), analysis_data)
        if data.project_types in ['all', 'nipt']:
            ana = MedMongo("pt_v2").get_one_bymore('sg_analysis_status', {'batch_id': ObjectId(board['_id']),
                                                                          'type': 'nipt'})
            if ana:
                analysis_data = {
                    "member_id": data.member_id,
                    "desc": '重新进行分析',
                    "is_show": "2"
                    # "end_counts": 0
                }
                MedMongo("pt_v2").update_table("sg_analysis_status", "batch_id", ObjectId(board['_id']), analysis_data)
        # 查到板号信息，开始拆分
        datasplit_batch_id = board['_id']
        # record_data = {
        #     "batch_id": datasplit_batch_id,
        #     "board_batch": board_batch_name,
        #     "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        #     "member_id": data.member_id,
        #     "desc": '数据下机开始自动拆分'
        # }
        # 添加一条记录
        # MedMongo("pt_v2").add_one('sg_datasplit_records', record_data)
        # 判断拆分类型
        board_batch_path = str(data.board_batch) if hasattr(data, "start_type") and str(data.start_type) == 'auto' \
            else board['board_path']
        run_parameters = os.path.join(board_batch_path, "RunParameters.xml")
        if not os.path.exists(run_parameters):
            info = {'success': False, 'info': '{}:文件不存在'.format(run_parameters)}
            return json.dumps(info)
        tree = ET.parse(run_parameters)
        root = tree.getroot()
        if root.find("IsPairedEnd").text == 'true':
            split_type = 'PE'
        elif root.find("IsPairedEnd").text == 'false':
            split_type = 'SE'
        else:
            info = {'success': False, 'info': '判断split_type出错'}
            return json.dumps(info)
        MedMongo("pt_v2").update_table('sg_datasplit', 'board_path', board_batch_path, {'split_type': split_type})
        # 调用WPM运行数据拆分
        task_name = 'medical.data_split'
        options = {
            "split_tab": board['split_info_path'],
            "batch_id": str(datasplit_batch_id),
            "board_batch": str(data.board_batch) if hasattr(data, "start_type") and str(data.start_type) == 'auto'
            else board['board_path'],
            "split_type": split_type,
            "member_id": data.member_id,
            "update_info": json.dumps({str(datasplit_batch_id): 'sg_datasplit'}),
            "project_types": data.project_types
        }
        self.set_sheet_data_(name=task_name, options=options, module_type="workflow", db_type='pt_v2',
                             analysis_name='med_datasplit')
        print options
        info = super(DataSplitAction, self).POST()
        return json.dumps(info)

    def check_can_split(self, board_batch):
        """
        用于检查板子是否下机了，不下机不能往后面运行
        add by HONGDONG 20180305
        :param board_batch:
        :return:
        """
        board_one = "/mnt/clustre/upload/nextseq1/" + board_batch
        board_two = "/mnt/clustre/upload/nextseq/" + board_batch
        if os.path.exists(board_one):
            if not os.path.exists(board_one + "/RTAComplete.txt") or not os.path.isfile(board_one + "/RTAComplete.txt"):
                info = {'success': False, 'info': '板子{}还没有完全下机不能进行拆分！'.format(board_batch)}
                return json.dumps(info)
        elif os.path.exists(board_two):
            if not os.path.exists(board_two + "/RTAComplete.txt") or not os.path.isfile(board_two + "/RTAComplete.txt"):
                info = {'success': False, 'info': '板子{}还没有完全下机不能进行拆分！'.format(board_batch)}
                return json.dumps(info)
        else:
            info = {'success': False, 'info': '板子{}在nextseq与nextseq1中不存在！'.format(board_batch)}
            return json.dumps(info)
