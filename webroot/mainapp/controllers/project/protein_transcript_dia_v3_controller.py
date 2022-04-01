# -*- coding: utf-8 -*-
from __future__ import print_function
import web
import json
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.itraq_and_tmt import ItraqTmt
from mainapp.models.mongo.dia import Dia
from mainapp.models.mongo.denovo_rna_v2 import DenovoRnaV2
from mainapp.models.mongo.ref_rna import RefRna
from mainapp.models.mongo.ref_rna_v2 import RefRnaV2
from mainapp.models.mongo.prok_rna import ProkRNA
from mainapp.models.mongo.core.base import Base
from mainapp.config.db import Config
from biocluster.file import download
from biocluster.file import getsize, exists
import os
import getpass


class ProteinTranscriptDiaV3Controller(Base):
    def __init__(self, bind_object=None, instant=False, ):
        super(ProteinTranscriptDiaV3Controller, self).__init__(bind_object)
        self._project_type = "dia"
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        self.itraq_tmt = ItraqTmt(bind_object=bind_object)
        self.dia = Dia(bind_object=bind_object)
        self.ref_rna_v2 = RefRnaV2(bind_object=bind_object)
        self.ref_rna = RefRna(bind_object=bind_object)
        self.denovo_rna_v2 = DenovoRnaV2(bind_object=bind_object)
        self.prok_rna = ProkRNA(bind_object=bind_object)

    @property
    def data(self):
        return self._post_data

    @property
    def return_msg(self):
        """
        获取Post数据
        """
        return self._return_msg

    @property
    def instant(self):
        """
        任务是否是即时计算
        """
        return self._instant

    @property
    def sheet_data(self):
        """
        获取运行流程所需的Json数据
        :return:
        """
        return self._sheet_data

    def _update_status_api(self):
        """
        根据client决定接口api为
        1. mbio.api.web.labelfree.update_status
        2. mbio.api.web.labelfree.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'dia.update_status'
        else:
            return 'dia.tupdate_status'

    def roll_back(self):
        """
        当任务投递失败时，如WPM服务出错时，主表写入start状态无法由API更新，此处进行更新

        :return:
        """
        print("INFO: 任务提交出错，尝试更新主表状态为failed。")
        try:
            update_info = json.loads(self.sheet_data['options']['update_info'])
            for i in update_info:
                if i == "batch_id":
                    continue
                self.dia.update_status_failed(update_info[i], i)
                print("INFO: 更新主表状态为failed成功: coll:{} _id:{}".format(update_info[i], i))
        except Exception as e:
            print('ERROR:尝试回滚主表状态为failed 失败:{}'.format(e))

    @check_sig
    def POST(self):
        workflow_client = Basic(data=self.sheet_data, instant=self.instant)
        try:
            run_info = workflow_client.run()
            print("打印出run_info： {}".format(run_info))
            try:
                info_msg = json.loads(run_info['info'])  # info嵌套info的情况
                if info_msg['code'] == 'R001' and info_msg['info'] != "":
                    del info_msg['code']
                run_info.update(info_msg)
            except:
                run_info['info'] = run_info['info']  # info只包含字符的情况
            self._return_msg = workflow_client.return_msg
            print("修改后的run_info： {}".format(run_info))
            if not run_info['success']:
                self.roll_back()
                return {"success": False, "info": "Failed to submit your task because too many tasks are being "
                                                  "calculated on the cluster. Please try again!"}
            else:
                return run_info
        except Exception as e:
            self.roll_back()
            return {"success": False, "info": "Failed to submit your task because too many tasks are being "
                                              "calculated on the cluster. Please try again!"}

    def set_sheet_data(self, name, options, main_table_name, task_id, project_sn,
                       module_type="workflow", params=None, to_file=None):
        """
        设置运行所需的Json文档

        :param name: workflow/module/tool相对路径
        :param module_type: workflow/module/tool
        :param main_table_name: 交互分析项主表名称
        :param options: workflow/module/tool参数
        :param params: 交互分析主表params字段
        :param to_file: workflow/module/tool mongo数据转文件
        :return:
        """
        self._post_data = web.input()
        new_task_id = self.dia.get_new_id(task_id)
        self._sheet_data = {
            'id': new_task_id,
            'stage_id': 0,
            'interaction': True,
            'name': name,  # 需要配置
            'type': module_type,  # 可以配置
            'client': self.data.client,
            'output': self._create_output_dir(task_id, main_table_name),
            'project_sn': project_sn,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(),
            'db_type': 'dia',  # 需要根据自己连接的哪个数据库来修改实际的名字
            'options': options,  # 需要配置
            'CLUSTER': getpass.getuser()
        }
        if self.instant:
            self._sheet_data["instant"] = True
        if params:
            self._sheet_data["params"] = params
        if to_file:
            self._sheet_data["to_file"] = to_file
        print('Sheet_Data: {}'.format(self._sheet_data))
        self.workflow_id = new_task_id
        return self._sheet_data

    def _create_output_dir(self, task_id, main_table_name):
        data = web.input()
        task_info = self.dia.get_task_info(task_id)
        config = Config()
        if config.RGW_ENABLE:
            target_dir = os.path.join(config.get_project_region_bucket(project_type="dia"), "files",
                                      str(task_info['member_id']), str(task_info['project_sn']),
                                      task_id, 'interaction_results', main_table_name)
        else:
            client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
            if client == 'client01':
                target_dir = 'sanger:'
            else:
                target_dir = 'tsanger:'

            target_dir += 'files/' + str(task_info['member_id']) + \
                          '/' + str(task_info['project_sn']) + '/' + \
                          task_id + '/interaction_results/' + main_table_name
        return target_dir

    def create_tmp_dir(self, task_id, main_table_name):
        tmp_dir = os.path.join(Config().get_work_dir(), "tmp", "tmp_s3", task_id, main_table_name)
        if os.path.exists(tmp_dir):
            pass
        else:
            os.makedirs(tmp_dir)
        return tmp_dir

    def use_s3(self, path):
        return path
        # if "/mnt/ilustre" in path:
        #     if os.path.exists(path):
        #         return path
        #     else:
        #         # 推测历史文件的buket
        #         return 's3://rerewrweset' + path.split('rerewrweset')[1]
        # elif "s3://" in path:
        #     return path
        # else:
        #     return "s3://" + path

    def download_from_s3(self, from_file, to_path="download/", inter_dir="", cover=True):
        base_name = os.path.basename(from_file)
        to_file = os.path.join(inter_dir, base_name)
        # print('from {} to {}'.format(from_file, to_file))

        if os.path.exists(to_file) and os.path.getsize(to_file) != 0:
            pass
        else:
            print('from {} to {}'.format(from_file, to_file))
            download(from_file, to_file)
        return to_file
