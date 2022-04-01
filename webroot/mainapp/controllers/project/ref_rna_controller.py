# -*- coding: utf-8 -*-
# __author__ = 'qin.danhua'
import web
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.ref_rna import RefRna
from meta_controller import MetaController
from mainapp.config.db import Config
import os
import json
import shutil
from biocluster.file import getsize, exists
from biocluster.file import download
import glob
from boto.s3.bucket import Bucket
import re
from biocluster.file import exists
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.file import download
import getpass
#from mainapp.models.mongo.meta import RefRna


class RefRnaController(MetaController):
    def __init__(self, instant=False):
        super(RefRnaController, self).__init__(instant)
        self.meta = RefRna()
        self.ref_rna = RefRna()

    def _update_status_api(self):
        """
        根据client决定接口api为ref_rna.update_status/ref_rna.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'ref_rna.update_status'
        else:
            return 'ref_rna.tupdate_status'

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
            return run_info
        except Exception, e:
            self.roll_back()
            return {"success": False, "info": "运行出错: %s" % e }

    def set_sheet_data(self, name, options, main_table_name, task_id, project_sn, module_type="workflow", params=None, to_file=None):
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

        #table_info = self.ref_rna.get_main_info(main_id=main_id,collection_name=collection_name)
        #task_id = table_info["project_sn"]
        #task_id = table_info["task_id"]
        new_task_id = self.get_new_id(task_id)
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
            'db_type': 'ref_rna',  # 特殊用途，仅用于basic中判断是哪个数据库
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

    def get_main_info(self, main_id, collection_name):
        express_info = self.ref_rna.get_main_info(main_id=main_id, collection_name=collection_name)
        return express_info


    def _create_output_dir(self, task_id, main_table_name):
        data = web.input()
        task_info = self.meta.get_task_info(task_id)
        config = Config()
        if config.RGW_ENABLE:
            target_dir = os.path.join(config.get_project_region_bucket(project_type="ref_rna"), "files",
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
        # if path.startswith("/"):
        #     if os.path.exists(path):
        #         return path
        #     else:
        #         raise Exception("file not exists {}".format(path))
        #         # 推测历史文件的buket
        #         # return 's3://rerewrweset' + path.split('rerewrweset')[1]
        # elif re.match(r'^\w+://\S+/.+$', path):
        #     if exists(path) or path.endswith("/"):
        #         return path
        #     else:
        #         raise Exception("file not exists {}".format(path))
        # else:
        #     raise Exception("file format wrong".format(path))
        #     # return "s3://" + path

    def download_from_s3(self, from_file, to_path="download/", inter_dir="", cover=True):
        base_name = os.path.basename(from_file)
        to_file = os.path.join(inter_dir, base_name)
        #print('from {} to {}'.format(from_file, to_file))

        if os.path.exists(to_file) and os.path.getsize(to_file)!= 0:
            pass
        else:
            ## 接口不允许多线程下载，20181217
            # transfer = MultiFileTransfer()
            # print('from {} to {}'.format(from_file, to_file))
            # transfer.add_download(os.path.dirname(from_file) + "/" , os.path.dirname(to_file) + "/")
            # transfer.perform()
            download(from_file, to_file)
        return to_file
