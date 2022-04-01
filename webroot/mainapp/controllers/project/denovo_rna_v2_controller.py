# -*- coding: utf-8 -*-

from __future__ import print_function
import web
import json
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.denovo_rna_v2 import DenovoRnaV2
from mainapp.models.mongo.core.base import Base
from mainapp.config.db import Config
from biocluster.file import download
import os
import getpass


class DenovoRnaV2Controller(Base):
    def __init__(self, bind_object=None, instant=False, ):
        super(DenovoRnaV2Controller, self).__init__(bind_object)
        self._project_type = "denovo_rna_v2"
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        # 下面的denovo_rna成了实例对象，拥有多种和数据库交互的函数。
        self.denovo_rna_v2 = DenovoRnaV2(bind_object=bind_object)

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
        1. mbio.api.web.denovo_rna_v2.update_status
        2. mbio.api.web.denovo_rna_v2.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'denovo_rna_v2.update_status'
        else:
            return 'denovo_rna_v2.tupdate_status'

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
                self.denovo_rna_v2.update_status_failed(update_info[i], i)
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
                       module_type="workflow", params=None, to_file=None, new_task_id=None):
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
        if not new_task_id:
            new_task_id = self.denovo_rna_v2.get_new_id(task_id)
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
            'db_type': 'denovo_rna_v2',  # 需要根据自己连接的哪个数据库来修改实际的名字
            'options': options,  # 需要配置
            'task_id': task_id,
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
        task_info = self.denovo_rna_v2.get_task_info(task_id)
        config = Config()
        if config.RGW_ENABLE:
            target_dir = os.path.join(config.get_project_region_bucket(project_type="denovo_rna_v2"), 'files',
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
        if "/mnt/ilustre" in path:
            if os.path.exists(path):
                return path
            else:
                # 推测历史文件的buket
                return 's3://rerewrweset' + path.split('rerewrweset')[1]
        elif "s3://" in path:
            return path
        elif "s3nb://" in path:
            return path
        elif "s3nb1://" in path:
            return path
        else:
            return "s3://" + path

    def download_from_s3(self, from_file, to_path='download/', inter_dir='', cover=True):
        base_name = os.path.basename(from_file)
        to_file = os.path.join(inter_dir, base_name)
        if os.path.exists(to_file) and os.path.getsize(to_file) != 0:
            pass
        else:
            print('from {} to {}'.format(from_file, to_file))
            download(from_file, to_file)
        return to_file

    '''
    def download_from_s3(self, from_file, to_path="download/", inter_dir="", cover=True):
        """
        从s3对象存储下载数据到本地, 为了避免堵塞进程，此功能应该放置在流程最后执行。
        :param from_file: 需要下载的文件路径或文件路径, 必须是类似s3region://bucket/key写法。
        因为对象存储中没有文件夹的概念，需要下载文件夹必须使用"/"结尾，以明确表明下载的是文件夹
        :param to_path: 下载文件相对于当前工作目录的存放目录。
        当路径为"/"结尾时，表示下载文件存放在此文件夹下，否者为下载完整路径。
        当from_file为文件夹时，此参数也必须以"/"结尾。目录层级与下载的s3目录层级结构相同。
        默认情况下放置在当前模块工作目录的download目录下。
        :param cover: 对已存在的文件是否覆盖
        :return:
        """
        if re.match(r"^/|^\.\.", to_path):
            raise Exception("不能使用绝对路径或切换到其他目录!")
        if os.path.basename(to_path) == ".":
            raise Exception("目标文件不能叫\".\"!")
        target_dir = False
        if re.match(r"/$", to_path):
            target_dir = True
        self.s3transfer = S3TransferManager()
        self.s3transfer.base_path = inter_dir
        self.s3transfer.overwrite = cover
        m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", from_file)
        if not m:
            raise Exception("下载路径%s格式不正确!" % from_file)
        else:
            region = m.group(1)
            bucket_name = m.group(2)
            key_name = m.group(3)
            if re.match(r"/$", key_name):
                if not target_dir:
                    raise Exception("下载文件为文件夹时，源路径%s也必须为文件夹,以\"/\"结尾!" % to_path)
                conn = self.s3transfer.config.get_rgw_conn(region, bucket_name)
                bucket = Bucket(connection=conn, name=bucket_name)
                for key in bucket.list(prefix=key_name):
                    source = os.path.join(from_file, key.name)
                    target = os.path.join(target_dir, os.path.relpath(key.name, key_name))
                    self.s3transfer.add(source, target)
            else:
                if not target_dir:  # 处理已存在文件的情况
                    target = os.path.join(inter_dir, to_path)
                    if os.path.exists(target):
                        if cover:
                            if os.path.isdir(target):
                                shutil.rmtree(target)
                            else:
                                os.remove(target)
                        else:
                            raise Exception("目标文件夹%s已经存在!" % target)
                else:
                    target = os.path.join(inter_dir, to_path, os.path.basename(key_name))
                self.s3transfer.add(from_file, target)
        self.s3transfer.wait()
        '''
