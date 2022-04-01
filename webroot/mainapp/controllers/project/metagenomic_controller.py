# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'
# last modified BY guhaidong, 2017-10-31

import web
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.metagenomic import Metagenomic
from meta_controller import MetaController
import os
from biocluster.config import Config
# from mainapp.config.db import Config
import pandas as pd
import json

class MetagenomicController(MetaController):
    def __init__(self, instant=False):
        super(MetagenomicController, self).__init__(instant)
        self.metagenomic = Metagenomic()
        self.db_path = os.path.join(Config().SOFTWARE_DIR, "database/metagenome")
        #self.db_path = "/mnt/ilustre/users/sanger-dev/app/database/metagenome"
        self.level_file = os.path.join(self.db_path, "level_to_file.xls")
        self.level_name_file = os.path.join(self.db_path, "database_id_to_file.xls")

    def _update_status_api(self):
        """
        根据client决定接口api为metagenomic.update_status/metagenomic.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'metagenomic.update_status'
        else:
            return 'metagenomic.tupdate_status'

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
        except Exception, e:
            self.roll_back()
            return {"success": False, "info": "ERROR: %s" % e }

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
        new_task_id = self.get_new_id(task_id)
        self._sheet_data = {
            'id': new_task_id,
            'interaction': True,
            'name': name,  # 需要配置
            'type': module_type,  # 可以配置
            'client': self.data.client,
            'output': self._create_output_dir(task_id, main_table_name),
            'project_sn': project_sn,
            'IMPORT_REPORT_DATA': True,
            'UPDATE_STATUS_API': self._update_status_api(),
            #'db_type': '_metagenomic',  # 特殊用途，仅用于basic中判断是哪个数据库
            'db_type': 'metagenomic',
            'options': options  # 需要配置
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
        task_info = self.metagenomic.get_task_info(task_id)
        if Config().RGW_ENABLE:
            part_dir = str(task_info['member_id']) + \
                          '/' + str(task_info['project_sn']) + '/' + \
                          task_id + '/interaction_results/' + main_table_name
            sanger_prefix = Config().get_project_region_bucket(project_type="metagenomic")
            target_dir = os.path.join(sanger_prefix, 'files/' + part_dir + '/')  # 需要配置bucket
        else:
            client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
            if client == 'client01':
                target_dir = 'sanger'
            else:
                target_dir = 'tsanger'
            target_dir += ':rerewrweset/files/' + str(task_info['member_id']) + \
                          '/' + str(task_info['project_sn']) + '/' + \
                          task_id + '/interaction_results/' + main_table_name
        # target_dir += ':rerewrweset/files/' + part_dir + ';' + sanger_prefix + 'rerewrweset/files/' + part_dir + '/'  # 同时上传到rere和s3
        return target_dir

    def use_s3(self, path):
        if path.startswith("rere"):
            return "s3://" + path  # 防止tsg_31796任务路径出错
        else:
            return path
        '''
        if "/mnt/ilustre" in path:
            return path
        elif "s3://" in path:
            return path
        else:
            print "s3://" + path
            return "s3://" + path
        '''

    def level_id(self, level, type=1):
        level_id_table = pd.read_table(self.level_file, sep='\t', header=0)
        # level_id_table.set_index("level_id")
        # act_level = level_id_table.ix[mylevel,"file_level"]
        # act_level = level_id_table[level_id_table["level_id"== mylevel]]
        if type == 1:
            mylevel = int(level)
            act_level1 = level_id_table[level_id_table["level_id"] == mylevel]["file_level"]
            act_level = "".join(act_level1.tolist())
            if act_level1.empty:
                raise Exception("level id错误或没有输入类型type!")
        elif type == 2:
            mylevel = str(level)
            act_level1 = level_id_table[level_id_table["file_level"] == mylevel]["level_id"]
            act_level = "".join(act_level1.tolist())
            if act_level1.empty:
                raise Exception("level_name错误!")
            act_level = int(act_level)
        else:
            raise Exception("输入type只能为1或2!")
        return act_level

    def level_convert(self, level_name, level_id, type=1 ):
        # type 1 : id to name, type 2 : name to id
        samename = [1,2,3,4,5,6,7,8,11,14,15,16,17,24,31,32,33,34,38,39,46,44,45,51,53,59,60,61,62,72]
        level_id = int(level_id)
        if not level_id in samename:
            level_table = pd.read_table(self.level_name_file, sep='\t', header=0)
            if type == 1:
                act_level_name1 = level_table[level_table["value"] == level_name]["file_name"]
                act_level_name = "".join(act_level_name1.tolist())
                if act_level_name1.empty:
                    raise Exception("level id错误或没有输入类型type!")
            elif type == 2:
                act_level_name1 = \
                level_table[(level_table["file_name"] == level_name) & (level_table["level_id"] == level_id)][
                    "value"]
                act_level_name = "".join(act_level_name1.tolist())
                if act_level_name1.empty:
                    raise Exception("level_name错误!")
            else:
                raise Exception("输入type只能为1或2!")
        else:
            act_level_name = level_name
        if int(level_id) == 10:
            act_level_name = act_level_name[0]
        return act_level_name

    def check_level_with_database(self,anno_type,level_id,database=None):
        """
        :param anno_type: 页面的数据库名称
        :param level_id:
        :param database:
        :return:
        """
        pass
