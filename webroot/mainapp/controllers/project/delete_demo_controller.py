# -*- coding: utf-8 -*-
from __future__ import print_function
import web
import json
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.delete_demo import DeleteDemo
from mainapp.models.mongo.core.base import Base
from mainapp.config.db import Config
from boto.s3.bucket import Bucket
import re
from biocluster.file import exists


class DeleteDemoController(Base):
    def __init__(self, bind_object=None, instant=False, ):
        super(DeleteDemoController, self).__init__(bind_object)
        self._project_type = "project"
        self._instant = instant
        self._post_data = None
        self._sheet_data = None
        self._return_msg = None
        # 下面的ref_rna成了实例对象，拥有多种和数据库交互的函数。
        self.delete_demo = DeleteDemo(bind_object=bind_object)

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
                self.ref_rna_v2.update_status_failed(update_info[i], i)
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
            return run_info
        except Exception as e:
            self.roll_back()
            return {"success": False, "info": "运行出错: {}".format(e)}
