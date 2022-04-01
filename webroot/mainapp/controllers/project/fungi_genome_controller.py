# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'
# modified by zouxuan
import web
from ..core.basic import Basic
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.fungi_genome import FungiGenome
from meta_controller import MetaController
from biocluster.config import Config
import json

class FungiGenomeController(MetaController):
    def __init__(self, instant=False):
        super(FungiGenomeController, self).__init__(instant)
        self._instant = instant
        self._post_data = None
        self._return_msg = None
        self._sheet_data = None
        self.fungi_genome = FungiGenome()

    def _update_status_api(self):
        """
        根据client决定接口api为fungi_genome.update_status/fungi_genome.tupdate_status
        """
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            return 'fungi_genome.update_status'
        else:
            return 'fungi_genome.tupdate_status'

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
            return {"success": False, "info": "ERROR: %s" % e}

    def set_sheet_data(self, name, options, main_table_name, task_id, project_sn, module_type="workflow", params=None,to_file=None):
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

        # table_info = self.ref_rna.get_main_info(main_id=main_id,collection_name=collection_name)
        # task_id = table_info["project_sn"]
        # task_id = table_info["task_id"]
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
            'db_type': 'fungigenome',  # 特殊用途，仅用于basic中判断是哪个数据库
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

    def get_main_info(self, main_id, collection_name):
        express_info = self.ref_rna.get_main_info(main_id=main_id, collection_name=collection_name)
        return express_info

    def _create_output_dir(self, task_id, main_table_name):
        data = web.input()
        task_info = self.fungi_genome.get_task_info(task_id)
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            target_dir = 'sanger'
        else:
            target_dir = 'tsanger'
        target_dir = Config().get_project_region_bucket(project_type="fungi")  #guanqing.zou 20180928
        target_dir += 'files/' + str(task_info['member_id']) + '/' + str(task_info['project_sn']) + '/' + task_id \
                      + '/interaction_results/' + main_table_name
        return target_dir