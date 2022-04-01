# -*- coding: utf-8 -*-
# __author__ = 'hesheng'
# last_modified guhaidong "add metagenomic target_member_type @ 20171201
import web
import json
from mainapp.libs.signature import check_sig
from mainapp.controllers.core.basic import Basic
from biocluster.core.function import filter_error_info
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.ref_rna import RefRna
from biocluster.config import Config
from mainapp.models.mongo.core.base import Base
from biocluster.wpm.client import *
import random
import datetime


class DemoMongodataCopy(Base):
    def __init__(self):
        super(DemoMongodataCopy, self).__init__()
        self._project_type = 'ref_rna'

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        requires = ['type', 'task_id', 'target_task_id', 'target_project_sn', 'target_member_id', 'target_member_type']
        for i in requires:
            if not (hasattr(data, i)):
                return json.dumps({"success": False, "info": "parameters missing:%s" % i})
        workflow_id = self.get_new_id(data.task_id)
        task_id1 = data.task_id
        if data.type == 'meta':
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }
            }
        elif data.type == "ref_rna":
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.refrna_copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id
                }
            }
            # mongodb = Config().mongo_client[Config().MONGODB + "_ref_rna"]
            # collection = mongodb['sg_task']
            collection = self.db['sg_task']
            nums1 = collection.count({"task_id": {"$regex": data.task_id}})
            nums = collection.count({"task_id": {"$regex": data.task_id}, "demo_status": "end"})
            if nums1:
                if nums1 <= 3:
                    id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
                    json_obj = {
                        "type": "workflow",
                        "id": id,
                        "name": "copy_demo.demo_backup",
                        "IMPORT_REPORT_DATA": True,
                        "IMPORT_REPORT_AFTER_END": False,
                        "options": {
                            "task_id": data.task_id,
                            "target_task_id": data.task_id + '_' + id,
                            "target_project_sn": "refrna_demo",
                            "target_member_id": "refrna_demo"
                        }
                    }
                    worker = worker_client()
                    worker.add_task(json_obj)
                if nums < 3:
                    info = {"success": False, "info": "demo数据正在准备中，请一段时间后再次进行拉取", 'code':'C2200801'}
                    return json.dumps(info)
            else:
                info = {"success": False, "info": "demo数据有问题，不能进行拉取", 'code':'C2200802'}
                return json.dumps(info)  # add by hongdongxuan 20171022
        elif data.type == "metagenomic":  # modified by guhaidong 20171215
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.metagenomic_copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }  # 增加target_member_type @ 20171201
            }
        elif data.type == "bacgenome":  # modified bygaohao 20180426
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.bac_copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }  # 增加target_member_type @ 20171201
            }
        elif data.type == "fungigenome":  # modified bygaohao 20180612
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.fungi_copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }  # 增加target_member_type @ 20171201
            }
        elif data.type == "metabolome":  # modified by shaohua.yuan 20180720
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.metabolome_copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }
            }
        elif data.type == "metagbin":  # modified by qingchen.zhang 20190104
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.metagbin_copy_demo",  # 需要配置
                'type': 'workflow',  # 可以配置
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }  # 增加target_member_type @ 20171201
            }
        elif data.type == "bac_assem":  # modified by haidong.gu 20190528
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.bacassem_copy_demo",
                'type': 'workflow',
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }
            }

        elif data.type == 'meta_asv':  #modify by qingchen.zhang @ 20200528
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.metaasv_copy_demo",
                'type': 'workflow',
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }
            }
        elif data.type == 'bac_comparative':  #modify by gao.hao @ 20200915
            data_json = {
                'id': workflow_id,
                'stage_id': 0,
                'name': "copy_demo.bac_comparative_demo",
                'type': 'workflow',
                'client': data.client,
                'project_sn': data.target_project_sn,
                'options': {
                    "task_id": data.task_id,
                    "target_task_id": data.target_task_id,
                    "target_project_sn": data.target_project_sn,
                    "target_member_id": data.target_member_id,
                    "target_member_type": data.target_member_type
                }
            }
        else:
            variables = []
            variables.append(data.type)
            info = {"success": False, "info": "项目类型type：{}不合法!".format(data.type), 'code':'C2200803', 'variables':variables}
            return json.dumps(info)
        workflow_client = Basic(data=data_json, instant=True)
        try:
            run_info = workflow_client.run()
            run_info['info'] = filter_error_info(run_info['info'])
            if data.type == "ref_rna": ## last modified by shicaiping 20171116
                time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                col = self.db["sg_task"]
                result = col.find_one({"task_id": {"$regex": task_id1 + "_.*_.*"}, "demo_status": "end"})
                demo_task_id = result["task_id"]
                run_info['demo_task_id'] = demo_task_id ## 返回给前端的run_info信息增加demo task id字段
                col.update({"task_id" : demo_task_id},{"$set":{'demo_status':"occupy"}}, upsert=False) ## 更新demo_status状态，防止重复拉取
                col.update({"task_id" : demo_task_id}, {"$set": {"is_demo": 2}})## 更新is_demo状态，用于获取序列信息
                col.update({"task_id" : demo_task_id}, {"$set": {"created_ts": time}})
                col.update({"task_id" : demo_task_id}, {"$set": {"member_id": data.target_member_id}}) ## 更新member_id
                col.update({"task_id" : demo_task_id}, {"$set": {"project_sn": data.target_project_sn}})## 更新project_sn
                col.update({"task_id" : demo_task_id}, {"$set": {"member_type": data.target_member_type}})## 更新member_type
            return json.dumps(run_info)
        except Exception as e:
            return json.dumps({"success": False, "info": "运行出错: %s" % filter_error_info(str(e))})

    def get_new_id(self, task_id, otu_id=None):
        """
        根据旧的ID生成新的workflowID，固定为旧的后面用“_”，添加两次随机数或者一次otu_id一次随机数
        """
        if otu_id:
            new_id = "{}_{}_{}".format(task_id, otu_id[-4:], random.randint(1, 10000))
        else:
            new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"),
                                       random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id, otu_id)
        return new_id
