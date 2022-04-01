# -*- coding: utf-8 -*-
# __author__ = 'xuting'

import web
import json
from mainapp.controllers.core.basic import Basic
import random
import re,os
from mainapp.libs.signature import check_sig
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.submit.sequence.sample_extract import SampleExtract as SE
from biocluster.config import Config
from biocluster.file import list_dir, exists
import datetime



class SampleExtract(object):
    """
    检测序列文件，获取序列中的样本信息，序列个类型可能是文件，也可能是文件夹; 序列的格式可能是fastq, 也可能是fasta
    """
    @check_sig
    def POST(self):
        """
        接口参数:
        file_path: 序列文件的路径
        task_id: 任务id
        type: 文件的类型， 值为file或者dir
        format: 序列的格式，值为fasta或者fastq
        """
        data = web.input()
        print data
        params = ["file_info", "format", "task_id", "client"]
        ## 上面这几个参数不一定是必传的，因为加入了一种新的情况：copy已有任务
        # for name in params:
        #     if not hasattr(data, name):
        #         variables = []
        #         variables.append(name)
        #         info = {"success": False, "info": "参数{}不存在".format(name), 'code':'C3000101', 'variables':variables}
        #         return json.dumps(info)
        if hasattr(data, "client"):
            if data.client not in ["client01", "client03"]:
                variables = []
                variables.append(data.client)
                info = {"success": False, "info": "未知的client：{}".format(data.client), 'code':'C3000102', 'variables':variables}
                return json.dumps(info)
            if data.client == "client01":
                pre_path = "sanger:"
                type_name = "sanger"
            elif data.client == "client03":
                pre_path = "tsanger:"
                type_name = "tsanger"
        if hasattr(data, "format"):
            if data.format not in ["sequence.fasta", "sequence.fastq", "sequence.fasta_dir", "sequence.fastq_dir"]:
                info = {"success": False, "info": "参数format的值不正确", 'code':'C3000103'}
                return json.dumps(info)
        if hasattr(data, "format"):
            if data.format in ["sequence,fasta_dir", "sequence.fastq_dir"] and not hasattr(data, "params_path"):
                info = {"success": False, "info": "params_path not exists!"}
                return json.dumps(info)

        file_info = json.loads(data.file_info)
        # if self.is_s3_path(file_info["path"]):
        #     if data.format in ["sequence.fasta", "sequence.fastq"] and not exists(file_info["path"]):
        #         variables = []
        #         variables.append(file_info["path"])
        #         info = {"success": False, "info": "文件{}不存在".format(file_info["path"]), 'code':'C3000104', 'variables':variables}
        #         return json.dumps(info)
        #     elif data.format in ["sequence.fasta_dir", "sequence.fastq_dir"] and not list_dir(file_info["path"]):
        #         variables = []
        #         variables.append(file_info["path"])
        #         info = {"success": False, "info": "文件夹{}不存在".format(file_info["path"]), 'code':'C3000105', 'variables':variables}
        #         return json.dumps(info)
        # else:
        #     suff_path = re.split(":", file_info["path"])[-1]
        #     config = Config().get_netdata_config(type_name)
        #     rel_path = os.path.join(config[type_name + "_path"], suff_path)
        #     print rel_path
        #     if data.format in ["sequence.fasta", "sequence.fastq"] and not exists(rel_path):
        #         variables = []
        #         variables.append(file_info["path"])
        #         print rel_path
        #         info = {"success": False, "info": "文件{}不存在".format(file_info["path"]), 'code':'C3000104', 'variables':variables}
        #         return json.dumps(info)
        #     elif data.format in ["sequence.fasta_dir", "sequence.fastq_dir"] and not list_dir(rel_path):
        #         variables = []
        #         variables.append(file_info["path"])
        #         info = {"success": False, "info": "文件夹{}不存在".format(file_info["path"]), 'code':'C3000105', 'variables':variables}
        #         return json.dumps(info)
        my_params = dict()
        if hasattr(data, "task_id"):
            my_params["task_id"] = data.task_id
        if hasattr(data, "file_info"):
            my_params["file_info"] = json.loads(data.file_info)
        if hasattr(data, "format"):
            my_params["format"] = data.format
        if hasattr(data, "query_id"):
            my_params["query_id"] = data.query_id
        if hasattr(data, "from_task"):
            my_params["from_task"] = data.from_task
        if hasattr(data, "file_name"):
            my_params["file_name"] = data.file_name
        params = json.dumps(my_params, sort_keys=True, separators=(',', ':'))
        main_id = SE().add_sg_seq_sample(data.task_id, data.file_info, params, data.query_id)
        if hasattr(data, "from_task"):
            table_id = SE().add_main_table(data.task_id, data.from_task, params, data.query_id)
        else:
            table_id = SE().add_sample_check(data.task_id, data.file_info, params, data.query_id)
        json_obj = dict()
        json_obj["name"] = "sequence.sample_extract"
        json_obj["id"] = self.get_new_id(data.task_id)
        # json_obj["id"] = data.task_id
        json_obj['type'] = "workflow"
        json_obj["IMPORT_REPORT_DATA"] = True
        json_obj["IMPORT_REPORT_AFTER_END"] = True
        json_obj["USE_DB"] = True
        json_obj['client'] = data.client
        json_obj["options"] = dict()
        if hasattr(data, "format"):
            if data.format == "sequence.fastq":
                if re.match(r"^([\w\-]+)://.*", file_info["path"]):
                    json_obj["options"]["in_fastq"] = "{}||{}".format(data.format, file_info["path"])
                else:
                    json_obj["options"]["in_fastq"] = "{}||{}{}".format(data.format, pre_path, file_info["path"])
            elif data.format == "sequence.fasta":
                if re.match(r"^([\w\-]+)://.*", file_info["path"]):
                    json_obj["options"]["in_fasta"] = "{}||{}".format(data.format, file_info["path"])
                else:
                    json_obj["options"]["in_fasta"] = "{}||{}{}".format(data.format, pre_path, file_info["path"])
            elif data.format == "sequence.fastq_dir":
                request_obj = {
                "format": data.format,
                "task_id": data.task_id,
                "client": data.client,
                "params_path": data.params_path,
                "mapping_file": True
                }
                json_obj["options"]["in_fastq"] = request_obj
            elif data.format == "sequence.fasta_dir":
                request_obj = {
                "format": data.format,
                "task_id": data.task_id,
                "client": data.client,
                "params_path": data.params_path,
                "mapping_file": True
                }
                json_obj["options"]["in_fasta"] = request_obj
            else:
                info = {"success": False, "info": "have unknown format"}
                print "unknown format %s" % data.format
                return json.dumps(info)

        if hasattr(data, "from_task"):
            json_obj["options"]["from_task"] = data.from_task
        if hasattr(data, "query_id"):
            json_obj["options"]["query_id"] = data.query_id
        if hasattr(data, "task_id"):
            json_obj["options"]["task_id"] = data.task_id
        if hasattr(data, "file_name"):
            json_obj["options"]["file_name"] = data.file_name
        # 给json文件的options字段指定输入的序列，这个option的字段可能是in_fastq或者是in_fasta
        # if self.is_s3_path(file_info['path']):
        #     if file_info["file_list"] in [None, "none", "None", "null", 'Null', '[]', '', []]:
        #         json_obj["options"]["in_fastq"] = "{}||{}".format(data.format, file_info["path"])
        #     else:
        #         json_obj["options"]["in_fastq"] = "{}||{}".format(data.format, json.dumps(file_info["file_list"]))
        # elif file_info["file_list"] in [None, "none", "None", "null", 'Null', '[]', '', []]:
        #     if data.format in ["sequence.fasta_dir", "sequence.fasta"]:
        #         json_obj["options"]["in_fasta"] = "{}||{}/{}".format(data.format, pre_path, suff_path)
        #     elif data.format in ["sequence.fastq_dir", "sequence.fastq"]:
        #         json_obj["options"]["in_fastq"] = "{}||{}/{}".format(data.format, pre_path, suff_path)
        # else:
        #     if data.format in ["sequence.fasta_dir", "sequence.fasta"]:
        #         json_obj["options"]["in_fasta"] = "{}||{}/{};;{}".format(data.format, pre_path, suff_path, file_info["file_list"])
        #     elif data.format in ["sequence.fastq_dir", "sequence.fastq"]:
        #         json_obj["options"]["in_fastq"] = "{}||{}/{};;{}".format(data.format, pre_path, suff_path, json.dumps(file_info["file_list"]))

        json_obj["options"]["table_id"] = str(table_id)
        json_obj["options"]["main_id"] = str(main_id)
        update_info = json.dumps({str(table_id): "sg_sample_check"})
        json_obj["options"]["update_info"] = update_info
        if data.client == "client01":
            update_api = "meta.update_status"
        elif data.client == "client03":
            update_api = "meta.tupdate_status"
        json_obj["UPDATE_STATUS_API"] = update_api
        """
        workflow_module = Workflow()

        workflow_module.add_record(insert_data)
        info = {"success": True, "info": "样本提取接口投递成功，正在计算中..."}
        return json.dumps(info)
        """
        print json_obj
        workflow_client = Basic(data=json_obj, instant= False)
        try:
            run_info = workflow_client.run()
            self._return_msg = workflow_client.return_msg
            return json.dumps(run_info)
        except Exception, e:
            variables = []
            variables.append(e)
            return json.dumps({"success": False, "info": "运行出错: %s" % e,  'code':'C3000106', 'variables':variables})

    def get_new_id(self, task_id):
        new_id = "{}_{}_{}".format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"), random.randint(1, 10000))
        workflow_module = Workflow()
        workflow_data = workflow_module.get_by_workflow_id(new_id)
        if len(workflow_data) > 0:
            return self.get_new_id(task_id)
        return new_id

    def is_s3_path(self, path):
        '''
        判断此路径是否是存储于对象存储中
        :param path:
        :return:
        '''
        if path.startswith("tsanger") or path.startswith("sanger") or path.startswith("rerewrweset"):
            return False
        else:
            return True
