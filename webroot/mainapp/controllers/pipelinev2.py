# -*- coding: utf-8 -*-
# __author__ = 'HD'
# last modified by hd@20200818 (传入数据由xml变为json)
import web
# from biocluster.config import Config
# from mainapp.models.mongo.core.base import Base
from mainapp.libs.signature import check_sig, CreateSignature, check_sig_v2
from web import form
from mainapp.libs.input_check import check_format, check_format_json
import json
# from mainapp.models.workflow import Workflow
import os
# from mainapp.libs.jsonencode import CJsonEncoder
# import xml.etree.ElementTree as ET
from mainapp.config.db import Config
from biocluster.wpm.client import worker_client
# import traceback
import re
# from mainapp.libs.getip import get_ip
# from mainapp.models.admin.userlog import UserlogModel


class Pipelinev2(object):

    def GET(self):
        path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '../views/'))
        render = web.template.render(path)
        return render.pipline(self.get_form())

    @check_sig_v2
    @check_format_json
    def POST(self):
        data = web.input()
        # print data
        client = data.client if hasattr(
            data, "client") else web.ctx.env.get('HTTP_CLIENT')
        # dbversion = int(data.db_version) if hasattr(data, "db_version") else 0
        dbversion = 0
        if hasattr(data, "db_version"):
            try:
                dbversion = int(data.db_version)
            except Exception, e:
                return json.dumps({"success": False, "info": "db_version传入的字段不合法:{}".format(str(e))})
        dydb = 0
        if hasattr(data, "dydb"):
            try:
                dydb = int(data.dydb)
            except Exception, e:
                return json.dumps({"success": False, "info": "dydb传入的字段不合法:{}".format(str(e))})
        contract_id = ""
        if hasattr(data, "contract_id"):
            contract_id = data.contract_id
        dydb_contract_ids = {}
        if hasattr(data, "dydb_contract_ids"):
            try:
                dydb_contract_ids = json.loads(data.dydb_contract_ids)
            except Exception as e:
                print(e)
                return json.dumps({"success": False, "info": "dydb_contract_ids传入的字段不合法:{}".format(str(e))})
            # dydb_contract_ids = data.dydb_contract_ids
        api = "sanger"
        try:
            if client == "client01" or client == "client03":
                if client != "client01":
                    api = "tsanger"
                json_obj = self.sanger_submit_json()
            else:
                json_obj = self.json_submit()
        except Exception, e:
            return json.dumps({"success": False, "info": str(e)})
        if "type" not in json_obj.keys() or "id" not in json_obj.keys():
            info = {"success": False, "info": "Json内容不正确!!"}
            return json.dumps(info)
        json_obj["UPDATE_STATUS_API"] = api
        json_obj['client'] = client
        json_obj['interaction'] = False
        if dbversion != 0:
            json_obj["DBVersion"] = dbversion
        # if dydb != 0:
        json_obj["dydb"] = dydb
        if dydb ==1:
            if contract_id =="":
                info = {"success": False, "info": "参数dydb为1时必须设置contract_id"}
                return json.dumps(info)
            json_obj['contract_id']=contract_id
        if dydb_contract_ids:
            json_obj['dydb_contract_ids']=dydb_contract_ids
        response = worker_client().add_task(json_obj)
        if "success" in response.keys() and response["success"]:
            info = {"success": True, "info": "添加任务到队列成功!", "workflow_id": json_obj['id']}
            return json.dumps(info)
        else:
            info = {"success": False, "info": "添加任务到队列失败：{}!".format(response["info"]),
                    "workflow_id": json_obj['id']}
            return json.dumps(info)

    @staticmethod
    def sanger_submit_json():
        """
        前端传入的数据变为json格式，{"params" {}, "basis": {}}
        params中存入的数据是分析需要的参数；basis中传入的是工作流程需要的参数。
        {
            params:{
                qc："true",
                min_contig:"1000",
                sofware_bin:"metabat",
                cdhit_coverage:"0.9",
                in_fastq:{"format":"sequence.fastq_dir","alias":"rawdata","dir":"true"},
                specimen_info:{"format":"meta_genomic.specimen_info","alias":"specimen_info.txt",
                "path":"s3nb://commonbucket/files/m_6192/6192_5ccffbfa64abc/3sample/b38f5e9edeec033ffc.txt"},
                rm_host:"true",
                ref_database:"Vertebrates,Homo_sapiens",
                assemble_dir:{"format":"metagbin.assemble_dir","alias":"assemble","dir":"true"},
                cdhit_identity:"0.95",
                mapping_file:{"format":"mapping_file","alias":"mapping_file.txt",
                "path":"s3nb://commonbucket/files/m_6873/6873_5f279482aa8bb/mapping_file_1596429581.txt"}
            },
            basis:{
                member_id:"m_6873",
                cmd_id: 200,
                project_sn:"6873_5f279482aa8bb",
                task_id:"majorbio_274215",
                name:"metagbin.metagbin",
                module_type:"workflow",
                member_type:"1",
                type:"sanger",
                id : 'cmd_200_1596429581'
            }
        }
        add by hd@20200818
        :return:
        """
        data = web.input()
        json_obj = {}
        basis = json.loads(data.basis)
        # print basis
        params = json.loads(data.params)
        # print params
        # print "7777:{}".format(basis["member_id"])
        client = data.client if hasattr(data, "client") else None
        mapping_file = None
        file_path_type = "sanger" if client == "client01" else "tsanger"
        if "member_id" in basis.keys():
            json_obj["member_id"] = str(basis["member_id"])
        if "cmd_id" in basis.keys():
            json_obj["cmd_id"] = int(basis["cmd_id"])
        if "project_sn" in basis.keys():
            json_obj["project_sn"] = basis["project_sn"]
        if "name" in basis.keys():
            json_obj["name"] = basis["name"]
        if "task_id" in basis.keys():
            json_obj["id"] = basis["task_id"]
        if "member_type" in basis.keys():
            json_obj["member_type"] = int(basis["member_type"])
        else:
            json_obj['member_type'] = 0
        if "id" in basis.keys():
            json_obj["stage_id"] = basis["id"]
        else:
            json_obj['stage_id'] = "0"
        if "IMPORT_REPORT_DATA" in basis.keys() and basis["IMPORT_REPORT_DATA"] in ["false", "False", "FALSE", "no",
                                                                                    "No", "NO", "N"]:
            json_obj["IMPORT_REPORT_DATA"] = False  # 更新报告数据
        else:
            json_obj["IMPORT_REPORT_DATA"] = True  # 更新报告数据
        if "IMPORT_REPORT_AFTER_END" in basis.keys() and \
                basis["IMPORT_REPORT_AFTER_END"] in ["false", "False", "FALSE", "no", "No", "NO", "N"]:
            json_obj["IMPORT_REPORT_AFTER_END"] = False
        else:
            json_obj["IMPORT_REPORT_AFTER_END"] = True
        if "module_type" in basis.keys():
            json_obj['type'] = basis["module_type"]
        if "mapping_file" in params.keys():
            mapping_file = params["mapping_file"]["path"]
        config = Config()
        out_file_path = config.get_bucket_from_path(json_obj['name'])
        json_obj['output'] = os.path.join(out_file_path, "files", json_obj["member_id"], json_obj['project_sn'],
                                          json_obj['id'], 'workflow_results/')
        json_obj['options'] = {}
        # print "1111:{}".format(json_obj)
        for key in params:
            # print "333:{}".format(key)
            if key == "mapping_file":
                continue
            if isinstance(params[key], dict):  # 判断该参数值file格式
                print "test is ok"
                if "dir" in params[key].keys() and params[key]['dir']:
                    if not mapping_file:
                        raise Exception("文件夹参数未找到mapping file!")
                    else:
                        json_obj['options'][key] = "filelist[%s](%s):%s" % (key, file_path_type, mapping_file)
                else:
                    # print str(params[key])
                    if re.match(r"^([\w\-]+)://.*", str(params[key]['path'])):
                        # print "s3 file go this way"
                        json_obj['options'][key] = "%s" % params[key]['path']
                    else:
                        # print "not s3 path go this way"
                        json_obj['options'][key] = "%s:%s" % (file_path_type, params[key]['path'])
                if "alias" in params[key].keys() and params[key]['alias']:
                    json_obj['options'][key] += "{%s}" % params[key]['alias']
                if "format" in params[key].keys() and params[key]["format"]:
                    json_obj['options'][key] = "%s||%s" % (params[key]["format"], json_obj['options'][key])
            else:
                json_obj['options'][key] = params[key]
                # print "test is ok2"
        # print "222:{}".format(json_obj)
        return json_obj

    @staticmethod
    def json_submit():
        data = web.input()
        json_obj = json.loads(data.json)
        return json_obj

    @staticmethod
    def get_form():
        sig_obj = CreateSignature("test")
        return form.Form(
            form.Hidden(name='client', value=sig_obj.client),
            form.Hidden(name='nonce', value=sig_obj.nonce),
            form.Hidden(name='timestamp', value=sig_obj.timestamp),
            form.Hidden(name='signature', value=sig_obj.signature),
            form.Textarea("json", description="Json", rows="20", cols="100"),
            form.Dropdown(name='cluster', args=['sanger', 'isanger'], value='isanger'),
            form.Button("submit", type="submit", description="提交")
        )


class Pipelinev2Stop(object):

    @check_sig_v2
    def POST(self):
        data = web.input()
        if not hasattr(data, "id") or data.id.strip() == "":
            info = {"success": False, "info": "停止任务需要参数: id!"}
            return json.dumps(info)
        json_data = {"id": data.id.strip(), "msg": "stop"}
        response = worker_client().run_cmd(json_data)
        return json.dumps(response)


class Pipelinev2Pause(object):

    @check_sig_v2
    def POST(self):
        data = web.input()
        print data
        if not hasattr(data, "id") or data.id.strip() == "":
            info = {"success": False, "info": "暂停任务需要参数: id!"}
            return json.dumps(info)
        json_data = {"id": data.id.strip(), "msg": "pause"}
        response = worker_client().run_cmd(json_data)
        return json.dumps(response)


class Pipelinev2StopPause(object):

    def POST(self):
        return self.GET()

    @check_sig_v2
    def GET(self):
        data = web.input()
        print data
        if not hasattr(data, "id") or data.id.strip() == "":
            info = {"success": False, "info": "重新开始暂停任务需要参数 id!"}
            return json.dumps(info)
        json_data = {"id": data.id.strip(), "msg": "continue"}
        response = worker_client().run_cmd(json_data)
        return json.dumps(response)
