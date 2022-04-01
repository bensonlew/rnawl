# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
# last modified by liulinmeng (add member_type)
import web
# from biocluster.config import Config
# from mainapp.models.mongo.core.base import Base
from mainapp.libs.signature import check_sig, CreateSignature
from web import form
from mainapp.libs.input_check import check_format
import json
# from mainapp.models.workflow import Workflow
import os
# from mainapp.libs.jsonencode import CJsonEncoder
import xml.etree.ElementTree as ET
# from mainapp.config.db import get_use_api_clients, get_api_type, get_mongo_client, Config
from mainapp.config.db import Config
from biocluster.wpm.client import worker_client
# import traceback
import re
# from mainapp.libs.getip import get_ip
# from mainapp.models.admin.userlog import UserlogModel


class Pipeline(object):

    def GET(self):
        path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '../views/'))
        render = web.template.render(path)
        return render.pipline(self.get_form())

    @check_sig
    @check_format
    def POST(self):
        data = web.input()
        client = data.client if hasattr(
            data, "client") else web.ctx.env.get('HTTP_CLIENT')
        dbversion = 0
        if hasattr(data, "db_version"):
            try:
                dbversion = int(data.db_version)
            except Exception, e:
                return json.dumps({"success": False, "info": "db_version传入的字段不合法:{}".format(str(e))})
        api = "sanger"
        try:
            if client == "client01" or client == "client03":
                if client != "client01":
                    api = "tsanger"
                json_obj = self.sanger_submit()
            else:
                json_obj = self.json_submit()
        except Exception, e:
            return json.dumps({"success": False, "info": str(e)})
        if "type" not in json_obj.keys() or "id" not in json_obj.keys():
            info = {"success": False, "info": "Json内容不正确!!"}
            return json.dumps(info)
        # if client in get_use_api_clients():
        #     api = get_api_type(client)
        #     if api:
        #         json_obj["UPDATE_STATUS_API"] = api
        json_obj["UPDATE_STATUS_API"] = api
        json_obj['client'] = client
        json_obj['interaction'] = False
        if dbversion != 0:
            json_obj['DBVersion'] = dbversion
        response = worker_client().add_task(json_obj)
        # session = web.config.get('_session')
        # if session and session.is_login:
        #     model = UserlogModel()
        #     model.add("添加Workflow: %s" % json_obj['id'])
        if "success" in response.keys() and response["success"]:
            info = {"success": True, "info": "添加任务到队列成功!", "workflow_id": json_obj['id']}
            return json.dumps(info)
        else:
            info = {"success": False, "info": "添加任务到队列失败：{}!".format(response["info"]),
                    "workflow_id": json_obj['id']}
            return json.dumps(info)

    @staticmethod
    def sanger_submit():
        data = web.input()
        xml_data = "".join(data.content)
        root = ET.fromstring(xml_data)
        json_obj = {}
        client = None
        new_version2 = False
        mapping_file = None
        if hasattr(data, "client"):
            client = data.client
        if client == "client01":
            file_path_type = "sanger"
        else:
            file_path_type = "tsanger"
        for child_of_root in root:
            if child_of_root.tag == "member_id":
                json_obj["member_id"] = child_of_root.text
            if child_of_root.tag == "cmd_id":
                json_obj["cmd_id"] = int(child_of_root.text)
            if child_of_root.tag == "project_sn":
                json_obj['project_sn'] = child_of_root.text
            if child_of_root.tag == "name":
                json_obj['name'] = child_of_root.text
            if child_of_root.tag == "task_id":
                json_obj['id'] = child_of_root.text
            if child_of_root.tag == "version" and int(child_of_root.text) == 2:
                new_version2 = True
            if child_of_root.tag == "member_type":
                try:
                    json_obj['member_type'] = int(child_of_root.text)
                except Exception:
                    json_obj['member_type'] = 0
            if child_of_root.tag == "mapping_file":
                mapping_file = child_of_root.text
            # if child_of_root.tag == "bucket":
            #     json_obj["bucket"] = child_of_root.text
            #     file_path = file_path_type + child_of_root.text
            # if child_of_root.tag == "region":
            #     json_obj['region'] = child_of_root.text
            if child_of_root.tag == "IMPORT_REPORT_DATA" and child_of_root.text in ["false", "False", "FALSE",
                                                                                    "no", "No", "NO", "N"]:
                json_obj["IMPORT_REPORT_DATA"] = False  # 更新报告数据
            else:
                json_obj["IMPORT_REPORT_DATA"] = True
            if child_of_root.tag == "IMPORT_REPORT_AFTER_END" and child_of_root.text in ["false", "False", "FALSE",
                                                                                         "no", "No", "NO", "N"]:
                json_obj["IMPORT_REPORT_AFTER_END"] = False
            else:
                json_obj["IMPORT_REPORT_AFTER_END"] = True
        first_stage = root.find("stage")
        json_obj['stage_id'] = first_stage.find("id").text
        if not new_version2:
            raise Exception("不兼容此接口版本!")
        if not json_obj['stage_id']:
            json_obj['stage_id'] = 0
        json_obj['type'] = first_stage.find("type").text
        json_obj['name'] = first_stage.find("name").text
        # json_obj['output'] = "%s/files/%s/%s/%s/%s" % (file_path, json_obj["member_id"], json_obj['project_sn'],
        #                                                json_obj['id'], json_obj['stage_id'])
        config = Config()
        out_file_path = config.get_bucket_from_path(json_obj['name'])
        json_obj['output'] = os.path.join(out_file_path, "files", json_obj["member_id"], json_obj['project_sn'],
                                          json_obj['id'], 'workflow_results/')
        option = first_stage.find("parameters")
        # print json_obj
        json_obj['options'] = {}
        for opt in option:
            if 'type' in opt.attrib.keys():
                if "dir" in opt.attrib.keys() and opt.attrib["dir"]:
                    if not mapping_file:
                        raise Exception("文件夹参数未找到mapping file!")
                    else:
                        json_obj['options'][opt.tag] = "filelist[%s](%s):%s" % (
                            opt.tag, file_path_type, mapping_file)
                else:
                    if re.match(r"^([\w\-]+)://.*", opt.text):
                        json_obj['options'][opt.tag] = "%s" % opt.text
                    else:
                        json_obj['options'][opt.tag] = "%s:%s" % (file_path_type, opt.text)
                if "alias" in opt.attrib.keys() and opt.attrib["alias"]:
                    json_obj['options'][opt.tag] += "{%s}" % opt.attrib["alias"]
                if "format" in opt.attrib.keys() and opt.attrib["format"]:
                    json_obj['options'][opt.tag] = "%s||%s" % (opt.attrib["format"], json_obj['options'][opt.tag])
            else:
                json_obj['options'][opt.tag] = opt.text
            # if 'fileList_file' in opt.attrib.keys():
            #     path = "filelist:%s" % opt.text
            #     if "format" in opt.attrib.keys():
            # else:
            #     if 'type' in opt.attrib.keys() and opt.attrib["type"]:
            #
            #             file_list = ''
            #             if "fileList" in opt.attrib.keys():
            #                 file_list = opt.attrib['fileList']
            #             tmp_list = [None, "none", "None",
            #                         "null", 'Null', '[]', '']
            #             if new_version2:
            #                 if file_list in tmp_list:
            #                     if m:
            #                         json_obj['options'][
            #                             opt.tag] = "%s||%s" % (opt.attrib["format"], opt.text)
            #                     else:
            #                         json_obj['options'][
            #                             opt.tag] = os.path.join("%s||%s" % (opt.attrib["format"], file_path_type), opt.text)
            #                     if "alias" in opt.attrib.keys() and opt.attrib['alias']:
            #                         json_obj['options'][opt.tag] += "{%s}" % opt.attrib['alias']
            #                 else:
            #                     new_file_list = []
            #                     _file_list_obj = json.loads(file_list)
            #                     for f in _file_list_obj:
            #                         if not re.match(r"^([\w\-]+)://.*", f["file_path"]):
            #                             f["file_path"] = "%s:%s" % (file_path_type, f["file_path"])
            #                         new_file_list.append(f)
            #                     json_obj['options'][opt.tag] = "%s||%s" % (opt.attrib["format"],
            #                                                                json.dumps(new_file_list))
            #             else:
            #                 if file_list in tmp_list:
            #                     if m:
            #                         json_obj['options'][
            #                             opt.tag] = "%s||%s" % (opt.attrib["format"], opt.text)
            #                     else:
            #
            #                         json_obj['options'][
            #                             opt.tag] = os.path.join("%s||%s" % (opt.attrib["format"], file_path), opt.text)
            #                 else:
            #                     # json_obj['options'][opt.tag] = "{}||{}/{};;{}".format(opt.attrib["format"], file_path, opt.text, file_list)
            #                     if m:
            #                         path = "{}||{}".format(opt.attrib["format"], opt.text)
            #                         json_obj['options'][opt.tag] = path + ";;{}".format(
            #                             file_list)
            #                     else:
            #                         path = "{}||{}".format(opt.attrib["format"], file_path)
            #                         json_obj['options'][opt.tag] = os.path.join(path, opt.text) + ";;{}".format(
            #                             file_list)
            #         else:
            #             # json_obj['options'][opt.tag] = "%s/%s" % (file_path, opt.text)
            #             if m:
            #                 json_obj['options'][opt.tag] = opt.text
            #             else:
            #                 if new_version2:
            #                     json_obj['options'][opt.tag] = file_path_type + opt.text
            #                 else:
            #                     json_obj['options'][opt.tag] = os.path.join(file_path, opt.text)
            #             if "alias" in opt.attrib.keys() and opt.attrib['alias']:
            #                 json_obj['options'][opt.tag] += "{%s}" % opt.attrib['alias']
            #         # else:
            #         #     if "format" in opt.attrib.keys():
            #         #         if m:
            #         #             json_obj['options'][opt.tag] = "%s||%s" % (opt.attrib["format"], opt.text)
            #         #         else:
            #         #             json_obj['options'][opt.tag] = "%s||%s:%s" % (opt.attrib["format"], opt.attrib["type"], opt.text)
            #         #     else:
            #         #         if m:
            #         #             json_obj['options'][opt.tag] = opt.text
            #         #         else:
            #         #             json_obj['options'][opt.tag] = "%s:%s" % (opt.attrib["type"], opt.text)
            #     else:
            #         if "format" in opt.attrib.keys():
            #             json_obj['options'][opt.tag] = "%s||%s" % (opt.attrib["format"], opt.text)
            #         else:
            #             json_obj['options'][opt.tag] = opt.text
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


class PipelineStop(object):

    @check_sig
    def POST(self):
        data = web.input()
        # client = data.client if hasattr(
        #     data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if not hasattr(data, "id") or data.id.strip() == "":
            info = {"success": False, "info": "停止任务需要参数: id!"}
            return json.dumps(info)
        json_data = {"id": data.id.strip(), "msg": "stop"}
        response = worker_client().run_cmd(json_data)
        return json.dumps(response)


class PipelinePause(object):

    @check_sig
    def POST(self):
        data = web.input()
        print data
        if not hasattr(data, "id") or data.id.strip() == "":
            info = {"success": False, "info": "暂停任务需要参数: id!"}
            return json.dumps(info)
        json_data = {"id": data.id.strip(), "msg": "pause"}
        response = worker_client().run_cmd(json_data)
        return json.dumps(response)


class PipelineStopPause(object):

    def POST(self):
        return self.GET()

    @check_sig
    def GET(self):
        data = web.input()
        print data
        if not hasattr(data, "id") or data.id.strip() == "":
            info = {"success": False, "info": "重新开始暂停任务需要参数 id!"}
            return json.dumps(info)
        json_data = {"id": data.id.strip(), "msg": "continue"}
        response = worker_client().run_cmd(json_data)
        return json.dumps(response)
