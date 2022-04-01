# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from biocluster.core.function import load_class_by_path
from biocluster.api.file.remote import RemoteFileManager
from mainapp.libs.signature import check_sig
import web
import json
from biocluster.core.exceptions import FileError
import os
import traceback
from biocluster.file import download

class FileCheck(object):

    @check_sig
    def POST(self):
        data = web.input()
        print "******"
        print data
        print "******"
        for i in ["format", "module", "type", "path"]:
            if not hasattr(data, i):
                msg = {"success": False, "info": "参数不全"}
                return json.dumps(msg)
            if getattr(data, i).strip() == "":
                msg = {"success": False, "info": "参数%s不能为空" % i}
                return json.dumps(msg)
        return self.check(data)

    def check(self, data):
        try:
            if data.format == "sequence.fastq_dir":  # sequence.fastq_dir
                info = {"success": True, "info": "检测通过"}
                return json.dumps(info)
            elif data.format == "sequence.fasta_dir":  # sequence.fasta_dir
                info = {"success": True, "info": "fasta_dir检测通过%s"}
                return json.dumps(info)
            file_obj = load_class_by_path(data.format, "File")()
            full_path = download(data.path)
            # file_manager = RemoteFileManager(data.path)
            # if file_manager.type != "local" and file_manager.type != "http":
            #     config = file_manager.config.get_netdata_config(file_manager.type)
            #     full_path = os.path.join(config[file_manager.type + "_path"], file_manager.path)
            # else:
            #     full_path = data.path
            file_obj.set_path(full_path)
            paths = data.module.split(".")
            function_name = "_".join(paths)
            if data.type.lower() == "tool":
                function_name += "_tool_check"
            elif data.type.lower() == "module":
                function_name += "_module_check"
            elif data.type.lower() == "workflow":
                function_name += "_workflow_check"
            else:
                msg = {"success": False, "info": "参数type不能为" % data.type}
                return json.dumps(msg)
            if hasattr(data, "check") and data.check:
                if hasattr(file_obj, data.check):
                    getattr(file_obj, data.check)()
                else:
                    info = {"success": False, "info": "文件类%s中未定义指定的检测函数%s!" % (data.format, data.check)}
                    return json.dumps(info)
            else:
                if hasattr(file_obj, function_name):
                    getattr(file_obj, function_name)()
                else:
                    getattr(file_obj, "check")()
        except ImportError, e:
            info = {"success": False, "info": "文件模块错误: %s" % e}
            return json.dumps(info)
        except FileError, e:
            info = {"success": False, "info": "文件检测错误: %s" % e}
            return json.dumps(info)
        except Exception, e:
            # exstr = traceback.format_exc()
            # print exstr
            info = {"success": False, "info": "错误: %s" % e}
            return json.dumps(info)
        else:
            info = {"success": True, "info": "检测通过"}
            return json.dumps(info)

    def check_group(self, data, file_list):
        file_obj = load_class_by_path(data.format, "File")()
        file_manager = RemoteFileManager(data.path)
        if file_manager.type != "local" and file_manager.type != "http":
            config = file_manager.config.get_netdata_config(file_manager.type)
            full_path = os.path.join(config[file_manager.type + "_path"], file_manager.path)
        else:
            full_path = data.path
        file_obj.set_path(full_path)
        try:
            file_obj.get_info()
        except:
            info = {"success": False, "info": "group文件检测错误"}
            return json.dumps(info)
        if not isinstance(file_list, (dict)):
            info = {"success": True, "info": "检测通过"}
            return json.dumps(info)
        else:
            try:
                # file_list = eval(file_list)
                sample_list = file_obj.prop["sample"]
                new_list = [file_list[x][0] for x in file_list.keys()]
                new = []
                for item in new_list:
                    item = item.decode("unicode_escape")
                    new.append(item)
                print "sample list: " + str(sample_list)
                print "new_list: " + str(new)
                for new_name in file_list.keys():
                    item = file_list[new_name][0].decode("unicode_escape")
                    if item not in sample_list:
                        raise FileError("分组文件中样本名与检测的样本信息不匹配")
                for sample in sample_list:
                    if sample not in new:
                        raise FileError("分组文件中样本名与检测的样本信息不匹配")
            except FileError as e:
                info = {"success": False, "info": "错误:%s" % e}
                print "group 文件的info：" + str(info)
                return json.dumps(info)
            else:
                info = {"success": True, "info": "检测通过"}
                return json.dumps(info)


class TestData(object):
    def __init__(self):
        self.format = None
        self.module = None
        self.type = None
        self.path = None
        self.check = None


class MultiFileCheck(object):

    def __init__(self):
        self.checker = FileCheck()

    @check_sig
    def POST(self):
        data = web.input()
        print "************"
        print data
        print "************"
        if not hasattr(data, "content"):
                msg = {"success": False, "info": "缺少参数content！"}
                return json.dumps(msg)
        try:
            json_obj = json.loads(data.content)
        except ValueError:
            info = {"success": False, "info": "content必须为Json格式!"}
            return json.dumps(info)
        else:
            for i in ["module", "type", "files", "file_list"]:
                if i not in json_obj.keys():
                    info = {"success": False, "info": "缺少参数:%s" % i}
                    return json.dumps(info)
            if not isinstance(json_obj["files"], list) or len(json_obj["files"]) == 0:
                info = {"success": False, "info": "必须至少有一个检测文件!"}
                return json.dumps(info)
            message = {"success": True, "files": []}
            for f in json_obj["files"]:
                d = TestData()
                d.module = json_obj["module"]
                d.type = json_obj["type"]
                d.format = f["format"]
                d.path = f["path"]
                if d.format == "meta.otu.group_table" and isinstance(json_obj["file_list"], (dict)):
                    result = json.loads(self.checker.check_group(d, json_obj["file_list"]))
                else:
                    result = json.loads(self.checker.check(d))
                if result["success"]:
                    x = {
                        "name": f["name"],
                        "pass": True
                    }
                else:
                    x = {
                        "name": f["name"],
                        "pass": False,
                        "info": result["info"]
                    }
                message["files"].append(x)
            return json.dumps(message)
