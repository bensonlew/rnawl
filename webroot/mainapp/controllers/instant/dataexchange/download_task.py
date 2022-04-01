# -*- coding: utf-8 -*-
# __author__ = 'xuting'
from __future__ import division
import web
import json
import os
import re
import datetime
from biocluster.config import Config
from mainapp.models.data_exchange import Download, Identity


class DownloadTask(object):
    def POST(self):
        data = web.input()
        ip = data.ip
        code = data.identity
        user = data.user
        mode = data.mode
        info = Identity().get_task_id(code)
        if not info["success"]:
            info["data"] = ""
            return json.dumps(info)
        else:
            # 获取主流程上的文件
            path = Download().get_path_by_workflow_id(info["task_id"])
            if path == "empty":
                info["info"] = "找不到任务 {}, 距该任务的运行时间可能太久，该任务的信息已经从数据库里被移除".format(info["task_id"])
                info["success"] = False
                return json.dumps(info)
            full_path = os.path.join(Config().WORK_DIR, path)
            file_list = list()
            for d in os.walk(full_path):
                # d 为一个元祖, 里面有4个元素， 第一个元素是字符串，代表上级目录，
                # 第二个元素是一个列表，里面的内容是这个上级目录下的文件夹，
                # 第三个元素是一个列表，里面的内容是这个上级文件目录下的文件
                # 第四个元素是一个数字，代表这个文件的类型，1代表是工作流(项目任务)的文件， 2代表是report上的文件
                for f in d[2]:
                    full_path = os.path.join(d[0], f)
                    rel_path = re.sub(Config().WORK_DIR, "", full_path)
                    rel_path = re.sub("^/*\d+/", "", rel_path)
                    rel_path = rel_path.lstrip("/")
                    file_size = self.get_size(full_path)
                    file_list.append([full_path, file_size, rel_path, 1])
            if len(file_list) == 0:
                info["success"] = False
                info["info"] = "任务{}的结果文件夹为空文件， 该任务可能没有正常运行".format(info["task_id"])

            # 获取report上的文件
            report_path = Download().get_report_by_workflow_id(info["task_id"])
            if report_path == "empty":
                info["info"] = "未能找到任务{}的运行记录，请联系相关人员解决该问题".format(info["task_id"])
                info["success"] = False
                return json.dumps(info)
            self.config = Config().get_netdata_config(mode)
            full_path = os.path.join(self.config[mode + "_path"], report_path)
            if os.path.isdir(full_path):
                for d in os.walk(full_path):
                    for f in d[2]:
                        full_path = os.path.join(d[0], f)
                        pattern = os.path.join(self.config[mode + "_path"], "rerewrweset/files/")
                        rel_path = re.sub(pattern, "", full_path)
                        rel_path = re.sub("^/*\w+/\d+/\w+/", "", rel_path)
                        rel_path = rel_path.lstrip("/")
                        file_size = self.get_size(full_path)
                        file_list.append([full_path, file_size, rel_path, 2])

            info["data"] = file_list
            download_info = dict()
            download_info["code"] = code
            download_info["request_time"] = datetime.datetime.now().strftime('%Y-%m-%d %X')
            download_info["ip"] = ip
            download_info["user"] = user
            Identity().add_download_record(download_info)
            return json.dumps(info)

    def get_size(self, path):
        """
        获取文件大小
        """
        b = os.path.getsize(path)
        if b / 1024 < 1:
            return "{}B".format(b)
        elif b >= 1024 and b < 1024 * 1024:
            return "{:.3f}KB".format(b / 1024)
        elif b >= 1024 * 1024 and b < 1024 * 1024 * 1024:
            return "{:.3f}MB".format(b / (1024 * 1024))
        else:
            return "{:.3f}GB".format(b / (1024 * 1024 * 1024))
