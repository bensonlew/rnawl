# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
# lasted modified by hd
import json
import os
import urllib
import random
import hashlib
import time
import sys
import traceback
from biocluster.wpm.log import Log
from biocluster.core.function import CJsonEncoder, filter_error_info
import copy


class Sanger(Log):

    def __init__(self, data):
        super(Sanger, self).__init__(data)
        self._client = "client01"
        self._key = "e0bb04dd111155b2e6bc6db26d0e1fef"
        # self._url = "http://api.sanger.com/task/add_task_log"
        self._binds_id = "5f4324f09b7900009300805c"
        self._interface_id = 57
        self._env_name = "online"
        self._url = "http://apicenter.lab.majorbio.com/index/in"
        # self._post_data = "%s&%s" % (self.get_sig(), self.post_data)

    def update(self):
        self.send()

    @property
    def post_data(self):
        my_content = self.data["sync_task_log"]
        my_data = dict()
        base_info = dict()
        if "task" in my_content.keys():
            base_info["basis"] = copy.deepcopy(my_content["task"])
        if "log" in my_content.keys():
            base_info["log"] = my_content["log"]
        if "base_path" in my_content.keys():
            base_info['basis']["base_path"] = my_content["base_path"]
            # noinspection PyBroadException
            try:
                project_sn = my_content["base_path"].split('/')[3]
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                print e
                sys.stdout.flush()
            else:
                base_info['basis']["project_sn"] = project_sn
        if "region" in my_content.keys():
            base_info['basis']["region"] = my_content["region"]
        if "files" in my_content.keys() and len(my_content["files"]) > 5000:
            files_len = len(my_content["files"])
            dirs = {}
            if "dirs" in my_content.keys():
                dirs = my_content["dirs"]
            start_index = 0
            if files_len % 1000 > 0:
                count = files_len / 1000 + 1
            else:
                count = files_len / 1000
            for i in xrange(count):
                end_index = start_index + 1000
                if end_index > files_len:
                    split_data = {
                        "files": my_content["files"][start_index:]
                    }
                else:
                    split_data = {
                        "files": my_content["files"][start_index:end_index]
                    }
                if dirs:
                    split_dirs = []
                    for f in split_data["files"]:
                        for d in dirs:
                            rel_path = os.path.relpath(f["path"], d["path"])
                            if not rel_path.startswith("../"):
                                split_dirs.append(d)
                    split_data["dirs"] = split_dirs
                start_index = end_index
                split_data.update(base_info)
                my_data["sync_task_log"] = json.dumps(split_data, cls=CJsonEncoder)
                yield urllib.urlencode(my_data)
        else:
            if "files" in my_content.keys():
                base_info["files"] = my_content["files"]
            if "dirs" in my_content.keys():
                base_info["dirs"] = my_content["dirs"]
            my_data["sync_task_log"] = json.dumps(base_info, cls=CJsonEncoder)
            yield urllib.urlencode(my_data)

    def get_sig(self):
        """
        增加bind_id等前端需要的验证
        :return:
        """
        nonce = str(random.randint(1000, 10000))
        timestamp = str(int(time.time()))
        x_list = [self._key, timestamp, nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        sig = sha1.hexdigest()
        signature = {
            "client": self._client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": sig,
            "binds_id": self._binds_id,
            "interface_id": self._interface_id,
            "env_name": self._env_name
        }
        return urllib.urlencode(signature)
