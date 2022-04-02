# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import json
from .lib.transfer import MultiFileTransfer
import logging
from biocluster.config import Config
import os
from biocluster.core.exceptions import RunningError
# import xattr
# import shutil
from boto.s3.bucket import Bucket
import re
import random
import time
import hashlib
import sys
from biocluster.file import exists
from biocluster.core.function import get_hostname
from biocluster.proto import filetrans_pb2
try:
    import urllib2
    from urllib import urlencode
except ImportError:
    from urllib import request as urllib2
    from urllib.parse import urlencode

PY3 = sys.version_info[0] == 3

API_URL = {
    # "sanger": "http://api.sanger.com/file/check_mapping_file",
    # "tsanger": "http://api.tsg.com/file/check_mapping_file",
    "sanger": "http://openapi.labsanger.sanger.com/file/check_mapping_file",
    "tsanger": "http://openapi.nsg.com/file/check_mapping_file",
}

CLIENT_KEY = {
    "client01": "1ZYw71APsQ",
    "client03": "hM4uZcGs9d"
}


class Fileapicache(object):
    def __init__(self, type_name, path):
        self.config = Config()
        self.type_name = type_name
        self._path = path
        self._task_id = None
        self._api_client = None
        self._data_type = None
        self._alias = None
        self._file_list = []
        self.transfer = MultiFileTransfer()
        self.logger = logging.getLogger("FileAPICache")
        logger = logging.getLogger("boto")
        logger.setLevel(logging.INFO)
        m = re.match(r"fileapi\[(.*)\]\((.*)\):(.*)$", path)
        if m:
            self._api_client = m.group(1)
            self._data_type = self.config.get_api_type(self._api_client)
            self._task_id = m.group(2)
            self._path = m.group(3)
            list_data = self.get_file_list()
            if list_data:
                self.logger.debug("获取到接口数据:%s" % list_data)
                try:
                    data = json.loads(list_data)
                except:
                    raise Exception("fileapi返回数据不正确！")
                else:
                    if not isinstance(data, dict):
                        raise Exception("fileapi返回数据不正确！")
                    elif data["success"] is True or data["success"] == "true":
                        self._file_list = data["d"]
                        if not isinstance(self._file_list, list):
                            raise Exception("fileapi返回数据格式不正确！")
                    else:
                        raise Exception("fileapi未能正确获取数据！")

            else:
                raise Exception("fileapi返回数据不正确！")
        else:
            raise Exception("fileapi格式不正确！")

    def get_file_list(self):
        http_handler = urllib2.HTTPHandler(debuglevel=1)
        https_handler = urllib2.HTTPSHandler(debuglevel=1)
        opener = urllib2.build_opener(http_handler, https_handler)
        urllib2.install_opener(opener)
        data = {
                "task_id": self._task_id,
                "params_path": self._path,
                "ceshi": 1
        }
        post_data = "%s&%s" % (self.get_sig(), urlencode(data))
        request = urllib2.Request(API_URL[self._data_type], post_data)
        response = urllib2.urlopen(request, timeout=60)
        response_data = response.read()
        response.close()
        return response_data

    def get_sig(self):
        nonce = str(random.randint(1000, 10000))
        timestamp = str(int(time.time()))
        x_list = [CLIENT_KEY[self._api_client], timestamp, nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        sig = sha1.hexdigest()
        signature = {
            "client": self._api_client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": sig
        }
        return urlencode(signature)

    def download(self, to_path):
        file_list = []
        for f in self._file_list:
            if re.match(r"^([\w\-]+)://.*", f["file_path"]):
                m1 = re.match(r"^([\w+\-]+)://([\w\-]+)/(.*)$", f["file_path"])
                region = m1.group(1)
                bucket = m1.group(2)
                key = m1.group(3)
                source = f["file_path"]
                conn = self.config.get_rgw_conn(region, bucket)
                bucket_obj = Bucket(connection=conn, name=bucket)
                key = bucket_obj.get_key(key)
                if not key:
                    # raise Exception("文件{}不存在".format(source))
                    raise RunningError("文件%s不存在", (source,), '008')
                # cache_path = _get_s3_cache_path(f["file_path"])
                # tag = None
                # try:
                #     tag = xattr.getxattr(cache_path, "user.etag")
                # except IOError:
                #     pass
                target_path = os.path.join(to_path, f["alias"])
                if not PY3:
                    if isinstance(source, unicode):
                        source = source.encode("utf-8")
                    if isinstance(target_path, unicode):
                        target_path = target_path.encode("utf-8")
                # if not (os.path.isfile(cache_path) and key and tag and key.etag and tag == key.etag.strip('"')):
                    # self.transfer.add(from_uri=source, to_uri=target_path)
                file_list.append([target_path, source])
                # else:
                #     self.logger.info("文件%s找到本地缓存文件%s，跳过下载..." % (source, cache_path))
                #     if os.path.exists(target_path):
                #         if os.path.islink(target_path):
                #             os.remove(target_path)
                #         elif os.path.isdir(target_path):
                #             shutil.rmtree(target_path)
                #         else:
                #             os.remove(target_path)
                #     dir_name = os.path.dirname(target_path)
                #     if not os.path.exists(dir_name):
                #         os.makedirs(dir_name)
                #     os.link(cache_path, target_path)
            else:
                path = "%s:%s" % (self._data_type, f["file_path"])
                source = self.config.convert_path_to_http(path)
                # info = _get_http_info(source)
                if not exists(source):
                    # raise Exception("文件{}不存在".format(source))
                    raise RunningError("文件%s不存在", (source,), '008')
                # cache_path = _get_http_cache_path(source)
                # tag = None
                # try:
                #     tag = xattr.getxattr(cache_path, "user.etag")
                # except IOError:
                #     pass
                target_path = os.path.join(to_path, f["alias"])
                if not PY3:
                    if isinstance(source, unicode):
                        source = source.encode("utf-8")
                    if isinstance(target_path, unicode):
                        target_path = target_path.encode("utf-8")
                # if not (os.path.isfile(cache_path) and tag and "etag" in info.keys() and info["etag"]
                #         and tag == info["etag"].strip('"')):
                    # self.transfer.add(from_uri=source, to_uri=target_path)
                file_list.append([target_path, source])
                # else:
                #     self.logger.info("文件%s找到本地缓存文件%s，跳过下载..." % (source, cache_path))
                #     if os.path.exists(target_path):
                #         if os.path.islink(target_path):
                #             os.remove(target_path)
                #         elif os.path.isdir(target_path):
                #             shutil.rmtree(target_path)
                #         else:
                #             os.remove(target_path)
                #     dir_name = os.path.dirname(target_path)
                #     if not os.path.exists(dir_name):
                #         os.makedirs(dir_name)
                #     os.link(cache_path, target_path)

        # self.transfer.wait()
        # for cache, target in file_list:
        #     if os.path.exists(cache):
        #         if os.path.islink(cache):
        #             os.remove(cache)
        #         elif os.path.isdir(cache):
        #             shutil.rmtree(cache)
        #         else:
        #             os.remove(cache)
        #     dir_name = os.path.dirname(cache)
        #     if not os.path.exists(dir_name):
        #         os.makedirs(dir_name)
        #     os.link(target, cache)
        if len(file_list) > 0:
            self.transfer.trans_file(self.file_iter(file_list))
            self.transfer.perform()
        return to_path

    def upload(self, from_path):
        raise Exception("文件列表类型不支持上传!")

    def file_iter(self, file_list, upload=False):
        for f in file_list:
            data = filetrans_pb2.FileTrans(
                workflow_id=self.config.current_workflow_id,
                tool_id=self.config.current_tool_id,
                host=get_hostname(),
                processid=int(os.getpid()),
                upload=upload,
                frompath=f[1],
                topath=f[0],
                basepath="",
                usecache=True
            )
            yield data
