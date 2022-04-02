# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from biocluster.file import exists
import json
import re
from .lib.transfer import MultiFileTransfer
import logging
from biocluster.config import Config
import os
from biocluster.core.exceptions import RunningError
# import xattr
# import shutil
import sys
from biocluster.core.function import get_hostname
from biocluster.proto import filetrans_pb2

PY3 = sys.version_info[0] == 3


class Httpcache(object):
    def __init__(self, type_name, path):
        self.config = Config()
        self.type_name = type_name
        # self._cache_base = None
        # self._cache_path = None
        self._file_list = None
        self._alias = None
        m = re.match(r"^(.*){(.*)}$", path)
        if m:
            self._path = self.config.convert_path_to_http(
                    os.path.join(self.config.get_netdata_config(type_name)[type_name + "_path"], m.group(1)))
            self._alias = m.group(2)
        else:
            tmp = re.split(";;", path)
            if len(tmp) > 1:
                self._file_list = json.loads(tmp[1])
            if self.type_name == "http":
                self._path = self.config.convert_path_to_http(tmp[0])
            else:
                self._path = self.config.convert_path_to_http(
                    os.path.join(self.config.get_netdata_config(type_name)[type_name + "_path"], tmp[0]))
        self.transfer = MultiFileTransfer()
        # self.transfer.base_path = self.cache_base
        self.logger = logging.getLogger("HttpCache")

    # @property
    # def cache_base(self):
    #     if not self._cache_base:
    #         self._cache_base = self.config.get_netdata_config("http")["http_cache_dir"]
    #     return self._cache_base

    def download(self, to_path):
        if os.path.isfile(to_path):
            os.remove(to_path)
        if not os.path.exists(to_path):
            os.makedirs(to_path)
        if self._file_list:
            file_list = []
            for my_dict in self._file_list:
                if not re.match("rerewrweset", my_dict["file_path"]):
                    my_dict["file_path"] = "rerewrweset/" + my_dict["file_path"]
                if self.type_name == "http":
                    source = os.path.join(self._path, my_dict["file_path"])
                else:
                    source = self.config.convert_path_to_http(os.path.join(
                        self.config.get_netdata_config(self.type_name)[self.type_name + "_path"], my_dict["file_path"]))
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
                target_path = os.path.join(to_path, my_dict["alias"])
                # if not (os.path.isfile(cache_path) and tag and tag == info["etag"].strip('"')):
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
        else:
            if re.search(r"/$", self._path):
                # url_list = HttpDir(self.config.convert_path_to_http(self._path)).get_files()
                # if len(url_list) == 0:
                #     # raise Exception("文件夹{}为空!".format(self.path))
                #     raise RunningError("文件夹%s为空!", (self._path,), '009')
                if self._alias:
                    target_path = os.path.join(to_path, self._alias)
                else:
                    target_path = os.path.join(to_path, os.path.basename(self._path.rstrip("/")))
                if not PY3:
                    if isinstance(self._path, unicode):
                        self._path = self._path.encode("utf-8")
                    if isinstance(target_path, unicode):
                        target_path = target_path.encode("utf-8")
                # file_list = []
                # for f in url_list:
                #     # cache_path = _get_http_cache_path(f)
                #     # tag = None
                #     # try:
                #     #     tag = xattr.getxattr(cache_path, "user.etag")
                #     # except IOError:
                #     #     pass
                #     # info = _get_http_info(f)
                #     target_path = os.path.join(to_path, os.path.relpath(f, os.path.dirname(self._path.rstrip("/"))))
                #     if not PY3:
                #         if isinstance(f, unicode):
                #             f = f.encode("utf-8")
                #         if isinstance(target_path, unicode):
                #             target_path = target_path.encode("utf-8")
                    # if os.path.isfile(cache_path) and tag:
                    #     self.logger.info("文件%s找到本地缓存文件%s，跳过下载..." % (f, cache_path))
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
                    # else:
                        # self.transfer.add(from_uri=f, to_uri=target_path)
                #    file_list.append([target_path, f])
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
                # if len(url_list) > 0:
                self.transfer.add_download(self._path, target_path, self._path)
                self.transfer.perform()
                return os.path.join(to_path, os.path.basename(self._path.rstrip("/")))
            else:
                # if not exists(self._path):
                #     # raise Exception("文件{}不存在".format(self.path))
                #     raise RunningError("文件%s不存在", (self._path,), '008')
                # cache_path = _get_http_cache_path(self._path)
                # tag = None
                # try:
                #     tag = xattr.getxattr(cache_path, "user.etag")
                # except IOError:
                #     pass
                # info = _get_http_info(self._path)
                if self._alias:
                    target_path = os.path.join(to_path, self._alias)
                else:
                    target_path = os.path.join(to_path, os.path.basename(self._path.rstrip("/")))
                if not PY3:
                    if isinstance(self._path, unicode):
                        self._path = self._path.encode("utf-8")
                    if isinstance(target_path, unicode):
                        target_path = target_path.encode("utf-8")
                    # if not (os.path.isfile(cache_path) and tag and tag == info["etag"].strip('"')):
                    #
                    # self.transfer.add(from_uri=self._path, to_uri=target_path)
                    # self.transfer.wait()
                    # if os.path.exists(cache_path):
                    #     if os.path.islink(cache_path):
                    #         os.remove(cache_path)
                    #     elif os.path.isdir(cache_path):
                    #         shutil.rmtree(cache_path)
                    #     else:
                    #         os.remove(cache_path)
                    # dir_name = os.path.dirname(cache_path)
                    # if not os.path.exists(dir_name):
                    #     os.makedirs(dir_name)
                    # os.link(target_path, cache_path)
                # file_list = list()
                # file_list.append([target_path, self._path])
                # self.transfer.trans_file(self.file_iter(file_list))
                # self.transfer.perform()
                # else:
                #     self.logger.info("文件%s找到本地缓存文件%s，跳过下载..." % (self._path, cache_path))
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
                file_list = list()
                file_list.append([target_path, self._path])
                self.transfer.trans_file(self.file_iter(file_list))
                self.transfer.perform()
                return target_path

    def upload(self, from_path):
        raise Exception("http类型不支持上传!")

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
