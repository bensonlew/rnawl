# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from .remote import RemoteFile
from .lib.transfer import MultiFileTransfer
from biocluster.config import Config
from boto.s3.bucket import Bucket
import re
import os
import json
# import xattr
# import shutil
import logging
from biocluster.core.exceptions import RunningError
from biocluster.core.function import get_hostname
import sys
from biocluster.proto import filetrans_pb2

PY3 = sys.version_info[0] == 3


class S3cache(RemoteFile):
    def __init__(self, type_name, path):
        super(S3cache, self).__init__(type_name, path)
        if not PY3:
            if isinstance(path, unicode):
                path = path.encode("utf-8")
        self.config = Config()
        self.type_name = type_name
        self._region = None
        self._bucket = None
        self._key = None
        # self._cache_base = None
        # self._cache_path = None
        self._file_list = None
        self._alias = None
        m = re.match(r'^(.*){(.*)}$', path)
        if m:
            self._path = m.group(1)
            self._alias = m.group(2)
        else:
            tmp = re.split(";;", path)
            if len(tmp) > 1:
                self._file_list = json.loads(tmp[1])
            self._path = tmp[0]
        self._type = type_name
        self._conn = None
        self._key_obj = None
        self._bucket_obj = None
        self.__parse_path()
        self.transfer = MultiFileTransfer()
        # self.transfer.base_path = self.cache_base
        self.logger = logging.getLogger("S3Cache")
        logger = logging.getLogger("boto")
        logger.setLevel(logging.INFO)

    def __parse_path(self):
        m = re.match(r"^([\w+\-]+)://([\w\-]+)/(.*)$", self._path)
        if not m:
            raise Exception("路径%s格式错误!" % self._path)
        self._region = m.group(1)
        self._bucket = m.group(2)
        self._key = m.group(3)

    # @property
    # def cache_base(self):
    #     if not self._cache_base:
    #         self._cache_base = self.config.get_netdata_config(self.type_name)["%s_cache_dir" % self.type_name]
    #     return self._cache_base

    # @property
    # def cache_path(self, ):
    #     if not self._cache_path:
    #         self._cache_path = os.path.join(self.cache_base, self._region, self._bucket, self._key)
    #     return self._cache_path

    @property
    def conn(self):
        if not self._conn:
            self._conn = self.config.get_rgw_conn(self._region, self._bucket)
        return self._conn

    @property
    def bucket(self):
        if not self._bucket_obj:
            self._bucket_obj = Bucket(connection=self.conn, name=self._bucket)
        return self._bucket_obj

    @property
    def key(self):
        if not self._key_obj:
            self._key_obj = self.bucket.get_key(self._key)
        return self._key_obj

    def download(self, to_path):
        if not PY3:
            if isinstance(to_path, unicode):
                to_path = to_path.encode("utf-8")
        if not self._key:
            raise RunningError("不能下载整个Bucket!", None, '007')
        if os.path.isfile(to_path):
            os.remove(to_path)
        if not os.path.exists(to_path):
            os.makedirs(to_path)
        if self._file_list:
            file_list = []
            for my_dict in self._file_list:
                m1 = re.match(r"^([\w+\-]+)://([\w\-]+)/(.*)$", my_dict["file_path"])
                if m1:
                    region = m1.group(1)
                    bucket = m1.group(2)
                    key = m1.group(3)
                    source = my_dict["file_path"]
                    conn = self.config.get_rgw_conn(region, bucket)
                    bucket_obj = Bucket(connection=conn, name=bucket)
                    key = bucket_obj.get_key(key)
                else:
                    if re.match(self._bucket, my_dict["file_path"]):
                        my_dict["file_path"] = os.path.relpath(my_dict["file_path"], self._bucket)
                    source = os.path.join("%s://%s" % (self._region, self._bucket), my_dict["file_path"])
                    key = self.bucket.get_key(my_dict["file_path"])
                if not key:
                    # raise Exception("文件{}不存在".format(source))
                    raise RunningError("文件%s不存在", (source,), '008')
                # cache_path = os.path.join(self.cache_base, self._region, self._bucket, my_dict["file_path"])
                target_path = os.path.join(to_path, my_dict["alias"])
                # tag = None
                # try:
                #     tag = xattr.getxattr(cache_path, "user.etag")
                # except IOError:
                #     pass
                # if not (os.path.isfile(cache_path) and key and tag and tag == key.etag.strip('"')):
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
            if len(file_list) > 0:
                self.transfer.trans_file(self.file_iter(file_list))
                self.transfer.perform()

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

            return to_path
        else:
            if re.search(r"/$", self._key):
                # key_list = []
                # for key in self.bucket.list(prefix=self._key):
                #     if re.search(r"/$", key.name):
                #         continue
                #     key_list.append(key)
                # if len(key_list) == 0:
                #     # raise Exception("文件夹{}为空!".format(self.path))
                #     raise RunningError("文件夹%s为空!", (self.path,), '009')
                from_uri = "%s://%s/%s" % (self._region, self._bucket, self._key)
                if self._alias:
                    target_path = os.path.join(to_path, self._alias)
                else:
                    target_path = os.path.join(to_path, os.path.basename(self._key.rstrip("/")))
                if not PY3:
                    if isinstance(from_uri, unicode):
                        from_uri = from_uri.encode("utf-8")
                    if isinstance(target_path, unicode):
                        target_path = target_path.encode("utf-8")
                # file_list = []
                # for key in key_list:
                #     # cache_path = os.path.join(self.cache_base, self._region, self._bucket, key.name)
                #     # tag = None
                #     # try:
                #     #     tag = xattr.getxattr(cache_path, "user.etag")
                #     # except IOError:
                #     #     pass
                #     from_uri = "%s://%s/%s" % (self._region, self._bucket, key.name)
                #     if self._alias:
                #         target_path = os.path.join(to_path, self._alias,
                #                                    os.path.relpath(key.name, self._key))
                #     else:
                #         target_path = os.path.join(to_path,
                #                                    os.path.relpath(key.name, os.path.dirname(self._key.rstrip("/"))))
                #     if not PY3:
                #         if isinstance(from_uri, unicode):
                #             from_uri = from_uri.encode("utf-8")
                #         if isinstance(target_path, unicode):
                #             target_path = target_path.encode("utf-8")
                #     # if not (os.path.isfile(cache_path) and key and tag and key.etag and tag == key.etag.strip('"')):
                #         # self.transfer.add(from_uri=from_uri, to_uri=target_path)
                #     file_list.append([target_path, from_uri])
                #     # else:
                #     #     self.logger.info("文件%s找到本地缓存文件%s，跳过下载..." % (from_uri, cache_path))
                #     #     if os.path.exists(target_path):
                #     #         if os.path.islink(target_path):
                #     #             os.remove(target_path)
                #     #         elif os.path.isdir(target_path):
                #     #             shutil.rmtree(target_path)
                #     #         else:
                #     #             os.remove(target_path)
                #     #     dir_name = os.path.dirname(target_path)
                #     #     if not os.path.exists(dir_name):
                #     #         os.makedirs(dir_name)
                #     #     os.link(cache_path, target_path)
                # # self.transfer.wait()
                # # for cache, target in file_list:
                # #     if os.path.exists(cache):
                # #         if os.path.islink(cache):
                # #             os.remove(target)
                # #         elif os.path.isdir(cache):
                # #             shutil.rmtree(cache)
                # #         else:
                # #             os.remove(cache)
                # #     dir_name = os.path.dirname(cache)
                # #     if not os.path.exists(dir_name):
                # #         os.makedirs(dir_name)
                # #     os.link(target, cache)
                # if len(key_list) > 0:
                self.transfer.add_download(from_uri, target_path, from_uri)
                self.transfer.perform()
                return target_path
            else:
                # if not self.key:
                #     # raise Exception("文件{}不存在".format(self.path))
                #     raise RunningError("文件%s不存在", (self.path,), '008')
                # tag = None
                # try:
                #     tag = xattr.getxattr(self.cache_path, "user.etag")
                # except IOError:
                #     pass
                if self._alias:
                    target_path = os.path.join(to_path, self._alias)
                else:
                    target_path = os.path.join(to_path, os.path.basename(self._key.rstrip("/")))
                if not PY3:
                    if isinstance(target_path, unicode):
                        target_path = target_path.encode("utf-8")
                    if isinstance(self.path, unicode):
                        self._path = self.path.encode("utf-8")
                # if not (os.path.isfile(self.cache_path) and self.key and tag and tag == self.key.etag.strip('"')):
                file_list = list()
                file_list.append([target_path, self.path])
                self.transfer.trans_file(self.file_iter(file_list))
                self.transfer.perform()
                # self.transfer.add(from_uri=self.path, to_uri=target_path)
                # self.transfer.wait()
                # if os.path.exists(self.cache_path):
                #     if os.path.islink(self.cache_path):
                #         os.remove(self.cache_path)
                #     elif os.path.isdir(self.cache_path):
                #         shutil.rmtree(self.cache_path)
                #     else:
                #         os.remove(self.cache_path)
                # dir_name = os.path.dirname(self.cache_path)
                # if not os.path.exists(dir_name):
                #     os.makedirs(dir_name)
                # os.link(target_path, self.cache_path)
                # else:
                #     self.logger.info("文件%s找到本地缓存文件%s，跳过下载..." % (self.path, self.cache_path))
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
                #     os.link(self.cache_path, target_path)
                return target_path

    def upload(self, from_path):
        if not PY3:
            if isinstance(from_path, unicode):
                from_path = from_path.encode("utf-8")
        if not os.path.exists(from_path):
            raise Exception("源文件%s不存存在" % from_path)
        if os.path.isdir(from_path):
            file_list = []
            for root, dirs, files in os.walk(from_path):
                rel_path = os.path.relpath(root, from_path)
                for i_file in files:
                    if rel_path == ".":
                        key_path = os.path.join(self._key, i_file)
                    else:
                        key_path = os.path.join(self._key, rel_path, i_file)
                    i_file_path = os.path.join(root, i_file)
                    if os.path.islink(i_file_path):
                        real_path = os.path.realpath(i_file_path)
                        if os.path.exists(real_path):
                            target = "%s://%s/%s" % (self._region, self._bucket, key_path)
                            # cache_path = os.path.join(self.cache_base, self._region, self._bucket, key_path)
                            # self.transfer.add(from_uri=real_path, to_uri=target)
                            file_list.append([target, real_path])
                    else:
                        target = "%s://%s/%s" % (self._region, self._bucket, key_path)
                        # cache_path = os.path.join(self.cache_base, self._region, self._bucket, key_path)
                        # self.transfer.add(from_uri=i_file_path, to_uri=target)
                        file_list.append([target, i_file_path])
            # self.transfer.wait()
            # for cache, source, key_path in file_list:
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
            #     if os.path.exists(cache):
            #         os.remove(cache)
            #     os.link(source, cache)
                # key = self.bucket.get_key(key_path)
                # if key:
                #     if key.etag:
                #         xattr.setxattr(cache, "user.etag", key.etag.strip('"'))
                #     key.close(fast=True)
            if len(file_list) > 0:
                self.transfer.trans_file(self.file_iter(file_list, True))
                self.transfer.perform()

        else:
            self._path = os.path.join(self._path, os.path.basename(from_path))
            self.__parse_path()
            real_path = from_path
            if os.path.islink(from_path):
                real_path = os.path.realpath(from_path)
                if not os.path.exists(real_path):
                    raise Exception("源文件%s是一个无效的软链接!" % from_path)
            # self.transfer.add(from_uri=real_path, to_uri=self.path)
            # self.transfer.wait()
            file_list = list()
            file_list.append([self.path, real_path])
            self.transfer.trans_file(self.file_iter(file_list, True))
            self.transfer.perform()
            # if os.path.exists(self.cache_path):
            #     if os.path.isdir(self.cache_path):
            #         shutil.rmtree(self.cache_path)
            #     else:
            #         os.remove(self.cache_path)
            # dir_name = os.path.basename(self.cache_path)
            # if not os.path.exists(dir_name):
            #     os.makedirs(dir_name)
            # if os.path.exists(self.cache_path):
            #     os.remove(self.cache_path)
            # os.link(real_path, self.cache_path)
            # if self.key:
            #     if self.key.etag:
            #         xattr.setxattr(self.cache_path, "user.etag", self.key.etag.strip('"'))
            #     self.key.close(fast=True)
        return self.path

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
