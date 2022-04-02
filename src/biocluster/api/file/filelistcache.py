# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from biocluster.file import exists, download
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
import sys
from biocluster.core.function import get_hostname
from biocluster.proto import filetrans_pb2

PY3 = sys.version_info[0] == 3


class Filelistcache(object):
    def __init__(self, type_name, path):
        self.config = Config()
        self.type_name = type_name
        self._path = path
        self._option_name = None
        self._olddata_type = None
        self._alias = None
        self._file_list = []
        self.transfer = MultiFileTransfer()
        self.logger = logging.getLogger("FilelistCache")
        self._mapping_file = None
        logger = logging.getLogger("boto")
        logger.setLevel(logging.INFO)
        m = re.match(r"filelist\[(.*)\]\((.*)\):(.*){(.*)}$", path)
        if m:
            self._option_name = m.group(1)
            self._olddata_type = m.group(2)
            mapping_file = m.group(3)
            self._alias = m.group(4)
            mapping_file_cache = download(mapping_file)
            self._mapping_file = mapping_file_cache
            with open(mapping_file_cache, "r") as f:
                data = json.load(f)
            if self._option_name in data.keys():
                if isinstance(data[self._option_name], list) and len(data[self._option_name]) > 0:
                    self._file_list = data[self._option_name]
                else:
                    raise Exception("文件夹列表文件%s内容不正确或者参数%s的列表内容为空!" %
                                    (mapping_file, self._option_name))
            else:
                raise Exception("文件夹列表文件%s内容中没有包含参数%s信息!" % (mapping_file, self._option_name))
        else:
            raise Exception("filelist格式不正确！")

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
                target_path = os.path.join(to_path, self._alias, f["alias"])
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
                # conn.close()

            else:
                path = "%s:%s" % (self._olddata_type, f["file_path"])
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
                target_path = os.path.join(to_path, self._alias, f["alias"])
                if not PY3:
                    if isinstance(source, unicode):
                        source = source.encode("utf-8")
                    if isinstance(target_path, unicode):
                        target_path = target_path.encode("utf-8")
                # if not (os.path.isfile(cache_path) and tag and "etag" in info.keys() and info["etag"]
                #         and info["etag"].strip('"')):
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
        if self._mapping_file and os.path.exists(self._mapping_file):
            if not os.path.exists(to_path):
                os.makedirs(to_path)
            p = os.path.join(to_path, "mapping_file.txt")
            if os.path.exists(p):
                os.remove(p)
            os.link(self._mapping_file, p)
        return os.path.join(to_path, self._alias)

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
