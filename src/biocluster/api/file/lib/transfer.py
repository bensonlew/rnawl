# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
# from biocluster.core.singleton import singleton
from biocluster.config import Config
import os
import socket
# from multiprocessing import Queue, Value, Lock
# from Queue import Empty
# import datetime
import sys
import gevent
from biocluster.core.exceptions import RunningError, TransferError
# from biocluster.core.watcher import Watcher
# import gipc
# import setproctitle
# import time
import traceback
# import threading
# from biocluster.core.function import friendly_size
# from .s3 import S3FileTransfer
# from .http import HttpFileTransfer
import re
# import glob
# from biocluster.file import _get_s3_dir_list, HttpDir, exists, _get_s3_cache_path, _get_http_cache_path,
# _get_s3_key, _get_http_info
import logging

import grpc
from biocluster.proto import filetrans_pb2_grpc, filetrans_pb2, public_pb2
from biocluster.core.function import get_hostname, friendly_size
from collections import Iterable


class MultiFileTransfer(object):

    """
    用于大文件和大批量文件的传输
    """
    def __init__(self):
        # self.manager = TransferManager()
        self._files = {}
        self.config = Config()
        self.logger = logging.getLogger("MultiFileTransfer")
        self._wait_times = 0
        self._add_file_times = 0

    # def set_uncover_mode(self):
    #     """
    #     设置为不允许覆盖模式，当目标文件存在时报错!
    #     :return:
    #     """
    #     self.manager.overwrite = False

    def add_upload(self, from_path, to_path, base_path=None):
        """
        批量添加上传文件
        :param from_path: 需要上传的文件或文件夹，可使用通配符, 如 "/root/123/456/789/*.txt"， 只支持本地文件
        :param to_path: 上传的目标文件夹，必须以"/"结尾,只支持对象存储路径
        :param base_path: 上传文件路径的相对基础路径。用于决定上传文件/文件夹的目录结构，如：
        add_upload("/root/123/456/789/test.txt", "s3://test/"): 不写相对基础路径时，
        上传的结果路径为"s3://test/test.txt", 如果有多个文件重名，则会相互覆盖
        add_upload("/root/123/456/789/test.txt", "s3://test/", "/root/123/"): 指定了相对基础路径时，
        上传的结果路径为"s3://test/456/789/test.txt"，保留了相对基础路径之后的目录结构
        :return:
        """
        if not re.match(r"^[\w\-]+://\S+/.+$", to_path):
            raise Exception("上传路径%s格式不正确!" % to_path)
        if not to_path.endswith("/"):
            raise Exception("文件夹路径%s必须以\"/\"号结尾!" % to_path)
        if base_path:
            if not from_path.startswith(base_path):
                raise Exception("base_path %s与from_path %s不能匹配!" % (base_path, from_path))
        if not os.path.exists(from_path):
            raise RunningError("文件%s不存在", (from_path,), code="008")
        if self.config.get_file_type(from_path) != "local":
            raise Exception("上传路径%s非本地路径，不支持上传!" % from_path)
        if self.config.get_file_type(to_path) != "s3":
            raise Exception("上传路径%s非s3路径，不支持上传!" % to_path)

        # for f in glob.glob(os.path.abspath(from_path)):
        #     if os.path.isdir(f):
        #         for root, dirs, files in os.walk(f):
        #             if base_path:
        #                 rel_path = os.path.relpath(root, base_path)
        #                 for i_file in files:
        #                     source = os.path.join(root, i_file)
        #                     target = os.path.join(to_path, rel_path, i_file)
        #                     self._files[source] = (target, _get_s3_cache_path(target))
        #             else:
        #                 for i_file in files:
        #                     source = os.path.join(root, i_file)
        #                     target = os.path.join(to_path, i_file)
        #                     self._files[source] = (target,  _get_s3_cache_path(target))
        #     else:
        #         if base_path:
        #             target = os.path.join(to_path, os.path.relpath(f, base_path))
        #         else:
        #             target = os.path.join(to_path, os.path.basename(f))
        #         self._files[f] = (target,  _get_s3_cache_path(target))
        env_dist = os.environ
        tool_id = ""
        if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
            tool_id = self.config.current_tool_id
        data = filetrans_pb2.FileTrans(
            workflow_id=self.config.current_workflow_id,
            tool_id=tool_id,
            host=get_hostname(),
            processid=int(os.getpid()),
            upload=True,
            frompath=from_path,
            topath=to_path,
            basepath=base_path if base_path else "",
            usecache=True
        )
        self._add_file_trans(data)

    def _add_file_trans(self, data):
        self._add_file_times += 1
        env_dist = os.environ
        port = self.config.wfm_port
        if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
            port = self.config.ntm_port
        try:
            # hostname = socket.gethostname()
            # addr = socket.gethostbyname(socket.gethostname())
            # self.logger.info("hostname:%s" % hostname)
            self.logger.info("port:%s" % port)
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = filetrans_pb2_grpc.FileTransServerStub(channel)
                responses = stub.Trans(data)
                if not responses.ok:
                    self.logger.error("添加上传下载失败,原因: %s" % responses.reason)
                    raise TransferError("添加上传下载失败,原因: %s" % responses.reason)
                self._add_file_times = 0
                self.logger.info("添加上传下载成功")
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._add_file_times > 3:
                self.logger.error("添加上传下载发生错误超过3次,退出运行: %s" % e)
                raise e
            else:
                self.logger.error("添加上传下载发生错误10秒后重试: %s" % e)
                gevent.sleep(20)
                self._add_file_trans(data)

    def add_download(self, from_path, to_path, base_path=None):
        """
        批量添加下载文件
        :param from_path: 支持对象存储和http路径，会自动进行路径转换, 不支持通配符，文件夹必须以"/"结尾，否则无法正常工作
        :param to_path: 下载文件的目标路径文件夹
        :param base_path: 下载文件的相对基础路径，必须与from_path有相同的开头,用于决定下载文件/文件夹的目录结构，如：
        add_download("s3://test/123/456/789/test.txt", "/root/"): 不写相对基础路径时，
        下载的结果路径为"/root/test.txt", 如果有多个文件重名，则会相互覆盖
        add_download("s3://test/123/456/789/test.txt", "/root/", "s3://test/123/"): 指定了相对基础路径时，
        下载的结果路径为"/root/456/789/test.txt"，保留了相对基础路径之后的目录结构
        :return:
        """
        type_name = self.config.get_file_type(from_path)
        if not (type_name == "s3" or type_name == "http"):
            raise Exception("下载路径%s不正确，只支持对象存储和http类型!" % from_path)
        # if not to_path.endswith("/"):
        #     raise Exception("文件夹路径%s必须以\"/\"号结尾!" % to_path)
        if not to_path.endswith("/"):
            to_path = "%s/" % to_path
        if base_path:
            if not from_path.startswith(base_path):
                raise Exception("base_path %s与from_path %s不能匹配!" % (base_path, from_path))
        env_dist = os.environ
        tool_id = ""
        if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
            tool_id = self.config.current_tool_id
        data = filetrans_pb2.FileTrans(
            workflow_id=self.config.current_workflow_id,
            tool_id=tool_id,
            host=get_hostname(),
            processid=int(os.getpid()),
            upload=False,
            frompath=from_path,
            topath=to_path,
            basepath=base_path if base_path else "",
            usecache=True
        )
        self.logger.info("datadata:{}-{}".format(self.config.current_workflow_id, tool_id))
        self._add_file_trans(data)
        # if from_path.endswith("/"):
        #     if type_name == "s3":
        #         file_list = _get_s3_dir_list(from_path)
        #         for f in file_list:
        #             cache_path = _get_s3_cache_path(f)
        #             if base_path:
        #                 target = os.path.join(to_path, os.path.relpath(f, base_path))
        #             else:
        #                 target = os.path.join(to_path, os.path.basename(f))
        #             key = _get_s3_key(f)
        #             tag = None
        #             try:
        #                 tag = xattr.getxattr(cache_path, "user.etag")
        #             except IOError:
        #                 pass
        #             if os.path.isfile(cache_path) and key and tag and tag == key.etag.strip('"'):
        #                 self.logger.info("文件%s站到本地缓存，跳过下载..." % f)
        #                 if os.path.exists(target):
        #                     if os.path.islink(target):
        #                         os.remove(target)
        #                     elif os.path.isdir(target):
        #                         shutil.rmtree(target)
        #                     else:
        #                         os.remove(target)
        #                 dir_name = os.path.dirname(target)
        #                 if not os.path.exists(dir_name):
        #                     os.makedirs(dir_name)
        #                 os.link(cache_path, target)
        #             else:
        #                 self._files[f] = (target, cache_path)
        #     else:
        #         file_list = HttpDir(self.config.convert_path_to_http(from_path)).get_files()
        #         for f in file_list:
        #             cache_path = _get_http_cache_path(f)
        #             if base_path:
        #                 base_path = self.config.convert_path_to_http(base_path)
        #                 target = os.path.join(to_path, os.path.relpath(f, base_path))
        #             else:
        #                 target = os.path.join(to_path, os.path.basename(f))
        #             info = _get_http_info(f)
        #             tag = None
        #             try:
        #                 tag = xattr.getxattr(cache_path, "user.etag")
        #             except IOError:
        #                 pass
        #             if os.path.isfile(cache_path) and tag and tag == info["etag"].strip('"'):
        #                 self.logger.info("文件%s站到本地缓存，跳过下载..." % f)
        #                 if os.path.exists(target):
        #                     if os.path.islink(target):
        #                         os.remove(target)
        #                     elif os.path.isdir(target):
        #                         shutil.rmtree(target)
        #                     else:
        #                         os.remove(target)
        #                 dir_name = os.path.dirname(target)
        #                 if not os.path.exists(dir_name):
        #                     os.makedirs(dir_name)
        #                 os.link(cache_path, target)
        #             else:
        #                 self._files[f] = (target, cache_path)
        # else:
        #     if not exists(from_path):
        #         raise RunningError("文件%s不存在", (from_path,), code="008")
        #     if base_path:
        #         target = os.path.join(to_path, os.path.relpath(from_path, base_path))
        #     else:
        #         target = os.path.join(to_path, os.path.basename(from_path))
        #     real_path = self.config.convert_real_path(from_path)
        #     if type_name == "s3":
        #         cache_path = _get_s3_cache_path(from_path)
        #         key = _get_s3_key(from_path)
        #         tag = None
        #         try:
        #             tag = xattr.getxattr(cache_path, "user.etag")
        #         except IOError:
        #             pass
        #         if os.path.isfile(cache_path) and key and tag and tag == key.etag.strip('"'):
        #             self.logger.info("文件%s站到本地缓存，跳过下载..." % from_path)
        #             if os.path.exists(target):
        #                 if os.path.islink(target):
        #                     os.remove(target)
        #                 elif os.path.isdir(target):
        #                     shutil.rmtree(target)
        #                 else:
        #                     os.remove(target)
        #             dir_name = os.path.dirname(target)
        #             if not os.path.exists(dir_name):
        #                 os.makedirs(dir_name)
        #             os.link(cache_path, target)
        #         else:
        #             self._files[real_path] = (target, _get_s3_cache_path(real_path))
        #     else:
        #         cache_path = _get_http_cache_path(real_path)
        #         info = _get_http_info(real_path)
        #         tag = None
        #         try:
        #             tag = xattr.getxattr(cache_path, "user.etag")
        #         except IOError:
        #             pass
        #         if os.path.isfile(cache_path) and tag and tag == info["etag"].strip('"'):
        #             self.logger.info("文件%s站到本地缓存，跳过下载..." % real_path)
        #             if os.path.exists(target):
        #                 if os.path.islink(target):
        #                     os.remove(target)
        #                 elif os.path.isdir(target):
        #                     shutil.rmtree(target)
        #                 else:
        #                     os.remove(target)
        #             dir_name = os.path.dirname(target)
        #             if not os.path.exists(dir_name):
        #                 os.makedirs(dir_name)
        #             os.link(cache_path, target)
        #         else:
        #             self._files[real_path] = (target, _get_http_cache_path(real_path))

    def trans_file(self, file_iter):
        if not isinstance(file_iter, Iterable):
            raise Exception("file_iter必须是Iterable类型!" % file_iter)
        self._add_file_times += 1
        env_dist = os.environ
        port = self.config.wfm_port
        if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
            port = self.config.ntm_port
        try:
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = filetrans_pb2_grpc.FileTransServerStub(channel)
                responses = stub.TransFile(file_iter)
                if not responses.ok:
                    self.logger.error("添加下载失败,原因: %s" % responses.reason)
                    raise TransferError("添加下载失败,原因: %s" % responses.reason)
                self._add_file_times = 0
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._add_file_times > 3:
                self.logger.error("添加下载发生错误超过3次,退出运行: %s" % e)
                raise e
            else:
                self.logger.error("添加上传下载发生错误10秒后重试: %s" % e)
                gevent.sleep(20)
                self.trans_file(file_iter)

    def perform(self):
        """
        开始传输文件队列，并等待传输完成
        :return:
        """
        # for k, v in self._files.items():
        #     self.manager.add(k, v[0])
        # self.manager.wait(end=True)
        # for k, v in self._files.items():
        #     if os.path.exists(v[1]):
        #         if os.path.islink(v[1]):
        #             os.remove(v[1])
        #         elif os.path.isdir(v[1]):
        #             shutil.rmtree(v[1])
        #         else:
        #             os.remove(v[1])
        #     dir_name = os.path.dirname(v[1])
        #     if not os.path.exists(dir_name):
        #         os.makedirs(dir_name)
        #     if os.path.exists(k):
        #         if k != v[1]:
        #             os.link(k, v[1])
        #     elif os.path.exists(v[0]):
        #         if v[0] != v[1]:
        #             os.link(v[0], v[1])
        # self._files = {}
        env_dist = os.environ
        tool_id = ""
        if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
            tool_id = self.config.current_tool_id
        data = filetrans_pb2.WorkflowOrTool(
            workflow_id=self.config.current_workflow_id,
            tool_id=tool_id,
            host=get_hostname(),
            processid=int(os.getpid()),
        )
        self._wait_grpc(data)

    def _wait_grpc(self, data):
        self._wait_times += 1
        env_dist = os.environ
        port = self.config.wfm_port
        if "current_mode" in env_dist.keys() and env_dist["current_mode"] == "tool":
            port = self.config.ntm_port
        try:
            with grpc.insecure_channel('localhost:%s' % port) as channel:
                stub = filetrans_pb2_grpc.FileTransServerStub(channel)
                responses = stub.Wait(data)
                for progres in responses:
                    self.logger.info("文件传输进度:共计 %s个(%s), 其中成功%s个(%s),缓存 %s个(%s) ,"
                                     "失败%s个(%s),速度: %s/s" %
                                     (progres.total, friendly_size(progres.total_size), progres.trans_number,
                                      friendly_size(progres.trans_size), progres.cache_number,
                                      friendly_size(progres.cache_size),
                                      progres.error_number, friendly_size(progres.error_size),
                                      friendly_size(progres.speed)))
                    if progres.end:
                        self.logger.info("文件传输完成")
                        if progres.total == 0:
                            self.logger.error("没有可传输的文件,请检查源路径是否为空: %s" % progres.info)
                        if progres.error_number > 0:
                            self.logger.error("文件传输错误详情: %s" % progres.info)
                self._wait_times = 0
        except grpc.RpcError as rpc_error:
            self.logger.error('RPC通信错误 %s: %s' % (rpc_error.code(), rpc_error))
            raise rpc_error
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            sys.stdout.flush()
            if self._wait_times > 3:
                self.logger.error("等待传输完成发生错误超过3次,退出运行: %s" % e)
                raise e
            else:
                self.logger.error("等待传输完成发生错误10秒后重试: %s" % e)
                gevent.sleep(20)
                self._wait_grpc(data)

# class TransferFailedError(Exception):
#     """
#     文件传输失败时触发
#     """
#     def __init__(self, value):
#         Exception.__init__(self, value)
#         self.value = value
#
#     def __str__(self):
#         return str(self.value)


# @singleton
# class TransferManager(object):
#     def __init__(self):
#         self.config = Config()
#         self._conn = None
#         self._max_threads = self.config.get_rgw_max_threads()
#         self.base_path = os.getcwd()
#         self.overwrite = True
#         self._start_time = None
#         self._count_bytes = 0
#         self.log = True
#         self.etag = True
#         self.process = None
#         self.queue = Queue()
#         self.failed_queue = Queue()
#         self.is_wait = Value("i", 0)
#         self.to_end = Value("i", 0)
#         self._threads = {}
#         self.files = {}
#         self.chunk_queue = Queue()
#         self.lock = Lock()
#         self.has_error = Value("i", 0)
#         self._has_start = False
#         self.logger = logging.getLogger("FileTransfer")
#         self.is_end = False
#         self.failed_retry = {}
#         self._error_files = {}
#         self.failed_files = []
#
#     def update(self, size):
#         with self.lock:
#             self._count_bytes += size
#
#     @property
#     def speed(self):
#         with self.lock:
#             if self._start_time:
#                 second = (datetime.datetime.now() - self._start_time).seconds
#                 if second == 0:
#                     second = 1
#                 data = "%s/s" % friendly_size(self._count_bytes*1.0/second)
#                 if second > 600:
#                     self._start_time = datetime.datetime.now()
#                     self._count_bytes = 0
#                 return data
#             else:
#                 return 0
#
#     def add(self, from_uri, to_uri=None):
#         self.is_end = False
#         has_warn = False
#         while self.is_wait.value == 1:
#             if not has_warn:
#                 self.logger.debug("Last wait has not finished, waiting ....")
#                 has_warn = True
#             gevent.sleep(1)
#         if isinstance(from_uri, unicode):
#             from_uri = from_uri.encode("utf-8")
#         if isinstance(to_uri, unicode):
#             to_uri = to_uri.encode("utf-8")
#
#         self.add_queue(from_uri, to_uri)
#
#     # def _join_process(self):
#     #     if self.is_end is True:
#     #         return "exit"
#     #     if self.process:
#     #         if not self.process.is_alive():
#     #             self.process.join()
#     #             self.process = None
#
#     def _check_queue(self):
#         if self.is_end is True:
#             self._has_start = False
#             return "exit"
#         with self.lock:
#             if self.process:
#                 if not self.process.is_alive():
#                     self.process.join()
#                     self.process = None
#                     if not self.queue.empty():
#                         self.process = gipc.start_process(self._start_process, args=())
#             else:
#                 if not self.queue.empty():
#                     self.process = gipc.start_process(self._start_process, args=())
#
#     def wait(self, end=False):
#         with self.lock:
#             if self.process:
#                 if not self.process.is_alive():
#                     self.process.join()
#                     self.process = None
#                     if not self.queue.empty():
#                         self.process = gipc.start_process(self._start_process, args=())
#             else:
#                 if not self.queue.empty():
#                     self.process = gipc.start_process(self._start_process, args=())
#             self.is_wait.value = 1
#             if end:
#                 self.to_end.value = 1
#         while self.is_wait.value == 1:
#             # if self.has_error.value == 1:
#             #     raise Exception("文件上传下载模块出现错误!")
#             if self.process:
#                 if not self.process.is_alive():
#                     self.process.join()
#                     self.process = None
#                     break
#             else:
#                 if self.queue.empty() and self.chunk_queue.empty():
#                     break
#                 else:
#                     self.process = gipc.start_process(self._start_process, args=())
#             gevent.sleep(1)
#         self.logger.debug("等待任务结束已完成")
#         if end:
#             while True:
#                 # if self.has_error.value == 1:
#                 #     raise Exception("文件上传下载模块出现错误!")
#                 if self.process and self.process.is_alive():
#                     gevent.sleep(1)
#                 else:
#                     if self.process:
#                         self.process.join()
#                         self.process = None
#                     if self.queue.empty() and self.chunk_queue.empty():
#                         self.is_end = True
#                         self._has_start = False
#                         self.logger.debug("所有任务完成，结束等待状态和传输进程...." )
#                         break
#                     else:
#                         self.process = gipc.start_process(self._start_process, args=())
#         self.is_wait.value = 0
#         failed_list = self.get_failed()
#         retry_list = []
#         for f in failed_list:
#             if f in self.failed_retry.keys():
#                 self.failed_retry[f] += 1
#                 if self.failed_retry[f] <= 3:
#                     retry_list.append(f)
#                 else:
#                     self.failed_files.append(f)
#             else:
#                 self.failed_retry[f] = 0
#                 retry_list.append(f)
#         if len(retry_list) > 0:
#             self.logger.debug("文件%s下载失败.开始重试...." % retry_list)
#             for f in retry_list:
#                 self.add(f, self._error_files[f])
#             if not (self.process and self.process.is_alive()):
#                 self.process = gipc.start_process(self._start_process, args=())
#             return self.wait(end)
#         if len(self.failed_files):
#             raise RunningError("文件传输失败:%s", (";\t".join(failed_list),), '010')
#         # if self.has_error.value == 1:
#         #     raise Exception("文件传输失败:%s" % failed_list)
#
#     def get_failed(self):
#         file_list = []
#         while True:
#             try:
#                 data = self.failed_queue.get_nowait()
#             except Empty:
#                 break
#             else:
#                 if data:
#                     file_list.append(data[0])
#                     self._error_files[data[0]] = data[1]
#         return list(set(file_list))
#
#     def add_queue(self, from_path, to_path=None):
#         if from_path.endswith("/"):
#             raise Exception("不支持文件夹的上传下载!")
#         if not exists(from_path):
#             raise RunningError("文件%s不存在", (from_path,), '008')
#         type1 = self.config.get_file_type(from_path)
#         if to_path:
#             type2 = self.config.get_file_type(to_path)
#             if type1 == "local":
#                 if type2 != "s3":
#                     if type2 == "http":
#                         raise Exception("目标地址%s不支持上传!" % to_path)
#                     else:
#                         raise Exception("不支持的传输类型: %s, %s！" % (from_path, to_path))
#             elif type1 == "s3" or type1 == "http":
#                 if type2 != "local":
#                     raise Exception("不支持的传输类型: %s, %s！" % (from_path, to_path))
#             else:
#                 raise Exception("不支持的传输类型: %s, %s！" % (from_path, to_path))
#         else:
#             if not (type1 == "s3" or type1 == "http"):
#                 raise Exception("不支持的传输类型: %s, %s！" % (from_path, to_path))
#         while self.is_wait.value == 1:
#             gevent.sleep(1)
#         self.is_wait.value = 0
#         if self.chunk_queue.empty():
#             self._start_time = datetime.datetime.now()
#             self._count_bytes = 0
#         while self.queue.full():
#             gevent.sleep(1)
#         self.queue.put((from_path, to_path))
#         self.logger.debug("添加文件%s到队列..." % from_path)
#         self.to_end.value = 0
#         self.is_end = False
#         self.run()
#
#     def run(self):
#         if self.process and self.process.is_alive():
#             return
#         if not self._has_start:
#             Watcher().add(self._check_queue, 3)
#             self._has_start = True
#         if self.process:
#             if self.process.is_alive():
#                 return
#             else:
#                 self.process.join()
#                 self.process = None
#                 if not self.queue.empty():
#                     self.process = gipc.start_process(self._start_process, args=())
#         else:
#             self.process = gipc.start_process(self._start_process, args=())
#
#     def _start_process(self):
#         try:
#             self.logger.debug("开始启动新进程%s ...." % os.getpid())
#             old_name = setproctitle.getproctitle()
#             setproctitle.setproctitle("FileTransfer %s" % old_name)
#             while True:
#                 for i, f in self.files.items():
#                     if f.finish:
#                         self.files.pop(i)
#                 try:
#                     data = self.queue.get_nowait()
#                 except Empty:
#                     # with self.lock:
#                     all_end_work = True
#                     for thread in self._threads.values():
#                         if thread.is_work:
#                             all_end_work = False
#                         if not thread.is_alive():
#                             thread.join()
#                             self._threads.pop(str(id(thread)))
#                     all_file_end = True
#                     for i, f in self.files.items():
#                         if f.is_end:
#                             if f.finish:
#                                 self.files.pop(i)
#                         else:
#                             all_file_end = False
#                     if self.is_wait.value == 1 and all_file_end is True and all_end_work is True:
#                         self.is_wait.value = 0
#                     # if all_end_work is True and all_file_end is False:
#                     #     self._start_transfer_thread(1)
#                     if len(self._threads.keys()) == 0 or all_end_work:
#                         if self.queue.empty() and self.chunk_queue.empty() and all_file_end is True:
#                             if self.to_end.value == 1:
#                                 self.logger.debug("所有文件队列传输完成，退出进程%s!!" % os.getpid())
#                                 self.is_wait.value = 0
#                                 break
#                             else:
#                                 self.logger.debug("所有传输线程结束请队列为空，进程%s等待30秒...!" % os.getpid())
#                                 time.sleep(30)
#                                 if self.queue.empty() and self.chunk_queue.empty():
#                                     self.logger.debug("等待结束后队列仍然为空，退出进程%s!" % os.getpid())
#                                     self.is_wait.value = 0
#                                     break
#                                 else:
#                                     self.logger.debug("发现队列有数据，进程%s继续运行...!" % os.getpid())
#                                     for thread in self._threads.values():
#                                         if not thread.is_alive():
#                                             thread.join()
#                                             self._threads.pop(str(id(thread)))
#                                     self._start_transfer_thread(self._max_threads - len(self._threads.keys()))
#                                     continue
#                         else:
#                             # if self.has_error.value == 1:
#                             #     raise Exception("传输过程发生错误!")
#                             if self.queue.empty() and self.chunk_queue.empty():
#                                 all_file_end = True
#                                 for i, f in self.files.items():
#                                     if not f.is_end:
#                                         all_file_end = False
#                                         self.failed_queue.put((f.from_path, f.to_path))
#                                 if all_file_end is False:
#                                     self.is_wait.value = 0
#                                     raise Exception("传输过程发生不可预计的错误!")
#                             else:
#                                 for thread in self._threads.values():
#                                     if not thread.is_alive():
#                                         thread.join()
#                                         self._threads.pop(str(id(thread)))
#                                 self._start_transfer_thread(self._max_threads - len(self._threads.keys()))
#
#                     time.sleep(1)
#                 else:
#                     if data:
#                         try:
#                             if len(self._threads.keys()) < self._max_threads:
#                                 # if self.has_error.value == 1:
#                                 #     raise Exception("传输过程发生错误!")
#                                 self._start_transfer_thread(self._max_threads - len(self._threads.keys()))
#                             if self.config.get_file_type(data[0]) == "local":
#                                 if self.config.get_file_type(data[1]) == "s3":
#                                     file_transfer = S3FileTransfer(self)
#                                 elif self.config.get_file_type(data[1]) == "http":
#                                     raise Exception("http文件类型不支持上传!")
#                                 else:
#                                     raise Exception("不支持的传输类型: %s, %s！" % (data[0], data[1]))
#                             elif self.config.get_file_type(data[1]) == "local":
#                                 if self.config.get_file_type(data[0]) == "s3":
#                                     file_transfer = S3FileTransfer(self)
#                                 elif self.config.get_file_type(data[0]) == "http":
#                                     file_transfer = HttpFileTransfer(self)
#                                 else:
#                                     raise Exception("不支持的传输类型: %s, %s！" % (data[0], data[1]))
#                             else:
#                                 raise Exception("不支持的传输类型: %s, %s！" % (data[0], data[1]) )
#                             self.files[str(id(file_transfer))] = file_transfer
#                             file_transfer.set_uri(data[0], data[1])
#                             self.logger.debug("开始传输文件%s 到 %s, 大小: %s, Etag: %s, chunks %s" %
#                                               (data[0], data[1], file_transfer.size, file_transfer.etag,
#                                                file_transfer.chunks))
#                             # file_transfer.conn.close()
#                             file_transfer.start()
#                             all_end_work = True
#                             for thread in self._threads.values():
#                                 if thread.is_work:
#                                     all_end_work = False
#                                 if not thread.is_alive():
#                                     thread.join()
#                                     self._threads.pop(str(id(thread)))
#                             if all_end_work:
#                                 self._start_transfer_thread(self._max_threads - len(self._threads.keys()))
#                             time.sleep(0.1)
#                         except Exception, e:
#                             exstr = traceback.format_exc()
#                             print exstr
#                             print e
#                             sys.stdout.flush()
#                 gc.collect()
#         except Exception, e:
#             exstr = traceback.format_exc()
#             print exstr
#             print e
#             # self.has_error.value = 1
#             sys.stdout.flush()
#
#     def _start_transfer_thread(self, number):
#         for i in xrange(number):
#             thread = LoopThread(self)
#             thread.start()
#             self._threads[str(id(thread))] = thread


# class LoopThread(threading.Thread):
#     def __init__(self,  manager):
#         super(LoopThread, self).__init__()
#         self.manager = manager
#         self.is_work = False
#         self._last_work = datetime.datetime.now()
#         self.logger = None
#
#     def run(self):
#         super(LoopThread, self).run()
#         self.logger = logging.getLogger("FileTransfer(thread %s)" % self.ident)
#         self.logger.debug("开始新线程运行...")
#         while True:
#             try:
#                 data = self.manager.chunk_queue.get_nowait()
#             except Empty:
#                 self.is_work = False
#                 if (datetime.datetime.now() - self._last_work).seconds > 30:
#                     self.logger.info("30秒钟没有检测到新任务，退出运行...")
#                     break
#                 if self.manager.to_end.value == 1 and self.manager.chunk_queue.empty() and self.manager.queue.empty():
#                     self.logger.info("队列为空，结束运行...")
#                     break
#                 time.sleep(1)
#             else:
#                 try:
#                     if data:
#                         self._last_work = datetime.datetime.now()
#                         self.is_work = True
#                         file_transfer = self.manager.files[str(data[0])]
#                         if file_transfer.type == 0:
#                             # thread = S3DownloadThread(file_transfer, data[1])
#                             file_transfer.download(data[1])
#                         else:
#                             # thread = S3UploadThread(file_transfer, data[1])
#                             file_transfer.upload(data[1])
#                 except Exception:
#                     exstr = traceback.format_exc()
#                     print exstr
#                     sys.stdout.flush()
