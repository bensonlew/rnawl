# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import importlib
import re
import os
import sys
import json
from datetime import datetime, date
import inspect
import ctypes
import setproctitle
from bson.objectid import ObjectId
# from ctypes import Structure, c_char, c_int, addressof, memmove,  POINTER, sizeof, byref
# import time
# from .exceptions import MaxLengthError
# import urllib2
# import urllib
from ..config import Config
# import random
import hashlib
# import gevent
# import types
# from .exceptions import OptionError, FileError, RunningError, CodeError


def get_from_friendly_size(size):
    pattern = re.compile(r'([\d.]+)([gmk])', re.I)
    match = re.match(pattern, size)
    normal_size = 0
    if match:
        unit = match.group(2)
        if unit.upper() == "G":
            normal_size = int(match.group(1)) * 1024 * 1024 * 1024
        elif unit.upper() == "M":
            normal_size = int(match.group(1)) * 1024 * 1024
        elif unit.upper() == "K":
            normal_size = int(match.group(1)) * 1024
    else:
        pattern = re.compile(r'([\d.]+)', re.I)
        match = re.match(pattern, size)
        if match:
            normal_size = int(match.group(1))
    return normal_size


def _async_raise(tid, exctype):
    """raises the exception, performs cleanup if needed"""
    if not inspect.isclass(exctype):
        exctype = type(exctype)
    res = ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, ctypes.py_object(exctype))
    if res == 0:
        raise ValueError("invalid thread id")
    elif res != 1:
        # """if it returns a number greater than one, you're in trouble,
        # and you should call it again with exc=NULL to revert the effect"""
        ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, None)
        raise SystemError("PyThreadState_SetAsyncExc failed")


def stop_thread(thread):
    _async_raise(thread.ident, SystemExit)


def get_clsname_form_path(path, tp="Agent"):
    """
    从path中获取ClassName

    :param path:  string 模块路径
    :param tp:  string 种类Agent Tool,Module,Workflow
    :return:  string
    """
    path = re.sub(r'[^a-zA-Z0-9._]', '', path)
    name = path.split(".").pop()
    ll = name.split("_")
    ll.append(tp)
    l1 = [el.capitalize() for el in ll]
    return "".join(l1)


def load_class_by_path(path, tp="Agent"):
    """
    根据Path找到对应的类，并返回类对象

    :param path: string 模块路径
    :param tp:  string 种类Agent Tool,Module,Workflow
    :return: 对应的class对象
    """
    dir_name = {
        "Agent": "mbio.tools.",
        "Tool": "mbio.tools.",
        "Module": "mbio.modules.",
        "Workflow": "mbio.workflows.",
        "Package": "mbio.packages.",
        "File": "mbio.files."
    }
    class_name = get_clsname_form_path(path, tp)
    module_name = dir_name[tp] + path
    # try:
    #    sys.modules[module_name]
    # except KeyError:
    if module_name in sys.modules.keys():
        imp = sys.modules[module_name]
    else:
        imp = importlib.import_module(module_name)
    # else:
    #     del sys.modules[module_name]
    #     imp = importlib.import_module(module_name)
    if tp == "Package":
        return imp
    return getattr(imp, class_name)


def get_classpath_by_object(obj):
    """
    通过对象获取类导入路径  包括自定义扩展的 Agent Tool Module  Workflow File子类属性
    :param obj:
    :return:
    """
    class_name = str(type(obj))
    m = re.match(r"<class\s\'(.*)\'>", class_name)
    class_name = m.group(1)
    paths = class_name.split(".")
    paths.pop()
    return ".".join(paths)


def daemonize(stdout='/dev/null', stderr='/dev/null'):

    try:
        pid = os.fork()
        if pid > 0:
            sys.exit(0)
    except OSError as e:
        sys.stderr.write("fork #1 failed: (%d) %s\n" % (e.errno, e.strerror))
        sys.exit(1)

    # 从母体环境脱离
    # os.chdir("/")
    # os.umask(0)
    os.setsid()
    # 执行第二次fork
    try:
        pid = os.fork()
        if pid > 0:
            sys.exit(0)  # second parent out
    except OSError as e:
        sys.stderr.write("fork #2 failed: (%d) %s]n" % (e.errno, e.strerror))
        sys.exit(1)

    # 进程已经是守护进程了，重定向标准文件描述符
    for f in sys.stdout, sys.stderr:
        f.flush()

    so = open(stdout, 'a+')
    se = open(stderr, 'a+', 0)
    os.dup2(so.fileno(), sys.stdout.fileno())
    os.dup2(se.fileno(), sys.stderr.fileno())
    # so.close()
    # se.close()


def get_hostname():
    sys_name = os.name
    if sys_name == 'nt':
        host_name = os.getenv('computername')
        return host_name.replace('.local', '')
    elif sys_name == 'posix':
        with os.popen('echo $HOSTNAME') as f:
            host_name = f.readline()
            return host_name.strip('\n').replace('.local', '')
    else:
        return 'Unkwon hostname'


hostname = get_hostname()


class CJsonEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(obj, date):
            return obj.strftime('%Y-%m-%d')
        elif isinstance(obj, ObjectId):
            return str(obj)
        else:
            return json.JSONEncoder.default(self, obj)


def change_process_tile(tile):
    """
    更改当前进程名称
    :param tile: String 进程名称
    :return:
    """
    setproctitle.setproctitle(tile)


def filter_error_info(info):
    """
    过滤返回前端的错误信息，避免出现系统工作路径等敏感信息

    :param info: String 错误信息
    :return:
    """
    suberror = re.compile(r'/.+/.+/.+/')
    info = suberror.sub('', str(info))
    info = info.replace('\"', ' ').replace('\'', ' ')
    return info


def friendly_size(size):
    gb = 1024 * 1024 * 1024.0
    mb = 1024 * 1024.0
    kb = 1024.0
    if size >= gb:
        return "%.2fG" % (size / gb)
    elif size >= mb:
        return "%.2fM" % (size / mb)
    elif size >= kb:
        return "%.2fk" % (size / kb)
    else:
        return "%s" % size


def pretty_date(diff):
    """
    Get a datetime object or a int() Epoch timestamp and return a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc
    """
    second_diff = diff.seconds
    day_diff = diff.days

    if day_diff < 0:
        return '刚刚'
    if day_diff == 0:
        if second_diff < 10:
            return "刚刚"
        if second_diff < 60:
            return str(second_diff) + "秒前"
        if second_diff < 120:
            return "1分钟前"
        if second_diff < 3600:
            return str(second_diff / 60) + "分钟前"
        if second_diff < 7200:
            return "1小时前"
        if second_diff < 86400:
            return str(second_diff / 3600) + "小时前"
    if day_diff == 1:
        return "昨天"
    if day_diff < 7:
        return str(day_diff) + "天前"
    if day_diff < 31:
        return str(day_diff / 7) + "星期前"
    if day_diff < 365:
        return str(day_diff / 30) + "月前"
    return str(day_diff / 365) + "年前"


def friendly_time(time_along):
    if time_along < 0:
        return ""
    elif time_along < 60:
        return str(time_along) + "秒"
    elif time_along < 3600:
        s = str(time_along/60) + "分"
        second = time_along % 60
        if second > 0:
            s += "%s秒" % second
        return s
    else:
        s = str(time_along / 3600) + "小时"
        second = time_along % 3600
        if second > 60:
            m = str(second / 60)
            if int(m) > 0:
                s += "%s分" % m
        return s


# class MessageData(Structure):
#     _fields_ = [('action', c_char*50), ('id', c_char*1024), ("data", c_char*2048)]
#
#     def __str__(self):
#         return "<MessageData: action=%s, id=%s, other=%s, addr=%ld>" % \
#                (self.action, self.id, self.other, addressof(self))
#
#
# class LogData(Structure):
#     _fields_ = [('action', c_char*50), ("data", c_char*20480), ("chunks", c_int)]
#
#     def __str__(self):
#         return "<MessageData: action=%s, data=%s, addr=%ld>" % \
#                (self.action,  self.data, addressof(self))
#
#
# class ActionData(Structure):
#     _fields_ = [('action', c_char*50), ('id', c_char*50), ("timestamp", c_int)]
#
#     def __str__(self):
#         return "<ActionData: action=%s, id=%s, timestamp=%s, addr=%ld>" % \
#                (self.action, self.data, self.timestamp, addressof(self))


# def copier_factory(typ):
#     def f(a, b):
#         memmove(a, b, sizeof(typ))
#     f.argtypes = [POINTER(typ), POINTER(typ)]
#     return f
#
#
# def copy_msg(a):
#     b = MessageData()
#     copy_func = copier_factory(MessageData)
#     copy_func(byref(b), byref(a))
#     return b
#
#
# def copy_action(a):
#     b = ActionData()
#     copy_func = copier_factory(ActionData)
#     copy_func(byref(b), byref(a))
#     return b
#
#
# def copy_log(a):
#     b = LogData()
#     copy_func = copier_factory(LogData)
#     copy_func(byref(b), byref(a))
#     return b
#
#
# def add_run_queue(queue, action, ids, msg=None):
#     # length = len(self.run_info_queue)
#     is_full = True
#     if len(action) > 50:
#         raise MaxLengthError("action %s 长度超过最长50个限制!" % action)
#     id_str = json.dumps(ids)
#     if len(id_str) > 1024:
#         raise MaxLengthError("ids %s 长度超过最长1024个限制!" % id_str)
#     msg_str = json.dumps(msg, cls=CJsonEncoder) if msg else "\"\""
#     if len(msg_str) > 2048:
#         raise MaxLengthError("msg %s 长度超过最长1024个限制!" % msg_str)
#     with queue.get_lock():
#         for i, v in enumerate(queue):
#             if v.action == "":
#                 is_full = False
#                 queue[i] = (action, id_str, msg_str)
#                 break
#     if is_full:
#         print "状态队列已满，请稍后...."
#         sys.stdout.flush()
#         time.sleep(1)
#         return add_run_queue(queue, action, id_str, msg_str)
#
#
# def add_log_queue(queue, action, data):
#     length = len(queue)
#     is_full = True
#     size = 20480
#     if len(action) > 50:
#         raise MaxLengthError("action %s 长度超过最长50个限制!" % action)
#     json_data = json.dumps(data, cls=CJsonEncoder) if data else '""'
#     json_len = len(json_data)
#     if json_len > size * 1000:
#         raise MaxLengthError("data 长度超过最长10M系统不支持!")
#     with queue.get_lock():
#         for i, v in enumerate(queue):
#             if queue[i].action == "":
#                 is_full = False
#                 if json_len > size * (length - i):
#                     print "日志队列已满，请稍后...."
#                     sys.stdout.flush()
#                     is_full = True
#                 elif json_len > size:
#                     if json_len % size > 0:
#                         chunks = json_len // size + 1
#                     else:
#                         chunks = json_len // size
#                     for m in xrange(chunks):
#                         queue[i+m] = (action, json_data[(m*size):(m*size+size)], chunks)
#                     break
#                 else:
#                     queue[i] = (action, json_data, 0)
#                     break
#     if is_full:
#         print "日志队列已满，请稍后...."
#         sys.stdout.flush()
#         time.sleep(1)
#         return add_log_queue(queue, action, data)
#
#
# def add_action_queue(queue, action, wid):
#     time_stamp = int(time.time())
#     is_full = True
#     if len(action) > 50:
#         raise MaxLengthError("action %s 长度超过最长50个限制!" % action)
#     if len(wid) > 50:
#         raise MaxLengthError("wid %s 长度超过最长50个限制!" % wid)
#     with queue.get_lock():
#         for i, v in enumerate(queue):
#             if v.id == wid and v.action == action:
#                 return
#         for i, v in enumerate(queue):
#             if queue[i].action == "":
#                 is_full = False
#                 queue[i] = (action, wid, time_stamp)
#                 break
#             if time_stamp - queue[i].timestamp > 200:
#                 is_full = False
#                 queue[i] = (action, wid, time_stamp)
#                 break
#     if is_full:
#         print "指令队列已满，请稍后...."
#         sys.stdout.flush()
#         time.sleep(1)
#         return add_action_queue(queue, action, wid)
#
#
# def get_run_queue(queue):
#     msg_list = []
#     with queue.get_lock():
#         for i, v in enumerate(queue):
#             if queue[i].action == "":
#                 break
#             msg_list.append(copy_msg(queue[i]))
#             queue[i] = ("", "", "")
#     data = []
#     for msg in msg_list:
#         data.append((msg.action, json.loads(msg.id), json.loads(msg.data)))
#     return data
#
#
# def get_log_queue(queue):
#     log_list = []
#     with queue.get_lock():
#         for i, v in enumerate(queue):
#             if queue[i].action == "":
#                 break
#             log_list.append(copy_log(queue[i]))
#             queue[i] = ("", "", 0)
#     data = []
#     x = 0
#     while x < len(log_list):
#
#         big_data = ""
#         if log_list[x].chunks > 1:
#             chunks = log_list[x].chunks
#             action = log_list[x].action
#             for m in xrange(chunks):
#                 big_data += log_list[x].data
#                 x += 1
#             data.append((action, json.loads(big_data)))
#         else:
#             data.append((log_list[x].action, json.loads(log_list[x].data)))
#             x += 1
#     return data
#
#
# def get_action_queue(queue, wid):
#     with queue.get_lock():
#         for i, v in enumerate(queue):
#             if v.id == wid:
#                 action = v.action
#                 queue[i] = ("", "", 0)
#                 return action
#     return None
#
#
# class StateData(Structure):
#     _fields_ = [("id", c_char*500), ('state', c_char*50), ('jobid', c_char*50), ('host', c_char*50),
#                 ("data", c_char*2048),  ("version", c_int), ("timestamp", c_int)]
#
#     def __str__(self):
#         return "< StateData: id=%s, state=%s, jobid=%s, host=%s, version=%s, data=%s, timestamp=%s, addr=%ld>" % \
#                (self.id, self.state,  self.jobid, self.host, self.version, self.data, self.timestamp, addressof(self))
#
#     def to_dict(self):
#         return {
#             "id": self.id,
#             "state": self.state,
#             "jobid": self.jobid,
#             "host": self.host,
#             "data": json.loads(self.data),
#             "version": self.version
#         }
#
#
# def copy_state(a):
#     b = StateData()
#     copy_func = copier_factory(StateData)
#     copy_func(byref(b), byref(a))
#     return b
#
#
# class CallbackData(Structure):
#     _fields_ = [("id", c_char*500), ('action', c_char*50), ("data", c_char*2048),  ("version", c_int)]
#
#     def __str__(self):
#         return "< CallbackData: id=%s, action=%s,  data=%s,  version=%s, addr=%ld>" % \
#                (self.id, self.action,  self.data, self.version, addressof(self))
#
#     def to_dict(self):
#         return {
#             "id": self.id,
#             "action": self.action,
#             "data": json.loads(self.data),
#             "version": self.version
#         }
#
#
# def copy_callback(a):
#     b = CallbackData()
#     copy_func = copier_factory(CallbackData)
#     copy_func(byref(b), byref(a))
#     return b


def link(source, target, exclude=None, include_source_dir=True):
    """
    将当期的文件路径建立硬链接到目标路径

    :param source: 需要建立硬链接的源路径，可以是文件或文件夹

    :param target: 目标路径文件夹，源路径下的所有文件将被链接到此目录下，如果此目录不存在，将自动创建。
    如果存在，此目录下冲突的文件/文件夹将被删除并替换。

    :param exclude: list列表，其元素必须为正则表达式,相对于于source路径，只要匹配上了就会被排除链接，默认为None

    :param include_source_dir: 当source为文件夹时，是否需要在target下创建source指向的文件夹名。默认: True

    :return: None
    """
    if exclude is not None:
        if not isinstance(exclude, list):
            raise Exception("exclude必须为list数组，每个元素为相对路径的正则表达式")
    source = os.path.abspath(source)
    if not os.path.exists(source):
        raise Exception("source文件%s不存在" % source)
    # basename = os.path.basename(from_path)
    target = os.path.abspath(target)
    if not os.path.exists(target):
        os.makedirs(target)
    else:
        if not os.path.isdir(target):
            raise Exception("target文件%s不是文件夹" % target)

    if os.path.isdir(source):
        for root, dirs, files in os.walk(source):
            find_exclude_dir = False
            for r_rule in exclude:
                r_path = os.path.relpath(root, source)
                pattern = re.compile(r_rule)
                match = pattern.match(r_path)
                if match:
                    find_exclude_dir = True
            if find_exclude_dir:
                continue   # 匹配到文件夹跳过
            if include_source_dir:
                rel_path = os.path.relpath(root, os.path.dirname(source))
            else:
                rel_path = os.path.relpath(root, source)
            for i_file in files:
                find_exclude = False
                for r_rule in exclude:
                    r_path = os.path.relpath(os.path.join(root, i_file), source)
                    pattern = re.compile(r_rule)
                    match = pattern.match(r_path)
                    if match:
                        find_exclude = True
                if find_exclude:
                    continue  # 匹配到文件跳过
                i_file_path = os.path.join(root, i_file)
                if os.path.islink(i_file_path):
                    real_path = os.path.realpath(i_file_path)
                    if os.path.exists(real_path):
                        dir_path = os.path.join(target, rel_path)
                        if not os.path.exists(dir_path):
                            os.makedirs(dir_path)
                        file_path = os.path.join(dir_path, i_file)
                        if os.path.exists(file_path):
                            os.remove(file_path)
                        os.link(real_path, file_path)
                else:
                    dir_path = os.path.join(target, rel_path)
                    if not os.path.exists(dir_path):
                        os.makedirs(dir_path)
                    file_path = os.path.join(dir_path, i_file)
                    if os.path.exists(file_path):
                        os.remove(file_path)
                    os.link(i_file_path, file_path)
    else:
        if os.path.islink(source):
            real_path = os.path.realpath(source)
            if not os.path.exists(real_path):
                raise Exception("源文件%s是一个无效的软链接!" % source)
            os.link(real_path, os.path.join(target, os.path.basename(source)))
        else:
            os.link(source, os.path.join(target, os.path.basename(source)))


def get_signature(nonce, timestamp):
    key = Config().WEB_INTERFACE_KEY
    x_list = [key, timestamp, nonce]
    x_list.sort()
    sha1 = hashlib.sha1()
    map(sha1.update, x_list)
    return sha1.hexdigest()


# def get_web_data(params, logger, method="post", api="batch"):
#     http_handler = urllib2.HTTPHandler()
#     https_handler = urllib2.HTTPSHandler()
#     opener = urllib2.build_opener(http_handler, https_handler)
#     urllib2.install_opener(opener)
#     data = urllib.urlencode(params)
#     config = Config()
#     nonce = str(random.randint(1000, 10000))
#     timestamp = str(int(time.time()))
#     signature = {
#         "client": config.WEB_INTERFACE_CLIENT,
#         "nonce": nonce,
#         "timestamp": timestamp,
#         "signature": get_signature(nonce, timestamp)
#     }
#     signature = urllib.urlencode(signature)
#     url = "http://%s/%s" % (config.WEB_INTERFACE_DOMAIN, api)
#     if method == "get":
#         if "?" in url:
#             url = "%s&%s&%s" % (url, signature, data)
#         else:
#             url = "%s?%s&%s" % (url, signature, data)
#         logger.debug("get data to url %s ...\n\n" % url)
#         request = urllib2.Request(url)
#     else:
#         if "?" in url:
#             url = "%s&%s" % (url, signature)
#         else:
#             url = "%s?%s" % (url, signature)
#         logger.debug("post data to url %s ...\n\n" % url)
#         request = urllib2.Request(url, data)
#     try:
#         response = urllib2.urlopen(request)
#         sys.stdout.flush()
#         sys.stderr.flush()
#     except urllib2.HTTPError, e:
#         logger.warning("Web API状态出错: %s , 30秒后重试..." % e)
#         gevent.sleep(30)
#         return get_web_data(params, logger, method, api)
#     else:
#         return response.read()


def get_dir_size(dir_path):
    size = 0
    count = 0
    for root, dirs, files in os.walk(dir_path, followlinks=True):
        count = len(files)
        size += sum([os.path.getsize(os.path.realpath(os.path.join(root, name))) for name in files])
    return count, size


def get_error_str(error):
    try:
        error_data = json.loads(error)
    except Exception:
        return error
    else:
        if isinstance(error_data, dict) and "error_type" in error_data.keys():
            if "info" not in error_data.keys():
                error_data["info"] = "系统默认错误，尚未编码！"
            if "variables" not in error_data.keys():
                error_data["variables"] = None
            if "code" not in error_data.keys():
                error_data["code"] = "001"
            if error_data["error_type"] == "running":
                s = "RunningError, ErrorCode: %s, " % error_data["code"] + error_data["info"]
            elif error_data["error_type"] == "option":
                if "option" in error_data.keys():
                    s = "OptionError, 参数: %s, ErrorCode: %s, " % (error_data["option"], error_data["code"])\
                        + error_data["info"]
                else:
                    s = "OptionError, ErrorCode: %s, " % error_data["code"] + error_data["info"]
            elif error_data["error_type"] == "file":
                if "option" in error_data.keys():
                    if "file" in error_data.keys():
                        s = "FileError, 文件:%s, 参数: %s , ErrorCode: %s, " % (error_data["file"], error_data["option"],
                                                                            error_data["code"]) + error_data["info"]
                    else:
                        s = "FileError, 参数: %s, ErrorCode: %s, " % (error_data["option"], error_data["code"]) \
                            + error_data["info"]
                else:
                    if "file" in error_data.keys():
                        s = "FileError, 文件: %s ErrorCode: %s, " % (error_data["file"], error_data["code"])\
                            + error_data["info"]
                    else:
                        s = "FileError, ErrorCode: %s, " % error_data["code"] + \
                            error_data["info"]
            else:
                s = "CodeError, ErrorCode: %s, " % error_data["code"] + error_data["info"]
            return s
        else:
            return error
