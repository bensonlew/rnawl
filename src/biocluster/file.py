# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
import pycurl
from .config import Config
import traceback
import gevent
from pyquery import PyQuery as pq
import os
import re
from boto.s3.bucket import Bucket
from .core.exceptions import RunningError
import xattr
import shutil
import logging
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import tempfile

PY3 = sys.version_info[0] == 3


_http_file_info = {}
config = Config()
logger = logging.getLogger("boto")
logger.setLevel(logging.INFO)
logger1 = logging.getLogger("requests.packages.urllib3.connectionpool")
logger.setLevel(logging.INFO)


def exists(file_path):
    """
    判断文件或文件夹是否存在,
    当类型为对象存储或http时，不支持文件夹检测，否则可能返回错误的结果
    :param file_path:文件路径，支持所有路径格式
    :return:
    """
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    type_name = config.get_file_type(file_path)
    if type_name == "s3":
        if file_path.endswith("/"):
            raise Exception("文件路径%s类型为对象存储，不支持文件夹检测!" % file_path)
        else:
            key = _get_s3_key(file_path)
            if key:
                # key.close()
                return True
    elif type_name == "http":
        if file_path.endswith("/"):
            raise Exception("文件路径%s类型为http，不支持文件夹检测!" % file_path)
        else:
            info = _get_http_info(config.convert_path_to_http(file_path))
            return info["exits"]
    else:
        return os.path.exists(config.convert_real_path(file_path))


def list_dir(dir_path):
    """
    列出文件夹下所有文件的完整路径，支持所有类型文件
    注意：请勿对子文件数量过多的文件夹执行此方法(1000个子文件以上)，否则会导致程序卡住，严重时会导致系统崩溃
    :param dir_path: 文件夹路径必须以"/"结尾，否则会报错, 支持所有路径格式
    :return: list数组,所有级别子文件的完整路径，不含子文件夹
    """
    if not PY3:
        if isinstance(dir_path, unicode):
            dir_path = dir_path.encode("utf-8")
    if not dir_path.endswith("/"):
        raise Exception("文件夹路径必须以'/'号结尾!")
    type_name = config.get_file_type(dir_path)
    if type_name == "s3":
        return _get_s3_dir_list(dir_path)
    elif type_name == "http":
        return HttpDir(config.convert_path_to_http(dir_path)).get_files()
    else:
        file_list = []
        for root, dirs, files in os.walk(os.path.abspath(config.convert_real_path(dir_path)), topdown=False):
            for name in files:
                file_list.append(os.path.join(root, name))
            for name in dirs:
                file_list.append(os.path.join(root, name))
        return file_list


def get_reader(file_path):
    """
    获取文件对象的只读句柄
    注意：当文件类型为对象存储或http时，此方法会自动下载文件至缓存，打开文件时会导致进程被长时间阻塞，并且可能下载中断，
    因此此方法只适用于较小的文件(10M以内)，对于大文件请使用get_top_lines或者get_last_line方法获取文件部分内容
    :param file_path: 文件路径，不支持文件夹, 支持所有路径格式
    :return:
    """
    if file_path.endswith("/"):
        raise Exception("不支持文件夹读取!")
    file_path = download(file_path)
    return open(file_path, "r")


def get_top_lines(file_path, lines=10):
    """
    获取文件最开始部分指定的行数，lines值不能大于1000。此方法只针对文本文件，对于二进制文件将无法获取正确的内容
    此方法可用用于大文件内容的检测
    :param file_path: 支持所有路径格式
    :param lines: 读取行数，默认值10 ，不能大于1000
    :return:
    """
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    count = 0
    buf = StringIO()
    line_size = 1000
    file_size = getsize(file_path)
    while True:
        off_set = count * line_size * lines
        length = line_size * lines
        if off_set + length > file_size:
            length = file_size - off_set
            buf.write(_get_file_part(file_path, off_set, length))
            data = buf.getvalue().split("\n")
            if len(data) < lines:
                return [l.strip() for l in data]
            else:
                return [l.strip() for l in data[:lines]]
        else:
            count += 1
            buf.write(_get_file_part(file_path, off_set, length))
            data = buf.getvalue().split("\n")
            if len(data) > lines:
                return [l.strip() for l in data[:lines]]


def get_last_lines(file_path, lines=10):
    """
    获取文件结尾部分指定的最后行数，lines值不能大于1000。此方法只针对文本文件，对于二进制文件将无法获取正确的内容
    此方法可用用于大文件内容的检测
    :param file_path: 支持所有路径格式
    :param lines: 读取行数，默认值10 ，不能大于1000
    :return:
    """
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    count = 1
    buf = StringIO()
    line_size = 1000
    file_size = getsize(file_path)
    while True:
        off_set = file_size - count * line_size * lines
        length = line_size * lines
        if off_set < 0:
            off_set = 0
            length = file_size
            buf.write(_get_file_part(file_path, off_set, length))
            data = buf.getvalue().split("\n")
            if len(data) < lines:
                return [l.strip() for l in data]
            else:
                return [l.strip() for l in data[-lines:]]
        else:
            count += 1
            buf.write(_get_file_part(file_path, off_set, length))
            data = buf.getvalue().split("\n")
            if len(data) > lines:
                return [l.strip() for l in data[-lines:]]


def download(from_path, to_path=None):
    """
    单线程下载文件
    此方法只适用于单个小文件的下载(100M以下)，并会堵塞当前进程，不能用于大文件,不支持文件夹，也不支持并行
    :param from_path:,支持所有格式
    :param to_path: 本地路径，默认下载到缓存路径
    :return: 返回下载文件本地地址
    """
    if not PY3:
        if isinstance(from_path, unicode):
            from_path = from_path.encode("utf-8")
        if isinstance(to_path, unicode):
            to_path = to_path.encode("utf-8")
    if from_path.endswith("/"):
        raise Exception("不支持文件夹下载!")
    if to_path:
        to_path = os.path.abspath(to_path)
    else:
        if not PY3:
            temp_dir = tempfile.mkdtemp()
        else:
            temp_dir = tempfile.TemporaryDirectory().name
        print "temp_path:{}".format(temp_dir)
        to_path = os.path.join(temp_dir, os.path.basename(from_path))
    if exists(from_path):
        type_name = config.get_file_type(from_path)
        if type_name == "s3":
            # cache_path = _get_s3_cache_path(from_path)
            # tag = None
            # if os.path.isfile(cache_path):
            #     tag = xattr.getxattr(cache_path, "user.etag")
            key = _get_s3_key(from_path)
            # if os.path.isfile(cache_path) and key and tag and tag == key.etag.strip('"'):
            #     print("文件%s找到本地缓存文件%s，跳过下载..." % (from_path, cache_path))
            #     sys.stdout.flush()
            # else:
            #     _rm(cache_path)
            #     dir_name = os.path.dirname(cache_path)
            #     if not os.path.exists(dir_name):
            #         os.makedirs(dir_name)
            #     key.get_contents_to_filename(cache_path)
            #     # if key.etag:
            #     #     xattr.setxattr(cache_path, "user.etag", key.etag.strip('"'))
            # if to_path:
            #     _rm(to_path)
            #     dir_name = os.path.dirname(to_path)
            #     if not os.path.exists(dir_name):
            #         os.makedirs(dir_name)
            #     os.link(cache_path, to_path)
            # else:
            #     to_path = cache_path
            _rm(to_path)
            dir_name = os.path.dirname(to_path)
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
            key.get_contents_to_filename(to_path)
            key.close()
        elif type_name == "http":
            http_path = config.convert_path_to_http(from_path)
            _rm(to_path)
            dir_name = os.path.dirname(to_path)
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
            with open(to_path, "w") as f:
                if not PY3:
                    if isinstance(http_path, unicode):
                        http_path = http_path.encode("utf-8")
                c = pycurl.Curl()
                c.setopt(c.URL, http_path)
                c.setopt(c.FOLLOWLOCATION, 1)
                c.setopt(c.WRITEDATA, f)
                c.perform()
                c.close()
            # cache_path = _get_http_cache_path(http_path)
            # tag = None
            # if os.path.isfile(cache_path):
            #     try:
            #         tag = xattr.getxattr(cache_path, "user.etag")
            #     except:
            #         pass
            # if os.path.isfile(cache_path) and tag == _http_file_info[http_path]["etag"].strip('"'):
            #     print("文件%s找到本地缓存文件%s，跳过下载..." % (from_path, cache_path))
            #     sys.stdout.flush()
            # else:
            #     _rm(cache_path)
            #     dir_name = os.path.dirname(cache_path)
            #     if not os.path.exists(dir_name):
            #         os.makedirs(dir_name)
            #     with open(cache_path, "w") as f:
            #         if not PY3:
            #             if isinstance(http_path, unicode):
            #                 http_path = http_path.encode("utf-8")
            #         c = pycurl.Curl()
            #         c.setopt(c.URL, http_path)
            #         c.setopt(c.FOLLOWLOCATION, 1)
            #         c.setopt(c.WRITEDATA, f)
            #         c.perform()
            #         c.close()
            #     if _http_file_info[http_path]["etag"]:
            #         xattr.setxattr(cache_path, "user.etag", _http_file_info[http_path]["etag"].strip('"'))
            # if to_path:
            #     _rm(to_path)
            #     dir_name = os.path.dirname(to_path)
            #     if not os.path.exists(dir_name):
            #         os.makedirs(dir_name)
            #     os.link(cache_path, to_path)
            # else:
            #     to_path = cache_path
        else:
            if to_path:
                _rm(to_path)
                dir_name = os.path.dirname(to_path)
                if not os.path.exists(dir_name):
                    os.makedirs(dir_name)
                if config.get_netdata_lib("type_name") == "lustre":
                    os.link(config.convert_real_path(from_path), to_path)
                else:
                    shutil.copy(config.convert_real_path(from_path), to_path)
            else:
                to_path = from_path
        return to_path
    else:
        raise RunningError("文件%s不存在", (from_path, ), code="008")


def getsize(file_path):
    """
    获取文件大小,单位为字节数，此方法只适用于文件，文件夹不能获取正确的值
    """
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    if file_path.endswith("/"):
        raise Exception("不支持文件夹获取大小!")
    if exists(file_path):
        type_name = config.get_file_type(file_path)
        if type_name == "s3":
            key = _get_s3_key(file_path)
            size = int(key.size)
            # key.close()
        elif type_name == "http":
            http_path = config.convert_path_to_http(file_path)
            info = _get_http_info(http_path)
            if info["is_file"]:
                size = int(info["length"])
            else:
                raise Exception("文件%s可能是文件夹，无法正常获取其大小!" % http_path)
        else:
            size = int(os.path.getsize(config.convert_real_path(file_path)))
        return size
    else:
        raise RunningError("文件%s不存在", (file_path,), code="008")


def _get_http_info(file_path):
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    headers = {}
    if file_path in _http_file_info.keys():
        return _http_file_info[file_path]
    else:

        def _header_function(header_line):
            header_line = header_line.decode('iso-8859-1')
            if ':' not in header_line:
                return
            name, value = header_line.split(':', 1)
            name = name.strip()
            value = value.strip()
            name = name.lower()
            headers[name] = value

        try:
            if file_path in _http_file_info.keys():
                _http_file_info[file_path]['try'] += 1
            else:
                _http_file_info[file_path] = {"try": 1}
            c = pycurl.Curl()
            c.setopt(c.URL, file_path)
            c.setopt(c.NOBODY, 1)  # header only, no body
            c.setopt(c.FOLLOWLOCATION, 1)
            c.setopt(c.MAXREDIRS, 5)
            c.setopt(c.NOSIGNAL, 1)
            c.setopt(c.CONNECTTIMEOUT, 60)
            c.setopt(c.TIMEOUT, 300)
            c.setopt(c.HEADERFUNCTION, _header_function)
            c.perform()
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            print(e)
            sys.stdout.flush()
            if _http_file_info[file_path]["try"] >= 3:
                raise Exception("获取%s信息超过3次仍然发生错误，请与系统管理员联系。" % file_path)
            else:
                print("获取%s信息发生错误,3秒后重试..." % file_path)
                sys.stdout.flush()
                gevent.sleep(3)
                return _get_http_info(file_path)
        else:
            code = c.getinfo(c.RESPONSE_CODE)
            if code == 200:
                _http_file_info[file_path]["exits"] = True
                real_url = c.getinfo(pycurl.EFFECTIVE_URL)
                _http_file_info[file_path]["real_path"] = real_url
                if real_url.endswith("/"):
                    _http_file_info[file_path]["is_file"] = False
                else:
                    _http_file_info[file_path]["is_file"] = True
                    _http_file_info[file_path]["length"] = headers["content-length"]
                    _http_file_info[file_path]["etag"] = headers["etag"]
            elif code == 404:
                _http_file_info[file_path]["exits"] = False
            else:
                if _http_file_info[file_path]["try"] >= 3:
                    raise Exception("获取%s信息出现服务器内部错误，请与系统管理员联系。" % file_path)
                else:
                    return _get_http_info(file_path)
            return _http_file_info[file_path]


def _get_s3_key(file_path):
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    m = re.match(r"^([\w+\-]+)://([\w+\-]+)/(.*)$", file_path)
    if m:
        region = m.group(1)
        bucket = m.group(2)
        key = m.group(3)
        conn = config.get_rgw_conn(region=region, bucket=bucket)
        return Bucket(connection=conn, name=bucket).get_key(key)
    else:
        raise Exception("文件路径%s不是对象存储路径格式。" % file_path)


def _get_s3_dir_list(dir_path):
    if not PY3:
        if isinstance(dir_path, unicode):
            dir_path = dir_path.encode("utf-8")
    if not dir_path.endswith("/"):
        raise Exception("文件夹路径必须以'/'号结尾!")
    m = re.match(r"^([\w+\-]+)://([\w+\-]+)/(.*)$", dir_path)
    file_list = []
    if m:
        region = m.group(1)
        bucket = m.group(2)
        key = m.group(3)
        conn = config.get_rgw_conn(region=region, bucket=bucket)
        bucket_obj = Bucket(connection=conn, name=bucket)
        for k in bucket_obj.list(prefix=key):
            if re.search(r"/$", k.name):
                continue
            file_list.append(os.path.join("%s://%s" % (region, bucket), k.name))
    return file_list


def _get_s3_cache_path(file_path):
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    m = re.match(r"^([\w+\-]+)://([\w+\-]+)/(.*)$", file_path)
    if m:
        region = m.group(1)
        bucket = m.group(2)
        key = m.group(3)
        return os.path.join(config.get_netdata_config("s3")["s3_cache_dir"], region, bucket, key)
    else:
        raise Exception("文件路径%s不符合对象存储路径!" % file_path)


def _get_http_cache_path(http_path):
    if not PY3:
        if isinstance(http_path, unicode):
            http_path = http_path.encode("utf-8")
    m = re.match(r"^(http://|^https://)([^/]+)/(.*)$", http_path)
    if m:
        domain = m.group(2)
        key = m.group(3)
        return os.path.join(config.get_netdata_config("http")["http_cache_dir"], domain, key)
    else:
        raise Exception("文件路径%s不符合http路径!" % http_path)


def _rm(file_path):
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    if os.path.exists(file_path):
        if os.path.islink(file_path):
            os.remove(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
        else:
            os.remove(file_path)


def _get_file_part(file_path, offset, length):
    if not PY3:
        if isinstance(file_path, unicode):
            file_path = file_path.encode("utf-8")
    if file_path.endswith("/"):
        raise Exception("不支持文件夹下载!")
    if exists(file_path):
        type_name = config.get_file_type(file_path)
        if type_name == "s3":
            m = re.match(r"^([\w+\-]+)://([\w+\-]+)/(.*)$", file_path)
            if m:
                region = m.group(1)
                bucket = m.group(2)
                key = m.group(3)
                conn = config.get_rgw_conn(region=region, bucket=bucket)
                resp = conn.make_request("GET", bucket, key, headers={"Range": "bytes=%d-%d" %
                                                                               (offset, offset + length)})
                data = resp.read(len)
                return data
            else:
                raise Exception("文件路径%s不是对象存储路径格式。" % file_path)
        elif type_name == "http":
            buf = StringIO()
            http_path = config.convert_path_to_http(file_path)
            if not PY3:
                if isinstance(http_path, unicode):
                    http_path = http_path.encode("utf-8")
            c = pycurl.Curl()
            c.setopt(c.URL, http_path)
            c.setopt(c.FOLLOWLOCATION, 1)
            c.setopt(c.RANGE, '%s-%s' % (offset, offset + length))
            c.setopt(c.WRITEDATA, buf)
            c.perform()
            c.close()
            code = c.getinfo(c.RESPONSE_CODE)
            if code == 200:
                data = buf.getvalue()
                buf.close()
                return data
            else:
                return _get_file_part(file_path, offset, length)

        else:
            with open(config.convert_real_path(file_path), "r") as f:
                f.seek(offset)
                data = f.read(length)
            return data


class HttpDir(object):
    """
    获取文件夹下所有子文文件路径, 文件夹必须以"/"结尾。否则无法获取信息, 只支持http路径
    注意：文件夹子文件数过多时可能会导致进程阻塞，严重时导致系统崩溃
    """
    def __init__(self, http_path):
        if not PY3:
            if isinstance(http_path, unicode):
                http_path = http_path.encode("utf-8")
        self.http_path = http_path
        self._file_list = []
        self._dir_list = []
        self._has_get = False

    def _get_files(self, dir_path):
        d = pq(url=dir_path)
        dir_list = []
        if not d("title").text().startswith("Index of"):
            raise Exception("http目录%s中含有index.html等网页文件，不能获取文件列表信息!" % dir_path)
        for a in d("pre>a:gt(0)").items():
            path = a.attr("href")
            url = os.path.join(dir_path, path)
            if path.endswith("/"):
                self._dir_list.append(url)
                dir_list.append(url)
            else:
                self._file_list.append(url)
        for dir_url in dir_list:
            self._get_files(dir_url)
        self._has_get = True

    def get_files(self):
        """
        获取所有级别子文件的完整url路径
        :return:
        """
        if not self._has_get:
            self._get_files(self.http_path)
        return self._file_list

    def get_dirs(self):
        if not self._has_get:
            self._get_files(self.http_path)
        return self._dir_list
