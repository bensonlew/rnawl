# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""文件基本模块"""

import os
import hashlib
from .core.exceptions import FileError
import re
import pickle
from .core.function import get_dir_size, load_class_by_path
import gzip
import gevent
import glob
import types
import sys

PY3 = sys.version_info[0] == 3


class FileBase(object):
    """
    文件对象基类模块
    """

    def __init__(self):
        """
        """
        self._properties = {}
        self._is_set = False
        self.option = None
        self.parent = None

    @property
    def prop(self):
        """
        返回文件对象所有属性，dict字典 key是属性名
        """
        return self._properties

    @property
    def is_set(self):
        """
        返回文件对象是否被设置过路径
        :return:
        """
        return self._is_set

    @property
    def path(self):
        """
        返回文件路径

        :return:
        """
        if self.is_set:
            return self._properties["path"]
        else:
            return None

    @property
    def size(self):
        """
        获取文件大小
        :return:
        """
        if "size" in self._properties.keys():
            return self._properties["size"]
        else:
            return 0

    @property
    def format(self):
        """
        返回文件类的Path路径

        :return:
        """
        class_name = str(type(self))
        m = re.match(r"<class\s\'(.*)\'>", class_name)
        class_name = m.group(1)
        paths = class_name.split(".")
        paths.pop()
        paths.pop(0)
        paths.pop(0)
        return ".".join(paths)

    def set_property(self, name, value):
        """
        添加文件对象属性

        :param name: string 属性名称
        :param value: 属性值
        """
        try:
            pickle.dumps(value)
        except Exception:
            raise Exception("属性%s的值%s不能被序列化，请重新设置!" % (name, value))
        self._properties[name] = value
        return self

    def get_info(self):
        """
        获取当前文件对象的所有属性,需在扩展子类中重写此方法
        """
        if self.is_set:
            self.set_property('dirname', os.path.dirname(self.prop['path']))
            self.set_property('basename', os.path.basename(self.prop['path']))
        else:
            raise Exception("请先调用set_path方法设置文件路径!!")

    def set_path(self, path):
        """
        设置当前文件对象的路径，并获取所有属性

        :param path: 路径
        """
        if not os.path.exists(path):
            raise Exception("%s不存在，请先设置正确的文件或文件夹路径路径!" % path)
        self.set_property("path", path)
        # self.get_info()
        self._is_set = True
        if self.option:
            self.option._is_set = True

    def check(self):
        """
        检测当前文件对象是否满足运行需求,需在扩展子类中重写
        """
        if not self.is_set:
            raise Exception("请先调用set_path方法设置文件路径!")
        return True

    def link(self, link_path=None):
        """
        为当前文件对象创建软连接

        :param link_path:  生成链接的目标文件夹或希望生成的链接文件路径  所在文件夹必须已经存在
        :return: string 生成的链接路径
        """
        if "path" in self.prop.keys() and os.path.exists(self.prop['path']):
            source_path = os.path.abspath(self.prop['path'])
            file_name = os.path.basename(source_path)
            if link_path:
                if os.path.isdir(link_path):   # 指定目标为文件夹
                    target_path = os.path.join(link_path, file_name)
                    os.symlink(source_path, target_path)
                    return target_path
                else:    # 指定目标为文件
                    os.symlink(source_path, link_path)
                    return link_path
            else:   # 未指定目标 默认链接到当前目录下
                os.symlink(source_path, file_name)
                return os.path.join(os.getcwd(), file_name)

    @property
    def api(self):
        """
        获取文件所属参数Workflow/Module/Agent/Tool的ApiManager对象
        :return:
        """
        parent = self
        while parent.parent is not None:
            parent = parent.parent

        if not parent.option:
            raise Exception("文件对象必须属于某个参数！")
        return parent.option.bind_obj.api


class File(FileBase):
    """
    单个文件对象
    """
    def __init__(self):
        super(File, self).__init__()
        self._is_gzip_mark = None

    def check(self):
        """
        检测文件是否正常,并可用于后续分析

        :return:
        """
        super(File, self).check()
        if not ('path' in self.prop.keys() and os.path.isfile(self.prop['path'])):
            raise FileError("文件路径不在正确或文件不存在", None, "002")
        return True

    def get_md5(self):
        """
        获取文件MD5值
        """
        md5 = hashlib.md5()
        with open(self.prop['path'], 'rb') as f:
            while 1:
                data = f.read(8096)
                if not data:
                    break
                md5.update(data)
        return md5.hexdigest()

    def get_size(self):
        """
        获取文件大小
        """
        if os.path.islink(self.path):
            path = os.path.abspath(os.path.realpath(self.path))
        else:
            path = self.path
        size = os.path.getsize(path)
        return size

    def get_reader(self):
        """
        获取文件读取句柄，只读,注意读完成后需要将句柄关闭。
        :return: 返回文件打开句柄，支持gzip压缩文件作为文本文件读取
        """
        if self.is_gzip:
            f = gzip.open(self.path, 'rb')
        else:
            f = open(self.path, "r")
        return f

    @property
    def is_gzip(self):
        """
        判断是否为压缩违建
        :return:
        """
        if self._is_gzip_mark is None:
            return self._is_gzip()
        else:
            return self._is_gzip_mark

    def _is_gzip(self):
        if self.is_set:
            try:
                with gzip.open(self.path, 'rb') as f:
                    f.readline()
            except (IOError, Exception):
                self._is_gzip_mark = False
                return False
            else:
                self._is_gzip_mark = True
                return True
        return False

    def ungzip(self):
        """
        将gzip文件解压成新的文件存储,对于大文件建议在tool中运行此方法
        :return:
        """
        if self.is_gzip:
            to_path = self._create_unzip_path()
            if os.path.exists(to_path) and os.path.isfile(to_path):
                os.remove(to_path)
            g_zip = GZipTool(1048576)
            g_zip.decompress(self.path, to_path)
            self.set_property("ungzip_path", to_path)
            return to_path

    def gzip(self):
        """
        将文件解压成gz文件存储
        :return:
        """
        if not self.is_gzip:
            to_path = "%s.gz" % self.path
            if os.path.exists(to_path) and os.path.isfile(to_path):
                os.remove(to_path)
            g_zip = GZipTool(1048576)
            g_zip.compress(self.path, to_path)
            self.set_property("gzip_path", to_path)
            return to_path

    @property
    def ungzip_path(self):
        """
        解压缩后的文件路径，只有执行ungzip方法后才有值
        :return:
        """
        if "size" in self._properties.keys():
            return self._properties["ungzip_path"]
        else:
            return None

    @property
    def gzip_path(self):
        """
        压缩后的文件路径，只有执行gzip方法后才有值
        :return:
        """
        if "size" in self._properties.keys():
            return self._properties["gzip_path"]
        else:
            return None

    def _create_unzip_path(self):
        m = re.match(r"^(.*)\.gz$", self.path)
        if m:
            return m.group(1)
        else:
            return "%s.ungzip" % self.path

    def get_first_lines(self, line=10):
        """
        获取文件的最开始几行，一般用于大文件格式是否正确的检测，注意不要将line值设置过大，否则可能超出内存和使运行时间过长
        :param line: 指定需要读取几行
        :return: list数组。每一行为一个数组元素，去除了最后的换行符
        """
        count = 0
        data = []
        with self.get_reader() as f:
            for l in f:  # 一行行的把数据从硬盘加载到内存里读出来
                if count < line:  # 读取前五行
                    data.append(l.strip())
                    count += 1
                else:
                    break
        return data

    def get_last_lines(self, line=10):
        """
        获取文件的最后几行，一般用于大文件格式是否正确的检测，注意不要将line值设置过大，否则可能超出内存和使运行时间过长
        :param line: 指定需要读取几行
        :return: list数组。每一行为一个数组元素，去除了最后的换行符
        """
        line_size = -100
        with self.get_reader() as f:
            while True:
                off_set = line_size * line
                try:
                    f.seek(off_set, 2)
                except IOError:
                    f.seek(0, 0)
                    lines = f.readlines()
                    data = [l.strip() for l in lines]
                    break
                else:
                    lines = f.readlines()
                    if len(lines) >= line + 1:
                        data = [l.strip() for l in lines[-line:]]
                        break
                    else:
                        line_size *= 2
        return data


class Directory(FileBase):
    """
    文件夹对象
    """
    def __init__(self):
        super(Directory, self).__init__()
        self._file_list = []

    def check(self):
        """

        :return:
        """
        super(Directory, self).check()
        if not('path' in self.prop.keys() and os.path.isdir(self.prop['path'])):
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!", None, "002")
        for f in self._file_list:
            f.check()
            gevent.sleep(0)
        return True

    def add_files(self, file_type, path="*"):
        """
        添加文件夹下属文件

        :param file_type: 添加的文件类型自动加载路径
        :param path:  需要添加文件的路径，可以使用通配符，但是只能添加当前文件夹下的子文件或文件夹，
        不包含子文件夹迭代,默认"*"
        :return:
        """
        if not self.is_set:
            raise Exception("必须先设置路径!")
        class_module = load_class_by_path(file_type, tp="File")
        module = class_module()
        is_file = False
        is_dir = False
        if isinstance(module, File):
            is_file = True
        elif isinstance(module, Directory):
            is_dir = True
        if PY3:
            check = isinstance(path, str)
        else:
            check = isinstance(path, types.StringTypes)
        if check:
            if re.match(r"^/|^\.\.", path):
                raise Exception("不能使用绝对路径或切换到其他目录!")
            for f in glob.glob(os.path.join(self.path, path)):
                if (os.path.isfile(f) and is_file) or (os.path.isdir(f) and is_dir):
                    module = class_module()
                    module.option = self.option
                    module.parent = self
                    module.set_path(f)
                    self._file_list.append(module)
        else:
            module = class_module()
            module.option = self.option
            module.parent = self
            module.set_path(path)
            self._file_list.append(module)

    @property
    def file_number(self):
        """
        获取下属文件对象数量,只包含直接下属，不包含间接下属
        :return:
        """
        return len(self._file_list)

    @property
    def files(self):
        """
        获取所有下属文件对象
        :return:
        """
        return self._file_list

    def get_info(self):
        super(Directory, self).get_info()
        self.set_property("file_number", self.file_number)

    def get_size(self):
        """
        获取文件大小
        """
        super(Directory, self).check()
        count, size = get_dir_size(self.path)
        return size


class GZipTool(object):
    """
    gzip文件压缩解压缩工具
    """

    def __init__(self, buf_size):
        self.buf_size = buf_size
        self.fin = None
        self.fout = None

    def compress(self, src, dst):
        self.fin = open(src, 'rb')
        self.fout = gzip.open(dst, 'wb')
        self.__in2out()

    def decompress(self, gz_file, dst):
        self.fin = gzip.open(gz_file, 'rb')
        self.fout = open(dst, 'wb')
        self.__in2out()

    def __in2out(self,):
        while True:
            buf = self.fin.read(self.buf_size)
            if len(buf) < 1:
                break
            self.fout.write(buf)
        self.fin.close()
        self.fout.close()
