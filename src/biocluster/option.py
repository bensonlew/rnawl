# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""参数类"""

from .iofile import FileBase
from .core.function import load_class_by_path, get_classpath_by_object
import os
from .core.exceptions import OptionError
import re
import traceback
from .core.exceptions import FileError, CodeError
from .api.file.remote import RemoteFileManager
import sys


class Option(object):
    """
    参数选项的属性库设定
    """
    def __init__(self, opt, bind_obj=None):
        """
        初始化一个参数
        参数使用
        """
        self._options = opt
        self._format = None
        self._check = None
        self._format_list = []
        self._check_list = []
        self.bind_obj = bind_obj
        if not isinstance(opt, dict):
            raise Exception("opt必须为一个字典")
        for attr in ('name', 'type'):
            if attr not in opt.keys():
                raise Exception("必须设置参数属性 {}".format(attr))
        self._name = opt['name']
        self._type = opt['type']
        self._required = opt['required'] if "required" in opt.keys() else False
        self._max = opt['max'] if "max" in opt.keys() else None
        self._min = opt['min'] if "min" in opt.keys() else None
        self._choose = opt['choose'] if "choose" in opt.keys() else None
        self._max_choose = opt['max_choose'] if "max_choose" in opt.keys() else None
        self._is_set = False
        if opt['type'] in {'outfile', 'infile'}:
            if 'format' not in opt.keys():
                raise Exception("必须设置参数属性 format")
            else:
                formats = re.split('\s*,\s*', opt['format'])
                if len(formats) > 1:
                    if 'check' in opt.keys():
                        checks = re.split('\s*,\s*', opt['check'])
                        if len(checks) == 1:
                            for index in range(len(formats)):
                                self._format_list.append(formats[index].strip())
                                self._check_list.append(checks[0].strip())
                        else:
                            for index in range(len(formats)):
                                if index + 1 > len(checks):
                                    self._format_list.append(formats[index].strip())
                                    self._check_list.append(False)
                                else:
                                    self._format_list[index] = formats[index].strip()
                                    check_str = checks[index].strip()
                                    self._check_list.append(check_str if check_str != "" else False)
                    else:
                        # print formats
                        for index in range(len(formats)):
                            # print index
                            self._format_list.append(formats[index].strip())
                            self._check_list.append(False)
                    self._value = False
                else:
                    self._format = opt['format'].strip()
                    self._check = opt['check'].strip() if 'check' in opt.keys() else False
                    self._value = load_class_by_path(self._format, "File")()
        else:
            if opt['type'] == "select":
                if not self._choose:
                    raise Exception("select参数必须设置choose选择范围")
            if self._choose:
                if not (isinstance(self._choose, list) or isinstance(self._choose, tuple)):
                    raise Exception("choose选择范围必须为数组或元组")

            if 'default' in opt.keys():
                self._value = opt['default']
            else:
                self.bind_obj.logger.warning("参数%s没有默认值，请确认确实不需要默认值？" % self._name)
                self._value = None

        if opt['type'] not in {"int", "float", "string", "select", "bool", "infile", "outfile"}:
            raise Exception("参数属性不在规范范围内type：{}".format(self._type))

    @property
    def name(self):
        """
        获取参数名

        :return: string
        """
        return self._name

    @property
    def type(self):
        """
        获取参数属性

        :return: String
        """
        return self._type

    @property
    def format_list(self):
        """
        多格式支持时返回格式列表
        :return:
        """
        return self._format_list

    @property
    def value(self):
        """
        获取参数值

        也可直接赋值,当参数为文件类参数时自动判断，当value为实际存在的文件路径时，调用set_path设置File对象路径。
        当value为File文件对象时，检测是否与自身文件类文件格式类型相同，如果相同则替换,不相同时报错

        :return:
        """
        return self._value

    @value.setter
    def value(self, value):
        """


        :param value:
        :return:
        """
        if self._type in {'outfile', 'infile'}:
            if isinstance(value, unicode) or isinstance(value, str):
                path_list = value.split("||")
                if len(path_list) > 1:
                    file_path = path_list[1]
                    file_format = path_list[0]
                else:
                    file_path = path_list[0]
                    file_format = None
                if self._type == "infile":  # 远程文件复制
                    remote_file = RemoteFileManager(file_path)
                    if remote_file.type != "local":
                        self.bind_obj.logger.info("发现参数%s为远程文件%s,开始复制..." % (self.name, value))
                        remote_file.download(os.path.join(self.bind_obj.work_dir, "remote_input", self.name))
                        file_path = remote_file.local_path
                if os.path.exists(file_path):
                    # if self.type == "infile":  # 检查输出文件是否满足要求
                        # class_obj = load_class_by_path(self._options[name].format, "File")
                    if file_format is not None:
                        if len(self._format_list) > 1:
                            if file_format not in self._format_list:
                                # e = "输入参数%s的文件格式%s必须在范围%s内!" % (self.name, file_format, self._format_list)
                                self._file_check_error("输入的文件格式%s不在范围%s内!",
                                                       (file_format, self._format_list), "010")
                            else:
                                self._format = file_format
                                self._check = self._check_list[self._format_list.index(file_format)]
                        else:
                            if file_format != self._format:
                                # e = "输入参数%s的文件格式必须为%s,不能为%s!" % (self.name, self._format, file_format)
                                self._file_check_error("输入的文件格式必须为%s,不能为%s!",
                                                       (self._format, file_format), "011")
                        success, file_obj_or_error = self._check_file(self._format, self._check, file_path)
                        if success and file_obj_or_error:
                            self._value = file_obj_or_error
                        else:
                            # e = "输入参数%s的文件格式不正确:%s!" % (self.name, file_obj_or_error)
                            self._file_check_error(file_obj_or_error)
                    else:
                        if len(self._format_list) > 1:
                            has_pass = False
                            error_info = ""
                            for index in range(len(self._format_list)):
                                format_path = self._format_list[index]
                                check = self._check_list[index]
                                success, file_obj_or_error = self._check_file(format_path, check, file_path, loop=True)
                                if success:
                                    has_pass = True
                                    self._value = file_obj_or_error
                                    self._format = format_path
                                    self._check = check
                                    break
                                else:
                                    error_info += "格式%s:%s " % (format_path, file_obj_or_error)
                            if not has_pass:
                                self.bind_obj.logger.error(error_info)
                                self._file_check_error("多格式文件检测未通过!", None, "012")
                        else:
                            success, file_obj_or_error = self._check_file(self._format, self._check, file_path)
                            if success and file_obj_or_error:
                                self._value = file_obj_or_error
                            else:
                                self._file_check_error(file_obj_or_error)
                else:
                    # raise OptionError("参数%s的输入文件%s不存在，请确认路径正确！", (self.name, file_path), "002")
                    # e = OptionError("文件%s不存在！", (os.path.basename(file_path)), 002)
                    # e.bind_object = self.bind_obj
                    # e.option = self.name
                    # self.bind_obj.set_error(e)
                    self._file_check_error("文件%s不存在！", (os.path.basename(file_path)), "002")
            else:
                if len(self._format_list) > 1:
                    has_pass = False
                    self._check_type(value)
                    for index in range(len(self._format_list)):
                        class_obj = load_class_by_path(self._format_list[index], "File")
                        if isinstance(value, class_obj):
                            self._value = value
                            self._format = self._format_list[index]
                            self._check = self._check_list[index]
                            has_pass = True
                            break
                    if not has_pass:
                        self._file_check_error("文件对象类型设置不正确!", None, "013")
                else:
                    self._check_type(value)
                    class_obj = load_class_by_path(self._format, "File")
                    if isinstance(value, class_obj):
                        self._value = value
                    else:
                        self._file_check_error("文件对象类型设置不正确!", None, "013")
            if isinstance(self._value, FileBase):
                self._value.option = self
                self._value.parent = None
        else:
            self._value = self._check_type(value)
        self._is_set = True

    @property
    def format(self):
        """
        获取文件类型参数
        """
        return self._format

    @property
    def check(self):
        """
        获取文件类型参数
        """
        return self._check

    @property
    def is_set(self):
        """
        返回参数对象是否被设置了值
        :return:
        """
        return self._is_set

    def _check_type(self, value):
        """
        检测值是否符合要求
        :param value:
        :return:
        """
        try:
            if self._type == "int":
                try:
                    value = int(value)
                except ValueError:
                    raise OptionError("参数值类型不符合%s:%s", (self._type, value), "003")
            if self._type == "float":
                try:
                    value = float(value)
                except ValueError:
                    raise OptionError("参数值类型不符合%s:%s", (self._type, value), "003")
            if self._type == "bool":
                try:
                    if isinstance(value, str) or isinstance(value, unicode):
                        if re.match(r"^[\-\+]?\d+$", value):
                            if int(value) > 0:
                                value = True
                            else:
                                value = False
                        elif re.match(r"true", value, re.I) or re.match(r"yes", value, re.I):
                            value = True
                        elif re.match(r"false", value, re.I) or re.match(r"no", value, re.I):
                            value = False
                    if isinstance(value, int) or isinstance(value, float):
                        if int(value) > 0:
                            value = True
                        else:
                            value = False
                    if value is None:
                        value = False
                except ValueError:
                    raise OptionError("参数值类型不符合%s:%s", (self._type, value), "003")
                if not isinstance(value, bool):
                    raise OptionError("参数值类型不符合%s:%s", (self._type, value), "003")
            if self._type == "string":
                if not (isinstance(value, unicode) or isinstance(value, str)):
                    raise OptionError("参数值类型不符合%s:%s", (self._type, value), "003")
            if self._type in {"infile", "outfile"}:
                if not isinstance(value, FileBase):
                    raise OptionError("参数值类型不符合%s:%s", (self._type, value), "003")
            if self._type in ["int", "float"]:
                if self._max and value > self._max:
                        raise OptionError("参数值不能大于%s", self._max, "004")
                if self._min and value < self._min:
                    raise OptionError("参数值不能小于%s", self._min, "005")
            if self._type in ["int", "float", "string"]:
                if self._choose:
                    if value not in self._choose:
                        raise OptionError("参数值必须在选择范围内!", None, "006")
            if self._type == "select":
                if isinstance(value, list) or isinstance(value, tuple):
                    if self._max_choose and len(value) > self._max_choose:
                        raise OptionError("选择参数个数不能超过最大限制数%s个!", self._max_choose, "007")
                    for v in value:
                        if v not in self._choose:
                            raise OptionError("参数值%s不在选择范围内!", value, "008")
                else:
                    raise OptionError("参数值必须为数组!", None, "009")
        except OptionError as e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            e.bind_object = self.bind_obj
            e.option = self.name
            self.bind_obj.set_error(e)
        return value

    def _check_file(self, format_path, check, path, loop=False):
        """

        :param format_path:
        :param check:
        :return:
        """
        class_name = self.bind_obj.__class__.__name__
        self_class_path = get_classpath_by_object(self.bind_obj)
        paths = self_class_path.split(".")[2:]
        function_name = "_".join(paths)
        if re.search(r"Agent$", class_name) or re.search(r"Tool$", class_name):
            function_name += "_agent_check"
        elif re.search(r"Module$", class_name):
            function_name += "_module_check"
        elif re.search(r"Workflow$", class_name):
            function_name += "_workflow_check"
        else:
            raise Exception("类名称不正确!")
        try:
            file_obj = load_class_by_path(format_path, "File")()
            file_obj.set_path(path)
            if check:
                if hasattr(file_obj, check):
                    getattr(file_obj, check)()
                else:
                    raise Exception("文件类%s中未定义指定的检测函数%s!" %
                                    (format_path, check))
            else:
                if hasattr(file_obj, function_name):
                    getattr(file_obj, function_name)()
                else:
                    getattr(file_obj, "check")()
        except FileError, e:
            e.bind_object = self.bind_obj
            e.option = self.name
            e.file_name = os.path.basename(path)
            exstr = traceback.format_exc()
            if loop:
                self.bind_obj.logger.debug("检测未通过(以下为调试信息，可忽略):\n%s" % exstr)
            else:
                print exstr
            sys.stdout.flush()
            return False, e
        except Exception, e:
            exstr = traceback.format_exc()
            print exstr
            sys.stdout.flush()
            # self._file_check_error(str(e))
            return False, e
        else:
            return True, file_obj

    def _file_check_error(self, value, variables=None, code="001"):
        """
        文件检测错误后的处理

        :param error:
        :return:
        """
        # class_name = self.bind_obj.__class__.__name__
        # if re.search(r"Tool$", class_name):
        #     self.bind_obj.set_error(error)
        # else:
        #     self.bind_obj.get_workflow().exit(exitcode=1, data=error)
        if isinstance(value, CodeError):
            e = value
        else:
            e = OptionError(str(value), variables, code)
            e.bind_object = self.bind_obj
            e.option = self.name
        self.bind_obj.set_error(e)

    def run_file_check(self, name, *args, **kwargs):
        if self._type in {"infile", "outfile"} and self.is_set:
            if hasattr(self.value, name):
                try:
                    getattr(self.value, name)(*args, **kwargs)
                except FileError, e:
                    exstr = traceback.format_exc()
                    print exstr
                    sys.stdout.flush()
                    e.bind_object = self.bind_obj
                    e.option = self.name
                    e.file_name = os.path.basename(self.value.path)
                    self.bind_obj.set_error(e)
            else:
                raise Exception("文件对象%s不存在检测方法%s" % (self.value.format, name))

