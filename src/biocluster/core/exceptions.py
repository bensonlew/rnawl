# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import json


class Error(Exception):
    def __init__(self, value):
        Exception.__init__(self, value)
        self.value = value

    def __str__(self):
        return str(self.value)

class ExitError(Error):
    """
    用于主动跟踪程序退出
    """
    pass




class MaxLengthError(Error):
    """
    参数错误
    """
    pass

class TransferError(Error):
    """
    文件传输失败时触发
    """
    pass


class CodeError(Exception):
    """
    带有编码的错误,用于web端多语言支持
    """
    def __init__(self, value, variables=None, code="001", prefix="S"):
        """

        :param value: 错误提示，变量使用%s代替
        :param variables: tuple元组或字符串变量, 数量与value中的%s对应
        :param code:  错误编码
        """
        if isinstance(value, Exception):
            self.value = str(value)
        else:
            self.value = value

        self.variables = tuple(variables) if isinstance(variables, list) else variables
        self._code = "%s%s" % (prefix, code)
        self.bind_object = None
        Exception.__init__(self, self.__str__())

    def __str__(self):
        s = ""
        if self.bind_object:
            s = "模块%s运行出错," % self.bind_object.id
        if self.variables:
            s += "ErrorCode: %s, " % self.code + self.value % self.variables
        else:
            s += "ErrorCode: %s, " % self.code + self.value
        return s

    @property
    def code(self):
        return self._code

    def json(self):
        pass

    @property
    def info(self):
        if self.variables:
            return self.value % self.variables
        else:
            return self.value

class OptionError(CodeError):
    """
    参数错误
    """
    pass   

class EventStopError(Error):
    """
    事件已停止异常，发生此异常时是由于事件已经停止监听或尚未启动监听时触发此事件
    """
    pass


class UnknownEventError(Error):
    """
    事件未定义异常，发生此异常时是由于事件尚未定义时触发此事件
    """

    pass


class OptionError(CodeError):
    """
    参数错误
    """
    def __init__(self, value, variables=None, code="001"):
        self.option = None
        super(OptionError, self).__init__(value, variables, code, "O")

    def __str__(self):
        s = ""
        if self.bind_object:
            s = "模块%s " % self.bind_object.id
        if self.option:
            s += "参数%s检测出错," % self.option
        if self.variables:
            s += "ErrorCode: %s, " % self.code + self.value % self.variables
        else:
            s += "ErrorCode: %s, " % self.code + self.value
        return s

    def json(self):
        data = {
            "error_type": "option",
            "name": self.bind_object.id if self.bind_object else None,
            "option": self.option,
            "code": self.code,
            "variables": self.variables,
            "info": self.info
        }
        return data


class FileError(CodeError):
    """
    文件错误
    """
    def __init__(self, value, variables=None, code="001"):
        self.option = None
        self.file_name = None
        super(FileError, self).__init__(value, variables, code, "F")

    def __str__(self):
        s = ""
        if self.bind_object:
            s = "模块%s " % self.bind_object.id
        if self.option:
            s += "参数%s," % self.option
        if self.file_name:
            s += "文件%s检测出错," % self.file_name
        if self.variables:
            s += "ErrorCode: %s, " % self.code + self.value % self.variables
        else:
            s += "ErrorCode: %s, " % self.code + self.value
        return s

    def json(self):
        data = {
            "error_type": "file",
            "name": self.bind_object.id if self.bind_object else None,
            "file_name": self.file_name,
            "option": self.option,
            "code": self.code,
            "variables": self.variables,
            "info": self.info
        }
        return data


class RunningError(CodeError):
    """
    运行错误
    """
    def __init__(self, value, variables=None, code="001"):
        super(RunningError, self).__init__(value, variables, code, "R")

    def json(self):
        data = {
            "error_type": "running",
            "name": self.bind_object.id if self.bind_object else None,
            "code": self.code,
            "variables": self.variables,
            "info": self.info
        }
        return data


class MaxLengthError(Error):
    """
    参数错误
    """
    pass
