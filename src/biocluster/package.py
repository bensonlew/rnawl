# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
"""
功能包基础类

"""
from multiprocessing import Process, Queue
from .core.function import load_class_by_path, get_clsname_form_path
import time
from .core.exceptions import RunningError
import traceback
import sys


class PackageManager(object):
    def __init__(self, parent):
        self._parent = parent
        self._packages = []

    @property
    def packages(self):
        return self._packages

    def add_package(self, pkg):
        if pkg not in self._packages:
            self._packages.append(pkg)

    def get(self, path, class_name=False):
        module = load_class_by_path(path, tp="Package")
        if class_name == "":
            class_name = None
        elif class_name is False:
            class_name = get_clsname_form_path(path, tp="Package")
        # for p in self._packages:
        #     if p.module == module and p.class_name == class_name:
        #         return p
        if class_name is None:
            pkg = Package(module, None, self._parent)
        else:
            if not hasattr(module, str(class_name)):
                raise Exception("模块%s中没有定义%s类" % (path, class_name))
            else:
                pkg = Package(module, class_name, self._parent)
        return pkg

    def wait(self):
        while True:
            has_run = False
            for package in self._packages:
                for p in package:
                    if p.is_alive():
                        has_run = True
                        time.sleep(1)
                    else:
                        if p.is_error():
                            if isinstance(p.error_info, dict) and "error_type" in p.error_info.keys():
                                self._parent.logger.error(p.error_info["info"])
                                self._parent.set_error(p.error_info["value"], p.error_info["variables"],
                                                       p.error_info["code"])
                            else:
                                self._parent.set_error(p.error_info)
            if not has_run:
                break


class Package(object):
    def __init__(self, module, class_name, tool):
        self._module = module
        self._class_name = class_name
        self._class_object = self._get_class()
        self._processes = None
        self._tool = tool

    @property
    def processes(self):
        return self._processes

    @property
    def module(self):
        return self._module

    @property
    def name(self):
        if self._class_name is None:
            return str(self._module)
        else:
            return str(self._class_object)

    @property
    def tool(self):
        return self._tool

    @property
    def class_name(self):
        return self._class_name

    def _get_class(self):
        try:
            if self._class_name is not None:
                cls = getattr(self._module, self._class_name)(self._tool)
                if not isinstance(cls, PackageBase):
                    raise Exception("模块%s中的类%s不是PacakgeBase的子类，不能自动加载!" %
                                    (self._module, self._class_name))
                return cls
            else:
                return False
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            print(e)
            sys.stdout.flush()
            self._tool.set_error(e)
            raise e

    def run(self, *args, **kwargs):
        try:
            self._tool.package.add_package(self)
            if self._class_name is not None:
                return self.run_method("run", *args, **kwargs)
            else:
                return self.run_function("main", *args, **kwargs)
        except Exception as e:
            exstr = traceback.format_exc()
            print(exstr)
            print(e)
            sys.stdout.flush()
            self._tool.set_error(e)

    def run_function(self, func_name, *args, **kwargs):
        if self._class_name is not None:
            raise Exception("本方法只适用于调用普通函数，请使用get(path, class_name=None)方法来获取不含类的包对象")
        else:
            if hasattr(self._module, func_name):
                self._tool.logger.info("开始运行模块%s，函数%s" % (self._module, func_name))
                return self._run_process(getattr(self._module, func_name), *args, **kwargs)
            else:
                raise Exception("模块%s中没有定义方法%s" % (self._module, func_name))

    def run_method(self, func_name, *args, **kwargs):
        if self._class_name is None:
            raise Exception("本方法只适用于定义了PackageBase子类模块的加载!")
        else:
            if hasattr(self._class_object, func_name):
                self._tool.logger.info("开始运行模块%s，类%s，方法%s" % (self._module, self.class_name, func_name))
                return self._run_process(getattr(self._class_object, func_name), *args, **kwargs)
            else:
                raise Exception("模块%s中没有定义方法%s" % (self._class_object, func_name))

    def _run_process(self, func, *args, **kwargs):
        self._has_run = True
        process = PackageProcess(self, target=func, args=args, kwargs=kwargs)
        self._processes.append(process)
        process.start()
        self._tool.logger.info("开始运行子进程PID:%s" % process.pid)
        return process

    def wait(self):
        while True:
            has_run = False
            for p in self._processes:
                if p.is_alive():
                    time.sleep(1)
                    has_run = True
                else:
                    if p.is_error():
                        if isinstance(p.error_info, dict) and "error_type" in p.error_info.keys():
                            self._tool.logger.error(p.error_info["info"])
                            self._tool.set_error(p.error_info["value"], p.error_info["variables"], p.error_info["code"])
                        else:
                            self._tool.set_error(p.error_info)
            if not has_run:
                break


class PackageProcess(Process):
    def __init__(self, package, *args, **kwargs):
        super(PackageProcess, self).__init__(*args, **kwargs)
        self._package = package
        self._return_queue = Queue()
        self._return_queue.cancel_join_thread()
        self._is_error = False
        self._error_ = None
        self._has_get = False

    @property
    def package(self):
        return self._package

    def run(self):
        if self._target:
            try:
                return_data = self._target(*self._args, **self._kwargs)
                self._return_queue.put(return_data)
            except RunningError as e:
                name = "%s(Package: %s)" % (self.package.tool.id, self.package.name)
                if e.variables:
                    info = "模块%s运行出错,ErrorCode: %s," % (name, e.code) + str(e.value) % e.variables
                else:
                    info = "模块%s运行出错,ErrorCode: %s," % (name, e.code) + str(e.value)
                data = {
                    "error_type": "running",
                    "name": name,
                    "value": e.value,
                    "code": e.code,
                    "variables": e.variables,
                    "info": info
                }
                self._return_queue.put({"_error_": data})
            except Exception as e:
                name = "%s(Package: %s)" % (self.package.tool.id, self.package.name)
                exstr = traceback.format_exc()
                print(exstr)
                self.package.tool.logger.error("%s运行出错:%s" % (name, e))
                self._return_queue.put({"_error_": "%s" % e})
                raise e

    def join(self, timeout=None):
        super(PackageProcess, self).join(timeout)
        if not self._has_get:
            data = self._return_queue.get()
            self._has_get = True
            if isinstance(data, dict) and "_error_" in data.keys():
                self._is_error = True
                self._error_ = data["_error_"]

    def is_alive(self):
        if super(PackageProcess, self).is_alive():
            return True
        else:
            if not self._has_get:
                data = self._return_queue.get()
                self._has_get = True
                if isinstance(data, dict) and "_error_" in data.keys():
                    self._is_error = True
                    self._error_ = data["_error_"]
            return False

    def is_error(self):
        return self._is_error

    @property
    def error_info(self):
        return self._error_


class PackageBase(object):

    """
    Package抽象类，所有Package应该为此类的子类，并实现run方法,并且初始化才是必须与此抽象类相同
    """

    def __init__(self, tool):
        self._tool = tool

    @property
    def tool(self):
        return self._tool

    def run(self):
        """
        抽象方法，必须在子类中实现此方法
        :return:
        """
        pass

    def set_error(self, value, variables=None, code="001"):
        e = RunningError(value, variables, code)
        e.bind_object = self._tool
        raise e
