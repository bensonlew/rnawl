# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import gevent
import types
from gevent.event import AsyncResult
import inspect
# from gevent.lock import BoundedSemaphore
# import socket
from .exceptions import EventStopError, UnknownEventError, CodeError
import traceback
from datetime import datetime
import sys


class EventHandler(object):
    """
    事件处理器，核心事件处理类
    """

    def __init__(self, name):
        """
        Constructor
        """
        self._name = name
        self._event = AsyncResult()
        self._func = []
        self._greenlets = []
        self._start = False
        self.parent = None
        # self.sem = BoundedSemaphore()

    @property
    def name(self):
        """
        只读属性，获取事件名称
        :return:str 返回事件名称
        """
        return self._name

    @property
    def is_start(self):
        """
        只读属性,返回是否事件开始监听
        """
        return self._start

    def __bind(self, func, bindobject, data):
        """
        将函数包装成事件触发函数
        :param func: 事件处理函数
        :param data: 需要传递给事件处理函数的数据
        :return: waiter函数
        """

        def waiter():
            para = self._event.get()
            start = datetime.now()
            if hasattr(bindobject, "logger"):
                bindobject.logger.debug("Event %s , Function %s 开始执行..." % (self.name, func.func_name))
            argspec = inspect.getargspec(func)
            args = argspec.args
            try:
                if self.is_start:
                    event_data = {"bind_object": bindobject, "event": self, "data": data}
                    if inspect.ismethod(func):
                        args.pop(0)
                    length = len(args)
                    if length == 0:
                        func()
                    if length == 1:
                        if args[0] == 'event':
                            func(event_data)
                        else:
                            func(para)
                    if length == 2:
                        if args[0] == 'event':
                            func(event_data, para)
                        else:
                            func(para, event_data)
                    if length > 2:
                        if inspect.ismethod(func):
                            raise Exception("指定的绑定函数%s为bound method,参数超过限制个数!" % func)
                        else:
                            raise Exception("指定的绑定函数%s为unbound function,参数超过限制个数!" % func)
            except Exception, e:
                exstr = traceback.format_exc()
                print exstr
                sys.stdout.flush()
                if isinstance(e, CodeError):
                    e.bind_object = bindobject
                if hasattr(bindobject, "set_error"):
                    bindobject.set_error(e)
            end = datetime.now()
            if hasattr(bindobject, "logger"):
                bindobject.logger.debug("Event %s, Function %s 执行完毕，耗时%s s。" % (self.name, func.func_name,
                                        (end - start).seconds))
        return waiter

    def bind(self, func, bindobject, data=None):
        """
        添加事件监听函数

        :param func: 为自定义unbound function 或bound method(对象方法或类方法),最多允许2个参数（bound method第一个参数除外,即self）,
                      bound method定义时第一个参数应该为self或cls,如:

        .. code-block:: python
            :linenos:

            class A(object):
                def b(self,c):
                    pass
            x=A()
            obj.bind(A.b)

        function可定义为:``def func(a,b)`` 或 ``def func(a)`` 或 ``def func()``

        1. 第一个参数运行时获取值为 **EventHandler.fire** 函数传递的参数

        #. 第二个参数运行时获取值为一个 **event_data** 字典

           * *event_data['bind_object']* 值是该方法绑定的EventObject对象

           * *event_data['event']* 值是该方法绑定的EventHandler对象

           * *event_data['data']* 是bind方法传递的data参数

        #. 当定义时指定参数名为 *event* 时(如 ``def func(event)`` )，此时*event*的值为 **event_data** 字典,不区分参数位置

        :param bindobject: 事件绑定的对象
        :param data: 可选参数 将通过 **event_data['data']** 形式传递给运行函数 默认值: **None**
        :return: self
        """
        if self.is_start:
            raise Exception("模块%s, %s事件已经启动监听，绑定事件处理函数应该在启动事件前进行!" %
                            (self.parent, self.name))
        else:
            waiter = self.__bind(func, bindobject, data)
            self._func.append(waiter)
        return self

    def fire(self, para=None):
        """
        触发此事件,执行后此事件绑定的方法将会执行

        :param para: 可选参数，触发事件时传递给事件绑定函数的参数
        """
        # with self.sem:
        if self.is_start:
            self._event.set(para)
        else:
            # print self.name, self
            raise EventStopError("模块%s, 事件%s尚未开始监听或已经运行，不能被触发!有子模块未完成而此模块先结束"
                                 "时可能会导致此问题，请搜索Warning日志" % (self.parent, self.name))

    def start(self):
        """
        开始启动事件监听,启动事件监听后将不能添加新的绑定事件
        事件触发前，必须先启动事件监听

        :return: self
        """
        for item in self._func:
            self._greenlets.append(gevent.spawn(item))
        self._start = True
        return self

    def stop(self):
        """
        停止事件监听,停止事件监听后事件将不能被再次触发::

            注意：不能在事件触发的处理函数中来停止当前事件监听

        :return: self
        """
        current = gevent.getcurrent()
        for gl in self._greenlets:
            if gl is current:
                raise Exception("模块%s, 不能在当前事件触发的函数中停止当前事件监听!" % self.parent)
        if self._event.ready():
            if len(self._greenlets) > 0:
                gevent.joinall(self._greenlets)
        else:
            if len(self._greenlets) > 0:
                gevent.killall(self._greenlets, block=False)
        self._start = False
        return self

    def restart(self):
        """
        重启事件监听,将已经stop的事件重新启动,并重置所有事件处理函数等待再次被触发

        :return: self
        """
        if self.is_start:
            # raise Exception("%s事件已经启动监听，需的调用EventHandler.stop方法停止后才能重新启动!" % self.name)
            self.stop()
        self._event = AsyncResult()
        self._greenlets = []
        self.start()
        return self


class LoopEventHandler(object):
    """
    可循环触发的事件处理器, 继承至 **EventHandler**
    """

    def __init__(self, name):
        """
        Constructor
        """
        self._name = name
        # self._raw_handler = EventHandler(name)
        self._loop_handlers = []
        self._start = False
        # self.sem = BoundedSemaphore()
        self._fire_count = 0
        self._bind_list = []
        self.parent = None

    @property
    def name(self):
        """
        只读属性，获取事件名称
        :return:str 返回事件名称
        """
        return self._name

    @property
    def is_start(self):
        """
        只读属性,返回是否事件开始监听
        """
        return self._start

    def bind(self, func, bindobject, data=None):
        """
        添加事件监听函数

        :param func: 为自定义unbound function 或bound method(对象方法或类方法),最多允许2个参数（bound method第一个参数除外,即self）,
                      bound method定义时第一个参数应该为self或cls,如:

        .. code-block:: python
            :linenos:

            class A(object):
                def b(self,c):
                    pass
            x=A()
            obj.bind(A.b)

        function可定义为:``def func(a,b)`` 或 ``def func(a)`` 或 ``def func()``

        1. 第一个参数运行时获取值为 **EventHandler.fire** 函数传递的参数

        #. 第二个参数运行时获取值为一个 **event_data** 字典

           * *event_data['bind_object']* 值是该方法绑定的EventObject对象

           * *event_data['event']* 值是该方法绑定的EventHandler对象

           * *event_data['data']* 是bind方法传递的data参数

        #. 当定义时指定参数名为 *event* 时(如 ``def func(event)`` )，此时*event*的值为 **event_data** 字典,不区分参数位置

        :param func:
        :param bindobject:
        :param data:
        :return:
        """
        if self.is_start:
            raise Exception("模块%s, %s事件已经启动监听，绑定事件处理函数应该在启动事件前进行!" %
                            (self.parent, self.name))
        else:
            self._bind_list.append([func, bindobject, data])
        return self

    def fire(self, para=None):
        """
        触发事件::

            注意：不能在事件绑定的处理函数中尝试再次触发本事件，否则将导致异常

        :param para: 可选参数，触发事件时传递给事件绑定函数的参数 默认值:None
        """
        if self.is_start:
            # handler = copy.deepcopy(self._raw_handler)
            self._fire_count += 1
            name = "%s[%s]" % (self.name, self._fire_count)
            handler = EventHandler(name)
            for bind_para in self._bind_list:
                handler.bind(*bind_para)
            self._loop_handlers.append(handler)
            handler.start()
            handler.fire(para)
        else:
            raise EventStopError("模块%s, 事件%s尚未开始监听或已经运行，不能被触发!" % (self.parent, self.name))

        # print "%s start fire ....." % self.name
        # with self.sem:
        #     if not self.is_start:
        #         # print "%s 没有启动 ....." % self.name
        #         raise EventStopError("%s事件没有启动!" % self.name)
        #     if self._event.ready():
        #         self.stop()
        #         self.restart()
        #     self._event.set(para)
        # # self._current_handler.fire(para)
        # # print "%s Fire ................" % self.name

    def start(self):
        """
        开始启动事件监听,启动事件监听后将不能添加新的绑定事件
        事件触发前，必须先启动事件监听

        :return: self
        """
        self._start = True
        return self

    def stop(self):
        """
        停止事件监听,停止事件监听后事件将不能被再次触发::

            注意：不能在事件触发的处理函数中来停止当前事件监听

        :return: self
        """
        for handler in self._loop_handlers:
            handler.stop()
        self._start = False
        return self

    def restart(self):
        """
        重启事件监听,将已经stop的事件重新启动,并重置所有事件处理函数等待再次被触发

        :return: self
        """
        if not self.is_start:
            self.start()
        return self


class EventObject(object):
    """
    核心事件对象，用于扩展出支持事件绑定与触发的对象。
    每个 **EventObject** 可以添加多个事件，每个事件可以绑定多个对应的函数，当事件被触发时，绑定的函数即开始执行

    使用举例:

    .. code-block:: python
            :linenos:

            from biocluster.core.event import EventObject


            class ExtendEventObject(EventObject):

                def __init__(self):
                    super(ExtendEventObject, self).__init__()
                    self.test1 = "I'm testfunc1"
                    self.test2 = "I'm testfunc2"
                    self.add_event("test1")         # 添加事件test1
                    self.on("test1", self.testfunc1)  # 绑定事件方法  bound function
                    self.add_event("test2")   # 添加事件test2
                    self.on("test2", self.testfunc2)  # 绑定事件方法  bound function

                def testfunc1(self):
                    print self.test1

                def testfunc2(self, x):
                    print "%s and %s" % (self.test2, x)


            def testfunc3(eventobj):
                print "wait function will exit"
                eventobj.set_end()   # 设置EventObject为完成状态，wait()函数即停止阻塞主线程

            testobj = ExtendEventObject()  # 实例化类
            testobj.add_event("test3")    # 添加事件test3
            testobj.on("test3", testfunc3)  # 绑定事件方法 unbound function
            testobj.start_listener()  # 启动事件监听
            testobj.fire("test1")   # 触发事件test1
            testobj.fire("test2", "I'm x")   # 触发事件test2
            testobj.fire("test3", testobj)   # 触发事件test3,并传递testobj对象
            testobj.wait()  # 阻塞主线程等待时间执行结束

    输出结果::

        I'm testfunc1
        I'm testfunc2 and I'm x
        wait function will exit
    """

    def __init__(self):
        self.events = {}
        self._start = False   # 判断是否开始事件监听
        self._end = False     # 判断是否结束事件监听
        # self.semaphore = BoundedSemaphore(1)

    def __str__(self):
        if hasattr(self, "id") and hasattr(self, "name"):
            return "%s(%s)" % (getattr(self, "name"), getattr(self, "id"))
        else:
            return self.__repr__()

    @property
    def is_start(self):
        """
        只读属性,返回是否开始事件监听
        """
        return self._start

    @property
    def is_end(self):
        """
        只读属性,返回是否事结束件监听
        """
        return self._end

    def add_event(self, event, loop=False):
        """
        为对象添加单个事件处理器

        :param event: string 事件名称
        :param loop: bool  是否为可多次触发的事件
        :return: self
        """
        # if self.is_start:
        #     raise Exception("%s已经开始运行，无法添加事件!" % self)
        # if self.is_end:
        #     raise Exception("%s已经运行结束，无法添加事件!" % self)
        # with self.semaphore:
        if event in self.events.keys():
            raise Exception("事件已经存在，请勿重复添加")
        if self.check_event_name(event):
            if loop:
                self.events[event] = LoopEventHandler(event)
            else:
                self.events[event] = EventHandler(event)
            self.events[event].parent = self
        return self

    @staticmethod
    def check_event_name(name):
        """
        静态方法，检测事件名称是否符合要求

        :param name: string 事件名
        :return: bool
        """
        if not isinstance(name, types.StringTypes):
            raise Exception("事件名称必须为字符串")
        elif not name.islower():
            raise Exception("事件名称必须都是小写字母！")
        else:
            return True

    def on(self, name, func, data=None):
        """
        为对象自身绑定事件处理函数

        :param name: string 事件名称

        :param func: 自定义unbound function 或bound method(对象方法或类方法),最多允许2个参数（bound method第一个参数除外,即self）,
                      bound method定义时第一个参数应该为self或cls,如:

            .. code-block:: python
                :linenos:

                class A(object):
                    def b(self,c):
                        pass
                x=A()
                obj.bind(x.b)

        function可定义为: ``def func(a,b)`` 或 ``def func(a)`` 或 ``def func()``

        #. 第一个参数运行时获取值为**EventHandler.fire**函数传递的参数

        #. 第二个参数运行时获取值为一个**event_data**字典

               * *event_data['bind_object']* 值是该方法绑定的EventObject对象

               * *event_data['event']* 值是该方法绑定的EventHandler对象

               * *event_data['data']* 是bind方法传递的data参数

        当定义时指定参数名为 *event* 时(如 ``def func(event)`` )，此时 *event* 的值为 **event_data** 字典,不区分参数位置

        :param data: 可选参数 将通过 **event_data['data']** 形式传递给运行函数
        :return:
        """
        # with self.semaphore:
        if name not in self.events.keys():
            raise UnknownEventError(name)
        e = self.events[name]
        e.bind(func, self, data)
        return self

    def fire(self, name, para=None):
        """
        触发事件

       :param name: string 需要触发的事件名称
       :param para:  function 需要传递给事件绑定函数的参数
       :return: none
        """
        # with self.semaphore:
        if name not in self.events.keys():
            raise UnknownEventError(name)
        e = self.events[name]
        e.fire(para)

    def start_listener(self):
        """
        开始事件监听，只有开始事件监听后才能触发事件

        :return: none
        """
        if self.is_start:
            return
        else:
            if self.is_end:
                raise Exception("已经运行结束，如需重新运行，请使用restartListener方法!")
            else:
                for eve in self.events.values():
                    eve.start()
        self._start = True

    def wait(self):
        """
        阻塞主线程，直至设置为该对象为完成状态

        :return: none
        """
        while True:
            if self.is_end:
                break
            else:
                gevent.sleep(1)

    def set_end(self):
        """
        设置对象为完成状态

        :return:none
        """
        self._start = False
        self._end = True

    def stop_listener(self):
        """
        停止事件监听，将对象设置为完成状态::

            注意：此函数不能使用在自己的事件处理函数中，否则将导致异常

        :return:none
        """
        for eve in self.events.values():
            eve.stop()
        self.set_end()

    def restart_listener(self):
        """
        重置所有事件，并重新开始监听,已经停止监听的对象才能重新启动::

            注意：此函数不能使用在自己的事件处理函数中，否则将导致异常

        :return: None
        """
        if self.is_start:
            raise Exception("已经启动事件监听，使用此方法前需使用setEnd或stopListener方法结束事件监听!")
        else:
            self.stop_listener()
            for eve in self.events.values():
                eve.restart()
            self._start = True
            self._end = False
