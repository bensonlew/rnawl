# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import gevent
from .singleton import singleton
import traceback
import sys
from .exceptions import CodeError


@singleton
class Watcher(object):
    """
    管理定时循环运行的监控任务，
    """
    def __init__(self):
        self._jobs = []

    def add(self, func,  interval=1):
        """
        添加需要定时重运行的函数，并设置运行时间间隔
        当函数返回值为exit时，退出循环运行

        :param func:  需要执行的函数
        :param interval:  运行时间间隔
        :return: None
        """
        def watcher():
            # print "add watcher %s " % func
            while True:
                # print "%s is start" % func
                # print "start watcher %s" % os.getpid()
                try:
                    r = func()
                except CodeError as e:
                    raise e
                except Exception as e:
                    ex_str = traceback.format_exc()
                    print("执行定时任务出错，即将重试...")
                    print(ex_str, e)
                    sys.stdout.flush()
                    sys.stderr.flush()
                else:
                    if r and r == "exit":
                        # print "%s exit" % func
                        return
                gevent.sleep(interval)
                # print "end watcher %s" % os.getpid()
        self._jobs.append(gevent.spawn(watcher))

    def stopall(self):
        """
        终止所有正在运行的监控函数运行

        :return:
        """
        if len(self._jobs) > 0:
            gevent.killall(self._jobs)
