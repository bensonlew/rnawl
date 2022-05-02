# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from multiprocessing.managers import BaseManager
from ..config import Config

config = Config()


class LogManager(BaseManager):
    pass


class WorkerManager(BaseManager):
    pass


def worker_client():
    WorkerManager.register("worker")
    wpm_manager = WorkerManager(address=config.wpm_listen, authkey=config.wpm_authkey)
    wpm_manager.connect()
    worker = wpm_manager.worker()
    return worker


def event_client(*wid):
    WorkerManager.register("get_event")
    wpm_manager = WorkerManager(address=config.wpm_listen, authkey=config.wpm_authkey)
    wpm_manager.connect()
    return wpm_manager.get_event(*wid)


def wait(wid, timeout=None):
    if isinstance(wid, list) or isinstance(wid, tuple):
        event = event_client(*wid)
    else:
        event = event_client(wid)
    if str(event) == "None":
        raise Exception("获取Event错误!")
    if timeout is None:
        timeout = config.wpm_instant_timeout
        return event.wait(timeout)
    elif timeout == 0:
        return event.wait()
    else:
        return event.wait(timeout)


def log_client():
    LogManager.register("apilog")
    m = LogManager(address=config.wpm_logger_listen, authkey=config.wpm_logger_authkey)
    m.connect()
    log = m.apilog()
    return log
