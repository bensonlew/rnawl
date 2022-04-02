# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""单例模式"""


def singleton(cls, *args, **kw):
    """

    定义单例模式

    :param cls:
    :param args:
    :param kw:
    :return:
    """
    instances = {}

    def _singleton(*args, **kw):
        if cls not in instances:
            instances[cls] = cls(*args, **kw)
        return instances[cls]
    return _singleton

