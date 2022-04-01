# -*- coding: utf-8 -*-

"""
该模块中保存不同产品controller中公用的一些函数方法
"""
import datetime
import time


def get_new_id(task_id, main_id=None, analysis_name=None):
    """
    按照时间去生成task_id
    :param task_id: 基础工作流的任务id
    :param main_id: 交互分析主表的id， 一般是update_status中的key
    :param analysis_name: 分析自定义的名字
    :return:
    """
    # time.time()的结果，保留精度到20位
    tmpid = "".join(("%.20f" % time.time()).split("."))
    if main_id:
        tmpid = main_id
    if analysis_name:
        new_id = '{}_{}_{}'.format(task_id, analysis_name + tmpid, datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
    else:
        new_id = '{}_{}_{}'.format(task_id, tmpid, datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
    return new_id
