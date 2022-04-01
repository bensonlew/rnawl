# -*- coding: utf-8 -*-
# __author__ = 'guoquan'


class Workflow(object):
    def __init__(self):
        """
        因为很多地方使用到该类中get_by_workflow_id方法去获取task_id并初始化为唯一的值，因为新框架中不支持去查找Postgresql。
        所以如果您发现你因为task_id重复导致的报错，可以尝试如下方法：
        方法一. new_id = '{}_{}'.format(task_id, datetime.datetime.now().strftime("%m%d%H%M%S%f"),
        random.randint(1000, 10000))  时间精确到毫秒，一般情况是不会有问题了
        方法二. 使用mainapp/controllers/core/base中的get_new_id方法去组建task_id
        """
        pass

    def get_by_workflow_id(self, wid):
        return []
