# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

from biocluster.wpm.client import worker_client, wait


class Instant(object):
    def __init__(self, data):
        self._json = data
        self._id = data["id"]
        self._return_msg = None

    @property
    def id(self):
        """
        获取运行任务的ID

        :return:
        """
        return self._id

    @property
    def return_msg(self):
        """
        获取运行任务的返回值

        :return:
        """
        return self._return_msg

    def run(self):
        worker = worker_client()
        info = worker.add_task(self._json)
        if info["success"]:
            end = wait(self._id)
            if end is True:
                self._return_msg = worker.get_msg(self._id)
            else:
                raise Exception("运行超时!")
        else:
            raise Exception("任务提交失败:%s" % info["info"])
