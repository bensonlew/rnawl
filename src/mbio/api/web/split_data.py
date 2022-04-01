# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from biocluster.wpm.log import Log


class SplitData(Log):

    def __init__(self, data):
        super(SplitData, self).__init__(data)
        # self._client = "client01"
        # self._key = "1ZYw71APsQ"
        self._url = "http://172.16.6.15:8080/api/split/receive_pipeline"
        # self._url = "http://172.16.6.96/html/code.php"
        self._post_data = self.post_data

    def update(self):
        self.send()