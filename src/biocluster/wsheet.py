# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""工作流表单"""

import json


class Sheet(object):
    """
    workflow 表单
    """
    def __init__(self, jsonfile="", data=None):
        """
        根据配置文件或者配置对象生成Sheet对象
        :param jsonfile: file path
        :param data:  Python object convert from json file,未指定jsonfile时生效
        :return:
        """
        self._data = {}
        if jsonfile:
            with open(jsonfile, 'r') as f:
                self._data = json.load(f)
        else:
            self._data = data

    # @property
    # def id(self):
    #     """
    #     任务ID
    #     """
    #     return self._data['id']
    #
    # @property
    # def name(self):
    #     """
    #     模块名称
    #     """
    #     if 'name' in self._data.keys():
    #         return self._data['name']
    #     else:
    #         return "link"
    #
    # @property
    # def type(self):
    #     """
    #     调用类型
    #     """
    #     if 'type' in self._data.keys():
    #         return self._data['type']
    #     else:
    #         return "link"
    #
    # @property
    # def output(self):
    #     """
    #     需要上传的远程路径
    #     :return:
    #     """
    #     if 'output' in self._data.keys():
    #         return self._data['output']
    #     else:
    #         return None

    def __getattr__(self, name):
        """
        通过下属步骤的名字直接访问下属步骤对象

        :param name:
        :return:
        """
        if name in self._data.keys():
            return self._data[name]
        else:
            return None

    @property
    def data(self):
        return self._data

    def option(self, name):
        """
        获取参数值

        :param name:  参数名
        """

        if self.type == "pipeline":
            raise Exception("pipeline类型没有参数")
        data = self._data['options']
        if name not in data.keys():
            raise Exception("没有参数%s" % name)
        return data[name]

    def set_option(self, name, value):
        """
        重新设置参数的值

        :param name:  参数名
        :param value: 参数值
        :return:
        """
        if name not in self._data['options'].keys():
            raise Exception("没有参数%s" % name)
        else:
            self._data['options'][name] = value

    def options(self):
        """
        获取所有Option

        :return: dict name/value
        """
        data = self._data['options']
        return data
