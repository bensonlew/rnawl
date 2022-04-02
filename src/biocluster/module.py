# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""
功能基本类
"""
from .core.function import load_class_by_path, get_classpath_by_object
from .basic import Basic
import time
import copy


class Module(Basic):
    """
    功能基本类
    """
    def __init__(self, parent):
        super(Module, self).__init__(parent)
        self._tools_report_data = []
        self._start_time = None
        self._end_time = None

    def add_tool_report(self, data):
        """
        添加Tool运行报告
        :return:
        """
        self._tools_report_data.append(data)

    def get_report_data(self):
        """
        获取运行报告
        :return:
        """

        data = {
            "run_id": self.id,
            "path": get_classpath_by_object(self),
            "work_dir": self.work_dir,
            "start_time": self._start_time if self._start_time else int(time.time()),
            "end_time": self._end_time if self._end_time else int(time.time()),
            "tool_num": len(self.children),
            "tools": self._tools_report_data,
            "modules": []
        }
        for c in self.children:
            if isinstance(c, Module):
                data["modules"].append(c.get_report_data())
        return data

    def add_tool(self, path):
        """
        添加下属工具，跳过模块

        :param path: String   Tool路径名称
        :return:  Agent 返回Tool对应的 :py:class:`biocluster.agent.Agent` 对象
        """
        agent = load_class_by_path(path, tp="Agent")(self)
        self.add_child(agent)
        return agent

    def find_tool_by_id(self, toolid):
        """
        通过id搜索所属Tool

        :param toolid:  :py:class:`biocluster.tool.Tool` 对象的ID
        """
        # ids = toolid.split(".")
        # if len(ids) < 3:
        #     return False
        # if (ids[0] + "." + ids[1]) != self.id:
        #     return False
        childs = copy.copy(self.children)
        for c in childs:
            if str(c.id) == str(toolid):
                return c
            else:
                if hasattr(c, "find_tool_by_id"):
                    tool = c.find_tool_by_id(toolid)
                    if tool:
                        return tool
        return False

    def add_module(self, path):
        """
        添加下属 :py:class:`biocluster.module.Module`

        :param path:  :py:class:`biocluster.module.Module` 对应的自动加载path路径，请参考教程中对应的说明

        """
        module = load_class_by_path(path, tp="Module")(self)
        self.add_child(module)
        return module

    def end(self):
        self._end_time = int(time.time())
        super(Module, self).end()

    def run(self):
        self._start_time = int(time.time())
        super(Module, self).run()

    def count_used(self):
        childs = self.children
        cpu_used = 0
        memory_used = 0
        for c in childs:
            (x, y) = c.count_used()
            cpu_used += x
            memory_used += y
        return cpu_used, memory_used
