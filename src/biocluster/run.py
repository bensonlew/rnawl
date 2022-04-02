# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""
功能基本类
"""
from .core.function import load_class_by_path, get_classpath_by_object
from .basic import Basic
import time
import copy
import argparse
import json


class Run(object):
    """
    功能基本类
    """
    def __init__(self, tp="Tool", name=None):
        self.tp = tp
        self.name = ""

    def run(self, obj_json={}):
        obj = load_class_by_path(self.name, self.tp)
        if self.tp == 'Tool':
            '''
            tool 使用agent的检查方法检查
            '''
            agent = load_class_by_path(self.name, "Agent")
            agent.check_options()

        obj.set_options(obj_json["options"])
        obj.run()

    def set_options(self, options):
        """
        批量设置参数值

        :param options: dict key是参数名,value是参数值
        :return:
        """
        if not isinstance(options, dict):
            raise Exception("参数格式错误!")
        for name, value in options.items():
            self.option(name, value)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-tp", type=str, required=True, 
                        help = "one in Workflow, Module, Tool")
    parser.add_argument("-name", type=str, required=True,
                        help = "path joined by ."
                        )
    parser.add_argument("-json", type=str, required=False,
                        help = "paramters for tools and so on"
                        )
    args = parser.parse_args()
    run = Run(tp=args.tp, name=args.name)
    run.run()
    

        



