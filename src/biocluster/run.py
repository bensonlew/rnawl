# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

"""
功能基本类
"""
import imp
from biocluster.core.function import load_class_by_path, get_classpath_by_object
from biocluster.basic import Basic
from biocluster.config import Config
import time
import copy
import argparse
import json
import os
import pickle
import sys
from biocluster.wsheet import Sheet
from biocluster.core.exceptions import ExitError

PY3 = False



class Run(Basic):
    """
    功能基本类
    """
    def __init__(self, tp="Tool", name=None, sheet=None):
        super(Run, self).__init__()
        self.tp = tp.capitalize()
        
        self.pathname = name
        self._full_name = "Run"
        # self.parent = None
        self.sheet = sheet
        self.sheet_update()
        self._work_dir = "./"
        self.debug = False
        
        if self.tp=="Tool":
            self.agent = load_class_by_path(self.pathname, "Agent")(self)
        elif self.tp=="Workflow":
            a  = load_class_by_path(self.pathname, "Workflow")
            self.workflow = a(self.sheet)
        else:
            self.obj = load_class_by_path(self.pathname, self.tp)(self)
            
    def sheet_update(self):
        self.sheet._data.update({
                "PACKAGE_DIR": "/home/liubinxu/work/rnawl/src/mbio/packages",
                "SOFTWARE_DIR": "/home/liubinxu/miniconda2/bin"

            }

        )
        self.sheet.work_dir  = "./"
        self.sheet.id = "test"

    def _update(self, error_type, error_str):
        pass
            
    def get_options(self):
        return self.obj._options

    def exit(self, exitcode=1, data=None, terminated=False):
        """
        立即退出当前流程

        :param exitcode:
        :param data:
        :param terminated:
        :return:
        """
        # exstr = traceback.format_exc()
        # print exstr
        sys.stdout.flush()
        print("data {}".format(data))
        # self.end_unfinish_job()
        if isinstance(data, dict) and "error_type" in data.keys() and "info" in data.keys():
            error_str = json.dumps(data)
            if PY3:
                error_str = bytes(json.dumps(data), encoding='utf-8')
            # self._save_report_data(data)
        else:
            error_str = str(data)
            # self._save_report_data()
        if terminated:
            self.step.terminated(data)
        else:
            self.step.failed(data)
        self.step.update()
        self._update("error", error_str)  # 退出流程不再给wfm发送error状态
        self.logger.error("Failed: %s " % error_str)
        # self.rpc_server.close()
        print("Failed: %s " % error_str)
        raise ExitError(error_str)

    @property
    def step(self):
        """
        主步骤

        :return:
        """
        return self._main_step

    def run(self):
        if self.tp == 'Tool':
            '''
            tool 使用agent的检查方法检查
            '''
            agent = load_class_by_path(self.pathname, "Agent")(self)
            print(agent)
            options = self.sheet.options()
            agent.set_options(options)
            agent.check_options()
            # path = os.path.join(self.workdir, self.name + ".pk")
            path = agent.save_config()
            print("path :%s", path)
            with open(path, "rb") as f:
                pickle_config = pickle.load(f)

            self.obj = load_class_by_path(self.pathname, "Tool")(pickle_config)
            self.obj.run()
        elif self.tp == 'Workflow':
            self.workflow.run()
        else:
            self.obj.set_options(options)
            self.obj.run()

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
    parser.add_argument("-tp", type=str, required=False, 
                        help = "one in Workflow, Module, Tool")
    parser.add_argument("-name", type=str, required=False,
                        help = "path joined by ."
                        )
    parser.add_argument("-par", type=bool, default=False)
    parser.add_argument("-json", type=str, required=False,
                        help = "paramters for tools and so on",
                        default = None
                        )
    parser.add_argument("-options", type=str, required=False,
                        help = "options like para1=para1 para2=para2",
                        default = None
                        )
    args = parser.parse_args()
    print(args)

    sheet = dict()

    if args.json:
        with open(args.json, "r") as f:
            sheet = json.load(f)


    ## 如果添加参数根据参数替换
    if args.tp and args.name:
        sheet["type"] = args.tp
        sheet["name"] = args.name



    if args.options:
            for opt in args.options.split(" "):
                for k,v in opt.split("="):
                    sheet["options"].update({k: v})
        
    run = Run(tp = sheet["type"], name = sheet["name"], sheet = Sheet(jsonfile = args.json))

    if args.par:
        opt = run.get_options()
        print(json.dumps(opt, indent=4))
    else:    
        run.run()
    

        



