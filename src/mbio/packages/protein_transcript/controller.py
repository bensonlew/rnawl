# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')

from run_commands import RunCommands
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import time


class Controller(object):
    def __init__(self):
        self.command_run_dict = defaultdict(list)
        self.func_run_dict = defaultdict(list)
        self.command_list = list()
        self.function_list = list()

    def end(self, normal=0, out = None):
        if not normal:
            print('your program finished normally')
        else:
            print('your program finished unexceptally')
        if out:
            print(out)
        os._exit(normal)

    def add_command(self, params = None, name = None):
        if not params:
            params = dict()
        if type(params) != dict:
            self.end(normal=1, out='the params of your command must be wrong')
        if not 'run_wd' in params:
            try:
                run_wd = self.output
            except:
                run_wd = os.getcwd()
            params.update({'run_wd': run_wd, 'name': name})
        else:
            params.update({'name': name})
        command = RunCommands(params)
        self.command_list.append(command)
        return command

    def command_run(self, c):
        if not c in self.command_list:
            self.end(normal=1,out='your command has not been added')
        if not os.path.exists(c.running) and not os.path.exists(c.judge):
            c.run()

    def function_run(self, func):
        if func in self.function_list:
            return
        if not hasattr(self, func.__name__):
            self.end(normal=1, out='your function %s has not been defined' %func)
        can = 1
        if func in self.func_run_dict:
            for c in self.func_run_dict[func]:
                if not os.path.exists(c.judge):
                    can = 0
                    break
        if can:
            self.function_list.append(func)
            func = getattr(self, func.__name__)
            try:
                print('%s is running'%func.__name__)
                func()
                print('%s runned' % func.__name__)
            except Exception as e:
                print(e)
            else:
                print('%s run successful'%func.__name__)


    def init_func(self, func):
        if not hasattr(self, func.__name__):
            self.end(normal=1, out='your function %s has not been defined' %func)
        while func in self.function_list:
            self.function_list.remove(func)

    def after_one(self, c, func):
        if not c in self.command_list:
            self.end(normal=1, out='your command %s has not been added' %c)
        if not hasattr(self, func.__name__):
            self.end(normal=1, out='your function %s has not been defined' %func)
        self.command_run_dict[c].append(func)
        self.func_run_dict[func].append(c)

    def after_some(self, cs, func):
        if type(cs) != list:
            cs = [cs]
        if not hasattr(self, func.__name__):
            self.end(normal=1, out='your function %s has not been defined' %func)
        for c in cs:
            if not c in self.command_list:
                self.end(normal=1, out='your command %s has not been added' % c)
            else:
                self.command_run_dict[c].append(func)
        self.func_run_dict[func] += cs

    def fire(self):
        while 1:
            for c in self.command_run_dict:
                if os.path.exists(c.judge):
                    with ThreadPoolExecutor(6) as pool:
                        pool.map(self.function_run, self.command_run_dict[c])
            time.sleep(20)

