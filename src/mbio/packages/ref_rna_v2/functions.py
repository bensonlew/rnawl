# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
from functools import wraps

def workfuncdeco(func):
    def wrapper(*args, **kwargs):
        args[0].logger.info('begin of the function ({}) at ({})'.format(func.__name__, func.__module__))
        result = func(*args, **kwargs)
        args[0].logger.info('final of the function ({}) at ({})'.format(func.__name__, func.__module__))
        return result
    return wrapper

def modlfuncdeco(func):
    def wrapper(*args, **kwargs):
        args[0].logger.info('begin of the function ({}) at ({})'.format(func.__name__, func.__module__))
        result = func(*args, **kwargs)
        args[0].logger.info('final of the function ({}) at ({})'.format(func.__name__, func.__module__))
        return result
    return wrapper

def toolfuncdeco(func):
    def wrapper(*args, **kwargs):
        args[0].logger.info('begin of the function ({}) at ({})'.format(func.__name__, func.__module__))
        result = func(*args, **kwargs)
        args[0].logger.info('final of the function ({}) at ({})'.format(func.__name__, func.__module__))
        return result
    return wrapper


def tryforgood(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper


def runcmd(tool, cmd_name, cmd, shell=False, block=True):
    if shell:
        cmd = tool.config.SOFTWARE_DIR + '/' + cmd
    command = tool.add_command(cmd_name, cmd, shell=shell)
    command.run()
    command.no_check = True
    if block:
        tool.wait()
        for name, command in tool.commands.items():
            if command.no_check:
                if command.return_code == command.default_return_code:
                    command.no_check = False
                    tool.logger.info('succeed in running {}'.format(name))
                else:
                    tool.set_error('fail to run {}, abord'.format(name))

def pkgsfuncdeco(func):
    def wrapper(*args, **kwargs):
        getime = lambda: datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print '{}\tINFO: begin of the function ({}) at ({})'.format(
            getime(), func.__name__, func.__module__
        )
        try:
            result = func(*args, **kwargs)
        except Exception as e:
            print '{}\tERROR: exception occurred when calling ({}) at ({})'.format(
                getime(), func.__name__, func.__module__
            )
            print '{}\tDEBUG: args -> {}'.format(getime(), args)
            print '{}\tDEBUG: kwargs -> {}'.format(getime(), kwargs)
            raise Exception(e)
        print '{}\tINFO: final of the function ({}) at ({})'.format(
            getime(), func.__name__, func.__module__
        )
        return result
    return wrapper
