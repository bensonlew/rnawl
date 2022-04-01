# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import logging

def reminder(func):
    def wrapper(*args, **kwargs):
        args[0].logger.info('begin of the function ({}) at ({})'.format(func.__name__, func.__module__))
        result = func(*args, **kwargs)
        args[0].logger.info('final of the function ({}) at ({})'.format(func.__name__, func.__module__))
        return result
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
                    tool.logger.info('succeed in running command ({})'.format(name))
                else:
                    tool.set_error('fail to run command ({}), abord'.format(name))

def watcher(func):
    logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    def wrapper(*args, **kwargs):
        logging.debug('args -> {}'.format(args))
        logging.debug('kwargs -> {}'.format(kwargs))
        logging.info('begin of the function ({}) at ({})'.format(func.__name__, func.__module__))
        result = func(*args, **kwargs)
        logging.info('final of the function ({}) at ({})'.format(func.__name__, func.__module__))
        return result
    return wrapper
