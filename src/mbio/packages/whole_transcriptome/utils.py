# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import os
from collections import OrderedDict


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


def runcmd(tool, cmd_name, cmd, shell=False, block=True, ignore_error=False):
    if shell:
        cmd = tool.config.SOFTWARE_DIR + '/' + cmd
    if ignore_error:
        command = tool.add_command(cmd_name, cmd, shell=shell, ignore_error=True)
    else:
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
        print '{}\tDEBUG: args -> {}'.format(getime(), args)
        print '{}\tDEBUG: kwargs -> {}'.format(getime(), kwargs)
        result = func(*args, **kwargs)
        print '{}\tINFO: final of the function ({}) at ({})'.format(
            getime(), func.__name__, func.__module__
        )
        return result

    return wrapper


def check_map_dict(map_dict):
    for k, v in map_dict.items():
        if not os.path.exists(v):
            raise Exception('can not find {} with label ({}) in map dict'.format(v, k))


def read_fastq_dir(fastq_dir):
    is_se = False
    fastqs_dict = OrderedDict()
    for line in open(os.path.join(fastq_dir, 'list.txt')):
        eles = line.strip().split('\t')
        fastq = os.path.join(fastq_dir, eles[0])
        sample, mate_type = eles[1:]
        if sample in fastqs_dict:
            if mate_type == 'l':
                fastqs_dict[sample].insert(0, fastq)
            elif mate_type == 'r':
                fastqs_dict[sample].append(fastq)
        else:
            fastqs_dict[sample] = [fastq]
    else:
        pe_sample_count = len(filter(lambda item: len(item[1]) > 1, fastqs_dict.items()))
        if pe_sample_count:
            if pe_sample_count == len(fastqs_dict):
                is_se = False
            else:
                raise Exception('mix mate type found in {} -> {}'.format(fastq_dir, fastqs_dict))
        else:
            is_se = True
        return is_se, fastqs_dict
