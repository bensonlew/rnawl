# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import pandas as pd
import os
import functools
import shutil
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base
import colors
from colors import darken, is_foreground_light, lighten
from pymongo import MongoClient
from biocluster.config import Config
import json

# 文件移动、运行时间相关方法


def get_old_mongo():
    '''
    固定的方法获取老mongo服务器
    '''
    client = MongoClient("mongodb://meta:v6m4t7w9y6x5@10.100.1.10/sanger?authMechanism=SCRAM-SHA-1",connect=False)
    db = client['sanger_metagenomic']
    return client, db


def get_mongo():
    client = Config().get_mongo_client(mtype="metagenomic")
    db = client[Config().get_mongo_dbname("metagenomic")]
    return client, db

def get_save_pdf_status(task_id=None):
    """
    查询任务save_pdf字段
    :param task_id: task_id
    :return:

    mongo, db = get_mongo()
    collection = db['sg_task']
    result = collection.find_one({"task_id": task_id})
    if result:
        if "save_pdf" in result and result["save_pdf"] == 1:
            save_pdf = True
        else:
            save_pdf = False
    else:
        save_pdf = False
    # return save_pdf
    """
    return True

def get_name(table_id=None, table=None, name="name"):
    """
    查询主表name字段，用于pdf存图
    """
    mongo, db = get_mongo()
    collection = db[table]
    result = collection.find_one({"_id": ObjectId(table_id)})
    return result[name]

def get_submit_loc(table_id=None, table=None):
    """
    查询主表submit_location字段，用于pdf存图
    """
    mongo, db = get_mongo()
    collection = db[table]
    result = collection.find_one({"_id": ObjectId(table_id)})
    params = result["params"]
    submit_location = json.loads(params)["submit_location"]
    return submit_location

def link_file(oldfile, newfile):
    """
    hard link file from oldfile to newfile
    :param oldfile:
    :param newfile:
    :return:
    """
    if not os.path.isfile(oldfile):
        raise Exception("不存在文件：%s" % oldfile)
    if os.path.exists(newfile):
        os.remove(newfile)
    os.link(oldfile, newfile)

def link_dir(olddir, newdir):
    """
    hard link directory from olddir to newdir
    :param olddir:
    :param newdir:
    :return:
    """
    if not os.path.isdir(olddir):
        raise Exception("不存在路径: %s" % olddir)
    allfiles = os.listdir(olddir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    # else:
    #     shutil.rmtree(newdir)
    #     os.makedirs(newdir)
    for file in allfiles:
        oldfile = os.path.join(olddir, file)
        newfile = os.path.join(newdir, file)
        if os.path.isdir(newfile):
            shutil.rmtree(newfile)
        if os.path.isfile(oldfile):
            link_file(oldfile, newfile)
        else:
            link_dir(oldfile, newfile)

def time_count(func):  # 统计函数运行时间，作为方法的装饰器
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run %s at %s' % (func_name, start_time))
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，共运行{}s".format(func_name, end - start))
    return wrapper

def wait_file(path, times=10, sleep=10):
    '''
    等待某个文件生成
    :param path: 文件名称
    :param wait_times: 最大等待次数
    :param sleep: 每次等待时间，单位秒
    :return:
    '''

    while times > 0:
        if not os.path.isfile(path):
            time.sleep(sleep)
            times -= 1
            self.wait_file(path, times=times, sleep=sleep)
        return path
    raise Exception("超过文件等待次数，需检查文件%s" % path)

# 代谢组处理相关方法

def check_info(info):
    '''
    检查原始表中的注释内容，将代表无注释的字符统一处理成'-'
    :param info: 注释信息
    :return:  str
    '''
    if info in ['_', '-', ' ', '', '\r']:
        return '-'
    else:
        return info

def check_command(bind_tool, command, success_list, memory_limit_list=None, success_fun=None, fail_fun=None, memory_limit_fun=None):
    """
    检查Tool里面调用的command返回值，根据返回值判断具体处理
    :param bind_tool: Tool
    :param command: Command
    :param success_list: app运行成功的return_code list
    :param memory_limit_list: app运行超内存时的return_code list
    :param success_fun: 运行成功时采取的行为
    :param fail_fun: 运行失败后采取的行为
    :param memory_limit_fun: 超出内存时采取的行为
    :return:
    """
    bind_tool.logger.info("<--START CHECK COMMAND %s-->" % command.name)
    if type(success_list) != type(list()):
        if success_list is None:
            success_list = list()
        elif success_list == 0:
            success_list = [0]
        else:
            bind_tool.set_error("command %s check_command success_list arg error: %s" % (command.name, success_list))
    if type(memory_limit_list) != type(list()):
        if memory_limit_list is None:
            memory_limit_list = list()
        else:
            bind_tool.set_error("command %s check_command memory_limit_list arg error: %s" % (command.name, memory_limit_list))
    if command.return_code in success_list:
        if success_fun:
            success_fun()
        else:
            bind_tool.logger.info("command %s success" % command.name)
    elif command.return_code in memory_limit_list:
        if memory_limit_fun:
            memory_limit_fun()
        else:
            bind_tool.add_state("memory_limit", "memory is low!")
    else:
        if fail_fun:
            fail_fun()
        else:
            bind_tool.set_error("command %s failed" % command.name)

    bind_tool.logger.info("<--START CHECK COMMAND END-->")


class Style(object):
    """basic style"""
    colors = (
        '#F44336',  # 0
        '#3F51B5',  # 4
        '#009688',  # 8
        '#FFC107',  # 13
        '#FF5722',  # 15
        '#9C27B0',  # 2
        '#03A9F4',  # 6
        '#8BC34A',  # 10
        '#FF9800',  # 14
        '#E91E63',  # 1
        '#2196F3',  # 5
        '#4CAF50',  # 9
        '#FFEB3B',  # 12
        '#673AB7',  # 3
        '#00BCD4',  # 7
        '#CDDC39',  # 11b
        '#9E9E9E',  # 17
        '#607D8B',  # 18
    )

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def get_colors(self, len_):
        colors = []
        if len(self.colors) < len_:
            missing = len_ - len(self.colors)
            cycles = 1 + missing // len(self.colors)
            for i in range(0, cycles + 1):
                for color_ in self.colors:
                    colors.append(darken(color_, 33 * i / cycles))
                    if len(colors) >= len_:
                        break
                else:
                    continue
                break
        else:
            colors = self.colors[:len_]
        return colors

DefaultStyle = Style

class DarkStyle(Style):
    """A dark style (old default)"""
    # colors = (
    #     '#ff5995', '#b6e354', '#feed6c', '#8cedff', '#9e6ffe', '#899ca1',
    #     '#f8f8f2', '#bf4646', '#516083', '#f92672', '#82b414', '#fd971f',
    #     '#56c2d6', '#808384', '#8c54fe', '#465457'
    # )
    # 以下为冯超确认的配色方案 20190122
    colors = (
        '#82b414','#ff5995',             '#feed6c', '#8cedff', '#9e6ffe', '#899ca1',
                   '#bf4646', '#516083', '#f92672',
        '#56c2d6', '#808384', '#8c54fe', '#465457'
    )

class LightStyle(Style):
    """A light style"""
    colors = (
        '#242424', '#9f6767', '#92ac68', '#d0d293', '#9aacc3', '#bb77a4',
        '#77bbb5', '#777777'
    )

class DarkSolarizedStyle(Style):
    """Dark solarized popular theme"""
    colors = (
        '#b58900', '#cb4b16', '#dc322f', '#d33682', '#6c71c4', '#268bd2',
        '#2aa198', '#859900'
    )

class RedBlueStyle(Style):
    """A red and blue theme"""
    colors = (
        '#d94e4c', '#e5884f', '#39929a', lighten('#d94e4c', 10),
        darken('#39929a', 15), lighten('#e5884f', 17), darken('#d94e4c', 10),
        '#234547'
    )

class LightColorizedStyle(Style):
    """A light colorized style"""
    colors = (
        '#fe9592', '#534f4c', '#3ac2c0', '#a2a7a1', darken('#fe9592', 15),
        lighten('#534f4c', 15), lighten('#3ac2c0', 15), lighten('#a2a7a1', 15),
        lighten('#fe9592', 15), darken('#3ac2c0', 10)
    )

class DarkColorizedStyle(Style):
    """A dark colorized style"""
    colors = (
        '#c900fe', '#01b8fe', '#59f500', '#ff00e4', '#f9fa00',
        darken('#c900fe', 20), darken('#01b8fe', 15), darken('#59f500', 20),
        darken('#ff00e4', 15), lighten('#f9fa00', 20)
    )

class TurquoiseStyle(Style):
    """A turquoise style"""
    colors = (
        '#93d2d9', '#ef940f', '#8C6243', '#fff', darken('#93d2d9', 20),
        lighten('#ef940f', 15), lighten('#8c6243', 15), '#1b8088'
    )

class LightGreenStyle(Style):
    """A light green style"""
    colors = (
        '#7dcf30', '#247fab', lighten('#7dcf30', 10), '#ccc',
        darken('#7dcf30', 15), '#ddd', lighten('#247fab', 10),
        darken('#247fab', 15)
    )

class DarkGreenStyle(Style):
    """A dark green style"""
    colors = (
        '#adde09', '#6e8c06', '#4a5e04', '#fcd202', '#C1E34D',
        lighten('#fcd202', 25)
    )

class DarkGreenBlueStyle(Style):
    """A dark green and blue style"""
    colors = (
        lighten('#34B8F7', 15), '#7dcf30', '#247fab', darken('#7dcf30', 10),
        lighten('#247fab', 10), lighten('#7dcf30', 10), darken('#247fab', 10),
        '#fff'
    )

class BlueStyle(Style):
    """A blue style"""
    colors = (
        '#00b2f0', '#43d9be', '#0662ab', darken('#00b2f0', 20),
        lighten('#43d9be', 20), lighten('#7dcf30', 10), darken('#0662ab', 15),
        '#ffd541', '#7dcf30', lighten('#00b2f0', 15), darken('#ffd541', 20)
    )

class SolidColorStyle(Style):
    """A light style with strong colors"""
    colors = (
        '#FF9900', '#DC3912', '#4674D1', '#109618', '#990099', '#0099C6',
        '#DD4477', '#74B217', '#B82E2E', '#316395', '#994499'
    )

class MetaDefine1(Style):
    """冯超挑选的颜色20190222"""
    colors = (
        '#82b414', '#FF0000', '#0000FF', '#993300', '#FF9900', '#9900FF',
        '#FF6600', '#FF00FF', '#00FFFF'
    )

styles = {
    'default': DefaultStyle,  # basic style
    'dark': DarkStyle,  # A dark style (old default)
    'light': LightStyle,  # A light style
    'light_red_blue': RedBlueStyle,  # A red and blue theme
    'dark_solarized': DarkSolarizedStyle,  # Dark solarized popular theme
    'dark_colorized': DarkColorizedStyle,  # A dark colorized style
    'light_colorized': LightColorizedStyle,  # A light colorized style
    'turquoise': TurquoiseStyle,  # A turquoise style
    'green': LightGreenStyle,  # A light green style
    'dark_green': DarkGreenStyle,  # A dark green style
    'dark_green_blue': DarkGreenBlueStyle,  # A dark green and blue style
    'blue': BlueStyle,  # A blue style
    'solid_color': SolidColorStyle,  # A light style with strong colors
    'define1': MetaDefine1  # ah ooh, pure color by fengchao
}

# class ParametricStyleBase(Style): pygal包的style包含的类，有需要再加，根据一个初始化的颜色生成后续颜色，基于style