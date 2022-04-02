# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
import os
import json


def json_to_file(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_input.json" % option_name)
    bind_obj.logger.debug("正在将参数%s转换为文件，路径:%s" % (option_name, file_path))
    with open(file_path, "w") as f:
        json.dump(data, f, indent=4)

    return file_path


def text_to_file(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_input.json" % option_name)
    bind_obj.logger.debug("正在将参数%s转换为文件，路径:%s" % (option_name, file_path))
    with open(file_path, "w") as f:
        f.write(data)
    return file_path
