# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last modified @ 20200421
import json
import os
from biocluster.file import exists
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.file import download


def params_check(toollab_params):
    if toollab_params["target_field"].split(",")[-1] in ["all", "mixed"]:
        if type(toollab_params["search_string"]) == dict:
            path = toollab_params["search_string"]["path"]
        elif type(toollab_params["search_string"]) == str:
            path = toollab_params["search_string"]
        else:
            info = {"success": False, "info": "文件格式错误"}
        id_file = download_from_s3(path)
        with open(id_file, 'r') as f:
            line_num = len(f.readlines())
            if line_num > 200:
                info = {"success": False, "info": "id 数量大于200时不可以选择\"For More Information\""}
                return json.dumps(info)
    return None


def download_from_s3(from_file, to_path="download/", inter_dir="", cover=True):
    base_name = os.path.basename(from_file)
    to_file = os.path.join(inter_dir, base_name)
    #print('from {} to {}'.format(from_file, to_file))
    if os.path.exists(to_file) and os.path.getsize(to_file) != 0:
        pass
    else:
        download(from_file, to_file)
    return to_file
