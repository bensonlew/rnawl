# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# last modified @ 20210316
import json


def params_check(toollab_params):
    if ('fasta' not in toollab_params) and ('fasta_dir' not in toollab_params):
        info = {"success": False, "info": "必须输入单个基因组或者多个基因组文件夹"}
        return json.dumps(info)
    if toollab_params["species"] not in ["Bacteria", "Archaea", "Viruses"]:
        info = {"success": False, "info": "species参数{}必须为Bacteria, Archaea or Viruses".format(toollab_params["species"])}
        return json.dumps(info)
    return None
