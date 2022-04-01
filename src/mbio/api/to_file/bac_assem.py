# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import os
import re
import json
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import defaultdict

client = Config().get_mongo_client(mtype="bac_assem")
db = client[Config().get_mongo_dbname("bac_assem")]


def export_log_by_gapid(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s.log" % option_name)
    collection = db["gap_fill_detail"]
    results = collection.find({"gap_id": data})
    print results.count()
    with open(file_path, "w") as f:
        for result in results:
            run_log = result["run_log"]
            run_log = run_log.split("\n")
            is_circled = "Circle" if result["scaf_status"] == "Circular" else "linear"
            log_list = []
            for log in run_log:
                if log.startswith("Scaffold"):
                    log_list.append(log)
            f.write(result["scaf_name"] + "\t" + ",".join(log_list) + "\t" + is_circled + "\n")
    return file_path

def export_log_by_sof(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s.log" % option_name)
    collection1 = db["draft"]
    print(bind_obj.sheet.option("task_id"))
    main = collection1.find_one({'task_id': bind_obj.sheet.option("task_id")})
    sample = bind_obj.sheet.option("samp")
    print(main)
    print(main['_id'])
    main_id = main['_id']
    collection = db["draft_seq_detail"]
    results = collection.find({"draft_id": main_id, "sof_type": data, "type": "scaffold", "samp": sample})
    print results.count()
    with open(file_path, "w") as f:
        for result in results:
            is_circled = "Circle" if result["status"] == "Circular" else "linear"
            f.write(result["seq_id"] + "\t" + result["seq_id"] + "\t" + is_circled + "\n")
    return file_path