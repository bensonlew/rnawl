#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__: konghualei 20170418

import json
import os
import sys
from collections import OrderedDict

def conf2json(conf):
    array_stack = list()
    clean = []
    a= open(conf, 'r')
    n = 0
    for line in a:
        n += 1
        line = line.strip()
        line_clean = line.replace("=>", ":")
        if "array(" in line_clean:
            if "column" in line_clean:
                line_clean = line_clean.replace("array(", "list([")
                array_stack.append("list")
            else:
                line_clean = line_clean.replace("array(", "dict({")
                array_stack.append("array(")

        # elif "(" in line_clean:
        #     array_stack.append("(")

        if ")," in line_clean or line_clean.endswith(")"):
            try:
                pair_type = array_stack.pop()
                if pair_type == "array(":
                    line_clean = line_clean.replace(")", "})")
                elif pair_type == "list":
                    line_clean = line_clean.replace(")", "])")
                else:
                    pass
            except Exception as e:
                print("line {} err".format(n))
                raise e
        clean.append(line_clean)
    with open(conf+".tmp", 'w') as f_w:
        f_w.write("{" + "\n".join(clean) + "}")
    conf_dict=eval("{" + "\n".join(clean) + "}")
    return conf_dict

one_touch_file = sys.argv[1]
deposit_file = sys.argv[2]
one_touch_dict = conf2json(one_touch_file)
deposit_dict = conf2json(deposit_file)
name2id = dict()
imgtype2id = dict()

deposit_dict_new = OrderedDict()
for k, v in deposit_dict.items():
    if "imgs" in v:
        for img, img_dict in v["imgs"].items():
            deposit_a = {
                "img_type": img_dict['method'],
                "task_id": "test",
                "img_url": "test_s3",
                "serial_no": 0,
                "submit_location": v.get("submit_location", ""),
                "tab_index": 0,
                "result_table_name": "undefined",
                "table_name": k,
                "type": "img",
                "result_table_id": "test",
                "with_table": False
            }
            name2id[img_dict['name']] = str(img)
            imgtype2id[img_dict['method']] = str(img)
            deposit_dict_new[str(img)] = deposit_a
with open("deposit.json", 'w') as f:
    f.write(json.dumps(deposit_dict_new, indent=4, ensure_ascii=False))

print json.dumps(name2id, indent=4, ensure_ascii=False)
print json.dumps(imgtype2id, indent=4, ensure_ascii=False)

one_touch_new = OrderedDict()
for k, v in one_touch_dict.items():
    onetouch_a = {
        "categroy_name": v["categroy_name"],
        "graphic_url": "/medical/show_assessmentcoverage_report.html",
        "task_id": "test",
        "result_table_id": "test",
        "img_url": "test_s3",
        "create_id": "219",
        "categroy_url": v["categroy_url"],
        "submit_location": v["submit_location"],
        "method": v["method"],
        "api": v["api"],
        "params": [],
        "table_deposit_id": "6107682acef2bb251e8b456c",
        "created_time": "2021-08-02 11:36:10",
        "delay_time": 10,
        "is_dir": "",
        "origin_table": v["origin_table"],
        "type": v["type"],
        "id": 1024,
        "result_table": v["result_table"],
        "name": v["name"]
    }

    if v["method"] in imgtype2id:
        _id = imgtype2id[v["method"]]
        print("method right {}".format(v["method"]))
    elif v["name"] in name2id:
        _id = name2id[v["name"]]
        print("right {}".format(v["name"]))
    else:
        print("error {}".format(v["name"]))
        _id = "insert{}".format(v['name'])
    one_touch_new[_id] = onetouch_a

with open("one_touch.json", 'w') as f:
    f.write(json.dumps(one_touch_new, indent=4, ensure_ascii=False))


config = OrderedDict()
for k,v in sorted(one_touch_new.items()):
    # print v
    config[k] = {
        "category_name": v["categroy_name"],
        "img": "",
        "name": v["name"],
        "main_coll": deposit_dict_new[k]["table_name"]
    }

print(json.dumps(config, indent=4, ensure_ascii=False))
