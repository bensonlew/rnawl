# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last modify 20180417

import os
import json
from bson.objectid import ObjectId
from biocluster.config import Config


db = Config().get_mongo_client(mtype="dna_wgs")[Config().get_mongo_dbname("dna_wgs")]


def get_venn_params(id):
    snp_coll = db['sg_snp_compare']
    result = snp_coll.find_one({"_id": ObjectId(id)})
    main_id = result["main_id"]
    type = result["type"]
    name = result["name"]
    params = json.loads(result["params"])
    if type == "sample":
        name1 = params["sample1"] + "_gt"
        name2 = params["sample2"] + "_gt"
    elif type == "two_group":
        group_names = params["group_dict"].keys()
        name1 = group_names[0] + "_frequence"
        name2 = group_names[1] + "_frequence"
    else:
        group_names = params["group_dict"].keys()
        name1 = group_names[0] + "_frequence"
        name2 = None
    return main_id, name, name1, name2


def export_venn_file(data, option_name, dir_path, bind_obj=None):
    # bind_obj.logger.debug("正在导出")
    snp_coll = db['sg_snp_compare']
    snp_diff_coll = db["sg_snp_compare_detail"]
    diff_dict = {}
    new_data = data
    all_list = []
    for id in data:
        new_data.remove(id)
        print new_data
        main_id, name, name1, name2 = get_venn_params(id)
        venn_path = os.path.join(dir_path, "%s_venn.list" % name)
        venn_file = open(venn_path, "w")
        diff_path = os.path.join(dir_path, "%s_venn_variant.xls" % name)
        diff_file = open(diff_path, "w")
        if name2:
            header = "#chrom\tpos\tref\talt\t" + name1 + "\t" + name2
        else:
            header = "#chrom\tpos\tref\talt\t" + name1
        diff_file.write(header + "\n")
        results = snp_diff_coll.find({"compare_id": main_id})
        for result in results:
            chr = result["chr"]
            pos = result["pos"]
            ref = result["ref"]
            venn = chr + "_" + str(pos) + "_" + ref
            venn_file.write(venn + "\n")
            if venn not in all_list:
                all_list.append(venn)
            line = chr + "\t" + str(pos) + "\t" + ref + "\t" + result[name1]
            if name2 in result.keys():
                line += "\t" + result[name2]
            diff_file.write(line + "\n")


data = ["5ad43a8ba4e1af305c6bfc0f", "5ad55d81a4e1af014bd55d30", "5ad57c4da4e1af247b66c499"]
option_name = ""
dir_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/venn"
export_venn_file(data, option_name, dir_path, bind_obj=None)
