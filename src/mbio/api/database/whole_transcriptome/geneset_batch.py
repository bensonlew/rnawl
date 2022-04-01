# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'


import datetime
import os
import json
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from bson.objectid import ObjectId


class GenesetBatch(ApiBase):
    def __init__(self, bind_object):
        super(GenesetBatch, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    def add_diff_genesets(self, task_id, geneset_file, category, geneset_type, geneset_suffix, diff_id):
        if not os.path.exists(geneset_file):
            raise Exception("文件不存在{}".format(geneset_file))
        main_info = self.db["diff"].find_one({"main_id": ObjectId(diff_id)})
        project_sn = main_info["project_sn"]
        diff_name = main_info["name"]
        try:
            group_id = main_info["group_id"]
        except:
            group_id = json.loads(main_info["params"])["group_id"]
        with open(geneset_file, "r") as f:
            head = f.readline()
            for line in f:
                items = line.strip().split("\t")
                if len(items) <= 2:
                    continue
                if not (items[0] or items[1]):
                    continue
                name = items[0]
                name = geneset_type + "_" + category + "_" + name + geneset_suffix
                seq_list = items[1].split(",")
                regulate_list = items[2].split(",")
                kind_list = items[3].split(",")
                category_list = items[4].split(",")
                type = ""
                if int(category_list.count("mRNA")) >= 1:
                    type += "mRNA:" + str(category_list.count("mRNA"))
                if int(category_list.count("lncRNA")) >= 1:
                    if type != "":
                        type += " lncRNA:" + str(category_list.count("lncRNA"))
                    else:
                        type += "lncRNA:" + str(category_list.count("lncRNA"))
                if int(category_list.count("miRNA")) >= 1:
                    if type != "":
                        type += " miRNA:" + str(category_list.count("miRNA"))
                    else:
                        type += "miRNA:" + str(category_list.count("miRNA"))
                if int(category_list.count("circRNA")) >= 1:
                    if type != "":
                        type += " circRNA:" + str(category_list.count("circRNA"))
                    else:
                        type += "circRNA:" + str(category_list.count("circRNA"))
                # prepare main table info
                desc = "batch create genesets of table {}".format(str(diff_name))
                if len(seq_list) == 0:
                    continue
                main_info = dict(
                    name=name,
                    task_id=task_id,
                    length=len(seq_list),
                    project_sn=project_sn,
                    level=geneset_type,
                    type=type,
                    desc=desc,
                    status="start",
                    source="diff_exp",
                    is_use=0,
                    group_id=group_id,
                    created_ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                )
                self.bind_object.logger.info(main_info)
                main_id = self.create_db_table('geneset', [main_info])
                # prepare detail table info
                row_dict_list = [{"seq_list": seq_list, "regulate_list": regulate_list, "kind_list": kind_list, "category_list": category_list}]
                # insert detail
                tag_dict = dict(geneset_id=main_id)
                self.create_db_table('geneset_detail', row_dict_list, tag_dict=tag_dict)
                self.update_db_record('geneset', main_id, status="end", main_id=main_id)
                # # delete origin genesets
                # query_info = dict(name=name, task_id=task_id, type=geneset_type)
                # result_ids = self.find_table_by_query_record('sg_geneset', query_info)
                # for result_id in result_ids:
                #     if ObjectId(result_id) != ObjectId(main_id):
                #         print(result_id)
                #         self.remove_db_record('sg_geneset', main_id=ObjectId(result_id))
                #         self.remove_db_record('sg_geneset_detail', geneset_id=ObjectId(result_id))