# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'


import datetime
import os
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from bson.objectid import ObjectId


class GenesetBatch(ApiBase):
    def __init__(self, bind_object):
        super(GenesetBatch, self).__init__(bind_object)
        self._project_type = 'denovo_rna_v2'

    def add_diff_genesets(self, task_id, geneset_file, geneset_type, geneset_suffix, diff_id):
        if not os.path.exists(geneset_file):
            raise Exception("文件不存在{}".format(geneset_file))
        main_info = self.db["sg_diff"].find_one({"main_id": ObjectId(diff_id)})
        project_sn = main_info["project_sn"]
        diff_name = main_info["name"]
        group_id = main_info["group_id"]
        with open(geneset_file, "r") as f:
            head = f.readline()
            for line in f:
                items = line.strip().split("\t")
                if len(items) <= 2:
                    continue
                if not (items[0] or items[1]):
                    continue
                name = items[0]
                name = name + "_" + geneset_type + geneset_suffix
                seq_list = items[1].split(",")
                regulate_list = items[2].split(",")
                # prepare main table info
                desc = "batch create genesets of table {}".format(str(diff_name))
                if len(seq_list) == 0:
                    continue
                main_info = dict(
                    name=name,
                    task_id=task_id,
                    gene_length=len(seq_list),
                    project_sn=project_sn,
                    type=geneset_type,
                    desc=desc,
                    status="start",
                    source="diff_exp",
                    is_use=0,
                    group_id=group_id,
                    created_ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                )
                self.bind_object.logger.info(main_info)
                main_id = self.create_db_table('sg_geneset', [main_info])
                # prepare detail table info
                row_dict_list = [{"seq_list": seq_list, "regulate_list": regulate_list}]
                # insert detail
                tag_dict = dict(geneset_id=main_id)
                self.create_db_table('sg_geneset_detail', row_dict_list, tag_dict=tag_dict)
                self.update_db_record('sg_geneset', main_id, status="end", main_id=main_id)
                # # delete origin genesets
                # query_info = dict(name=name, task_id=task_id, type=geneset_type)
                # result_ids = self.find_table_by_query_record('sg_geneset', query_info)
                # for result_id in result_ids:
                #     if ObjectId(result_id) != ObjectId(main_id):
                #         print(result_id)
                #         self.remove_db_record('sg_geneset', main_id=ObjectId(result_id))
                #         self.remove_db_record('sg_geneset_detail', geneset_id=ObjectId(result_id))