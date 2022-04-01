# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from mbio.api.database.tool_lab.api_base import ApiBase
from biocluster.api.database.base import report_check
from bson.objectid import ObjectId
from bson.son import SON
import os
import math


class Enrichment(ApiBase):
    def __init__(self, bind_object):
        super(Enrichment, self).__init__(bind_object)

    @report_check
    def add_enrich_detail(self, result, main_id):
        insert_data = list()
        if not os.path.exists(result):
            self.bind_object.set_error("结果文件不存在")
        with open(result, "r") as f:
            head = f.readline()
            for line in f:
                items = line.strip().split("\t")
                gene_ids = list(items[7].split("/"))
                data = [
                    ("enrich_id", ObjectId(main_id)),
                    ("id", items[0]),
                    ("description", items[1]),
                    ("gene_ratio", items[2]),
                    ("gene_count", int(items[2].split("/")[0])),
                    ("bg_ratio", items[3]),
                    ("bg_count", int(items[3].split("/")[0])),
                    ("pvalue", float(items[4])),
                    ("padjust", float(items[5])),
                    ("qvalue", float(items[6])),
                    ("gene_ids", gene_ids),
                    ("log10padjust", -(math.log10(float(items[5])))),
                    ("rich_factor", round(float(items[2].split("/")[0]) / float(items[3].split("/")[0]), 4)),
                ]
                insert_data.append(SON(data))
        try:
            collection = self.db['enrichment_detail']
            collection.insert_many(insert_data)
        except Exception, e:
            self.bind_object.set_error("导入富集分析结果失败, {}".format(e))
        else:
            self.bind_object.logger.info("导入富集分析结果成功")
        try:
            main_table = self.db['enrichment']
            table_data = "{\"column\":[{\"field\":\"id\", \"title\":\"ID\", \"filter\":\"false\", \"sort\":\"false\", \"type\": \"string\"}, {\"field\":\"description\", \"title\":\"Description\", \"filter\":\"false\", \"sort\":\"false\", \"type\": \"string\"}, {\"field\":\"gene_ratio\", \"title\":\"GeneRatio\", \"filter\":\"false\", \"sort\":\"false\", \"type\": \"string\"}, {\"field\":\"bg_ratio\", \"title\":\"BgRatio\", \"filter\":\"false\", \"sort\":\"false\", \"type\": \"string\"}, {\"field\":\"pvalue\", \"title\":\"P_value\", \"filter\":\"true\", \"sort\":\"true\", \"type\": \"float\"}, {\"field\":\"padjust\", \"title\":\"P_adjust\", \"filter\":\"true\", \"sort\":\"true\", \"type\": \"float\"}, {\"field\":\"qvalue\", \"title\":\"Q_value\", \"filter\":\"true\", \"sort\":\"true\", \"type\": \"float\"}, {\"field\":\"gene_ids\", \"title\":\"GeneIDs\", \"filter\":\"false\", \"sort\":\"false\", \"type\": \"string\"}, {\"field\":\"gene_count\", \"title\":\"Count\", \"filter\":\"false\", \"sort\":\"false\", \"type\": \"int\"}]}"
            venn_upset_data = "{\"name\": \"description\", \"data\": \"gene_ids\"}"
            bubble_data = "{\"name\": \"description\", \"x\": \"rich_factor\", \"y\": \"description\", \"size\": \"gene_count\", \"color\": \"padjust\", \"category\": \"category\"}"
            line_data = "{\"name\": \"name\", \"y\": \"gene_count\", \"x\": \"description\"}"
            column_data = "{\"name\": \"description\", \"data\": \"log10padjust\"}"
            main_table.update({'_id': ObjectId(main_id)}, {
                '$set': {'table_data': table_data, 'venn_upset_data': venn_upset_data, 'bubble_data': bubble_data,
                         'line_data': line_data, 'column_data': column_data}}, upsert=True)
        except Exception, e:
            self.bind_object.set_error("更新富集分析主表失败, {}".format(e))
        else:
            self.bind_object.logger.info("更新富集分析主表成功")
