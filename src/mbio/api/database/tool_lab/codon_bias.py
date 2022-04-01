# -*- coding: utf-8 -*-



import os
from bson.objectid import ObjectId
import pandas as pd
import glob
from api_base import ApiBase
import json

class CodonBias(ApiBase):
    def __init__(self, bind_object):
        super(CodonBias, self).__init__(bind_object)

    def add_detail(self,infile,main_id):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)

        insert_data = []
        with open(infile) as fr:
            fr.readline()
            for line in fr:
                spline = line.strip().split('\t')
                tmp = {
                    'aa' :spline[0],
                    'codon' : spline[1].replace('NAN',''),
                    'num' : float(spline[2]),
                    'codon_id' : main_id
                }
                insert_data.append(tmp)
        try:
            self.db['codon_bias_detail'].insert_many(insert_data)
            tmp = {"column":[{"filter": "false", "field": "aa", "type": "string", "sort": "false", "title": "aa"},
                             {"filter": "false", "field": "num", "type": "string", "sort":"false", "title": "Sample"},
                             {"filter": "false", "field": "codon", "type": "string", "sort": "false", "title": "Codon"}],
                   "condition": {}}
            updata_info= {"data_table": json.dumps(tmp) }
            self.db['codon_bias'].update({"_id":main_id},{"$set":updata_info})
        except Exception as e:
            self.bind_object.logger.info("导入codon_bias_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入codon_bias_detail数据成功")

    def add_stats(self, infile,main_id,pic_path):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        insert_data = []
        with open(infile) as fr:
            fr.readline()
            for line in fr:
                spline = line.strip().split('\t')
                tmp = {
                    'codon_id' : main_id,
                    'title' : spline[0].strip(),
                    't3s': float(spline[1]),
                    'c3s': float(spline[2]),
                    'a3s': float(spline[3]),
                    'g3s': float(spline[4]),
                    'cai': float(spline[5]),
                    'cbi': float(spline[6]),
                    'fop': float(spline[7]),
                    'nc': float(spline[8]),
                    'gc3s': float(spline[9]),
                    'gc': float(spline[10]),
                    'l_sym': float(spline[11]),
                    'l_aa': float(spline[12]),
                    'gravy': float(spline[13]),
                    'aromo': float(spline[14]),
                }
                insert_data.append(tmp)
        try:
            self.db['codon_bias_stats'].insert_many(insert_data)
            tmp = {"column":[{"filter": "false", "field": "a3s", "type": "string", "sort": "false", "title": "A3S"},
                             {"filter": "false", "field": "gravy", "type": "string", "sort": "false", "title": "Gravy"},
                             {"filter": "false", "field": "cbi", "type": "string", "sort": "false", "title": "CBI"},
                             {"filter": "false", "field": "gc", "type": "string", "sort": "false", "title": "GC"},
                             {"filter": "false", "field": "fop", "type": "string", "sort": "false", "title": "FOP"},
                             {"filter": "false", "field": "g3s", "type": "string", "sort": "false", "title": "G3S"},
                             {"filter": "false", "field": "gc3s", "type": "string", "sort": "false", "title": "GC3S"},
                             {"filter": "false", "field": "title", "type": "string", "sort": "false", "title": "title"},
                             {"filter": "false", "field": "nc", "type": "string", "sort": "false", "title": "NC"},
                             {"filter": "false", "field": "cai", "type": "string", "sort": "false", "title": "CAI"},
                             {"filter": "false", "field": "l_sym", "type": "string", "sort": "false", "title": "L_sym"},
                             {"filter": "false", "field": "t3s", "type": "string", "sort": "false", "title": "T3s"},
                             {"filter": "false", "field": "c3s", "type": "string", "sort": "false", "title": "C3S"},
                             {"filter": "false", "field": "l_aa", "type": "string", "sort": "false", "title": "L_aa"}],
                   "condition": {}}

            self.db['codon_bias'].update({"_id":main_id},{"$set":{"data_table2": json.dumps(tmp), "pic_path": pic_path}})

        except Exception as e:
            self.bind_object.logger.info("导入codon_bias_stats数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入codon_bias_stats数据成功")




