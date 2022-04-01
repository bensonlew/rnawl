# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modify:20180605
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
import math
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metabolome.common import check_metab
import pandas as pd


class Enrich(Base):
    # 对应表格链接：http://git.majorbio.com/liu.linmeng/metabolome/wikis/collection/metab_set/ipath
    def __init__(self, bind_object):
        super(Enrich, self).__init__(bind_object)
        self._project_type = "metabolome"
        self.metab_map = {}  # map metab_id to metab name

    @report_check
    def add_enrich(self, name, params,metab_set=None):
        """
        工作流对metab_table主表进行导表
        :param name: 工作流需要导两张主表，第一张对应Raw_，为原始数据表；第二张对应Table_，为预处理后的结果
        :param params: 工作流选择的运行参数
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'MetabsetEnrich_Origin',
            'name': name,
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'created_ts': created_ts,
            'status': 'end',
            'has_topo' : 1
        }
        if metab_set:
            if not isinstance(metab_set,ObjectId):
                metab_set = ObjectId(metab_set)
            insert_data['metab_set'] = metab_set

        collection = self.db['metabset_kegg_enrich']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    @report_check
    def add_enrich_detail(self, main_id, overview_table, enrich_table):
        """
        :param main_id: 主表id
        :param enrich_table: 富集分析结果表
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54700601")
        if not os.path.exists(enrich_table):
            self.bind_object.set_error("enrich_table的路径不存在，检查：%s" , variables=(enrich_table), code="54700602")
        self.metab_init(overview_table)
        data_list = list()
        kegg_type1 = list()
        with open(enrich_table, "r") as f1:
            lines = f1.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    'kegg_enrich_id': main_id,
                    'study_count': int(line[0]),
                    'description':line[1],
                    'database': line[2],
                    'pathway_id': line[3].split("path:")[1] if "path:" in line[3] else line[3],
                    'ratio_in_study': line[4],
                    'ratio_in_pop': line[5],
                    'bg_number': line[5].split("/")[1],
                    'pvalue': round(float(line[6]), 4),
                    'qvalue': line[7] if line[7] == "None" else round(float(line[7]), 4),
                    'metab_id': line[8],
                    'metab_list': self.get_metab_list(line[8]),
                    'hyperlink': line[9],
                    'kegg_type': "".join([x[0] for x in line[11].split(' ')]),
                    'log_p': self.neg_log(line[6]),
                    'log_q': self.neg_log(line[7]),
                }
                try:
                    data['enrich_factor'] = float(line[4].split("/")[0])/float(line[5].split("/")[0])
                except:
                    data['enrich_factor'] = 0.0
                kegg_type1.append("".join([x[0] for x in line[11].split(' ')]))
                data_list.append(data)
        try:
            if data_list:
                detail_coll = self.db["metabset_kegg_enrich_detail"]
                detail_coll.insert_many(data_list)
                main_col = self.db['metabset_kegg_enrich']
                kegg_type1 = list(set(kegg_type1))
                self.bind_object.logger.info("main_id: %s\ncategories: %s" % (main_id, kegg_type1))
                main_col.update({'_id':main_id}, {'$set': {'categories': kegg_type1}})
            else:
                self.bind_object.logger.info("no result with the metab")
        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(overview_table,e), code="54700603")
        else:
            self.bind_object.logger.info("导入kegg_enrich_detail表成功")

    def add_topology(self,main_id, topo_table,remote_path):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        if not os.path.exists(topo_table):
            self.bind_object.logger.info("%s 不存在。不导metabset_kegg_enrich_topo"%topo_table)
            return
        my_data = pd.read_table(topo_table,sep='\t',header=0)
        insert_list = []
        insert_list_has_pic = []
        dir_path = os.path.dirname(topo_table)

        for i in range(len(my_data)):
            insert_data = {
                "kegg_enrich_id":main_id,
                "pathway_id": my_data['ko'][i],
                "description": my_data['pathway_name'][i],
                "match_status": my_data['match_status'][i],
                "pvalue": my_data['p_value'][i],
                "d_pvalue": my_data['FDR'][i],
                "impact": my_data['impact_value'][i],
                "metab_ids": my_data['metab_ids'][i],
                "compound_ids" : my_data['compound_ids'][i]
            }
            if os.path.exists(dir_path + '/'+my_data['ko'][i]+'.network.png') :
                insert_data['pic'] = os.path.join(remote_path,my_data['ko'][i]+'.network.png')
                insert_list_has_pic.append(insert_data)
            else:
                insert_list.append(insert_data)
        try:
            collection = self.db['metabset_kegg_enrich_topo']
            collection.insert_many(insert_list_has_pic+insert_list)  ## 有图片的先导库
        except Exception:
            self.bind_object.set_error('metabset_kegg_enrich_topo 导表错误： %s '%topo_table)
        else:
            self.bind_object.logger.info('metabset_kegg_enrich_topo 导表成功')



    def metab_init(self, overview_table):
        """
        通过查询数据库，获取到metab_id与代谢物名称的对应关系
        :return: self.metab_map
        """
        f1 = open(overview_table).readlines()
        data_dic = {}
        for line in f1:
            line = line.strip().split("\t")
            data_dic[line[0]] = line[1]
        self.metab_map = data_dic

    def get_metab_list(self, metab_id_str):
        """
        根据metab_id，生成对应的代谢物名称
        :param metab_id_str:
        :return:
        """
        if len(self.metab_map) == 0:
            self.bind_object.set_error("self.metab_map is none", code="54700604")
        metab_id_list = metab_id_str.split("|")
        metab_name_list = []
        for metab_id in metab_id_list:
            metab_name = self.metab_map[metab_id]
            metab_name = check_metab(metab_name)
            metab_name_list.append(metab_name)
        metab_name_str = '|'.join(metab_name_list)
        return metab_name_str

    def neg_log(self, value):
        value = float(value)
        if value > 0 :
            log10value = -math.log10(value)
        else:
            log10value = -math.log10(0.0001)
        return  log10value



