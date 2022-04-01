# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last_modify:20180605
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
from contextlib import nested
from mbio.packages.metabolome.common import check_metab


class Ipath(Base):
    # 对应表格链接：http://git.majorbio.com/liu.linmeng/metabolome/wikis/collection/metab_set/ipath
    def __init__(self, bind_object):
        super(Ipath, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_ipath(self, name, rere_path, params=None):
        """
        工作流对metab_table主表进行导表
        :param name: 工作流需要导两张主表，第一张对应Raw_，为原始数据表；第二张对应Table_，为预处理后的结果
        :param type: 流程的类型, {GC|LC}
        :param params: 工作流选择的运行参数
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'MetabsetKeggIpath_Origin',
            'name': name,
            'created_ts': created_ts,
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')) if params else "",
            'status': 'end',
            "pathways": ["Metabolism.svg", "Antibiotics.svg", "Microbial_metabolism.svg","Secondary_metabolites.svg"] ,   #["Metabolic_pathways.svg", "Regulatory_pathways.svg","Biosynthesis_of_secondary_metabolities.svg"],
            "result_dir": rere_path
        }

        collection = self.db['metabset_ipath']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    @report_check
    def add_ipath_detail(self, main_id, result_dir, metab_set=None):
        """
        对阴离子，阳离子和合成表分别进行导表
        :param main_id: 主表id
        :param table_path: 三种类型数据表各自的路径，里面对应存储代谢物明细表和丰度表，逗号分割
        :param metab_set: 代谢集表格
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54700701")
        main_col = self.db['metabset_ipath']
        insert_file = os.path.join(result_dir, 'gene_ipath_input.xls')
        data_list = list()
        with open(insert_file, "r") as f1:
            lines = f1.readlines()
            for line in lines:
                line = line.strip().split("\t")
                data = [
                    ("ipath_id", main_id),
                    ("metabolite", check_metab(line[0])),
                    ("compound_id", line[1]),
                    ("color", line[2]),
                    ("width", line[3])
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["metabset_ipath_detail"]
        try:
            if len(data_list) != 0:
                detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(insert_file,e), code="54700702")
        else:
            self.bind_object.logger.info("导入stat信息成功")
        if metab_set:
            doc_keys = set()
            with open(metab_set,"r") as f2:
                for l in f2.readlines():
                    name = l.strip().split('\t')[0]
                    doc_keys.add(name)
            metab_set_name = list(doc_keys)
            # main_col.update_one({'_id': main_id}, {'$set': {'result_dir': result_dir, 'table_columns': metab_set_name}})
        main_col.update_one({'_id': main_id}, {'$set': {'table_columns': metab_set_name}})

    def add_relation_ipath(self, main_id, set_to_ipath):
        gene_set = []
        metab_set = []
        pathways = ["Antibiotics.svg", "Metabolism.svg", "Microbial_metabolism.svg", "Secondary_metabolites.svg"]
        with open(set_to_ipath, 'r') as r:
             for l in r:
                 l = l.strip().split('\t')
                 if len(l) < 4:
                     continue
                 if l[0] == "metab_set":
                     metab_set.append(l[3].split(';'))
                 else:
                     gene_set.append(l[3].split(';'))
        main_col = self.db['relation_ipath']
        main_id = ObjectId(main_id)
        main_col.update_one({'_id': main_id}, {'$set': {'pathways': pathways, 'metab_list': metab_set, 'gene_list': gene_set}})

    def add_relation_detail(self, main_id, result_file):
        detail_col = self.db['relation_ipath_detail']
        data_list = []
        with open(result_file, 'r') as r:
            for l in r :
                l = l.strip('\n').split('\t')
                data = {
                    "ipath_id": ObjectId(main_id),
                    "item": l[0],
                    "annot": l[1],
                    "color": l[2],
                    "width": l[3],
                }
                data_list.append(data)
        try:
            detail_col.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入{}信息出错： {}".format(result_file, e))
        else:
            self.bind_object.logger.info("导入stat信息成功")
