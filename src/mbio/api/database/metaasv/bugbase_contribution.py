# -*- coding: utf-8 -*-
# __author__ = zzg

from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name
import json


class BugbaseContribution(Base):
    def __init__(self, bind_object):
        super(BugbaseContribution, self).__init__(bind_object)
        self._project_type = "metaasv"

    #@report_check
    def add_contribute(self, geneset_id, fun_tax_file, tax_species_old, main=True, main_table_id=None, name=None,
                       params=None):
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52800501")
        if main:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            if params != None:
                params = params
                print params
            else:
                params = ""
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "Contribute_Origin",
                'params': params,
                'status': 'end',
            }
            try:
                collection = self.db['bugbase_contribution ']
                contribute_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入bugbase_contribution 主表异常:{}'.format(e))
                self.bind_object.set_error("导入bugbase_contribution 主表异常")
            print contribute_id
        return contribute_id

    #@report_check
    def add_contribute_detail(self, contribute_id, fun_tax_file,tax_species_old, project_sn='bugbase',task_id='bugbase',params='bugbase',update_main=True):
        collection = self.db["bugbase_contribution"]
        if (contribute_id is None) or (contribute_id == ""):
            name = "Bugbase" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Bugbase',
                params=params if params else "null",
                status="start",
                group="C,N,S"
            )
            # main_id = self.create_db_table('sg_bugbase', [main_info])
            main_id = collection.insert_one(main_info).inserted_id
        else:
            main_id = ObjectId(contribute_id)

        data_list = []
        print main_id
        self.bind_object.logger.error("开始导表")
        self.bind_object.logger.error("开始fun_tax_file")
        tax_file = {}
        with open(tax_species_old, 'rb') as f2:
            data = f2.readlines()
            for i in data[1:]:
                tax_file[i.split("\t")[0]] = i.strip().split("\t")[-1]
        fun_list = []
        all_phenotypes = []
        with open(fun_tax_file, 'rb') as f1:
            data = f1.readlines()
            for i in data[1:]:
                if i.strip().split("\t")[0] not in all_phenotypes:
                    all_phenotypes.append(i.strip().split("\t")[0])
                total = 0
                for x in range(0, len(data[0].strip().split("\t")) - 3):
                    total += float(i.strip().split("\t")[x + 2])
                if total != 0:
                    insert_data = {
                        'bugbase_id': main_id,
                        'taxon': i.strip().split("\t")[1],
                        'phenotypes': i.strip().split("\t")[0],
                        'tax_abundance': float(tax_file[i.strip().split("\t")[1]]) if i.strip().split("\t")[
                                                                                          1] in tax_file else 0
                    }
                    for x in range(0, len(data[0].strip().split("\t")) - 3):
                        insert_data[data[0].strip().split("\t")[x + 2]] = i.strip().split("\t")[x + 2]
                    total = 0
                    data_list.append(insert_data)
        all_name = ["phenotypes", "taxon"]
        for i in data[0].strip().split("\t")[2:-1]:
            all_name.append(i)
        table_dict = {"table_data": all_name}
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
        column_dict = {"column_data": {"name":"phenotypes","data":data[0].strip().split("\t")[2:-1],"category": ""}}
        column_info = json.dumps(column_dict, sort_keys=True, separators=(',', ':'))

        try:
            collection1 = self.db['bugbase_contribution_detail']
            collection1.insert_many(data_list)
            collection.update_one({'_id': main_id}, {"$set": {"main_id": main_id, "table_data": table_info, "column_data": column_info,"all_function":",".join(all_phenotypes)}})
        except Exception as e:
            self.bind_object.logger.error("导入bugbase_contribution_detail表格信息出错:%s" % (e))
            self.bind_object.set_error("导入bugbase_contribution_detail表信息出错", code="52800506")
        else:
            self.bind_object.logger.info("导入bugbase_contribution_detail表格信息成功!")
