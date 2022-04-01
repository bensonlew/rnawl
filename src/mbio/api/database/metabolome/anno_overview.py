# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base,report_check
from bson import ObjectId
import json
import types
import datetime
from bson.son import SON
import pandas as pd
import re


class AnnoOverview(Base):
    def __init__(self, bind_object):
        super(AnnoOverview, self).__init__(bind_object)
        self._project_type = 'metabolome'

    """
    ##################################################
    EXPORT TEST
    ##################################################

    def export_test(self, table_path):
        params = {
            "set_id": "5b0e4349a4e1af498abb645a"
        }
        main = self.add_overview(params=params)
        self.add_overview_detail(main, table_path)
        self.add_overview_ko(main, table_path + "_ko.xls")
    """

    @report_check
    def add_overview(self, params, name="Overview_Origin"):
        # params{set_id}用来保存代谢集主表记录
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn" : project_sn,
            "task_id": task_id,
            "desc": "代谢物总览",
            "created_ts": created_ts,
            "status": "end",
            "name": name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "version" : "2.0"
        }
        collection = self.db['anno_overview']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    @report_check
    def add_overview_detail(self, main_id, overview_table,origin_exp_file=None):
        if origin_exp_file:  #用来标记没有表达量的代谢物
            exp_data = pd.read_table(origin_exp_file,sep='\t',header=0)
            has_exp_list = exp_data['metab_id'].tolist()

        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54700201")
        data_list = list()
        mdata = pd.read_table(overview_table,header=0,sep='\t')
        #with open(overview_table, "r") as f1:
        #    file = f1.readlines()
        #    for line in file[1:]:
        #        line = line.strip().split("\t")
        add_all_heads = ['super_class','class','sub_class','kegg_id_desc', 'ID', 'm/z', 'Retention time','RT (min)', 'Adducts','hmdb_id',
                         'cas_id','formula','hmdb_name','element','mass','Similarity', "Theoretical Fragmentation Score", "Fragmentation Score", "Score"]
        mongo_k = ['superclass','class','subclass','kegg_id_desc', 'o_id', 'mz', 'rt','rt', 'adducts','hmdb_id',
                   'cas_id','formula','hmdb_name','element','mass','similarity', 'theo_frag', 'frag', 'score']  # 去'hmdb_name','element','mass'
        map_mk = dict(zip(add_all_heads,mongo_k))

        c_add = []
        c_add_mk = []
        for k in add_all_heads:
            if k in mdata.columns:
                c_add.append(k)
                c_add_mk.append(map_mk[k])


        for i in range(len(mdata)):

            data = [
                ("overview_id", main_id),
                ("metab_id", mdata['metab_id'][i]),
                ("metab", mdata['metab'][i]),
                ("mode", mdata['mode'][i]),
                #("hmdb_id", mdata['hmdb_id'][i]),
                #("hmdb_name", mdata['hmdb_name'][i]),
                #("mass", mdata['mass'][i]),
                #("cas_id", mdata['cas_id'][i]),
                #("formula", mdata['formula'][i]),
                #("element", mdata['element'][i]),
                ("compound_id", mdata['compound_id'][i]),
                ("compound_name", mdata['compound_name'][i]),
                ("c1_category", mdata['compound_first_category'][i]),
                ("c2_category", mdata['compound_second_category'][i]),
                ("pathway_id", mdata['pathway_id'][i]),
                ("description", mdata['description'][i]),
                ("p1_category", mdata['kegg_first_category'][i]),
                ("p2_category", mdata['kegg_second_category'][i])
            ]
            if re.match('metab_\d*$',mdata['metab'][i]):    # 标记没有代谢物名称的数据
                data.append(('no_name',1))
            if origin_exp_file:
                if mdata['metab_id'][i] not in has_exp_list: # 标记没有表达量的代谢物
                    data.append(('no_exp',1))
                else:
                    data.append(('has_exp',1))

            for k in c_add:  #zouguanqing 20190617
                data.append((map_mk[k],mdata[k][i]))
            # if 'super_class' in mdata.columns:
            #     data.append(("superclass", mdata['super_class'][i]))
            #     data.append(("class", mdata['class'][i]))
            #     data.append(("subclass", mdata['sub_class'][i]))
            data = SON(data)
            data_list.append(data)
        detail_coll = self.db['anno_overview_detail']
        try:
            detail_coll.insert_many(data_list)
            var_head = 'metab,mode'
            var_head += ',c1_category,c2_category,pathway_id,p1_category,p2_category,description'
            for rm_k in ['hmdb_name','element','mass']:
                if rm_k in c_add_mk:
                    c_add_mk.remove(rm_k)
            var_head += ',' +','.join(c_add_mk)
            main_collection = self.db['anno_overview']
            var_head = ','.join(set(var_head.split(',')))
            main_collection.update({"_id":main_id},{"$set":{"var_head":var_head}})

        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(overview_table,e), code="54700202")
        else:
            self.bind_object.logger.info("导入stat信息成功")

    @report_check
    def add_overview_ko(self, main_id, ko_table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54700203")
        data_list = list()
        with open(ko_table, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("overview_id", main_id),
                    ('pathway_id', line[0]),
                    ("description", line[1]),
                    ("first_category", line[2]),
                    ("second_category", line[3]),
                    ("compound_id", line[4]),
                    ("metab_id", line[5]),
                    ("hyperlink", line[6])
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db['anno_overview_ko']
        try:
            detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(overview_table,e), code="54700204")
        else:
            self.bind_object.logger.info("导入stat信息成功")