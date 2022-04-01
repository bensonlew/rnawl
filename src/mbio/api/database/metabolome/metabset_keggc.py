# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base
from bson import ObjectId
import json
import types
import datetime
import pandas as pd
from bson.son import SON

class MetabsetKeggc(Base):
    def __init__(self, bind_object):
        super(MetabsetKeggc, self).__init__(bind_object)
        self._project_type = "metabolome"

    def add_metabsetc(self, params, name="KEGGC_Origin"):
        # params{set_id}用来保存代谢集主表记录
        # project_sn, task_id, name, status, desc, created_ts,
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "desc": "代谢集KEGG化合物分类",
            "created_ts": created_ts,
            "status": "end",
            "name": name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':'))
        }
        collection = self.db['metabset_keggc']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    def add_metabsetc_detail(self, main_id, table):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54701301")
        data_list = list()
        in_data = pd.read_table(table,sep="\t",header=0)
        first_key = 'fist_category'
        if first_key not in in_data.columns:
            first_key = 'first_category'
        for i in range(len(in_data)):
            data = [
                ("kegg_id", main_id),
                ("first_category", in_data[first_key][i]),
                ("second_category", in_data['second_category'][i]),
                ("metab_id", in_data['metab_id'][i]),
                ("count", in_data['count'][i])
            ]
            data = SON(data)
            data_list.append(data)
        detail_coll = self.db["metabset_keggc_detail"]
        try:
            if data_list:
                detail_coll.insert_many(data_list)
            main_table = self.db['metabset_keggc']
            main_data = main_table.find_one({"_id":main_id})
            if main_data:
                params = main_data['params']
                metab_set_id = eval(params)['metabset']
                main_table.update({"_id":main_id},{"$set":{"metab_set":ObjectId(metab_set_id)}})
        except Exception,e:
            self.bind_object.set_error("导入%s信息出错： %s" , variables=(table,e), code="54701302")
        else:
            self.bind_object.logger.info("导入stat信息成功")
