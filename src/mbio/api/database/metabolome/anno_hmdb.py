# -*- coding: utf-8 -*-
from biocluster.config import Config
from biocluster.api.database.base import Base,report_check
from bson import ObjectId
import json,os
import types
import datetime
from bson.son import SON
import datetime
import pandas as pd

class AnnoHmdb(Base):
    def __init__(self, bind_object):
        super(AnnoHmdb, self).__init__(bind_object)
        self._project_type = 'metabolome'

    def add_anno_hmdb_main(self, name, params, anno_hmdb=None, database="annohmdb",table_type=None,metab_set=None):
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "desc": "",
            "name": name,
            "created_ts": created_ts,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "status": "end",
            "anno_hmdb": anno_hmdb
        }
        if table_type:
            insert_data['table_type'] = table_type
        if metab_set:
            insert_data['metab_set'] = metab_set

        if database == "annohmdb":
            collection = self.db['anno_hmdb']
        elif database == "metabsethmdb":
            collection = self.db['metabset_hmdb']
        else:
            self.bind_object.set_error("databse name is error")
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    def add_level_detail(self, main_id, level_path, database):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id must be ObjectID type or corresponding string！')
        data_list = list()
        pdata = pd.read_table(level_path,header=0,sep='\t')

        #with open(level_path, "r") as f1:
        #    file = f1.readlines()
        var_head = []
        cur_var_head = []
        mk_map = {'ID':'o_id'}
        if 'ID' in pdata.columns:
            cur_var_head.append('ID')
            var_head.append(mk_map['ID'])
        for i in range(len(pdata)):
            # for line in file[1:]:
            #     line = line.strip().split("\t")
            data = [
                ('hmdb_id', main_id),
                ('metab', pdata['Metabolite'][i]),
                ('metab_id', pdata['Metab_id'][i]),
                ('kingdom', pdata['Kingdom'][i]),
                ('superclass', pdata['Superclass'][i]),
                ('class', pdata['Class'][i]),
                ('subclass', pdata['Subclass'][i])
            ]
            for k in cur_var_head:
                data.append((mk_map[k],pdata[k][i]))
            data = SON(data)
            data_list.append(data)
        if database == "annohmdb":
            level_collection = self.db['anno_hmdb_level']
            main_collection = self.db['anno_hmdb']

        elif database == "metabsethmdb":
            level_collection = self.db['metabset_hmdb_level']
            main_collection = self.db['metabset_hmdb']
        else:
            self.bind_object.set_error("database name {} is error".format(database))
        try:
            if data_list:
                level_collection.insert_many(data_list)
                if var_head != []:
                    main_collection.update({"_id":main_id},{"$set":{"var_head":','.join(var_head)}})  #zouguanqing 20190617
            else:
                self.bind_object.logger.info("level表为空")
        except Exception, e:
            self.bind_object.set_error("%s import error：%s" , variables=(level_path, e))
        else:
            self.bind_object.logger.info("import level succeed")

    def add_stat_detail(self, main_id, stat_dir, database):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id must be ObjectID type or corresponding string")
        data_list = list()
        file_dict = {"superclass": "HmdbSuperclass.xls","class": "HmdbClass.xls","subclass":"HmdbSubclass.xls"}
        for each in file_dict.keys():
            mytype = each
            file = os.path.join(stat_dir, file_dict[each])
            if not os.path.exists(file):
                self.bind_object.set_error("% file not exists".format(file))
            with open(file, "r") as f:
                table = f.readlines()
                for line in table[1:]:
                    line = line.strip().split("\t")
                    data = [
                        ("hmdb_id", main_id),
                        ("name", line[0]),
                        ("number", int(line[1])),
                        ("metab_id", line[2]),
                        ("type", mytype),
                    ]
                    data = SON(data)
                    data_list.append(data)
        if database == "annohmdb":
            stat_coll = self.db['anno_hmdb_stat']
        elif database == "metabsethmdb":
            stat_coll = self.db['metabset_hmdb_stat']
        else:
            self.bind_object.set_error("databse name is error")
        try:
            if data_list:
                stat_coll.insert_many(data_list)
            else:
                self.bind_object.logger.info("stat table is empty")
        except Exception, e:
            self.bind_object.set_error("%s import error：%s" , variables=(level_path, e))
        else:
            self.bind_object.logger.info("import stat succeed")

