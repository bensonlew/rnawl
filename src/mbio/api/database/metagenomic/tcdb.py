# -*- coding:utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import datetime
import types
from mbio.packages.metagenomic.id_convert import name2id
import json

class Tcdb(Base):
    def __init__(self, bind_object):
        super(Tcdb, self).__init__(bind_object)
        self._object_type = 'metagenomic'

    @report_check
    def add_tcdb(self, geneset_id, specimen, anno_file, name=None, params=None):
        geneset_id_str = ''
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                geneset_id = ObjectId(geneset_id)
            else:
                self.bind_object.set_error('geneset_id必须为ObjectId对象或者其对应的字符串！', code="52804801")

        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "正在计算",
            "name": name if name else 'AnnoTcdb_Origin',
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "anno_file" : anno_file,
            "geneset_id" : geneset_id,
            "specimen": specimen,
            "is_origin" : 1,
            "lowest_level" : "TCDB ID",
            "settled_params": json.dumps({"version": "tcdb_v20200917"})  # by zhaozhigang 20200929
        }

        main = self.db['anno_tcdb']
        main_id = main.insert_one(data).inserted_id
        return main_id

    @report_check
    def add_tcdb_detail(self, path, tcdb_id):
        if not isinstance(tcdb_id, ObjectId):
            tcdb_id = ObjectId(tcdb_id)
        with open(path) as f:
            f.readline()
            insert_list = []
            put_list = []
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                if len(spline) < 6:
                    continue
                if spline[1] in put_list:
                    continue
                put_list.append(spline[1])
                data = {
                    'tcdb_id': tcdb_id,
                    't_name' : spline[1],
                    'desc' : spline[2].replace('/','&'),
                    'family' : spline[3],
                    'subclass' : spline[4],
                    'class' : spline[5],
                }
                insert_list.append(data)
        try:
            collection = self.db['anno_tcdb_detail']
            collection.insert_many(insert_list)
            self.db['anno_tcdb'].update_one({'_id': tcdb_id}, {'$set':{'main_id':tcdb_id,"settled_params": json.dumps({"version": "tcdb_v20200917"})}})
        except Exception as e:
            self.bind_object.logger.info('导入anno_tcdb_detail 出错 %s' % e )
            self.bind_object.set_error('import anno_tcdb_detail ERROR', code="52804802")
        else:
            self.bind_object.logger.info('导入 anno_tcdb_detail 成功')


    @report_check
    def add_tcdb_abund(self, path, tcdb_id, level):
        if not isinstance(tcdb_id, ObjectId):
            tcdb_id = ObjectId(tcdb_id)
        result = self.db['anno_tcdb'].find_one({'_id': tcdb_id})
        task_id = result['task_id']
        samples_dic = name2id(task_id, type="task")
        with open(path) as f:
            head = f.readline().strip()
            sp_head = head.split('\t')
            insert_list = []
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                if len(spline) < 2:
                    continue
                data= {
                    'tcdb_id' : tcdb_id,
                    'level' : level,
                    'name' : spline[0],
                    'desc' : spline[-1].replace('/','&')
                }

                for i,e in enumerate(sp_head[1:-1],1):
                    if e != 'Total':
                        e = samples_dic[e]
                    data[e] = int(spline[i])
                insert_list.append(data)

            try:
                collection = self.db['anno_tcdb_abund']
                collection.insert_many(insert_list)
            except Exception as e:
                self.bind_object.logger.info('导入 anno_tcdb_abund 失败')
                self.bind_object.set_error('import anno_tcdb_abund ERROR', code="52804803")
            else:
                self.bind_object.logger.info('导入 anno_tcdb_abund 成功')









