# -*- coding:utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import datetime
import types
from mbio.packages.metagenomic.id_convert import name2id


class Mvirdb(Base):
    def __init__(self, bind_object):
        super(Mvirdb, self).__init__(bind_object)
        self._object_type = 'metagenomic'

    @report_check
    def add_mvirdb(self,geneset_id, specimen, anno_file, name=None, params=None):
        geneset_id_str = ''
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                geneset_id = ObjectId(geneset_id)
            else:
                self.bind_object.set_error('geneset_id必须为ObjectId对象或者其对应的字符串！', code="52804101")
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "正在计算",
            "name": name if name else 'AnnoMvirdb_Origin',
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "type" : 1,
            "anno_file": anno_file,
            "geneset_id" : geneset_id,
            "specimen": specimen,
            "is_origin" : 1,
            "lowest_level": "Virulence Factor ID"
        }

        main = self.db['anno_mvir']
        main_id = main.insert_one(data).inserted_id
        return main_id

    @report_check
    def add_mvirdb_detail(self, path, mvir_id):
        if not isinstance(mvir_id, ObjectId):
            mvir_id = ObjectId(mvir_id)
        with open(path) as f:
            f.readline()
            insert_list = []
            put_list = []
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                if len(spline) < 8:
                    continue
                if spline[1] in put_list:
                    continue
                put_list.append(spline[1])
                data = {
                    'mvir_id': mvir_id,
                    'factor_id' : spline[1],
                    'source' : spline[2],
                    'designate' : spline[3],
                    'gene' : spline[4],
                    'vir_type' : spline[6],
                    'desc' : spline[5],
                    'status' : spline[7]
                }
                insert_list.append(data)
        try:
            collection = self.db['anno_mvir_detail']
            collection.insert_many(insert_list)
        except Exception as e:
            self.bind_object.logger.info('导入anno_mvir_detail 出错 %s' % e )
            self.bind_object.set_error('import anno_mvir_detail ERROR', code="52804102")
        else:
            self.bind_object.logger.info('导入 anno_mvir_detail 成功')


    @report_check
    def add_mvirdb_abund(self, path, mvir_id, level):
        if not isinstance(mvir_id, ObjectId):
            mvir_id = ObjectId(mvir_id)
        result = self.db['anno_mvir'].find_one({'_id':mvir_id})
        task_id = result['task_id']
        sample_dic = name2id(task_id, 'task')
        with open(path) as f:
            head = f.readline().strip()
            heads = head.split('\t')
            if level == 'factor':
                sp_head = heads[1:-1]
            else:
                sp_head = heads[1:]
            insert_list = []
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                if len(spline) !=  len(heads):
                    continue
                data= {
                    'mvir_id' : mvir_id,
                    'level' : level,
                    'name' : spline[0]
                }
                for i,e in enumerate(sp_head,1):
                    if e != 'Total':
                        e = sample_dic[e]
                    data[e] = int(spline[i])
                insert_list.append(data)

            try:
                collection = self.db['anno_mvir_abund']
                collection.insert_many(insert_list)
                self.db['anno_mvir'].update_one({'_id': mvir_id}, {'$set':{'main_id':mvir_id}})
            except Exception as e:
                self.bind_object.logger.info('导入 anno_mvir_abund 失败')
                self.bind_object.set_error('import anno_mvir_abund ERROR', code="52804103")
            else:
                self.bind_object.logger.info('导入 anno_mvir_abund 成功')









