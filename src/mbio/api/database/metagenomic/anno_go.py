# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify:20181114
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


class AnnoGo(Base):
    def __init__(self, bind_object):
        super(AnnoGo, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_anno_go(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None, name=None,
                      params=None, group_id=None, group_detail=None,is_origin =None):
        if not isinstance(geneset_id, ObjectId):  # 检查传入的geneset_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52804301")
        else:
            geneset_id_str = str(geneset_id)
        if main:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            if params == None:
                if group_id != None:
                    if not isinstance(group_id, types.StringTypes):
                        if isinstance(group_id, ObjectId):
                            group_id = str(group_id)
                        else:
                            self.bind_object.set_error('geneset_id must be string！', code="52804302")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52804303")
                params = {
                    "database": "GO",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "submit_location": "anno_go",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "GO_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'level4'
            }
            if is_origin == 1:
                insert_data["is_origin"] =int(1)
            try:
                collection = self.db['anno_go']
                go_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.set_error('导入anno_go主表异常:%s', variables=(e), code="52804304")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52804305")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52804306")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                go_id = main_table_id
            except Exception as e:
                self.bind_object.set_error('更新anno_go主表anno_file_path出错:%s', variables=(e), code="52804307")
        return go_id

    @report_check
    def add_go_func(self,go_id,file,level,main=True):
        data_list = []
        result = self.db['anno_go'].find_one({'_id': ObjectId(go_id)})
        samples_dic = {}
        if not result:
            self.bind_object.set_error('not find anno_go file id !', code="52804308")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        with open(file, "r") as f:
            head = f.next()
            sams = []
            if "GO (Lev1)" in head:
                heads = head.strip().split("\t")
                if level == 59:
                    sams = heads[1:len(heads)]
                elif level == 60:
                    sams = heads[3:len(heads)]
                elif level == 61:
                    sams = heads[5:len(heads)]
                elif level == 62:
                    sams = heads[7:len(heads)]
            else:
                self.bind_object.set_error('go_level_profile file error！', code="52804309")
            n = 1
            for line in f:
                lin = line.strip('\r\n').split('\t')
                insert_data = {
                    'go_id': ObjectId(go_id),
                    'goterm': lin[0],
                    }
                if level == 59:
                    insert_data['level'] = 59
                    for i in range(0, len(sams)):
                        if not sams[i] == "Total":
                            sample_id = samples_dic[sams[i]]
                        else:
                            sample_id = sams[i]
                        insert_data[sample_id] = float(lin[(i + 1)])
                elif level == 60:
                    insert_data['goterm_2'] = lin[1]
                    insert_data['goid_2'] = lin[2]
                    insert_data['level'] = 60
                    if main:
                        insert_data['goterm_2_name'] = 'L2_' + str(n)
                    for i in range(0, len(sams)):
                        if not sams[i] == "Total":
                            sample_id = samples_dic[sams[i]]
                        else:
                            sample_id = sams[i]
                        insert_data[sample_id] = float(lin[(i + 3)])
                elif level == 61:
                    insert_data['goterm_2'] = lin[1]
                    insert_data['goid_2'] = lin[2]
                    insert_data['goterm_3'] = lin[3]
                    insert_data['goid_3'] = lin[4]
                    insert_data['level'] = 61
                    if main:
                        insert_data['goterm_3_name'] = 'L3_' + str(n)
                    for i in range(0, len(sams)):
                        if not sams[i] == "Total":
                            sample_id = samples_dic[sams[i]]
                        else:
                            sample_id = sams[i]
                        insert_data[sample_id] = float(lin[(i + 5)])
                elif level == 62:
                    insert_data['goterm_2'] = lin[1]
                    insert_data['goid_2'] = lin[2]
                    insert_data['goterm_3'] = lin[3]
                    insert_data['goid_3'] = lin[4]
                    insert_data['goterm_4'] = lin[5]
                    insert_data['goid_4'] = lin[6]
                    insert_data['level'] = 62
                    if main:
                        insert_data['goterm_4_name'] = 'L4_' + str(n)
                    for i in range(0, len(sams)):
                        if not sams[i] == "Total":
                            sample_id = samples_dic[sams[i]]
                        else:
                            sample_id = sams[i]
                        insert_data[sample_id] = float(lin[(i + 7)])
                data_list.append(insert_data)
                n += 1
        try:
            collection = self.db['anno_go_detail']
            collection.insert_many(data_list)
            self.db['anno_go'].update_one({'_id':ObjectId(go_id)},{'$set':{'main_id':ObjectId(go_id)}})
        except Exception as e:
            self.bind_object.logger.error("导入表格%s信息出错:%s" % (file, e))
        else:
            self.bind_object.logger.info("导入表格%s信息成功!" % file)

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_go'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})