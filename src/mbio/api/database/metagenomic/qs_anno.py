# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify:20180914
from biocluster.api.database.base import Base, report_check
import os,re
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name
import json


class QsAnno(Base):
    def __init__(self, bind_object):
        super(QsAnno, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_anno_qs(self, geneset_id, specimen, anno_file_path,main=True, main_table_id=None, name=None,
                      params=None, group_id=None, group_detail=None,is_origin =None):
        if not isinstance(geneset_id, ObjectId):  # 检查传入的geneset_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52804901")
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
                            self.bind_object.set_error('geneset_id must be string！', code="52804902")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52804903")
                params = {
                    "database": "qs",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "identity": 0,
                    "align_length": 0,
                    "submit_location": "annoqs",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': 'QS群体感应分析',
                'created_ts': created_ts,
                'name': name if name else "QS_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'QS_id'
            }
            if is_origin == 1:
                insert_data["is_origin"] =int(1)
            try:
                collection = self.db['anno_qs']
                qs_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.set_error('导入anno_qs主表异常:%s', variables=(e), code="52804904")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52804905")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52804906")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                qs_id = main_table_id
            except Exception as e:
                self.bind_object.set_error('更新anno_qs主表anno_file_path出错:%s', variables=(e), code="52804907")
        return qs_id

    @report_check
    def add_anno_qs_class(self, qs_id, qs_dir):
        qs_profile = qs_dir + "/qs_class_profile.xls"
        if not isinstance(qs_id, ObjectId):
            if isinstance(qs_id, types.StringTypes):
                qs_id = ObjectId(qs_id)
            else:
                self.bind_object.set_error('qs_id必须为ObjectId对象或其对应的字符串！', code="52804908")
        if not os.path.exists(qs_profile):
            self.bind_object.set_error('qs_profile所指定的路径不存在，请检查！', code="52804909")
        data_list = []
        result = self.db['anno_qs'].find_one({'_id': qs_id})
        if not result:
            self.bind_object.set_error('找不到anno_qs对应的主表id', code="52804910")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        with open(qs_profile, 'rb') as f:
            head = f.next()
            if "#Class" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('qs_class_profile.xls文件错误！', code="52804911")
            for line in f:
                line = line.strip().split('\t')
                qsclass = line[0]
                insert_data = {
                    'qs_id': qs_id,
                    'class': qsclass,
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_qs_detail']
                collection.insert_many(data_list)
                self.db['anno_qs'].update_one({'_id': qs_id}, {'$set':{'main_id':qs_id}})
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (qs_profile, e))
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % qs_profile)

    def add_qs_graph(self,qs_id,file):
        data_list = []
        with open(file, 'rb') as f:
            for line in f:
                line = line.strip('\r\n').split('\t')
                des =line[2].split(',')
                des = map(float, des)
                insert_data = {
                    'qs_id': ObjectId(qs_id),
                    'group_label': line[0],
                    'func':line[1],
                    'data':des
                }
                if not re.search('-', line[3]):
                    dess =''
                    if re.search(',', line[3]):
                        dess = line[3].split(',')
                        dess = map(float, dess)
                    else:
                        dess=[]
                        dess.append(float(line[3]))
                        dess = map(float, dess)
                    insert_data['outlier'] = dess
                data_list.append(insert_data)
        try:
            collection = self.db['anno_qs_graph']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入表格%s信息出错:%s" % (file, e))
        else:
            self.bind_object.logger.info("导入表格%s信息成功!" % file)



    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_qs'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
