# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20170923
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


class MgAnnoArdb(Base):
    def __init__(self, bind_object):
        super(MgAnnoArdb, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_anno_ardb(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None,
                      name=None, params=None, group_id=None, group_detail=None, software_ver={}):
        # 主表, 所有的函数名称以add开头，里面可以加需要导入数据库而表格里没有的信息作为参数
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_ardb_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes): # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52801701")
        else:
            geneset_id_str = str(geneset_id)
        #if not os.path.exists(anno_file_path):  # 调整为上传的永久路径，由于先导表，暂不检查
            #raise Exception('anno_file_path所指定的路径不存在，请检查！')
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
                            self.bind_object.set_error('geneset_id必须为字符串类型！', code="52801702")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52801703")
                params = {
                    "database": "ardb",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "identity": 0,
                    "align_length": 0,
                    "submit_location": "annoardb",
                    "task_type": 2
                }
                print params
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "ARDB_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'ARG',
                "settled_params": json.dumps({"version": "ARDB_v1.1"})
            }
            if software_ver:
                insert_data.update(software_ver)
            try:
                collection = self.db['anno_ardb']
                anno_ardb_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_ardb主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_ardb主表异常", code="52801704")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52801705")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52801706")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_ardb_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_ardb主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_ardb主表anno_file_path出错", code="52801707")
        return anno_ardb_id

    @report_check
    def add_anno_ardb_arg(self, anno_ardb_id, ardb_profile_dir, update_main=True):
        ardb_profile = ardb_profile_dir + "/ardb_ARG_profile.xls"
        if not isinstance(anno_ardb_id, ObjectId):
            if isinstance(anno_ardb_id, types.StringTypes):
                anno_ardb_id = ObjectId(anno_ardb_id)
            else:
                self.bind_object.set_error('anno_ardb_id必须为ObjectId对象或其对应的字符串！', code="52801708")
        if not os.path.exists(ardb_profile):
            self.bind_object.set_error('ardb_profile所指定的路径不存在，请检查！', code="52801709")
        data_list = []
        result = self.db['anno_ardb'].find_one({'_id': anno_ardb_id})
        if not result:
            self.bind_object.set_error('找不到ardb_ARG_profile对应的主表id', code="52801710")
        else:
            task_id = result['task_id']
            #task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
            #print samples_dic
        with open(ardb_profile, 'rb') as f:
            head = f.next()
            if "#ARG" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads)]
            else:
                self.bind_object.set_error('ardb_ARG_profile.xls文件错误！', code="52801711")
            for line in f:
                line = line.strip().split('\t')
                arg = line[0]
                insert_data = {
                    'ardb_id': anno_ardb_id,
                    'arg': arg
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_ardb_arg']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (ardb_profile, e))
                self.bind_object.set_error("导入ardb_profile信息出错", code="52801712")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % ardb_profile)
            if update_main:
                main_table = self.db['anno_ardb']
                if "Total" in sams:
                    sams.remove("Total")
                specimen = ",".join([samples_dic[i] for i in sams])
                # specimen = ",".join(sams)
                main_table.update_one({'_id': ObjectId(anno_ardb_id)}, {'$set': {'specimen': specimen,
                                                                                 "settled_params": json.dumps({"version": "ARDB_v1.1"})}})

    @report_check
    def add_anno_ardb_type(self, anno_ardb_id, ardb_profile_dir):
        ardb_profile = ardb_profile_dir + "/ardb_type_profile.xls"
        if not isinstance(anno_ardb_id, ObjectId):
            if isinstance(anno_ardb_id, types.StringTypes):
                anno_ardb_id = ObjectId(anno_ardb_id)
            else:
                self.bind_object.set_error('anno_ardb_id必须为ObjectId对象或其对应的字符串！', code="52801713")
        if not os.path.exists(ardb_profile):
            self.bind_object.set_error('ardb_profile所指定的路径不存在，请检查！', code="52801714")
        data_list = []
        result = self.db['anno_ardb'].find_one({'_id': anno_ardb_id})
        if not result:
            self.bind_object.set_error('找不到ardb_type_profile对应的主表id', code="52801715")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        with open(ardb_profile, 'rb') as f:
            head = f.next()
            if "#Type" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            for line in f:
                line = line.strip().split('\t')
                type = line[0]
                antibiotic_type = line[-1]
                insert_data = {
                    'ardb_id': anno_ardb_id,
                    'type': type,
                    'antibiotic_type': antibiotic_type
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_ardb_type']
                collection.insert_many(data_list)
                main_table = self.db['anno_ardb']
                main_table.update_one({'_id': ObjectId(anno_ardb_id)}, {'$set': {"settled_params": json.dumps({"version": "ARDB_v1.1"})}})
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (ardb_profile, e))
                self.bind_object.set_error("导入表格ardb_profile出错", code="52801716")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % ardb_profile)

    @report_check
    def add_anno_ardb_class(self, anno_ardb_id, ardb_profile_dir):
        ardb_profile = ardb_profile_dir + "/ardb_class_profile.xls"
        if not isinstance(anno_ardb_id, ObjectId):
            if isinstance(anno_ardb_id, types.StringTypes):
                anno_ardb_id = ObjectId(anno_ardb_id)
            else:
                self.bind_object.set_error('anno_ardb_id必须为ObjectId对象或其对应的字符串！', code="52801717")
        if not os.path.exists(ardb_profile):
            self.bind_object.set_error('ardb_profile所指定的路径不存在，请检查！', code="52801718")
        data_list = []
        result = self.db['anno_ardb'].find_one({'_id': anno_ardb_id})
        if not result:
            self.bind_object.set_error('找不到ardb_type_profile对应的主表id', code="52801719")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(ardb_profile, 'rb') as f:
            head = f.next()
            if "#Class" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            for line in f:
                line = line.strip().split('\t')
                classname = line[0]
                class_des = line[-1]
                insert_data = {
                    'ardb_id': anno_ardb_id,
                    'class': classname,
                    'class_des': class_des
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_ardb_class']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (ardb_profile, e))
                self.bind_object.set_error("导入ardb_profile信息出错", code="52801720")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % ardb_profile)

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_ardb'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
