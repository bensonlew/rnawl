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


class MgAnnoCazy(Base):
    def __init__(self, bind_object):
        super(MgAnnoCazy, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_anno_cazy(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None, name=None,
                      params=None, group_id=None, group_detail=None, software_ver={}):
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_cazy_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52801901")
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
                            self.bind_object.set_error('geneset_id必须为字符串类型！', code="52801902")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52801903")
                params = {
                    "database": "cazy",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "identity": 0,
                    "align_length": 0,
                    "submit_location": "annocazy",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "CAZy_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'Family',
                "settled_params": json.dumps({"version": "cazy_v20200408"})
            }
            if software_ver:
                insert_data.update(software_ver)
            try:
                collection = self.db['anno_cazy']
                anno_cazy_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_cazy主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_cazy主表异常", code="52801904")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52801905")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52801906")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_cazy_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_cazy主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_file_path出错", code="52801907")
        return anno_cazy_id

    @report_check
    def add_anno_cazy_family(self, anno_cazy_id, cazy_profile_dir, update_main=True):
        cazy_profile = cazy_profile_dir + "/cazy_family_profile.xls"
        if not isinstance(anno_cazy_id, ObjectId):
            if isinstance(anno_cazy_id, types.StringTypes):
                anno_cazy_id = ObjectId(anno_cazy_id)
            else:
                self.bind_object.set_error('anno_cazy_id必须为ObjectId对象或其对应的字符串！', code="52801908")
        if not os.path.exists(cazy_profile):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('cazy_profile所指定的路径不存在，请检查！', code="52801909")
        data_list = []
        result = self.db['anno_cazy'].find_one({'_id': anno_cazy_id})
        if not result:
            self.bind_object.set_error('找不到anno_cazy_nog对应的主表id', code="52801910")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(cazy_profile, 'rb') as f:
            head = f.next()
            if "#Family" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('cazy_family_profile.xls文件错误！', code="52801911")
            for line in f:
                line = line.strip().split('\t')
                family = line[0]
                description = line[-1]
                insert_data = {
                    'cazy_id': anno_cazy_id,
                    'family': family,
                    'description': description
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_cazy_family']
                collection.insert_many(data_list)
                main_table = self.db['anno_cazy']
                main_table.update_one({'_id': ObjectId(anno_cazy_id)}, {'$set': {"settled_params": json.dumps({"version": "cazy_v20200408"})}})
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (cazy_profile, e))
                self.bind_object.set_error("导入cazy_profile信息出错", code="52801912")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % cazy_profile)
            if update_main:
                main_table = self.db['anno_cazy']
                if "Total" in sams:
                    sams.remove("Total")
                specimen = ",".join([samples_dic[i] for i in sams])
                # specimen = ",".join(sams)
                main_table.update_one({'_id': ObjectId(anno_cazy_id)}, {'$set': {'specimen': specimen,
                                                                                 "settled_params": json.dumps({"version": "cazy_v20200408"})}})

    @report_check
    def add_anno_cazy_class(self, anno_cazy_id, cazy_profile_dir):
        cazy_profile = cazy_profile_dir + "/cazy_class_profile.xls"
        if not isinstance(anno_cazy_id, ObjectId):
            if isinstance(anno_cazy_id, types.StringTypes):
                anno_cazy_id = ObjectId(anno_cazy_id)
            else:
                self.bind_object.set_error('anno_cazy_id必须为ObjectId对象或其对应的字符串！', code="52801908")
        if not os.path.exists(cazy_profile):
            self.bind_object.set_error('cazy_profile 所指定的路径不存在，请检查！', code="52801909")
        data_list = []
        result = self.db['anno_cazy'].find_one({'_id': anno_cazy_id})
        if not result:
            self.bind_object.set_error('找不到anno_cazy_nog对应的主表id', code="521801910")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(cazy_profile_dir + "/cazy_class_profile.xls", 'rb') as f:
            head = f.next()
            if "#Class" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('cazy_class_profile.xls文件错误！', code="52801913")
            for line in f:
                line = line.strip().split('\t')
                classes = line[0]
                description = line[-1]
                insert_data = {
                    'cazy_id': anno_cazy_id,
                    'class': classes,
                    'description': description
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_cazy_class']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (cazy_profile, e))
                self.bind_object.set_error("导入cazy_profile信息出错", code="52801914")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % cazy_profile)

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_cazy'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
