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

class MgAnnoVfdb(Base):
    def __init__(self, bind_object):
        super(MgAnnoVfdb, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_anno_vfdb(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None, name=None,
                      params=None, group_id=None, group_detail=None, software_ver={}):
        # 主表, 所有的函数名称以add开头，里面可以加需要导入数据库而表格里没有的信息作为参数
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_vfdb_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52802401")
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
                            self.bind_object.set_error('group_id必须为字符串类型！: 内容为%s\n类型为%s',
                                                       variables=(group_id,type(group_id)), code="52802402")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52802403")
                params = {
                    "database": "vfdb",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "align_database" : "All",
                    "identity": 0,
                    "align_length": 0,
                    "submit_location": "annovfdb",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "VFDB_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'VFs',
                "group_name": specimen,
                "settled_params": json.dumps({"version": "vfdb_v20200703"})
            }
            if software_ver:
                insert_data.update(software_ver)
            try:
                collection = self.db['anno_vfdb']
                anno_vfdb_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_vfdb主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_vfdb主表异常", code="52802404")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52802405")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52802406")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_vfdb_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_vfdb主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_vfdb主表anno_file_path出错", code="52802407")
        return anno_vfdb_id

    @report_check
    def add_anno_vfdb_vfs(self, anno_vfdb_id, vfdb_profile_dir, database, update_main=True):
        if not isinstance(anno_vfdb_id, ObjectId):
            if isinstance(anno_vfdb_id, types.StringTypes):
                anno_vfdb_id = ObjectId(anno_vfdb_id)
            else:
                self.bind_object.set_error('anno_vfdb_id必须为ObjectId对象或其对应的字符串！', code="52802408")
        if database == "core":
            vfdb_profile = vfdb_profile_dir + "/vfdb_core_VF_profile.xls"
        elif database == "predict":
            vfdb_profile = vfdb_profile_dir + "/vfdb_predict_VF_profile.xls"
        #elif database == "all":
            #vfdb_profile = vfdb_profile_dir + "/vfdb_all_VF_profile.xls"
        if not os.path.exists(vfdb_profile):
            self.bind_object.set_error('vfdb_profile所指定的路径不存在，请检查！', code="52802409")
        data_list = []
        result = self.db['anno_vfdb'].find_one({'_id': anno_vfdb_id})
        if not result:
            self.bind_object.set_error('找不到vfdb_profile对应的主表id', code="52802410")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(vfdb_profile, 'rb') as f:
            head = f.next()
            if "#VFs" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 4]
            for line in f:
                line = line.strip().split('\t')
                vfs = line[0]
                species = line[-4]
                function = line[-3]
                level1 = line[-2]
                level2 = line[-1]
                insert_data = {
                    'vfdb_id': anno_vfdb_id,
                    'vfs': vfs,
                    'species': species,
                    'function': function,
                    'level1': level1,
                    'level2': level2,
                    'data_type': database
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_vfdb_vfs']
                collection.insert_many(data_list)
                main_table = self.db['anno_vfdb']
                main_table.update_one({'_id': ObjectId(anno_vfdb_id)}, {'$set': {"settled_params": json.dumps({"version": "vfdb_v20200703"})}})
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (vfdb_profile, e))
                self.bind_object.set_error("导入vfdb_profile信息出错", code="52802411")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % vfdb_profile)
            if update_main:
                main_table = self.db['anno_vfdb']
                if "Total" in sams:
                    sams.remove("Total")
                specimen = ",".join([samples_dic[i] for i in sams])
                # specimen = ",".join(sams)
                main_table.update_one({'_id': ObjectId(anno_vfdb_id)}, {'$set': {'specimen': specimen,
                                                                                 "settled_params": json.dumps({"version": "vfdb_v20200703"})}})

    @report_check
    def add_anno_vfdb_pie(self, anno_vfdb_id, vfdb_profile_dir):
        vfdb_profile = vfdb_profile_dir + "/vfdb_level_pie.xls"
        if not isinstance(anno_vfdb_id, ObjectId):
            if isinstance(anno_vfdb_id, types.StringTypes):
                anno_vfdb_id = ObjectId(anno_vfdb_id)
            else:
                self.bind_object.set_error('anno_vfdb_id必须为ObjectId对象或其对应的字符串！', code="52802408")
        if not os.path.exists(vfdb_profile):
            self.bind_object.set_error('vfdb_profile所指定的路径不存在，请检查！', code="52402409")
        data_list = []
        result = self.db['anno_vfdb'].find_one({'_id': anno_vfdb_id})
        if not result:
            self.bind_object.set_error('找不到vfdb_profile对应的主表id', code="52802410")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        with open(vfdb_profile, 'rb') as f:
            lines = f.readlines()
            heads = lines[0].strip().split("\t")
            #sam_num = len(heads) / 2 - 1
            sams = heads[2:len(heads)]
            #sams = heads[2:(sam_num + 2)]
            # sam_pers = heads[(sam_num + 2):len(heads)]
            for line in lines[1:]:
                line = line.strip().split('\t')
                level1 = line[0]
                level2 = line[1]
                # abu = line[2]
                # percent = line[3]
                insert_data = {
                    'vfdb_id': anno_vfdb_id,
                    'level1': level1,
                    'level2': level2,
                }
                for i in range(0, len(sams)):
                    if samples_dic.has_key(sams[i]):
                        sample_id = samples_dic[sams[i]]
                        #sam_pers_id = sample_id + "_percent"
                    else:
                        sample_id = sams[i]
                        #sam_pers_id = sample_id + "_percent"
                    insert_data[sample_id] = float(line[i + 2])
                    #insert_data[sam_pers_id] = float(line[(i + 2 + sam_num)])
                data_list.append(insert_data)
            if len(lines) > 1:
                try:
                    collection = self.db['anno_vfdb_pie']
                    collection.insert_many(data_list)
                except Exception as e:
                    self.bind_object.logger.error("导入表格%s信息出错:%s" % (vfdb_profile, e))
                    self.bind_object.set_error("导入vfdb_profile信息出错", code="52802411")
                else:
                    self.bind_object.logger.info("导入表格%s信息成功!" % vfdb_profile)
            else:
                self.bind_object.logger.info("注释结果pie图数据为空")
    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_vfdb'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
