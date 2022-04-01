# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20171114
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
import json
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name


class MgAnnoCog(Base):
    def __init__(self, bind_object):
        super(MgAnnoCog, self).__init__(bind_object)
        # self._db_name = Config().MONGODB + '_metagenomic'
        self._project_type = "metagenomic"

    @report_check
    def add_anno_cog(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None, name=None,
                     params=None, group_id=None, group_detail=None, software_ver={}):
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_cog_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52802001")
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
                            self.bind_object.set_error('geneset_id必须为字符串类型！', code="52802002")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52802003")
                params = {
                    "database": "cog",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "identity": 0,
                    "align_length": 0,
                    "submit_location": "annocog",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "COG_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'NOG',
                "settled_params": json.dumps({"version": "eggNOG_v4.5.1"})
            }
            if software_ver:
                insert_data.update(software_ver)
            try:
                collection = self.db['anno_cog']
                anno_cog_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_cog主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_cog主表异常", code="52802004")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52802005")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52802006")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_cog_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_cog主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_cog主表anno_file_path出错", code="52802007")
        return anno_cog_id

    @report_check
    def add_anno_cog_nog(self, anno_cog_id, cog_profile_dir, update_main=True):
        cog_profile = cog_profile_dir + "/cog_nog_profile.xls"
        if not isinstance(anno_cog_id, ObjectId):  # 检查传入的anno_cog_id是否符合ObjectId类型
            if isinstance(anno_cog_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                anno_cog_id = ObjectId(anno_cog_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_cog_id必须为ObjectId对象或其对应的字符串！', code="52802008")
        if not os.path.exists(cog_profile):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('cog_nog_profile所指定的路径不存在，请检查！', code="52802009")
        data_list = []
        result = self.db['anno_cog'].find_one({'_id': anno_cog_id})
        if not result:
            self.bind_object.set_error('找不到anno_cog_nog对应的主表id', code="52802010")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(cog_profile, 'rb') as f:
            head = f.next()
            if "#NOG" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('cog_nog_profile.xls文件错误！', code="52802011")
            for line in f:
                line = line.strip().split('\t')
                nog = line[0]
                des = line[-1]
                insert_data = {
                    'cog_id': anno_cog_id,
                    'nog': nog,
                    'description': des
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_cog_nog']
                collection.insert_many(data_list)
                main_table = self.db['anno_cog']
                main_table.update_one({'_id': ObjectId(anno_cog_id)}, {'$set': {"settled_params": json.dumps({"version": "eggNOG_v4.5.1"})}})
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (cog_profile, e))
                self.bind_object.set_error("导入cog_profile信息出错", code="52802012")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % cog_profile)
            if update_main:
                main_table = self.db['anno_cog']
                if "Total" in sams:
                    sams.remove("Total")
                specimen = ",".join([samples_dic[i] for i in sams])
                # specimen = ",".join(sams)
                main_table.update_one({'_id': ObjectId(anno_cog_id)}, {'$set': {'specimen': specimen,
                                                                                "settled_params": json.dumps({"version": "eggNOG_v4.5.1"})}})

    @report_check
    def add_anno_cog_function(self, anno_cog_id, cog_profile_dir):
        cog_profile = cog_profile_dir + "/cog_function_profile.xls"
        if not isinstance(anno_cog_id, ObjectId):  # 检查传入的anno_cog_id是否符合ObjectId类型
            if isinstance(anno_cog_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                anno_cog_id = ObjectId(anno_cog_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_cog_id必须为ObjectId对象或其对应的字符串！', code="52802008")
        if not os.path.exists(cog_profile):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('cog_function_profile所指定的路径不存在，请检查！', code="52802013")
        data_list = []
        result = self.db['anno_cog'].find_one({'_id': anno_cog_id})
        if not result:
            self.bind_object.set_error('找不到anno_cog_nog对应的主表id', code="52802010")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(cog_profile, 'rb') as f:
            head = f.next()
            if "#Function" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('cog_function_profile文件错误！', code="52802014")
            for line in f:
                line = line.strip().split('\t')
                function = line[0]
                des = line[-1]
                insert_data = {
                    'cog_id': anno_cog_id,
                    'function': function,
                    'description': des
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_cog_function']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (cog_profile, e))
                self.bind_object.set_error("导入cog_profile信息出错", code="52802015")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % cog_profile)

    @report_check
    def add_anno_cog_category(self, anno_cog_id, cog_profile_dir):
        cog_profile = cog_profile_dir + "/cog_category_profile.xls"
        if not isinstance(anno_cog_id, ObjectId):  # 检查传入的anno_cog_id是否符合ObjectId类型
            if isinstance(anno_cog_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                anno_cog_id = ObjectId(anno_cog_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_cog_id必须为ObjectId对象或其对应的字符串！', code="52802008")
        if not os.path.exists(cog_profile):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('cog_category_profile所指定的路径不存在，请检查！', code="52802016")
        data_list = []
        result = self.db['anno_cog'].find_one({'_id': anno_cog_id})
        if not result:
            self.bind_object.set_error('找不到anno_cog_nog对应的主表id', code="52802010")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(cog_profile, 'rb') as f:
            head = f.next()
            if "#Category" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads)]
            else:
                self.bind_object.set_error('cog_category_profile文件错误！', code="52802017")
            for line in f:
                line = line.strip().split('\t')
                category = line[0]
                insert_data = {
                    'cog_id': anno_cog_id,
                    'category': category
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_cog_category']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (cog_profile, e))
                self.bind_object.set_error("导入cog_profile信息出错", code="52802018")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % cog_profile)

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_cog'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
