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
import pymongo

class MgAnnoCard(Base):
    def __init__(self, bind_object):
        super(MgAnnoCard, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_anno_card(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None, name=None,
                      params=None, group_id=None, group_detail=None, software_ver={}):
        if not isinstance(geneset_id, ObjectId):  # 检查传入的geneset_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52801801")
        else:
            geneset_id_str = str(geneset_id)
        #if not os.path.exists(anno_file_path):  # # 调整为上传的永久路径，由于先导表，暂不检查
            #raise Exception('anno_file_path所指定的路径不存在，请检查！')
        if main:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            if params == None:
                if group_id != None:
                    if not isinstance(group_id, types.StringTypes):
                        if isinstance(group_id,ObjectId):
                            group_id = str(group_id)
                        else:
                            self.bind_object.set_error('geneset_id必须为字符串类型！', code="52801802")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52801803")
                params = {
                    "database": "card",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "identity": 0,
                    "align_length": 0,
                    "submit_location": "annocard",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "CARD_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'ARO',
                "settled_params": json.dumps({"version": "card_v3.0.9"})
            }
            if software_ver:
                insert_data.update(software_ver)
            try:
                collection = self.db['anno_card']
                anno_card_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_card主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_card主表异常", code="52801804")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52801805")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52801806")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_card_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_card主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_file_path出错", code="52801807")
        return anno_card_id

    @report_check
    def add_anno_card_aro(self, anno_card_id, card_profile_dir, update_main=True):
        task_id = "_".join(self.bind_object.sheet.id.split("_")[0:2])
        self.old_task = 0
        self.bind_object.logger.info("task_id : {}".format(task_id))
        sg_task_result = self.db['sg_task'].find_one({"task_id": task_id})
        if sg_task_result:
            if "database" in sg_task_result:
                self.bind_object.logger.info("sg_task_result : {}".format(sg_task_result['database']))
                card_profile = card_profile_dir + "/card_ARO_profile_all.xls"
            else:
                card_profile = card_profile_dir + "/card_ARO_profile.xls"
                self.old_task = 1
        else:
            card_profile = card_profile_dir + "/card_ARO_profile.xls"
            self.old_task = 1
        if not isinstance(anno_card_id, ObjectId):
            if isinstance(anno_card_id, types.StringTypes):
                anno_card_id = ObjectId(anno_card_id)
            else:
                self.bind_object.set_error('anno_card_id必须为ObjectId对象或其对应的字符串！', code="52801808")
        self.bind_object.logger.info("self.old_task: {}".format(self.old_task))
        if not os.path.exists(card_profile):
            self.bind_object.set_error('card_profile所指定的路径不存在，请检查！', code="52801809")
        data_list = []
        result = self.db['anno_card'].find_one({'_id': anno_card_id})
        if not result:
            self.bind_object.set_error('找不到anno_card_nog对应的主表id', code="52801810")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(card_profile, 'rb') as f:
            head = f.next()
            if "#ARO" in head:
                heads = head.strip().split("\t")
                # sams = heads[1:len(heads) - 2]
                if self.old_task == 1 :
                    sams = heads[2:len(heads) - 1]
                else:
                    sams = heads[7:len(heads)]
            else:
                self.bind_object.set_error('card_ARO_profile.xls文件错误！', code="52801811")
            for line in f:
                line = line.strip().split('\t')
                if self.old_task == 1 : ## 兼容老的数据库版本 add by qingchen.zhang
                    ARO = line[0]
                    aro_name = line[1]
                    description = line[-1]
                    insert_data = {
                        'card_id': anno_card_id,
                        'aro': ARO,
                        'aro_name': aro_name,
                        'description': description
                    }
                    for i in range(0, len(sams)):
                        if not sams[i] == "Total":
                            sample_id = samples_dic[sams[i]]
                        else:
                            sample_id = sams[i]
                        insert_data[sample_id] = float(line[i + 2])

                else:
                    ARO = line[0]
                    aro_name = line[1]
                    description = line[2]
                    amr_gene_family = line[3]
                    drug_class = line[4]
                    type = line[5]
                    resistance = line[6]
                    insert_data = {
                        'card_id': anno_card_id,
                        'aro': ARO,
                        'aro_name': aro_name,
                        'description': description,
                        'amr_gene_family':amr_gene_family,
                        'drug_class':drug_class,
                        'type':type,
                        'resistance':resistance
                    }
                    for i in range(0, len(sams)):
                        if not sams[i] == "Total":
                            sample_id = samples_dic[sams[i]]
                        else:
                            sample_id = sams[i]
                        insert_data[sample_id] = float(line[i + 7])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_card_aro']
                collection.insert_many(data_list)
                main_table = self.db['anno_card']
                if self.old_task == 0:
                    main_table.update_one({'_id': ObjectId(anno_card_id)}, {'$set': {"settled_params": json.dumps({"version": "card_v3.0.9"})}})
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (card_profile, e))
                self.bind_object.set_error("导入card_profile信息出错", code="52801812")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % card_profile)
                # collection.ensure_index('nog', unique=False)
        with open(card_profile_dir + "/card_ARO_gene_number.xls", 'rb') as f1:
            head = f1.next()
            for line in f1:
                line = line.strip().split('\t')
                ARO = line[0]
                gene_number = line[1]
                try:
                    collection.update({"card_id": anno_card_id, 'aro': ARO}, {'$set': {'gene_number': gene_number}})
                except Exception as e:
                    self.bind_object.logger.error("更新%s信息gene_number出错:%s" % (card_profile, e))
                    self.bind_object.set_error("更新card_profile信息出错", code="52801813")
                else:
                    self.bind_object.logger.info("更新%s信息gene_number成功!" % card_profile)
            if update_main:
                main_table = self.db['anno_card']
                if "Total" in sams:
                    sams.remove("Total")
                specimen = ",".join([samples_dic[i] for i in sams])
                # specimen = ",".join(sams)
                main_table.update_one({'_id': ObjectId(anno_card_id)}, {'$set': {'specimen': specimen}})

    @report_check
    def add_anno_card_class(self, anno_card_id, card_profile_dir):
        card_profile = card_profile_dir + "/card_class_profile.xls"
        if not isinstance(anno_card_id, ObjectId):
            if isinstance(anno_card_id, types.StringTypes):
                anno_card_id = ObjectId(anno_card_id)
            else:
                self.bind_object.set_error('anno_card_id必须为ObjectId对象或其对应的字符串！', code="52801814")
        if not os.path.exists(card_profile):
            self.bind_object.set_error('card_profile所指定的路径不存在，请检查！', code="52801815")
        data_list = []
        result = self.db['anno_card'].find_one({'_id': anno_card_id})
        if not result:
            self.bind_object.set_error('找不到anno_card_nog对应的主表id', code="52801816")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
            print samples_dic
        with open(card_profile, 'rb') as f:
            head = f.next()
            if "#Class" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            for line in f:
                line = line.strip().split('\t')
                classes = line[0]
                description = line[-1]
                insert_data = {
                    'card_id': anno_card_id,
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
                collection = self.db['anno_card_class']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (card_profile, e))
                self.bind_object.set_error("导入card_profile信息出错", code="52801817")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % card_profile)

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_card'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
