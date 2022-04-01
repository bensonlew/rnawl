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


class Probiotics(Base):
    def __init__(self, bind_object):
        super(Probiotics, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_anno_probio(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None,
            name=None, params=None, group_id=None, group_detail=None,nr_method = "best_hit",is_origin=None, task_id=None):
        # 主表, 所有的函数名称以add开头，里面可以加需要导入数据库而表格里没有的信息作为参数
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_probio_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes): # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52803401")
        else:
            geneset_id_str = str(geneset_id)
        #if not os.path.exists(anno_file_path):  # 调整为上传的永久路径，由于先导表，暂不检查
            #raise Exception('anno_file_path所指定的路径不存在，请检查！')
        if main:
            if not task_id:
                task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "Probio_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'Probiotic_name',
                "is_origin": is_origin,
                "nr_method": nr_method
            }
            try:
                collection = self.db['anno_probio']
                anno_probio_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_ardb主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_ardb主表异常", code="52803402")
            else:
                self.bind_object.logger.info('更新anno_probio_id主表成功:{}'.format(anno_probio_id))
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52803403")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52803404")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_probio_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_probio_id主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_probio_id主表anno_file_path出错", code="52803405")
        return anno_probio_id

    @report_check
    def add_anno_probio_detail(self, anno_probio_id, anno_profile, update_main=False):
        if not isinstance(anno_probio_id, ObjectId):
            if isinstance(anno_probio_id, types.StringTypes):
                anno_probio_id = ObjectId(anno_probio_id)
            else:
                self.bind_object.set_error('anno_probio_id必须为ObjectId对象或其对应的字符串！', code="52803406")
        if not os.path.exists(anno_profile):
            self.bind_object.set_error('ardb_profile所指定的路径不存在，请检查！', code="52803407")
        data_list = []
        result = self.db['anno_probio'].find_one({'_id': anno_probio_id})
        if not result:
            self.bind_object.set_error('找不到probio_anno.xls对应的主表id', code="52803408")
        else:
            task_id = result['task_id']
            #task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
            #print samples_dic
        with open(anno_profile, 'rb') as f:
            head = f.next()
            if "#Probiotic" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads)]
            else:
                self.bind_object.set_error('probio_anno.xls文件错误！', code="52803409")
            for line in f:
                line = line.strip().split('\t')
                link_id = line[-1]
                probiotics = line[0]
                brand = line[1]
                strain = line[2]
                dev_stage = line[4]
                use_in = line[5]
                effect = line[6]
                dis_class = line[7]
                disease = line[8]
                taxon= line[9]
                insert_data = {
                    'probio_id': anno_probio_id,
                    'probiotics': probiotics,
                    'brand': brand,
                    'strain': strain,
                    'dev_stage': dev_stage,
                    'use_in': use_in,
                    'dis_class': dis_class,
                    'effect': effect,
                    'disease': disease,
                    'link_id': link_id,
                    'taxon': taxon,
                }
                data_list.append(insert_data)
            try:
                collection = self.db['anno_probio_detail']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入ardb_profile信息出错", code="52803410")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % anno_profile)
            if update_main:
                main_table = self.db['anno_probio']
                if "Total" in sams:
                    sams.remove("Total")
                specimen = ",".join([samples_dic[i] for i in sams])
                main_table.update_one({'_id': ObjectId(anno_probio_id)}, {'$set': {'specimen': specimen}})

    @report_check
    def add_anno_probio_abun(self, anno_probio_id, abu_profile):
        if not isinstance(anno_probio_id, ObjectId):
            if isinstance(anno_probio_id, types.StringTypes):
                anno_probio_id = ObjectId(anno_probio_id)
            else:
                self.bind_object.set_error('anno_probio_id必须为ObjectId对象或其对应的字符串！', code="52803411")
        if not os.path.exists(abu_profile):
            self.bind_object.set_error('ardb_profile所指定的路径不存在，请检查！', code="52803412")
        data_list = []
        result = self.db['anno_probio'].find_one({'_id': anno_probio_id})
        if not result:
            self.bind_object.set_error('找不到ardb_type_profile对应的主表id', code="52803413")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        with open(abu_profile, 'rb') as f:
            head = f.next()
            if "#Probiotics" in head:
                heads = head.strip().split("\t")
                sams = heads[2:len(heads) - 1]
                self.bind_object.logger.info(sams)
            for line in f:
                line = line.strip().split('\t')
                probiotics = line[0]
                genus = line[1]
                insert_data = {
                    'probio_id': anno_probio_id,
                    'probiotics': probiotics,
                    'genus': genus,
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 2])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_probio_abun']
                collection.insert_many(data_list)
                self.db['anno_probio'].update_one({'_id': anno_probio_id}, {'$set':{'main_id':anno_probio_id}})
            except Exception as e:
                self.bind_object.set_error("导入表格ardb_profile出错", code="52803414")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % abu_profile)

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_probio'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="52803415")
        return object_id
