# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180329
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


class Contribute(Base):
    def __init__(self, bind_object):
        super(Contribute, self).__init__(bind_object)
        # self._db_name = Config().MONGODB + '_metagenomic'
        self._project_type = "metagenomic"

    @report_check
    def add_contribute(self, geneset_id, specimen, anno_type, group, main=True, main_table_id=None, name=None,
                       params=None):
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52800501")
        if main:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            if params != None:
                params = params
                print params
            else:
                params = ""
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "Contribute_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_type': anno_type,
                'group': group,
            }
            try:
                collection = self.db['contribute']
                contribute_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入contribute主表异常:{}'.format(e))
                self.bind_object.set_error("导入contribute主表异常", code="52800502")
            print contribute_id
        return contribute_id

    @report_check
    def add_contribute_detail(self, contribute_id, tax_fun_file, fun_tax_file, update_main=True):
        if not isinstance(contribute_id, ObjectId):  # 检查传入的Contribute_id是否符合ObjectId类型
            if isinstance(contribute_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                contribute_id = ObjectId(contribute_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('contribute_id必须为ObjectId对象或其对应的字符串！', code="52800503")
        if not os.path.exists(tax_fun_file):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(tax_fun_file), code="52800504")
        if not os.path.exists(fun_tax_file):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(fun_tax_file), code="52800504")
        data_list = []
        print contribute_id
        result = self.db['contribute'].find_one({'_id': contribute_id})
        groups = result["group"].split(",")
        len_groups = len(groups)
        if not result:
            self.bind_object.set_error('找不到corr_network_by_cut.txt对应的主表id', code="52800505")
        else:
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
        self.bind_object.logger.error("开始tax_fun_file")
        with open(tax_fun_file, 'rb') as f:
            head = f.next()
            heads = head.strip().split("\t")
            abu = heads[2:len(heads)]
            self.bind_object.logger.info("len_abu:%s" % (len(abu)))
            #self.bind_object.logger.error(abu)
            for line in f:
                line = line.strip().split('\t')
                taxon = line[0]
                function = line[1]
                insert_data = {
                    'contribute_id': contribute_id,
                    'taxon': taxon,
                    'function': function,
                    'type': 2
                }
                for i in range(0, len(abu)-len_groups-1):
                    if samples_dic.has_key(abu[i]):
                        sample_id = samples_dic[abu[i]]
                    else:
                        sample_id = abu[i]
                    insert_data[sample_id] = float(line[i + 2])
                ### modification when sample name is name with group name
                for j in range(len(abu)-len_groups,len(abu)):
                    sample_id = abu[j]
                    splits = sample_id.split("_")
                    sample_id = "_".join(splits[1:len(splits)])
                    insert_data[sample_id] = float(line[j + 2])
                data_list.append(insert_data)
        self.bind_object.logger.error("开始fun_tax_file")
        with open(fun_tax_file, 'rb') as f1:
            head = f1.next()
            heads = head.strip().split("\t")
            abu = heads[2:len(heads)]
            for line in f1:
                line = line.strip().split('\t')
                taxon = line[1]
                function = line[0]
                insert_data = {
                    'contribute_id': contribute_id,
                    'taxon': taxon,
                    'function': function,
                    'type': 1
                }
                for i in range(0, len(abu)-len_groups-1):
                    if samples_dic.has_key(abu[i]):
                        sample_id = samples_dic[abu[i]]
                    else:
                        sample_id = abu[i]
                    insert_data[sample_id] = float(line[i + 2])
                ### modification when sample name is name with group name
                for j in range(len(abu)-len_groups,len(abu)):
                    sample_id = abu[j]
                    splits = sample_id.split("_")
                    sample_id = "_".join(splits[1:len(splits)])
                    insert_data[sample_id] = float(line[j + 2])
                data_list.append(insert_data)
            try:
                collection = self.db['contribute_detail']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入contribute_detail表格信息出错:%s" % (e))
                self.bind_object.set_error("导入contribute_detail表信息出错", code="52800506")
            else:
                self.bind_object.logger.info("导入contribute_detail表格信息成功!")
