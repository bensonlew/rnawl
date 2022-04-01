# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify:20170922
from biocluster.api.database.base import Base, report_check
import os
import re
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
import pandas as pd

class DataStat(Base):
    def __init__(self, bind_object):
        super(DataStat, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.specimen2id = ""
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_data_stat(self, data_type, stat_path, dir_path, raw_data_stat_id, params=None, software_ver={}):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if data_type == "unigene":
            tmp_data_type = "raw"
        else:
            tmp_data_type = data_type
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '数据质控统计主表',
            'name': 'data_stat',
            'created_ts': created_ts,
            'params': "1",  # 符合存pdf逻辑，无实际意义
            # 'params': params,
            'status': 'end',
            'type': tmp_data_type,
        }
        if software_ver:  # 增加软件版本字段 键为软件，值版本
            insert_data.update(software_ver)
        collection = self.db['data_stat']
        # 将主表名称写在这里
        data_stat_id = collection.insert_one(insert_data).inserted_id
        # 将导表数据通过insert_one函数导入数据库，将此条记录生成的_id作为返回值，给detail表做输入参数
        if data_type in ["raw"]:
            self.add_data_stat_detail(stat_path, data_stat_id, data_type, raw_data_stat_id=None)
            if not dir_path == "null":
                if os.path.isdir(dir_path):
                    self.add_specimen_graphic(dir_path, task_id)
                else:
                    self.new_add_specimen_graphic(dir_path, task_id)
        if data_type in ["clean", "optimised"]:
            self.add_data_stat_detail(stat_path, data_stat_id, data_type, raw_data_stat_id)
        if data_type in ["unigene"]:
            self.add_data_stat_detail(stat_path, data_stat_id, "unigene", raw_data_stat_id=None)
            self.add_specimen_graphic(dir_path, task_id)
        return data_stat_id

    @report_check
    def add_data_stat_detail(self, stat_path, data_stat_id, data_type, raw_data_stat_id=None):
        if not isinstance(data_stat_id, ObjectId):
            if isinstance(data_stat_id, types.StringTypes):
                data_stat_id = ObjectId(data_stat_id)
            else:
                self.bind_object.set_error('data_stat_id必须为ObjectId对象或其对应的字符串！', code="52800601")
        if not os.path.exists(stat_path):
            self.bind_object.set_error('stat_path所指定的路径不存在，请检查！', code="52800602")
        data_list = []  # 存入表格中的信息，然后用insert_many批量导入
        if data_type in ["clean", "optimised"]:
            mydb = self.db['data_stat_detail']
            raws = mydb.find({"data_stat_id": ObjectId(raw_data_stat_id)})
            self.specimen2id = name2id(raw_data_stat_id, type="raw")
            # raws= mydb.find({"data_stat_id": ObjectId("59cde2d3a4e1af116ca21c32")})
            if raws is None:
                self.bind_object.set_error('没有找到样品集数据', code="52800603")
            specimen = {}
            for raw in raws:
                specimen[raw['specimen_name']] = {"raw_read_num": raw['raw_read_num']}
                specimen[raw['specimen_name']]["raw_base"] = raw['raw_base']
        with open(stat_path, 'rb') as f:
            lines = f.readlines()
            for line in lines[1:]:  # 从第二行记录信息，因为第一行通常是表头文件，忽略掉
                line = line.strip().split('\t')
                data = {
                    "data_stat_id": data_stat_id,
                    #"specimen_name": line[0],
                }
                if data_type == "raw":
                    data['specimen_name'] = line[0]
                    data['specimen_source'] = line[1]
                    data['insert_size'] = int(line[2])
                    data['raw_read_len'] = int(line[3])
                    data['raw_read_num'] = int(line[4])
                    data['raw_base'] = int(line[5])
                    data['new_name'] = line[0]  #zouguanqing 20181126 宏基因组升级
                    data['desc'] = '--'      #zouguanqing 20181126 宏基因组升级
                if data_type == "unigene": # ysh 非冗余基因集上传
                    data['specimen_name'] = line[0]
                    data['raw_read_num'] = int(line[1])
                    data['raw_base'] = int(line[2])
                    data['raw_read_len'] = float(line[3])
                    data['insert_size'] = int(line[4])
                    data['new_name'] = line[0]  #zouguanqing 20181126 宏基因组升级
                    data['desc'] = '--'      #zouguanqing 20181126 宏基因组升级
                if data_type == "clean":
                    data['specimen_name'] = self.specimen2id[line[0]]
                    read_num = specimen[line[0]]["raw_read_num"]
                    base = specimen[line[0]]["raw_base"]
                    clean_ratio = float(line[1]) / read_num
                    clean_base_ratio = float(line[2]) / base
                    data['clean_read_num'] = int(line[1])
                    data['clean_base'] = int(line[2])
                    data['clean_ratio'] = clean_ratio
                    data['clean_base_ratio'] = clean_base_ratio
                if data_type == "optimised":
                    data['specimen_name'] = self.specimen2id[line[0]]
                    read_num = specimen[line[0]]["raw_read_num"]
                    base = specimen[line[0]]["raw_base"]
                    opt_ratio = float(line[1]) / read_num
                    opt_base_ratio = float(line[2]) / base
                    data['opt_read_num'] = int(line[1])
                    data['opt_base'] = int(line[2])
                    data['opt_ratio'] = opt_ratio
                    data['opt_base_ratio'] = opt_base_ratio
                data_list.append(data)
        try:
            collection = self.db['data_stat_detail']
            insert_ids = collection.insert_many(data_list).inserted_ids  # 用insert_many批量导入数据库，insert_one一次只能导入一条记录
            if data_type == 'raw' or data_type == 'unigene':  # add by guhaidong @ 20171204 对raw类型的detail表进行更新，增加origin_id字段
                for insert_id in insert_ids:
                    collection.update({'data_stat_id': data_stat_id, '_id': insert_id}, {'$set': {'origin_id': insert_id}})
            # self.update_raw_detail(data_stat_id, data_type)  # add by guhaidong @ 20171204
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (stat_path, e))
            self.bind_object.set_error("导入stat_path信息出错", code="52800604")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % stat_path)

    '''
    def update_raw_detail(self, data_stat_id, data_type):
        """
        added by guhaidong @ 20171204 用来增加raw表中的origin_id字段
        :param data_stat_id: 主表_id
        :param data_type: detail表类型, 只用为'raw'时执行此程序
        :return:
        """
        if data_type != 'raw':
            return
        collection = self.db['data_stat_detail']
        result = collection.find({'data_stat_id': data_stat_id})
        for one in result:
            origin_id = one['_id']
            collection.update({'_id': origin_id}, {'$set': {'origin_id': origin_id}})
    '''

    @report_check
    def new_add_specimen_graphic(self, file_path, task_id):
        df = pd.read_csv(file_path, sep='\t', chunksize=50000)
        collection = self.db['specimen_graphic']
        for one in df:
            one["task_id"] = task_id
            insert_data = one.to_dict('records')
            try:
                collection.insert_many(insert_data)
            except Exception, e:
                self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
                self.bind_object.set_error("导入specimen_graph信息出错", code="52800606")
        self.bind_object.logger.info("导入%s信息成功!" % file_path)

    @report_check
    def add_specimen_graphic(self, dir_path, task_id):
        self.specimen2id = name2id(task_id, type="task")
        filelist = os.listdir(dir_path.rstrip('/'))
        for i in filelist:
            line = i.strip().split('.', 1)
            specimen = self.specimen2id[line[0]]
            if re.search(r'1\.(fastq|fq)\.', line[1]):
                reads_direct = "left"
            elif re.search(r'2\.(fastq|fq)\.', line[1]):
                reads_direct = "right"
            elif re.search(r's\.(fastq|fq)\.', line[1]): # uniuque gene
                continue
            else:
                self.bind_object.logger.error(line[1])
                self.bind_object.set_error("没有匹配到fastq序列", code="52800605")
            data_list = []
            specimen_graphic = dir_path.rstrip('/') + '/' + i
            with open(specimen_graphic, 'rb') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.strip().split('\t')
                    data = {
                        "task_id": task_id,
                        "specimen_name": specimen,
                        "type": reads_direct,
                        "column": line[0],
                        "min": line[10],
                        "max": line[11],
                        "q1": line[6],
                        "median": line[7],
                        "q3": line[8],
                        "A": line[12],
                        "C": line[13],
                        "G": line[14],
                        "T": line[15],
                        "N": line[16],
                    }
                    data_list.append(data)
            try:
                collection = self.db['specimen_graphic']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.error("导入%s信息出错:%s" % (specimen_graphic, e))
                self.bind_object.set_error("导入specimen_graph信息出错", code="52800606")
            else:
                self.bind_object.logger.info("导入%s信息成功!" % specimen_graphic)

    @report_check
    def add_data_stat_detail_simple(self, stat_path, data_stat_id, data_type, raw_data_stat_id=None):
        if not isinstance(data_stat_id, ObjectId):
            if isinstance(data_stat_id, types.StringTypes):
                data_stat_id = ObjectId(data_stat_id)
            else:
                self.bind_object.set_error('data_stat_id必须为ObjectId对象或其对应的字符串！', code="52800601")
        self.bind_object.logger.info(stat_path)
        if not os.path.exists(stat_path):
            self.bind_object.set_error('stat_path所指定的路径不存在，请检查！', code="52800602")
        data_list = []  # 存入表格中的信息，然后用insert_many批量导入
        with open(stat_path, 'rb') as f:
            for line in f:
                line = line.strip().split('\t')
                data = {
                    "data_stat_id": data_stat_id,
                    #"specimen_name": line[0],
                }
                if data_type == "raw":
                    data['specimen_name'] = line[0]
                    data['insert_size'] = int(line[1])
                    data['new_name'] = line[0]  #zouguanqing 20181126 宏基因组升级
                    data['desc'] = '--'      #zouguanqing 20181126 宏基因组升级
                data_list.append(data)
        try:
            collection = self.db['data_stat_detail']
            insert_ids = collection.insert_many(data_list).inserted_ids  # 用insert_many批量导入数据库，insert_one一次只能导入一条记录
            if data_type == 'raw':  # add by guhaidong @ 20171204 对raw类型的detail表进行更新，增加origin_id字段
                for insert_id in insert_ids:
                    collection.update({'data_stat_id': data_stat_id, '_id': insert_id}, {'$set': {'origin_id': insert_id}})
            # self.update_raw_detail(data_stat_id, data_type)  # add by guhaidong @ 20171204
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (stat_path, e))
            self.bind_object.set_error("导入stat_path信息出错", code="52800604")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % stat_path)
