# -*- coding: utf-8 -*-
#__author__ = "zhangqingchen"

from biocluster.api.database.base import Base,report_check
from bson.objectid import ObjectId
from bson.son import SON
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name
import datetime
import json
import os
import types

class MgAnnoCyps(Base):
    def __init__(self, bind_object):
        super(MgAnnoCyps, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_anno_cyps(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None,
                      name=None, params=None, group_id=None, group_detail=None,is_origin=None):
        """
        cyps注释分析的结果主表
        params根据接口设计调整
        geneset_id同其他主表
        anno_file_path路径同其他主表
        """
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:
                self.bind_object.set_error('geneset_id必须为ObjectId对象或者其对应的字符串！', code="52804501")
        else:
            geneset_id_str = str(geneset_id)
        if main:
            task_id = self.bind_object.sheet.id
            #task_id = "tsg_32241"
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            if params == None:
                if group_id != None:
                    if not isinstance(group_id, types.StringTypes):
                        if isinstance(group_id, ObjectId):
                            group_id = str(group_id)
                        else:
                            self.bind_object.set_error('geneset_id必须为字符串类型', code="52804502")
                    if group_detail == None:
                        self.bind_object.set_error('传入geneset_id时必须输入group_detail！', code="52804503")
                params = {
                    'database': 'cyps',
                    'group_detail': group_detail,
                    'group_id': group_id,
                    'geneset_id': geneset_id_str,
                    'identity': 0,
                    'align_length': 0,
                    'submit_location': 'annocyps',
                    'task_typs': 2
                }
                print params
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "desc": "p450注释",
                "created_ts": created_ts,
                "status": "end",
                "name": name if name else "CYPS_Origin",
                "params": params,
                "geneset_id": geneset_id,
                "specimen": specimen,
                "anno_file": anno_file_path,
                "lowest_level": "Sid",
                "is_origin": is_origin 
            }
            try:
                collection = self.db["anno_cyps"]
                anno_cyps_id = collection.insert_one(insert_data).inserted_id
            except Exception,e:
                self.bind_object.logger.error('导入的anno_cyps主表异常：{}'.format(e))
                self.bind_object.set_error("导入anno_cyps主表异常", code="52804504")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为FALSE时需要提供main_table_id！", code="52804505")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52804506")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_cyps_id = main_table_id
            except Exception, e:
                self.bind_object.logger.error('更新anno_cyps主表anno_file出错：{}'.format(e))
                self.bind_object.set_error("更新anno_cyps主表anno_file出错", code="52804507")
        return anno_cyps_id

    @report_check
    def add_anno_cyps_detail(self, anno_cyps_id, cyps_profile_dir):
        cyps_profile = cyps_profile_dir + "/cyps_anno_stat.xls"
        if not isinstance(anno_cyps_id, ObjectId):
            if isinstance(anno_cyps_id, types.StringTypes):
                anno_cyps_id = ObjectId(anno_cyps_id)
            else:
                self.bind_object.set_error("anno_cyps_id必须为ObjectId对象或其对应的字符串!", code="52804508")
        if not os.path.exists(cyps_profile):
            self.bind_object.set_error('cyps_profile所指定的路径不存在，请检查！', code="52804509")
        data_list = list()
        result = self.db['anno_cyps'].find_one({'_id': ObjectId(anno_cyps_id)})
        if not result:
            self.bind_object.set_error('找不到cyps_anno_stat对应的主表id', code="52804510")
        else:
            pass
        with open(cyps_profile, "r") as f:
            head = f.next()
            if "#" in head:
                for line in f:
                    line = line.strip().split("\t")
                    data = [
                        ("cyps_id", anno_cyps_id),
                        ("sfam", line[1]),
                        ("hfam", line[2]),
                        ("nr", line[0]),
                        ("species", line[3]),
                        ("homo", line[4]),
                        ("super", line[5])
                    ]
                    data_son =SON(data)
                    data_list.append(data_son)
            else:
                self.bind_object.set_error('cyps_anno_stat.xls不是以#开头', code="52804511")
            try:
                collection = self.db["anno_cyps_detail"]
                collection.insert_many(data_list)
            except Exception,e:
                self.bind_object.logger.error("导入表格%s信息出错：%s" % (cyps_profile, e))
                self.bind_object.set_error("导入cyps_profile信息出错", code="52804512")
            else:
                self.bind_object.logger.info("导入表格%s信息成功！" % cyps_profile)

    @report_check
    def add_anno_cyps_abu(self, anno_cyps_id, cyps_profile_dir, type, update_main=True):
        if not isinstance(anno_cyps_id, ObjectId):
            if isinstance(anno_cyps_id, types.StringTypes):
                anno_cyps_id = ObjectId(anno_cyps_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串!", code="52804513")
        if type == "homo":
            cyps_profile = cyps_profile_dir + "/cyps_homo_family_profile.xls"
        elif type == "super":
            cyps_profile = cyps_profile_dir + "/cyps_super_family_profile.xls"
        else:
            self.bind_object.set_error("缺少typs字段！", code="52804514")
        if not os.path.exists(cyps_profile):
            self.bind_object.set_error('cyps_profile所指定的路径不存在，请检查！', code="52804515")
        result = self.db['anno_cyps'].find_one({'_id': ObjectId(anno_cyps_id)})
        if not result:
            self.bind_object.set_error('找不到cyps_profile对应的主表id', code="52804516")
        else:
            task_id = result['task_id']
            #task_id = "tsg_32241"
            self.bind_object.logger.info("从主表中取出的task_id:%s" % task_id)
            samples_dic = name2id(task_id, type="task")
        data_list = []
        with open(cyps_profile, "r") as f:
            head = f.next()
            if "#" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads)-1]
            else:
                self.bind_object.set_error('生成的cyps_profile文件错误！', code="52804517")
            for line in f:
                line = line.strip().split("\t")
                family = line[0]
                data = {
                    "cyps_id": anno_cyps_id,
                    "family": family,
                    "type": type,
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    #self.bind_object.logger.info("打印样本id：%s" % sample_id)
                    data[sample_id] = float(line[i + 1])
                data = SON(data)
                data_list.append(data)
            try:
                collection = self.db["anno_cyps_abu"]
                collection.insert_many(data_list)
                self.db['anno_cyps'].update_one({'_id': anno_cyps_id}, {'$set':{'main_id':anno_cyps_id}})
            except Exception,e:
                main_table = self.db['anno_cyps']
                if "Total" in sams:
                    sams.remove('Total')
                specimen = ','.join([samples_dic[i] for i in sams])
                main_table.update_one({'_id': ObjectId(anno_cyps_id)}, {'$set': {'specimen':specimen}})

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_cyps'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
