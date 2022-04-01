# -*- coding: utf-8 -*-
#__author__ = "zhangqingchen"

from biocluster.api.database.base import Base,report_check
from bson import ObjectId
from bson.son import SON
from mbio.packages.metagenomic.id_convert import name2id
import datetime
import json
import os
import types


class MgAnnoPfam(Base):
    def __init__(self, bind_object):
        super(MgAnnoPfam, self).__init__(bind_object)
        self.project_type = "metagenomic"

    @report_check
    def add_anno_pfam(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None,
            name=None, params=None, group_id=None, group_detail=None, is_origin=None):
        """
        pfam注释分析的结果主表
        """
        if not isinstance(geneset_id, ObjectId):
            if isinstance(geneset_id, types.StringTypes):
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:
                self.bind_object.set_error('geneset_id必须为ObjecId对象或者其对应的字符串！', code="52803601")
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
                            self.bind_object.set_error('geneset_id必须为字符串类型', code="52803602")
                    if group_detail == None:
                        self.bind_object.set_error('传入geneset_id时必须输入group_detail！', code="52803603")
                params = {
                    'database': 'pfam',
                    'group_detail': group_detail,
                    'group_id': group_id,
                    'geneset_id': geneset_id_str,
                    'identity': 0,
                    'align_length': 0,
                    'submit_location': 'annopfam',
                    'task_typs': 2
                }
                print params
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "desc": "pfam注释",
                "created_ts": created_ts,
                "status": "end",
                "name": name if name else "PFAM_Origin",
                "params": params,
                "geneset_id": geneset_id,
                "specimen": specimen,
                "anno_file": anno_file_path,
                "lowest_level": "Pfam ID",
                "is_origin": is_origin,
                "settled_params": json.dumps({"version": "pfam_v33.1"})
            }
            try:
                collection = self.db["anno_pfam"]
                anno_pfam_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入的anno_pfam主表异常：{}'.format(e))
                self.bind_object.set_error("导入的anno_pfam主表异常", code="52803604")
        else:
            if main_table_id is None:
                self.bind_object.set_error('main为FALSE时需要提供main_table_id！', coede="", code="52803605")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52803606")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_pfam_id = main_table_id
            except Exception, e:
                self.bind_object.logger.error('更新anno_pfam主表anno_file出错：{}'.format(e))
                self.bind_object.set_error('更新anno_pfam主表anno_file出错', code="52803607")
        return anno_pfam_id

    @report_check
    def add_anno_pfam_stat(self, anno_pfam_id, pfam_profile_dir):
        pfam_profile = pfam_profile_dir + "/gene_pfam_anno_stat.xls"
        if not isinstance(anno_pfam_id, ObjectId):
            if isinstance(anno_pfam_id, types.StringTypes):
                anno_pfam_id= ObjectId(anno_pfam_id)
            else:
                self.bind_object.set_error("anno_pfam_id必须为ObjectId对象或其对应的字符串！", code="52803608")
        if not os.path.exists(pfam_profile):
            self.bind_object.set_error('pfam_profile所指定的路径不存在，请检查！', code="52803609")
        data_list = list()
        result = self.db['anno_pfam'].find_one({'_id': ObjectId(anno_pfam_id)})
        if not result:
            self.bind_object.set_error('找不到pfam_anno_stat对应的主表id', code="52803610")
        else:
            pass
        with open(pfam_profile, "r") as f:
            head = f.next()
            if "#" in head:
                for line in f:
                    line = line.strip().split("\t")
                    if len(line) == 5:
                        data = [
                            ("pfam_id", anno_pfam_id),
                            ("pfam", line[0]),
                            ("pfam_name", line[1]),
                            ("pfam_desc", line[2]),
                            ("type", line[3]),
                            ("clan", line[4])
                        ]
                    if len(line) == 4:
                        data = [
                            ("pfam_id", anno_pfam_id),
                            ("pfam", line[0]),
                            ("pfam_name", line[1]),
                            ("pfam_desc", line[2]),
                            ("type", line[3]),
                        ]
                    data_son =SON(data)
                    data_list.append(data_son)
            else:
                self.bind_object.set_error('pfam_anno_stat.xls不是以#开头', code="52803611")
            try:
                collection = self.db["anno_pfam_stat"]
                collection.insert_many(data_list)
            except Exception,e:
                self.bind_object.logger.error("导入表格%s信息出错：%s" %(pfam_profile, e))
                self.bind_object.set_error("导入pfam_profile信息出错", code="52803612")
            else:
                self.bind_object.logger.info('导入表格信息成功: %s' % pfam_profile)

    @report_check
    def add_anno_pfam_detail(self, anno_pfam_id, pfam_profile_dir, ana_type, update_main=True):
        if not isinstance(anno_pfam_id, ObjectId):
            if isinstance(anno_pfam_id, types.StringTypes):
                anno_pfam_id = ObjectId(anno_pfam_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803613")
        if ana_type == "pfam":
            pfam_profile = pfam_profile_dir + "/pfam_acc_profile.xls"
        elif ana_type == "type":
            pfam_profile = pfam_profile_dir + "/pfam_type_profile.xls"
        elif ana_type == "clan":
            pfam_profile = pfam_profile_dir + "/pfam_clan_profile.xls"
        else:
            self.bind_object.set_error("缺少type字段！", code="52803614")
        if not os.path.exists(pfam_profile):
            self.bind_object.set_error('pfam_profile所指定的路径不存在，请检查', code="52803615")
        result = self.db['anno_pfam'].find_one({'_id': ObjectId(anno_pfam_id)})
        if not result:
            self.bind_object.set_error('找不到pfam_profile对应的主表id', code="52803616")
        else:
            task_id = result['task_id']
            #task_id = "tsg_32241"
            samples_dic = name2id(task_id, type="task")
        data_list = list()
        with open(pfam_profile, "r") as f:
            head = f.next()
            if "#" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads)-1]
            else:
                self.bind_object.set_error('生成的pfam_profile文件错误', code="52803617")
            for line in f:
                line = line.strip().split("\t")
                family_id = line[0]
                data = {
                    "family_id": family_id,
                    "ana_type": ana_type,
                    "pfam_id": anno_pfam_id,
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    data[sample_id] = float(line[i + 1])
                data = SON(data)
                data_list.append(data)
            try:
                collection = self.db["anno_pfam_detail"]
                collection.insert_many(data_list)
                self.db['anno_pfam'].update_one({'_id': anno_pfam_id}, {'$set':{'main_id':anno_pfam_id,
                                                                                "settled_params": json.dumps({"version": "pfam_v33.1"})}})
            except Exception, e:
                self.bind_object.logger.error("导入表格%s信息出错：%s" % (pfam_profile, e))
                self.bind_object.set_error("导入pfam_profile信息出错", code="52803618")
            else:
                self.bind_object.logger.info("导入anno_pfam_detail成功: %s" % pfam_profile)
            if update_main:
                main_table = self.db['anno_pfam']
                if "Total" in sams:
                    sams.remove('Total')
                specimen = ','.join([samples_dic[i] for i in sams])
                main_table.update_one({'_id': ObjectId(anno_pfam_id)}, {'$set': {'specimen': specimen}})

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_pfam'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})
