# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20171114
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

class MgAnnoNr(Base):
    def __init__(self, bind_object):
        super(MgAnnoNr, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_anno_nr(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None, name=None, params=None,
                    group_id=None, group_detail=None, task_id=None, software_ver={}):
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_nr_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52802201")
        else:
            geneset_id_str = str(geneset_id)
        #if not os.path.exists(anno_file_path): # 调整为上传的永久路径，由于先导表，暂不检查
            #raise Exception('anno_file_path所指定的路径不存在，请检查！')
        if main:
            if not task_id:
                task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            if params == None:
                if group_id != None:
                    if not isinstance(group_id, types.StringTypes):
                        if isinstance(group_id, ObjectId):
                            group_id = str(group_id)
                        else:
                            self.bind_object.set_error('geneset_id必须为字符串类型！', code="52802202")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52802203")
                if name and name == 'NR_Origin_LCA': # add by qingchen.zhang@20190924工作流用于增加nr_method字段
                    nr_method = 'lca'
                elif name and name == 'NR_Origin_Deunclassified':
                    nr_method = 'de_unclassified'
                else:
                    nr_method = 'best_hit'
                params = {
                    "database": "nr",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "identity": 0,
                    "nr_method": nr_method, # add by qingchen.zhang@20190924工作流用于增加nr_method字段
                    "align_length": 0,
                    "submit_location": "annonr",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': "",
                'created_ts': created_ts,
                'name': name if name else "NR_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                "settled_params": json.dumps({"version": "nr_v20200604"})
            }
            if software_ver:
                insert_data.update(software_ver)
            try:
                collection = self.db['anno_nr']
                anno_nr_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_nr主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_nr异常", code="52802204")
        else:
            if main_table_id is None:
                self.bind_object.logger.error("main为False时需提供main_table_id!", code="52802205")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52802206")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_nr_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_nr主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_nr主表anno_file_path出错", code="52802207")
        return anno_nr_id

    @report_check
    def add_anno_nr_detail(self, anno_nr_id, nr_profile_dir, update_main=True):
        self.bind_object.logger.info("start nr detail")
        if not isinstance(anno_nr_id, ObjectId):  # 检查传入的anno_nr_id是否符合ObjectId类型
            if isinstance(anno_nr_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                anno_nr_id = ObjectId(anno_nr_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_nr_id必须为ObjectId对象或其对应的字符串！', code="52802208")
        if not os.path.exists(nr_profile_dir):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('nr_profile_dir所指定的路径不存在，请检查！', code="52802209")
        # data_list = list()  # 存入表格中的信息，然后用insert_many批量导入
        result = self.db['anno_nr'].find_one({'_id': anno_nr_id})
        if not result:
            self.bind_object.set_error('找不到anno_nr_detail对应的主表id', code="52802210")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        f_d = os.path.join(nr_profile_dir, "tax_d.xls")
        f_k = os.path.join(nr_profile_dir, "tax_k.xls")
        f_p = os.path.join(nr_profile_dir, "tax_p.xls")
        f_c = os.path.join(nr_profile_dir, "tax_c.xls")
        f_o = os.path.join(nr_profile_dir, "tax_o.xls")
        f_f = os.path.join(nr_profile_dir, "tax_f.xls")
        f_g = os.path.join(nr_profile_dir, "tax_g.xls")
        f_s = os.path.join(nr_profile_dir, "tax_s.xls")
        file_list = [f_d, f_k, f_p, f_c, f_o, f_f, f_g, f_s]
        level_list = [1, 2, 3, 4, 5, 6, 7, 8]
        i = 0
        d_id = k_id = p_id = c_id = o_id = f_id = g_id = s_id = 0
        name = {}
        for each in file_list:
            self.bind_object.logger.info("start {} profile".format(each))
            if not os.path.exists(each):
                self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(each), code="52802211")
            else:
                i += 1
                level_id = level_list[i - 1]
                data_list = []
                with open(each, 'rb') as f:
                    head = f.next()
                    if "#Taxonomy" in head:
                        sams = head.strip().split("\t")[1:len(head)]
                    else:
                        self.bind_object.set_error("丰度文件错误或者丰度文件head错误!", code="52802212")
                    for line in f:
                        line = line.strip().split('\t')
                        tax = line[0].split(";")
                        insert_data = {
                            'level_id': level_id,
                            'nr_id': anno_nr_id,
                        }
                        for j in range(0, len(sams)):
                            if not sams[j] == "Total":
                                sample_id = samples_dic[sams[j]]
                            else:
                                sample_id = sams[j]
                            insert_data[sample_id] = float(line[j + 1])
                        if len(tax) >= 1:
                            insert_data['d__'] = tax[0]
                            if len(tax) == 8:
                                if not name.has_key(tax[0]):
                                    d_id += 1
                                    name[tax[0]] = "d" + str(d_id)
                                insert_data['d_id'] = name[tax[0]]
                        if len(tax) >= 2:
                            insert_data['k__'] = tax[1]
                            if len(tax) == 8:
                                if not name.has_key(tax[1]):
                                    k_id += 1
                                    name[tax[1]] = "k" + str(k_id)
                                insert_data['k_id'] = name[tax[1]]
                        if len(tax) >= 3:
                            insert_data['p__'] = tax[2]
                            if len(tax) == 8:
                                if not name.has_key(tax[2]):
                                    p_id += 1
                                    name[tax[2]] = "p" + str(p_id)
                                insert_data['p_id'] = name[tax[2]]
                        if len(tax) >= 4:
                            insert_data['c__'] = tax[3]
                            if len(tax) == 8:
                                if not name.has_key(tax[3]):
                                    c_id += 1
                                    name[tax[3]] = "c" + str(c_id)
                                insert_data['c_id'] = name[tax[3]]
                        if len(tax) >= 5:
                            insert_data['o__'] = tax[4]
                            if len(tax) == 8:
                                if not name.has_key(tax[4]):
                                    o_id += 1
                                    name[tax[4]] = "o" + str(o_id)
                                insert_data['o_id'] = name[tax[4]]
                        if len(tax) >= 6:
                            insert_data['f__'] = tax[5]
                            if len(tax) == 8:
                                if not name.has_key(tax[5]):
                                    f_id += 1
                                    name[tax[5]] = "f" + str(f_id)
                                insert_data['f_id'] = name[tax[5]]
                        if len(tax) >= 7:
                            insert_data['g__'] = tax[6]
                            if len(tax) == 8:
                                if not name.has_key(tax[6]):
                                    g_id += 1
                                    name[tax[6]] = "g" + str(g_id)
                                insert_data['g_id'] = name[tax[6]]
                        if len(tax) == 8:
                            insert_data['s__'] = tax[7]
                            if len(tax) == 8:
                                if not name.has_key(tax[7]):
                                    s_id += 1
                                    name[tax[7]] = "s" + str(s_id)
                                insert_data['s_id'] = name[tax[7]]
                        # data = SON(data)
                        data_list.append(insert_data)
                    try:
                        collection = self.db['anno_nr_detail']
                        collection.insert_many(data_list)
                        self.db['anno_nr'].update_one({'_id':anno_nr_id},{'$set':{'main_id':anno_nr_id}})
                    except Exception as e:
                        self.bind_object.logger.error("导入tax_profile表格%s信息出错:%s" % (each, e))
                        self.bind_object.set_error("导入tax_profile表格信息出错", code="52802213")
                    else:
                        self.bind_object.logger.info("导入tax_profile表格%s信息成功!" % each)
        if update_main:
            main_table = self.db['anno_nr']
            if "Total" in sams:
                    sams.remove("Total")
            specimen = ",".join([samples_dic[i] for i in sams])
            # specimen = ",".join(sams)
            main_table.update_one({'_id': ObjectId(anno_nr_id)}, {'$set': {'specimen': specimen}})
        main_table = self.db['anno_nr']
        main_table.update_one({'_id': ObjectId(anno_nr_id)}, {'$set': {'status': 'end',
                                                                       "settled_params": json.dumps({"version": "nr_v20200604"})}})  #important


    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_nr'].update_one({'_id': ObjectId(main_table_id)}, \
                                      {'$set': {'anno_file': anno_file}})
