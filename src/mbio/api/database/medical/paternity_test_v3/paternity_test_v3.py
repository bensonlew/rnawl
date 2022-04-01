# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
import re
from biocluster.api.database.base import Base, report_check
from bson import ObjectId
from types import StringTypes
import datetime
import phpserialize
from dateutil.relativedelta import relativedelta
import requests
import hashlib
import json
import time
import os


class PaternityTestV3(Base):
    """
    亲子V3版的各种导表
    """
    def __init__(self, bind_object):
        super(PaternityTestV3, self).__init__(bind_object)
        self._project_type = 'pt_v3'
        self.err_list = [5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        self.err_list_homohybrid = [5, 6, 7, 8, 9, 10, 15, 20]

    def get_case_sampls(self, case_name):
        """
        根据case号到sg_sample_info中查找case中的所有样本名
        ('sg_sample_info').find({"sample_id": {"$regex": /^WS(.*)-(\d)-(M|F|S)(.*)/}})
        """
        result_list = []
        if re.match("(WQ.*)-(\d)$", case_name):
            result = self.db['sg_sample_info'].find({'case_name': case_name.split("-")[0]})
            result1 = self.db['sg_sample_info'].find({'case_name': case_name})
            result_list.extend(list(result))
            result_list.extend(list(result1))
        else:
            result = self.db['sg_sample_info'].find({'case_name': case_name})
            result1 = self.db['sg_sample_info'].find({'case_name': {"$regex": "^{}-(\d)$".format(case_name)}})
            result_list.extend(list(result))
            result_list.extend(list(result1))
        return result_list

    def make_familys(self, case_name):
        """
        更根据case name 找到这个家系中所有的父本与母本的组合
        :param case_name:
        :return:
        """
        family_list = []
        dad_id = []
        mom_id = []
        result = self.get_case_sampls(case_name)
        if result:
            for m in result:
                if "-M" in m['sample_id']:
                    mom_id.append(m['sample_id'])
                elif "-F" in m['sample_id']:
                    dad_id.append(m['sample_id'])
                else:
                    pass
            if len(list(set(dad_id))) > 0 and len(list(set(mom_id))) > 0:
                for dad in dad_id:
                    for mom in mom_id:
                        family_list.append([dad, mom])
        self.bind_object.logger.info("重组后的家系列表为：{}".format(family_list))
        return family_list

    def update_datasplit_html(self, path, batch_id):
        """
        在sg_datasplit表中添加split_result_html路径，用于页面点击查看
        add by yuguo @20171214
        """
        batch_id = ObjectId(batch_id)
        collection = self.db['sg_datasplit']
        try:
            collection.update({'_id': batch_id}, {"$set": {'split_result_html': path}}, upsert=True)
        except Exception, e:
            raise Exception('更新sg_datasplit中的htmlresult出错：{}'.format(e))
        else:
            self.bind_object.logger.info('更新sg_datasplit中的htmlresult成功！')

    def update_datasplit_status(self, batch_id):
        """
        更新拆分的状态为end
        add by yuguo @20171214
        """
        batch_id = ObjectId(batch_id)
        collection = self.db['sg_datasplit']
        try:
            collection.update({'_id': batch_id}, {"$set": {'status': "end", 'desc': "运行结束"}})
        except Exception, e:
            self.bind_object.logger.info('更新sg_datasplit中的status失败！{}'.format(e))  # 这里错误不用终止
        else:
            self.bind_object.logger.info('更新sg_datasplit中的status成功！')

    def add_sg_father_err(self, batch_id, family_id, father_id, err_min, types, name, member_id, postype='hybrid'):
        """
        用于添加不同错配的father主表 add by hongdong @20171121
        :param batch_id: 拆分表的主表id
        :param family_id: 家系的id
        :param father_id:
        :param err_min:
        :param types: 为1代表是家系分析的主表，不然就是自由交互的主表
        :param name: 这个家系的父本+母本+胎儿+error
        :param member_id: 会员id
        :param postype: 位点类型,默认杂合（hybrid），homohybrid
        :return:
        """
        if types == '1':
            if not isinstance(batch_id, ObjectId):
                if isinstance(batch_id, StringTypes):
                    batch_id = ObjectId(batch_id)
                else:
                    raise Exception("batch_id必须为ObjectId对象或其对应的字符串!")
            if not isinstance(family_id, ObjectId):
                if isinstance(family_id, StringTypes):
                    family_id = ObjectId(family_id)
                else:
                    raise Exception("family_id必须为ObjectId对象或其对应的字符串!")
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        insert_data = {
            "batch_id": batch_id,
            "family_id": family_id,
            "father_id": father_id,
            "err_min": err_min,
            "types": types,
            "name": name,
            "member_id": member_id,
            "status": "start",
            "is_filter": "no",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        try:
            if postype == 'homohybrid':
                collection = self.db["sg_father_err_homohybrid"]
            else:
                collection = self.db["sg_father_err"]
            father_err_id = collection.insert_one(insert_data).inserted_id
        except Exception as e:
            raise Exception('导入sg_father_err主表{}出错：{}-{}'.format(err_min, e, postype))
        else:
            self.bind_object.logger.info("导入sg_father_err主表{}成功--{}".format(err_min, postype))
        return father_err_id

    def update_father_err_status(self, father_err_id, postype='hybrid'):
        """
        用于更新father_err主表的状态为end, add by hongdong@20171121
        :param father_err_id:
        :param postype:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        try:
            if postype == 'homohybrid':
                self.db['sg_father_err_homohybrid'].update({"_id": father_err_id}, {"$set": {"status": "end"}},
                                                           multi=True)
            else:
                self.db['sg_father_err'].update({"_id": father_err_id}, {"$set": {"status": "end"}}, multi=True)
        except Exception as e:
            raise Exception('更新sg_father_err主表状态为end出错：{}-{}'.format(e, postype))
        else:
            self.bind_object.logger.info("更新sg_father_err主表状态为end成功")

    def add_sg_pt_father_detail(self, file_path, father_err_id):
        """
        用于添加调试的表格, add by hongdong@20171121
        :param file_path:
        :param father_err_id:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                self.bind_object.logger.info({"father_err_id:{}".format(father_err_id)})
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "chrom":
                    continue
                if line[8] == 'NA':
                    dad_rf = 'NA'
                else:
                    dad_rf = str(round(float(line[8]), 8))    # 保留小数点8位然后转字符串
                if line[17] == 'NA':
                    preg_rf = 'NA'
                else:
                    preg_rf = str(round(float(line[17]), 8))
                if line[26] == 'NA':
                    mom_rf = 'NA'
                else:
                    mom_rf = str(round(float(line[26]), 8))

                insert_data = {
                    "father_err_id": father_err_id,
                    "chrom": line[0],
                    "pos": line[1],
                    "dad_id": line[2],
                    "dad_ref": line[3],
                    "dad_alt": line[4],
                    "dad_dp": line[5],
                    "dad_ref_dp": line[6],
                    "dad_alt_dp": line[7],
                    "dad_rf": dad_rf,
                    "dad_geno": line[9],
                    "dad_geno_bases": line[10],
                    "preg_id": line[11],
                    "preg_ref": line[12],
                    "preg_alt": line[13],
                    "preg_dp": line[14],
                    "preg_ref_dp": line[15],
                    "preg_alt_dp": line[16],
                    "preg_rf": preg_rf,
                    "preg_geno": line[18],
                    "preg_geno_bases": line[19],
                    "mom_id": line[20],
                    "mom_ref": line[21],
                    "mom_alt": line[22],
                    "mom_dp": line[23],
                    "mom_ref_dp": line[24],
                    "mom_alt_dp": line[25],
                    "mom_rf": mom_rf,
                    "mom_geno": line[27],
                    "mom_geno_bases": line[28],
                    "reg": line[29],
                    "from": line[30],
                    "to": line[31],
                    "rs": line[32],
                    "hapmap_rf": line[33],
                    "hapmap_geno": line[34],
                    "n": line[35],
                    "mj_ref": line[36],
                    "pA": line[37],
                    "pG": line[38],
                    "pC": line[39],
                    "pT": line[40],
                    "mj_dp": line[41],
                    "mj_gene": line[42],
                    "is_test": line[43],
                    "is_mis": line[44],
                    "mustbe": line[45],
                    "mustnotbe": line[46],
                    "good": line[47],
                    "pi": line[48]
                }
                sg_pt_family_detail.append(insert_data)
            try:
                collection = self.db['sg_father_err_debug']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入调试表格出错：{}'.format(e))
                raise Exception('导入调试表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入调试表格成功")

    def add_info_detail(self, file_path, father_err_id, dad_id=None, mom_id=None, son_id=None):
        """
        bed.preg.id     dp_preg percent error   s_signal        bed.mom.id      dp_mom  mom_preg
        WQ17114249-S-1  74.32   10.08   0.5     41.9512 WQ17114249-M-1  269.13  99.0148
        add by hongdong@20171122
        :param file_path:
        :param father_err_id:
        :param dad_id:
        :param mom_id:
        :param son_id:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                self.bind_object.logger.info({"father_err_id:{}".format(father_err_id)})
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        pregnant_week = '--'
        if dad_id and mom_id and son_id:
            if "-T" in dad_id:
                if dad_id.startswith("TEST"):
                    dad_id_ = re.match('(TESTWQ.*)-(T.*)', dad_id).group(1)
                else:
                    dad_id_ = re.match('(WQ.*)-(T.*)', dad_id).group(1)
            elif '-L' in dad_id:
                if dad_id.startswith("TEST"):
                    dad_id_ = re.match('(TESTWQ.*)-(L.*)', dad_id).group(1)
                else:
                    dad_id_ = re.match('(WQ.*)-(L.*)', dad_id).group(1)
            else:
                dad_id_ = dad_id
            if '-T' in mom_id:
                if mom_id.startswith("TEST"):
                    mom_id_ = re.match('(TESTWQ.*)-(T.*)', mom_id).group(1)
                else:
                    mom_id_ = re.match('(WQ.*)-(T.*)', mom_id).group(1)
            elif '-L' in mom_id:
                if mom_id.startswith("TEST"):
                    mom_id_ = re.match('(TESTWQ.*)-(L.*)', mom_id).group(1)
                else:
                    mom_id_ = re.match('(WQ.*)-(L.*)', mom_id).group(1)
            else:
                mom_id_ = mom_id
            son_id_new = son_id.split('_')[0]
            if '-T' in son_id_new:
                if son_id_new.startswith("TEST"):
                    son_id_ = re.match('(TESTWQ.*)-(T.*)', son_id_new).group(1)
                else:
                    son_id_ = re.match('(WQ.*)-(T.*)', son_id_new).group(1)
            elif '-L' in son_id_new:
                if son_id_new.startswith("TEST"):
                    son_id_ = re.match('(TESTWQ.*)-(L.*)', son_id_new).group(1)
                else:
                    son_id_ = re.match('(WQ.*)-(L.*)', son_id_new).group(1)
            else:
                son_id_ = son_id_new
            gest_info = self.gest_info_post(dad_id_, mom_id_, son_id_)
            if 'pregnant_week' in gest_info.keys():
                pregnant_week = gest_info['pregnant_week']
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "bed.preg.id":
                    continue
                result = 'no'
                if str(line[7]) != 'NA' and eval(line[7]) >= 95:
                    result = 'yes'
                insert_data = {
                    "father_err_id": father_err_id,
                    "preg_id": line[0],
                    "dp_preg": line[1],
                    "percent": line[2],
                    "error": line[3],
                    "s_signal": line[4],
                    "mom_id": line[5],
                    "dp_mom": line[6],
                    "mom_preg": line[7],
                    "result": result,
                    'pregnant_week': pregnant_week
                }
                sg_pt_family_detail.append(insert_data)
            try:
                collection = self.db['sg_father_err_ms_match']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入母本与胎儿匹配信息出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入母本与胎儿匹配信息成功")

    def import_dedup_data(self, file_path, father_dedup_id=None, father_err_id=None):
        """
        导入查重后的表格, 注这里少了样本的建库时间与抽提时间, add by hongdong@20171122
        :param file_path:
        :param father_dedup_id:
        :param father_err_id:有该字段的时候 就是自由交互中的查重表格
        :return:
        """
        if father_dedup_id != 'no':
            if not isinstance(father_dedup_id, ObjectId):
                if isinstance(father_dedup_id, StringTypes):
                    father_dedup_id = ObjectId(father_dedup_id)
                else:
                    raise Exception("father_dedup_id必须为ObjectId对象或其对应的字符串!")
        if father_err_id:
            if not isinstance(father_err_id, ObjectId):
                if isinstance(father_err_id, StringTypes):
                    father_err_id = ObjectId(father_err_id)
                else:
                    raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        dad_id = []
        sg_pt_family_detail = list()
        no_analyzed_sample_results = dict()

        no_analyzed_samples = set()
        no_analyzed = self.db['sg_father_dedup_detail'].find({"father_dedup_id": father_dedup_id, "result": '样本未上机'})
        for i in no_analyzed:
            if "dad_id" in i:
                no_analyzed_samples.add(i["dad_id"])

        with open(file_path, 'r') as f:
            data = f.readlines()[1:]
            for line in data:
                line = line.strip().split('\t')
                temp_fp = eval(line[4])
                rcp = float(temp_fp) / (float(temp_fp) + 1)
                if rcp > 0.5:
                    rcp_result = ">99.99%"
                else:
                    rcp_result = "<0.01%"
                if str(line[3]) == "NA":
                    err_rate = 'NA'
                else:
                    err_rate = float('%.2f' % float(line[3]))
                extract_batch, library_batch, deal_batch, receive_batch, sj_batch_no, zj_batch_no = \
                    self.find_father_batch(line[0])
                dad_time_ = '--'
                dad_time_result = self.db['sg_sample_qc'].find_one({"sample_id": line[0]})
                if dad_time_result and "date" in dad_time_result.keys():
                    dad_time_ = dad_time_result['date']
                insert_data = {
                    "father_dedup_id": father_dedup_id,
                    "dad_id": line[0],
                    "test_pos_n": line[1],
                    "err_pos_n": line[2],
                    "test_pos_n_int": int(line[1]),
                    "err_pos_n_int": int(line[2]),
                    "err_rate": err_rate,
                    "err_rate_str": line[3],
                    "fq": format(eval(line[4]), '.2e'),  # 科学计数法保留两位
                    "dp": line[5],
                    "eff_rate": line[6],
                    "ineff_rate": line[7],
                    "result": line[8],
                    "rcp": rcp_result,
                    "dad_time": dad_time_,
                    "extract_batch": extract_batch,
                    "library_batch": library_batch,
                    "deal_batch": deal_batch,
                    "receive_batch": receive_batch,
                    "sj_batch_no": sj_batch_no,
                    "zj_batch_no": zj_batch_no,
                    "eff_ineff_rate": str(float(line[6].split(":")[0].strip()[:-1]) - float(line[7].split(":")
                                                                                            [0].strip()[:-1])) + "%"
                }
                if father_err_id:
                    insert_data.update({"father_err_id": father_err_id})
                if line[0] in no_analyzed_samples:
                    no_analyzed_sample_results[line[0]] = insert_data
                else:
                    sg_pt_family_detail.append(insert_data)
                dad_id.append(line[0])
            collection = self.db['sg_father_dedup_detail']

            for sample in no_analyzed_sample_results:
                try:
                    collection.update({"father_dedup_id": father_dedup_id, "result": '样本未上机', "dad_id": sample},
                                      {"$set": no_analyzed_sample_results[sample]})
                except Exception as e:
                    self.bind_object.logger.error('更新批次中未上机样品出错：{}'.format(e))
                else:
                    self.bind_object.logger.info("更新批次中未上机样品成功".format())
            try:
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                self.bind_object.logger.error('导入查重表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入查重表格成功".format())

    def add_dedup_main(self, err_min, mom_id, son_id):
        """
        添加不同错配的查重 add by hongdong @20171123
        :param err_min:
        :param mom_id:
        :param son_id:
        :return:
        """
        father_dedup_id = ''
        insert_data = {
            'm_s_id': mom_id + '_' + son_id,
            'err_min': err_min,
            'create_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            # 'ready_dad_id': []
        }
        try:
            father_dedup_id = self.db['sg_father_dedup'].insert_one(insert_data).inserted_id
        except Exception as e:
            self.bind_object.logger.error('导入{}与{}查重{}主表出错：{}'.format(mom_id, son_id, err_min, e))
        else:
            self.bind_object.logger.info("导入{}与{}查重{}主表成功".format(mom_id, son_id, err_min))
        return father_dedup_id

    def find_dedup_err_id(self, err, mom_id, son_id):
        """
        用于查找同一母本+胎儿的主表，不同错配，id, add by hongdong@20171122
        :param err:
        :param mom_id:
        :param son_id:
        :return:
        """
        _id = ''
        m_s_id = mom_id + '_' + son_id
        try:
            result = self.db['sg_father_dedup'].find_one({'m_s_id': m_s_id, 'err_min': err})
            if result:
                _id = result['_id']
            else:
                return False
        except Exception as e:
            self.bind_object.logger.error('查找{}与{}的查重主表出错：{}'.format(mom_id, son_id, e))
        else:
            self.bind_object.logger.info("查找{}与{}查重主表成功".format(mom_id, son_id))
        return _id

    def update_dedup_father(self, err_min, mom_id, son_id, father_id):
        """
        用于更新样本的查重批次的一些信息，更新的字段主要在sg_father中，逻辑有点复杂，单独拎出来进行整理，整个逻辑是查找该
        父本的所有的批次信息，然后整合好，去查重库中进行比较，如果所有父本都在的时候，ok，显示全部已经查重，如果有部分不在的时候，、
        就要判断是否已经下机还是样本数据量低，如果没有下机的也要在sg_sample_qc中导入一份数据
        :param err_min:
        :param mom_id:
        :param son_id:
        :param father_id:
        :return:
        """
        if father_id:
            if not isinstance(father_id, ObjectId):
                if isinstance(father_id, StringTypes):
                    father_id = ObjectId(father_id)
                else:
                    raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        library_batch = '--'
        extract_batch = '--'
        deal_batch = "--"
        receive_batch = "--"
        sj_batch_no = "--"
        zj_batch_no = "--"

        dad_list = []  # 存储该父母本的查重批次中所有的其它父本
        dad_no_dedup = []  # 存储没有没有进行查重到的父本
        dad_dp_err = []  # 没有上机了，但是深度低
        dad_no_analysis = []  # 没有上机
        fa_collection = self.db['sg_father']
        fa_result = fa_collection.find_one({"_id": father_id})
        if fa_result:
            library_batch = fa_result['library_batch'] if 'library_batch' in fa_result and fa_result['library_batch'] \
                else "not used to retrieve samples"
            extract_batch = fa_result['extract_batch'] if 'extract_batch' in fa_result and fa_result['extract_batch'] \
                else "not used to retrieve samples"
            deal_batch = fa_result["deal_batch"] if "deal_batch" in fa_result and fa_result["deal_batch"] \
                else "not used to retrieve samples"
            receive_batch = fa_result["receive_batch"] if "receive_batch" in fa_result and fa_result["receive_batch"] \
                else "not used to retrieve samples"
            sj_batch_no = fa_result["sj_batch_no"] if "sj_batch_no" in fa_result and fa_result["sj_batch_no"] \
                else "not used to retrieve samples"
            zj_batch_no = fa_result["zj_batch_no"] if "zj_batch_no" in fa_result and fa_result["zj_batch_no"] \
                else "not used to retrieve samples"

        batch_result = self.db['sg_sample_batch'].find({"$or": [{"library_batch": library_batch},
                                                                {"extract_batch": extract_batch},
                                                                {"deal_batch": deal_batch},
                                                                {"receive_batch": receive_batch},
                                                                {"sj_batch_no": sj_batch_no},
                                                                {"zj_batch_no": zj_batch_no}]})
        if batch_result:
            for m in batch_result:
                if '-F' in str(m['sample_id']) or "-M" in str(m['sample_id']):
                    dad_list.append(m['sample_id'])
        dad_list = list(set(dad_list))                                              # 去重
        father_dedup_id = self.find_dedup_err_id(int(err_min), mom_id, son_id)
        if father_dedup_id:
            dedup_dad_list = self.find_dedup_father_id(father_dedup_id)
        else:
            dedup_dad_list = []
        for i in dad_list:
            if i not in set(dedup_dad_list):
                dad_no_dedup.append(i)
        if dad_no_dedup:
            for n in dad_no_dedup:
                if self.db['sg_sample_qc'].find_one({"sample_id": n, "is_problem": '1'}):
                    dad_dp_err.append(n)
                else:
                    dad_no_analysis.append(n)
            batch_dedup = '{}/{}'.format(len(dad_list) - len(dad_no_dedup), len(dad_list))
            problem_sample_num = str(len(dad_no_dedup))
            problem_sample_name = {"1": ','.join(dad_dp_err), '2': ','.join(dad_no_analysis)}  # 1代表的是深度低的样本，
            # 2代表的是没有上机的样本
        else:
            batch_dedup = '{}/{}'.format(len(dad_list), len(dad_list))
            problem_sample_num = '0'
            problem_sample_name = {'1': '', '2': ''}
        if int(err_min) == 5:                              # 将err_min为5的结果更新到sg_family中
            try:
                fa_collection.update({"_id": father_id}, {"$set": {"batch_dedup": batch_dedup,
                                                                   "problem_sample_num": problem_sample_num,
                                                                   "problem_sample_name": problem_sample_name}})
            except Exception as e:
                self.bind_object.logger.error("更新样本的排查批次以及问题样本信息到sg_father失败！{}".format(e))
        if not dad_no_dedup:
            return
        dp_err_sample_list = []
        if dad_dp_err:
            for pro_dad in dad_dp_err:
                extract_batch_, library_batch_, deal_batch_, receive_batch_, sj_batch_no_, zj_batch_no_ = \
                    self.find_father_batch(pro_dad)
                insert_data = {
                    "father_dedup_id": father_dedup_id,
                    "dad_id": pro_dad,
                    "test_pos_n": '',
                    "err_pos_n": '',
                    "err_rate": '',
                    "fq": '',
                    "dp": '',
                    "eff_rate": '',
                    "ineff_rate": '',
                    "result": '样本深度低',
                    "rcp": '',
                    "dad_time": '',
                    "extract_batch": extract_batch_,
                    "library_batch": library_batch_,
                    "deal_batch": deal_batch_,
                    "receive_batch": receive_batch_,
                    "sj_batch_no": sj_batch_no_,
                    "zj_batch_no": zj_batch_no_,
                    "eff_ineff_rate": ""
                }
                # 匹配到抽提批次或建库批次就会在页面上显示，如果有能显示的数据就不导表
                if self.db['sg_father_dedup_detail'].find_one({"father_dedup_id": father_dedup_id,
                                                               "dad_id": pro_dad,
                                                               "result": '样本深度低',
                                                               "$or": [{"extract_batch": extract_batch_},
                                                                       {"library_batch": library_batch_},
                                                                       {"deal_batch": deal_batch_},
                                                                       {"receive_batch": receive_batch_},
                                                                       {"sj_batch_no": sj_batch_no_},
                                                                       {"zj_batch_no": zj_batch_no_}]}):
                    continue
                dp_err_sample_list.append(insert_data)
        if dad_no_analysis:
            for no_ana_dad in dad_no_analysis:
                extract_batch_, library_batch_, deal_batch_, receive_batch_, sj_batch_no_, zj_batch_no_ = \
                    self.find_father_batch(no_ana_dad)
                insert_data = {
                    "father_dedup_id": father_dedup_id,
                    "dad_id": no_ana_dad,
                    "test_pos_n": '',
                    "err_pos_n": '',
                    "err_rate": '',
                    "fq": '',
                    "dp": '',
                    "eff_rate": '',
                    "ineff_rate": '',
                    "result": '样本未上机',
                    "rcp": '',
                    "dad_time": '',
                    "extract_batch": extract_batch_,
                    "library_batch": library_batch_,
                    "deal_batch": deal_batch_,
                    "receive_batch": receive_batch_,
                    "sj_batch_no": sj_batch_no_,
                    "zj_batch_no": zj_batch_no_,
                    "eff_ineff_rate": ""
                }
                if self.db['sg_father_dedup_detail'].find_one({"father_dedup_id": father_dedup_id, "dad_id": no_ana_dad,
                                                               "result": '样本未上机',
                                                               "$or": [{"extract_batch": extract_batch_},
                                                                       {"library_batch": library_batch_},
                                                                       {"deal_batch": deal_batch_},
                                                                       {"receive_batch": receive_batch_},
                                                                       {"sj_batch_no": sj_batch_no_},
                                                                       {"zj_batch_no": zj_batch_no_}]
                                                               }):
                    continue
                dp_err_sample_list.append(insert_data)
        if dp_err_sample_list:
            try:
                collection = self.db['sg_father_dedup_detail']
                collection.insert_many(dp_err_sample_list)
            except Exception as e:
                self.bind_object.logger.error('导入异常样本到查重表格出错：{}'.format(e))
            else:
                self.bind_object.logger.info("导入异常样本到查重表格成功")
        else:
            self.bind_object.logger.info("没有需要导入到查重详细表中的异常样本！")
        self.add_sample_pt(dad_no_analysis, 'no')  # 用于导入没有上机的样本数据到sg_sample_qc表中

    def find_father_batch(self, dad_id):
        """
        用于根据dad_id去获取到父本的建库批次与抽提批次
        :param dad_id:
        :return:extract_batch, library_batch
        """
        result = self.db['sg_sample_batch'].find_one({"sample_id": dad_id})
        if result:
            if "extract_batch" in result:
                extract_batch = result['extract_batch']
            else:
                extract_batch = '--'
            if "library_batch" in result:
                library_batch = result['library_batch']
            else:
                library_batch = '--'
            if "deal_batch" in result:
                deal_batch = result['deal_batch']
            else:
                deal_batch = '--'
            if 'receive_batch' in result:
                receive_batch = result["receive_batch"]
            else:
                receive_batch = "--"
            if "sj_batch_no" in result:
                sj_batch_no = result["sj_batch_no"]
            else:
                sj_batch_no = "--"
            if "zj_batch_no" in result:
                zj_batch_no = result["zj_batch_no"]
            else:
                zj_batch_no = "--"
        else:
            extract_batch = '--'
            library_batch = '--'
            deal_batch = "--"
            receive_batch = "--"
            sj_batch_no = "--"
            zj_batch_no = "--"
        return extract_batch, library_batch, deal_batch, receive_batch, sj_batch_no, zj_batch_no

    def get_all_batch_sample(self, dad_id, mom_id, son_id):
        sample_list = []
        dad_list = self.find_batch_ids(dad_id)
        mom_list = self.find_batch_ids(mom_id)
        son_list = self.find_batch_ids(son_id)
        dad_list.extend(mom_list)
        dad_list.extend(son_list)
        for m in list(set(dad_list)):
            if "-M" in m or "-F" in m:
                sample_list.append(m)
        return sample_list

    def find_batch_ids(self, sample_id_):
        """
        获取某一个样本的建库，抽提，处理批次的相关所有的样本
        如果某个批次为空，不拿空字符串去检索样本
        :param sample_id_:
        :return:
        """
        sample_id = []
        extract_batch, library_batch, deal_batch, receive_batch, sj_batch_no, zj_batch_no = \
            self.find_father_batch(sample_id_)
        extract_batch = extract_batch if extract_batch else "not used to retrieve samples"
        library_batch = library_batch if library_batch else "not used to retrieve samples"
        deal_batch = deal_batch if deal_batch else "not used to retrieve samples"
        receive_batch = receive_batch if receive_batch else "not used to retrieve samples"
        sj_batch_no = sj_batch_no if sj_batch_no else "not used to retrieve samples"
        zj_batch_no = zj_batch_no if zj_batch_no else "not used to retrieve samples"
        batch_result = self.db['sg_sample_batch'].find({"$or": [{"library_batch": library_batch},
                                                                {"extract_batch": extract_batch},
                                                                {"deal_batch": deal_batch},
                                                                {"receive_batch": receive_batch},
                                                                {"sj_batch_no": sj_batch_no},
                                                                {"zj_batch_no": zj_batch_no}]})
        if batch_result.count() != 0:
            for m in batch_result:
                sample_id.append(m['sample_id'])
        return sample_id

    def update_father_tab_batch(self, father_id, dad_id, types=None):
        """
        用于将家系父本id的父本的建库批次与抽提批次更新到sg_father主表中
        :param father_id:
        :param dad_id:
        :param types:为free的时候是更新自由交互的主表，当为空的时候就直接更新正常的家系分析结果
        :return:
        """
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        extract_batch, library_batch, deal_batch, receive_batch, sj_batch_no, zj_batch_no = \
            self.find_father_batch(dad_id)
        if types and types == 'free':
            collection = "sg_free_combine"
        else:
            collection = "sg_father"
        try:
            self.db[collection].update({"_id": father_id}, {"$set": {"extract_batch": extract_batch,
                                                                     "library_batch": library_batch,
                                                                     "deal_batch": deal_batch,
                                                                     "receive_batch": receive_batch,
                                                                     "sj_batch_no": sj_batch_no,
                                                                     "zj_batch_no": zj_batch_no}}, multi=True,
                                       upsert=True)
        except Exception as e:
            raise Exception("更新父本的建库批次与抽提批次到sg_father中失败{}".format(e))
        else:
            self.bind_object.logger.info("更新父本的建库批次与抽提批次到sg_father中成功！")

    def find_dedup_father_id(self, father_dedup_id):
        """
        用于查找某个father_dedup_id下面的所有的父本 add by hongdong @20171123
        注：设置father_dedup_id为索引
        :param father_dedup_id:
        :return:
        """
        if not isinstance(father_dedup_id, ObjectId):
            if isinstance(father_dedup_id, StringTypes):
                father_dedup_id = ObjectId(father_dedup_id)
            else:
                self.bind_object.logger.info({"father_err_id:{}".format(father_dedup_id)})
                raise Exception("father_dedup_id必须为ObjectId对象或其对应的字符串!")
        dad_list = []
        try:
            results = self.db['sg_father_dedup_detail'].find({"father_dedup_id": father_dedup_id})
            if results.count() != 0:
                for m in results:
                    if str(m['result']) in ['Yes', 'No']:
                        dad_list.append(m['dad_id'])
        except Exception as e:
            self.bind_object.logger.error("查找查重父本的list出错！{}".format(e))
        else:
            self.bind_object.logger.info('查找查重父本的list成功！')
        return dad_list

    def update_to_father_err(self, father_dedup_id, father_err_id):
        """
        用于更新查重的不同err的主表id到sg_father_err表中 add by hongdong @20171123
        :param father_dedup_id:
        :param father_err_id:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        if not isinstance(father_dedup_id, ObjectId):
            if isinstance(father_dedup_id, StringTypes):
                father_dedup_id = ObjectId(father_dedup_id)
            else:
                raise Exception("father_dedup_id必须为ObjectId对象或其对应的字符串!")

        try:
            self.db['sg_father_err'].update({"_id": father_err_id}, {"$set": {"father_dedup_id": father_dedup_id}},
                                            multi=True)
        except Exception as e:
            self.bind_object.logger.error('更新dedup_id到sg_father_err主表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("更新dedup_id到sg_father_err主表成功")

    def add_result_pic(self, file_path, family_name, father_id, father_err_id, postype='hybrid'):
        """
        用于添加不同错配的结果图片 add by hongdong @20171123
        :param file_path:
        :param family_name:
        :param father_id:
        :param father_err_id:
        :param postype:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        snp_pic = file_path + family_name + "_fig1.png"
        err_min_pic_2 = file_path + family_name + "_fig2.png"
        err_min_pic = file_path + family_name + "_family_match.png"
        family_pic = file_path + family_name + "_family.png"
        ss_nongdu = file_path + family_name + "_density_SS.png"
        ms_nongdu = file_path + family_name + "_density_MS.png"
        if postype == 'homohybrid':
            insert_data = {
                'father_err_id': father_err_id,
                'father_id': father_id,
                'snp_pic': file_path + family_name + "_fig1_update.png",
                'err_min_pic': file_path + family_name + "_family_match_update.png"
            }
        else:
            insert_data = {
                'father_err_id': father_err_id,
                'father_id': father_id,
                'snp_pic': snp_pic,
                'err_min_pic': err_min_pic,
                'err_min_pic_2': err_min_pic_2,
                'family_pic': family_pic,
                'ms_nongdu': ms_nongdu,
                'ss_nongdu': ss_nongdu
            }
        try:
            if postype == 'homohybrid':
                self.db['sg_father_err_pic_homohybrid'].insert_one(insert_data)
            else:
                self.db['sg_father_err_pic'].insert_one(insert_data)
        except Exception as e:
            raise Exception('导入sg_father_err_pic出错{}-{}'.format(e, postype))
        else:
            self.bind_object.logger.info('导入sg_father_err_pic成功！-{}'.format(postype))

    def add_test_pos(self, file_path, father_err_id, postype='hybrid'):
        """
        用于导入测试位点信息 add by hongdong @20171123
        :param file_path:
        :param father_err_id:
        :param postype:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "检测位点编号":
                    continue
                if line[5] == 'Mis':
                    mis = '错配'
                else:
                    mis = '-'
                insert_data = {
                    "father_err_id": father_err_id,
                    "test_no": int(line[0]),
                    "chrom": line[1],
                    "dad_geno": line[2],
                    "mom_geno": line[3],
                    "preg_geno": line[4],
                    "is_mis": mis
                }
                sg_pt_family_detail.append(insert_data)
        if len(sg_pt_family_detail) == 0:
            self.bind_object.logger.info("位点信息表格为空")
            pass
        else:
            try:
                if postype == 'homohybrid':
                    self.db['sg_father_err_test_position_homohybrid'].insert_many(sg_pt_family_detail)
                else:
                    self.db['sg_father_err_test_position'].insert_many(sg_pt_family_detail)
            except Exception as e:
                raise Exception('导入位点信息表格出错：{}-{}'.format(e, postype))
            else:
                self.bind_object.logger.info("导入位点信息表格成功!-{}".format(postype))

    def err_postions_list(self, father_err_id, postype='hybrid'):
        """
        用于获取对应的错配位点, err_postions_num是从调试表中找到错配位点的编号，然后更具编号去debug表中去找具体的位点pos
        :param father_err_id:
        :return:
        """
        err_postions_num = []
        err_postions = []
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        if postype == 'homohybrid':
            result = self.db['sg_father_err_test_position_homohybrid'].find({"father_err_id": father_err_id})
            if result.count() != 0:
                for n in result:
                    if n['is_mis'] != '-':
                        err_postions_num.append(str(n['test_no']))
            else:
                raise Exception("{}：在测试位点的表中没有找到相关记录！".format(father_err_id))
            num_pos_map = self.db['sg_father_err_num_pos_map'].find({"father_err_id": father_err_id})
            if num_pos_map.count() != 0:
                for m in num_pos_map:
                    if str(m['4num']) in err_postions_num:
                        err_postions.append(str(m['pos']))
            else:
                raise Exception("{}：在sg_father_err_num_pos_map表中没有找到相关记录！".format(father_err_id))
        else:
            result = self.db['sg_father_err_test_position'].find({"father_err_id": father_err_id})
            if result.count() != 0:
                for n in result:
                    if n['is_mis'] != '-':
                        err_postions_num.append(str(n['test_no']))
            else:
                raise Exception("{}：在测试位点的表中没有找到相关记录！".format(father_err_id))
            result_debug = self.db['sg_father_err_debug'].find({"father_err_id": father_err_id})
            if result_debug.count() != 0:
                for m in result_debug:
                    if str(m['n']) in err_postions_num:
                        err_postions.append(str(m['pos']))
            else:
                raise Exception("{}：在调试表中没有找到相关记录！".format(father_err_id))
        return err_postions

    def is_father(self, father_id, err_max, postype='hybrid'):
        """
        用于检测该家系是否是真父，将err_min的所有的结果都进行检查，如果都为no的时候就return Fales，否则返回True
        :param father_id:
        :param err_max:
        :return:
        """
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        err_analysis_result = []
        if postype == 'homohybrid':
            err_list_valid = [i for i in self.err_list_homohybrid if i <= err_max]
            for n in err_list_valid:
                result = self.db['sg_father_err_homohybrid'].find_one({"father_id": father_id, 'err_min': n,
                                                                       "is_filter": "no"})
                if result:
                    analysis_result = self.db['sg_father_err_result_homohybrid'].find_one(
                        {"father_err_id": result['_id']})
                    if analysis_result:
                        err_analysis_result.append(str(analysis_result['result']))
                    else:
                        raise Exception("father_err_id:{}在father_err_result表格中没有检索到-{}"
                                        .format(result['_id'], postype))
                else:
                    raise Exception("father_id:{}在father_err表格中没有检索到-{}".format(father_id, postype))
        else:
            err_list_valid = [i for i in self.err_list if i <= err_max]
            for n in err_list_valid:
                result = self.db['sg_father_err'].find_one({"father_id": father_id, 'err_min': n, "is_filter": "no"})
                if result:
                    analysis_result = self.db['sg_father_err_result'].find_one({"father_err_id": result['_id']})
                    if analysis_result:
                        err_analysis_result.append(str(analysis_result['result']))
                    else:
                        raise Exception("father_err_id:{}在father_err_result表格中没有检索到".format(result['_id']))
                else:
                    raise Exception("father_id:{}在father_err表格中没有检索到".format(father_id))
        if 'Yes' in err_analysis_result:
            return True
        else:
            return False

    def find_err_max(self, father_id, postype='hybrid'):
        """
        获取错配最大值
        :return:
        """
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        max_err = 5
        if postype == 'homohybrid':
            result = self.db['sg_father_err_homohybrid'].find({"father_id": father_id, "is_filter": "no"})
        else:
            result = self.db['sg_father_err'].find({"father_id": father_id, "is_filter": "no"})
        if result:
            for m in result:
                if int(m['err_min']) > max_err:
                    max_err = int(m['err_min'])
            return max_err
        else:
            raise Exception("father_id:{}在father_err表格中没有检索到".format(father_id))

    def add_analysis_tab(self, file_path, father_err_id, father_id, err_min, postype='hybrid'):
        """
        添加家系分析的匹配结果表 add by hongdong @20171123
        :param file_path:
        :param father_err_id:
        :param father_id:
        :param err_min:
        :param postype:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        sg_pt_family_detail = list()
        match = 'No'
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == "dad.id":
                    continue
                temp_fp = eval(line[4])
                rcp = temp_fp / (temp_fp + 1)
                if rcp > 0.5:
                    rcp_result = "> 99.99%"
                else:
                    rcp_result = "< 0.01%"
                if line[8] == 'Yes':
                    match = line[8]
                insert_data = {
                    "father_err_id": father_err_id,
                    "father_id": father_id,
                    "dad_id": line[0],
                    "test_pos_n": line[1],
                    "err_pos_n": line[2],
                    "err_rate": line[3],
                    "fq": line[4],
                    "dp": line[5],
                    "eff_rate": line[6],
                    "ineff_rate": line[7],
                    "result": line[8],
                    "rcp": rcp_result,
                    "err_min": err_min,
                    "eff_ineff_rate": str(float(line[6].split(":")[0].strip()[:-1]) - float(line[7].split(":")
                                                                                            [0].strip()[:-1])) + "%"
                }
                sg_pt_family_detail.append(insert_data)
            try:
                if postype == 'homohybrid':
                    collection = self.db['sg_father_err_result_homohybrid']
                else:
                    collection = self.db['sg_father_err_result']
                collection.insert_many(sg_pt_family_detail)
            except Exception as e:
                raise Exception('导入sg_father_err_result表格出错：{}-{}'.format(e, postype))
            else:
                self.bind_object.logger.info("导入sg_father_err_result表格成功-{}".format(postype))
            if err_min == 5 and postype != 'homohybrid':
                try:
                    self.db['sg_father'].update({"_id": father_id}, {"$set": {"result": match}}, multi=True)
                except Exception as e:
                    raise Exception('错配为5-更新sg_father中父本与胎儿匹配失败：{}'.format(e))

    def export_tab_file(self, sample, out_dir, new_sample=None, is_ms=False):
        """
        导出tab file,生成sample_id命名的tab文件 add by hongdong @ 20171122
        :param sample:
        :param out_dir:
        :param new_sample:
        :param is_ms:
        :return:
        """
        if not is_ms:
            collection = self.db['sg_sample_tab']
        else:
            collection = self.db['sg_sample_tab_ms']
        if new_sample:
            sample_tab = new_sample + '.tab'
            sample_name = new_sample
        else:
            sample_tab = sample + '.tab'
            sample_name = sample
        file_path = os.path.join(out_dir, sample_tab)

        if os.path.exists(file_path):
            pass
        else:
            search_result = collection.find({"sample_id": sample})  # 读出来是个地址
            if search_result.count() != 0:  # 判断是否找到了相应的结果
                with open(file_path, 'w+') as f:
                    for i in search_result:
                        f.write(sample_name + '\t' + i['chrom'] + '\t' + i['pos'] + '\t' + i['ref'] + '\t' +
                                i['alt'] + '\t' + i['dp'] + '\t' + i['ref_dp'] + '\t' + i['alt_dp'] + '\n')
                self.bind_object.logger.info("导出样本{}的tab文件成功！".format(sample))
            else:
                raise Exception('没有在数据库中搜到样本{}'.format(sample))
            if os.path.getsize(file_path):
                return file_path
            else:
                raise Exception('样本数据{}的tab文件为空，可能还未下机'.format(sample))

    def attempt_to_export_tab_file(self, sample, out_dir, new_sample=None, err_postions=None):
        """
        按照sg_sample_tab、sg_sample_tab_ms的优先级尝试导出tab文件，一旦成功不再尝试，全部无法导出则报错
        :param sample:
        :param out_dir:
        :param new_sample:
        :param err_postions:
        :return:
        """
        collection = self.db['sg_sample_tab']
        collection_ms = self.db['sg_sample_tab_ms']
        if new_sample:
            sample_tab = new_sample + '.tab'
            sample_name = new_sample
        else:
            sample_tab = sample + '.tab'
            sample_name = sample
        file_path = os.path.join(out_dir, sample_tab)

        if self.tab_exist(sample):
            pass
        else:
            raise Exception("导出tab文件前样品{}在sg_sample_name表中验证失败".format(sample))

        if os.path.exists(file_path):
            pass
        else:
            dedup_output = set()
            search_result = collection.find({"sample_id": sample})  # 读出来是个地址
            search_result_ms = collection_ms.find({"sample_id": sample})
            if search_result.count() != 0:  # 判断是否找到了相应的结果
                for i in search_result:
                    if err_postions:
                        if str(i['pos']) not in err_postions:
                            dedup_output.add(sample_name + '\t' + i['chrom'] + '\t' + i['pos'] + '\t' + i['ref']
                                             + '\t' + i['alt'] + '\t' + i['dp'] + '\t' + i['ref_dp'] + '\t'
                                             + i['alt_dp'] + '\n')
                    else:
                        dedup_output.add(sample_name + '\t' + i['chrom'] + '\t' + i['pos'] + '\t' + i['ref'] + '\t'
                                         + i['alt'] + '\t' + i['dp'] + '\t' + i['ref_dp'] + '\t' + i['alt_dp'] + '\n')
                with open(file_path, 'w+') as f:
                    for i in sorted(dedup_output):
                        f.write(i)
                    self.bind_object.logger.info("导出样本{}的tab文件成功！".format(sample))
            elif search_result_ms.count() != 0:
                for i in search_result_ms:
                    if err_postions:
                        if str(i['pos']) not in err_postions:
                            dedup_output.add(sample_name + '\t' + i['chrom'] + '\t' + i['pos'] + '\t' + i['ref'] +
                                             '\t' + i['alt'] + '\t' + i['dp'] + '\t' + i['ref_dp'] + '\t'
                                             + i['alt_dp'] + '\n')
                    else:
                        dedup_output.add(sample_name + '\t' + i['chrom'] + '\t' + i['pos'] + '\t' + i['ref'] + '\t'
                                         + i['alt'] + '\t' + i['dp'] + '\t' + i['ref_dp'] + '\t' + i['alt_dp'] + '\n')
                with open(file_path, 'w+') as f:
                    for i in sorted(dedup_output):
                        f.write(i)
                    self.bind_object.logger.info("导出样本{}的tab文件成功！".format(sample))
            else:
                raise Exception('没有在数据库中搜到样本{}'.format(sample))
            if os.path.getsize(file_path):
                return file_path
            else:
                raise Exception('样本数据{}的tab文件为空，可能还未下机'.format(sample))

    def get_board_name(self, datasplit_id):
        """
        获取拆分的板子名字
        add by hongdong @ 20171126
        :param datasplit_id:
        :return:
        """
        if not isinstance(datasplit_id, ObjectId):
            if isinstance(datasplit_id, StringTypes):
                datasplit_id = ObjectId(datasplit_id)
            else:
                raise Exception("datasplit_id必须为ObjectId对象或其对应的字符串!")
        result = self.db['sg_datasplit'].find_one({"_id": datasplit_id})
        if result:
            if result['split_type'] not in ['PE', 'SE']:
                raise Exception("split_type{},必须为PE 或者 SE".format(result['split_type']))
            return result['board_batch'], result['split_type']
        else:
            raise Exception("{}在sg_datasplit中没有找到board_batch or split_type字段！".format(datasplit_id))

    def get_sample_list(self, board_batch):
        """
        用于根据板号去找到这个板子中对应的样本id, 查找逻辑，首先查找到板子中所有的样本名，然后检查这个样本名是不是亲子的，
        是的话再检查是不是tab文件已经有了，没有的话，才会添加到运行队列中
        add by hongdong @ 20171126
        :param board_batch:
        :return:
        """
        sample_list = []
        result = self.db['sg_sample_info'].find({"board_batch": board_batch})
        if result:
            for m in result:
                name = "-".join([m['case_name'], m['sample_id']])
                if re.match('WQ([0-9]*)-.*', name):
                    if self.db['sg_sample_tab'].find_one({"sample_id": name}):
                        self.bind_object.logger.info("样本{}已经存在于库中，不添加到运行队列中")
                    else:
                        self.bind_object.logger.info("样本{}不存在，将添加到运行队列中")
                        sample_list.append(name)
        else:
            raise Exception("sg_sample_info表中不存在{}的信息！".format(board_batch))

    def get_family_list(self, samples):
        """
        根据样本列表获取到家系组合后的家系列表[['WQ123', 'WQ1234', 'WQ234'], ['WQ111', 'WQ222', 'WQ333']]
        add by hongdong @ 20171126
        :param samples:
        :return:
        """
        family_id_list = []
        family_list = []
        main_collection = self.db['sg_sample_name']
        for i in samples:
            m = re.match('WQ([0-9]{2,})-(M|F|S)(.*)', i)
            family_id = 'WQ' + m.group(1)
            if family_id in family_id_list:
                continue
            else:
                family_id_list.append(family_id)
            sample_dad = main_collection.find({"sample_id": {"$regex": family_id + "-F.*"}})
            sample_mom = main_collection.find({"sample_id": {"$regex": family_id + "-M.*"}})
            sample_son = main_collection.find({"sample_id": {"$regex": family_id + "-S.*"}})
            if sample_dad and sample_mom and sample_son:
                dad_id = []
                mom_id = []
                son_id = []
                for d in sample_dad:
                    dad_id.append(d['sample_id'])
                    dad_id = list(set(dad_id))
                for m in sample_mom:
                    mom_id.append(m['sample_id'])
                    mom_id = list(set(mom_id))
                for s in sample_son:
                    son_id.append(s['sample_id'])
                    son_id = list(set(son_id))
                for dad in dad_id:
                    for mom in mom_id:
                        for son in son_id:
                            family_member = list()
                            family_member.append(dad)
                            family_member.append(mom)
                            family_member.append(son)
                            family_list.append(family_member)
            else:
                self.bind_object.logger.info("样本 {} 没有组建成家系".format(i))
                continue
        return family_list

    def set_new_sample_to_empty(self, samples):
        collection = self.db["sg_family"]
        casenames = list()
        for sample in samples:
            # m = re.match('WQ(.*)-(M|F|S)(.*)', sample)
            # casename = 'WQ' + m.group(1)
            casenames.append(sample.split("-")[0])
        for case_name in list(set(casenames)):
            result = collection.find_one({"case_name": case_name})
            if result and "new_sample" in result.keys():
                try:
                    collection.update({"case_name": case_name}, {"$set": {"new_sample": ""}}, multi=True)
                except Exception as e:
                    raise Exception("清除{}家系的new_sample出错{}".format(case_name, e))
                else:
                    self.bind_object.logger.info("清除{}家系的new_sample成功".format(case_name))

    def add_sg_family(self, samples):
        """
        该函数用于添加sg_family表格,或者更新样本到sg_family, 注后面的时间统计还要重新写
        特加急：\u7279\u52a0\u6025， 加急：\u52a0\u6025， 普通：\u666e\u901a
        add by hongdong @ 20171126
        :param samples:
        :return:
        """
        self.set_new_sample_to_empty(samples)
        collection = self.db["sg_family"]
        case_name_valid = []                                           # 需要更新受理日期、报告截止日期的case_name
        for sample in samples:
            # sample : WQ1234-2-F1, WQ1234-F1
            if sample.startswith("TEST"):                              # 跳过测试样品
                continue
            n = re.match('WQ(.*)-(M|F|S)(.*)', sample)
            if n:
                id_ = n.group(2) + n.group(3)
            else:
                raise Exception("亲子样本名{}命令不规范--流程无法继续！".format(sample))
            # case_name = 'WQ' + m.group(1)
            case_name = sample.split("-")[0]
            if "-F" in sample or "-M" in sample or "-SC" in sample:     # 只下机了S样品时，不更新受理日期、报告截止日期
                case_name_valid.append(case_name)
            batch_info = self.db['sg_file_info'].find_one({"sample_id": sample})
            report_status = '1'
            if batch_info and "urgence" in batch_info.keys():
                if batch_info['urgence'] == u'\u52a0\u6025\u0031\u5929':  # 加急1天
                    report_status = '4'
                elif batch_info['urgence'] == u'\u52a0\u6025\u0032\u5929':  # 加急2天
                    report_status = '3'
                elif batch_info['urgence'] == u'\u52a0\u6025\u0033\u5929':  # 加急3天
                    report_status = '2'
                self.bind_object.logger.info(batch_info['urgence'])
            self.bind_object.logger.info(report_status)
            # accept_time, report_time = self.get_report_time(sample, report_status)  # 获取受理日期与报告截止日期
            result = collection.find_one({"case_name": case_name})
            if result:
                self.bind_object.logger.info("case_name {} 已经在sg_family存在，后面将更"
                                             "新{}相关信息！".format(case_name, sample))
                if str(result['new_sample']) != '':
                    if id_ not in str(result['new_sample']):
                        new_sample = str(result['new_sample']) + ',' + id_
                    else:
                        new_sample = str(result['new_sample'])
                else:
                    new_sample = n.group(2) + n.group(3)
                if '-S' in sample:
                    id_s = id_
                    if str(result['preg_ids']) != '':
                        if id_s not in str(result['preg_ids']):
                            preg_ids = str(result['preg_ids']) + ',' + id_s
                        else:
                            preg_ids = str(result['preg_ids'])
                    else:
                        preg_ids = id_s
                    update_data = {
                        "preg_ids": preg_ids,
                        "new_sample": new_sample,
                        "update_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "report_status": report_status,
                        # "accept_time": accept_time,
                        # "report_time": str(report_time).split(" ")[0]
                    }
                elif '-M' in sample:
                    id_m = id_
                    if str(result['mom_ids']) != '':
                        if id_m not in str(result['mom_ids']):
                            mom_ids = str(result['mom_ids']) + ',' + id_m
                        else:
                            mom_ids = str(result['mom_ids'])
                    else:
                        mom_ids = id_m
                    update_data = {
                        "mom_ids": mom_ids,
                        "new_sample": new_sample,
                        "update_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "report_status": report_status,
                        # "accept_time": accept_time,
                        # "report_time": str(report_time).split(" ")[0]
                    }
                elif '-F' in sample:
                    id_f = id_
                    if str(result['dad_ids']) != '':
                        if id_f not in str(result['dad_ids']):
                            dad_ids = str(result['dad_ids']) + ',' + id_f
                        else:
                            dad_ids = str(result['dad_ids'])
                    else:
                        dad_ids = id_f
                    report_undone = result['report_undone'].strip().split(",")
                    if "" in report_undone:
                        report_undone.remove("")
                    report_done = result['report_done'].strip().split(",")
                    if id_f not in report_undone and id_f not in report_done:
                        report_undone.append(id_f)
                    update_data = {
                        "dad_ids": dad_ids,
                        "new_sample": new_sample,
                        "update_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "report_undone": ",".join(report_undone),
                        "report_status": report_status,
                        # "accept_time": accept_time,
                        # "report_time": str(report_time).split(" ")[0]
                        # "is_report": '2'   # 有新父本更新进去的时候，出报告的状态改为2，没有出报告
                    }
                else:
                    raise Exception("样本{}命名不规范！".format(sample))
                # if "preg_ids" in result.keys() and "dad_ids" in result.keys() and "mom_ids" in result.keys():
                #     if result["preg_ids"] and result["dad_ids"] and result["mom_ids"]:
                #         update_data.update({"analysis_status": "2"})
                # else:
                #     raise Exception("preg_ids or dad_ids or mom_ids 不在sg_family中！")
                try:
                    collection.update({"case_name": case_name}, {"$set": update_data}, multi=True, upsert=True)
                except Exception as e:
                    raise Exception("更新样本{}到sg_family表格出错{}".format(sample, e))
                else:
                    self.bind_object.logger.info("更新样本{}到sg_family表格成功".format(sample))
            else:
                self.bind_object.logger.info("case_name {} 在sg_family不存在，后面将插入{}相关信"
                                             "息！".format(case_name, sample))
                insert_data = {
                    "case_name": case_name,
                    "new_sample": id_,
                    "update_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    # "accept_time": accept_time,
                    # "report_time": str(report_time).split(" ")[0],
                    "report_done": "",
                    "report_undone": id_ if '-F' in sample else "",
                    "remarks": "",
                    "report_status": report_status,
                    "analysis_status": "1",  # 用于区分该家系中是否有报告进行分析，3 表示有家系并且能够进行出报告，2 有家系但是报告还在分析中，1 缺样本
                    "is_report": "2"  # 用于区分family中父本是否都出报告了，2默认为没有出报告
                }
                if "-S" in sample:
                    insert_data.update({"preg_ids": id_, "dad_ids": "", "mom_ids": ""})
                elif '-M' in sample:
                    insert_data.update({"preg_ids": "", "dad_ids": "", "mom_ids": id_})
                elif '-F' in sample:
                    insert_data.update({"preg_ids": "", "dad_ids": id_, "mom_ids": ""})
                else:
                    raise Exception("样本{}命名不规范！".format(sample))
                try:
                    collection.insert_one(insert_data)
                except Exception as e:
                    raise Exception("导入样本{}到sg_family表格出错{}".format(sample, e))
                else:
                    self.bind_object.logger.info("导入样本{}到sg_family表格成功".format(sample))
        for case in list(set(case_name_valid)):
            self.update_family_accept_time_v3(case)

    def get_report_time(self, sample, report_status):
        """
        用于获取报告的受理日期与计算出报告的截止日期
        :param sample:  样本名称
        :param report_status:  报告的加急状态， 1,2,3
        :return:
        """
        if "-T" in sample:  # add by hongdong 20180102 重上机的样本客户信息获取
            sample_ = re.match('(WQ.*)-(T.*)', sample).group(1)
            guest_info = self.db['sg_guest_info'].find_one({"name": sample_})
        elif "-L" in sample:
            sample_ = re.match('(WQ.*)-(L.*)', sample).group(1)
            guest_info = self.db['sg_guest_info'].find_one({"name": sample_})
        else:
            guest_info = self.db['sg_guest_info'].find_one({"name": sample})
        accept_time = '--'
        report_time = '--'
        if guest_info and "accept_time" in guest_info.keys():
            accept_time = guest_info['accept_time']
            time_ = accept_time.strip().split("-")
            d1 = datetime.datetime(int(time_[0]), int(time_[1]), int(time_[2]))
            if report_status == '1':
                report_time = d1 + datetime.timedelta(days=7)
            elif report_status == '2':
                report_time = d1 + datetime.timedelta(days=3)
            elif report_status == '3':
                report_time = d1 + datetime.timedelta(days=2)
        return accept_time, report_time

    def set_T_L(self, m):
        if "-T" in m:
            _id = re.match('(WQ.*)-(T.*)', m).group(1)
        elif "-L" in m:
            _id = re.match('(WQ.*)-(L.*)', m).group(1)
        else:
            _id = m
        return _id

    def get_case_all_sample(self, case_name):
        """
        根据case号 获取所有的样本
        :return:
        """
        samples = []
        f_id = []
        m_id = []
        s_id = []
        result = self.get_case_sampls(case_name)
        if result:
            samples = [temp['sample_id'] for temp in result]
        for m in list(set(samples)):
            if "-F" in m:
                f_id.append(self.set_T_L(m))
            elif "-S" in m:
                s_id.append(self.set_T_L(m))
            elif "-M" in m:
                m_id.append(self.set_T_L(m))
        return f_id, m_id, s_id,

    def update_family_accept_time_v3(self, case_name):
        """
        用于更新family表中的受理日期，family表中的受理日期是能够组建成家系中的样本最晚受理日期，sg_father中的受理日期是三个样本
        的最晚的受理日期,要请求意见系统的一个接口，组建一个父本，母本，请求接口后会返回一个这个case号中的所有的受理日期
        :param case_name:
        :return:
        """
        report_time = ''
        res_fa = self.db['sg_family'].find_one({"case_name": case_name})
        acc_time = []
        f_ids, m_ids, s_ids = self.get_case_all_sample(case_name)
        if len(f_ids) == 0 or len(m_ids) == 0:
            self.bind_object.logger.info("缺少样本，流程正常结束！")
            return
        for f in f_ids:
            for m in m_ids:
                res = self.gest_info_post(f, m)
                if res and "audit_date" in res.keys() and res['audit_date']:
                    acc_time.append(res['audit_date'])
                else:
                    self.bind_object.logger.info("{}_{}没有找到客户信息表，--正常退出！".format(f, m))
        lasted_acc_time = self.get_last_acc_time(acc_time)
        if res_fa and lasted_acc_time:
            report_status = res_fa['report_status']
            lasted_acc_time__ = datetime.datetime(int(lasted_acc_time[0]), int(lasted_acc_time[1]),
                                                  int(lasted_acc_time[2]))
            if report_status == '4':
                report_time = lasted_acc_time__ + datetime.timedelta(days=1)
            elif report_status == '3':
                report_time = lasted_acc_time__ + datetime.timedelta(days=2)
            elif report_status == '2':
                report_time = lasted_acc_time__ + datetime.timedelta(days=3)
            elif report_status == '1':
                report_time = lasted_acc_time__ + datetime.timedelta(days=7)
            else:
                self.bind_object.logger.info("样本的加急信息不正确--警告！")
        else:
            self.bind_object.logger.info("整理出报告时间失败--程序正常退出！")
            return
        if report_time:
            try:
                self.db['sg_family'].update({"case_name": case_name}, {"$set": {"accept_time": str(lasted_acc_time__),
                                                                                "report_time": str(report_time),
                                                                                }},
                                            multi=True, upsert=True)
            except Exception as e:
                raise Exception("更新准确的样本受理时间与出报告时间失败{}".format(e))
            else:
                self.bind_object.logger.info("更新准确的样本受理时间与出报告时间成功！")

    def get_last_acc_time(self, acc_time):
        """
        输入一组时间列表， 获取到最新的一个时间
        :return:
        """
        lasted_acc_time = ''
        if acc_time:
            lasted_acc_time = acc_time[0].strip().split("-")
            lasted_acc_time_ = datetime.datetime(int(lasted_acc_time[0]), int(lasted_acc_time[1]),
                                                 int(lasted_acc_time[2]))
            for m in acc_time[1:]:
                time_ = m.strip().split("-")
                time1 = datetime.datetime(int(time_[0]), int(time_[1]), int(time_[2]))
                if (time1 - lasted_acc_time_).days > 0:
                    lasted_acc_time = time_
        # else:
        #     self.bind_object.logger.info('受理日期列表为空，不进行添加')
        return lasted_acc_time

    def update_family_accept_time(self, case_name):

        """
        用于更新family表中的受理日期，family表中的受理日期是能够组建成家系中的样本最晚受理日期，sg_father中的受理日期是三个样本
        的最晚的受理日期,逻辑有点复杂,最后如果运行了该家系的话，那么就要将sg_family表中的is_report状态改为2，未出报告
        :param case_name:
        :return:
        """
        report_time = ''
        res_fa = self.db['sg_family'].find_one({"case_name": case_name})
        sample_list = []
        fm_list = []
        s_list = []
        acc_time = []
        if res_fa and "dad_ids" in res_fa and "mom_ids" in res_fa and "preg_ids" in res_fa and "new_sample" in res_fa:
            sample_list.extend(res_fa["dad_ids"].strip().split(","))
            sample_list.extend(res_fa["mom_ids"].strip().split(","))
            sample_list.extend(res_fa["preg_ids"].strip().split(","))
            sample_list.extend(res_fa["new_sample"].strip().split(","))
        for sample in list(set(sample_list)):
            if "F" in sample or "M" in sample:
                fm_list.append(case_name + "-" + sample)
            if "S" in sample:
                s_list.append(case_name + "-" + sample)

        m = re.compile(r'(WQ\d+-SC\d*)$')
        n = re.compile(r'(WQ\d+-SC\d*)-')
        for son_id in s_list:
            if m.match(son_id) or n.match(son_id):     # 如果胎儿是SC样品，将对应MC样品的受理时间加入列表中进行比较
                son_id_ = m.match(son_id).group(1) if m.match(son_id) else n.match(son_id).group(1)
                mc_id = son_id_.replace("SC", "MC")
                mc = self.db['sg_guest_info'].find_one({"name": {"$regex": "{}$|{}-".format(mc_id, mc_id)}})
                if mc and "accept_time" in mc.keys():
                    acc_time.append(mc['accept_time'])
                else:
                    self.bind_object.logger.info("没有在sg_guest_info中查到{}的接受时间的信息".format(mc_id))

        for fm_id in fm_list:
            if "-T" in fm_id:                       # 重新上机样本，报告组反馈
                fm_id = re.match('(WQ.*)-(T.*)', fm_id).group(1)
            if "-L" in fm_id:                       # 实验端自己反馈
                fm_id = re.match('(WQ.*)-(L.*)', fm_id).group(1)
            sample_gu = self.db['sg_guest_info'].find_one({"name": fm_id})
            if sample_gu and "accept_time" in sample_gu.keys():
                acc_time.append(sample_gu["accept_time"])
            else:
                self.bind_object.logger.info("没有在sg_guest_info中查到{}的接受时间的信息".format(fm_id))
        if acc_time:
            lasted_acc_time = acc_time[0].strip().split("-")
            lasted_acc_time_ = datetime.datetime(int(lasted_acc_time[0]), int(lasted_acc_time[1]),
                                                 int(lasted_acc_time[2]))
            for m in acc_time[1:]:
                time_ = m.strip().split("-")
                time = datetime.datetime(int(time_[0]), int(time_[1]), int(time_[2]))
                if (time - lasted_acc_time_).days > 0:
                    lasted_acc_time = time_
        else:
            self.bind_object.logger.info("最晚受理日期查找失败--程序正常退出！")
            return
        if res_fa:
            report_status = res_fa['report_status']
            lasted_acc_time__ = datetime.datetime(int(lasted_acc_time[0]), int(lasted_acc_time[1]),
                                                  int(lasted_acc_time[2]))
            if report_status == '4':
                report_time = lasted_acc_time__ + datetime.timedelta(days=1)
            elif report_status == '3':
                report_time = lasted_acc_time__ + datetime.timedelta(days=2)
            elif report_status == '2':
                report_time = lasted_acc_time__ + datetime.timedelta(days=3)
            elif report_status == '1':
                report_time = lasted_acc_time__ + datetime.timedelta(days=7)
            else:
                self.bind_object.logger.info("样本的加急信息不正确--警告！")
        else:
            self.bind_object.logger.info("整理出报告时间失败--程序正常退出！")
            return
        if report_time:
            try:
                self.db['sg_family'].update({"case_name": case_name}, {"$set": {"accept_time": str(lasted_acc_time__),
                                                                                "report_time": str(report_time),
                                                                                }},
                                            multi=True, upsert=True)
            except Exception as e:
                raise Exception("更新准确的样本受理时间与出报告时间失败{}".format(e))
            else:
                self.bind_object.logger.info("更新准确的样本受理时间与出报告时间成功！")

    def update_report_status(self, dad_id):
        """
        如果运行了该家系的话，那么就要将sg_family表中的is_report状态改为2，未出报告
        :param dad_id:
        :return:
        """
        # m = re.match('WQ([0-9]{2,})-(M|F|S)(.*)', dad_id)
        m = re.match('WQ([0-9]{2,})-(.*)', dad_id)
        case_name = 'WQ' + m.group(1)
        try:
            self.db['sg_family'].update({"case_name": case_name}, {"$set": {'is_report': "2"}}, multi=True, upsert=True)
        except Exception as e:
            raise Exception("更新报告状态失败：{}".format(e))
        else:
            self.bind_object.logger.info("更新报告状态成功！")

    def add_sg_family_info(self, father_id, dad_id, mom_id, is_free=None):
        """
        出报告的时候需要的所有的家系信息，以及后面报告组需要的查看出了哪些报告以及截至时间等--这个不用了
        :param father_id:
        :param dad_id:
        :param mom_id:
        :param is_free:
        :return:
        """
        ask_time_list = []
        accept_time_list = []
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        if self.db['sg_family_info'].find_one({"name": "_".join([dad_id, mom_id])}):
            self.db['sg_family_info'].remove({"name": "_".join([dad_id, mom_id])})
        son = self.db['sg_father'].find_one({"_id": father_id})
        if son and "preg_id" in son.keys():
            preg_id = son["preg_id"]
            m = re.compile(r'(WQ\d+-SC\d*)$')
            n = re.compile(r'(WQ\d+-SC\d*)-')
            if m.match(preg_id) or n.match(preg_id):
                preg_id_ = m.match(preg_id).group(1) if m.match(preg_id) else n.match(preg_id).group(1)
                mc_id = preg_id_.replace("SC", "MC")
                mc = self.db['sg_guest_info'].find_one({"name": {"$regex": "{}$|{}-".format(mc_id, mc_id)}})
                if mc and "ask_time" in mc.keys() and "accept_time" in mc.keys():
                    ask_time_list.append(mc['ask_time'])
                    accept_time_list.append(mc['accept_time'])
        if '-T' in dad_id:  # 重新上机样本，报告组反馈
            dad_id_ = re.match('(WQ.*)-(T.*)', dad_id).group(1)
            dad = self.db['sg_guest_info'].find_one({"name": dad_id_})
        elif "-L" in dad_id:   # 实验端自己反馈
            dad_id_ = re.match('(WQ.*)-(L.*)', dad_id).group(1)
            dad = self.db['sg_guest_info'].find_one({"name": dad_id_})
        else:
            dad = self.db['sg_guest_info'].find_one({"name": dad_id})
        if '-T' in mom_id:
            mom_id_ = re.match('(WQ.*)-(T.*)', mom_id).group(1)
            mom = self.db['sg_guest_info'].find_one({"name": mom_id_})
        elif "-L" in mom_id:
            mom_id_ = re.match('(WQ.*)-(L.*)', mom_id).group(1)
            mom = self.db['sg_guest_info'].find_one({"name": mom_id_})
        else:
            mom = self.db['sg_guest_info'].find_one({"name": mom_id})
        if dad and "ask_time" in dad.keys() and "accept_time" in dad.keys():
            ask_time_list.append(dad['ask_time'])
            accept_time_list.append(dad['accept_time'])
        if mom and "ask_time" in mom.keys() and "accept_time" in mom.keys():
            ask_time_list.append(mom['ask_time'])
            accept_time_list.append(mom['accept_time'])
        ask_time = "--"
        if ask_time_list:
            ask_time = ask_time_list[0]
            for i in range(1, len(ask_time_list)):
                time_1 = ask_time.strip().split("-")
                time1 = datetime.datetime(int(time_1[0]), int(time_1[1]), int(time_1[2]))
                time_2 = ask_time_list[i].strip().split("-")
                time2 = datetime.datetime(int(time_2[0]), int(time_2[1]), int(time_2[2]))
                if (time2 - time1).days < 0:
                    ask_time = time2.strftime("%Y-%m-%d")
        accept_time = '--'
        if accept_time_list:
            accept_time = accept_time_list[0]
            for j in range(1, len(accept_time_list)):
                time_1 = accept_time.strip().split("-")
                time1 = datetime.datetime(int(time_1[0]), int(time_1[1]), int(time_1[2]))
                time_2 = accept_time_list[j].strip().split("-")
                time2 = datetime.datetime(int(time_2[0]), int(time_2[1]), int(time_2[2]))
                if (time2 - time1).days > 0:
                    accept_time = time2.strftime("%Y-%m-%d")
        insert_data = {
            "case_id": re.match('(WQ.*)-F(.*)', dad_id).group(1),
            "name": "_".join([dad_id, mom_id]),
            "ask_person": dad['ask_person'] if dad and "ask_person" in dad.keys() else "--",
            "mom_name": mom['sample_name'] if mom and "sample_name" in mom.keys() else "--",
            "mom_type": mom['sample_type'] if mom and "sample_type" in mom.keys() else "--",
            "mom_number": mom['sample_number'] if mom and "sample_number" in mom.keys() else "--",
            "mom_id": re.match('(.*)-(M.*)', mom_id).group(2),
            "dad_name": dad['sample_name'] if dad and "sample_name" in dad.keys() else "--",
            "dad_type": dad['sample_type'] if dad and "sample_type" in dad.keys() else "--",
            "dad_number": dad['sample_number'] if dad and "sample_number" in dad.keys() else "--",
            "dad_id": re.match('(.*)-(F.*)', dad_id).group(2),
            "ask_time": ask_time,
            "accept_time": accept_time,
            "F_accept_time": dad['accept_time'] if dad and 'accept_time' in dad.keys() else "--",
            "F_ask_time": dad['ask_time'] if dad and 'ask_time' in dad.keys() else "--",
            "F_result_time": dad['accept_time'] if dad and 'accept_time' in dad.keys() else "--",
            "M_accept_time": mom['accept_time'] if mom and 'accept_time' in mom.keys() else "--",
            "M_ask_time": mom['ask_time'] if mom and 'ask_time' in mom.keys() else "--",
            "M_result_time": mom['accept_time'] if mom and 'accept_time' in mom.keys() else "--",
            "report_time": "",
            "gestation_week": mom['gestation_week'] if mom and 'gestation_week' in mom.keys() else "--",
            "company_name": mom['company_name'] if mom and 'company_name' in mom.keys() else "--",
            "contacts_people": mom['contacts_people'] if mom and 'contacts_people' in mom.keys() else "--",
            "father_id": father_id,
            "insert_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        if is_free and is_free == 'free':
            insert_data.update({"is_free": 'yes'})   # 报告汇总的表格中添加is_free的标记用于区分自由交互的报告
        try:
            self.db['sg_family_info'].insert_one(insert_data)
        except Exception as e:
            raise Exception("sg_family_info导入失败{}".format(e))
        else:
            self.bind_object.logger.info("sg_family_info导入成功!")

    def find_params_result(self, params):
        """
        用于查找家系分析的主表参数，返回值为true or false,找到相关参数的时候返回true 否则返回false
        添加查询条件"is_show": "1"， 当为1的时候肯定是家系分析结束的，当为2的时候，存在是因为缺少样本，流程终止，
        但是运行状态还是为end
        :param params:
        :return:
        """
        reslut = self.db['sg_father'].find({"params": params, 'status': {'$in': ['end']}, "is_show": "1"})
        for m in reslut:
            self.bind_object.logger.info("result[params]:{}".format(m['params']))
        if reslut.count() == 0:
            self.bind_object.logger.info("{}参数比对没有找到相关结果, 将继续后面的家系分析！".format(params))
            return False
        else:
            return True

    def tab_exist(self, sample_id):
        """
        用于判断该样本的tab文件是否存在
        add by hongdong
        :param sample_id:
        :return:
        """
        result = self.db['sg_sample_name'].find_one({"sample_id": sample_id, "is_tab": "yes"})
        if result:
            return True
        else:
            return False

    def check_end(self, id_, collection):
        """
        用于更新程序运行的总进度 add by hongdong @ 20171122
        :param id_:
        :param collection:
        :return:
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        result = self.db[collection].find_one({"_id": id_})
        if result:
            return result['status']
        else:
            raise Exception("查找{}的状态status失败！".format(collection))

    def check_sample_in_board(self, batch_id, sample_name):
        """
        用于检查样本是不是在这个板子中 add by hongdong @ 20171122
        :param batch_id:
        :param sample_name:
        :return:
        """
        if not isinstance(batch_id, ObjectId):
            if isinstance(batch_id, StringTypes):
                batch_id = ObjectId(batch_id)
            else:
                raise Exception("batch_id必须为ObjectId对象或其对应的字符串!")
        try:
            result = self.db['sg_datasplit'].find_one({"_id": batch_id})
            if result:
                if sample_name in result['wq_sample']:
                    self.bind_object.logger.info("样本{}在当前板子中，即将进行call snp！".format(sample_name))
                    return True
            else:
                self.bind_object.logger.info("sg_datasplit表格中不存在batch_id为{}的记录，流程终止！".format(batch_id))
                return False
        except Exception as e:
            raise Exception("检查样本{}是否在{}板子中出错{}".format(sample_name, batch_id, e))

    def export_tab_file_select(self, sample, out_dir, err_postions=None, new_sample=None, tab_name="sg_sample_tab"):
        """
        按照位点进行过滤，导出tab file,生成sample_id命名的tab文件，add by hongdong @ 20171122
        :param sample:
        :param out_dir:
        :param err_postions:
        :param new_sample:
        :param tab_name:
        :return:
        """
        collection = self.db[tab_name]
        if new_sample:
            sample_tab = new_sample + '.tab'
            sample_name = new_sample
        else:
            sample_tab = sample + '.tab'
            sample_name = sample
        file_path = os.path.join(out_dir, sample_tab)

        if os.path.exists(file_path):
            pass
        else:
            search_result = collection.find({"sample_id": sample})  # 读出来是个地址
            if search_result.count() != 0:  # 判断是否找到了相应的结果
                with open(file_path, 'w+') as f:
                    for i in search_result:
                        if err_postions:
                            if str(i['pos']) not in err_postions:
                                f.write(sample_name + '\t' + i['chrom'] + '\t' + i['pos'] + '\t' + i['ref'] + '\t' +
                                        i['alt'] + '\t' + i['dp'] + '\t' + i['ref_dp'] + '\t' + i['alt_dp'] + '\n')
                        else:
                            f.write(sample_name + '\t' + i['chrom'] + '\t' + i['pos'] + '\t' + i['ref'] + '\t' +
                                    i['alt'] + '\t' + i['dp'] + '\t' + i['ref_dp'] + '\t' + i['alt_dp'] + '\n')
                self.bind_object.logger.info("导出样本{}的tab文件成功！".format(sample))
            else:
                raise Exception('没有在数据库中搜到样本{}'.format(sample))
            if os.path.getsize(file_path):
                return file_path
            else:
                raise Exception('样本数据{}的tab文件为空，可能还未下机'.format(sample))

    def check_report_result(self, name):
        """
        用于检查该家系是否已经存入报告了,检测到存入报告了，如果该次传进来的father_err_id是之前已经传入的，那么就是一模一样的数据，
        然后就更新该组数据的is_report为no
        add by hongdong @ 20171122
        :param name:
        :return:
        """
        # if not isinstance(father_err_id, ObjectId):
        #     if isinstance(father_err_id, StringTypes):
        #         father_err_id = ObjectId(father_err_id)
        #     else:
        #         raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        try:
            # rr = self.db['sg_father_summary'].find_one({"father_err_id": father_err_id})
            # if rr:
            #     return True
            result = self.db['sg_father_summary'].find_one({"name": name})
            if result:
                try:
                    self.db['sg_father_summary'].update({"name": name}, {"$set": {'is_report': 'no'}},
                                                        multi=True, upsert=True)
                except Exception as e:
                    raise Exception("更新存入报告的is_report为no失败{}".format(e))
                else:
                    self.bind_object.logger.info("更新存入报告的is_report为yes成功！")
                return True
            else:
                return False
        except Exception as e:
            raise Exception('存入报告模块在查找sg_father_summary失败{}'.format(e))

    def import_report_summary_result(self, dad_id, mom_id, son_id, father_err_id, file_path, err_min, report_time,
                                     name, report_id, postype='hybrid'):
        """
        用于整合导入报告的汇总记录, 当有father_summary_id的时候该函数只是用于更新
        add by hongdong @ 20171122
        :param dad_id:
        :param mom_id:
        :param son_id:
        :param father_err_id:
        :param file_path:整合的这条记录会存入一下到工作流的output_dir中
        :param err_min:
        :param report_time:
        :param name:
        :param report_id:
        :param postype:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        if not isinstance(report_id, ObjectId):
            if isinstance(report_id, StringTypes):
                report_id = ObjectId(report_id)
            else:
                raise Exception("report_id必须为ObjectId对象或其对应的字符串!")
        try:
            if postype == 'homohybrid':
                result_fq = self.db['sg_father_err_result_homohybrid'].find_one({"father_err_id": father_err_id})
            else:
                result_fq = self.db['sg_father_err_result'].find_one({"father_err_id": father_err_id})
            if "-T" in dad_id:
                dad_id_ = re.match('(WQ.*)-(T.*)', dad_id).group(1)
            elif '-L' in dad_id:
                dad_id_ = re.match('(WQ.*)-(L.*)', dad_id).group(1)
            else:
                dad_id_ = dad_id
            if '-T' in mom_id:
                mom_id_ = re.match('(WQ.*)-(T.*)', mom_id).group(1)
            elif '-L' in mom_id:
                mom_id_ = re.match('(WQ.*)-(L.*)', mom_id).group(1)
            else:
                mom_id_ = mom_id
            son_id_new = son_id.split('_')[0]
            if '-T' in son_id_new:
                son_id_ = re.match('(WQ.*)-(T.*)', son_id_new).group(1)
            elif '-L' in son_id_new:
                son_id_ = re.match('(WQ.*)-(L.*)', son_id_new).group(1)
            else:
                son_id_ = son_id_new
            gest_info = self.gest_info_post(dad_id_, mom_id_, son_id_)
            if not (result_fq and gest_info):
                raise Exception("存入报告的时候father_err_result与客户信息表中没有查找到相关信息！")
            if postype == 'homohybrid':
                # 因为胎儿纯合与杂合综合过滤的结果与单纯胎儿杂合的胎儿浓度等值是一样的，所以胎儿杂合与纯合的过滤结果，
                # 没有重读计算胎儿的浓度等值，所以在进行数据整合的时候，获取不到胎儿浓度等的值，这里做了兼容处理，
                # 首先根据father_err_id sg_father_err_result_homohybrid找到errmin与fatherid值，有了这些值就可以
                # 从sg_father_err中找到胎儿杂合的结果也就是father_err_id，用这个id去sg_father_err_ms_match中找胎儿浓度等值
                father_result = self.db['sg_father_err_result_homohybrid'].find_one({"father_err_id": father_err_id})
                if father_result:
                    err_ = father_result['err_min']
                    father_id_ = father_result['father_id']
                    err_id_result = self.db['sg_father_err'].find_one({"father_id": father_id_, 'err_min': err_})
                    if err_id_result:
                        err_id = err_id_result['_id']
                        ms_match = self.db['sg_father_err_ms_match'].find_one({"father_err_id": err_id})
                    else:
                        ms_match = None
                else:
                    ms_match = None
            else:
                ms_match = self.db['sg_father_err_ms_match'].find_one({"father_err_id": father_err_id})
            if ms_match:
                error = ms_match['error']
                percent = ms_match['percent']
            else:
                percent = ''
                error = ''
            urgency = self.db['sg_file_info'].find_one({"sample_id": dad_id})
            if urgency:
                urgency_degree = urgency['urgence']
            else:
                urgency_degree = ''
            summary_data = {
                "case_name": dad_id.split("-")[0],
                "ask_person": gest_info['customer_name'],  # 申请人
                "mom_name": gest_info['mom_sample_name'],  # 母本名字
                "mom_type": gest_info['mom_sample_type'],  # 模本样本类型
                "mom_id": mom_id,
                "son_id": son_id,
                # "dad_id": re.match("(.*)-(F.*)", dad_id).group(2),
                # "mom_id": re.match("(.*)-(M.*)", mom_id).group(2),
                # "son_id": re.match("(.*)-(S.*)", son_id).group(2),
                "dad_name": gest_info['dad_sample_name'],  # 父本名字
                "dad_type": gest_info['dad_sample_type'],  # 父本样本类型
                "dad_id": dad_id,
                "ask_time": gest_info["entrust_date"],
                "accept_time": gest_info["audit_date"],
                "received_time": gest_info["receive_sample_date"],
                "err_min": str(err_min),  # 允许错配
                "test_pos_n": result_fq['test_pos_n'],  # 测试位点数
                "err_pos_n": result_fq['err_pos_n'],  # 错配为点数
                "err_rate": result_fq['err_rate'],  # 错配率
                "cpi": result_fq['fq'],
                "rcp": "< 0.0001" if str(result_fq['rcp']) == "< 0.01%" else "> 0.9999",
                "show_rcp": result_fq['rcp'],
                "result": "不支持" if str(result_fq['result']) == 'No' else "支持",
                # "ask_time": gest_info_m['ask_time'],  # 委托日期, 申请日期
                # "accept_time": gest_info_m['accept_time'],  # 受理日期
                "quality_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),  # 质控日期,插入该条记录的时间
                "report_time": report_time,  # 报告日期, 存入报告的时间
                "urgency_degree": urgency_degree,  # 加急信息
                "gestation_week": gest_info['pregnant_week'],  # 孕周
                "percent": percent,  # 胎儿浓度
                "error": error,  # 测序错误
                "is_report": "yes",
                "name": name,
                "father_err_id": father_err_id,
                "report_id": report_id,
                "num": "",
                "is_delay": "",
                "delay_reason": "",
                "remark": "",
                "success_material": "",
                "failed_material": "",
                "report_material": "",
                "mom_problem": self.get_sample_problem(mom_id_),
                "mom_frequency": self.get_sample_frequency(mom_id_),
                "dad_problem": self.get_sample_problem(dad_id_),
                "dad_frequency": self.get_sample_frequency(dad_id_),
                "son_problem": self.get_sample_problem(son_id),
                "son_frequency": self.get_sample_frequency(son_id),
                "report_period": "",
                "report_people": self.get_report_people_email(father_err_id),
                "report_check": "",
                "report_first_peopel": "",
                "report_final_people": "",
                'report_result': self.get_report_result(father_err_id)  # 报告组认定的该家系的判定结果
            }
            self.db['sg_father_summary'].insert_one(summary_data)
            with open(file_path + '/summary_data.xls', 'w') as w:
                w.write('家系编号\t父本\t母本\t胎儿\t申请人\t父本姓名\t父本样本类型\t母本姓名\t母本样本类型\t胎儿'
                        '信息位点筛选\t检测为点数\t错配位点数\t错配率\tCPI\tRCP\tRCP显示\t判定\t申请日期\t受理日期\t质'
                        '控日期\t报告日期\t加急信息\t孕周\t胎儿浓度\n')
                w.write(summary_data['case_name'] + '\t' + summary_data['dad_id'] + '\t' +
                        summary_data['mom_id'] + '\t' + summary_data['son_id'] + '\t' +
                        summary_data['ask_person'] + '\t' + summary_data['dad_name'] + '\t' +
                        summary_data['dad_type'] + '\t' + summary_data['mom_name'] + '\t' +
                        summary_data['mom_type'] + '\t' + summary_data['err_min'] + '\t' +
                        summary_data['test_pos_n'] + '\t' + summary_data['err_pos_n'] + '\t' +
                        summary_data['err_rate'] + '\t' + summary_data['cpi'] + '\t' + summary_data['rcp'] + '\t' +
                        summary_data['show_rcp'] + '\t' + summary_data['result'] + '\t' +
                        str(summary_data['ask_time']) + '\t' + str(summary_data['accept_time']) + '\t' +
                        str(summary_data['quality_time']) + '\t' + str(summary_data['report_time']) + '\t' +
                        summary_data['urgency_degree'] + '\t' + summary_data['gestation_week'] + '\t' +
                        summary_data['percent'] + '\n')
        except Exception as e:
            raise Exception("报告模块在在查表整合信息的时候出错！{}".format(e))
        else:
            self.bind_object.logger.info("报告模块查表整合信息成功！")

    def get_sample_frequency(self, sample_id):
        """
        用于获取每个样本出现的次数
        :param sample_id: WQ1234-S-T1_M
        :return:
        """
        if '-S' in sample_id:
            sample_id = str(sample_id).split('_')[0]
            if "-T" in sample_id:
                n = re.match('(WQ.*)-(T.*)', sample_id)
                if n:
                    sample_id = n.group(1)
            elif '-L' in sample_id:
                n = re.match('(WQ.*)-(L.*)', sample_id)
                if n:
                    sample_id = n.group(1)
        return self.db['sg_sample_pt'].find({'sample_id': {'$regex': "^{}.*".format(sample_id)},
                                             "is_analysis": "yes"}).count()

    def get_report_result(self, father_err_id):
        """
        获取报告组对该报告的认定结果
        :return:
        """
        report_result = '--'
        result = self.db['sg_father_err_homohybrid'].find_one({"_id": father_err_id})
        if result:
            father_result = self.db['sg_father'].find_one({"_id": result['father_id']})
            if father_result and 'report_result' in father_result.keys():
                report_result = father_result['report_result']
            else:
                free_result = self.db['sg_free_combine'].find_one({"_id": result['father_id']})
                if free_result and 'report_result' in free_result.keys():
                    report_result = free_result['report_result']
        return report_result

    def get_sample_problem(self, sample_id):
        """
        获取问题样本的具体问题如果有多个异常描述就，分割
        :param sample_id:
        :return:
        """
        if '-S' in sample_id:
            sample_id = str(sample_id).split('_')[0]
            if "-T" in sample_id:
                n = re.match('(WQ.*)-(T.*)', sample_id)
                if n:
                    sample_id = n.group(1)
            elif '-L' in sample_id:
                n = re.match('(WQ.*)-(L.*)', sample_id)
                if n:
                    sample_id = n.group(1)
        self.bind_object.logger.info(sample_id)
        err_list = []
        result = self.db['sg_sample_pt'].find({'sample_id': {'$regex': "^{}.*".format(sample_id)},
                                               "is_analysis": "yes", "is_problem": "1"})  # 这里只查找异常库中的样本
        for m in result:
            if m['remarks']:
                data = phpserialize.loads(str(m['remarks']))
                for key in sorted(data.keys()):
                    if "desc" in data[key].keys():
                        err_list.append(data[key]['desc'])
        dd = ';'.join(err_list) if err_list else '--'
        return dd

    def get_report_people_email(self, father_err_id):
        """
        获取出报告的人的邮箱
        :param father_err_id:
        :return:
        """
        result = self.db['sg_report_result'].find_one({"table_id": str(father_err_id)})
        if result:
            member_email = result['member_email']
        else:
            member_email = '--'
        return member_email

    def update_report_done(self, dad_id, father_err_id, postype='hybrid'):
        """
        更新sg_family中出报告的情况, 首先更新已出报告与未出报告，然后更新这个家系已有的父本是否都已经出报告了
        add by hongdong @ 20171122
        :param dad_id:
        :param father_err_id:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        report_undone = []
        # case_name = re.match("(.*)-F.*", dad_id).group(1)
        case_name = dad_id.split("-")[0]
        sample_id = "-".join(dad_id.split('-')[1:])
        result = self.db['sg_family'].find_one({"case_name": case_name})
        if not result:
            raise Exception("在sg_family中没有查找到{}的case号码！".format(dad_id))
        report_done_ = result['report_done']
        report_undone_ = result['report_undone']
        if sample_id not in report_done_.strip().split(","):
            report_done_ += ',{}'.format(sample_id) if report_done_.strip() != '' else sample_id
            for m in report_undone_.strip().split(','):
                if str(m.strip()) != str(sample_id):
                    report_undone.append(m.strip())
            self.bind_object.logger.info("111{}{}".format(report_done_, report_undone))
            try:
                self.db['sg_family'].update({"case_name": case_name},
                                            {"$set": {'report_done': report_done_,
                                                      "report_undone": ','.join(report_undone)}},
                                            multi=True, upsert=True)
            except Exception as e:
                raise Exception("更新sg_family中的报告状态失败！{}".format(e))
            else:
                self.bind_object.logger.info("更新sg_family中的报告状态成功！")
        else:
            self.bind_object.logger.info("样本{}已经在已出报告中了，不需要再次更新！".format(dad_id))
        if len(report_undone) == 0:
            try:
                self.db['sg_family'].update({"case_name": case_name}, {"$set": {'is_report': '1'}}, multi=True,
                                            upsert=True)
            except Exception as e:
                raise Exception("更新sg_family中所有父本已经出报告状态失败！{}".format(e))
            else:
                self.bind_object.logger.info("更新sg_family中所有父本已经出报告状态成功！")
        if postype == 'hybrid':
            result = self.db['sg_father_err'].find_one({"_id": father_err_id})  # 下载报告的接口运行失败 没有types字段
            if result and result['types'] == "1":
                try:
                    self.db['sg_father'].update({"_id": result['father_id']}, {"$set": {'is_report': '1'}}, multi=True,
                                                upsert=True)
                except Exception as e:
                    raise Exception("更新sg_father中父本已经出报告状态失败！{}".format(e))
                else:
                    self.bind_object.logger.info("更新sg_father中父本已经出报告状态成功！")
        else:
            result = self.db['sg_father_err_homohybrid'].find_one({"_id": father_err_id, 'types': "1"})
            if result:
                self.db['sg_father'].update({"_id": result['father_id']}, {"$set": {'is_report': '1'}}, multi=True,
                                            upsert=True)

    def add_sample_pt(self, samples_list, is_analysis=None):
        """
        用于整合样本的批次信息以及样本的qc信息 add by hongdong @ 20171122
        ot 代表文库质量， ot_dedup代表杂交效率
        position_num 代表每个样本的tab文件行数
        :param samples_list:
        :param is_analysis:  用于没有上机的样本整合信息
        :return:
        """
        info_list = []
        if len(samples_list) == 0:
            return
        for sample_id in samples_list:
            if re.match('(WQ.*)-(S.*)_(M.*)', sample_id):   # 过滤掉S_M的样本
                continue
            sample_result = self.db['sg_sample_pt'].find_one({"sample_id": sample_id})
            if sample_result:
                if sample_result["is_analysis"] == 'no':  # 如果找到该样本没有上机，删除后重新导入
                    self.db['sg_sample_pt'].remove({"sample_id": sample_id})
                else:
                    continue
            sample_info = self.db['sg_file_info'].find_one({"sample_id": sample_id})
            sample_qc = self.db['sg_sample_qc'].find_one({"sample_id": sample_id})
            ot = '/'
            ot_dedup = '/'
            if sample_qc:
                if '-S' in sample_id and sample_info['analysis_type'] == "pt":
                    ot = sample_qc['ot_dedup']  # 文库质量
                    ot_dedup = sample_qc['ot']  # 杂交效率
                else:
                    ot = sample_qc['ot']  # 文库质量
            insert_data = {
                "sample_id": sample_id,   # 样本编号
                "board_num": sample_info['s_sequence_num'] if sample_info and "s_sequence_num" in sample_info.keys()
                else "--",  # 上机测序量
                "num": sample_qc['num'] if sample_qc and 'num' in sample_qc.keys() else "--",  # 数据量
                "dp": sample_qc['dp'] if sample_qc and 'dp' in sample_qc.keys() else "--",   # 测序深度
                "ot": ot,   # 文库质量
                "ot_dedup": ot_dedup,  # 杂交效率
                "0Xcoveragerate": sample_qc['0Xcoveragerate'] if sample_qc and '0Xcoveragerate' in sample_qc.keys()
                else "--",
                "15Xcoveragerate": sample_qc['15Xcoveragerate'] if sample_qc and '15Xcoveragerate' in sample_qc.keys()
                else "--",
                "50Xcoveragerate": sample_qc['50Xcoveragerate'] if sample_qc and '50Xcoveragerate' in sample_qc.keys()
                else "--",
                "100Xcoveragerate": sample_qc['100Xcoveragerate'] if sample_qc and '100Xcoveragerate' in
                                                                                   sample_qc.keys() else "--",
                "500Xcoveragerate": sample_qc['500Xcoveragerate'] if sample_qc and '500Xcoveragerate' in
                                                                                   sample_qc.keys() else "--",
                "1000Xcoveragerate": sample_qc['1000Xcoveragerate'] if sample_qc and '1000Xcoveragerate' in
                                                                                     sample_qc.keys() else "--",
                "num_chrY": sample_qc['num_chrY'] if sample_qc and 'num_chrY' in sample_qc.keys() else "--",  # y位点数
                # num_chrY_int ：排序用
                "num_chrY_int": int(sample_qc['num_chrY']) if sample_qc and 'num_chrY' in sample_qc.keys() else -1,
                "pollution": sample_qc['pollution'] if sample_qc and 'pollution' in sample_qc.keys() else "--",  # 污染程度
                "deal_batch": sample_info['deal_batch'] if sample_info and "deal_batch" in sample_info.keys()
                else "--",  # 处理批次
                "extract_batch": sample_info['extract_batch'] if sample_info and "extract_batch" in sample_info.keys()
                else "--",  # 抽提批次
                "library_batch": sample_info['library_batch'] if sample_info and "library_batch" in sample_info.keys()
                else '--',  # 建库批次
                "board_batch": sample_info['board_batch'] if sample_info and "board_batch" in sample_info.keys()
                else '--',   # 上机板号，上级批次
                "receive_batch": sample_info["receive_batch"] if sample_info and "receive_batch" in sample_info
                else "--",
                "sj_batch_no": sample_info["sj_batch_no"] if sample_info and "sj_batch_no" in sample_info else "--",
                "zj_batch_no": sample_info["zj_batch_no"] if sample_info and "zj_batch_no" in sample_info else "--",
                "is_problem": sample_qc['is_problem'] if sample_qc and "is_problem" in sample_qc.keys() else "--",
                "finish_time": sample_qc['date'] if sample_qc and 'date' in sample_qc.keys() else "--",
                "is_analysis": "no" if is_analysis else 'yes',
                "index": sample_info['index'] if sample_info and "index" in sample_info.keys() else "--",
                "index2": sample_info['index2'] if sample_info and "index2" in sample_info.keys() else "--",
                "position_num": sample_qc['position_num'] if sample_qc and 'position_num' in sample_qc.keys() else "--",
                "remarks": ""
            }
            info_list.append(insert_data)
        if len(info_list) == 0:
            self.bind_object.logger.info("所有的样本的批次信息与qc信息已经存在了--不再进行更新！")
        else:
            try:
                self.db['sg_sample_pt'].insert_many(info_list)
            except Exception as e:
                raise Exception("整合样本的批次信息以及样本的qc信息出错！{}".format(e))
            else:
                self.bind_object.logger.info("整合样本的批次信息以及样本的qc信息ok！")

    def find_father_id(self, case_id, dad_id):
        """
        add by  hoongdong 20171208
        模糊匹配寻找父本id （亲子鉴定自由交互需要查找）
        "sample_id": {"$regex": "WQ1708.*M.*"}
        :param case_id: case号，也就是家系号
        :param dad_id:
        :return: 返回一个父本列表
        """
        dad_list = []
        col = self.db['sg_sample_name']
        sample_dad = col.find({"$or": [{"sample_id": {"$regex": case_id + '.*' + dad_id + '.*'}, "is_tab": "yes"},
                                       {"sample_id": {"$regex": "^{}-(\d)-{}.*".format(case_id, dad_id)},
                                        "is_tab": "yes"}]})
        if sample_dad:
            for i in sample_dad:
                dad_list.append(i["sample_id"])
            return dad_list
        else:
            raise Exception('按照{}与{}进行模糊匹配没有查找到相应的样本，'
                            '请更换查询规则，重新运行!'.format(case_id, dad_id))

    def set_sample_id(self, sample_id):
        m = re.match('(WQ.*)-(.*)', sample_id)
        if m:
            case_name = m.group(1)
            id_name = m.group(2)
        else:
            case_name = "--"
            id_name = "--"
        if self.db['sg_sample_name'].find_one({"sample_id": sample_id}):
            return sample_id
        elif self.db['sg_sample_name'].find_one({"sample_id": {"$regex": "^{}-(\d)-{}.*".format(case_name, id_name)}}):
            return self.db['sg_sample_name'].find_one({"sample_id": {"$regex": "^{}-(\d)-{}.*"
                                                      .format(case_name, id_name)}})['sample_id']
        else:
            return sample_id

    def add_datasplit_info(self, board_batch, fastq_dir):
        """
        用于在拆分计算好了之后，拆分无误后，添加拆分成功的字段，用于避免重复运算
        :param board_batch:
        :param fastq_dir:
        :return:
        """
        try:
            self.db['sg_datasplit'].update({"board_batch": board_batch}, {"$set": {'is_ok': 'yes',
                                                                                   "fastq_dir": fastq_dir}},
                                           multi=True, upsert=True)
        except Exception as e:
            raise Exception("更新拆分结果信息到拆分主表失败！{}".format(e))
        else:
            self.bind_object.logger.info("更新拆分结果信息到拆分主表成功")

    def find_fastq_info(self, board_batch):
        """
        获取对应板子的相关信息，主要是这块板子有没有拆分过以及拆分的状态
        :param board_batch:
        :return:
        """
        result = self.db['sg_datasplit'].find_one({"board_batch": board_batch})
        if "wq_sample" not in result.keys():
            raise Exception("sg_datasplit表格中缺少了wq_sample字段，请在文件检查部分进行确认！")
        if result and "is_ok" in result.keys() and "fastq_dir" in result.keys():
            if result['fastq_dir'] and result['is_ok'] == "yes":
                self.bind_object.logger.info("检测到该板子的数据已经拆分完成了，不再进行拆分分析！")
                return True, result['fastq_dir'], result['wq_sample']
            else:
                self.bind_object.logger.info("检测到该板子的数据没有拆分完成，将进行拆分分析！")
                return False, '', result['wq_sample']
        else:
            self.bind_object.logger.info("检测到该板子的数据没有拆分完成，将进行拆分分析！")
            return False, '', result['wq_sample']

    def add_pt_guest(self, data_list, types):
        """
        type:为insert时执行插入，为update时执行更新
        auther = hongdong
        :param data_list:
        :param types:
        :return:
        """
        collection = self.db["sg_guest_info"]
        if str(types) == "insert":
            try:
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入客户信息表出错:%s" % e)
            else:
                self.bind_object.logger.info("导入客户信息表成功!")
        elif str(types) == 'update':
            for m in data_list:
                try:
                    collection.find_one_and_update({"name": m['case_name'] + m['sample_id']}, {'$set': m}, upsert=True)
                except Exception, e:
                    raise Exception("更新客户信息表出错:%s" % e)
                else:
                    self.bind_object.logger.info("更新客户{}信息表成功!".format(m['case_name'] + m['sample_id']))

    def find_sample_id(self, sample_id):
        """hongdong"""
        collection = self.db["sg_guest_info"]
        result = collection.find_one({"name": sample_id})
        if result:
            self.bind_object.logger.info("数据库中已经有该样本，将进行更新操作！")
            return True
        else:
            return False

    def import_batch_table(self, file_path):
        """
        导入批次表
        内部订单编号,抽提批次,建库批次,加急
        20171021-1_1566677.xls
        :param file_path:
        :return:
        """
        file_name = os.path.basename(file_path)
        data_id = file_name.split('_')[0]
        collection = self.db["sg_sample_batch"]
        try:
            if collection.find_one({'table_name': data_id}):
                collection.remove({'table_name': data_id})
        except Exception as e:
            self.bind_object.logger.error('删除批次表{}信息出错：{}'.format(data_id, e))
            raise Exception('删除批次表{}信息出错：{}'.format(data_id, e))
        data_list = []
        sample_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for line in data:
                line = line.strip().split("\t")
                sample_list.append(line[0])
                insert_data = {
                    'sample_id': line[0],
                    'extract_batch': line[1],
                    'library_batch': line[2],
                    'urgency_degree': line[3],
                    'table_name': data_id
                }
                data_list.append(insert_data)
        try:
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error('导入批次表信息出错：{}'.format(e))
            raise Exception('导入批次表信息出错：{}'.format(e))
        else:
            self.bind_object.logger.info("导入批次表信息成功！")

    def add_to_mongo(self, data_list, table_name, types, query_tuple=()):
        """
        chy
        :param data_list:
        :param types:
        :param table_name:
        :param query_tuple:
        :return:
        """
        collection = self.db[table_name]
        if str(types) == 'insert':
            try:
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception('导入{}出错：{}'.format(table_name, e))
            else:
                self.bind_object.logger.info('导入{}成功！'.format(table_name))
        elif str(types) == 'update':
            if len(query_tuple) < 2:
                raise Exception('update模式必须提供一个长度大于等于2的元组来定义查询比较的方法，'
                                '1定义mongo中的字段，2以后定义插入数据中的字段')
            for m in data_list:
                mongo_spec = query_tuple[0]
                query_value = ''
                for i in range(1, len(query_tuple)):
                    query_value += m[query_tuple[i]]

                try:
                    collection.update({mongo_spec: query_value}, {'$set': m}, upsert=True)
                except Exception, e:
                    raise Exception('更新{}出错:{}'.format(query_value, e))
                else:
                    self.bind_object.logger.info('更新{}成功！'.format(query_value))

    def check_existence(self, table_name, query):
        """
        chy
        :param table_name:
        :param query:
        :return:
        """
        if not isinstance(query, dict):
            raise Exception('调用paternity_test_v2.check_existence方法时query须为字典')
        if len(query) == 0:
            raise Exception('调用paternity_test_v2.check_existence方法时query为空')
        collection = self.db[table_name]
        result = collection.find_one(query)
        if result:
            self.bind_object.logger.info("数据库中已经有该样本！")
            return True
        else:
            return False

    def add_tab(self, table_path, table_name):
        """
        多重、杂捕导tab表
        :param table_path:
        :param table_name:
        :return:
        """
        collection = self.db[table_name]
        count = 0
        data_list_insert = list()
        sample_id = ''                                             # 测试用
        with open(table_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n').split('\t')
                insert_data = {
                    "sample_id": line[0],
                    "chrom": line[1],
                    "pos": line[2],
                    "ref": line[3],
                    "alt": line[4],
                    "dp": line[5],
                    "ref_dp": line[6],
                    "alt_dp": line[7]
                }
                count += 1
                data_list_insert.append(insert_data)
                sample_id = line[0]
        if self.check_existence(table_name, {'sample_id': sample_id}):
            self.bind_object.logger.info("{}样品tab文件在{}表中已经存在，不再重复导入".format(sample_id, table_name))
            return count
        try:
            collection.insert_many(data_list_insert)
        except Exception, e:
            raise Exception('导入{}出错：{}'.format(table_name, e))
        else:
            self.bind_object.logger.info('导入{}成功！'.format(table_name))
        return count

    def update_sample_tab(self, sample_id, board_batch, batch_id, sample_type):
        collection = self.db['sg_sample_name']
        update_data = {
            "board_batch": board_batch,
            "type": sample_type,
            "sample_id": sample_id,
            "batch_id": ObjectId(batch_id),
            "storage_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "is_tab": "yes"
        }
        try:
            collection.update({'sample_id': sample_id}, {'$set': update_data}, upsert=True)
        except Exception, e:
            raise Exception('更新sg_sample_name出错:{}'.format(e))
        else:
            self.bind_object.logger.info('更新sg_sample_name成功！')

    def pt_qc(self, qc_path, sample_id, line_count, board_batch, table_has_content=True):
        """
        杂捕导入qc信息
        :param qc_path:
        :param sample_id:
        :param line_count:
        :param board_batch:
        :param table_has_content:
        :return:
        """
        collection = self.db['sg_sample_qc']
        qc = dict()
        with open(qc_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n').split(':')
                qc[line[0]] = line[1]
        if table_has_content:
            is_problem = '2' if float(qc['dp']) > 10 else '1'
            desc = '/' if float(qc['dp']) > 10 else '样本深度小于阈值'
            sample_type = '1'
            ot = round(int(qc['n_hit']) / float(qc['num']), 4) if float(qc['num']) > 0 else '/'
            ot_dedup = round(int(qc['n_hit_dedup']) / float(qc['n_hit']), 4) if float(qc['n_hit']) > 0 else '/'
        else:
            is_problem = '1'
            desc = 'tab文件大小为0'
            sample_type = '1'
            ot = round(int(qc['n_hit']) / float(qc['num']), 4) if float(qc['num']) > 0 else '/'
            ot_dedup = round(int(qc['n_hit_dedup']) / float(qc['n_hit']), 4) if float(qc['n_hit']) > 0 else '/'
        update_data = {
            'sample_id': sample_id,
            'num': qc['num'].strip(),
            'n_mapped': qc['n_mapped'].strip(),
            'n_hit': qc['n_hit'].strip(),
            'ot': str(ot),
            'dp': qc['dp'].strip(),
            'properly_paired': qc['properly_paired'].strip(),
            'pcr_s': qc['pcr_s'].strip(),
            '0Xcoveragerate': qc['0Xcoveragerate'].strip(),
            '15Xcoveragerate': qc['15Xcoveragerate'].strip(),
            '50Xcoveragerate': qc['50Xcoveragerate'].strip(),
            '100Xcoveragerate': '/',
            'num_chrY': qc['num_chrY'].strip(),
            'is_problem': is_problem,
            'desc': desc,
            'type': sample_type,
            'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'n_dedup': qc['n_dedup'].strip(),
            'n_mapped_dedup': qc['n_mapped_dedup'].strip(),
            'n_hit_dedup': qc['n_hit_dedup'].strip(),
            'ot_dedup': str(ot_dedup),
            'position_num': str(line_count),
            'board_batch': board_batch
        }
        try:
            collection.update({'sample_id': sample_id}, {'$set': update_data}, upsert=True)
        except Exception, e:
            raise Exception('更新sg_sample_qc出错:{}'.format(e))
        else:
            self.bind_object.logger.info('更新sg_sample_qc成功！')
        self.update_sg_family(sample_id)                                 # 更新sg_family中的qc_update_time字段

    def dcpt_qc(self, qc_path, sample_id, line_count, board_batch, table_has_content=True):
        """
        多重导qc信息
        :param qc_path:
        :param sample_id:
        :param line_count:
        :param board_batch:
        :param table_has_content:
        :return:
        """
        collection = self.db['sg_sample_qc']
        qc = dict()
        with open(qc_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n').split(':')
                qc[line[0]] = line[1]
        if table_has_content:
            is_problem = '2' if float(qc['dp1']) > 10 else '1'
            desc = '/' if float(qc['dp1']) > 10 else '样本深度小于阈值'
            sample_type = '2'
            ot = round(int(qc['n_hit']) / float(qc['num']), 4) if float(qc['num']) > 0 else '/'
        else:
            is_problem = '1'
            desc = 'tab文件大小为0'
            sample_type = '2'
            ot = round(int(qc['n_hit']) / float(qc['num']), 4) if float(qc['num']) > 0 else '/'
        file_info = self.db['sg_file_info'].find_one({"sample_id": sample_id})
        update_data = {
            'sample_id': sample_id,
            'num': qc['num'].strip(),
            'n_mapped': qc['n_mapped'].strip(),
            'n_hit': qc['n_hit'].strip(),
            'ot': str(ot),
            'dp': qc['dp1'].strip(),
            'properly_paired': qc['properly_paired'].strip(),
            'pcr_s': qc['pcr_s'].strip(),
            '0Xcoveragerate': qc['0Xcoveragerate'].strip(),
            '15Xcoveragerate': qc['15Xcoveragerate'].strip(),
            '50Xcoveragerate': qc['50Xcoveragerate'].strip(),
            '100Xcoveragerate': qc['100Xcoveragerate'].strip(),
            '500Xcoveragerate': qc['500Xcoveragerate'].strip() if "1000Xcoveragerate" in qc else "/",
            '1000Xcoveragerate': qc['1000Xcoveragerate'].strip() if "1000Xcoveragerate" in qc else "/",
            'num_chrY': qc['num_chrY'].strip(),
            'is_problem': is_problem,
            'desc': desc,
            'type': sample_type,
            'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'n_dedup': '/',
            'n_mapped_dedup': '/',
            'n_hit_dedup': '/',
            'ot_dedup': '/',
            'position_num': str(line_count),
            'board_batch': board_batch,
            "extract_batch": file_info["extract_batch"] if file_info and "extract_batch" in file_info else "--",
            "library_batch": file_info["library_batch"] if file_info and "library_batch" in file_info else "--",
            "deal_batch": file_info["deal_batch"] if file_info and "deal_batch" in file_info else "--",
            "receive_batch": file_info["receive_batch"] if file_info and "receive_batch" in file_info else "--",
            "sj_batch_no": file_info["sj_batch_no"] if file_info and "sj_batch_no" in file_info else "--",
            "zj_batch_no": file_info["zj_batch_no"] if file_info and "zj_batch_no" in file_info else "--",
        }
        try:
            collection.update({'sample_id': sample_id}, {'$set': update_data}, upsert=True)
        except Exception, e:
            raise Exception('更新sg_sample_qc出错:{}'.format(e))
        else:
            self.bind_object.logger.info('更新sg_sample_qc成功！')
        self.update_sg_family(sample_id)                              # 更新sg_family中的qc_update_time字段

    def update_analysis_status(self, batch_id, types, is_update=None):
        """
        更新亲子任务进度
        :param batch_id:
        :param is_update:
        :param types  用于区分是call snp 还是family
        :return:
        """
        if not isinstance(batch_id, ObjectId):
            if isinstance(batch_id, StringTypes):
                batch_id = ObjectId(batch_id)
            else:
                raise Exception("batch_id必须为ObjectId对象或其对应的字符串!")
        if types not in ['snp', 'family']:
            raise Exception("分析类型{}不正确！必须为snp or family".format(types))
        collection = self.db['sg_analysis_status']
        results = collection.find({'_id': batch_id})
        if not results:
            raise Exception('sg_analysis_status查询结果为空')
        try:
            if types == 'snp':
                if is_update == "true":
                    collection.update({'_id': batch_id, 'type': "pt", "is_show": "1"}, {'$inc': {'snp_end_counts': 1}},
                                      upsert=True, multi=True)
                else:
                    return
            elif types == "family":
                result_ = collection.find({'batch_id': batch_id, 'type': "pt"}).sort([("created_ts", -1)])
                if result_:
                    _id = result_[0]['_id']
                    end_count, all_count = self.get_family_ana_num(batch_id)
                    collection.update({'_id': _id}, {'$set': {'family_end_counts': end_count,
                                                              "family_all_count": all_count}}, upsert=True, multi=True)
                else:
                    raise Exception("sg_analysis_status查询结果为空11")
            else:
                pass
        except Exception, e:
            raise Exception('更新sg_analysis_status出错:{}'.format(e))
        else:
            self.bind_object.logger.info('更新sg_analysis_status成功！')

    def get_family_ana_num(self, batch_id):
        """
        用于获取已经该批次中所有家系的个数，以及完成分析的个数
        :param batch_id:
        :return:
        """
        family_all_count = self.db['sg_father'].find({"batch_id": batch_id}).count()
        family_end_counts = self.db['sg_father'].find({"batch_id": batch_id, 'status': {'$in': ['end', 'fai'
                                                                                                       'led']}}).count()
        return family_end_counts, family_all_count

    def update_snp_status(self, batch_id):
        """
        用于更新snp的运行状态，当snp模块不进行分析的时候，状态也要更新为10/10
        :param batch_id:
        :return:
        """
        batch_id = ObjectId(batch_id)
        results = self.db['sg_analysis_status'].find({'batch_id': batch_id, 'type': 'pt'}).sort([("created_ts", -1)])
        if not results:
            raise Exception('查询sg_analysis_status结果为空!')
        _id = results[0]['_id']
        all_counts = results[0]['snp_all_counts']
        try:
            self.db['sg_analysis_status'].update({"_id": _id}, {'$set': {'snp_end_counts': all_counts}}, upsert=True)
        except Exception as e:
            raise Exception('更新sg_analysis_status进度条出错:{}'.format(e))
        else:
            self.bind_object.logger.info('更新sg_analysis_status进度条成功！')

    def update_sg_family(self, sample_id):
        """
        更新sg_family中的qc_update_time字段
        :param sample_id:
        :return:
        """
        update_data = {'qc_update_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
        collection = self.db['sg_family']
        if isinstance(sample_id, StringTypes):
            if sample_id != '':
                casename = sample_id.split('-')[0]
                if self.check_existence('sg_family', {'case_name': casename}):           # 存在casename才更新
                    try:
                        collection.update({'case_name': casename}, {'$set': update_data}, upsert=True)
                    except Exception, e:
                        raise Exception('更新sg_family: {}: qc_update_time出错:{}'.format(casename, e))
                    else:
                        self.bind_object.logger.info('更新sg_family: {}: qc_update_time成功！'.format(casename))
            else:
                raise Exception('sample_id 为空')
        else:
            raise Exception('sample_id类型有误')

    def query_dp(self, sample_id):
        collection = self.db['sg_sample_qc']
        search_result = collection.find_one({"sample_id": sample_id})
        if search_result:
            dp = float(search_result['dp'])
            return dp
        else:
            raise Exception('在sg_sample_qc中未找到{}样品'.format(sample_id))

    def find_split_type(self, batch_id):
        result = self.db['sg_datasplit'].find_one({"_id": ObjectId(batch_id)})
        if result and "split_type" in result.keys():
            if result['split_type'] in ['PE', 'SE']:
                return result['split_type']
            else:
                raise Exception("split_type{},必须为PE 或者 SE".format(result['split_type']))
        else:
            raise Exception("{}在sg_datasplit中没有对应信息！".format(batch_id))

    def analysis_status(self, dad_id, types):
        """
        用于加急信息推荐状态的更新
        :param dad_id:
        :param types: type为2 有家系 但是还没有分析完成，3 有家系，并且分析完成
        :return:
        """
        if types not in ['1', '2', '3']:
            raise Exception("加急分析状态类型不正确！")
        collection = self.db["sg_family"]
        case_name = dad_id.split('-')[0]
        family_result = collection.find_one({"case_name": case_name})
        if family_result and 'analysis_status' in family_result.keys():
            if family_result['analysis_status'] == '3':  # 该family表中已经有能够展示的报告就不进行重复更新了
                return
        try:
            collection.update({"case_name": case_name}, {"$set": {"analysis_status": types}}, multi=True, upsert=True)
        except Exception as e:
            raise Exception("更新sg_family中analysis_status为{}失败！{}".format(types, e))

    def find_ms_id(self, time_past):
        """
        用于查询某一天之后查重主表中的ms_id
        :param time_past:
        :return:
        """
        ms_ids = []
        collection = self.db['sg_father_dedup']
        result = collection.find({"create_time": {"$gt": time_past}})
        if result:
            for i in result:
                if 'm_s_id' in i:
                    ms_ids.append(i['m_s_id'])
            return list(set(ms_ids))
        else:
            raise Exception("查询{}之后的m_s_id无结果".format(time_past))

    def search_father_id(self, mom_id, son_id):
        father_ids = []
        collection = self.db["sg_father"]
        results = collection.find({"mom_id": mom_id, "preg_id": son_id, "status": "end"})
        for i in results:
            if "library_batch" in i and "extract_batch" in i:
                father_ids.append(i["_id"])
        return father_ids

    def get_sample_analysis_type(self, sample_id):
        """
        获取样本对应的建库类型--pt dcpt wqcf
        :param sample_id:
        :return:
        """
        result = self.db["sg_file_info"].find_one({"sample_id": sample_id})
        if result:
            analysis_type = result['analysis_type']
            return analysis_type
        else:
            self.bind_object.logger.info("样本{}没有找到对应的分析类型！--默认设置为wqcf".format(sample_id))
            return "wqcf"

    def find_sample_batch(self, sample_id, sample_type="ALL"):
        """
        用于查找一个样本的处理批次，建库批次与抽提批次，之间是or的关系
        :return:
        sample_type:默认是所有样本包含了父本母本胎儿，但是如果是M的时候， 只保留M的样本，如果是F的时候，保留F的样本，
        如果是S就保留S的样本
        """
        samples = []
        # noinspection PyBroadException
        extract_batch, library_batch, deal_batch, receive_batch, sj_batch_no, zj_batch_no = \
            self.find_father_batch(sample_id)
        extract_batch = extract_batch if extract_batch not in ["", "--"] else "not used to retrieve samples"
        library_batch = library_batch if library_batch not in ["", "--"] else "not used to retrieve samples"
        deal_batch = deal_batch if deal_batch not in ["", "--"] else "not used to retrieve samples"
        receive_batch = receive_batch if receive_batch not in ["", "--"] else "not used to retrieve samples"
        sj_batch_no = sj_batch_no if sj_batch_no not in ["", "--"] else "not used to retrieve samples"
        zj_batch_no = zj_batch_no if zj_batch_no not in ["", "--"] else "not used to retrieve samples"
        batch_result = self.db['sg_sample_batch'].find({"$or": [{"library_batch": library_batch},
                                                                {"extract_batch": extract_batch},
                                                                {"deal_batch": deal_batch},
                                                                {"receive_batch": receive_batch},
                                                                {"sj_batch_no": sj_batch_no},
                                                                {"zj_batch_no": zj_batch_no}]})
        if batch_result:
            for m in batch_result:
                if re.match('WQ.*', m['sample_id']):
                    if sample_type == "ALL":
                        samples.append(m['sample_id'])
                    elif sample_type == 'M' and "-M" in m['sample_id']:
                        samples.append(m['sample_id'])
                    elif sample_type == "F" and '-F' in m['sample_id']:
                        samples.append(m['sample_id'])
                    elif sample_type == "S" and "-S" in m['sample_id']:
                        samples.append(m['sample_id'])
                    else:
                        pass
        return samples

    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        :param id_:
        :return:
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        :param file_path:
        :return:
        """
        if not os.path.exists(file_path):
            return False
        else:
            return True

    def check_is_match(self, file_path):
        """
        :param file_path:
        :return:
        """
        with open(file_path, "r") as r:
            data = r.readlines()[1:]
            temp = data[0].strip().split("\t")
            if float(temp[7]) >= 95:
                return True
            else:
                return False

    def update_db_record(self, collection, query_dict, update_dict, is_show_log="true", upsert=True, multi=True):
        """
        用于更新表格的字段，暂时只是简单的更新一个表的对应字段，后面有需求再进行完善
        :param collection: 集合的名称
        :param query_dict: 查询的字段，可以有多个字段进行联合查询, 是一个字典
        :param update_dict: 表格中要更新进去的字段
        :param is_show_log: 用于设定是否要显示出更新成功的日志信息, 为false的时候就不显示成功的日志信息
        :param upsert: 表中没有查找到对应条件的记录，确认是否要更新 默认为True
        :param multi: 表中按照对应条件查找到多条记录，是否要全部更新 默认为True
        :return:
        """
        try:
            self.db[collection].update(query_dict, {'$set': update_dict}, upsert=upsert, multi=multi)
        except Exception as e:
            raise Exception("更新表格{}失败！{}".format(collection, e))
        else:
            if is_show_log == "true":
                self.bind_object.logger.info("更新{}到表格{}成功！".format(update_dict, collection))

    def col_find_one(self, collection, query_dic):
        """
        用于mongo查询一条记录，查询字段是一个字典,支持一个字段的查询，也支持多个字段的查询
        example: {"batch_id": ObjectID("5a782577a4e1af477e1ac081"), "type": "pt"}
        :param collection:
        :param query_dic:
        :return:
        """
        collection = self.db[collection]
        result = collection.find_one(query_dic)
        # if not result:
        #     print ("没有找到{}表中的{}结果，请检查".format(collection, query_dic.keys()))
        return result

    def col_find(self, collection, query_dic):
        """
        用于mongo查询多条记录，查询字段是一个字典,支持一个字段的查询，也支持多个字段的查询
        example: {"batch_id": ObjectID("5a782577a4e1af477e1ac081"), "type": "pt"}
        :param collection:
        :param query_dic:
        :return:
        """
        collection = self.db[collection]
        return collection.find(query_dic)

    def col_insert_data(self, collection, data_list, is_show_log="true"):
        """
        插入多条记录，data_list样例格式[{"1":"2"},{"2":"3"},{"3": "4"}]
        这里只是抽象出，进行大规模的mongo插表,这里兼容了插入一条记录与多条记录
        :param collection:
        :param data_list:
        :param is_show_log: 用于设定是否要显示出插入成功的日志信息, 为false的时候就不显示成功的日志信息
        :return:
        """
        start = time.time()
        if not data_list:
            raise Exception("列表为空，不能进行后面的插表操作")
        table_id = None
        con = self.db[collection]
        record_num = len(data_list)
        try:
            if record_num > 5000000:
                for i in range(0, record_num, 4000000):
                    temp = data_list[i: i + 4000000]
                    con.insert_many(temp)
            else:
                if record_num >= 2:
                    con.insert_many(data_list)
                else:
                    table_id = con.insert_one(data_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                raise Exception("往{}插入多条记录失败{}".format(collection, e))
            else:
                raise Exception("往{}插入一条记录失败{}".format(collection, e))
        else:
            if is_show_log == 'true':
                if record_num >= 2:
                    print "往{}插入多条记录成功".format(collection)
                    # self.bind_object.logger.info("往{}插入多条记录成功".format(collection))
                else:
                    print "往{}插入一条记录成功".format(collection)
                    # self.bind_object.logger.info("往{}插入一条记录成功".format(collection))
                end = time.time()
                # self.bind_object.logger.info("文件导入mongo中花费的时间:{}".format(end - start))
                print ("文件导入mongo中花费的时间:{}".format(end - start))
            return table_id

    def check_sg_family(self, case_name):
        """
        用于查询sg_family表中特定case_name下的已有样本
        :param case_name:
        :return:
        """
        dad_ids, mom_ids, preg_ids = "", "", ""

        collection = self.db["sg_family"]
        result = collection.find_one({"case_name": case_name})
        if result:
            if "dad_ids" in result:
                dad_ids = result["dad_ids"]
            if "mom_ids" in result:
                mom_ids = result["mom_ids"]
            if "preg_ids" in result:
                preg_ids = result["preg_ids"]
        return dad_ids, mom_ids, preg_ids

    def add_sg_father_msmatch(self, file_path, father_id, son_id):
        """
        MS匹配情况结果查询
        :return:
        """
        father_id = self.check_objectid(father_id)
        with open(file_path, "r") as r:
            data = r.readlines()
            inser_data = {
                "father_id": father_id,
                "son_id": son_id,
                "data_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            if len(data) <= 1:
                inser_data.update({
                    "match_m": []
                })
            else:
                match_m = []
                for line in data[1:]:
                    temp = line.strip().split("\t")
                    match_m.append(temp[0])
                inser_data.update({
                    "match_m": match_m
                })
            self.col_insert_data("sg_father_msmatch", [inser_data])

    def gest_info_post(self, dad_id, mom_id, son_id=None):
        """
        根据父本 母本 请求http://beta.med.com/sync/order2-fetch获取对应的客户信息
        {
            u'pregnant_week': u'6+',
            u'son_sample_no': u'WQ1802490-S',
            u'audit_date': u'',
            u'order_no': u'WQ1802490',
            u'entrust_date': u'',
            u'mom_sample_no': u'WQ1802490-M',
            u'order_product_no': u'WQ1802490',
            u'urgent_status': 1,
            u'dad_sample_name': u'谷天明',
            u'dad_sample_type': u'女方静脉血',
            u'mom_sample_name': u'余楠',
            u'dad_sample_no': u'WQ1802490-F',
            u'receive_sample_date': u'2018-08-16',
            u'mom_sample_type': u'女方静脉血',
            u'customer_name': u'谷天明'
        }
        :param dad_id:
        :param mom_id:
        :param son_id:
        :return:
        """
        # self.bind_object.logger.info("F:{}, M:{}".format(dad_id, mom_id))
        # url = "http://beta.med.com/sync/order2-fetch"
        url = "http://med.i-sanger.cn/sync/order2-fetch"
        key = "pW4mbT92qrB3"
        timestamp = str(int(time.time()))
        sign = hashlib.md5(key + timestamp + key).hexdigest()
        params = {
            "dad_sample_no": dad_id,
            "mom_sample_no": mom_id,
            "timestamp": timestamp,
            "sign": sign
        }
        if son_id:
            params['son_sample_no'] = son_id
        dataa = json.dumps(params)
        headers = {
            'Content-Type': 'application/json',
            'Content-Length': str(len(dataa))}
        d = requests.post(url, data=dataa, headers=headers)
        # noinspection PyBroadException
        try:
            result = json.loads(d.text)
        except:
            self.bind_object.logger.info("请求接口返回的信息为：{}，无法进行json loads, "
                                         "对应的params:{}".format(d.text, params))
            result = {'state': -1}
        # print result
        if result['state'] == -1:
            return {}
        else:
            return result['data']

    def batch_info_post(self, batch_numbers):
        """
        根据批次号（多个值以|分隔） 请求http://med.i-sanger.cn/sync/sample2-fetch 获取满足批次的所有样品的信息

        :param batch_numbers:
        :return:
        """
        url = "http://med.i-sanger.cn/sync/sample2-fetch"
        key = "pW4mbT92qrB3"
        timestamp = str(int(time.time()))
        sign = hashlib.md5(key + timestamp + key).hexdigest()
        if not isinstance(batch_numbers, list) and not isinstance(batch_numbers, tuple) and not isinstance(
                batch_numbers, set):
            raise Exception("输入参数须为list、tuple或set类型")
        params = {
            "batch_no": "|".join(batch_numbers),
            "timestamp": timestamp,
            "sign": sign
        }
        dataa = json.dumps(params)
        headers = {
            'Content-Type': 'application/json',
            'Content-Length': str(len(dataa))}
        d = requests.post(url, data=dataa, headers=headers)
        try:
            result = json.loads(d.text)
        except Exception as e:
            self.bind_object.logger.info("请求接口返回的信息为：{}，无法进行json loads, "
                                         "对应的params:{}".format(d.text, params))
            self.bind_object.logger.info("抓取更多批次信息时出错，错误信息：{}".format(e))
            raise Exception("抓取更多批次信息时出错，错误信息：{}".format(e))

        return result['data']

    def get_ref_list(self, month_back=4):
        """
        根据参数查询指定月份之内的所有样品ID，默认为四个月
        :param month_back:
        :return:
        """
        time_now = datetime.datetime.now()                                # 当前时间戳
        time_past = time_now - relativedelta(months=month_back)           # 4个月之前的时间戳
        time_past_ = time_past.strftime('%Y-%m-%d %H:%M:%S')              # 转化为字符串形式
        ref_list = []
        collection = self.db['sg_sample_name']
        result = collection.find({"storage_time": {"$gt": time_past_}})   # 查询日期比time_past_新的文档
        if result:
            for i in result:
                if 'sample_id' in i:
                    ref_list.append(i['sample_id'])
        return list(set(ref_list))

    @staticmethod
    def s3_enable_flag():
        # 添加项目层面S3启用控制
        # 当self.config.RGW_ENABLE和此方法的返回值都为True时对象存储才启用
        return False

    def add_num_pos_map(self, file_path, father_err_id, father_id):
        """
        将转换后的n值与pos进行匹配关联，为了后面错配位点过滤做参考，该问题也是因为太胎儿是纯合与杂合的综合结果，
        不用进行调试表的导表
        :param file_path:
        :param father_err_id:
        :param father_id:
        :return:
        """
        if not isinstance(father_err_id, ObjectId):
            if isinstance(father_err_id, StringTypes):
                father_err_id = ObjectId(father_err_id)
            else:
                raise Exception("father_err_id必须为ObjectId对象或其对应的字符串!")
        if not isinstance(father_id, ObjectId):
            if isinstance(father_id, StringTypes):
                father_id = ObjectId(father_id)
            else:
                raise Exception("father_id必须为ObjectId对象或其对应的字符串!")
        sg_pt_family_detail = list()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[0] == "new_num":
                    continue
                insert_data = {
                    "father_err_id": father_err_id,
                    'father_id': father_id,
                    "num": str(line[2]),
                    '4num': line[0],
                    'pos': line[1]
                }
                sg_pt_family_detail.append(insert_data)
        if len(sg_pt_family_detail) == 0:
            self.bind_object.logger.info("位点信息表格为空")
            pass
        else:
            try:
                self.db['sg_father_err_num_pos_map'].insert_many(sg_pt_family_detail)
            except Exception as e:
                raise Exception('导入num pos匹配表失败：{}'.format(e))
            else:
                self.bind_object.logger.info("导入num pos匹配表成功!")


if __name__ == "__main__":
    pt = PaternityTestV3(None)
    print pt.gest_info_post("WQ1802334-F2", "WQ1802334-M")
    print "yes"
