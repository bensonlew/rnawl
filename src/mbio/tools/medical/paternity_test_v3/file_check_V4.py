#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
作者：kefei.huang
时间：20171226
# last modify 20180822 by hongyu.chen

检查样本信息

"""
from __future__ import print_function
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from bson.objectid import ObjectId
import os
import re
import codecs
import MySQLdb as mdb
import datetime
import pandas as pd
from collections import defaultdict
import paramiko
import random
import subprocess


class FileCheckV4Agent(Agent):
    def __init__(self, parent):
        super(FileCheckV4Agent, self).__init__(parent)
        self.logger.info("run")
        options = [
            {"name": "up", "type": "infile", "format": "medical.common"},
            {"name": "date", "type": "string"},
            {"name": "flowcell", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": 'indextype', 'type': 'string'}
        ]
        self.add_option(options)
        self.logger.info("options finished")
        self.queue = "gada"

    def check_options(self):
        if not self.option("up"):
            raise OptionError("请输入线上表格")
        if not self.option("date"):
            raise OptionError("请输入结果文件名称，以时间为单位，格式请参照20171023")
        if not self.option("flowcell"):
            raise OptionError("需要输入板号")
        if not self.option("main_id"):
            raise OptionError("需要输入主表ID")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(FileCheckV4Agent, self).end()


class FileCheckV4Tool(Tool):
    def __init__(self, config):
        super(FileCheckV4Tool, self).__init__(config)
        self.up = self.option("up").prop["path"]
        self.date = self.option("date")
        self.flowcell = self.option("flowcell")
        self.update_info = self.option("main_id")
        self.split = self.date + ".split.csv"
        self.errorlog = self.date + ".error.log"
        self.indextype = self.option("indextype")

        self.api_pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')

        self.sample_info = pd.read_excel(self.up).fillna(value="")           # 读取V3版样品信息表
        self.sample_info = self.sample_info.astype(str)
        self.data_convert()

        self.libraryType = {
                                "NIPT文库构建实验流程":"WS",
                                "WQ gDNA建库实验流程_多重":"WQ",
                                "WQ cfDNA建库实验流程_杂捕":"WQ",
                                "FA gDNA建库实验流程_杂捕":"FA",
                                "FA ctDNA建库实验流程_杂捕":"FA",
                                "FA gDNA建库实验流程":"FA",
                                "FA cfDNA建库实验流程_杂捕":"FA",
                                "YCZL 建库实验流程": "YCZL",
                                "QA 建库实验流程":"QA",
                                "CT 建库实验流程":"CT",
                                "HLY 建库实验流程": "HLY,ZL",
                                "RXA 建库实验流程":"RXA",
                                "snp 建库实验流程":"SNP",
                                "XA RNA建库实验流程":"XA",
                                "XA DNA建库实验流程":"XA",
                                "甲基化 建库实验流程":"METH",
                                "转录组文库 建库实验流程":"TRSC",
                                "单细胞DNA 建库实验流程": "scWGS,sc,MOS",
                                "small RNA建库实验流程": "small",
                                "YGY 建库实验流程": "YG",
                                "全外白细胞建库": "QBC",
                                "WQCF gDNA建库实验流程": "WQ",
                                "WQCF cfDNA建库实验流程": "WQ",
                                "XLT gDNA建库实验流程": "XLT"
                            }
        self.sampletype = {
                                "白细胞": "BXB",
                                "全血": "QX",
                                "血浆": "XJ",
                                "蜡块": "SL",
                                "石蜡": "SL",
                                "石蜡切片": "SL",
                                "石蜡切片(白片)": "SL,FP",
                                "胸腹水": "XS",
                                "手术标本": "XZ,SS",
                                "穿刺标本": "XZ",
                                "穿刺样本": "XZ,CC",  # add by HD 20180131 添加CC
                                "组织标本": "XZ",
                                "新鲜组织": "XZ",
                                "蜡卷": "SL,LJ",
                                "精斑": "JB",
                                "亲子父本全血":"QQ",
                                "指甲": "ZJ"
                           }
        self.analysistype = {
                                "NIPT文库构建实验流程": "nipt",
                                "WQ gDNA建库实验流程_多重": "dcpt",
                                "WQ cfDNA建库实验流程_杂捕": "pt",
                                "FA gDNA建库实验流程_杂捕": "ctdna",
                                "FA ctDNA建库实验流程_杂捕": "ctdna",
                                "FA gDNA建库实验流程": "ctdna",
                                "FA cfDNA建库实验流程_杂捕": "ctdna",
                                "YCZL 建库实验流程": "genetic_tumor",
                                "QA 建库实验流程": "QA_str",
                                "CT 建库实验流程": "ct_str",
                                "HLY 建库实验流程": "dc_str",
                                "XA RNA建库实验流程": "blood_cancer",
                                "XA DNA建库实验流程": "blood_cancer",
                                "甲基化 建库实验流程": "",
                                "RXA 建库实验流程": "",
                                "snp 建库实验流程": "",
                                "转录组文库 建库实验流程": "",
                                "单细胞DNA 建库实验流程": "",
                                "small RNA建库实验流程": "",
                                "YGY 建库实验流程": "",
                                "全外白细胞建库": "",
                                "WQCF gDNA建库实验流程": "wqcf",
                                "WQCF cfDNA建库实验流程": "wqcf",
                                "XLT gDNA建库实验流程": "XLT",
                                "CJSwailaidata": 'pass'
                             }
        print("initialization finished")

    def run(self):
        super(FileCheckV4Tool, self).run()
        self.logger.info("tools start run")

        Checkbarcode_Error = self.Checkbarcode()
        NamePattern_Error = self.NamePattern()
        self.GenerateSplitTable()
        ###判断是否可以导表
        self.logger.info(Checkbarcode_Error)
        self.logger.info(NamePattern_Error)
        Error_sum = Checkbarcode_Error + NamePattern_Error
        tempout = "Checkbarcode_Error = {}\nNamePattern_Error = {}\n".format(Checkbarcode_Error, NamePattern_Error)
        print (tempout)
        MAINID = self.update_info
        if Error_sum == 0:
            self.dump_error("end")
            self.dump_experiment_batch()
            self.dump_sample_info()
            self.dump_split()
            self.dump_sg_customer(MAINID)
            self.add_file_family(MAINID)
        else:
            self.dump_error("failed")
        self.trigger_file_check_api()                         # 触发线下文件检查
        self.end()

    def data_convert(self):
        """
        处理各种转换（适应样本信息表的变化）
        1、将加急状态中的加急24小时、加急48小时转换为加急1天、加急2天
        2、更改列名： 接收批次 -> 处理批次
        :return:
        """
        for i in range(0, len(self.sample_info)):
            if self.sample_info.loc[i, u'加急状态'] == u"加急24小时":
                self.sample_info.loc[i, u'加急状态'] = u"加急1天"
            elif self.sample_info.loc[i, u'加急状态'] == u"加急48小时":
                self.sample_info.loc[i, u'加急状态'] = u"加急2天"
            elif self.sample_info.loc[i, u'加急状态'] == u"1" or self.sample_info.loc[i, u'加急状态'] == u"1.0":
                self.sample_info.loc[i, u'加急状态'] = u"不加急"
            elif self.sample_info.loc[i, u'加急状态'] == u"2" or self.sample_info.loc[i, u'加急状态'] == u"2.0":
                self.sample_info.loc[i, u'加急状态'] = u"加急3天"
            elif self.sample_info.loc[i, u'加急状态'] == u"3" or self.sample_info.loc[i, u'加急状态'] == u"3.0":
                self.sample_info.loc[i, u'加急状态'] = u"加急2天"
            elif self.sample_info.loc[i, u'加急状态'] == u"4" or self.sample_info.loc[i, u'加急状态'] == u"4.0":
                self.sample_info.loc[i, u'加急状态'] = u"加急1天"

    def Checkbarcode(self):
        """
            这一部分要处理简单的查重
            1. 首先检查barcode列是不是对的。
            2. 其次检查里面有没有重复
            3. 如果有重复，就要检查样品名是不是一致，如果一致就跳过。否则就要报错
            4. 如果是双index的时候，这里要增加barcode2的判断，但是具体怎么判断还需要和试验端确定。
        """
        error_count = 0
        barcode_set = set()
        sample_set = set()
        ERROR = codecs.open(self.errorlog, "w", "utf-8")

        for i in range(0, len(self.sample_info)):
            barcode = self.sample_info.loc[i, u'index序列'].encode("utf-8")
            sample_id = self.sample_info.loc[i, u'样本编号'].encode("utf-8")
            if re.match(r'^[ATCG]+$', barcode) is None:
                errorlog = "barcode consist of ATCG only,so barcode column may be error,check your file"
                print (errorlog, file=ERROR)
                error_count += 1
            # 检查barcode有没有重复
            if barcode not in barcode_set:
                barcode_set.add(barcode)
            else:
                errorlog = "barcode {} repeat\n".format(barcode)
                error_count += 1
                ERROR.write(errorlog)
            # 检查样品名有没有重复
            if sample_id not in sample_set:
                sample_set.add(sample_id)
            else:
                errorlog = "sample {} repeat\n".format(sample_id)
                error_count += 1
                ERROR.write(errorlog)
        ERROR.close()
        return error_count

    def NamePattern(self):
        """
            这一部分就是处理杂捕，多重之类的命名规则。
            字典放置在程序的最上层。考虑到以后很多东西都会重新命名
        """
        pt_db = Config().get_mongo_client(mtype="pt_v3")[Config().get_mongo_dbname("pt_v3")]
        nipt_db = Config().get_mongo_client(mtype="nipt_v2")[Config().get_mongo_dbname("nipt_v2")]
        error_count = 0
        ERROR = codecs.open(self.errorlog, "a", "utf-8")

        for i in range(0, len(self.sample_info)):

            sample_id = self.sample_info.loc[i, u'样本编号'].encode("utf-8")
            jk_type = self.sample_info.loc[i, u'建库类型'].encode("utf-8")
            use_type = self.sample_info.loc[i, u'分析类型'].encode("utf-8")

            ###检查样品名称中是否有多余的"."符号、空格符号
            if '.' in sample_id:
                error_count += 1
                errorlog = '{} has . in name!'.format(sample_id)
                print(errorlog, file=ERROR)
            if ' ' in sample_id:
                error_count += 1
                errorlog = '{} has space symbol in name!'.format(sample_id)
                print(errorlog, file=ERROR)

            # 检查测试样品是否以TEST起始
            if use_type == "测试" and not sample_id.startswith("TEST"):
                error_count += 1
                errorlog = '{}: test sample must start with "TEST"'.format(sample_id)
                print(errorlog, file=ERROR)
            if use_type == "生产" and sample_id.startswith("TEST"):
                error_count += 1
                errorlog = '{}: product sample must NOT start with "TEST"'.format(sample_id)
                print(errorlog, file=ERROR)

            # 检查亲子样本是否有亲本信息
            if "WQ" in sample_id and not re.search(r'-S|-M|-F', sample_id):
                error_count += 1
                errorlog = '{}：亲子样品名中没有-F、-M或-S！'.format(sample_id)
                print(errorlog, file=ERROR)

            ###检查样品是否已经分析过
            if "WQ" in sample_id:
                tab_result = pt_db["sg_sample_tab"].find_one({"sample_id": sample_id})
                tab_result_problem = pt_db["sg_sample_tab_problem"].find_one({"sample_id": sample_id})
                if tab_result or tab_result_problem:
                    error_count += 1
                    errorlog = '{}：样品已分析过，请检查是否命名错误！'.format(sample_id)
                    print(errorlog, file=ERROR)
            elif "WS" in sample_id:
                bed_result = nipt_db["sg_bed"].find_one({"sample_id": sample_id})
                if bed_result:
                    error_count += 1
                    errorlog = '{}：样品已分析过，请检查是否命名错误！'.format(sample_id)
                    print(errorlog, file=ERROR)
        ERROR.close()
        return error_count

    def GenerateSplitTable(self):
        if self.indextype == 'single':
            out_title = ["sample_name", "analysis_type", "emergency", "department", "M reads", "index", "length", "zj_group"]
        else:
            out_title = ["sample_name", "analysis_type", "emergency", "department", "M reads", "index", "index2", "length", "zj_group"]
        for i in range(0, len(self.sample_info)):
            self.sample_info.loc[i, u"任务类型"] = self.analysistype[self.sample_info.loc[i, u"建库类型"].encode("utf-8")].decode("utf-8")
        self.sample_info[u"产品线"] = u"MED"
        if self.indextype == 'single':
            columns = [u"样本编号", u"任务类型", u"加急状态", u"产品线", u"上机数据量（M）", u"index序列", u"插入片段（bp）", u"杂交组号"]
        else:
            columns = [u"样本编号", u"任务类型", u"加急状态", u"产品线", u"上机数据量（M）", u"index序列", u"index序列2", u"插入片段（bp）", u"杂交组号"]
        self.sample_info.to_csv(self.split, columns=columns, header=out_title, index=False, encoding='utf-8')

    def dump_error(self, status):
        with codecs.open(self.errorlog, "r", "utf-8") as ERROR:
            errorlog = ERROR.readlines()
            errorlog = "".join(errorlog)

        db = Config().get_mongo_client(mtype="pt_v3")[Config().get_mongo_dbname("pt_v3")]

        insert_data = {
            "log": errorlog
        }
        db['sg_file_check'].update({"_id": ObjectId(self.option("main_id"))}, {"$set": insert_data}, upsert=True)
        if status == 'failed':
            self.set_error("error")  # 这里error不能变，是用于后面更新的
            raise Exception("文件检查中发现错误，流程终止！")
        print ("add sg_file_check!")

    def dump_experiment_batch(self):
        db = Config().get_mongo_client(mtype="pt_v3")[Config().get_mongo_dbname("pt_v3")]
        collection = db['sg_sample_batch']
        # 选取样本信息表中的批次信息
        col_n = [u'样本编号', u'抽提批次', u'建库批次', u'处理批次', u"上机批次", u"接收批次", u"杂交批次"]
        batch_info = pd.DataFrame(self.sample_info, columns=col_n)
        # 根据样本信息表中的批次信息 向医检系统查询批次相同的样本的信息
        batch_numbers = list()
        for col in [u'抽提批次', u'建库批次', u'处理批次', u"上机批次", u"接收批次", u"杂交批次"]:
            batch_numbers.extend(batch_info[col].values.tolist())
        batch_numbers = [i.encode("utf-8") for i in batch_numbers]
        batch_numbers = set(batch_numbers)

        query_results = self.api_pt.batch_info_post(batch_numbers)
        mapping = {
            u'sample_no': u'样本编号',
            u'ct_batch_no': u'抽提批次',
            u'jk_batch_no': u'建库批次',
            u'cl_batch_no': u'处理批次',
            u"sj_batch_no": u"上机批次",
            u"js_batch_no": u"接收批次",
            u"zj_batch_no": u"杂交批次"
        }
        query_info = pd.DataFrame.from_records(query_results).rename(columns=mapping).fillna(value="")
        # 合并样本信息表中的批次信息和查询到的批次信息，去重
        batch_info = pd.concat([batch_info, query_info], axis=0, join="outer",
                               ignore_index=True).drop_duplicates().reset_index(drop=True)
        batch_info.to_excel("batch_info.xlsx", index=False)
        # 导表
        for i in range(0, len(batch_info)):
            sample_id = batch_info.loc[i, u'样本编号'].encode("utf-8")
            if "WQ" in sample_id or "WS" in sample_id:
                insert_data = {
                    'sample_id': sample_id,
                    'extract_batch': batch_info.loc[i, u'抽提批次'].encode("utf-8"),
                    'library_batch': batch_info.loc[i, u'建库批次'].encode("utf-8"),
                    'deal_batch': batch_info.loc[i, u'处理批次'].encode("utf-8"),
                    'receive_batch': batch_info.loc[i, u'接收批次'].encode("utf-8"),
                    'sj_batch_no': batch_info.loc[i, u'上机批次'].encode("utf-8"),
                    'zj_batch_no': batch_info.loc[i, u'杂交批次'].encode("utf-8"),
                    'insert_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                }
                collection.update({"sample_id": sample_id}, {'$set': insert_data}, upsert=True, multi=True)

    def dump_sample_info(self):

        db = Config().get_mongo_client(mtype="pt_v3")[Config().get_mongo_dbname("pt_v3")]
        collection = db['sg_sample_info']
        count = 1
        for i in range(0, len(self.sample_info)):
            sample_id = self.sample_info.loc[i, u'样本编号'].encode("utf-8")
            case_name = sample_id.split("-")[0]
            sample = "-".join(sample_id.split("-")[1:]) if len(sample_id.split("-")) > 1 else ""  # 亲本
            insert_table = {
                'number': count,
                'sample_accept_time': self.sample_info.loc[i, u'收样日期'].encode("utf-8"),
                'library_type': self.sample_info.loc[i, u'建库类型'].encode("utf-8"),
                'analysis_type': self.sample_info.loc[i, u'任务类型'].encode("utf-8"),
                'product_type': self.sample_info.loc[i, u'产品线'].encode("utf-8"),
                'case_name': case_name,
                'sample_id': sample_id,   # modified by hd 20180111 统一sample_id 与余果姐后面家系组建衔接起来
                'sample': sample,
                'sample_type': self.sample_info.loc[i, u'样本类型'].encode("utf-8"),
                's_sequence_num': self.sample_info.loc[i, u'上机数据量（M）'].encode("utf-8"),
                'index': self.sample_info.loc[i, u'index序列'].encode("utf-8"),
                'sequence_size': self.sample_info.loc[i, u'插入片段（bp）'].encode("utf-8"),
                'use_type': self.sample_info.loc[i, u'分析类型'].encode("utf-8"),         # 生产|测试
                'board_batch': self.flowcell,
                'pattern': sample_id,
                'index2': self.sample_info.loc[i, u'index序列2'].encode("utf-8") if self.indextype == 'double' else '--',
            }
            count += 1
            collection.update({'sample_id': sample_id}, {'$set': insert_table}, upsert=True)

    def dump_split(self):

        db = Config().get_mongo_client(mtype="pt_v3")[Config().get_mongo_dbname("pt_v3")]
        collection = db['sg_datasplit']

        WQ_sample = []
        WS_sample = []

        for i in range(0, len(self.sample_info)):
            sample_id = self.sample_info.loc[i, u'样本编号'].encode("utf-8")
            if "WS" in sample_id:
                WS_sample.append(sample_id)
            elif "WQ" in sample_id:
                WQ_sample.append(sample_id)

        WQ_num = len(WQ_sample)
        WS_num = len(WS_sample)

        split_abs = os.getcwd() + "/" + self.split
        insert_table = {
            'board_batch': self.flowcell,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),  # modified by hongdong
            #  20171230 create_time -> created_ts
            "wq_sample": WQ_sample,
            "ws_sample": WS_sample,
            "wq_num": WQ_num,
            "ws_num": WS_num,
            "start_time": "",
            "status": "produced" if re.match(r".*_wailaidata_.*", self.flowcell) else "unproduced",
            "desc": "",
            "end_time": "",
            'split_type': 'PE' if re.match(r".*_wailaidata_.*", self.flowcell) else "",
            "split_info_path": split_abs,
            'board_path': self.get_board_path(self.flowcell),
            'indextype': self.indextype
        }
        collection.update({'board_batch': self.flowcell}, {'$set': insert_table}, upsert=True)

    def get_board_path(self, board_batch):
        """
        用于添加具体的板子所在路径，目前是nextseq 与nextseq1  add by hongdong 20180108
        :param board_batch:
        :return:
        """
        if os.path.exists("/mnt/ilustre/upload/nextseq1/" + board_batch):
            board_path = "/mnt/ilustre/upload/nextseq1/" + board_batch
        elif os.path.exists("/mnt/ilustre/upload/nextseq/" + board_batch):
            board_path = "/mnt/ilustre/upload/nextseq/" + board_batch
        elif os.path.exists("/mnt/ilustre/users/yixuezhuanhua/raw-data/" + board_batch):
            board_path = "/mnt/ilustre/users/yixuezhuanhua/raw-data/" + board_batch
        elif os.path.exists("/mnt/clustre/upload/nextseq1/" + board_batch):
            board_path = "/mnt/clustre/upload/nextseq1/" + board_batch
        elif os.path.exists("/mnt/clustre/upload/nextseq/" + board_batch):
            board_path = "/mnt/clustre/upload/nextseq/" + board_batch
        elif os.path.exists("/mnt/clustre/users/yixuezhuanhua/raw-data/" + board_batch):
            board_path = "/mnt/clustre/users/yixuezhuanhua/raw-data/" + board_batch
        else:
            # add by hd 20200810 增加外来数据样本路径信息
            if re.match(r".*_wailaidata_.*", self.flowcell):
                board_path = "/mnt/ilustre/users/yixuezhuanhua/raw-data/" + board_batch
            else:
                board_path = ''
                self.set_error("{}不在下机路径中".format(board_batch))
        self.logger.info("board_path{}".format(board_path))
        return board_path

    def dump_sg_customer(self, MAINID):
        pt_db = Config().get_mongo_client(mtype="pt_v3")[Config().get_mongo_dbname("pt_v3")]
        pt_collection = pt_db['sg_file_info']
        nipt_count = 1    # modified by hongdong 20171229 747-749 删除每行后的； 因为python 语法不对
        pt_count = 1
        count = 0

        for i in range(0, len(self.sample_info)):
            sample_id = self.sample_info.loc[i, u'样本编号'].encode("utf-8")

            if "WQ" in sample_id:
                count = pt_count
                pt_count += 1
            elif "WS" in sample_id:
                count = nipt_count
                nipt_count += 1

            insert_table = {
                "number": count,
                "sample_id": sample_id,
                "extract_batch": self.sample_info.loc[i, u'抽提批次'].encode("utf-8"),
                "library_batch": self.sample_info.loc[i, u'建库批次'].encode("utf-8"),
                "deal_batch": self.sample_info.loc[i, u'处理批次'].encode("utf-8"),
                'receive_batch': self.sample_info.loc[i, u'接收批次'].encode("utf-8"),
                'sj_batch_no': self.sample_info.loc[i, u'上机批次'].encode("utf-8"),
                'zj_batch_no': self.sample_info.loc[i, u'杂交批次'].encode("utf-8"),
                "library_type": self.sample_info.loc[i, u'建库类型'].encode("utf-8"),
                "index": self.sample_info.loc[i, u'index序列'].encode("utf-8"),
                "s_sequence_num": self.sample_info.loc[i, u'上机数据量（M）'].encode("utf-8"),
                "sequence_size": self.sample_info.loc[i, u'插入片段（bp）'].encode("utf-8"),
                "urgence": self.sample_info.loc[i, u'加急状态'].encode("utf-8"),
                "use_type": self.sample_info.loc[i, u'分析类型'].encode("utf-8"),         # 生产|测试
                "upload_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "board_batch": self.flowcell,
                "check_id": ObjectId(MAINID),
                'analysis_type': self.sample_info.loc[i, u'任务类型'].encode("utf-8"),
                'index2': self.sample_info.loc[i, u'index序列2'].encode("utf-8") if self.indextype == 'double' else '--',
            }

            if "WQ" in sample_id:
                pt_collection.update({'sample_id': sample_id}, {'$set': insert_table}, upsert=True)
            elif "WS" in sample_id:
                pt_collection.update({'sample_id': sample_id}, {'$set': insert_table}, upsert=True)

    def add_file_family(self, MAINID):
        """
        用于导入下机信息查看的记录（V3）
        :param MAINID:
        :return:
        """
        cases = defaultdict(list)
        for i in range(0, len(self.sample_info)):
            sample_id = self.sample_info.loc[i, u'样本编号'].encode("utf-8")
            case_name = sample_id.split("-")[0]
            parent = "-".join(sample_id.split("-")[1:]) if len(sample_id.split("-")) > 1 else ""  # 亲本
            cases[case_name].append(parent)

        data_list_insert = list()
        api = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        for case_name in cases:
            results = list()                           # 保存此case_name下缺少的样本
            dad_ids, mom_ids, preg_ids = api.check_sg_family(case_name)
            samples_in_board = ",".join(cases[case_name])
            if "F" in samples_in_board or dad_ids:
                self.logger.info("{}: F样品存在".format(case_name))
            else:
                results.append("F")
            if "M" in samples_in_board or mom_ids:
                self.logger.info("{}: M样品存在".format(case_name))
            else:
                results.append("M")
            if "S" in samples_in_board or preg_ids:
                self.logger.info("{}: S样品存在".format(case_name))
            else:
                results.append("S")
            not_ok_sample = ",".join(results)
            for parent in cases[case_name]:
                sample_id = case_name + "-" + parent
                insert_data = {
                    "board_batch": self.flowcell,
                    "sample_id": sample_id,
                    "case_id": case_name,
                    "not_ok_sample": not_ok_sample,
                    "file_id": MAINID
                }
                data_list_insert.append(insert_data)

        api.add_to_mongo(data_list_insert, "sg_file_family", "insert")             # 导表

    def trigger_file_check_api(self):
        # 将样本信息表用SCP传到线下并触发线下文件检查API
        ip = "192.168.12.46"
        api = "http://192.168.12.46:5000/file_check"
        user = "hongyu.chen"
        key = "mj123456"
        port = 22
        remote_dir = "/mnt/ilustre/users/hongyu.chen/web/uploads"
        random_num = datetime.datetime.now().strftime("%H%M%S") + str(random.randint(1000, 10000))
        filename = "{}_{}.xlsx".format(self.flowcell, random_num)
        remote_file = os.path.join(remote_dir, filename)
        # SFTP
        try:
            scp = paramiko.Transport(ip, port)
            scp.connect(username=user, password=key)
            sftp = paramiko.SFTPClient.from_transport(scp)
            sftp.put(self.up, remote_file)                    # 将样本信息表用SCP传到线下
            scp.close()
        except Exception, e:
            self.logger.info(e)
        # trigger file check api
        self.logger.info("开始触发线下文件检查")
        cmd = 'curl {} -d "sample_info={}&board_name={}" -X POST -v'.format(api, remote_file, self.flowcell)
        self.logger.info("运行： {}".format(cmd))
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        return_code = p.returncode
        self.logger.info(stdout)
        if return_code != 0:
            self.logger.info(stderr)
        else:
            self.logger.info("线下文件检查成功触发")
