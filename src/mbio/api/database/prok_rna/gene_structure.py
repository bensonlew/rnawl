# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import csv
import json
import pandas as pd
import unittest
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from api_base import ApiBase

class GeneStructure(ApiBase):
    def __init__(self, bind_object):
        super(GeneStructure, self).__init__(bind_object)
        self.task_id = self.bind_object.sheet.id

        self.rock_params = json.dumps({"software": "Rockhopper"}, sort_keys=True, separators=(',', ':'))
        self.promote_params = json.dumps({"software": "PromPredict"}, sort_keys=True, separators=(',', ':'))
        self.terminator_params = json.dumps({"software": "TranstermHP"}, sort_keys=True, separators=(',', ':')) # added @20210622

    def add_structure_all(self, rock_dir=None, promote=None, antisense=None):
        """
        基因结构导表函数
        """
        if rock_dir:
            self.add_rock_structure(rock_dir)
        if promote:
            self.add_promote_structure(promote)
        if antisense:
            self.add_antisense(antisense)

    def add_rock_structure(self, rock_dir):
        """
        导入 rockhopper 得到的opera, utr tss 等信息
        """
        for rock_file in ['operon.xls', 'UTR.xls', 'TSS_and_TTS.xls']:
            if os.path.exists(rock_dir + "/" + rock_file):
                pass
            else:
                raise Exception('基因结构对应的结果文件{} 不存在，请检查'.format(rock_dir + "/" + rock_file))

        task_id = self.task_id
        self.remove_table_by_main_record(main_table='sg_structure_operon', task_id=task_id, detail_table=['sg_structure_operon_detail', 'sg_structure_operon_len', 'sg_structure_operon_num'], detail_table_key='operon_id')
        operon_id = self.add_operon()
        self.add_operon_detail(operon_id=operon_id, operon_path=rock_dir + "/" + 'operon.xls')
        self.update_db_record('sg_structure_operon', operon_id, status="end", main_id=operon_id)

        self.remove_table_by_main_record(main_table='sg_structure_utr', task_id=task_id, detail_table=['sg_structure_utr_detail', 'sg_structure_utr_len'], detail_table_key='utr_id')
        utr_id = self.add_utr()
        self.add_utr_detail(utr_id=utr_id, utr_path=rock_dir + "/" + 'UTR.xls')
        self.update_db_record('sg_structure_utr', utr_id, status="end", main_id=utr_id)


        self.remove_table_by_main_record(main_table='sg_structure_tsstts', task_id=task_id, detail_table=['sg_structure_tsstts_detail', 'sg_structure_tsstts_len'], detail_table_key='tsstts_id')
        tsstts_id = self.add_tsstts()
        self.add_tsstts_detail(tsstts_id=tsstts_id, tsstts_path=rock_dir + "/" + 'TSS_and_TTS.xls')
        self.update_db_record('sg_structure_tsstts', tsstts_id, status="end", main_id=tsstts_id)


    def add_promote_structure(self, promote):
        """
        导入 rockhopper 得到的opera, utr tss 等信息
        """
        if os.path.exists(promote):
             pass
        else:
             raise Exception('基因结构对应的结果文件{} 不存在，请检查'.format(promote))

        task_id = self.task_id
        self.remove_table_by_main_record(main_table='sg_structure_promote', task_id=task_id, detail_table=['sg_structure_promote_detail'], detail_table_key='promote_id')
        promote_id = self.add_promote()
        self.add_promote_detail(promote_id=promote_id, promote_path=promote)
        self.update_db_record('sg_structure_promote', promote_id, status="end", main_id=promote_id)
  
    def add_novelgene(self, novelgene):
        """
        导入novelgene 信息
        """
        task_id = self.task_id
        rename_dict = {
            "Gene ID": "gene_id",
            "Location": "location",
            "Start": 'start',
            "End": "end",
            "Strand": "strand",
            "Length": "length"
        }
        self.remove_table_by_main_record(main_table='sg_novelgene', task_id=task_id, detail_table=['sg_novelgene_detail'], detail_table_key='novelgene_id')
        novelgene_id = self.add_novelgene_main()
        self.add_novelgene_detail(novelgene_id=novelgene_id, novelgene_path = novelgene, rename_dict=rename_dict)
        self.update_db_record('sg_novelgene', novelgene_id, status="end", main_id=novelgene_id)

    @report_check
    def add_novelgene_main(self, name=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'novelGene_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': self.promote_params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': 'novelgene预测结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_novelgene']
        novelgene_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_novelgene!")
        return novelgene_id

    def add_novelgene_detail(self, novelgene_id, novelgene_path, rename_dict = {}):
        if not isinstance(novelgene_id, ObjectId):
            if isinstance(novelgene_id, types.StringTypes):
                novelgene_id = ObjectId(novelgene_id)
            else:
                raise Exception('novelgene_id必须为ObjectId对象或其对应的字符串！')

        data_list = []
        df = pd.read_table(novelgene_path, header=0)
        df.rename(columns=rename_dict)
        df['novelgene_id'] = novelgene_id
        data_list = df.to_dict('records')

        if data_list:
            try:
                collection = self.db['sg_novelgene_detail']
                collection.insert_many(data_list)
            except Exception as e:
                raise Exception("导入novelgene详情信息:%s失败！" % novelgene_path)
            else:
                self.bind_object.logger.info("导入novelgene信息:%s成功" % novelgene_path)
    
    
    def add_antisense(self, antisense):
        """
        导入antisense 信息
        """
        task_id = self.task_id
        rename_dict = {
            "GeneID(+)": "gene_id1",
            "Start(+)": "start1",
            "End(+)" : "end1",
            "Description(+)": "description1",
            "GeneID(-)": "gene_id2",
            "Start(-)": "start2",
            "End(-)": "end2",
            "Description(-)": "description2",
            "Overlap_start": "overlap_start",
            "Overlap_end": "overlap_end",
            "Overlap_length": "overlap_length",
            "Type": "type"
        }
        self.remove_table_by_main_record(main_table='sg_antisense', task_id=task_id, detail_table=['sg_antisense_detail'], detail_table_key='antisense_id')
        antisense_id = self.add_antisense_main()
        self.add_antisense_detail(antisense_id=antisense_id, antisense_path = antisense, rename_dict = rename_dict)
        self.update_db_record('sg_antisense', antisense_id, status="end", main_id=antisense_id)

    @report_check
    def add_antisense_main(self, name=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'Antisense_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': self.promote_params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': 'antisense预测结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_structure_antisense']
        antisense_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_antisense!")
        return antisense_id

    @report_check
    def add_antisense_detail(self, antisense_id, antisense_path, rename_dict = {}):
        if not isinstance(antisense_id, ObjectId):
            if isinstance(antisense_id, types.StringTypes):
                antisense_id = ObjectId(antisense_id)
            else:
                raise Exception('antisense_id必须为ObjectId对象或其对应的字符串！')

        df = pd.read_table(antisense_path, header=0)
        df.rename(columns=rename_dict)
        df['antisense_id'] = antisense_id
        data_list = df.to_dict('records')

        if data_list:
            try:
                collection = self.db['sg_structure_antisense_detail']
                collection.insert_many(data_list)
            except Exception as e:
                self.logger.info(e)
                raise Exception("导入antisense详情信息:%s失败！" % antisense_path)
            else:
                self.bind_object.logger.info("导入antisense信息:%s成功" % antisense_path)
    

    @report_check
    def add_promote(self, name=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GeneStructurePromotor_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': self.promote_params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': 'promote预测结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_structure_promote']
        promote_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_structure_promote!")
        return promote_id

    @report_check
    def add_promote_detail(self, promote_id, promote_path):
        if not isinstance(promote_id, ObjectId):
            if isinstance(promote_id, types.StringTypes):
                promote_id = ObjectId(promote_id)
            else:
                raise Exception('promote_id必须为ObjectId对象或其对应的字符串！')
        data_list = []
        with open(promote_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip("\n").split("\t")
                genes_detail = line[6].split(";")
                genes_list = [x.split("(")[0] for x in genes_detail]

                data = [
                    ('promote_id', promote_id),
                    ("gene_id", line[0]),
                    ("location", line[1]),
                    ("upstream", int(line[3])),
                    ("len", int(line[4])),
                    ("seq", line[5]),
                    ("lsp", int(line[6]) if line[6] else None),
                    ("lspe", float(line[7]) if line[7] else None),
                    ("dmaxp", int(line[8]) if line[8] else None),
                    ("dmx", float(line[9]) if line[9] else None),
                    ("dave", float(line[10]) if line[10] else None)
                ]
                data = SON(data)
                data_list.append(data)

        if data_list:
            try:
                collection = self.db['sg_structure_promote_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入promote详情信息:%s失败！" % promote_path)
            else:
                self.bind_object.logger.info("导入promote信息:%s成功" % promote_path)

    def add_terminator_structure(self, terminator):
        """
        导入 TranstermHP 结果
        """
        if os.path.exists(terminator):
             pass
        else:
             raise Exception('基因结构对应的结果文件{} 不存在，请检查'.format(terminator))

        task_id = self.task_id
        self.remove_table_by_main_record(main_table='sg_structure_terminator', task_id=task_id, detail_table=['sg_structure_terminator_detail'], detail_table_key='terminator_id')
        terminator_id = self.add_terminator()
        self.add_terminator_detail(terminator_id=terminator_id, terminator_path=terminator)
        self.update_db_record('sg_structure_terminator', terminator_id, status="end", main_id=terminator_id)


    @report_check
    def add_terminator(self, name=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'Terminator_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': self.terminator_params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': 'terminator预测结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_structure_terminator']
        terminator_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_structure_terminator!")
        return terminator_id

    @report_check
    def add_terminator_detail(self, terminator_id, terminator_path):
        data_list = []
        with open(terminator_path, "r") as f:
            f.readline()
            for line in f:
                line = line.strip("\n").split("\t")
                data = [
                    ('terminator_id', ObjectId(terminator_id)),
                    ("gene_id", line[0]),
                    ("location", line[1]),
                    ("start", int(line[2])),
                    ("end", int(line[3])),
                    ("strand", line[4]),
                    ("5tail", line[5]),
                    ("5stem", line[6]),
                    ("loop", line[7]),
                    ("3stem", line[8]),
                    ("3tail", line[9]),
                    ("score", int(line[10])),
                    ("downstream", int(line[11])),
                ]
                data = SON(data)
                data_list.append(data)

        if data_list:
            try:
                collection = self.db['sg_structure_terminator_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入terminator详情信息:%s失败！" % terminator_path)
            else:
                self.bind_object.logger.info("导入terminator信息:%s成功" % terminator_path)

    @report_check
    def add_operon(self, name=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GeneStructureOpera_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': self.rock_params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': 'operon预测结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_structure_operon']
        operon_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_structure_operon!")
        return operon_id

    def stat_by_range(self, num_list, step):
        """
        :param num_list: 数值列表
        :param step: 统计步长
        根据区间统计数字列表, 最小值为 1 区间为 [1, step], [step+1, step*2]...
]       """
        stat_num = (max(num_list) - 1)/step + 1
        stat_result = [0 for i in range(stat_num)]
        for lens in num_list:
            stat_result[(lens - 1)/step] += 1
        return stat_result

    @report_check
    def add_tsstts_detail(self, tsstts_id, tsstts_path):
        if not isinstance(tsstts_id, ObjectId):
            if isinstance(tsstts_id, types.StringTypes):
                tsstts_id = ObjectId(tsstts_id)
            else:
                raise Exception('tsstts_id必须为ObjectId对象或其对应的字符串！')
        data_list = []
        with open(tsstts_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                genes_detail = line[6].split(";")
                genes_list = [x.split("(")[0] for x in genes_detail]
                try:
                    tts = line[8]
                except:
                    tts = ''
                try:
                    tss = line[5]
                except:
                    tss = ''
                try:
                    css = line[6]
                except:
                    css = ''
                try:
                    cts = line[7]
                except:
                    cts = ''
                data = [
                    ('tsstts_id', tsstts_id),
                    ('gene_id', line[0]),
                    ('gene_name', line[1]),
                    ('description', line[2]),
                    ('location', line[3]),
                    ('tss', tss),
                    ('tts', tts),
                    ('strand', line[4]),
                    ('css', css),
                    ('cts', cts),
                ]
                data = SON(data)
                data_list.append(data)

        if data_list:
            try:
                collection = self.db['sg_structure_tsstts_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入tsstts详情信息:%s失败！" % tsstts_path)
            else:
                self.bind_object.logger.info("导入tsstts信息:%s成功" % tsstts_path)


    @report_check
    def add_operon_detail(self, operon_id, operon_path):
        if not isinstance(operon_id, ObjectId):
            if isinstance(operon_id, types.StringTypes):
                operon_id = ObjectId(operon_id)
            else:
                raise Exception('operon_id必须为ObjectId对象或其对应的字符串！')
        data_list = []
        operon_length_list = []
        operon_genenum_list = []
        with open(operon_path, "r") as f:
            lines = f.readlines()
            operonid = 1
            for line in lines[1:]:
                line = line.strip().split("\t")
                genes_detail = line[6].split(";")
                genes_list = [x.split("(")[0] for x in genes_detail]
                genes_name_list = [x.split("(")[1] for x in genes_detail]
                gene_names = [x.split("|")[0] for x in genes_name_list]
                data = [
                    ('operon_id', operon_id),
                    ('operon', "Operon" + str(operonid).zfill(4)),
                    ('location', line[0]),
                    ('start', line[1]),
                    ('stop', line[2]),
                    ('strand', line[3]),
                    ('gene_num', int(line[4])),
                    ('gene_list', ";".join(genes_list)),
                    ('gene_names', ";".join(gene_names))
                ]
                operon_length_list.append(abs(int(line[2]) - int(line[1])) + 1)
                operon_genenum_list.append(int(line[4]))
                operonid += 1

                data = SON(data)
                data_list.append(data)

        if data_list:
            try:
                collection = self.db['sg_structure_operon_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入operon详情信息:%s失败！" % operon_path)
            else:
                self.bind_object.logger.info("导入operon信息:%s成功" % operon_path)

            step = 2000
            operon_length_stat = self.stat_by_range(operon_length_list, step)
            num_data = [{ "{}~{}".format(i*step +1 , (i+1)*step) : num } for i,num in enumerate(operon_length_stat)]

            data = [
                ('operon_id', operon_id),
                ('step', 2000),
                ('len_data', num_data),
            ]
            data = SON(data)
            try:
                collection = self.db['sg_structure_operon_len']
                collection.insert_one(data)
            except Exception, e:
                self.bind_object.set_error("导入operon长度统计信息出错")
            else:
                self.bind_object.logger.info("导入operon长度统计信息成功")

            step = 1
            operon_num_stat = self.stat_by_range(operon_genenum_list, step)
            num_data = [{ "{}".format(i+1) : num } for i,num in enumerate(operon_num_stat)]
            data = [
                ('operon_id', operon_id),
                ('step', 1),
                ('num_data', num_data),
            ]
            data = SON(data)
            try:
                collection = self.db['sg_structure_operon_num']
                collection.insert_one(data)
            except Exception, e:
                self.bind_object.set_error("导入operon基因数统计信息出错")
            else:
                self.bind_object.logger.info("导入operon基因数统计信息成功")

    @report_check
    def add_utr(self, name=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GeneStructureUtr_'  + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': self.rock_params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': 'utr预测结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_structure_utr']
        utr_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_structure_utr!")
        return utr_id

    @report_check
    def add_utr_detail(self, utr_id, utr_path):
        if not isinstance(utr_id, ObjectId):
            if isinstance(utr_id, types.StringTypes):
                utr_id = ObjectId(utr_id)
            else:
                raise Exception('utr_id必须为ObjectId对象或其对应的字符串！')
        data_list = []
        utr5_length_list = []
        utr3_length_list = []
        utr_genenum_list = []
        with open(utr_path, "r") as f:
            lines = f.readlines()
            utrid = 1
            for line in lines[1:]:
                line = line.strip().split("\t")
                genes_detail = line[6].split(";")
                genes_list = [x.split("(")[0] for x in genes_detail]
                utr_type = "5'UTR"
                if line[6] == "UTR3":
                    utr_type = "3'UTR"
                data = [
                    ('utr_id', utr_id),
                    ('location', line[0]),
                    ('start', line[1]),
                    ('end', line[2]),
                    ('gene_id', line[3]),
                    ('strand', line[5]),
                    ('gene_name', line[4]),
                    ('type', utr_type),
                    ('description', line[7]),
                    ('length', abs(int(line[2]) - int(line[1])) + 1),
                ]
                if line[6] == "UTR5":
                    utr5_length_list.append(abs(int(line[2]) - int(line[1])) + 1)
                else:
                    utr3_length_list.append(abs(int(line[2]) - int(line[1])) + 1)

                data = SON(data)
                data_list.append(data)

        if data_list:
            try:
                collection = self.db['sg_structure_utr_detail']
                collection.insert_many(data_list)
            except Exception, e:
                raise Exception("导入utr详情信息:%s失败！" % utr_path)
            else:
                self.bind_object.logger.info("导入utr信息:%s成功" % utr_path)

            step = 50
            num5_data = []
            if utr5_length_list:
                utr5_length = [x for x in utr5_length_list if x <= 1000]
                utr5_length_large = [x for x in utr5_length_list if x > 1000]
                if utr5_length: # modified by zhangyitong on 20210909, for empty utr5_length
                    utr5_length_stat = self.stat_by_range(utr5_length, step)
                    num5_data += [{"{}~{}".format(i * step + 1, (i + 1) * step): num} for i, num in enumerate(utr5_length_stat)]
                if len(utr5_length_large) > 0:
                    num5_data += [{">1000" : len(utr5_length_large)}]
            # else:
            #     num5_data = []

            num3_data = []
            if utr3_length_list:
                utr3_length = [x for x in utr3_length_list if x<=1000]
                utr3_length_large = [x for x in utr3_length_list if x>1000]
                if utr3_length:
                    utr3_length_stat = self.stat_by_range(utr3_length, step)
                    num3_data += [{ "{}~{}".format(i*step +1 , (i+1)*step) : num } for i,num in enumerate(utr3_length_stat)]
                if len(utr3_length_large) > 0:
                    num3_data += [{">1000" : len(utr3_length_large)}]
            # else:
            #     num3_data = []
            data = [
                ('utr_id', utr_id),
                ('step', 50),
                ('len5_data', num5_data),
                ('len3_data', num3_data),
            ]
            data = SON(data)
            try:
                collection = self.db['sg_structure_utr_len']
                collection.insert_one(data)
            except Exception, e:
                self.bind_object.set_error("导入utr长度统计信息出错")
            else:
                self.bind_object.logger.info("导入utr长度统计信息成功")


    @report_check
    def add_tsstts(self, name=None):
        task_id = self.task_id
        project_sn = self.bind_object.sheet.project_sn

        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GeneStructureTsstts_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'params': self.rock_params,
            # 'result_dir': result_dir,
            'status': 'start',
            'desc': 'tsstts预测结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        collection = self.db['sg_structure_tsstts']
        tsstts_id = collection.insert_one(insert_data).inserted_id
        self.bind_object.logger.info("add sg_structure_tsstts!")
        return tsstts_id

class TestFunction(unittest.TestCase):
    """
    测试导表函数

    """
    def test_mongo(test):
        from mbio.workflows.itraq_and_tmt.itraq_test_api import ItraqApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "prok_rna_srna",
            #+ str(random.randint(1,10000)),
            #"id": "denovo_rna_v2",
            "project_sn": "prok_rna_srna",
            #+ str(random.randint(1,10000)),
            "type": "workflow",
            "name": "prok_rna.prokrna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = ItraqApiWorkflow(wsheet)

        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/api_test/'

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("prok_rna.gene_structure")
        wf.test_api.add_structure_all(rock_dir=test_dir )

if __name__ == '__main__':
    unittest.main()
