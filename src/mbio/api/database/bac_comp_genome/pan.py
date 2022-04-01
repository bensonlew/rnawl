# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20190902

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
import pandas as pd
import numpy as np


class Pan(Base):
    def __init__(self, bind_object):
        super(Pan, self).__init__(bind_object)
        self._project_type = "bac_comparative"
        #self.id = 'tsg_123'
        #self.project_sn = '188_5b5acb3018'

    @report_check
    def add_pan(self, params=None, name=None):
        """
        pan_genomes的主表；params中记录了group_detail等字段
        :param params:
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        # task_id = self.id
        # project_sn = self.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "泛基因组分析主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Pan_Origin",
        }
        collection = self.db["pan"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    def add_pan_detail(self, inserted_id, cluster_data, anno_file, cluster_path=None):
        """
        导表的cluster表
        :param inserted_id: 主表id
        :param cluster: cluster表
        :param stat: 变异分析的统计表
        :param anno_file: 注释结果表
        :param cluster_path: s3路径
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        self.bind_object.logger.info("cluster_data:path %s" % (cluster_data))
        cluster_dict = {}
        # with open(stat, 'r') as f2:
        #     variation_lines = f2.readlines()
        #     for line in variation_lines[1:]:
        #         line = line.strip().split("\t")
        #         cluster = line[0]
        #         dn = int(line[2])
        #         ds = int(line[3])
        #         if ds == 0 or dn ==0:
        #             stat = 0
        #         else:
        #             stat = float(dn) / ds
        #         cluster_dict[cluster] = stat
        cog_cluster = {}
        with open(anno_file, 'r') as anno:
            lins = anno.readlines()
            for lin in lins:
                new_li = lin.strip().split("\t")
                cog_cluster[new_li[0]] = new_li

        with open(cluster_data, "r") as f:
            data_list = []
            header = f.readline().strip().split("\t")
            sample_list = header[3:]
            core = "no"
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "cluster_id": line[0],
                    "sample_num": int(line[1]),
                    "gene_num": int(line[2]),
                    "pan_id": new_inserted_id,
                }
                if int(line[1]) == len(sample_list):
                    core = "yes"
                for i in range(3,len(sample_list)+3):
                    sample = sample_list[i-3].replace('.', '_')
                    new_gene_list = []
                    if line[i] != '-':
                        gene_list = line[i].split(",")
                        for j in gene_list:
                            if j.split("|")[1] not in new_gene_list:
                                new_gene_list.append(j.split("|")[1])
                        data[sample] = ",".join(new_gene_list)
                    else:
                        data[sample] = '-'
                if line[0] in cluster_dict.keys():
                    data["dn_ds"] = cluster_dict[line[0]]
                else:
                    data["dn_ds"] = 0
                if line[0] in cog_cluster.keys():
                    anno_list = cog_cluster[line[0]]
                    data['cog_id'] = anno_list[2]
                    data['cog_des']= anno_list[3]
                    data['ko_id'] = anno_list[4]
                    data['ko_des'] = anno_list[5]
                else:
                    data['cog_id'] = "-"
                    data['cog_des']= "-"
                    data['ko_id'] = "-"
                    data['ko_des'] = "-"
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["pan_detail"]
            main_collection = self.db["pan"]
            collection.insert_many(data_list)
            self.bind_object.logger.info("开始更新主表的cluster路径")
            main_collection.update({'_id': new_inserted_id},{'$set':{'cluster_path': cluster_path, 'spe_list': ','.join(sample_list), 'main_id':new_inserted_id,
                                                                     "core": core}})
        except Exception, e:
            self.bind_object.logger.info("导入cluster结果表出错:%s" % (e))

    def add_pan_genome(self, inserted_id, type, genomes, genome_data):
        """
        导表计算pan、core、newgenes的数据和公式
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        data = pd.read_table(genome_data, sep="\t", header=0)###画图的10条数据
        data1 = data.T
        sample_list = list(data.columns)
        index_num = len(data1.index)
        data2 = {
            "pan_id": new_inserted_id,
                }
        for i in range(index_num):
            index_name = str(sample_list[i]).strip("g")
            new_list = list(data1.iloc[i])
            # num_list = map(lambda x:str(i), new_list)
            num_list =  [str(x) for x in new_list]
            num_string = ",".join(num_list)
            # data2[index_name] = float(np.array(num_list).mean())
            if type in ["core", "pan"]:
                data2[index_name] = num_string
            elif type in ["newgene"]:
                num_list =  [float(x) for x in num_list]
                data2[index_name] = float(np.array(num_list).mean())
        coresize = "-"
        pansize = "-"
        if type in ["pan"]:
            pan_number = list(data1.iloc[index_num - 1])
            pan_number = [float(x) for x in pan_number]
            pansize = float(pan_number[0])

        elif type in ["core"]:
            core_number = list(data1.iloc[index_num - 1])
            core_number = [float(x) for x in core_number]
            coresize = float(core_number[0])

        data2["type"] = type
        data_son = SON(data2)
        data_list.append(data_son)

        with open(genomes, 'r') as f: ###画图的公式和划图线的一条数据
            data = {}
            lines = f.readlines()
            if type in ["newgene"]:
                if len(lines) >= 4:
                    formula = lines[2].strip()
                    formula = re.sub("\"", "", formula)
                    lin = lines[5].strip().split("\t")
                    # self.bind_object.logger.info("cluster_data:path %s" % (lin))
                    for i in range(len(sample_list)):
                        sample = sample_list[i].replace('g', '')
                        data[sample] = str(float(lin[i]))
                    data["pan_id"] = new_inserted_id
                    data['mean'] = 'mean'
                    data["type"] = type
                    data_son = SON(data)
                    data_list.append(data_son)
            elif type in ["pan", "core"]:
                if len(lines) >= 4:
                    if not lines[2].startswith("#"):
                        formula = lines[2].strip()
                        formula = re.sub("\"", "", formula)
                    else:
                        formula = lines[3].strip()
                        formula = re.sub("\"", "", formula)
                    for li in lines[4:]:
                        if li.startswith("#"):
                            continue
                        new_lin = li.strip().split("\t")
                        geno_name = str(new_lin[0])
                        gene_num = new_lin[1]
                        data[geno_name] = str(float(gene_num))
                    data["pan_id"] = new_inserted_id
                    data['mean'] = 'mean'
                    data["type"] = type
                    data_son = SON(data)
                    data_list.append(data_son)
            try:
                collection = self.db["pan_graph"]
                collection.insert_many(data_list)
                main_collection = self.db["pan"]
                if len(lines) >= 4:
                    if type in ["pan"]:
                        self.bind_object.logger.info("开始更新主表的公式")
                        main_collection.update({"_id": new_inserted_id}, {'$set': {"pan_genomes": formula,
                                                                                       "pansize": int(pansize),
                                                                                   'main_id':new_inserted_id}})
                    elif type in ["core"]:
                        self.bind_object.logger.info("开始更新主表的公式")
                        main_collection.update({"_id": new_inserted_id}, {'$set': {"core_genomes": formula,
                                                                                       "coresize": int(coresize),
                                                                                   'main_id':new_inserted_id}})
                    elif type in ["newgene"]:
                        self.bind_object.logger.info("开始更新主表的公式")
                        main_collection.update({"_id": new_inserted_id}, {'$set': {"newgene": formula,
                                                                                   'main_id':new_inserted_id}})
            except Exception, e:
                self.bind_object.logger.info("导入pangenomes结果表出错:%s" % (e))
            else:
                self.bind_object.logger.info("导入pangenomes结果表成功!")

    def add_pan_variation(self, inserted_id, variation):
        """
        导入变异分析结果表
        :param inserted_id: 主表ID
        :param variation: 导入变异分析的表格
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        with open(variation, 'r') as f:
            f.readline()
            new_insert = 0
            for log_index,line in enumerate(f):
                if log_index % 200000 == 0 and log_index > 0:
                    new_insert = 1
                    self.bind_object.logger.info("have done %s lines" % log_index)
                line = line.strip().split('\t')
                data = {
                    "cluster_id": line[0],
                    "spe_num": float(line[1]),
                    "gene_num": float(line[2]),
                    "pan_id": new_inserted_id,
                    "pos": line[3],
                    "aa_type": line[4],
                    "nt_type": line[5],
                    "nt_seq": line[6],
                    "variation_type": line[7]
                    }
                data_list.append(data)
                if new_insert == 1 and log_index > 0:
                    self.bind_object.logger.info("开始导入数据库前%s lines" % log_index)
                    try:
                        collection = self.db['pan_variation']
                        collection.insert_many(data_list,ordered=True)
                        new_insert = 0
                        data_list = []
                    except Exception,e:
                        self.bind_object.logger.error("导入前%slines失败:%s" % (log_index,e))
                        self.bind_object.set_error("导入数据失败")
            self.bind_object.logger.info("导入剩余数据")
            try:
                collection = self.db["pan_variation"]
                collection.insert_many(data_list)
                main_collection = self.db["pan"]
                main_collection.update({'_id': new_inserted_id},{'$set':{'main_id': new_inserted_id}})
            except Exception, e:
                self.bind_object.logger.info("导入variation结果表出错:%s" % (e))
            else:
                self.bind_object.logger.info("导入variation结果表成功!")

    def add_pan_formula(self, inserted_id, genome_data):
        """
        添加主表的公式字段
        :param inserted_id: 主表id
        :param genome_data: pgap计算的结果
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        formula = ''
        core_formula = ""
        if genome_data:
            with open(genome_data, "r") as f:
                lines = f.readlines()
                for line in lines[6:10]:
                    line = line.strip()
                    if re.search(r"y =", line):
                        formula = line
                for line in lines[20:]:
                    line = line.strip()
                    if re.search(r"y =", line):
                        core_formula = line.strip().split("R-square")[0].strip()
            try:
                collection = self.db["pan"]
                collection.update({'_id': new_inserted_id},{'$set':{'pgap_formula': formula,"pgap_core": core_formula}})
                self.bind_object.logger.info("更新pan主表成功!")
            except Exception, e:
                self.bind_object.logger.info("更新pan主表出错:%s" % (e))

    def add_pan_variation_stat(self, inserted_id, variation):
        """
        导入变异分析统计表
        :param inserted_id: 主表ID
        :param variation: 导入变异分析统计的表格
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        with open(variation, 'r') as f:
            f.readline()
            new_insert = 0
            for log_index,line in enumerate(f):
                if log_index % 200000 == 0 and log_index > 0:
                    new_insert = 1
                    self.bind_object.logger.info("have done %s lines" % log_index)
                line = line.strip().split('\t')
                data = {
                    "cluster_id": line[0],
                    "pan_id": new_inserted_id,
                    "indel_num": int(line[1]),
                    "nonsy_num": int(line[2]),
                    "sy_num": int(line[3]),
                    }
                data_list.append(data)
                if new_insert == 1 and log_index > 0:
                    self.bind_object.logger.info("开始导入数据库前%s lines" % log_index)
                    try:
                        collection = self.db['pan_variation']
                        collection.insert_many(data_list,ordered=True)
                        new_insert = 0
                        data_list = []
                    except Exception,e:
                        self.bind_object.logger.error("导入前%slines失败:%s" % (log_index,e))
                        self.bind_object.set_error("导入数据失败")
            self.bind_object.logger.info("导入剩余数据")
            try:
                collection = self.db["pan_variation"]
                collection.insert_many(data_list)
                main_collection = self.db["pan"]
                main_collection.update({'_id': new_inserted_id},{'$set':{'main_id': new_inserted_id}})
            except Exception, e:
                self.bind_object.logger.info("导入variation结果表出错:%s" % (e))
            else:
                self.bind_object.logger.info("导入variation结果表成功!")

    def find_specimen(self, pan_id):
        """
        根据pan_id查到对应的样本名称
        :param pan_id:
        :return:
        """
        collection = self.db["pan"]
        result = collection.find_one({"main_id": ObjectId(pan_id)})
        if result:
            specimen_list = str(result['spe_list']).split(",")
            return specimen_list
        else:
            self.bind_object.set_error("未能查到对应主表的样本")

    def add_update_pan_formula(self, pan_id):
        """
        更新主表
        :param pan_id:
        :return:
        """
        try:
            new_inserted_id = ObjectId(pan_id)
            collection = self.db["pan_formula"]
            collection.update({'_id': new_inserted_id},{'$set':{'main_id': new_inserted_id}})
            self.bind_object.logger.info("更新pan_formula主表成功!")
        except Exception, e:
            self.bind_object.logger.info("更新pan_formula主表出错:%s" % (e))

