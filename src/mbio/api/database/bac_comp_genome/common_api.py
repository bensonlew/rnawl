# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# 20190916

from biocluster.api.database.base import Base, report_check
import datetime,json
from bson.objectid import ObjectId
import pandas as pd
from biocluster.config import Config
from Bio import SeqIO
import os,re
import copy
from bson.son import SON

class CommonApi(Base):
    def __init__(self, bind_object):
        super(CommonApi, self).__init__(bind_object)
        self._project_type = "bac_comparative"

    @report_check
    def add_main(self,table, name=None, params=None, others=None,desc=None, task_id=None):
        if task_id == None:
            task_id = self.bind_object.sheet.id
        else:
            task_id = task_id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": desc if desc else "Job has been finished",
            "name": name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        if others:
            for k in others.keys():      #others {mongo字段名: mongo字段值}
                insert_data[k] = others[k]
        collection = self.db[table]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id}, {'$set': {"main_id": ObjectId(inserted_id)}})
        return inserted_id

    @report_check
    def add_main_detail(self, infile, detail_table, main_id, mongo_key, has_head =True, main_name='main_id',
                        main_table=None, update_dic=None, convert_dic=None, insert_d = None):
        """
        :param infile:  需要导表的数据文件
        :param detail_table: detail表名称
        :param main_id: detail表对应的主表_id
        :param mongo_key: detail表中的字段名称,逗号分割
        :param has_head: 数据文件是否有表头，有则删除表头
        :param main_name: detail表中对应主表id的字段名称
        :param main_table: detail表对应的主表名称
        :param update_dic: 需要更新主表中的内容
        :param convert_dic: 更新detail表时数据格式的转换方法，需自定义函数{“字段1”： 转换函数1，“字段2”： 转换函数2}
                            不需要转换则不用写进参数变量中
        :return:
        """
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        key_list = mongo_key.split(',')
        num = len(key_list)
        with open(infile, 'rb') as r:
            if has_head:
                r.readline()
            for line in r:
                insert_data = {main_name: main_id}
                if insert_d:
                    for d in insert_d.keys():
                        insert_data[d] = insert_d[d]
                spline = line.strip("\n\t").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] != '':
                            insert_data[key_list[i]] = self.convert_format(spline[i], key_list[i], convert_dic)
                else:
                    self.bind_object.set_error("data incomplete: %s ", variables=(line))
                data_list.append(insert_data)
        try:
            collection = self.db[detail_table]
            self.bind_object.logger.info("开始导入%s数据")
            collection.insert_many(data_list)
            if update_dic:
                main_table = self.db[main_table]
                main_table.update({'_id': main_id}, {'$set': update_dic})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR")
        else:
            self.bind_object.logger.info("%s导入CommonApi数据成功" % infile)

    def convert_format(self, value, key_name, convert_dic):
        if convert_dic and convert_dic.has_key(key_name):
            try:
                value = convert_dic[key_name](value)
                return value
            except:
                self.bind_object.logger.error("字段%s,数据值%s,用方法%s转换不成功" % (key_name, value, convert_dic[key_name].__name__))
                return value
        else:
            return value

    @report_check
    def add_main_detail2(self, infile, detail_table, main_id, mongo_key, has_head=True, main_name='main_id',
                        main_table=None, convert_dic=None, insert_d =None, names_convert=None, names_co='True'):
        """
        :param infile:  需要导表的数据文件
        :param detail_table: detail表名称
        :param main_id: detail表对应的主表_id
        :param mongo_key: detail表中的字段名称,逗号分割
        :param has_head: 数据文件是否有表头，有则删除表头
        :param main_name: detail表中对应主表id的字段名称
        :param main_table: detail表对应的主表名称
        :return:
        """
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        key_list = mongo_key.split(',')
        num = len(key_list)
        with open(infile, 'rb') as r:
            lines =[]
            if has_head:
                lines = r.readlines()
            header = lines[0].strip("\n\r").split("\t")
            names = header[num:]
            for line in lines[1:]:
                insert_data = {main_name: main_id}
                if insert_d:
                    for d in insert_d.keys():
                        insert_data[d] = insert_d[d]
                spline = line.strip("\n").split("\t")
                if len(spline) == num + len(names):
                    for i in range(num):
                        if key_list[i] != '':
                            insert_data[key_list[i]] = self.convert_format(spline[i], key_list[i], convert_dic)
                else:
                    self.bind_object.set_error("data incomplete: %s ", variables=(line))
                if names_convert:
                    for name in names:
                        if names_co in ["True","true"]:
                            insert_data[names_convert[name]] = float(spline[header.index(name)])
                        elif names_co in ["False","false"]:
                            insert_data[names_convert[name]] = spline[header.index(name)]
                data_list.append(insert_data)
        try:
            collection = self.db[detail_table]
            self.bind_object.logger.info("开始导入%s数据")
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR")
        else:
            self.bind_object.logger.info("%s导入CommonApi数据成功" % infile)

    @report_check
    def add_main_detail3(self, infile, infile2, detail_table, main_id, mongo_key, type, has_head=True, main_name='main_id',
                        main_table=None, update_dic=None, convert_dic=None, insert_d= None):
        """
        :param infile:  需要导表的数据文件
        :param detail_table: detail表名称
        :param main_id: detail表对应的主表_id
        :param mongo_key: detail表中的字段名称,逗号分割
        :param has_head: 数据文件是否有表头，有则删除表头
        :param main_name: detail表中对应主表id的字段名称
        :param main_table: detail表对应的主表名称
        :param update_dic: 需要更新主表中的内容
        :param convert_dic: 更新detail表时数据格式的转换方法，需自定义函数{“字段1”： 转换函数1，“字段2”： 转换函数2}
                            不需要转换则不用写进参数变量中
        :return:
        """
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        key_list = mongo_key.split(',')
        num = len(key_list)
        with open(infile, 'rb') as r:
            if has_head:
                r.readline()
            for line in r:
                insert_data = {main_name:main_id}
                if insert_d:
                    for d in insert_d.keys():
                        insert_data[d] = insert_d[d]
                spline = line.strip("\n\t").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] != '':
                            insert_data[key_list[i]] = self.convert_format(spline[i], key_list[i], convert_dic)
                else:
                    self.bind_object.set_error("data incomplete: %s ", variables=(line))
                data_list.append(insert_data)
        try:
            collection = self.db[detail_table]
            self.bind_object.logger.info("开始导入%s数据")
            collection.insert_many(data_list)
            if update_dic:
                main_table = self.db[main_table]
                main_table.update({'_id':main_id},{'$set':update_dic})
            uniprot_iterator = SeqIO.parse(infile2, "fasta")
            records = list(uniprot_iterator)
            for i in records:
                id, sample = i.id.split("__")
                detail_table2 = self.db[detail_table]
                detail_table2.update({main_name: main_id, type: id, "specimen_id": sample}, {'$set': {"seq": "{}".format(i.seq)}})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR")
        else:
            self.bind_object.logger.info("%s导入CommonApi数据成功" % infile)

    @report_check
    def add_main_tree(self, infile, main_id, infile2=None, update_dic=None, insert_d=None):
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        with open(infile, 'rb') as r:
            line = r.readline()
            insert_data = {"tree_id": main_id, "value": line}
            if insert_d:
                for d in insert_d.keys():
                    insert_data[d] = insert_d[d]
            data_list.append(insert_data)
        try:
            collection = self.db["tree_detail"]
            self.bind_object.logger.info("开始导入%s数据")
            collection.insert_many(data_list)
            if update_dic:
                main_table = self.db["tree"]
                main_table.update({'_id': main_id}, {'$set': update_dic})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR")
        else:
            self.bind_object.logger.info("%s导入CommonApi数据成功" % infile)
        if infile2:
            data_list2 = []
            with open(infile2, 'rb') as r:
                lines = r.readlines()
                for line in lines[1:]:
                    lin = line.strip().split("\t")
                    insert_data2 = {
                        "tree_id": main_id,
                        "model": lin[0],
                        "aic": lin[2],
                        "aicc": lin[4],
                        "bic": lin[6]
                    }
                    data_list2.append(insert_data2)
            try:
                collection = self.db["modul_detail"]
                self.bind_object.logger.info("开始导入%s数据" % infile2)
                collection.insert_many(data_list2)
            except Exception, e:
                self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
                self.bind_object.set_error("import CommonApi data ERROR")
            else:
                self.bind_object.logger.info("%s导入CommonApi数据成功" % infile2)

    @report_check
    def add_main_tree2(self, infile, table_name, main_id, update_dic=None):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        with open(infile, 'rb') as r:
            line = r.readline()
            main_table = self.db[table_name]
            main_table.update({'_id': main_id}, {'$set': {"tree": line,"main_id": main_id}})
        try:
            if update_dic:
                main_table = self.db[table_name]
                main_table.update({'_id': main_id}, {'$set': update_dic})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR")
        else:
            self.bind_object.logger.info("%s导入CommonApi数据成功" % infile)


    @report_check
    def add_data_specimen(self,table_name, mongo_key, main_id ,main_name):
        """
        table_name:指导表的名称
        mongo_key：dict的原始样品：修改后的样品名称
        main_id：指主表id
        main_name：主表到detail的字段名称
        :return:
        """
        data_list = []
        for i in mongo_key.keys():
            insert_data = {
                 main_name: main_id,
                "specimen_id": i,
                "specimen_m": mongo_key[i],
                "new_spe_name": i,
            }
            data_list.append(insert_data)
        try:
            collection = self.db[table_name]
            self.bind_object.logger.info("开始导入%s数据" % table_name)
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入add_data_specimen数据出错:%s" % e)
            self.bind_object.set_error("importadd_data_specimen data ERROR")
        else:
            self.bind_object.logger.info("%s导入add_data_specimen数据成功" % table_name)

    @report_check
    def add_data_gene(self, database, samples, table_name, mongo_key, main_id, main_name, insert_da=None, convert_dic=None):
        """
        从majorbiodb查出数据导入新数据库
        samples:查询的样品名称
        table_name:指需要导detail表的名称
        mongo_key：指是字典形式{"存数据库的key":"取数据库的key"}
        main_id：指主表id
        main_name：主表到detail的字段名称
        :return:
        """
        #mj_db = Config().biodb_mongo_client.sanger_biodb
        dbclient = Config().get_mongo_client(mtype="metagenomic", ref=True)
        mj_db = dbclient[Config().get_mongo_dbname("metagenomic", ref=True)]
        db = mj_db[database]
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        for sample in samples:
            result = db.find_one({"ass_assession": sample})
            if result:
                insert_data = {"specimen_id": sample}
                insert_data[main_name] = main_id
                if insert_da:
                    for j in insert_da.keys():
                        insert_data[j] = insert_da[j]
                if convert_dic:
                    for i in mongo_key.keys():
                        insert_data[i] = self.convert_format(result[mongo_key[i]], i, convert_dic)
                else:
                    for i in mongo_key.keys():
                        insert_data[i] =result[mongo_key[i]]
                data_list.append(insert_data)
            else:
                insert_data = {"specimen_id": sample}
                insert_data[main_name] = main_id
                insert_data["collection_date"] = "-"
                insert_data["iso_source"] = "-"
                insert_data["env_biome"] = "-"
                insert_data["sample_type"] = "-"
                insert_data["from_refdb"] = "T"
                insert_data["geo_loc_name"] = "-"
                insert_data["host"] = "-"
                insert_data["lat_lon"] = "-"
                data_list.append(insert_data)
                self.bind_object.logger.info("MajorbioDB数据库不存在{}的数据库".format(sample))
        try:
            collection = self.db[table_name]
            self.bind_object.logger.info("开始导入%s数据" % table_name)
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入add_data_gene数据出错:%s" % e)
            self.bind_object.set_error("import add_data_gene data ERROR")
        else:
            self.bind_object.logger.info("%s导入add_data_gene数据成功" % table_name)

    @report_check
    def add_data_gene2(self, file, table_name, mongo_key, main_id, main_name, has_head=True, insert_da=None, dict_da=None, up_date=None):
        """
        samples:查询的样品名称
        table_name:指需要导detail表的名称
        mongo_key：指是字典形式{"存数据库的key":"取数据库的key"}
        main_id：指主表id
        main_name：主表到detail的字段名称
        name：
        has_head :是否有表头
        :return:
        """
        data_list = []
        key_list = mongo_key.split(',')
        num = len(key_list)
        with open(file, 'rb') as r:
            if has_head:
                r.readline()
            for line in r:
                insert_data = {main_name: main_id}
                spline = line.strip("\n\t").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] != '':
                            if insert_da:
                                if key_list[i] in insert_da.keys():
                                    for j in insert_da[key_list[i]]:
                                        insert_data[j] = spline[i]
                                else:
                                    insert_data[key_list[i]] = spline[i]
                            else:
                                insert_data[key_list[i]] = spline[i]
                else:
                    self.bind_object.set_error("data incomplete: %s ", variables=(line))
                if up_date:
                    if spline[0] in up_date.keys():
                        insert_data['strain'] = up_date[spline[0]]
                if dict_da:
                    for j in dict_da.keys():
                        insert_data[j] = dict_da[j]
                data_list.append(insert_data)
        try:
            collection = self.db[table_name]
            self.bind_object.logger.info("开始导入%s数据" % table_name)
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入add_data_gene数据出错:%s" % e)
            self.bind_object.set_error("import add_data_gene data ERROR")
        else:
            self.bind_object.logger.info("%s导入add_data_gene数据成功" % table_name)

    @report_check
    def add_num_update(self, dir, type_spli, table_name, main_id, main_name, type, sample_list):
        """
        dir :导入表的文件目录
        type_spli:  分割文件字段
        table_name：导入表的mongo的detail表名
        main_id：主表名称的id
        main_name：主表字段名称
        type：导入表中的字段名
        :return:
        """
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        files = os.listdir(dir)
        list1 = copy.copy(sample_list)
        if len(files) >0:
            for file in files:
                if file.endswith(type_spli):
                    sample = file.split(type_spli)[0]
                    self.bind_object.logger.info("aaa:%s" % sample)
                    list1.remove(sample)
                    with open(dir + "/" + file, "r") as f:
                        lines = f.readlines()
                    num = len(lines) - 1
                    collection = self.db[table_name]
                    self.bind_object.logger.info("导入add_num_update数据:%s" % sample)
                    collection.update({main_name: main_id, "specimen_id": sample}, {'$set': {type: int(num)}})
        if len(list1) >0:
            for sample in list1:
                collection = self.db[table_name]
                self.bind_object.logger.info("导入add_num_update数据不存在:%s" % sample)
                collection.update({main_name: main_id, "specimen_id": sample}, {'$set': {type: 0}})

    @report_check
    def add_num_anno_update(self, file, gene_id, sample_name, table_name, main_id, main_name, sample_list, type):
        """
        file :导入表的文件是所有样品的总文件，考虑部分样品没有结果
        gene_id:  根据gene_id字段统计个数
        sample_name :不同文件的样品名称的列名
        table_name：导入表的mongo的detail表名
        main_id：主表名称的id
        main_name：主表字段名称
        sample_list:所有样品的
        type：导入表中的字段名
        :return:
        """
        a = pd.read_table(file, sep='\t', header=0, dtype={'fullVisitorId': 'str'})
        a = dict(a.groupby([sample_name])[gene_id].count())
        collection = self.db[table_name]
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        for sample in sample_list:
            if sample in a:
                self.bind_object.logger.info("导入add_num_anno_update数据:%s" % sample)
                collection.update({main_name: main_id, "specimen_id": sample}, {'$set': {type: int(a[sample])}})
            else:
                self.bind_object.logger.info("导入add_num_anno_update数据不存在:%s" % sample)
                collection.update({main_name: main_id, "specimen_id": sample}, {'$set': {type: 0}})

    @report_check
    def add_gff_detail(self, file1, file2, main_id, table_name, path):
        """
        file :导入表的文件是所有样品的总文件，考虑部分样品没有结果
        gene_id:  根据gene_id字段统计个数
        sample_name :不同文件的样品名称的列名
        table_name：导入表的mongo的detail表名
        main_id：主表名称的id
        main_name：主表字段名称
        sample_list:所有样品的
        type：导入表中的字段名
        :return:
        """
        collection = self.db[table_name]
        collection.delete_many({"sample_id": ObjectId(main_id),"status":{"$exists": False}})
        dict1 = {}
        with open(file1, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                if (int(lin[1]) == 0) or (lin[4] in ["false","False"]):
                    dict1[lin[0]] = {"cds": lin[1], "trna": lin[2], 'rrna':lin[3], "gff_status": "false"}
                else:
                    dict1[lin[0]] = {"cds": lin[1], "trna": lin[2], 'rrna': lin[3], "gff_status": "true"}
        self.bind_object.logger.info(dict1)
        with open(file2, "r") as g:
            lines = g.readlines()
            data_list = []
            for line in lines[1:]:
                lin = line.strip().split("\t")
                cds = []
                rrna = []
                trna = []
                gff_status = []
                if not re.search("-",lin[6]):
                    if re.search(";", lin[6]):
                        for i in lin[6].split(";"):
                            self.bind_object.logger.info(i)
                            self.bind_object.logger.info(dict1)
                            self.bind_object.logger.info(dict1[i])
                            cds.append(dict1[i]["cds"])
                            rrna.append(dict1[i]["rrna"])
                            trna.append(dict1[i]["trna"])
                            gff_status.append(dict1[i]["gff_status"])
                    else:
                        cds.append(dict1[lin[6]]["cds"])
                        rrna.append(dict1[lin[6]]["rrna"])
                        trna.append(dict1[lin[6]]["trna"])
                        gff_status.append(dict1[lin[6]]["gff_status"])
                    if "false" in gff_status:
                        status = 'false'
                    else:
                        status = 'true'
                    data = {
                        "sample_id": ObjectId(main_id),
                        "sample": lin[0],
                        "species": lin[1],
                        "strain": lin[2],
                        "seq_status": lin[3],
                        "g_location": lin[4],
                        "seq_file": lin[5],
                        "gff_file": lin[6],
                        "seq_path": lin[7],
                        "gff_path": lin[8],
                        "cds_no": ";".join(cds),
                        "trna_no": ";".join(trna),
                        "rrna_no": ";".join(rrna),
                        "gff_status": status,
                    }
                else:
                    cds.append("-")
                    rrna.append("-")
                    trna.append("-")
                    data = {
                        "sample_id": ObjectId(main_id),
                        "sample": lin[0],
                        "species": lin[1],
                        "strain": lin[2],
                        "seq_status": lin[3],
                        "g_location": lin[4],
                        "seq_file": lin[5],
                        "gff_file": lin[6],
                        "seq_path": lin[7],
                        "gff_path": lin[8],
                        "cds_no": ";".join(cds),
                        "trna_no": ";".join(trna),
                        "rrna_no": ";".join(rrna),
                    }
                data_son = SON(data)
                data_list.append(data_son)
            try:
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file2, e))
            collection1 = self.db["sample"]
            collection1.update({"_id": ObjectId(main_id)}, {'$set': {"gff_path": path, "gff_status": 'true'}})

    @report_check
    def add_sample_detail(self, infile, detail_table, main_id, mongo_key, has_head =True, main_name='main_id',
                        main_table=None, update_dic=None, convert_dic=None, insert_d = None):
        """
        :param infile:  需要导表的数据文件
        :param detail_table: detail表名称
        :param main_id: detail表对应的主表_id
        :param mongo_key: detail表中的字段名称,逗号分割
        :param has_head: 数据文件是否有表头，有则删除表头
        :param main_name: detail表中对应主表id的字段名称
        :param main_table: detail表对应的主表名称
        :param update_dic: 需要更新主表中的内容
        :param convert_dic: 更新detail表时数据格式的转换方法，需自定义函数{“字段1”： 转换函数1，“字段2”： 转换函数2}
                            不需要转换则不用写进参数变量中
        :return:
        """
        collection = self.db[detail_table]
        collection.delete_many({"sample_id": ObjectId(main_id), "status": {"$exists": False}})
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        key_list = mongo_key.split(',')
        num = len(key_list)
        with open(infile, 'rb') as r:
            if has_head:
                r.readline()
            for line in r:
                insert_data = {main_name: main_id}
                if insert_d:
                    for d in insert_d.keys():
                        insert_data[d] = insert_d[d]
                spline = line.strip("\n\t").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] != '':
                            insert_data[key_list[i]] = self.convert_format(spline[i], key_list[i], convert_dic)
                else:
                    self.bind_object.set_error("data incomplete: %s ", variables=(line))
                data_list.append(insert_data)
        try:
            collection = self.db[detail_table]
            self.bind_object.logger.info("开始导入%s数据")
            collection.insert_many(data_list)
            if update_dic:
                main_table = self.db[main_table]
                main_table.update({'_id': main_id},{'$set': update_dic})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR")
        else:
            self.bind_object.logger.info("%s导入CommonApi数据成功" % infile)

    @report_check
    def update_data_gene(self, main_id, sample_dict):
        """
        目的是将data_gene中的gff的样本的new_gene_name更新为空，原因是让老师看到展示正确的gene_id
       :param main_id:
        :param sample_dict:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        if len(sample_dict)>0:
            for sample in sample_dict.keys():
                genome_list = sample_dict[sample]
                collection = self.db['data_gene']
                all_genome = ",".join(genome_list)
                result = collection.find_one({"data_id": main_id,"specimen_id":sample, "location": all_genome})
                zero_list = []
                for genome in genome_list:
                    zero_list.append("")
                new_gene_name = ",".join(zero_list)
                if result:
                    data_id = result['_id']
                    try:
                        collection.update_one({"_id": data_id}, {"$set": {"new_gene_name": new_gene_name}})
                    except:
                        self.bind_object.logger.info("更新data_gene表的基因前缀出错:%s" % e)
                        self.bind_object.set_error("update daa_gene data ERROR!")
                else:
                    self.bind_object.logger.info("此sample和location查不到对应信息，请检查是否已导入此样本信息:%s" % e)
                    self.bind_object.set_error("update daa_gene data ERROR!")
            self.bind_object.logger.info("更新data_gene表的基因前缀完成！")
        else:
            self.bind_object.logger.info("%s没有要修改的基因前缀" % sample_dict)

