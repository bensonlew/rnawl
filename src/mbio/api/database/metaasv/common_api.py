# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import datetime,json
from bson.objectid import ObjectId
import pandas as pd
import re
from Bio import SeqIO
from bson.son import SON
from types import DictType, StringTypes
from mainapp.libs.param_pack import group_detail_sort, param_pack


class CommonApi(Base):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(CommonApi, self).__init__(bind_object)
        self._project_type = "metaasv"

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
                        main_table=None, update_dic=None, convert_dic=None):
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
                spline = line.strip("\n\t").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] != '':
                            insert_data[key_list[i]] = self.convert_format(spline[i], key_list[i], convert_dic)
                else:
                    self.bind_object.set_error("data incomplete: %s ", variables=(line), code="52803302")

                data_list.append(insert_data)
        try:
            collection = self.db[detail_table]
            self.bind_object.logger.info("开始导入%s数据")
            collection.insert_many(data_list)
            if update_dic:
                main_table = self.db[main_table]
                main_table.update({'_id':main_id},{'$set':update_dic})
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

    def update_genome_taxon(self, genome_id, detail_table, task_id=None):
        if not task_id:
            task_id = "_".join(self.bind_object.sheet.id.split("_")[:2])
        result = self.db["genome"].find_one({"task_id": task_id, "genome_id": genome_id})
        data = pd.read_table(detail_table)
        data.sort_values("Identity(%)", ascending=False, inplace=True)
        taxon = data["Taxonomy"][0]
        self.db["genome"].update({"_id": result["_id"]}, {"$set": {"taxon": taxon}})

    def add_sg_status(self, db_name, desc=None, submit_location=None, params=None, table_id=None, table_name=None, genome_id=None, type_name=None):
        task_id ='_'.join(self.bind_object.sheet.id.split("_")[0:2])
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": desc if desc else "Job has been finished",
            "submit_location": submit_location,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "table_id": ObjectId(table_id),
            "table_name": table_name,
            "genome_id": genome_id,
            "type_name": type_name
        }
        collection = self.db[db_name]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def update_status_genome_id(self, table_id, genome_id):
        if not isinstance(table_id, ObjectId):
            table_id = ObjectId(table_id)
        collection = self.db["sg_status"]
        collection.update({"table_id": table_id}, {"$set": {"genome_id": genome_id}})

    def insert_ellipse_table(self, infile, main_id, analysis_type):
        ### 导入置信椭圆数据
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        insert_data = []
        name_list = []
        group_list = []
        with open(infile) as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                name = line[0]
                group = line[1]
                if name not in name_list:
                    name_list.append(name)
                if group not in group_list:
                    group_list.append(group)

            self.bind_object.logger.info("name_list : {}".format(name_list))
            for name2 in name_list:
                tmp = {}
                self.bind_object.logger.info("name2 : {}".format(name2))
                for line in lines[1:]:
                    line = line.strip()
                    spline = line.split('\t')
                    new_name = spline[0].strip()
                    if new_name in [name2]:
                        x = new_name.split("_")[0]
                        y = new_name.split("_")[1]
                        if analysis_type in ['pca', 'pcoa', 'nmds', 'dbrda']:
                            multi_analysis_id = str(analysis_type) + "_id"
                        elif analysis_type in ['rda_cca']:
                            multi_analysis_id = "rda_id"
                        else:
                            multi_analysis_id = "multi_analysis_id"
                        if tmp == {}:
                            tmp = {'name': new_name,
                                    'method': 'ci_circle',
                                    'type': 'ellipse',
                                    spline[1]: ','.join(spline[2:]),
                                    '{}'.format(multi_analysis_id): main_id,
                                    "x" : x,
                                    "y" : y}
                        else:
                            tmp = tmp
                            if new_name in [name2]:
                                tmp.update({
                                    spline[1]: ','.join(spline[2:])
                                })

                insert_data.append(tmp)
        try:
            if analysis_type in ['rda_cca', 'dbrda']:
                collection_name = str(analysis_type) + "_plugin"
                collection = self.db[collection_name]
                collection.insert_many(insert_data)
            elif analysis_type in ['pca', 'pcoa', 'nmds']:
                collection_name = str(analysis_type) + "_scatter_plugin"
                collection = self.db[collection_name]
                collection.insert_many(insert_data)
            main_collection = self.db[analysis_type]
            find_results = main_collection.find_one({"_id": main_id})
            if find_results:
                condition_dict = {}
                ellipse_data = json.loads(find_results["ellipse_data"])
                condition_list = ellipse_data["ellipse_data"]["condition"]['type']
                if len(condition_list) == 1:
                    condition_list.append("ci_circle")
                    condition_dict['type'] = condition_list
                    ellipse_data["condition"] = condition_dict
                    ellipse_data_json = json.dumps(ellipse_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": main_id}, {"$set": {"ellipse_data": ellipse_data_json,
                                                                       "main_id": main_id}})
            else:
                ellipse_data = {
                    "ellipse_data": {"data": group_list,
                    "data_option" : group_list,
                    "condition": {"type": ["ci_circle"]}}
                }
                ellipse_data_json = json.dumps(ellipse_data, sort_keys=True, separators=(',', ':'))
                main_collection.update_one({"_id": main_id}, {"$set": {"ellipse_data": ellipse_data_json,
                                                                       "main_id": main_id}})
        except Exception as e:
            self.bind_object.logger.error("导入detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入detail表格成功")

    def insert_anosim_detail(self, file_path, main_id, diff_type, correlation_key="pca_id",coll_name="pca_table", main_coll='pca'):
        ## 添加pca/pcoa/nmds组间差异检验结果
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        insert_data = {}
        with open(file_path, 'r') as f:
            all_lines = f.readlines()
            for line in all_lines[1:]:
                values = re.split('\t', line.strip())
                if diff_type == "anosim":
                    insert_data = {
                        correlation_key : main_id,
                        "method": values[0],
                        "statistic": values[1],
                        "p_value": values[2],
                        "permutation": values[3],
                        "type": diff_type
                    }
                elif diff_type == "adonis":
                    insert_data = {
                        correlation_key : main_id,
                        "method": values[0],
                        "df": values[1],
                        "sums_of_sqs": values[2],
                        "meansqs": values[3],
                        "f_model": values[4],
                        "r2": values[5],
                        "p_value": values[6],
                        "type": diff_type
                    }
                if insert_data["method"] in ["anosim"]:
                    insert_data["method"] = "ANOSIM"
                elif insert_data["method"] in ["adonis"]:
                    insert_data["method"] = "Adonis"
                data_list.append(insert_data)
            try:
                self.bind_object.logger.info("main_id：{}".format(main_id))
                collection = self.db[coll_name]
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("组间差异检验结果导入{}表格信息出错:{}".format(coll_name,e))
            else:
                self.bind_object.logger.info("组间差异检验结果导入{}表格成功".format(coll_name))
            try:
                main_collection = self.db[main_coll]
                if diff_type == "anosim":
                    adnosim_table_data = {
                            "table_data": ["method", "statistic", "p_value", "permutation"],
                            "condition": {"type": "anosim"}
                        }
                    text_data = {
                            "text_data": ["statistic", "p_value"],
                            "condition": {"type": "anosim"}
                        }
                    adnosim_table_data_json = json.dumps(adnosim_table_data, sort_keys=True, separators=(',', ':'))
                    text_data_json = json.dumps(text_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": main_id}, {"$set": {"anosim_table_data": adnosim_table_data_json,
                                                                           "text_data": text_data_json}})
                elif diff_type == "adonis":
                    adonis_table_data = {
                            "table_data": ["method", "df", "sums_of_sqs","f_model", "r2", "p_value"],
                            "condition": {"type": "adnois"}
                        }
                    text_data = {
                            "text_data": ["r2","p_value"],
                            "condition": {"type": "adnois"}
                        }
                    adonis_table_data_json = json.dumps(adonis_table_data, sort_keys=True, separators=(',', ':'))
                    text_data_json = json.dumps(text_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update_one({"_id": main_id}, {"$set": {"adonis_table_data": adonis_table_data_json,
                                                                           "text_data": text_data_json}})
                self.bind_object.logger.info('导入%s成功'%coll_name)
            except Exception as e:
                self.bind_object.logger.error("组间差异检验结果导入{}表格信息出错:{}".format(coll_name,e))
            else:
                self.bind_object.logger.info("组间差异检验结果导入{}表格成功".format(coll_name))

    def add_otu_table(self, file_path, major=False, otu_id=None, from_out_table=0, rep_path=None, task_id=None,
                      name=None, params=None, spname_spid=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        if from_out_table != 0 and not isinstance(from_out_table, ObjectId):
            if isinstance(from_out_table, StringTypes):
                from_out_table = ObjectId(from_out_table)
            else:
                self.bind_object.set_error("from_out_table必须为ObjectId对象或其对应的字符串!")
        if major:
            if spname_spid and params:
                group_detail = {'All': [str(i) for i in spname_spid.values()]}
                params['group_detail'] = group_detail_sort(group_detail)
            if task_id is None:
                task_id = self.bind_object.sheet.id
            insert_data = {
                "project_sn": self.bind_object.sheet.project_sn,
                "task_id": task_id,
                "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "ASV_Taxon_Origin",
                "status": "end",
                "level_id": 9,
                "set_id": "",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            collection = self.db["asv"]
            otu_id = collection.insert_one(insert_data).inserted_id
            params['asv_id'] = str(otu_id)
            insert_data["asv_id"] = otu_id
            new_params = param_pack(params)
            insert_data["params"] = new_params
            try:
                collection.update({"_id": otu_id}, {'$set': insert_data})
            except Exception as e:
                self.bind_object.logger.error('更新ASV表params出错:{}'.format(e))
                self.bind_object.set_error("ASV表更新出错")
        else:
            if otu_id is None:
                self.bind_object.set_error("major为False时需提供asv_id!")
        data_list = []
        # 读代表序列
        otu_reps = {}
        if rep_path:
            for record in SeqIO.parse(rep_path, 'fasta'):
                seq_id = record.id
                sequence = str(record.seq)
                if seq_id not in otu_reps.keys():
                    otu_reps[seq_id] = sequence
        # self.bind_object.logger.info("otu_reps : {}".format(otu_reps.keys()))
        with open(file_path, 'r') as f:
            l = f.readline()
            sample_list = l.strip().split("\t")[1:-1]
            sample_data = []
            for sample in sample_list:
                sample_data.append({"asv_id": otu_id, "specimen_id": spname_spid[sample]})
            collection = self.db["asv_specimen"]
            collection.insert_many(sample_data)
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                classify = line_data.pop()
                classify_list = re.split(r"\s*;\s*", classify)
                otu_list = [("task_id", task_id), ("asv_id", otu_id)]
                if rep_path:
                    otu_list.append(("asv_rep", otu_reps[line_data[0]]))
                for cf in classify_list:
                    if cf != "":
                        otu_list.append((cf[0:3], cf))
                i = 0
                otu_list.append(("asv", line_data[0]))
                for sample in sample_list:
                    i += 1
                    otu_list.append((sample, int(line_data[i])))
                data = SON(otu_list)
                data_list.append(data)
        try:
            collection = self.db["asv_detail"]
            collection.insert_many(data_list)
            collection_main = self.db["asv"]
            #collection_main.update({"_id": otu_id}, {"$set": {"main_id" : otu_id}})

        except Exception, e:
            self.bind_object.logger.error("导入ASV表格%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入ASV表格%s信息成功!" % file_path)
        return otu_id

    def add_meta_status(self, table_id=None, table_name=None, type_name=None, status='end', task_id=None, isnew='new', desc=None, submit_location=None, params=None):
        collection = self.db["sg_status"]
        if not task_id:
            task_id = self.bind_object.sheet.id
        else:
            if params is not None:
                if not isinstance(params, dict):
                    self.bind_object.set_error('提供的参数params必须为字典')
        if isinstance(table_id, StringTypes):
            try:
                table_id = ObjectId(table_id)
            except Exception as e:
                self.bind_object.logger.error('参数table_id不是mongo表ID类型:{}'.format(e))
                self.bind_object.set_error("table_id不是ObjectId类型")
        elif isinstance(table_id, ObjectId):
            pass
        else:
            self.bind_object.set_error('提供的table_id不是字符串或者ObjectId类型')
        if not type_name:
            self.bind_object.set_error('必须提供type_name!!!')
        if not table_name or not params:
            temp_collection = self.db[type_name]
            tempfind = temp_collection.find_one({'_id': table_id})
            if not tempfind:
                self.bind_object.logger.error('提供的ID:%s无法在表:%s中找到' % (table_id, type_name))
                self.bind_object.set_error("找不到相应记录")
            if 'name' in tempfind:
                table_name = tempfind['name']
            else:
                self.bind_object.logger.error('表：%s中没有name字段' % type_name)
                self.bind_object.set_error("表中没有name字段")
            if 'params' in tempfind:
                params = tempfind['params']
        if not submit_location:
            if isinstance(params, StringTypes):
                if 'submit_location' in params:
                    temp_params = json.loads(params)
                    submit_location = temp_params['submit_location']
            elif isinstance(params, DictType):
                submit_location = params['submit_location'] if 'submit_location' in params else None
        insert_data = {
            "table_id": table_id,
            "table_name": table_name,
            "task_id": task_id,
            "type_name": type_name,
            "status": status,
            "is_new": isnew,
            "desc": desc,
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, DictType) else params,
            "submit_location": submit_location,
            "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        insert_data = SON(insert_data)
        return_id = collection.insert_one(insert_data).inserted_id
        if return_id:
            return return_id
        else:
            self.bind_object.set_error('插入sg_status表出错')

    @report_check
    def add_database(self, task_id, database):
        try:
            collection = self.db["sg_task"]
            tempfind = collection.find_one({'task_id': task_id})
            collection.update_one({'_id': tempfind["_id"]}, {'$set': {'database': database}})
        except Exception as e:
            self.bind_object.logger.error("database导入sg_task表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("database导入sg_task表格成功")

    @report_check
    def find_group_name(self, task_id):
        try:
            collection = self.db["specimen_group"]
            tempfind = collection.find_one({'task_id': task_id})
            group_detail = {}
            if tempfind:
                group_names = tempfind['category_names']
                specimen_names = tempfind["specimen_names"]
                n_group_names = range(len(group_names))
                for i in n_group_names:
                    group = group_names[i]
                    if group not in group_detail:
                        group_detail[group] = ",".join([str(x) for x in specimen_names[i].keys()])
                # main_id = tempfind['_id']
                return group_detail
        except Exception as e:
            self.bind_object.logger.error("查找group名称出错:{}".format(e))
        else:
            self.bind_object.logger.info("查找group名称成功！")
