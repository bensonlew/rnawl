# -*- coding: utf-8 -*-
# __author__ = 'xuting'

from biocluster.api.database.base import Base, report_check
import re
import os
import json
import shutil
import pandas as pd
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes
from mbio.packages.meta.otu.export_otu import export_otu_table_by_level
# from biocluster.config import Config


class Venn(Base):
    def __init__(self, bind_object):
        super(Venn, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.new_otu_id = list()
        self.single = dict()
        self.num = dict()

    @report_check
    def create_venn_table(self, params, group_id, level_id, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51007801")
        if not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!", code="51007802")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("sg_otu表里找不到相应的记录", code="51007803")
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "正在计算venn表格..."
        if name:
            table_name = name
        else:
            table_name = self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "venn表格"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "level_id": level_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "status": "end",
            "desc": desc,
            "name": table_name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_otu_venn"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_venn_detail(self, venn_path, venn_id, otu_id, level,group_path=None):
        sub_otu_dir = os.path.join(self.bind_object.work_dir, "sub_otu")
        if os.path.exists(sub_otu_dir):
            shutil.rmtree(sub_otu_dir)
        os.mkdir(sub_otu_dir)
        self._find_info(otu_id)
        self.level = level
        self.all_otu = os.path.join(sub_otu_dir, "all.otu")
        export_otu_table_by_level(otu_id, 9, self.all_otu)
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                self.bind_object.set_error("venn_id必须为ObjectId对象或其对应的字符串!", code="51007804")
        venn_json = self._get_venn_json(venn_path,group_path)
        data_list = []
        collection = self.db['sg_otu_venn_detail']
        with open(venn_path, 'rb') as r:
            count = 0
            for line in r:
                line = line.strip('\n')
                line = re.split("\t", line)
                ##qingchen.zhang@20200509注释掉以下三行数据，主要是不再去sg_otu_detail中导数据
                # new_otu_id = self._add_sg_otu(otu_id, line[0])
                # self.new_otu_id.append(new_otu_id)
                # self._add_sg_otu_detail(line[2], otu_id, new_otu_id, line[0])
                tmp_name = re.split(',', line[2])
                name_list = list()
                for cla_info in tmp_name:
                    cla_info = re.split('; ', cla_info)
                    my_name = cla_info[-1]
                    name_list.append(my_name)
                display_name = ",".join(name_list)
                # collection = self.db['sg_otu_venn_detail']
                tmp_name = re.sub("only", "", line[0])
                tmp_name = re.sub(" ", "", tmp_name)
                tmp_list = re.split("&", tmp_name)
                insert_data = {
                    'otu_venn_id': venn_id,
                    'otu_id': ObjectId(otu_id),##由new_otu_id改为ObjectId(otu_id)，qingchen.zhang@20200509
                    'category_name': line[0],
                    # 'species_name': line[2],
                    #'species_name_display': display_name,
                    'display_count': int(line[1]),
                    # 'display_count': int(self.single[tuple(tmp_list)]),
                    'venn_json': json.dumps(venn_json[count])
                }
                data_list.append(insert_data)
                # collection.insert_one(insert_data)
                count += 1
        collection.insert_many(data_list)
        try:
            main_collection=self.db['sg_otu_venn']
            main_collection.update({"main_id": venn_id}, {"$set": {"origin": "no"}})
        except:
            self.bind_object.set_error("更新主表失败！")
        # 由于需求的变更，要在原来基础上再把几个分组的all的导入sg_otu和sg_otu_species
        # 所以要先获取到all的otu，再导入, 但是这几个的venn_detail表不用再进行导入
        """
        otu_list = defaultdict(list)
        with open(venn_path, "rb") as r:
            for line in r:
                tmp_list = list()
                line = line.rstrip().split("\t")
                if re.search("only", line[0]):
                    name = re.sub(" ", "", line[0])
                    name = re.sub("only", "", name)
                    if line[1] not in [0, "0"]:  # 当值为0的时候，列表少一个元素， line[2]不存在
                        otu_list[name].extend(re.split(",", line[2]))
                if re.search("&", line[0]):
                    name = re.sub(" ", "", line[0])
                    tmp_list = re.split("&", name)
                    for g in tmp_list:
                        if line[1] not in [0, "0"]:
                            otu_list[g].extend(re.split(",", line[2]))
        for g in otu_list:
            otu_list[g] = list(set(otu_list[g]))
        for g in otu_list:
            new_otu_id = self._add_sg_otu(otu_id, g + "__all")
            self._add_sg_otu_detail(",".join(otu_list[g]), otu_id, new_otu_id, g + "__all")
        """

    def _get_venn_json(self, venn_path,group_path=None):
        num = defaultdict(int)
        single = defaultdict(int)
        only = dict()
        gp = list()
        if group_path:
            pass
        else:
            group_path = self.bind_object.option("group_table").prop['path']
        with open(group_path, "rb") as r:
            r.next()
            for line in r:
                line = re.split('\t', line)
                if line[1] not in gp:
                    gp.append(line[1])
        gp_len = len(gp)
        sum_len = defaultdict(int)
        single_len = dict()
        # sum_len{分组的数目} : 该数目之下的和  例如该方案下有四个分组a,b,c,d 那么sum_len[3] 就应该是 abc， abd, bcd， acd 数目之和。
        # singel_len[(所包含的分组)]
        # num字典中记录了取并关系的各个分组的数目
        # single字典中记录了取交欢喜的各个分组的数目, 由于多次代码更新的原因，这个可能与only字典中的部分值有冗余
        with open(venn_path, 'rb') as r:
            for line in r:
                line = line.strip('\r\n')
                line = re.split("\t", line)
                if re.search("only", line[0]):
                    name = re.split('\s+', line[0])
                    only[(name[0],)] = int(line[1])
                    single[(name[0],)] = int(line[1])
                if re.search('&', line[0]):
                    line[0] = re.sub('\s+', '', line[0])
                    name = re.split('&', line[0])
                    single_len[tuple(name)] = int(line[1])
                    sum_len[len(name)] += int(line[1])
                    num[tuple(name)] += int(line[1])
        # print num
        # print only
        # 因为程序运行出来的venn表当中只给了某个组别的only的值, 因此需要计算相应的all的值
        # 计算某一个组取并的时候的数目和，例如如果有A,B,C三组，那就是计算所有的A当中有多少的OTU，B当中和C当中有多少的OTU
        for my_only in only:
            num[my_only] = only[my_only]
            for i in range(2, gp_len + 1):
                num[my_only] += ((-1) ** i) * sum_len[i]
                for s_len in single_len:
                    if len(s_len) == i and my_only[0] not in s_len:
                        num[my_only] = num[my_only] + ((-1) ** (i + 1) * single_len[s_len])

        # 计算所有组别取交的值
        for set_name in num:
            single[set_name] = num[set_name]
            c = 1
            for i in range(len(set_name) + 1, gp_len + 1):
                single[set_name] += ((-1) ** c) * sum_len[i]
                for new_set in num:
                    if len(new_set) == i and (not set(set_name) < set(new_set)):
                        single[set_name] = single[set_name] + ((-1) ** (c + 1) * num[new_set])
                c += 1

        # 为了Venn图美观，平均化单个的大小， 对其他部分的大小进行缩小
        avg = 0
        c = 0
        print num
        print single
        self.num = num
        self.single = single
        for name in num:
            if len(name) == 1:
                avg += num[name]
                c += 1
        avg = avg / c

        for name in num:
            num[name] = avg / len(name)

        # print num
        venn_json = list()
        with open(venn_path, 'rb') as r:
            for line in r:
                line = line.strip('\n')
                strline = line
                line = re.split("\t", line)
                tmp_label = line[0]
                tmp_list = list()
                if re.search("only", line[0]):
                    tmp_label = re.sub("only", "", tmp_label)
                    tmp_label = re.sub(" ", "", tmp_label)
                    name = re.split('\s+', line[0])[0]
                    """
                    sets = {"sets": [name]}
                    size = {"size": num[name]}
                    label = {"label": tmp_label}
                    tmp_list = [sets, size, label]
                    """
                    tmp_list = {"sets": [name], "size": num[(name,)], "label": tmp_label}
                elif re.search("&", line[0]):
                    line[0] = re.sub('\s+', '', line[0])
                    name = re.split('&', line[0])
                    """
                    sets = {'sets': name}
                    size = {"size": num[tuple(name)]}
                    label = {"label": tmp_label}
                    tmp_list = [sets, size, label]
                    """
                    tmp_list = {"sets": name, "size": num[tuple(name)], "label": ""}
                else:
                    raise Exception("Venn 表格中行{}无法解析".format(strline))
                venn_json.append(tmp_list)
        return venn_json

    def _add_sg_otu_detail(self, info, from_otu_id, new_otu_id, title):
        """
        对otu表进行筛选，当它符合venn_table里的结果时，将他输出到sub_otu_path中，形成一张otu子表，读取子表
        删除值全部为0的列，形成no_zero_path中的otu表，最后将这张表导入数据库中(sg_otu_detail, sg_speciem)
        """
        title = re.sub(r'\s', '', title)
        sub_otu_dir = os.path.join(self.bind_object.work_dir, "sub_otu")
        sub_otu_path = os.path.join(sub_otu_dir, title + ".sub")
        no_zero_path = os.path.join(sub_otu_dir, title + ".no_zero")
        selected_clas = re.split(',', info)
        o_otu_path = self.all_otu
        sample_num = defaultdict(int)
        level = self.bind_object.option("level")
        with open(o_otu_path, 'rb') as r, open(sub_otu_path, 'wb') as w:
            head = r.next()
            w.write(head)
            head = head.strip('\n\r')
            head = re.split('\t', head)
            for line in r:
                line = line.strip('\r\n')
                line = re.split('\t', line)
                full_classify = re.split("; ", line[0])
                venn_classsify = full_classify[0:level]
                str_ = "; ".join(venn_classsify)
                if str_ in selected_clas:
                    new_line = "\t".join(line)
                    w.write(new_line + "\n")
                    for i in range(1, len(line)):
                        sample_num[i] += int(line[i])
        # print sample_num
        new_head = self._del_zero_column(sub_otu_path, no_zero_path, sample_num, head)
        self._table_to_sg_otu_detail(no_zero_path, new_otu_id)
        self._add_sg_otu_specimen(from_otu_id, new_otu_id, new_head)

    def _del_zero_column(self, sub_otu_path, no_zero_path, sample_num, head):
        index_list = list()
        new_head = list()  # 也就是样本名的列表（除去了head[0]和为0的列head）
        for i in sample_num:
            if sample_num[i]:
                index_list.append(i)
                new_head.append(head[i])
        with open(sub_otu_path, 'rb') as r, open(no_zero_path, 'wb') as w:
            w.write("OTU ID\t" + "\t".join(new_head) + "\n")
            line = r.next()
            for line in r:
                line = line.strip('\n\r')
                line = re.split('\t', line)
                w.write(line[0] + '\t')
                tmp_line = list()
                for i in index_list:
                    tmp_line.append(line[i])
                w.write('\t'.join(tmp_line) + "\n")
        return new_head

    def _table_to_sg_otu_detail(self, no_zero_path, new_otu_id):
        data_list = list()
        with open(no_zero_path, 'rb') as r:
            line = r.next().strip('\r\n')
            head = re.split('\t', line)
            for line in r:
                insert_data = dict()
                line = line.strip('\r\n')
                line = re.split('\t', line)
                classify_list = re.split(r"\s*;\s*", line[0])
                for c in classify_list:
                    insert_data[c[0:3]] = c
                insert_data["otu_id"] = new_otu_id
                insert_data["task_id"] = self.task_id
                for i in range(1, len(line)):
                    insert_data[head[i]] = line[i]
                data_list.append(insert_data)
        collection = self.db['sg_otu_detail']
        if len(data_list) > 0:
            collection.insert_many(data_list)

    def _add_sg_otu(self, otu_id, name):
        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        if re.search(r'only', name):
            name = re.sub(' ', "_", name)
        if re.search(r'&', name):
            name = re.sub(' ', '', name)
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "from_id": str(otu_id),
            "name": "venn_otu_" + name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "status": "end",
            "show": 0,
            "type": "otu_venn",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db['sg_otu']
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def _find_info(self, otu_id):
        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        collection = self.db['sg_otu']
        result = collection.find_one({'_id': otu_id})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(otu_id))
            self.bind_object.set_error("sg_otu表中找不到相应的记录", code="51007803")
        self.project_sn = result['project_sn']
        self.task_id = result['task_id']

    def _add_sg_otu_specimen(self, otu_id, new_otu_id, samples):
        """
        添加sg_otu_specimen记录
        """
        matched_sample_id = list()
        data_list = list()
        collection = self.db['sg_otu_specimen']
        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        results = collection.find({"otu_id": otu_id})
        for result in results:
            my_collection = self.db["sg_specimen"]
            try:
                my_result = my_collection.find_one({"_id": result['specimen_id']})
            except:
                self.bind_object.logger.error("样本id:{}在sg_specimen表里未找到对应的记录".format(result['specimen_id']))
                self.bind_object.set_error("sg_specimen表里未找到对应记录", code="51007805")
            if my_result["specimen_name"] in samples:
                matched_sample_id.append(result['specimen_id'])
        # print samples
        for m_id in matched_sample_id:
            insert_data = {
                "otu_id": new_otu_id,
                "specimen_id": m_id
            }
            data_list.append(insert_data)
        # print data_list
        if len(data_list) > 0:
            collection.insert_many(data_list)

    @report_check
    def add_venn_graph(self, venn_graph_path, venn_id, project='meta'):
        data_list = []
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                self.bind_object.set_error("venn_id必须为ObjectId对象或其对应的字符串!", code="51007806")
        with open(venn_graph_path, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                insert_data = {
                    'venn_id': venn_id,
                    # 'otu_id': ObjectId(otu_id),
                    'category_name': line[0]
                }
                if project == 'meta':
                    insert_data['otu_names'] = line[1]
                if project == 'denovo':
                    insert_data['gene_list'] = line[1]
                data_list.append(insert_data)
        try:
            if project == 'meta':
                collection = self.db["sg_otu_venn_graph"]
                main = self.db["sg_otu"]
                #main.update_one({'_id': venn_id},{"$set": {"main_id": venn_id}})
            if project == 'denovo':
                collection = self.db["sg_denovo_venn_graph"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Venn画图数据出错:%s" % e)
            self.bind_object.set_error("导入Venn画图数据出错", code="51007807")
        else:
            self.bind_object.logger.error("导入Venn画图数据成功")

    @report_check
    def add_denovo_venn(self, express_id, venn_table=None, venn_graph_path=None, params=None):
        # self._db_name = Config().MONGODB + '_rna'
        self._project_type = 'refrna'
        if not isinstance(express_id, ObjectId):
            express_id = ObjectId(express_id)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "express_id": str(express_id),
            "name": "Venn_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "params": params
        }
        collection = self.db['sg_denovo_venn']
        inserted_id = collection.insert_one(insert_data).inserted_id
        if venn_table:
            self.add_denovo_venn_detail(venn_table, inserted_id)
        if venn_graph_path:
            self.add_venn_graph(venn_graph_path=venn_graph_path, venn_id=inserted_id, project='denovo')
        return inserted_id

    @report_check
    def add_denovo_venn_detail(self, venn_table, venn_id):
        data_list = []
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                self.bind_object.set_error("venn_id必须为ObjectId对象或其对应的字符串!", code="51007808")
        with open(venn_table, "r") as f:
            for line in f:
                line = line.strip('\n').split("\t")
                insert_data = {
                    'venn_id': venn_id,
                    'label': line[0],
                    'num': line[1],
                    'gene_ids': line[2]
                }
                data_list.append(insert_data)
        try:
            collection = self.db["sg_denovo_venn_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Venn数据出错:%s" % e)
            self.bind_object.set_error("导入Venn数据出错", code="51007809")
        else:
            self.bind_object.logger.error("导入Venn数据成功")

    @report_check
    def add_venn_pie(self, venn_table, venn_id, asv_table=None, group_table=None):
        """
        需要计算每个标签所含有的物种数，
        :param venn_table:
        :param venn_id:
        :param asv_table:
        :param group_table:
        :return:
        """
        data_list = []
        if not isinstance(venn_id, ObjectId):
            if isinstance(venn_id, StringTypes):
                venn_id = ObjectId(venn_id)
            else:
                self.bind_object.set_error("venn_id必须为ObjectId对象或其对应的字符串!")
        data = pd.read_table(asv_table, sep="\t", header=0)##分组名称
        columns = list(data.columns)
        data["level"] = (data["OTU ID"].str.split(";").str)[-1].str.strip()
        new_columns = ["level"] + columns[1:]
        all_data = data[new_columns]
        all_data = all_data.set_index("level")
        with open(venn_table, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                label = line[0].split(" only")[0].strip()
                label_list = label.split("&")
                label_list = [x.strip() for x in label_list]
                label_list = [x.split(" only")[0].strip() for x in label_list]
                try:
                    tmp_name = re.split(',', line[2])
                except:
                    tmp_name = ""
                name_list = []
                for cla_info in tmp_name:
                    cla_info = re.split('; ', cla_info)
                    my_name = cla_info[-1].strip()
                    if my_name not in name_list:
                        name_list.append(my_name)
                if tmp_name != "":
                    species_dict = {}
                    all_species_number = 0
                    for sp_name in name_list:
                        spe_number = 0
                        for group in label_list:
                            try:
                                sp_num = all_data[group].loc[sp_name] ###可能会出现某个水平上的值是有相同名称的（这种是bug，taxon可能是错误的，没有加上级水平的名称）
                                spe_number += int(sp_num)
                            except:
                                sp_num = all_data[group].loc[sp_name] ###可能会出现某个水平上的值是有相同名称的（这种是bug，taxon可能是错误的，没有加上级水平的名称）
                                spe_number += int(max(list(sp_num)))
                        all_species_number += spe_number
                        if sp_name not in species_dict:
                            species_dict[sp_name] = spe_number
                    others_number = 0
                    for sp_name in name_list:
                        spe_number = species_dict[sp_name]
                        percent = float(spe_number) / all_species_number
                        if int(spe_number) != 0 and (percent >= 0.0001):
                            insert_data = {
                                'venn_id': venn_id,
                                'specimen': label,
                                'species_name': sp_name,
                                'species_abu':spe_number
                            }
                            data_list.append(insert_data)
                        else:
                            others_number += spe_number
                    insert_data = {
                        'venn_id': venn_id,
                        'specimen': label,
                        'species_name': "others",
                        'species_abu':others_number
                    }
                    data_list.append(insert_data)
                else:
                    insert_data = {
                            'venn_id': venn_id,
                            'specimen': label,
                            'species_name': "",
                            'species_abu':0
                        }
                    data_list.append(insert_data)
        try:
            collection = self.db["sg_otu_venn_pie"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Venn数据出错:%s" % e)
            self.bind_object.set_error("导入Venn数据出错")
        else:
            self.bind_object.logger.error("导入Venn数据成功")
