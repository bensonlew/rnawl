# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base
from bson import ObjectId
import json
import types
import datetime
from bson.son import SON
import pandas as pd
import copy
import re


class AnnoOverview(Base):
    def __init__(self, work_type, organism=False):
        super(AnnoOverview, self).__init__()
        self._project_type = 'metabolome'
        self.unkown_list = ["", " ", "-", "_", "\r"]
        self.norm_unkown_list = ["-"]
        self.work_type = work_type
        self.table_info_dic = {}
        self.pathway_level_info = {}
        self.organism = organism

         #输入的表头，和输出文件的表头的对应关系
        self.name_map = {
            'ID':'ID',
            'm/z' :'m/z',
            'RT (min)':'Retention time',
            'Retention time':'Retention time',
            'Adducts':'Adducts',
            'Formula' : 'formula',
            'HMDB_ID':'hmdb_id',
            'Library ID':'hmdb_id',
            'CAS number' : 'cas_id',
            'CAS ID' : 'cas_id',
            'Similarity' : "similarity", # gc add score v3
            'Fragmentation Score' : 'Fragmentation Score',
            'Theoretical Fragmentation Score' : 'Theoretical Fragmentation Score',
            'Score' : 'Score',
        }
        self.mongo_map = self.name_map

    def init_table_info(self, metab_table, set_file, from_mongo=True):
        # metab_table can be metab_table_id or table_file.In condition one ,set from_mongo=True,
        # otherwise set the from_mongo=False.If the work_type is "LC", table_file should be three files seperated
        # with "," and the order must be pos/neg/mix
        if from_mongo:
            if not isinstance(metab_table, ObjectId):
                if isinstance(metab_table, types.StringTypes):
                    metab_table = ObjectId(metab_table)
                else:
                    raise Exception('metab_table必须为ObjectId对象或其对应的字符串！')
            self.init_table_info_from_mongo(metab_table)
            # 目前from_mongo读mix表，而原始表是没有mix表的
        elif self.work_type == "LC":
            table_list = metab_table.split(",")
            if len(table_list) != 2:
                raise Exception("LC分析传入的数据表详情需要两个，并且逗号分割，目前文件个数：%s" % len(table_list))
            self.init_lc_table_info_from_table(table_list, set_file)
        elif self.work_type == "GC":
            table_list = metab_table.split(",")
            if len(table_list) != 1:
                raise Exception("GC分析传入的数据表详情需要一个，目前文件个数：%s" % len(table_list))
            self.init_gc_table_info_from_table(table_list)
        else:
            raise Exception("work_type 必须为LC/GC,不可以为: %s" % self.work_type)

    def init_table_info_from_mongo(self, metab_table):
        # get information from mongo, disabled now
        if self.work_type == "LC":
            # collection = self.db['metab_table_mix']
            collection = self.db['exp_mix']
        elif self.work_type == "GC":
            # collection = self.db['metab_table_pos']
            collection = self.db['exp_pos']
        result = collection.find({"exp_id": metab_table})
        for one in result:
            # self.update_cas_id(collection, one['metab_id'])  # 暂时更新
            self.table_info_dic[one['exp_id']] = {
                "metab": one['metab'],
                "hmdb_id": one["hmdb_id"],
                "cas_id": one["cas_id"],  # 暂时更新
                "formula": one["formula"],
                "element": one["element"],
                "compound_id": one['kegg_id']
            }
            if self.work_type == "GC":
                self.table_info_dic[one["exp_id"]]["mode"] = one['mode']
            elif self.work_type == "LC":
                self.table_info_dic[one["exp_id"]]["mode"] = self.treat_mode(metab_table, one['exp_id'])
    """
    def update_cas_id(self, mix_coll, metab_id):  # 暂时更新
        check_result = mix_coll.find_one({"metab_id": metab_id})
        cas_id = check_result['cas_id']
        if cas_id != "-":
            cas_id = cas_id.replace("\r","")
            mix_coll.update({"metab_id": metab_id}, {"$set": {"cas_id": cas_id}},multi=True)
            pos_coll = self.db['metab_table_pos']
            neg_coll = self.db['metab_table_neg']
            pos_coll.update({"metab_id": metab_id}, {"$set": {"cas_id": cas_id}},multi=True)
            neg_coll.update({"metab_id": metab_id}, {"$set": {"cas_id": cas_id}},multi=True)
    """

    def treat_mode(self, metab_table, metab_id):
        # treat mode info , disabled now
        pos_col = self.db['exp_pos']
        neg_col = self.db['exp_neg']
        pos_result = pos_col.find_one({"exp_id": metab_table, "metab_id": metab_id})
        neg_result = neg_col.find_one({"exp_id": metab_table, "metab_id": metab_id})
        if pos_result and neg_result:
            return "pos,neg"
        elif pos_result:
            return "pos"
        elif neg_result:
            return "neg"
        else:
            raise Exception("数据表中找不到metab_id为%s的代谢物，主表：%s" % (metab_id, metab_table))

    def init_lc_table_info_from_table(self, table_list, set_file):
        # get information from table, in lc condition
        data = pd.read_table(set_file, header=None, index_col=0)
        pos_data = pd.read_table(table_list[0], index_col=['metab_id'])
        pos_index = pos_data.index
        neg_data = pd.read_table(table_list[1], index_col=['metab_id'])
        neg_index = neg_data.index

        var_heads = ['ID','m/z','RT (min)','Retention time','Adducts','Formula','HMDB_ID','Library ID','CAS number','CAS ID',
                     "Theoretical Fragmentation Score", "Fragmentation Score"]
        current_v = []  # 阴阳离子表 表头是一样的，故 current_p_v 和 current_n_v 统一成 current_v

        for v in var_heads:
            if v in pos_data.columns:
                current_v.append(v)
        self.current_v = current_v
        for metab_id in data.index:

            if metab_id in pos_index:
                 self.table_info_dic[metab_id] = {
                    "metab": pos_data["Metabolite"][metab_id],
                    #"hmdb_id": pos_data[h_name][metab_id],
                    #"cas_id": pos_data[cas_name][metab_id],
                    #"formula": pos_data["Formula"][metab_id],
                    #"element": self.element(pos_data["Formula"][metab_id]),
                    "compound_id": pos_data["KEGG Compound ID"][metab_id],
                    "mode": "pos"
                }
                 for k in self.current_v:
                     mk = self.mongo_map[k]
                     self.table_info_dic[metab_id][mk] = pos_data[k][metab_id]
                     if k == 'Formula':
                         self.table_info_dic[metab_id]['element'] = self.element(pos_data[k][metab_id])

            elif metab_id in neg_index:
                 self.table_info_dic[metab_id] = {
                    "metab": neg_data["Metabolite"][metab_id],
                    #"hmdb_id": neg_data[h_name][metab_id],
                    #"cas_id": neg_data[cas_name][metab_id],
                    #"formula": neg_data["Formula"][metab_id],
                    #"element": self.element(neg_data["Formula"][metab_id]),
                    "compound_id": neg_data["KEGG Compound ID"][metab_id],
                    "mode": "neg"
                }
                 for k in self.current_v:
                     mk = self.mongo_map[k]
                     self.table_info_dic[metab_id][mk] = neg_data[k][metab_id]
                     if k == 'Formula':
                         self.table_info_dic[metab_id]['element'] = self.element(neg_data[k][metab_id])

            else:
                raise Exception("数据表中找不到metab_id为%s的代谢物，表格路径：%s" % (metab_id, table_list))


    def init_gc_table_info_from_table(self, table_list):
        # get information from table in gc condition
        if self.work_type == "GC":
            data = pd.read_table(table_list[0], index_col=['metab_id']).fillna("-")
            var_heads = ['ID','m/z','RT (min)','Retention time','Adducts','Formula','HMDB_ID','Library ID','CAS number','CAS ID','Similarity',
                         "Score"]
            current_v = []
            for v in var_heads:
                if v in data.columns:
                    current_v.append(v)    #20190617 zouguanqing
            self.current_v = current_v
            for metab_id in data.index:
                self.table_info_dic[metab_id] = {
                    "metab": data["Metabolite"][metab_id],
                    #"hmdb_id": data["HMDB_ID"][metab_id],
                    #"cas_id": data["CAS number"][metab_id],
                    #"formula": data["Formula"][metab_id],
                    #"element": self.element(data["Formula"][metab_id]),
                    "compound_id": data["KEGG Compound ID"][metab_id],
                    "mode": data['Mode'][metab_id]
                }
                for v in self.current_v:
                    mv = self.mongo_map[v]
                    self.table_info_dic[metab_id][mv] = data[v][metab_id]
                    if v == 'Formula':
                         self.table_info_dic[metab_id]['element'] = self.element(data[v][metab_id])


    def check_info(self, info):
        if info in self.unkown_list:
            return "-"
        else:
            return info

    def element(self, formula):
        formula = self.check_info(formula)
        if formula == "-":
            return formula
        element_list = []
        for i in formula:
            if i in "+.[]-(); ":
                forward_letter = i
                continue
            if i == 'n':
                if forward_letter in [']', ')']:
                    continue
                elif forward_letter.isupper():
                    element_list[-1] += i
                else:
                    raise Exception("化合物中含有未处理的情况，需注意：%s" % formula)
            try:
                int(i)  # 去掉数字
                forward_letter = i
            except:
                forward_letter = i  # 用于判断n是在]后还是在元素后
                if i.isupper():
                    element_list.append(i)
                elif i.islower():
                    element_list[-1] += i
                else:
                    raise Exception('分子式%s中发现特殊字符%s' % (formula, i))
        element_set = set(element_list)
        return_str = ':'.join(list(element_set))
        return_str = ':' + return_str + ':'
        return return_str

    def get_table_info(self, metab_id):
        one_metab_dic = self.table_info_dic[metab_id]
        return_list = []
        for one_key in ["metab", "mode", "compound_id"]:  #, "hmdb_id", "cas_id", "formula", "element",
            info = one_metab_dic[one_key]
            info = self.check_info(info)
            return_list.append(info)
            # return_list.append(one_metab_dic[one_key])
        return return_list

    def get_var_head_info(self,metab_id):   # zouguanqing 20190617
        one_metab_dic = self.table_info_dic[metab_id]
        return_list = []
        for v in self.current_v:
            mv = self.mongo_map[v]
            return_list.append([mv,one_metab_dic[mv]])
            if mv == 'formula':
                return_list.append(['element',one_metab_dic['element']])
        return return_list


    def get_hmdb_info(self, hmdb_id):
        if hmdb_id in self.unkown_list:
            return self.norm_unkown_list * 2
        hmdb_db = self.ref_db['hmdb']
        try:
            hmdb_id + ";"
        except:
            raise Exception("check: %s" % hmdb_id)
        hmdb_ids = hmdb_id.split(';')
        hmdb_names = []
        masss = []
        for hmdb_id in hmdb_ids:
            if not  re.match('^HMDB\d+$',hmdb_id):
                continue
            result = hmdb_db.find_one({"hmdb_id": hmdb_id})  # 改为查询hmdb_id字段 by ghd @20191016
            if not result:
                result = hmdb_db.find_one({"second_id": {"$regex": hmdb_id + ";"}})
                if not result:
                    print "需注意hmdb_id%s不存在" % hmdb_id
                    #return "-","-"
                    hmdb_name = '-'
                    mass = '-'
                else:
                    print "需注意hmdb_id%s是二级分类" % hmdb_id
                    continue

            else:
                hmdb_name = result['name']
                mass = result['exact_mass']  # 用质量值
            hmdb_names.append(hmdb_name)
            masss.append(mass)

        return ';'.join(hmdb_names), ';'.join(masss)

    def get_compound_info(self, compound_id):
        if compound_id in self.unkown_list:
            return self.norm_unkown_list * 4
        compound_db = self.ref_db['kegg_compound']
        result = compound_db.find_one({"entry": compound_id})
        if not result:
            return self.norm_unkown_list * 4
        compound_name = result['name']
        compound_first_category = result['br1_category']
        compound_second_category = result['br2_category']
        pathway_id = result['pathway']
        pathway_id = self.check_pathway(pathway_id)
        return compound_name, compound_first_category, compound_second_category, pathway_id

    def check_pathway(self, pathway_id):
        if self.organism in [None, "False", "All"]:
            return pathway_id
        pathway_list = set(pathway_id.split(";"))
        organism_pathway_list = self.get_org_pathway_list()
        new_set = set(pathway_list) & set(organism_pathway_list)
        new_pathway_id = ';'.join(list(new_set))
        return new_pathway_id

    def get_org_pathway_list(self):
        """
        根据kegg_organism查询相关物种的通路信息
        :return:
        """
        tmp_list = self.organism.split(";")
        ref_org_db = self.ref_db["kegg_organisms"]
        ko_list = []
        if len(tmp_list) == 1:
            ko_set = set()
            results = ref_org_db.find({"first_category": tmp_list[0]})
            for result in results:
                tmp_ko_list = result['map_list'].split(";")
                ko_set = ko_set | set(tmp_ko_list)
            ko_list = list(ko_set)
        elif len(tmp_list) == 2:
            if tmp_list[0] == 'All':
                result = ref_org_db.find_one({"second_category": tmp_list[1]})
            else:
                result = ref_org_db.find_one({"first_category": tmp_list[0], "second_category": tmp_list[1]})
            ko_list = result['map_list'].split(";")
        ko_list = map(lambda x: 'map' + x, ko_list)
        return ko_list

    def get_kegg_info(self, pathway_id):
        if pathway_id in self.unkown_list:
            return self.norm_unkown_list * 3
        kegg_path_db = self.ref_db['kegg_pathway_level']
        pathway_id = pathway_id.split(";")
        description, kegg_first_category, kegg_second_category = [[], [], []]
        for one_path in pathway_id:
            result = kegg_path_db.find_one({"pathway_id": one_path})
            description.append(result['discription'])
            kegg_first_category.append(result["first_category"])
            kegg_second_category.append(result["second_category"])
            if self.pathway_level_info.has_key(one_path):
                continue
            self.pathway_level_info[one_path] = {
                "description": result['discription'],
                "first_category": result['first_category'],
                "second_category": result['second_category']
            }
        description = ';'.join(description)
        kegg_first_category = ';'.join(kegg_first_category)
        kegg_second_category = ';'.join(kegg_second_category)
        return description, kegg_first_category, kegg_second_category

    def check_element(self, category):
        # 用于检测分类名称中是否是一样的，如果是一样的，则归为一种 by ghd @20191016
        check_catogory = set(category)
        if len(check_catogory) == 1:
            return list(check_catogory)
        else:
            return category

    def make_table(self, set_file, output, ko_output):

        table_head = ["metab_id", "metab", "mode", "hmdb_id", "hmdb_name", "mass", "cas_id", "formula", "element",
                      "compound_id", "compound_name", "compound_first_category", "compound_second_category",
                      "pathway_id", "description", "kegg_first_category", "kegg_second_category", "similarity"]

        table_head.extend(['ID','m/z','Retention time','Adducts', "Theoretical Fragmentation Score", "Fragmentation Score", "Score"])
        current_table_head = copy.deepcopy(table_head)


        if 'RT (min)' not in self.current_v and 'Retention time' not in  self.current_v:
            if 'Retention time' in current_table_head:
                current_table_head.remove('Retention time')
        if  'CAS number' not in self.current_v and 'CAS ID' not in  self.current_v:
            if  'cas_id' in current_table_head:
                current_table_head.remove('cas_id')

        if 'Library ID' not in self.current_v and 'HMDB_ID' not in self.current_v:
            if 'hmdb_id' in current_table_head:
                current_table_head.remove('hmdb_id')

        for name in self.name_map:
            if name in ['RT (min)','Retention time','CAS number','CAS ID','Library ID','HMDB_ID']:
                continue
            if name not in self.current_v:
                if  self.name_map[name] in current_table_head:
                    current_table_head.remove(self.name_map[name])

        if 'hmdb_id' not in current_table_head:
            if 'hmdb_name' in current_table_head:
                current_table_head.remove('hmdb_name')
            if 'mass' in current_table_head:
                current_table_head.remove('mass')
        if 'formula' not in current_table_head:
            if 'element' in current_table_head:
                current_table_head.remove('element')

        table_dic = {}
        for v in current_table_head:
            table_dic[v] = []
        print table_dic.keys()

        set_data = pd.read_table(set_file, header=None)
        set_list = set_data[0].tolist()
        for metab_id in set_list:
            var_info_list = self.get_var_head_info(metab_id)
            for vk, vv in var_info_list:
                table_dic[vk].append(vv)
                if vk == 'hmdb_id':
                    hmdb_name, mass = self.get_hmdb_info(vv)  # 查参考库
                    table_dic["hmdb_name"].append(hmdb_name)
                    table_dic["mass"].append(mass)

            metab, mode, compound_id_list_str = self.get_table_info(metab_id)  #hmdb_id, cas_id, formula, element,
            compound_name_list = []
            compound_first_category_list = []
            compound_second_category_list = []
            pathway_id_set = set()
            pathway_info = {}
            for compound_id in compound_id_list_str.split(";"):
                compound_name, compound_first_category, compound_second_category, pathway_id = self.get_compound_info(
                    compound_id)
                compound_name_list.append(compound_name)
                pathway_info[compound_id] = pathway_id
                compound_first_category_list.append(compound_first_category)
                compound_second_category_list.append(compound_second_category)
                if pathway_id == "":
                    pathway_id = "-"
                if pathway_id != '-':
                    pathway_id_set.update(pathway_id.split(';'))
            if len(pathway_id_set) == 0:
                pathway_id = "-"
            else:
                pathway_id = ";".join(list(pathway_id_set))
            description, kegg_first_category, kegg_second_category = self.get_kegg_info(pathway_id)
            for c_id in pathway_info.keys():
                for one_path in pathway_info[c_id].split(";"):
                    if self.pathway_level_info.has_key(one_path):
                        if self.pathway_level_info[one_path].has_key("metab_id"):
                            self.pathway_level_info[one_path]["metab_id"].append(metab_id)
                            self.pathway_level_info[one_path]["compound_id"].append(c_id)
                        else:
                            self.pathway_level_info[one_path]["metab_id"] = [metab_id]
                            self.pathway_level_info[one_path]["compound_id"] = [c_id]
            table_dic['metab_id'].append(metab_id)
            table_dic["metab"].append(metab)
            table_dic["mode"].append(mode)
            table_dic["compound_id"].append(compound_id_list_str)  # 改成多个
            table_dic["compound_name"].append(";".join(compound_name_list)) # 改成多个
            table_dic["compound_first_category"].append(";".join(self.check_element(compound_first_category_list))) # 改成多个
            table_dic["compound_second_category"].append(";".join(self.check_element(compound_second_category_list))) # 改成多个
            table_dic["pathway_id"].append(pathway_id)
            table_dic["description"].append(description)
            table_dic["kegg_first_category"].append(kegg_first_category)
            table_dic["kegg_second_category"].append(kegg_second_category)

        print len(current_table_head)
        print len(table_dic.keys())
        #print table_dic
        for k in table_dic.keys():
            print k
            print len(table_dic[k])

        output_data = pd.DataFrame(data=table_dic, columns=current_table_head)
        output_data.to_csv(output, index=False, sep="\t", quoting=3)
        # ko_output = output + "_ko.xls"
        output_ko = open(ko_output, "w")
        output_ko.write("pathway_id\tdescription\tfirst_category\tsecond_category\tcompound_id\tmetab_id\thyperlink\n")
        for ko in sorted(self.pathway_level_info.keys()):
            metab_id = ';'.join(self.pathway_level_info[ko]["metab_id"])
            compound_id = ";".join(self.pathway_level_info[ko]["compound_id"])
            hyperlink = "http://www.genome.jp/dbget-bin/show_pathway?" + ko
            for c in set(self.pathway_level_info[ko]["compound_id"]):
                hyperlink += "+" + c + "%09green"
            output_ko.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ko, self.pathway_level_info[ko]['description'], self.pathway_level_info[ko]['first_category'], self.pathway_level_info[ko]['second_category'], compound_id, metab_id, hyperlink))
        output_ko.close()

    def run(self, metab_table, set_file, output, ko_output, from_mongo=True):
        self.init_table_info(metab_table, set_file, from_mongo=from_mongo)
        self.make_table(set_file, output, ko_output)
