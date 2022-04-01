# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20171223


import sys
import argparse
import os
import pandas as pd
from pymongo import MongoClient
#from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from biocluster.config import Config
from collections import Counter
import time,re
from mbio.packages.statistical.large_data import table_parser


class Overview(object):
    def __init__(self, db_version=None):
        Config().DBVersion = db_version
        self.client = Config().get_mongo_client(mtype="metagenomic")
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic")]
        self.ref_client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.ref_mongodb = self.ref_client[Config().get_mongo_dbname("metagenomic", ref=True)]
        #self.ref_mongodb = Config().mongo_client.sanger_biodb
        self.overview_main = self.mongodb.anno_overview
        self.overview = self.mongodb.anno_overview_detail
        self.db_path = os.path.join(Config().SOFTWARE_DIR, "database/metagenome")
        self.level_name_file = os.path.join(self.db_path, "database_id_to_file.xls")
        self.anno_nr = self.mongodb.anno_nr
        self.anno_nr_detail = self.mongodb.anno_nr_detail
        self.level_correspond = {
            "nr_1": "Domain",
            "nr_2": "Kingdom",
            "nr_3": "Phylum",
            "nr_4": "Class",
            "nr_5": "Order",
            "nr_6": "Family",
            "nr_7": "Genus",
            "nr_8": "Species",
            "query": "#Query",
            "geneset_id": "#Query",
            "cog_9": "COG_Category",
            "cog_10": "COG_Function",
            "cog_11": "NOG",
            "kegg_12": "KEGG_Level1",
            "kegg_13": "KEGG_Level2",
            "kegg_14": "KEGG_Pathway",
            "kegg_15": "KEGG_Modules",
            "kegg_16": "KEGG_Enzyme",
            "kegg_17": "KO",
            "cazy_18": "CAZY_Cl_description",
            "cazy_19": "CAZY_Family",
            "ardb_20": "ARDB_Class",
            "ardb_21": "ARDB_Type",
            "ardb_22": "ARDB_Antibiotic_type",
            "ardb_23": "ARDB_ARG",
            "card_24": "CARD_Class",
            "card_25": "CARD_ARO",
            "vfdb_26": "VFDB_Level1",
            "vfdb_27": "VFDB_Level2",
            "vfdb_28": "VFDB_VFs",
            #####
            "lca_nr_1": "lca_Domain",
            "lca_nr_2": "lca_Kingdom",
            "lca_nr_3": "lca_Phylum",
            "lca_nr_4": "lca_Class",
            "lca_nr_5": "lca_Order",
            "lca_nr_6": "lca_Family",
            "lca_nr_7": "lca_Genus",
            "lca_nr_8": "lca_Species",
            "duc_nr_1": "duc_Domain",
            "duc_nr_2": "duc_Kingdom",
            "duc_nr_3": "duc_Phylum",
            "duc_nr_4": "duc_Class",
            "duc_nr_5": "duc_Order",
            "duc_nr_6": "duc_Family",
            "duc_nr_7": "duc_Genus",
            "duc_nr_8": "duc_Species",
            "pfam_29" : "pfam_Domain",
            "pfam_30" : "pfam_Type",
            "p450_33" : "p450_homo_family",
            "p450_34" : "p450_super_family",
            "mvir_37" : "Mvirdb_Virulence_Factor_Type",
            "mvir_36" : "Mvirdb_Short_Description",  #"Gene Description（species）",  !!!
            "phi_40" : "phi_Pathogen_Species",
            "phi_41" : "phi_Phenotype",
            "phi_41" : "phi_PHI_ID",
            "phi_47" : "phi_Host",
            "tcdb_48" : "TCDB_Family",
            "tcdb_49" : "TCDB_Subclass",
            "tcdb_50" : "TCDB_Class",
            "tcdb_51" : "TCDB_ID",
            "qs_52" : "qs_Class",
            "probio_53" : "probio_Probiotic_name",
            "probio_54" : "probio_Genus",
            "probio_55" : "probio_Use_in",
            "probio_66" : "probio_Probiotic_Effect",
            "go_62" : "go_GO_Term_(Lev4)",
            "go_61" : "go_GO_Term_(Lev3)",
            "go_60" : "go_GO_Term_(Lev2)",
            "go_59" : "go_GO_(Lev1)",
            "secpro_63" : "sec_Gram_neg",
            "secpro_64" : "sec_Gram_pos",
            "secpro_65" : "sec_Fungi",
            "card_70" : "CARD_Resistance_Mechanism",
            "card_71" : "CARD_Antibiotic_class",
            "card_72" : "CARD_ARO_name",
            "ardb_73" : "ARDB_Antibiotic_class",
            #"t3ss_67" : ""

        }

    def run_select_overview(self, select, database, of, outDir):
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        select = select.replace("regex", "$regex").replace("\"flags\":\"i\"", "").replace("\$", "$")
        select = eval(select)
        select["anno_overview_id"] = ObjectId(select["anno_overview_id"])
        overview = select["anno_overview_id"]
        databases = database.split(",")
        self.all_names = self.convert_id_to_name()
        if of != "":
            print "from overview"
            #self.file_from_overview(select, databases, of, outDir)  # modify by shaohua.yuan 20180424
            # self.select_file(select, databases, of, outDir,overview)
            of = pd.read_csv(of, sep='\t', header=0, chunksize=100000, dtype='string')
            out = []
            for each in of:
                each_out = self.select_file(select, databases, each, outDir, overview)
                out.append(each_out)
            self.check_and_out(pd.concat(out))

    def check_and_out(self, df):
        df = df.drop_duplicates('#Query')
        anno_cols = len(df.columns) - 2
        if not df.empty:
            df["count"] = df.apply(lambda x: sum(x == "-"), axis=1)
            df = df[df["count"] != anno_cols]
            df = df.drop('count', 1)
            df.to_csv(outDir + "/anno_overview.xls", sep="\t", index=False)
            gene_list = pd.DataFrame(df["#Query"])
            gene_list.columns = ["GeneID"]
            gene_list = gene_list.drop_duplicates()
            gene_list.to_csv(outDir + "/gene_list.xls", sep="\t", index=False)
        else:
            raise Exception("该筛选条件未找到基因集！")

    def select_file(self, select, databases, of, outDir,anno_overview_id):
        print select
        # all_names = self.convert_id_to_name()
        if "CARD_Drug_Class" in of.columns:
            self.new_card = True
        else:
            self.new_card = False
        all_names = self.all_names
        if "anno_overview_id" in select:
            anno_overview_id = select["anno_overview_id"]
            print anno_overview_id
        else:
            raise Exception("未传入anno_overview主表id！")
        if len(select.keys()) == 1:
            select_file = of
        else:
            select_file = of
            if "length" in select:
                length_filter = select["length"]
                length = length_filter["$gte"]
                select_file = select_file[select_file["Length"] >= float(length)]
                select.pop("length")
            if "$and" in select:
                and_all = select["$and"]
                for each in and_all:
                    select_file = self.filter_level(all_names, select_file, each, overview_id=anno_overview_id)
            else:
                select.pop("anno_overview_id")
                if select:
                    # print "select:",select
                    select_file = self.filter_level(all_names, select_file, select, overview_id=anno_overview_id)
        data_cloumns = self.database_to_cloumns(databases)
        print(data_cloumns)
        print(select_file.columns)
        if "KEGG name" not in select_file.columns and "KEGG name" in data_cloumns:
            data_cloumns.remove("KEGG name")
        if not len(Counter(data_cloumns) - Counter(select_file.columns)) == 0:
            raise Exception("所选展示database列中有名称不包含于overview总表中")
        else:
            final = select_file[data_cloumns].drop_duplicates('#Query')
        return final

    def run_select_from_overviewfile(self, overviewfile, select, outDir, task_id, select_type="2"):
        select = eval(select)
        print select
        all_names = self.convert_id_to_name()
        print overviewfile
        of = pd.read_csv(overviewfile, sep='\t', header=0, chunksize=100000, dtype='string')
        out = []
        for select_file in of:
            print('select_file.head(): {}'.format(select_file.head()))
            tmp = None
            if "length" in select:
                length_filter = select["length"]
                length = length_filter["$gte"]
                select_file = select_file[select_file["Length"] >= float(length)]
                select.pop("length")
            if select_type == '1':
                tmp_list = []
                for each in select:
                    print ">>>>>>>select:", select
                    print ">>>each:", each
                    select_file1 = self.decompose_level(select_type, all_names, select_file, each, select[each], task_id=task_id)
                    tmp_list.append(select_file1)
                tmp = pd.concat(tmp_list)
            elif select_type == '2':
                for each in select:
                    print ">>>>>>>select:", select
                    print ">>>each:", each
                    select_file = self.decompose_level(select_type, all_names, select_file, each, select[each], task_id=task_id)
                    if select_file.empty:
                        break
                tmp = select_file
            else:
                raise Exception("type类型必须为1或2")
            out.append(tmp)
        self.check_and_out(pd.concat(out))

    def decompose_level(self, select_type, all_names, profile, database,level_dict, task_id=None):
        level_name_list = []
        merge_profile = pd.DataFrame(columns=["#Query"])
        if database == "nr":
            nr_method = level_dict.keys()[0]
            level = level_dict[nr_method].keys()
        elif database == "geneset_id":
            level = level_dict
        else:
            nr_method = "best_hit"
            level = level_dict.keys()
        print ">>>>>level%s" % level
        total_gene_id = []
        all_select_index = set()
        for i in level:
            print i
            if database == "nr":
                level_name_ids = level_dict[nr_method][i]
                level_total_name = database + "_" + i
            elif database == "geneset_id":
                level_name_ids = i
                level_total_name = "geneset_id"
            else:
                level_name_ids = level_dict[i]
                level_total_name = database + "_" + i
            act_name = self.level_correspond[level_total_name]
            if act_name == "#Query":
                level_name_list.append(i)
            else:
                for level_name_id in level_name_ids:
                    act_level_name = self.id_covert(all_names, i, level_name_id, task_id=task_id, nr_method=nr_method)
                    level_name_list.append(act_level_name)
            profile2 = profile[act_name].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename(
                act_name)
            profile2 = pd.DataFrame(profile2)
            print profile2.head()
            if select_type == "2":
                if "all" in level_name_list:
                    select_index = set(profile2[profile2[act_name] != "-"].index)
                else:
                    print ">>>>>>level_name_list:", level_name_list
                    tmp_index = []
                    for each_level_name in level_name_list:
                        if act_name in ['KEGG_Pathway', 'KEGG_Level2', 'KEGG_Level1']:
                            tmp_index.append(set(profile2[profile2[act_name].str.contains(each_level_name)].index))  # 每一轮筛选得到的index暂时存到一个列表中 by xieshichang 20200527
                        else:
                            tmp_index.append(set(profile2[profile2[act_name] == each_level_name].index))  # 每一轮筛选得到的index暂时存到一个列表中 by xieshichang 20200527

                    select_index = reduce(lambda x, y: x & y, tmp_index)  # select_type == 2时, 取交集
                    print ">>>select_index", select_index
                if not all_select_index:
                    all_select_index = select_index
                else:
                    all_select_index = all_select_index & select_index
            elif select_type == "1":
                if "all" in level_name_list:
                    level_name_list.pop["all"]
                if act_name == "#Query":
                    #tmp_profile = profile[profile["#Query"] == i]
                    total_gene_id.append(i)
                else:
                    if act_name in ['KEGG_Pathway', 'KEGG_Level2', 'KEGG_Level1']:
                        for x in level_name_list:
                            tmp_profile = profile[profile[act_name].str.contains(x)]
                            level_list = tmp_profile['#Query']
                            for j in level_list:
                                if j not in total_gene_id:
                                    total_gene_id.append(j)
                    else:
                        # tmp_profile = profile[profile[act_name].isin(level_name_list)]
                        # level_list = tmp_profile['#Query']
                        # for j in level_list:
                        #    total_gene_id.append(j)
                        # 改为在profile2中获得满足条件的index，然后在原始表profile中获得gene list
                        tmp_index = profile2[profile2[act_name].isin(level_name_list)].index
                        new_pro = profile.reset_index()
                        level_list = new_pro[new_pro['index'].isin(tmp_index)]['#Query']
                        total_gene_id.extend(list(level_list))
            else:
                raise Exception("select_type number must be 1 or 2！")
        if select_type == "2":
            select = profile.loc[all_select_index]
        elif select_type == "1":
            total_gene_id2 = set(total_gene_id)
            select = profile.loc[profile['#Query'].isin(total_gene_id2)]
        else:
            raise Exception("select_type number must be 1 or 2！")
        #print select.index
        #print select.head()
        return select

    def filter_level(self, all_names, profile, level_dict, overview_id=None):
        level_name_list = []
        level = level_dict.keys()
        print(type(profile))
        print(profile.head())
        for i in level:
            level_name_id = level_dict[i]
            act_name = self.level_correspond[i]
            if act_name == "#Query":
                act_level_name = level_name_id
            else:
                act_level_name = self.id_covert(all_names, i, level_name_id, overview=overview_id)
            level_name_list.append(act_level_name)
            profile2 = profile[act_name].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename(
                act_name)
            profile2 = pd.DataFrame(profile2)
            if "all" in level_name_list:
                select_index = profile2[profile2[act_name] != "-"].index
            else:
                select_index = profile2[profile2[act_name].isin(level_name_list)].index
        select = profile.loc[select_index]
        return select

    def id_covert(self, all_names, level_id, level_name_id, overview=None, task_id=None, nr_method="best_hit"):
        if level_id in ["nr_1", "nr_2", "nr_3", "nr_4", "nr_5", "nr_6", "nr_7", "nr_8"]:
            nr_dict1 = {"nr_1": "d_id", "nr_2": "k_id", "nr_3": "p_id", "nr_4": "c_id", "nr_5": "o_id", "nr_6": "f_id",
                        "nr_7": "g_id", "nr_8": "s_id"}
            nr_dict2 = {"nr_1": "d__", "nr_2": "k__", "nr_3": "p__", "nr_4": "c__", "nr_5": "o__", "nr_6": "f__",
                        "nr_7": "g__", "nr_8": "s__"}
            print nr_method
            if not task_id:
                task_id = self.overview_main.find_one({"_id": overview})["task_id"]
            if nr_method=="best_hit":
                origin_name = "NR_Origin"
            elif nr_method=="lca":
                origin_name = "NR_Origin_LCA"
            elif nr_method=="de_unclassified":
                origin_name = "NR_Origin_deunclassied"
            anno_nr_detail = self.anno_nr.find_one({"task_id": task_id, "status": "end", "name": origin_name})["_id"]
            nr_id = nr_dict1[level_id]
            nr_name = nr_dict2[level_id]
            name = self.anno_nr_detail.find_one({"nr_id": anno_nr_detail, "level_id": 8, nr_id: level_name_id})
            if name == None:
                raise Exception("id转化失败！")
            else:
                name = name[nr_name]
        else:
            if isinstance(level_name_id,dict):
                level_name_id = level_name_id.values()[0]
                string = re.compile("\.\*")
                level_name_id = string.sub('',level_name_id)
            if all_names.has_key(level_name_id):
                name = all_names[level_name_id]
            else:
                name = level_name_id
        return name

    def convert_id_to_name(self):
        name_convert = {}
        with open(self.level_name_file, "rb") as input:
            input.next()
            for line in input:
                line = line.strip().split("\t")
                id = line[0]
                name = line[1]
                if int(line[2]) == 10:
                    name = name.split(':')[0]
                elif int(line[2]) == 24:
                    name = line[0]
                if not name_convert.has_key(id):
                    name_convert[id] = name
        return name_convert

    def database_to_cloumns(self, databases):
        new_columns = ["#Query", "Length"]
        new_columns1 = []
        new_columns2 = []  # NR
        data_col = {}
        data_col["cog"] = ["NOG", "COG_description", "COG_Function", "COG_Fun_description", "COG_Category"]
        data_col["kegg"] = ["KEGG_Gene", "KO", "KO_description", "KEGG_Pathway", "KEGG_Path_description", "KEGG_Enzyme",
                            "KEGG_Modules", "KEGG_Level1", "KEGG_Level2"]
        data_col["cazy"] = ["CAZY_Family", "CAZY_Class", "CAZY_Cl_description"]
        data_col["ardb"] = ["ARDB_ARG", "ARDB_Type", "ARDB_Antibiotic_type", "ARDB_Class", "ARDB_Cl_description"]
        if self.new_card:
            data_col["card"] = ["CARD_ARO", "CARD_ARO_name", "CARD_ARO_description", "CARD_AMR_Gene_Family", "CARD_Drug_Class", "CARD_Antibiotic_class", "CARD_Resistance_Mechanism"]
        else:
            data_col["card"] = ["CARD_ARO", "CARD_ARO_name", "CARD_ARO_description", "CARD_ARO_category", "CARD_Class"]
        data_col["vfdb"] = ["VFDB_VFs", "VFDB_VF_Function", "VFDB_Species", "VFDB_Level1", "VFDB_Level2"]
        nr = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        nr2 = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        for each in databases:
            each = each.lower()
            if each in ["cog", "kegg", "cazy", "ardb", "card", "vfdb"]:
                new_list = data_col[each]
                new_columns1.extend(new_list)
            else:
                nr2.remove(each.capitalize())
        new_columns2 = [item for item in nr if item not in nr2]
        new_columns = new_columns + new_columns2 + new_columns1
        return new_columns

    def run_select_kegg(self, select, kegg, outDir):
        select = eval(select)
        pathway_str = select["pathway_id"]
        pathway_ids = pathway_str.split(",")
        reader5 = pd.read_table(kegg, sep='\t', header=0, iterator=True)    #modify by zhangqingchen 20180831
        annotable = table_parser(reader5, chunk_size=100000, ignore_index=True).drop_duplicates('#Query')     #modify by zhangqingchen 20180831
        kegg = pd.read_csv(kegg, sep='\t', header=0, chunksize=100000)
        has_header = True
        if os.path.exists(outDir + "/anno_kegg.xls"):
            os.path.remove(outDir + "/anno_kegg.xls")
            os.path.remove(outDir + "/gene_list.xls")
        for annotable in kegg:
            try:
                select = annotable[
                    annotable['Pathway'].apply((lambda x: len(set(x.split(";")).intersection(set(pathway_ids))) > 0))]
                select_genes = pd.DataFrame(select["#Query"])
                select_genes.columns = ["GeneID"]
                select_genes.to_csv(outDir + "/gene_list.xls", sep="\t", index=False, header=has_header, mode='a+')
                select.to_csv(outDir + "/anno_kegg.xls", sep="\t", index=False, header=has_header, mode='a+')
                has_header = False
            except Exception as e:
                raise Exception("kegg创建基因丰度表失败——{}".format(e))

    def get_id(self, id):
        if not isinstance(id, ObjectId):
            if isinstance(id, types.StringTypes):
                id = ObjectId(id)
        return id


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='[select]', required=True, help='Input seletc condition')
    parser.add_argument('-type', metavar='[reset type]', required=True, help='1 for overview, 2 for kegg pathway')
    parser.add_argument('-o', metavar='[output Directory]', required=True, help='output Directory name')
    parser.add_argument('-kegg', metavar='[gene_kegg_anno]', help='Input kegg gene_anno_file')
    parser.add_argument('-database', metavar='[database_list]', help='Input database_list result to output')
    parser.add_argument('-of', metavar='[overviewfile]', help='Get local overviewfile', default="")
    parser.add_argument('-id', metavar='[task_id]', help='Get task id for nr detail',default="")
    parser.add_argument('-st', metavar='[select type]', help='Get select type',default="2") # 1: or ; 2:and
    parser.add_argument('-mongo', help='version of mongo client', type=int, default=None)
    args = parser.parse_args()
    select = args.s
    mytype = args.type
    outDir = args.o
    of = args.of
    select_type = str(args.st)
    run_overview = Overview(args.mongo)
    task_id = args.id
    if mytype == 1 or mytype == "1":
        database = args.database
        run_overview.run_select_overview(select, database, of, outDir)
    elif mytype == 2 or mytype == "2":
        kegg = args.kegg
        run_overview.run_select_kegg(select, kegg, outDir)
    elif mytype == 3 or mytype == "3":
        #database = args.database
        #run_overview.run_select_overview(select, database, of, outDir, mytype=3)
        run_overview.run_select_from_overviewfile(of, select, outDir, task_id, select_type)
    else:
        raise Exception("type类型必须为1或2")
