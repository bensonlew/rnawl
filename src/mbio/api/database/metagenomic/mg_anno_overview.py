# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20170930
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from biocluster.config import Config
import copy


class MgAnnoOverview(Base):
    def __init__(self, bind_object):
        super(MgAnnoOverview, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'
        self.db_path = os.path.join(Config().SOFTWARE_DIR, "database/metagenome")
        self.level_name_file = os.path.join(self.db_path, "database_id_to_file.xls")
        self.all_convert_names = {}
        self.nr_dict = {}

    @report_check
    def add_anno_overview(self, geneset_id, name=None, task_id=None):
        # 主表, 所有的函数名称以add开头，里面可以加需要导入数据库而表格里没有的信息作为参数
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_overview_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52802301")
        if not task_id:
            task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '',
            'created_ts': created_ts,
            'name': name if name else "Overview_Origin",
            'params': 'null',
            'status': 'end',
            'geneset_id': geneset_id,
            'main_id': ''
        }
        collection = self.db['anno_overview']
        # 将主表名称写在这里
        anno_overview_id = collection.insert_one(insert_data).inserted_id
        # 将导表数据通过insert_one函数导入数据库，将此条记录生成的_id作为返回值，给detail表做输入参数
        return anno_overview_id

    @report_check
    def add_anno_overview_detail(self, anno_overview_id, overview_file):
        if not isinstance(anno_overview_id, ObjectId):  # 检查传入的anno_overview_id是否符合ObjectId类型
            if isinstance(anno_overview_id, str):  # 如果是string类型，则转化为ObjectId
                anno_overview_id = ObjectId(anno_overview_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_overview_id必须为ObjectId对象或其对应的字符串！', code="52802302")
        id_convert = MetagenomicController()
        result = self.db['anno_overview'].find_one({'_id': anno_overview_id})
        main_collection = self.db['anno_overview']
        main_collection.update({'_id': anno_overview_id}, {'$set': {'main_id': anno_overview_id}})
        self.all_convert_names = self.convert_name_to_id()
        data_list = []
        if not result:
            self.bind_object.set_error('找不到anno_overview对应的主表id', code="52802303")
        else:
            task_id = result['task_id']
            self.task_id = task_id
            # samples_dic = id2name(task_id, type="task")
        if not os.path.exists(overview_file):
            self.bind_object.logger.error('overview_file所指定的路径{}不存在，请检查！'.format(overview_file))
            self.bind_object.set_error("overview_file所指定的路径不存在", code="52802304")
        with open(overview_file, "rb") as f:
            self.bind_object.logger.info("start overview_file")
            head = f.next().strip().split("\t")
            this_map = self.get_this_personal_map(head)  # guanqing.zou
            if "Domain" in head:
                nr_mongo = self.db['anno_nr'].find_one({"task_id": task_id, "name": "NR_Origin", "status": "end"})
                nr_id = nr_mongo["_id"]
                self.bind_object.logger.info("start nr")
                result = self.db["anno_nr_detail"].find({"nr_id": nr_id, "level_id": 8})
                print result.count()
                print nr_id
                nr_dict = {}
                if result:
                    for each in result:
                        level_names = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
                        level_ids = ["d_id", "k_id", "p_id", "c_id", "o_id", "f_id", "g_id", "s_id"]
                        for n in range(0, len(level_names)):
                            nr_name = each[level_names[n]]
                            nr_id = each[level_ids[n]]
                            if not self.nr_dict.has_key(nr_name):
                                self.nr_dict[nr_name] = nr_id
                self.bind_object.logger.info("nr level ok")


            new_insert = 0
            for log_index,line in enumerate(f):
                if log_index % 200000 == 0 and log_index > 0:
                    new_insert = 1
                    self.bind_object.logger.info("have done %s lines" % log_index)
                line = line.strip().split('\t')
                gene = line[0]
                length = float(line[1])
                insert_data = {
                    'anno_overview_id': anno_overview_id,
                    'query': gene,
                    'length': length,
                }
                # NR
                if "Domain" in head:
                    nr_d = line[head.index("Domain")]
                    nr_k = line[head.index("Kingdom")]
                    nr_p = line[head.index("Phylum")]
                    nr_c = line[head.index("Class")]
                    nr_o = line[head.index("Order")]
                    nr_f = line[head.index("Family")]
                    nr_g = line[head.index("Genus")]
                    nr_s = line[head.index("Species")]
                    insert_data['nr_1'] = self.nr_covert(nr_d)
                    insert_data['nr_2'] = self.nr_covert(nr_k)
                    insert_data['nr_3'] = self.nr_covert(nr_p)
                    insert_data['nr_4'] = self.nr_covert(nr_c)
                    insert_data['nr_5'] = self.nr_covert(nr_o)
                    insert_data['nr_6'] = self.nr_covert(nr_f)
                    insert_data['nr_7'] = self.nr_covert(nr_g)
                    insert_data['nr_8'] = self.nr_covert(nr_s)
                # COG
                if "NOG" in head:
                    cog_id = line[head.index("NOG")]
                    Category = line[head.index("COG_Category")]
                    Function = line[head.index("COG_Function")]
                    cog_category = self.each_covert(Category, 9)
                    cog_function = self.each_covert(Function, 10)
                    insert_data['cog_9'] = cog_category
                    insert_data['cog_10'] = cog_function
                    insert_data['cog_11'] = cog_id
                # KEGG
                if "KO" in head:
                    kegg_gene = line[head.index("KEGG_Gene")]
                    kegg_orthology = line[head.index("KO")]
                    kegg_pathway = line[head.index("KEGG_Pathway")]
                    kegg_enzyme = line[head.index("KEGG_Enzyme")]
                    kegg_module = line[head.index("KEGG_Modules")]
                    l1 = line[head.index("KEGG_Level1")]
                    l2 = line[head.index("KEGG_Level2")]
                    kegg_l1 = self.each_covert(l1, 12)
                    kegg_l2 = self.each_covert(l2, 13)
                    insert_data['kegg_gene'] = kegg_gene
                    insert_data['kegg_12'] = kegg_l1
                    insert_data['kegg_13'] = kegg_l2
                    insert_data['kegg_14'] = kegg_pathway
                    insert_data['kegg_15'] = kegg_module
                    insert_data['kegg_16'] = kegg_enzyme
                    insert_data['kegg_17'] = kegg_orthology
                # cazy
                if "CAZY_Family" in head:
                    family = line[head.index("CAZY_Family")]
                    if "CAZY_Cl_description" in head:
                        cazy_classes = line[head.index("CAZY_Cl_description")]
                    else:
                        cazy_classes = line[head.index("CAZY_Class")]
                    cazy_family = self.each_covert(family, 19)
                    cazy_class = self.each_covert(cazy_classes, 18)
                    insert_data['cazy_19'] = cazy_family
                    insert_data['cazy_18'] = cazy_class
                # ARDB
                if "ARDB_ARG" in head:
                    ardb_arg = line[head.index("ARDB_ARG")]
                    types = line[head.index("ARDB_Type")]
                    ardb_arg = self.each_covert(ardb_arg, 23)
                    ardb_type = self.each_covert(types, 21)
                    resistance = line[head.index("ARDB_Antibiotic_type")]
                    ardb_resistance = self.each_covert(resistance, 22)
                    ardb_classes = line[head.index("ARDB_Class")]
                    ardb_class = self.each_covert(ardb_classes, 20)
                    if "ARDB_Antibiotic_class" in head:##兼容个性化注释报错，因为个性化注释可能是老表
                        ardb_anti_classes = line[head.index("ARDB_Antibiotic_class")]
                        ardb_anti_class = self.each_covert(ardb_anti_classes, 73)
                        insert_data['ardb_73'] = ardb_anti_class
                    insert_data['ardb_20'] = ardb_class
                    insert_data['ardb_21'] = ardb_type
                    insert_data['ardb_22'] = ardb_resistance
                    insert_data['ardb_23'] = ardb_arg
                # CARD
                if "CARD_ARO" in head:
                    aro = line[head.index("CARD_ARO")]
                    if "CARD_ARO_name" in head:
                        card_aro_name = line[head.index("CARD_ARO_name")]
                        aro_name = self.each_covert(card_aro_name, 72)
                        insert_data['card_72'] = aro_name
                    # card_amr = line[head.index("CARD_AMR_Gene_Family")]
                    # card_drug = line[head.index("CARD_Drug_Class")]
                    if "CARD_Antibiotic_class" in head:
                        card_type = line[head.index("CARD_Antibiotic_class")]
                        type = self.each_covert(card_type, 71)
                        insert_data['card_71'] = type
                    if "CARD_Resistance_Mechanism" in head:
                        card_resistance = line[head.index("CARD_Resistance_Mechanism")]
                        resistance = self.each_covert(card_resistance, 70)
                        insert_data['card_70'] = resistance
                    # card_class = self.each_covert(card_amr, 24)
                    card_aro = self.each_covert(aro, 25)
                    # amr = self.each_covert(card_amr, 68)
                    # drug= self.each_covert(card_drug, 69)
                    insert_data['card_25'] = card_aro
                    # insert_data['card_24'] = card_class
                    # insert_data['card_68'] = amr ## fix by qingchen.zhang @202011 原因是level存在多个导入的数据是错误的，不再导由前端直接从参考库中调取
                    # insert_data['card_69'] = drug
                # VFDB
                if "VFDB_VFs" in head:
                    if "VFDB_VFs" in head:
                        vfs = line[head.index("VFDB_VFs")]
                        vfdb_vfs = self.each_covert(vfs, 28)
                        insert_data['vfdb_28'] = vfdb_vfs
                    VFDB_Species = line[head.index("VFDB_Species")]
                    if "VFDB_Level1" in head:
                        vf_l1 = line[head.index("VFDB_Level1")]
                        insert_data['vfdb_26'] = vf_l1
                    if "VFDB_Level2" in head:
                        vf_l2 = line[head.index("VFDB_Level2")]
                        insert_data['vfdb_27'] = vf_l2
                    # vfdb_l1 = self.each_covert(vf_l1, 26)## fix by qingchen.zhang @202011 原因是level存在多个导入的数据是错误的，不再转level的id
                    # vfdb_l2 = self.each_covert(vf_l2, 27)
                    insert_data['vfdb_species'] = VFDB_Species

                ## zouguanqing
                insert_data = self.insert_peronal_data(this_map, line,insert_data)  #guanqing.zou

                data_list.append(insert_data)
                if new_insert == 1 and log_index > 0:
                    self.bind_object.logger.info("开始导入数据库前%s lines" % log_index)
                    try:
                        collection = self.db['anno_overview_detail']
                        collection.insert_many(data_list,ordered=True)
                        new_insert = 0
                        data_list = []
                    except Exception,e:
                        self.bind_object.logger.error("导入前%slines失败:%s" % (log_index,e))
                        self.bind_object.set_error("导入数据失败", code="52802305")
            self.bind_object.logger.info("导入剩余数据")
            try:
                collection = self.db['anno_overview_detail']
                collection.insert_many(data_list,ordered=True)
                self.db['anno_overview'].update_one({'_id':anno_overview_id},{'$set':{'main_id':anno_overview_id}})
            except Exception,e:
                self.bind_object.logger.error("导入总览表失败,%s" % e)
                self.bind_object.set_error("导入总览表失败", code="52802306")
            self.bind_object.logger.info("总览表成功！")

    @report_check
    def convert_name_to_id(self):
        ## id = 9,10,12,13,18,19,20,21,22,23,24,25,26,27,28
        name_convert = {}
        with open(self.level_name_file, "rb") as input:
            input.next()
            for line in input:
                line = line.strip().split("\t")
                id = line[0]
                name = line[1]
                level_id = line[2]
                #print id,name,level_id
                if level_id == '10':
                    name = name[0]
                key_name = level_id + "_" + name
                if not name_convert.has_key(key_name):
                    name_convert[key_name] = id
        return name_convert

    @report_check
    def each_covert(self, levels_name, level_id):
        list = []
        if levels_name != "-":
            levels_list = levels_name.split(";")
            for each in levels_list:
                each = each.strip()  # zouguanqing
                if self.all_convert_names.has_key(str(level_id) + "_" + each):
                    levels = self.all_convert_names[str(level_id) + "_" + each]
                else:
                    levels = '-'
                    self.bind_object.logger.info('wwww: '+ str(level_id) + ':' + each)


                    #self.bind_object.set_error("%s id error", variables=(str(level_id) + "_" + each), code="52802307")
                list.append(levels)
            levels_name = ";".join(list)
        else:
            levels_name = "-"
        return levels_name

    @report_check
    def nr_covert(self, nr_name):
        if nr_name != "-":
            if self.nr_dict.has_key(nr_name):
                nr_id = self.nr_dict[nr_name]
            else:
                # self.bind_object.set_error("%s nr name error", variables=(nr_name), code="52802308")
                nr_id = '-'  # 不报错
        else:
            nr_id = "-"
        return nr_id


    #zouguanqing >>>
    def personal_map_all(self):
        self.personal_maps = {
            'lca_Domain':{'id':'lca', 'mongo_key':'lca_nr_1'},
            'lca_Kingdom':{'id':'lca', 'mongo_key':'lca_nr_2'},
            'lca_Phylum':{'id':'lca', 'mongo_key':'lca_nr_3'},
            'lca_Class':{'id':'lca', 'mongo_key':'lca_nr_4'},
            'lca_Order':{'id':'lca', 'mongo_key':'lca_nr_5'},
            'lca_Family':{'id':'lca', 'mongo_key':'lca_nr_6'},
            'lca_Genus':{'id':'lca', 'mongo_key':'lca_nr_7'},
            'lca_Species':{'id':'lca', 'mongo_key':'lca_nr_8'},
            'duc_Domain':{'id': 'duc', 'mongo_key':'duc_nr_1'},
            'duc_Kingdom':{'id': 'duc', 'mongo_key':'duc_nr_2'},
            'duc_Phylum':{'id': 'duc', 'mongo_key':'duc_nr_3'},
            'duc_Class':{'id': 'duc', 'mongo_key':'duc_nr_4'},
            'duc_Order':{'id': 'duc', 'mongo_key':'duc_nr_5'},
            'duc_Family':{'id': 'duc', 'mongo_key':'duc_nr_6'},
            'duc_Genus':{'id': 'duc', 'mongo_key':'duc_nr_7'},
            'duc_Species':{'id': 'duc', 'mongo_key':'duc_nr_8'},
            'pfam_Domain':{'id':'-', 'mongo_key':'pfam_29'},  #Domain
            'pfam_Type':{'id':30, 'mongo_key':'pfam_30'}, # Type
            'p450_homo_family':{'id':33, 'mongo_key':'p450_33'}, #Homologous Family
            'p450_super_family':{'id':34, 'mongo_key':'p450_34'},  #Superfamily
            'Mvirdb_Virulence_Factor_Type':{'id':37, 'mongo_key':'mvir_37'},  #Virulence Factor Type
            'Mvirdb_Short_Description':{'id':36, 'mongo_key':'mvir_36'},  # Gene Description（species）
            'phi_Pathogen_Species':{'id':40, 'mongo_key':'phi_40'}, # PHI Species
            'phi_Phenotype':{'id':41, 'mongo_key':'phi_41'}, # Phenotype
            'phi_Experimental_Host_Species':{'id':42, 'mongo_key':'phi_42'}, # Experimental_Host_Species
            #'phi_HOST':{'id':'-', 'mongo_key':'phi_42'},
            'TCDB_Family':{'id':50, 'mongo_key':'tcdb_50'}, # Family
            'TCDB_Subclass':{'id':49, 'mongo_key':'tcdb_49'}, # Subclass
            'TCDB_Class':{'id':48, 'mongo_key':'tcdb_48'}, #  Class
            'qs_Class':{'id':52, 'mongo_key':'qs_52'}, # Class
            'probio_Probiotic_name':{'id':53, 'mongo_key':'probio_53'}, #Probiotic name
            'probio_Genus':{'id':54, 'mongo_key':'probio_54'}, #Probiotics genus
            'probio_Use_in':{'id':55, 'mongo_key':'probio_55'}, #Use in
            'probio_Probiotic_Effect':{'id':66, 'mongo_key':'probio_66'}, #Probiotic   Effect
            'go_GO_Term_(Lev4)':{'id':'-', 'mongo_key':'go_62'}, #Level4
            'go_GO_Term_(Lev3)':{'id':'-', 'mongo_key':'go_61'}, #Level3
            'go_GO_Term_(Lev2)':{'id':'-', 'mongo_key':'go_60'}, #Level2
            'go_GO_(Lev1)':{'id':'-', 'mongo_key':'go_59'}, #Level1
            'sec_Gram_neg':{'id':'-', 'mongo_key':'secpro_63'}, #Gram neg ,True
            'sec_Gram_pos':{'id':'-', 'mongo_key':'secpro_64'}, #Gram pos ,True
            'sec_Euk':{'id':'-', 'mongo_key':'secpro_65'}, #Fungi  ,True
            't3ss_score':{'id':'-', 'mongo_key':'t3ss_67'}, # ,True

        }


    def get_this_personal_map(self,header):
        self.personal_map_all()
        this_map = copy.deepcopy(self.personal_maps)
        for h in self.personal_maps.keys():
            if h not in header:
                this_map.pop(h)
            else:
                this_map[h]['index'] = header.index(h)

        self.personal_nr_dic = {}
        if 'lca_Domain' in header:
            self.bind_object.logger.info("start get personal lca nr dic")
            self.personal_nr_dic['lca'] = self.get_personal_nr_dic(self.task_id, 'NR_Origin_LCA')
        if 'duc_Domain' in header:
            self.bind_object.logger.info("start get personal nr duc dic")
            self.personal_nr_dic['duc'] = self.get_personal_nr_dic(self.task_id, 'NR_Origin_Deunclassified')

        return this_map

    def insert_peronal_data(self,this_map, line, insert_data):
        for k in this_map.keys():
            id = this_map[k]['id']
            mongo_key = this_map[k]['mongo_key']
            line_index = this_map[k]['index']
            value = line[line_index].strip()
            if id == '-':
                new_value = value
            elif id == 'lca' or id == 'duc' :
                new_value = self.personal_nr_covert(self.personal_nr_dic[id], value)
                #self.bind_object.logger.info(new_value)
            else:
                new_value = self.each_covert(value, id)
            insert_data[mongo_key] = new_value
        return insert_data


    def get_personal_nr_dic(self, task_id, name):
        nr_mongo = self.db['anno_nr'].find_one({"task_id": task_id, "name": name, "status": "end"})
        nr_id = nr_mongo["_id"]
        self.bind_object.logger.info("start get personal nr dic")
        result = self.db["anno_nr_detail"].find({"nr_id": nr_id, "level_id": 8})

        nr_dict = {}
        if result:
            level_names = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
            level_ids = ["d_id", "k_id", "p_id", "c_id", "o_id", "f_id", "g_id", "s_id"]
            for each in result:
                for n in range(0, len(level_names)):
                    nr_name = each[level_names[n]]
                    nr_id = each[level_ids[n]]
                    if not nr_dict.has_key(nr_name):
                        nr_dict[nr_name] = nr_id
        return nr_dict

    def personal_nr_covert(self, nr_dict, nr_name):
        if nr_name != "-":
            if nr_dict.has_key(nr_name):
                nr_id = nr_dict[nr_name]
            else:
                self.bind_object.set_error("%s nr name error", variables=(nr_name), code="52802308")
        else:
            nr_id = "-"
        return nr_id

     #zouguanqing <<<<

    def personal_anno_dic(self): ##add by qingchen.zhang@20191209
        """
        没有个性化注释的结果
        :return:
        """
        self.level_correspond = {
            "query": "#Query",
            "length" : "Length",
            "nr_1": "Domain",
            "nr_2": "Kingdom",
            "nr_3": "Phylum",
            "nr_4": "Class",
            "nr_5": "Order",
            "nr_6": "Family",
            "nr_7": "Genus",
            "nr_8": "Species",
            "cog_9": "COG_Category",
            "cog_10": "COG_Function",
            "cog_11": "NOG",
            "kegg_12": "KEGG_Level1",
            "kegg_13": "KEGG_Level2",
            "kegg_14": "KEGG_Pathway",
            "kegg_15": "KEGG_Modules",
            "kegg_16": "KEGG_Enzyme",
            "kegg_17": "KO",
            "cazy_18": "CAZY_Class",
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
        }

    ###宏基因组兼容老项目用,add by qingchen.zhang@20191209
    def from_mongo_to_file(self, task_id):
        """
        从MongoDB数据库中将overview的detail表导出
        :param main_id:
        :return:
        """
        self.personal_anno_dic()
        anno_overview_path = os.path.join(self.bind_object.work_dir, "anno_overview.xls")
        collection = self.db['anno_overview_detail']
        results = collection.find_one({"anno_overview_id": task_id})
        with open(anno_overview_path, "w") as w:
            header_name = self.level_correspond.values()
            w.write("\t".join(header_name) + "\n")
            for result in results:
                result_dict = {}
                level_list = []
                for key in self.level_correspond.keys():
                    if result.get(key) != "-":
                        result_dict[self.level_correspond[key]] = result.get(key)
                    else:
                        result_dict[self.level_correspond[key]] = '-'
                for new_key in header_name:
                    level_list.append(new_key)
                w.write("{}\n".format("\t".join(level_list)))
        return anno_overview_path
