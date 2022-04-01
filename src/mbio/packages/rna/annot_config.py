# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

"""获取注释相关文件目录，文件"""

import ConfigParser
import random
import os
import re
import sys
import json
from biocluster.config import Config
import glob
from biocluster.api.database.base import Base, report_check
import types
from bson.objectid import ObjectId
from collections import OrderedDict
import warnings

class Mydict2(dict):
    def __init__(self, input_data=None):
        if input_data:
            super(Mydict2, self).__init__(input_data)
        else:
            super(Mydict2, self).__init__({})

    def __getitem__(self, key):
        if key in super(Mydict2, self).keys():
            return super(Mydict2, self).__getitem__(key)
        else:
            return None

class Mydict(dict):
    def __init__(self, input_data=None):
        if input_data:
            super(Mydict, self).__init__(input_data)
        else:
            super(Mydict, self).__init__({})

    def __getitem__(self, key):
        if key in super(Mydict, self).keys():
            return super(Mydict, self).__getitem__(key)
        else:
            return Mydict2({})

class ApiBase(Base):
    def __init__(self, bind_object, project_type="ref_rna_v2"):
        super(ApiBase, self).__init__(bind_object)
        self._project_type = project_type

    @property
    def db(self):
        if self._db is None:
            self._db = self._ref_client[self._config.get_mongo_dbname(self._project_type, True)]
        return self._db

    @property
    def db_p(self):
        self._db_p = Config().get_mongo_client(mtype=self._project_type, db_version=1)[Config().get_mongo_dbname(self._project_type, db_version=1)]
        return self._db_p

    def create_db_table(self, table_name, content_dict_list, tag_dict=None):
        '''
        Create main/detail table in database system.
        :param table_name: table name
        :param content_dict_list: list with dict as elements
        :param tag_dict: a dict to be added into each record in content_dict_list.
        :return: None or main table id
        '''
        table_id = None
        conn = self.db[table_name]
        if tag_dict:
            for row_dict in content_dict_list:
                row_dict.update(tag_dict)
        record_num = len(content_dict_list)
        try:
            if record_num > 5000:
                for i in range(0, record_num, 3000):
                    tmp_list = content_dict_list[i: i + 3000]
                    conn.insert_many(tmp_list)
            else:
                if record_num >= 2:
                    conn.insert_many(content_dict_list)
                else:
                    table_id = conn.insert_one(content_dict_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                warnings.warn('fail to insert records into table {} -> ({})'.format(table_name, e))
            else:
                warnings.warn('fail to insert record into table {} -> ({})'.format(table_name, e))

        return table_id

    def create_db_table2(self, table_name, content_dict_list, tag_dict=None):
        '''
        Create main/detail table in database system.
        :param table_name: table name
        :param content_dict_list: list with dict as elements
        :param tag_dict: a dict to be added into each record in content_dict_list.
        :return: None or main table id
        '''
        table_id = None
        conn = self.db_p[table_name]
        if tag_dict:
            for row_dict in content_dict_list:
                row_dict.update(tag_dict)
        record_num = len(content_dict_list)
        try:
            if record_num > 5000:
                for i in range(0, record_num, 3000):
                    tmp_list = content_dict_list[i: i + 3000]
                    conn.insert_many(tmp_list)
            else:
                if record_num >= 2:
                    conn.insert_many(content_dict_list)
                else:
                    table_id = conn.insert_one(content_dict_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                warnings.warn('fail to insert records into table {} -> ({})'.format(table_name, e))
            else:
                warnings.warn('fail to insert record into table {} -> ({})'.format(table_name, e))

        return table_id

    def find_db_and_update(self, table_name, update_func):
        conn = self.db[table_name]
        results = conn.find({})
        for result in results:
            conn.update({'_id': result['_id']}, {'$set': update_func(result)}, upsert=True)

class AnnotConfig(object):
    def __init__(self, config_path=None, group_config_path=None):
        """
        config_path: 配置文件 默认脚本当前目录 annot.conf
        group_config_path: 配置文件 默认脚本当前目录 annot_group.conf
        """
        self.main_conf = Config()
        self.annot_config = self.import_annot_config(config_path)
        self.group_config = self.import_group_config(group_config_path)

    def import_annot_config(self, config_path=None):
        if not config_path:
            config_path =  os.path.dirname(os.path.realpath(__file__)) + '/annot.conf'
        annot_config = ConfigParser.RawConfigParser()
        annot_config.read(config_path)
        return annot_config

    def import_group_config(self, config_path=None):
        if not config_path:
            config_path =  os.path.dirname(os.path.realpath(__file__)) + '/annot_group.conf'
        group_config = ConfigParser.RawConfigParser()
        group_config.read(config_path)
        return group_config

    def list_group(self):
        """
        显示有哪些 group
        """
        return self.group_config.sections()

    def list_group_option(self, section=None):
        """
        显示指定section的 配置
        section: 指定section
        """
        if section:
            return self.group_config.options(section)
        else:
            print "section needed"

    def get_group_option_detail(self, section=None):
        """
        显示指定section的 指向的数据库信息
        section: 指定section
        """
        # 获取该group所有的属性信息
        goptions = self.list_group_option(section)
        gsection_dict = Mydict()
        for goption in goptions:
            if goption in ['version']:
                continue
            gvalue = self.group_config.get(section, goption)
            value2 = Mydict2({})
            try:
                value2 = Mydict2(self.get_option_by_section(section=gvalue))
            except Exception as e:
                print "{}, 检查配置文件".format(e)
            value2.update({"section": gvalue})
            gsection_dict[goption] = value2

        return gsection_dict


    def list_sections(self):
        """
        显示有哪些 数据库
        """
        return self.annot_config.sections()

    def list_options(self, section=None):
        """
        显示指定section的 指向的数据库信息
        section: 指定section
        """
        if section:
            return self.annot_config.options(section)
        else:
            print "section needed"

    def get_option_by_section(self, section=None, option=None):
        """
        指定section的 的信息
        section: 指定section
        option: 指定则返回指定的值， 否则返回所有值的字典
        """
        if option:
            value = self.annot_config.get(section, option)
            return value
        else:
            options_value = dict()
            for option in self.annot_config.options(section):
                options_value.update({option: self.annot_config.get(section, option)})
            return options_value

    def get_dmnd_path(self, db_name=None, version=None, diamond_version="v0.9.24.125", **kwargs):
        """
        获取diamond 索引路径
        db_name: 数据库名称如nr, kegg, Animal
        version: diamond 合并的版本号
        **kwargs: 其它指定各自数据库的版本号，如 nr_version="202006" 如果找到相关数据库会返回该指定的索引路径
        """
        db_path = None
        try:
            db_path = self.get_file_path(
                file=db_name,
                version=version,
                soft="diamond",
                db="diamond",
                soft_version=diamond_version)
        except:
            print "not found dmnd path"

        for k,v in kwargs.items():
            if v:
                try:
                    db_path = AnnotConfig().get_file_path(
                        file=db_name,
                        db=k.split("_")[0],
                        version=v,
                        soft="diamond",
                        soft_version=diamond_version)
                except:
                    pass
        if db_path:
            print "db is {}".format(db_path)
        else:
            raise Exception("not found diamand db")
        return db_path


    def get_option_by_option(self, **kwargs):
        '''
        根据参数选择数据库，如果匹配多条选择第一个
        '''
        section_list = list()
        for section in self.annot_config.sections():
            option_dict = self.get_option_by_section(section)
            option_select = {k: v for k,v in kwargs.items() if k not in ['section', 'file']}
            # print option_dict, kwargs
            if set(option_select.items()).issubset(set(option_dict.items())):
                print kwargs.items()
                print option_dict
                section_list.append((section, option_dict))
        if len(section_list) == 0:
            raise Exception("error not found database")
        else:
            if len(section_list) > 1:
                print "warning found multi database use first one, need check configure"
                return section_list[0][1]
            else:
                return section_list[0][1]

    def abs_path(self, path=None):
        return os.path.join(self.main_conf.SOFTWARE_DIR, path)

    def get_path(self, section=None, **kwargs):
        '''
        根据section或者其它参数获取数据库路径，如果结果有多条返回第一条
        '''
        if section:
            path = self.get_option_by_section(section, option="dir")
        else:
            option_dict = self.get_option_by_option(**kwargs)
            path =  option_dict['dir']

        path = self.abs_path(path)
        return path

    def get_file_dict(self, section=None, **kwargs):
        '''
        根据section或者其它参数获取数据库路径，获取相关数据库的文件，返回所有文件字典
        '''
        if section:
            option_dict = self.get_option_by_section(section)
        else:
            option_dict = self.get_option_by_option(**kwargs)

        path =  option_dict['dir']
        path = self.abs_path(path)

        files_path = option_dict['files_path']
        files = option_dict['files'].split(",")

        file_dict = dict()


        print "find database file is {} db attrbute is {}".format(path, option_dict)
        for file in files:
            file_dict.update({file: files_path.format(dir=path, file=file)})
        if "file_alias" in option_dict:
            alias_dict = {ali.split(":")[0]:ali.split(":")[1] for ali in option_dict['file_alias'].strip(",").split(",")}
            for file_alias in alias_dict:
                file_dict.update({file_alias: files_path.format(dir=path, file=alias_dict[file_alias])})
        print "file_dict is {}".format(file_dict)
        return file_dict

    def get_file_path(self, section=None, file=None, **kwargs):
        '''
        根据section或者其它参数获取数据库路径，获取相关指定的文件
        '''
        if section:
            option_dict = self.get_option_by_section(section)
        else:
            option_dict = self.get_option_by_option(**kwargs)

        print option_dict
        path = option_dict['dir']
        path = self.abs_path(path)

        files_path = option_dict['files_path']
        files = option_dict['files'].split(",")

        file_dict = dict()
        if file:
            if file in files:
                return files_path.format(dir=path, file=file)
            else:
                if "file_alias" in option_dict:
                    alias_dict = {ali.split(":")[0]:ali.split(":")[1] for ali in option_dict['file_alias'].split(",")}
                    if file in alias_dict:
                        return files_path.format(dir=path, file=alias_dict[file])
                    else:
                        raise Exception("error not found database {}".format(file))
                else:
                    raise Exception("error not found database {}".format(file))
        else:
            raise Exception("please using get_file_dict".format(file))


    def get_kegg_genome2dict(self, genome_file):
        name2id = dict()
        with open(genome_file, 'r') as f:
            for line in f:
                genome_id = line.split("\t")[0]
                abr = line.split("\t")[1].split(",")[0]
                species = line.split("\t")[1].split(";")[-1].strip()
                name2id[abr] = (genome_id, abr)
                name2id[species] = (genome_id, abr)
        return name2id

    def get_gene2path(self, gene_path):
        '''
        获取gene pathway对应关系
        '''
        print gene_path
        species_abr = os.path.basename(gene_path)
        path_confs = glob.glob(gene_path + "/*.conf")
        gene2path = dict()
        # print path_confs
        for path_conf in path_confs:
            pathway = os.path.basename(path_conf).split(".")[0]
            with open(path_conf, 'r') as f:
                for line in f:
                    # print line
                    cols = line.strip("\n").split("\t")
                    if len(cols) >= 3:
                        genes = [g_str.split(" ")[0] for g_str in cols[2].split(",")]
                        for gene in genes:
                            gene = species_abr + ":" + gene
                            if gene in gene2path:
                                gene2path[gene].add(pathway)
                            else:
                                gene2path[gene] = set([pathway])

        return gene2path

    def get_gene2path_new(self, gene_path_new=None, abr=None):
        '''
        获取gene pathway对应关系 2017版后
        '''
        gene2path = dict()
        with open(gene_path_new, 'r') as f:
            f.readline()
            for line in f:
                paths = line.strip().split("\t")[-1].replace("ko", abr).split(";")
                gene2path[line.strip().split("\t")[0]] = set(paths)
        return gene2path

    def get_kegg_gene2path(self, genome_abr=None, version=None):
        '''
        根据section单物种fasta 描述关系字典

        '''
        if version >= "202007":
            kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=version)
            gene_path = os.path.join(kegg_files_dict['data'], genome_abr + '_allGene_ko.xls')
            gene2path_dict = self.get_gene2path_new(gene_path, genome_abr)

        else:
            kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=version)
            gene_path = kegg_files_dict["species_path"] + "/{}".format(genome_abr)
            gene2path_dict = self.get_gene2path(gene_path)
        return gene2path_dict

    def get_kegg_ko2path(self, version=None, remove_global=True):
        '''
        获取ko pathway 对应关系
        '''
        kegg_pathway = self.get_file_path(file="pathway", db="kegg", db_type="file", version=version)
        ko2path = dict()
        global_path = ["map01100", "map01110", "map01120", "map01130", "map01200", "map01210", "map01212", "map01230", "map01220"]
        with open(kegg_pathway, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                path = cols[0].split(":")[1]
                ko = cols[1].split(":")[1]
                if path in global_path and remove_global:
                    continue
                if ko in ko2path:
                    if path in ko2path[ko]:
                        pass
                    else:
                        ko2path[ko].add(path)
                else:
                    ko2path[ko] = set([path])
        self.ko2path = ko2path
        return ko2path

    def get_kegg_gene2ko(self, version):
        '''
        获取gene ko对应关系
        '''
        gene2ko_file = self.get_file_path(file="ko_genes.list", db="kegg", db_type="file", version=version)
        gene2ko = dict()
        with open(gene2ko_file, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                ko = cols[0].split(":")[1]
                gene = cols[1]
                if gene in gene2ko:
                    gene2ko[gene].add(ko)
                else:
                    gene2ko[gene] = set([ko])
        gene2ko = gene2ko
        return gene2ko

    def get_kegg_ko2name(self, version=None):
        '''
        获取ko 描述对应关系
        '''
        if not version:
            raise Exception("kegg_version needed")
        else:
            kegg_des_file = self.get_file_path(file="ko_des", db="kegg", db_type="file", version=version)
            with open(kegg_des_file, 'r') as kegg_des_f:
                ko2des = [(line.strip().split("\t")[0][3:], line.strip().split("\t")[1].split(";")[0].strip())
                          for line in kegg_des_f]
            return dict(ko2des)

    def get_kegg_fa(self, version=None, species=None):
        '''
        获取kegg 单物种fasta 描述关系字典
        '''

        if version >= "202007":
            kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=version)
            name_dict = self.get_kegg_genome2dict(kegg_files_dict['genome2.xls'])
            (genome_id, genome_abr) = name_dict[species]
            fa_path = os.path.join(kegg_files_dict['data'], genome_abr + '.faa')

        else:
            kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version="2017")
            name_dict = self.get_kegg_genome2dict(kegg_files_dict['genome2.xls'])
            (genome_id, genome_abr) = name_dict[species]
            fa_path = os.path.join(kegg_files_dict['species_gene'], genome_abr, genome_id + '.pep')
        return fa_path

    def get_kegg_spe_path_dir(self, version=None, species=None):
        '''
        获取kegg 单物种fasta 描述关系字典
        '''

        if version >= "202007":
            kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=version)
            name_dict = self.get_kegg_genome2dict(kegg_files_dict['genome2.xls'])
            (genome_id, genome_abr) = name_dict[species]
            png_path = os.path.join(kegg_files_dict['html_database'], 'png_down_finished')

        else:
            kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version="2017")
            name_dict = self.get_kegg_genome2dict(kegg_files_dict['genome2.xls'])
            (genome_id, genome_abr) = name_dict[species]
            png_path = os.path.join(kegg_files_dict['species_pathway'], genome_abr)
        return genome_abr,png_path

    def get_kegg_ko2des(self, version=None):
        '''
        获取kegg ko 描述关系字典
        '''
        if not version:
            raise Exception("kegg_version needed")
        else:
            kegg_des_file = self.get_file_path(file="ko_des", db="kegg", db_type="file", version=version)
            with open(kegg_des_file, 'r') as kegg_des_f:
                ko2des = [(line.strip().split("\t")[0][3:], line.strip().split("\t")[1].split(";")[-1].strip())
                          for line in kegg_des_f.readlines()]
        return dict(ko2des)

    def get_keggdb_paths(self, version=None):
        kegg_json = self.get_file_path(file="br08901.json", db="kegg", db_type="file", version=version)
        with open(kegg_json, "rb") as f:
            root = json.load(f)
        classI = root['children']
        classII = []
        for i in classI:
            classII.extend(i['children'])
        classIII = []
        for i in classII:
            classIII.extend(i['children'])
        db_paths = ["map" + str(i['name']).split(" ")[0] for i in classIII]
        return db_paths

    def get_kegg_map2class(self, version=None):
        '''
        获取kegg 通路与分类字典
        '''
        if not version:
            raise Exception("kegg_version needed")
        map2class = dict()
        kegg_json = self.get_file_path(file="br08901.json", db="kegg", db_type="file", version=version)
        with open(kegg_json, "rb") as f:
            root = json.load(f)
        classI = root['children']
        for class1 in classI:
            class1_name = class1['name']
            for class2 in class1['children']:
                class2_name = class2['name']
                paths = class2['children']
                for path in paths:
                    path_name = "map" + str(path['name']).split(" ")[0]
                    des = str(path['name']).split("  ")[1]
                    map2class[path_name] = (class1_name, class2_name, des)
        return map2class

    def update_db_group(self, cmd_command=None, annot_group=None, type=None):
        '''
        更新mongo库
        '''

        options = self.get_group_option_detail(section=annot_group)
        params = {k: v['version'] for k, v in options.items() if 'version' in v}
        group_dict = {
            "cmd_command": cmd_command,
            "annot_group": annot_group,
            "type": type,
            "parameter": params
        }
        api = ApiBase(None)

        print group_dict
        api.create_db_table("annot_group", [group_dict])


    def update_database_software(self, group_list=None):
        '''
        更新mongo库
        '''

        default_list = [
            {
                "software_database" : "KEGG 数据库",
                "source" : "http://www.genome.jp/kegg/",
                "version" : "Version 2020.03",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "eggNOG 数据库",
                "source" : "http://eggnogdb.embl.de/#/app/home",
                "version" : "Version 5.0",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "Pfam 数据库",
                "source" : "http://pfam.xfam.org/",
                "version" : "Version v32.0",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "Swiss-prot 数据库",
                "source" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                "version" : "Version 2019.7.1",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "GO 数据库",
                "source" : "http://www.geneontology.org/",
                "version" : "Version 2019.7.1",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "NR 数据库",
                "source" : "ftp://ftp.ncbi.nlm.nih.gov/blast/db/",
                "version" : "Version 2019.6.26",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "Rfam 数据库",
                "source" : "http://rfam.janelia.org/",
                "version" : "Version Rfam v14.1",
                "usage" : "核糖体比例评估"
            }
        ]
        db_des = {
            "rfam" : "Rfam 数据库",
            "pfam" : "Pfam 数据库",
            "kegg" : "KEGG 数据库",
            "eggnog" : "eggNOG 数据库",
            "swissprot" : "Swiss-prot 数据库",
            "uniprot" : "Uniprot 数据库",
            "pir" : "PIR idmapping 数据库",
            "ncbi_taxonomy" : "NCBI 物种分类数据库",
            "rfam" : "Rfam 数据库",
            "nr" : "NR 数据库",
            "go" : "GO 数据库",
        }

        db_source = {
            "rfam" : "http://rfam.janelia.org/",
            "pfam" : "http://pfam.xfam.org/",
            "kegg" : "http://www.genome.jp/kegg/",
            "eggnog" : "http://eggnogdb.embl.de/#/app/home",
            "swissprot" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            "uniprot" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release",
            "pir" : "ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz",
            "ncbi_taxonomy" : "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
            "rfam" : "http://rfam.janelia.org/",
            "nr" : "https://www.ncbi.nlm.nih.gov/public/",
            "go" : "http://www.geneontology.org/",
        }

        db_useage = {
            "rfam" : "http://rfam.janelia.org/",
            "pfam" : "基因注释",
            "kegg" : "基因注释",
            "eggnog" : "基因注释",
            "swissprot" : "基因注释",
            "uniprot" : "基因注释",
            "pir" : "基因注释",
            "ncbi_taxonomy" : "基因注释",
            "rfam" : "核糖体比例评估",
            "nr" : "基因注释",
            "go" : "基因注释",
        }


        def convert_func(rec_dict):
            annot_group = rec_dict["annot_group"]
            options = self.get_group_option_detail(section=annot_group)
            db_show = list()

            print options
            for k, v in options.items():
                version = v["version"]
                if k == "string":
                    version = "11.0"

                if k == "cog":
                    version = "20"

                # print("k {} v {}".format(k, v))
                if version.startswith("20"):
                    version = version[0:4] + "." + version[4:]

                if k in db_des:
                    db_show.append({
                        "software_database" : db_des[k],
                        "source": db_source[k],
                        "useage": db_useage[k],
                        "version": "Version {}".format(version)
                    })

            return {"software_database": db_show}

        groups = []
        with open(group_list, 'r') as f:
            for line in f:
                [annot_group, cmd_command, type] = line.strip().split("\t")
                options = self.get_group_option_detail(section=annot_group)
                params = {k: v['version'] for k, v in options.items() if 'version' in v}
                group_dict = {
                    "cmd_command": cmd_command,
                    "annot_group": annot_group,
                    "collection" : "sg_task",
                    "condition": {"annot_group": annot_group}
                }
                group_dict.update(convert_func(group_dict))

                groups.append(group_dict)

                if type == "old":
                    group_dict_default = {
                        "cmd_command": cmd_command,
                        "annot_group": annot_group,
                        "collection" : "sg_task",
                        "condition" : {"annot_group": {"exists": False}}
                    }
                    group_dict_default.update(convert_func(group_dict_default))

                    groups.append(group_dict_default)
        # group_dict_default = {
        #     "cmd_command": cmd_command,
        #     "collection":  "sg_task",
        #     "condition" : {"annot_group": {"$exists": False}},
        #     "software_database_list" : default_list
        # }
        api = ApiBase(None)
        api.create_db_table("annot_database_software", groups)

    def update_database_software2(self, group_list=None,  project_type="ref_rna_v2"):
        '''
        更新mongo库
        '''

        default_list = [
            {
                "software_database" : "KEGG 数据库",
                "source" : "http://www.genome.jp/kegg/",
                "version" : "Version 2020.03",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "eggNOG 数据库",
                "source" : "http://eggnogdb.embl.de/#/app/home",
                "version" : "Version 5.0",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "Pfam 数据库",
                "source" : "http://pfam.xfam.org/",
                "version" : "Version v32.0",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "Swiss-prot 数据库",
                "source" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                "version" : "Version 2019.7.1",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "GO 数据库",
                "source" : "http://www.geneontology.org/",
                "version" : "Version 2019.7.1",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "NR 数据库",
                "source" : "ftp://ftp.ncbi.nlm.nih.gov/blast/db/",
                "version" : "Version 2019.6.26",
                "usage" : "参考基因组注释"
            },
            {
                "software_database" : "Rfam 数据库",
                "source" : "http://rfam.janelia.org/",
                "version" : "Version Rfam v14.1",
                "usage" : "核糖体比例评估"
            }
        ]
        db_des = {
            "rfam" : "Rfam 数据库",
            "pfam" : "Pfam 数据库",
            "kegg" : "KEGG 数据库",
            "eggnog" : "eggNOG 数据库",
            "swissprot" : "Swiss-prot 数据库",
            "uniprot" : "Uniprot 数据库",
            "pir" : "PIR idmapping 数据库",
            "ncbi_taxonomy" : "NCBI 物种分类数据库",
            "rfam" : "Rfam 数据库",
            "nr" : "NR 数据库",
            "go" : "GO 数据库",
        }

        db_source = {
            "rfam" : "http://rfam.janelia.org/",
            "pfam" : "http://pfam.xfam.org/",
            "kegg" : "http://www.genome.jp/kegg/",
            "eggnog" : "http://eggnogdb.embl.de/#/app/home",
            "swissprot" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            "uniprot" : "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release",
            "pir" : "ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz",
            "ncbi_taxonomy" : "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
            "rfam" : "http://rfam.janelia.org/",
            "nr" : "https://www.ncbi.nlm.nih.gov/public/",
            "go" : "http://www.geneontology.org/",
        }

        db_useage = {
            "rfam" : "http://rfam.janelia.org/",
            "pfam" : "基因注释",
            "kegg" : "基因注释",
            "eggnog" : "基因注释",
            "swissprot" : "基因注释",
            "uniprot" : "基因注释",
            "pir" : "基因注释",
            "ncbi_taxonomy" : "基因注释",
            "rfam" : "核糖体比例评估",
            "nr" : "基因注释",
            "go" : "基因注释",
        }


        def convert_func(rec_dict, db_annot_group=None):
            annot_group = rec_dict["annot_group"]
            options = self.get_group_option_detail(section=annot_group)
            db_show = list()

            print options
            for k, v in options.items():
                version = v["version"]
                if k == "string":
                    version = "11.0"

                if k == "cog":
                    version = "20"

                # print("k {} v {}".format(k, v))
                if version.startswith("20"):
                    version = version[0:4] + "." + version[4:]

                if k in db_des:
                    db_show.append({
                        "software_database" : db_des[k],
                        "source": db_source[k],
                        "usage": db_useage[k],
                        "annot_group": db_annot_group,
                        "version": "Version {}".format(version)
                    })

            return db_show

        recs = []
        with open(group_list, 'r') as f:
            for line in f:
                [annot_group, cmd_command, type] = line.strip().split("\t")
                if cmd_command.split(".")[0] == project_type:
                    options = self.get_group_option_detail(section=annot_group)
                    params = {k: v['version'] for k, v in options.items() if 'version' in v}
                    group_dict = {
                        "annot_group": annot_group,
                    }
                    rec = convert_func(group_dict, annot_group)

                    recs.extend(rec)

                    if type == "old":
                        group_dict_default = {
                            "annot_group": annot_group,
                        }
                        rec = convert_func(group_dict, "false")
                        recs.extend(rec)
        # group_dict_default = {
        #     "cmd_command": cmd_command,
        #     "collection":  "sg_task",
        #     "condition" : {"annot_group": {"$exists": False}},
        #     "software_database_list" : default_list
        # }
        api = ApiBase(None, project_type=project_type)
        api.create_db_table2("sg_software_database", recs)



    def update_db_group2(self, group_list=None):
        '''
        更新mongo库
        '''

        db_des = {
            "rfam" : "Rfam 数据库",
            "pfam" : "Pfam 数据库",
            "kegg" : "KEGG 数据库",
            "eggnog" : "eggNOG 数据库",
            "string" : "STRING 数据库",
            "swissprot" : "Swiss-prot 数据库",
            "uniprot" : "Uniprot 数据库",
            "pir" : "PIR idmapping 数据库",
            "ncbi_taxonomy" : "NCBI 物种分类数据库",
            "rfam" : "Rfam 数据库",
            "nr" : "NR 数据库",
            "go" : "GO 数据库",
        }


        def convert_func(rec_dict):
            annot_group = rec_dict["annot_group"]
            options = self.get_group_option_detail(section=annot_group)
            db_show = dict()

            for k, v in options.items():
                version = v["version"]
                if k == "string":
                    version = "11.0"

                if k == "cog":
                    version = "20"

                if version.startswith("2019") or version.startswith("2020") or version.startswith("2021"):
                    version = version[0:4] + "." + version[4:]

                if k in db_des:
                    db_show.update({db_des[k] : "Version {}".format(version)})

            return {"db_show": db_show}

        groups = []
        with open(group_list, 'r') as f:
            for line in f:
                [annot_group, cmd_command, type] = line.strip().split("\t")
                options = self.get_group_option_detail(section=annot_group)
                params = {k: v['version'] for k, v in options.items() if 'version' in v}
                group_dict = {
                    "cmd_command": cmd_command,
                    "annot_group": annot_group,
                    "type": type,
                    "parameter": params
                }
                group_dict.update(convert_func(group_dict))
                groups.append(group_dict)
        print groups
        api = ApiBase(None)
        api.create_db_table("annot_group", groups)

    def update_db_group_show_base(self, group_list=None):
        '''
        更新mongo库
        '''
        groups = []
        db_des = {
            "rfam" : "Rfam 数据库",
            "pfam" : "Pfam 数据库",
            "kegg" : "KEGG 数据库",
            "eggnog" : "eggNOG 数据库",
            "string" : "STRING 数据库",
            "swissprot" : "Swiss-prot 数据库",
            "uniprot" : "Uniprot 数据库",
            "pir" : "PIR idmapping 数据库",
            "ncbi_taxonomy" : "NCBI 物种分类数据库",
            "rfam" : "Rfam 数据库",
            "nr" : "NR 数据库",
            "go" : "GO 数据库",
        }


        def convert_func(rec_dict):
            annot_group = rec_dict["annot_group"]
            options = self.get_group_option_detail(section=annot_group)
            db_show = dict()

            for k, v in options.items():
                version = v["version"]
                if k == "string":
                    version = "11.0"

                if version.startswith("2019") or version.startswith("2020") or version.startswith("2021"):
                    version = version[0:4] + "." + version[4:]

                if k in db_des:
                    db_show.update({db_des[k] : "Version {}".format(version)})

            return {"db_show": db_show}

        api = ApiBase(None)
        api.find_db_and_update(table_name="annot_group", update_func=convert_func)


if __name__ == '__main__':
    cf = AnnotConfig()
    os.environ["current_mode"]="workflow"
    os.environ["NTM_PORT"]="7322"
    os.environ["WFM_PORT"]="7321"

    function_list = [fun for fun in dir(cf) if fun.startswith("get") or fun.startswith("list") or fun.startswith("update")]
    if len(sys.argv) == 1 or sys.argv[1] in ["-h", "-help", "--h", "--help"]:
        if len(sys.argv) >= 3:
            help(getattr(cf, sys.argv[2]))
        else:
            print function_list

    elif len(sys.argv) == 2:
        function = getattr(cf, sys.argv[1])
        print "AnnotConfig().{}()".format(sys.argv[1])
        print function()

    else:
        params_dict = dict()
        for par in sys.argv[2:]:
            params_dict.update({par.split("=")[0]: "=".join(par.split("=")[1:])})
        function = getattr(cf, sys.argv[1])
        print "AnnotConfig().{}({})".format(sys.argv[1], ",".join(sys.argv[2:]))
        print function(**params_dict)
