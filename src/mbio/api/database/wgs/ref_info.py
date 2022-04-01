# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'

from api_base import ApiBase
from bson.objectid import ObjectId
import json
import argparse


class RefInfo(ApiBase):  # 父类
    def __init__(self, bind_object):
        """
        __author__ = HONGDONG
        __last_modify__ = 20181112
        :param bind_object:
        """
        super(RefInfo, self).__init__(bind_object)
        self._project_type = "dna_wgs"
        self.json_dict = ""

    def add_sg_genome(self, json_path):
        """
        添加基因组数据，独立与整个流程的, 这里没有使用字典的循环，是为了确保json文件中下面字段全部有
        :return:
        """
        data_list = []
        self.json_dict = self.get_json(json_path)
        for keys in self.json_dict.keys():
            result = self.col_find_one("sg_species", {"species": keys})
            if not result:
                insert = {
                    "species": keys
                }
                main_id = self.db["sg_species"].insert_one(insert).inserted_id  # inserted_id &&
                self.update_db_record("sg_species", {"_id": main_id}, {"main_id": main_id})
            else:
                main_id = result["main_id"]  # &&
            for keys_ in self.json_dict[keys]:
                # if self.col_find_one("sg_genome", {"species": keys, "genome_version": keys_}):
                if self.col_find_one("sg_species_version", {"species": keys, "genome_version": keys_}):
                    continue
                else:
                    insert_data = {
                        "species_id": main_id,
                        "species": keys,
                        "genome_version": keys_,
                        "strain_characteristic": self.json_dict[keys][keys_]['strain_characteristic'],
                        "edition": self.json_dict[keys][keys_]['edition'],
                        "release_date": self.json_dict[keys][keys_]['release_date'],
                        "genome_size": self.json_dict[keys][keys_]['genome_size'],
                        "n50": self.json_dict[keys][keys_]['n50'],
                        "assembly": self.json_dict[keys][keys_]['assembly'],
                        "link": self.json_dict[keys][keys_]['link'],
                        "ref": self.json_dict[keys][keys_]['ref'],
                        "gff": self.json_dict[keys][keys_]['gff'],
                        "anno": self.json_dict[keys][keys_]['anno'],
                        "change_log": self.json_dict[keys][keys_]['change_log'],
                        "info_log": self.json_dict[keys][keys_]['info_log'],
                        "ref_fa_amb": self.json_dict[keys][keys_]['ref_fa_amb'],
                        "ref_fa_ann": self.json_dict[keys][keys_]['ref_fa_ann'],
                        "ref_fa_bwt": self.json_dict[keys][keys_]['ref_fa_bwt'],
                        "ref_fa_fai": self.json_dict[keys][keys_]['ref_fa_fai'],
                        "ref_fa_pac": self.json_dict[keys][keys_]['ref_fa_pac'],
                        "ref_fa_sa": self.json_dict[keys][keys_]['ref_fa_sa'],
                        "ref_dict": self.json_dict[keys][keys_]['ref_dict'],
                        "ref_chrlist": self.json_dict[keys][keys_]['ref_chrlist'],
                        "total_chrlist": self.json_dict[keys][keys_]['total_chrlist'],
                        "snpeff_path": self.json_dict[keys][keys_]['snpeff_path'],
                        "ssr_path": self.json_dict[keys][keys_]['ssr_path'],
                        "web_dir": self.json_dict[keys][keys_]['web_dir']
                    }
                    data_list.append(insert_data)
        if data_list:
            # self.col_insert_data("sg_genome", data_list)
            self.col_insert_data("sg_species_version", data_list)
        else:
            print "参考组信息库中已经有了不再进行插入！"

    def get_json(self, json_path):
        """
        用于解析出json文件，然后可以从中获取参考基因组数据
        :return:
        """
        f = open(json_path, "r")
        json_dict = json.loads(f.read())
        return json_dict

    def update_secret(self, genome_version_id):
        """
        用于区分是不是客户的不愿公开的基因组文件
        :param genome_version_id:
        :return:
        """
        self.update_db_record("sg_species_version", {"_id": genome_version_id}, {"is_secret": "true"})

    def find_genome_id(self, species, genome_version):
        return self.db['sg_species_version'].find_one({"species": species, "genome_version": genome_version})['_id']


if __name__ == "__main__":   # main 自己
    parser = argparse.ArgumentParser(description="手动导入物种参考基因组")
    parser.add_argument("-j", "--json_file", type=str, help="json文件")
    parser.add_argument("-p", "--project_type", type=str, help="project_type", default='dna_wgs')
    args = parser.parse_args()
    t = RefInfo(None)
    if args.project_type != 'dna_wgs':
        t._project_type = args.project_type
    t.add_sg_genome(args.json_file)
