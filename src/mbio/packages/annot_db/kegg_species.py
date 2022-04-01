# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

from biocluster.config import Config
import xml.etree.ElementTree as ET
import re
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis_t import KGMLCanvas
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from reportlab.lib import colors
import collections
import json
from itertools import islice
import subprocess
import gridfs
import os
import sys
from api_base import ApiBase
import unittest

class KeggSpecies(ApiBase):
    def __init__(self, bind_object):
        '''
        设置数据库，连接到mongod数据库，涉及kegg_ko，kegg_gene，kegg_pathway_png三个collections
        '''
        super(KeggSpecies, self).__init__(bind_object)


        self.species_path_dict = dict()
        self.translate = dict()

    def import_species_info(self, species_json, pathway_dir, translate):
        self.species_json = species_json
        self.pathway_dir = pathway_dir
        self.get_spe_pathways(pathway_dir)
        self.get_spe_transtate(translate)
        self.get_kegg_class()

    def get_spe_pathways(self, pathway_dir):
        for spe in os.listdir(pathway_dir):
            with open(os.path.join(pathway_dir, spe), 'r') as f:
                path_list = [line.split("\t")[0][-5:] for line in f]
            self.species_path_dict[spe] = path_list


    def get_spe_transtate(self, translate):
        with open(translate, 'r') as f:
            f.readline()
            for line in f:
                for col in line.strip().split("\t"):
                    if " (" in col:
                        self.translate[col.split(" (")[0]] = col.split(" (")[1].rstrip(")")
        jsonout = open("translate.json",'w')
        jsonout.write(json.dumps(self.translate, indent=4))


    def get_kegg_class(self):
        '''
        获取kegg分类注释
        '''
        map2class = dict()
        with open(self.species_json) as f:
            root = json.load(f)

        species_list = list()
        euk = root["children"][0]
        classIs = euk['children']

        kegg_species_map = list()

        data_list = list()
        for classI in classIs:
            classI["spe_list"] = list()
            classIIs = classI['children']
            for classII in classIIs:
                classII["spe_list"] = list()
                classIIIs = classII['children']
                for classIII in classIIIs:
                    print classIII['name']
                    classIII["spe_list"] = list()
                    if 'children' in classIII:
                        species = classIII['children']
                        for specie in species:
                            species_abr = specie['name'].split("  ")[0]
                            classI["spe_list"].append(species_abr)
                            classII["spe_list"].append(species_abr)
                            classIII["spe_list"].append(species_abr)
                            data = {
                                "abr" : species_abr,
                                "level" : "species",
                                "name" : specie['name'].split("  ")[1].split(" (")[0],
                                "map_list":  self.species_path_dict[species_abr] if species_abr in self.species_path_dict else []
                            }
                        data_list.append(data)
                        list3 = list()
                        [list3.extend(self.species_path_dict[spe]) for spe in classIII["spe_list"] if spe in self.species_path_dict]
                        classIII["map_list"] = list(set(list3))

                        data = {
                            "abr" : classIII['name'].split(" (")[0].lower().replace(" ", "_"),
                            "level" : "classIII",
                            "name" : classIII['name'].split(" (")[0],
                            "map_list": classIII["map_list"]
                        }
                        data_list.append(data)
                    else:
                        species_abr = classIII['name'].split("  ")[0]
                        classI["spe_list"].append(species_abr)
                        classII["spe_list"].append(species_abr)
                        data = {
                            "abr" : species_abr,
                            "level" : "species",
                            "name" : specie['name'].split("  ")[1].split(" (")[0],
                            "map_list": self.species_path_dict[species_abr]
                        }
                        data_list.append(data)

                list2 = list()
                [list2.extend(self.species_path_dict[spe]) for spe in classII["spe_list"] if spe in self.species_path_dict]
                classII["map_list"] = list(set(list2))
                data = {
                    "abr" : classII['name'].split(" (")[0].lower().replace(" ", "_"),
                    "level" : "classII",
                    "name" : classII['name'].split(" (")[0],
                    "map_list": classII["map_list"]
                }
                data_list.append(data)
            list1 = list()
            [list1.extend(self.species_path_dict[spe]) for spe in classI["spe_list"] if spe in self.species_path_dict]
            classI["map_list"] = list(set(list1))
            data = {
                "abr" : classI['name'].split(" (")[0].lower().replace(" ", "_"),
                "level" : "classI",
                "name" : classI['name'].split(" (")[0],
                "map_list": classI["map_list"]
            }
            data_list.append(data)
        self.create_db_table('kegg_species_map', data_list)

        jsonout = open("euk_pathways.json",'w')
        jsonout.write(json.dumps(euk, indent=4))

        # print self.translate


        data_stat = list()
        for classI in classIs:
            stat_dict = dict()
            stat_dict['name'] = classI['name'].split(" (")[0]
            stat_dict['id'] = stat_dict['name'].lower().replace(" ", "_")
            stat_dict['translate'] = self.translate[stat_dict['name']]
            stat_dict['children'] = list()

            classIIs = classI['children']
            for classII in classIIs:
                stat2_dict = dict()
                stat2_dict['name'] = classII['name'].split(" (")[0]
                print stat2_dict['name']
                stat2_dict['id'] = stat2_dict['name'].lower().replace(" ", "_")
                if stat2_dict['name'] in self.translate:
                    stat2_dict['translate'] = self.translate[stat2_dict['name']]
                    stat2_dict['children'] = list()
                    classIIIs = classII['children']
                    for classIII in classIIIs:
                        stat3_dict = dict()
                        stat3_dict['name'] = classIII['name'].split(" (")[0]
                        stat3_dict['id'] = stat3_dict['name'].lower().replace(" ", "_")
                        if stat3_dict['name'] in self.translate:
                            stat3_dict['translate'] = self.translate[stat3_dict['name']]
                            stat2_dict['children'].append(stat3_dict)
                    stat_dict['children'].append(stat2_dict)
            data_stat.append(stat_dict)
        self.create_db_table('kegg_species', data_stat)

        detail_list  = list()
        for classI in classIs:
            detail1_dict = dict()
            detail1_dict['cl1_name'] = classI['name'].split(" (")[0]
            detail1_dict['cl1_id'] = detail1_dict['cl1_name'].lower().replace(" ", "_")
            detail1_dict['cl1_translate'] = self.translate.get(detail1_dict['cl1_name'], "")


            classIIs = classI['children']
            for classII in classIIs:
                detail2_dict = dict()
                detail2_dict['cl2_name'] = classII['name'].split(" (")[0]
                detail2_dict['cl2_id'] = detail2_dict['cl2_name'].lower().replace(" ", "_")
                detail2_dict['cl2_translate'] = self.translate.get(detail2_dict['cl2_name'], "")
                classIIIs = classII['children']
                for classIII in classIIIs:
                    detail3_dict = dict()
                    if 'children' in classIII:
                        detail3_dict['cl3_name'] = classIII['name'].split(" (")[0]
                        detail3_dict['cl3_id'] = detail3_dict['cl3_name'].lower().replace(" ", "_")
                        detail3_dict['cl3_translate'] = self.translate.get(detail3_dict['cl3_name'], "")
                        species = classIII['children']
                        for spe in species:
                            data_dict = dict({"version": "2020"})
                            data_dict.update(detail1_dict)
                            data_dict.update(detail2_dict)
                            data_dict.update(detail3_dict)
                            data_dict['spe_abr'] = spe['name'].split("  ")[0]
                            data_dict['spe_name'] = spe['name'].split("  ")[1].split(" (")[0].replace(" ", "_")

                            detail_list.append(data_dict)
                    else:
                        data_dict = dict({"version": "2020"})
                        data_dict.update(detail1_dict)
                        data_dict.update(detail2_dict)

                        data_dict['spe_abr'] = classIII['name'].split("  ")[0]
                        data_dict['spe_name'] = classIII['name'].split("  ")[1].split(" (")[0].replace(" ", "_")
                        detail_list.append(data_dict)

        self.create_db_table('kegg_species_detail', detail_list)




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def __init__(self, method_name, params_dict=None):
        self.params_dict = params_dict
        super(TestFunction, self).__init__(methodName=method_name)
        self.toolbox = KeggSpecies(None)

    def import_species_info(self):
        self.toolbox.import_species_info(**self.params_dict)


if __name__ == '__main__':
    if sys.argv[1] in ["-h", "-help", "--h", "--help"]:
        print "\n".join(["import_species_info"])
        if len(sys.argv) == 3:
            help(getattr(KeggSpecies, sys.argv[2]))

    # /mnt/ilustre/users/sanger-dev/app/database/Annotation/download2020/kegg_20200318$ python  /mnt/ilustre/users/sanger-dev/sg-users/liubinxu/work/SangerBiocluster2/SangerBiocluster/src/mbio/packages/annot_db/kegg_species.py  import_species_info species_json=br08610.json pathway_dir=pathway translate=species.tran.xls

    elif len(sys.argv) >= 3:
        suite = unittest.TestSuite()
        params_dict = dict()
        for par in sys.argv[2:]:
            params_dict.update({par.split("=")[0]: "=".join(par.split("=")[1:])})

        suite.addTest(TestFunction(sys.argv[1], params_dict))
        unittest.TextTestRunner(verbosity=2).run(suite)
