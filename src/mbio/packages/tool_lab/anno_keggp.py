# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base
import pandas as pd
import os


class AnnoKeggp(Base):
    def __init__(self):
        super(AnnoKeggp, self).__init__()
        self._project_type = 'metabolome'
        self.stat_dic = {}
        self.detail_dict = {
            "metab_name": [],
            "first_category": [],
            "second_category": [],
            "description": [],
            "compound_id": [],
            'pathway_id':[]
        }


    def run(self, metab_names_str,detail_path,stat_path):
        metab_names = metab_names_str.split('\n')  ##;;;
        for metab_name in metab_names:
            if metab_name.strip() == '':
                continue
            name_info = self.find_name_compound(metab_name)
            if name_info:
                compound_id = name_info['entry']
                metab_ko_list = name_info['pathway'].split(';')
                self.detail_dict['metab_name'].append(metab_name)
                self.detail_dict['compound_id'].append(compound_id)
            else:
                self.detail_dict['metab_name'].append(metab_name)
                self.detail_dict['compound_id'].append('-')
                self.detail_dict['pathway_id'].append('-')
                self.detail_dict['first_category'].append('-')
                self.detail_dict['second_category'].append('-')
                self.detail_dict['description'].append('-')
                continue

            organisms = "All"
            if organisms not in  ["False", "All"]:
                ko_list = self.get_ko_list(organisms)
                ko_set = set(ko_list) & set(metab_ko_list)

            else:
                ko_set = set(metab_ko_list)

            ko_result = self.find_kegg(list(ko_set))
            # level_data = pd.DataFrame(columns=["pathway_id", "first_category", "second_category", "description", "compound_id", "metab_id", "count", "hyperlink"])

            tmp = {
                'pathway_id': [],
                'first_category' : [],
                'second_category' : [],
                'description' : []
            }

            for one in ko_result:
                tmp['pathway_id'].append(one['pathway_id'])
                tmp['first_category'].append(one['first_category'])
                tmp['second_category'].append(one['second_category'])
                tmp['description'].append(one['discription'])


                if self.stat_dic.has_key(one['first_category']):
                    if self.stat_dic[one['first_category']].has_key(one['second_category']):
                        self.stat_dic[one['first_category']][one['second_category']].append(compound_id)
                    else:
                        self.stat_dic[one['first_category']][one['second_category']] = [compound_id]
                else:
                    self.stat_dic[one['first_category']] = {
                        one['second_category']: [compound_id]
                    }

            if tmp['pathway_id']:
                self.detail_dict['pathway_id'].append(';'.join(tmp['pathway_id']))
            else:
                self.detail_dict['pathway_id'].append('-')
            if tmp['first_category']:
                self.detail_dict['first_category'].append(';'.join(tmp['first_category']))
            else:
                self.detail_dict['first_category'].append('-')

            if tmp['second_category']:
                self.detail_dict['second_category'].append(';'.join(tmp['second_category']))
            else:
                self.detail_dict['second_category'].append('-')

            if tmp['description']:
                self.detail_dict['description'].append(';'.join(tmp['description']))
            else:
                self.detail_dict['description'].append('-')


        detail_data = pd.DataFrame(data=self.detail_dict,
                        columns=["metab_name", "compound_id", "pathway_id","description", "first_category", "second_category"])

        detail_data.to_csv(detail_path, index=False, sep='\t')

        stat_file = open(stat_path, "w")
        stat_file.write("first_category\tsecond_category\tcompound_id\tcount\n")
        for fc in sorted(self.stat_dic.keys()):
            for sc in sorted(self.stat_dic[fc].keys()):
                if sc == "Global and overview maps":
                    continue  # 不对此二级分类作图 modified by GHD @ 20181205
                compound_list = self.stat_dic[fc][sc]
                comp_set = set(compound_list)
                count = len(comp_set)
                ids = ';'.join(comp_set)
                stat_file.write("%s\t%s\t%s\t%s\n" % (fc, sc, ids, count))
        stat_file.close()



    def get_ko_list(self, organisms):
        """
        根据kegg_organisms查询相关物种的通路信息
        :param organisms: 物种名称 XXX;xxx
        :return:
        """
        tmp_list = organisms.split(";")
        ref_org_db = self.ref_db["kegg_v94.2_organisms"]
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



    def find_name_compound(self, metab_name):
        database_name = 'kegg_v202007_br08001'
        query_field = 'name'
        regex = "(?:;{0};|;{0}$|^{0};|^{0}$)".format(metab_name)
        ref_result = self.ref_db[database_name].find_one({query_field: {"$regex": regex}})
        if ref_result:
            return ref_result
        else:
            return False

    def find_kegg(self, ko_list=None):
        """
        根据ko列表，查询kegg_pathway_level数据库
        :param ko_list:
        :return:
        """
        ref_kegg_db = self.ref_db["kegg_v94.2_pathway_level"]
        result = ref_kegg_db.find({"pathway_id": {"$in": ko_list}}, no_cursor_timeout=True)
        return result




