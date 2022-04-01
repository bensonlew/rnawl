# -*- coding: utf-8 -*-

import pandas as pd

class Metabsetc(object):
    def __init__(self):
        super(Metabsetc, self).__init__()
        self.category_dic = {}

    def run(self, set_table, overview_table, out_table):
        set_data = pd.read_table(set_table, header=None)
        set_list = set_data[0].tolist()
        overview_data = pd.read_table(overview_table)  #index_col=['metab_id']
        overview_data.drop_duplicates(['metab_id'],keep='first',inplace=True)   #hmdb和 Lib 分成了2行导致有重复metab_id.
        overview_data = overview_data.set_index(['metab_id'])
        overview_data = overview_data[overview_data['compound_second_category'] != '-']
        overview_list = overview_data.index.tolist()
        table_list = list(set(set_list) & set(overview_list))
        for indexs in table_list:
            metab_id = indexs
            first_category = overview_data["compound_first_category"][indexs]
            second_category = overview_data["compound_second_category"][indexs]
            self.count_category(metab_id, first_category, second_category)
        output_file = open(out_table, "w")
        table_head = ["first_category", "second_category", "metab_id", "count"]
        output_file.write("\t".join(table_head) + '\n')
        for fc in sorted(self.category_dic.keys()):
            for sc in sorted(self.category_dic[fc].keys()):
                metab_set = set(self.category_dic[fc][sc])
                count = len(metab_set)
                new_metab_id = ';'.join(list(metab_set))
                output_file.write("%s\t%s\t%s\t%s\n" % (fc, sc, new_metab_id, count))
        output_file.close()

    def count_category(self, metab_id, first_category, second_category):
        tmp_list = first_category.split(";")
        if len(tmp_list) == 1:
            self.count_first_category(metab_id, first_category, second_category)
        else:
            self.count_mul_first_category(metab_id, first_category, second_category)

    def count_first_category(self, metab_id, first_category, second_category):
        tmp_list = second_category.split(";")
        if self.category_dic.has_key(first_category):
            for one_second in tmp_list:
                if self.category_dic[first_category].has_key(one_second):
                    self.category_dic[first_category][one_second].append(metab_id)
                else:
                    self.category_dic[first_category][one_second] = [metab_id]
        else:
            for one_second in tmp_list:
                self.category_dic[first_category] = {one_second: [metab_id]}

    def count_mul_first_category(self, metab_id, first_category, second_category):
        tmp_list1 = first_category.split(";")
        tmp_list2 = second_category.split(";")
        for index,one_fc in enumerate(tmp_list1):
            if self.category_dic.has_key(one_fc):
                if self.category_dic[one_fc].has_key(tmp_list2[index]):
                    self.category_dic[one_fc][tmp_list2[index]].append(metab_id)
                else:
                    self.category_dic[one_fc][tmp_list2[index]] = [metab_id]
            else:
                self.category_dic[one_fc] = {tmp_list2[index]: [metab_id]}