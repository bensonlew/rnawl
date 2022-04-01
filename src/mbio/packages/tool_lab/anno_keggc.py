# -*- coding: utf-8 -*-
import pandas as pd
from biocluster.api.database.base import Base

class AnnoKeggc(Base):
    def __init__(self):
        super(AnnoKeggc, self).__init__()
        self._project_type = 'metabolome'
        self.ct_dic = {}
        self.level_table = ""

    def run(self, table, level_path, stat_path,database_name="kegg_compound",query_field='entry'):
        """
        查看相关代谢物信息的主程序
        :param table: 详情表
        :param level_path: 输出层级分类文件
        :param stat_path: 输出二级分类统计文件
        :param  query_field: entry （对应compound id）； cas_id； formula；name
        :return:
        """
        data = pd.read_table(table)
        data.rename(columns={data.columns[0]:"Name"},inplace=True)
        #data['Name'] =data[data.columns[0]]
        self.level_table = open(level_path, "w")
        # stat_table = open(stat_path, "w")
        self.level_table.write("compound_id\tfirst_category\tsecond_category\tname\tmetab_id\thyperlink\n")
        data.apply(self.process, axis=1, **{"database_name":database_name, "query_field": query_field})
        self.level_table.close()
        with open(stat_path, 'w') as file:
            file.write("fist_category\tsecond_category\thyperlink\tmetab_id\tcount\n")
            for i in sorted(self.ct_dic.keys()):
                for j in sorted(self.ct_dic[i].keys()):
                    metab_list = sorted(list(self.ct_dic[i][j]["metab_id"]))
                    metab_str = ';'.join(metab_list)
                    metab_count = len(metab_list)
                    file.write("%s\t%s\t%s\t%s\t%s\n" % (i, j, self.ct_dic[i][j]["hyperlink"], metab_str, metab_count))

    def process(self, one,database_name='kegg_compound', query_field='entry'):
        """
        根据表格中的一行信息，查询kegg_compound参考库，获得全部化合物相关的信息
        此方法兼容了metab含多个compound_id的处理 @20191016
        :param one: 详情表的一行数据
        :return:
        """

        c = one['Name']
        if c in ["_",'-']:
            return
        if query_field == 'name':
            regex = "(?:;{0};|;{0}$|^{0};|^{0}$)".format(c)
            ref_result = self.ref_db[database_name].find_one({query_field: {"$regex": regex}})
        else:
            ref_result = self.ref_db[database_name].find_one({query_field: c})
        if not ref_result:
            return


        sec_link = ref_result['br2_fig_str']
        sec_link_list = sec_link.split(';')
        compound_id = ref_result['entry']
        name = ref_result['name']
        link = "http://www.kegg.jp/entry/" + compound_id

        if len(sec_link_list) > 1:
            fir_list = ref_result['br1_category'].split(';')
            sec_list = ref_result['br2_category'].split(';')
            for i in range(len(sec_link_list)):
                self.store(compound_id, fir_list[i], sec_list[i], name, c, link, sec_link_list[i])
            # self.store(compound_id, fir_list[0], sec_list[0], name, metab_id, link, sec_link_list[0])
            # self.store(compound_id, fir_list[1], sec_list[1], name, metab_id, link, sec_link_list[1])
        else:
            fir_c = ref_result['br1_category']
            sec_c = ref_result['br2_category']
            self.store(compound_id, fir_c, sec_c, name, c, link, sec_link)

    def store(self, compound_id, first_category, second_category, name, metab_id, link, sec_link):
        """
        将数据存储在diction self.ct_dic 中
        :param compound_id:
        :param first_category:
        :param second_category:
        :param name:
        :param metab_id:
        :param link:
        :param sec_link:
        :return:
        """
        if sec_link != "-":
            sec_link = "http://www.kegg.jp" + sec_link
        self.level_table.write(
            "%s\t%s\t%s\t%s\t%s\t%s\n" % (compound_id, first_category, second_category, name, metab_id, link))
        if first_category == "-":
            return
        if first_category not in self.ct_dic.keys():
            self.ct_dic[first_category] = {second_category: {"metab_id": set(), "hyperlink": sec_link}}
        elif second_category not in self.ct_dic[first_category].keys():
            self.ct_dic[first_category][second_category] = {"metab_id": set(), "hyperlink": sec_link}
        self.ct_dic[first_category][second_category]['metab_id'].add(metab_id)

if  __name__ == "__main__":
    import sys
    table = sys.argv[1]
    level_path = 'level.xls'
    stat_path = 'stat.xls'
    obj = AnnoKeggc()
    obj.run(table, level_path, stat_path)