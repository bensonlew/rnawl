# -*- coding: utf-8 -*-
import pandas as pd
from biocluster.api.database.base import Base

class AnnoKeggc(Base):
    def __init__(self):
        super(AnnoKeggc, self).__init__()
        self._project_type = 'metabolome'
        self.ct_dic = {}
        self.level_table = ""

    def run(self, table, level_path, stat_path,database_name="kegg_compound"):
        """
        查看相关代谢物信息的主程序
        :param table: 详情表
        :param level_path: 输出层级分类文件
        :param stat_path: 输出二级分类统计文件
        :return:
        """
        data = pd.read_table(table)
        self.level_table = open(level_path, "w")
        # stat_table = open(stat_path, "w")
        self.level_table.write("compound_id\tfirst_category\tsecond_category\tname\tmetab_id\thyperlink\n")
        data.apply(self.process, axis=1, **{"database_name":database_name})
        self.level_table.close()
        with open(stat_path, 'w') as file:
            file.write("fist_category\tsecond_category\thyperlink\tmetab_id\tcount\n")
            for i in sorted(self.ct_dic.keys()):
                for j in sorted(self.ct_dic[i].keys()):
                    metab_list = sorted(list(self.ct_dic[i][j]["metab_id"]))
                    metab_str = ';'.join(metab_list)
                    metab_count = len(metab_list)
                    file.write("%s\t%s\t%s\t%s\t%s\n" % (i, j, self.ct_dic[i][j]["hyperlink"], metab_str, metab_count))

    def process(self, one,database_name='kegg_compound'):
        """
        根据表格中的一行信息，查询kegg_compound参考库，获得全部化合物相关的信息
        此方法兼容了metab含多个compound_id的处理 @20191016
        :param one: 详情表的一行数据
        :return:
        """
        if one['Metabolite'].startswith("pos") or one['Metabolite'].startswith("neg"):
            return  # 不要不带代谢物名称的化合物
        if one['KEGG Compound ID'] in ["_", "-"]:
            return
        # ref_result = self.ref_db[database_name].find_one({"entry": one['KEGG Compound ID']})
        for c in one["KEGG Compound ID"].split(";"):  # 允许一个代谢物有多个compound id 出现 by ghd @ 20191014
            ref_result = self.ref_db[database_name].find_one({"entry": c})
            if not ref_result:
                continue
            # if not ref_result:
            #     return
            sec_link = ref_result['br2_fig_str']
            sec_link_list = sec_link.split(';')
            compound_id = ref_result['entry']
            name = ref_result['name']
            link = "http://www.kegg.jp/entry/" + compound_id
            metab_id = one['metab_id']
            if len(sec_link_list) > 1:
                fir_list = ref_result['br1_category'].split(';')
                sec_list = ref_result['br2_category'].split(';')
                for i in range(len(sec_link_list)):
                    self.store(compound_id, fir_list[i], sec_list[i], name, metab_id, link, sec_link_list[i])
                # self.store(compound_id, fir_list[0], sec_list[0], name, metab_id, link, sec_link_list[0])
                # self.store(compound_id, fir_list[1], sec_list[1], name, metab_id, link, sec_link_list[1])
            else:
                fir_c = ref_result['br1_category']
                sec_c = ref_result['br2_category']
                self.store(compound_id, fir_c, sec_c, name, metab_id, link, sec_link)

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
