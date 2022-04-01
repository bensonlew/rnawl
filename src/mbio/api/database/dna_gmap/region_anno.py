# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
# modified 2018.02.26

from api_base import ApiBase
from collections import defaultdict
import os
import re
import datetime


class RegionAnno(ApiBase):
    """
    BSA关联区域注释导表
    """
    def __init__(self, bind_object):
        super(RegionAnno, self).__init__(bind_object)
        self._project_type = "dna_gmap"

    def add_sg_region_anno(self, region_id, trait, name=None, params=None):
        """
        关联区域注释主表
        """
        region_id = self.check_objectid(region_id)  # 检查id是否是OBjectID
        result = self.col_find_one("sg_qtl", {"main_id": region_id})
        task_id = result["task_id"]
        project_sn = result["project_sn"]
        name = result["name"]
        name = trait + "_" + name.split("QtlAnalysis_")[-1]
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "name": name if name else "origin_region_annotation",
            "params": params if params else "null",
            "desc": "关联区域基因注释主表",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        main_id = self.db['sg_region_anno'].insert_one(insert_data).inserted_id
        self.update_db_record("sg_region_anno", {"_id": main_id}, {"main_id": main_id, "region_id": region_id})
        return main_id

    def sg_region_anno_detail(self, anno_id, gene_total_path):
        """
        关联区域基因注释细节表
        region.threshold.gene.total
        """
        anno_id = self.check_objectid(anno_id)   # 检查id是否是OBjectID
        self.check_exists(gene_total_path)   # 检查文件是否存在
        data_list = []
        with open(gene_total_path, 'r') as r:
            lines = r.readlines()[1:]
            for line in lines:
                if line.startswith("#") or line.startswith("@"):
                    pass
                else:
                    tmp = line.strip().split('\t')
                    insert_data = {
                        "anno_id": anno_id,
                        "gene_name": tmp[0],
                        "gene_id": tmp[1],
                        "chr": tmp[4],
                        "start": int(tmp[5]),
                        "end": int(tmp[6]),
                        "high": int(tmp[7]),
                        "moderate": int(tmp[8]),
                        "low": int(tmp[9]),
                        "modifier": int(tmp[10]),
                        "nr_id": tmp[11],
                        "nr_desc": tmp[12],
                        "uniprot_id": tmp[13],
                        "uniprot_desc": tmp[14],
                        "ko_id": tmp[15],
                        "ko_desc": tmp[16],
                        "go_term": tmp[17],
                        "go_desc": tmp[18],
                        "egg": tmp[19],
                        "egg_desc": tmp[20]
                    }
                    data_list.append(insert_data)
        if len(data_list) == 0:  # add by hongdong 20180409
            self.bind_object.logger.info("gene_total注释列表为空！{}".format(gene_total_path))
        else:
            self.col_insert_data("sg_region_anno_detail", data_list)

    def sg_region_anno_go_stat(self, anno_id, go_stat_detail):
        """
        关联区域go分类统计表
        region.threshold.gene.go.final.stat.detail
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(go_stat_detail)  # 检查文件是否存在
        result = self.col_find_one("sg_region_anno", {"main_id": anno_id})
        task_id = result["task_id"]
        data_list = []
        categories = []  # x轴
        value_list = []
        with open(go_stat_detail, 'r') as r:
            lines = r.readlines()
            if len(lines) == 0:  # add by wzy 20180330  如果为空，则只导主表
                pass
            else:
                for line in lines[1:]:
                    tmp = line.strip().split('\t')
                    go_id = tmp[0].strip("\"")
                    des = tmp[1].strip("\"")
                    insert_data = {
                        "anno_id": anno_id,
                        "go_id": go_id,
                        "des": des,
                        "eff_variant": int(tmp[2]),
                        "total_variant": int(tmp[3])
                    }
                    data_list.append(insert_data)
                    value_list.append((go_id + ':' + des, int(tmp[2]), go_id))
                if len(value_list) == 0 or len(data_list) == 0:
                    self.bind_object.logger.info("anno_go注释列表为空！{}".format(go_stat_detail))
                else:
                    self.col_insert_data("sg_region_anno_go_stat", data_list)
                    sort_value = sorted(value_list, key=lambda value: value[1], reverse=True)  # 将数据按照value值倒序排列
                    # for tup in sort_value:
                    #     categories.append(tup[0])
                    bar_id = self.sg_bar(task_id, anno_id, 'GO分类统计柱形图', "", 2, location='region_anno_go_class')
                    for tup in sort_value:  # 按倒序顺序存入数据库，方便前端取数据，categories对应
                        self.sg_bar_detail(bar_id, tup[2], tup[1], "false", "2", tooltip=tup[0])

    def sg_region_anno_kegg_stat(self, anno_id, kegg_stat_detail, pathway_dir):
        """
        关联区域kegg分类统计表
        region.threshold.gene.kegg.final.stat.detail
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(kegg_stat_detail)  # 检查文件是否存在
        result = self.col_find_one("sg_region_anno", {"main_id": anno_id})
        task_id = result["task_id"]
        # if pathway_dir.startswith("/rerewrweset"):
        #     pathway_dir = "/rerewrweset" + pathway_dir.split("rerewrweset")[-1]
        data_list = []
        categories = []  # x轴
        value_list = []
        with open(kegg_stat_detail, 'r') as r:
            lines = r.readlines()
            if len(lines) == 0:  # add by wzy 20180330  如果为空，则只导主表
                pass
            else:
                for line in lines[1:]:
                    tmp = line.strip().split('\t')
                    ko_id = tmp[0].strip("\"")
                    des = tmp[1].strip("\"")
                    pdf_path = os.path.join(pathway_dir, ko_id + '.pdf')
                    png_path = os.path.join(pathway_dir, ko_id + '.png')
                    insert_data = {
                        "anno_id": anno_id,
                        "ko_id": ko_id,
                        "des": des,
                        "eff_variant": int(tmp[2]),
                        "total_variant": int(tmp[3]),
                        "pdf_path": pdf_path,
                        "png_path": png_path
                    }
                    data_list.append(insert_data)
                    value_list.append((ko_id + ':' + des, int(tmp[2]), ko_id))
                if len(value_list) == 0 or len(data_list) == 0:
                    self.bind_object.logger.info("anno_kegg注释列表为空！{}".format(kegg_stat_detail))
                else:
                    self.col_insert_data("sg_region_anno_kegg_stat", data_list)
                    sort_value = sorted(value_list, key=lambda value: value[1], reverse=True)  # 将数据按照value值倒序排列
                    # for tup in sort_value:
                    #     categories.append(tup[0])
                    bar_id = self.sg_bar(task_id, anno_id, 'KEGG分类统计柱形图', "", 2, location='region_anno_kegg_class')
                    for tup in sort_value:  # 按倒序顺序存入数据库，方便前端取数据，categories对应
                        self.sg_bar_detail(bar_id, tup[2], tup[1], "false", '2', tooltip=tup[0])

    def sg_region_anno_eggnog_stat(self, anno_id, eggnog_stat_detail):
        """
        关联区域eggnog分类统计表
        region.threshold.gene.eggnog.final.stat.detail
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(eggnog_stat_detail)  # 检查文件是否存在
        result = self.col_find_one("sg_region_anno", {"main_id": anno_id})
        task_id = result["task_id"]
        data_list = []
        categories = []   # x轴
        value_dict = defaultdict(list)
        with open(eggnog_stat_detail, 'r') as r:
            lines = r.readlines()
            if len(lines) == 0:  # add by wzy 20180330  如果为空，则只导主表
                pass
            else:
                for line in lines[1:]:
                    tmp = line.strip().split('\t')
                    eggnog_id = tmp[0].strip("\"")
                    categories.append(eggnog_id)
                    des = tmp[1].strip("\"")
                    value_dict[eggnog_id] = [eggnog_id + ':' + des, int(tmp[2])]   # 每个图例对应的值
                    insert_data = {
                        "anno_id": anno_id,
                        "eggnog_id": eggnog_id,
                        "des": des,
                        "eff_variant": int(tmp[2]),
                        "total_variant": int(tmp[3])
                    }
                    data_list.append(insert_data)
                if len(data_list) == 0:
                    self.bind_object.logger.info("anno_eggnog注释列表为空！{}".format(eggnog_stat_detail))
                else:
                    self.col_insert_data("sg_region_anno_eggnog_stat", data_list)
                    categories.sort()  # 将X轴按字母顺序排序
                    bar_id = self.sg_bar(task_id, anno_id,  'EggNog分类统计柱形图', "", 3,
                                         location='region_anno_eggnog_class')
                    for x in categories:   # 将value按照X轴排序后的顺序依次存放
                        self.sg_bar_detail(bar_id, value_dict[x][0], value_dict[x][1], "false", "3", categories=x)


if __name__ == "__main__":
    a = RegionAnno(None)
    project_sn = '188_5b48676c687ce'
    task_id = 'tsg_31253'
    region_id = "5b63af56a4e1af13fceff8c3"
    gene_total_path = "/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b48676c687ce/tsg_31236/interaction_results/QtlAnalysis/QtlAnalysis_20180803_092646/TGW.gene.total"
    # anno_id = a.add_sg_region_anno(region_id)
    anno_id = "5b642e1da4e1af4a1026de43"
    a.sg_region_anno_detail(anno_id, gene_total_path)
    # go_path = '/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b48676c687ce/tsg_31236/interaction_results/QtlAnalysis/QtlAnalysis_20180803_092646/TGW.go.final.stat.detail'
    # kegg_path = '/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b48676c687ce/tsg_31236/interaction_results/QtlAnalysis/QtlAnalysis_20180803_092646/TGW'
    # eggnog_path = '/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b48676c687ce/tsg_31236/interaction_results/QtlAnalysis/QtlAnalysis_20180803_092646/TGW'
    # pathway_dir = '/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b48676c687ce/tsg_31236/interaction_results/QtlAnalysis/QtlAnalysis_20180803_092646/pathway_dir'
    # a.sg_region_anno_go_stat(anno_id, go_path)
    # a.sg_region_anno_kegg_stat(anno_id, kegg_path, pathway_dir)
    # a.sg_region_anno_eggnog_stat(anno_id, eggnog_path)
