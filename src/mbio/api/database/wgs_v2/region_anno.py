# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# modified 20190325

import os
import json
import datetime
from api_base import ApiBase
from collections import defaultdict


class RegionAnno(ApiBase):
    """
    基因注释接口导表
    """
    def __init__(self, bind_object):
        super(RegionAnno, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.task_id

    def add_sg_region_anno(self, project_sn, task_id, params=None, name=None):
        """
        基因注释主表
        """
        params_json = {
            "task_id": task_id,
            "task_type": 2,
            "submit_location": "regionanno",
            "anno_id": "all",
            "region": "all",
            "region_type": "allregion",
            "chongmingming_result": ""
        }
        params_json = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": json.dumps(params, sort_keys=True, separators=(",", ":")) if params else params_json,
            "name": name if name else "origin_region_anno",
            "desc": "基因注释",
        }
        main_id = self.db["sg_region_anno"].insert_one(insert_data).inserted_id
        self.update_db_record("sg_region_anno", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_region_anno_detail(self, anno_id, pop_summary_path, species_version_id):
        """
        基因注释细节表
        pop.summary, 第三列：gene_name|GeneID|Genbank|transcript_id|protein
        """
        anno_id = self.check_objectid(anno_id)   # 检查id是否是OBjectID
        genome_id_ = self.check_objectid(species_version_id)  # 检查id是否是OBjectID
        self.check_exists(pop_summary_path)   # 检查文件是否存在
        data_list, gene_detail = [], []
        gene_detail_dict = {}
        with open(pop_summary_path, "r") as fr:
            lines = fr.readlines()[1:]
            for line in lines:
                if line.startswith("#") or line.startswith("@"):
                    continue
                item = line.strip().split("\t")
                gene_data_list = item[2].strip().split("|")
                if item[1] in gene_detail_dict.keys():
                    if item[0] == "" or item[0] == "--" or gene_data_list[3] == "" or gene_data_list[3] == "--":
                        transcrit_id_list = [item[1]]
                        new_transcrit_id_list = gene_detail_dict[item[1]]["transcripts"] + transcrit_id_list
                    else:
                        transcrit_id_list = gene_data_list[3].strip().split(";")
                        new_transcrit_id_list = gene_detail_dict[item[1]]["transcripts"] + transcrit_id_list
                    new_start = int(item[5]) if int(item[5]) < gene_detail_dict[item[1]]["start"] else gene_detail_dict[item[1]]["start"]
                    new_end = int(item[6]) if int(item[6]) > gene_detail_dict[item[1]]["end"] else gene_detail_dict[item[1]]["end"]
                    gene_detail_dict[item[1]] = {"start": new_start, "end": new_end, "transcripts": new_transcrit_id_list}
                else:
                    if len(gene_data_list) != 5:  # 解决奇葩的gff文件导致的gene_data_list不为5列的报错
                        print line
                        continue
                    if item[0] == "" or item[0] == "--" or gene_data_list[3] == "" or gene_data_list[3] == "--":
                        transcrit_id_list = [item[1]]
                    else:
                        transcrit_id_list = gene_data_list[3].strip().split(";")
                    gene_detail_dict[item[1]] = {"start": int(item[5]),"end": int(item[6]),"transcripts": transcrit_id_list}
        with open(pop_summary_path, "r") as r:
            lines = r.readlines()[1:]
            for line in lines:
                if line.startswith("#") or line.startswith("@"):
                    continue
                item = line.strip().split("\t")
                gene_data_list = item[2].strip().split("|")
                if len(gene_data_list) != 5:
                    continue
                if gene_data_list[3] == "" or gene_data_list[3] == "--":
                    transcrit_id_list = [item[1]]
                else:
                    transcrit_id_list = gene_data_list[3].strip().split(";")
                for transcrit_id in transcrit_id_list:
                    try:
                        insert_data = {
                            "anno_id": anno_id,
                            "gene_name": item[0],
                            "gene_id": item[1],
                            "gene_name_id": item[1] + ":" + item[0],
                            "transcrit_id": transcrit_id,
                            "chr": item[4],
                            "start": int(item[5]),
                            "end": int(item[6]),
                            "region": item[0] + ":" + item[5] + "-" + item[6],
                            "high": int(item[7]),
                            "moderate": int(item[8]),
                            "low": int(item[9]),
                            "modifier": int(item[10]),
                            "nr_id": item[11],
                            "nr_desc": item[12],
                            "uniprot_id": item[13],
                            "uniprot_desc": item[14],
                            "ko_id": item[15],
                            "ko_desc": item[16],
                            "go_term": item[17],
                            "go_desc": item[18],
                            "egg": item[19],
                            "egg_desc": item[20],
                            "pfam_acc": item[21],
                            "pfam_anno": item[22],
                            "interpro_acc": item[23],
                            "interpro_anno": item[24]
                        }
                    except:
                        insert_data = {
                            "anno_id": anno_id,
                            "gene_name": item[0],
                            "gene_id": item[1],
                            "gene_name_id": item[1] + ":" + item[0],
                            "transcrit_id": transcrit_id,
                            "chr": item[4],
                            "start": int(item[5]),
                            "end": int(item[6]),
                            "region": item[0] + ":" + item[5] + "-" + item[6],
                            "high": int(item[7]),
                            "moderate": int(item[8]),
                            "low": int(item[9]),
                            "modifier": int(item[10]),
                            "nr_id": item[11],
                            "nr_desc": item[12],
                            "uniprot_id": item[13],
                            "uniprot_desc": item[14],
                            "ko_id": item[15],
                            "ko_desc": item[16],
                            "go_term": item[17],
                            "go_desc": item[18],
                            "egg": item[19],
                            "egg_desc": item[20],
                            "pfam_acc": "--",
                            "pfam_anno": "--",
                            "interpro_acc": "--",
                            "interpro_anno": "--"
                        }
                    data_list.append(insert_data)
                gene_name, gen_id, gbank, tran_id, prote = self.check_split_info(item[2])
                gene_detail_ = {
                    "genome_version_id": genome_id_,
                    "gene_id": item[1],
                    "gene_name": item[0],
                    "genebank": gbank,
                    "chr": item[4],
                    "start": gene_detail_dict[item[1]]["start"],
                    "end": gene_detail_dict[item[1]]["end"],
                    "protein": prote,
                    "gene_length": int(item[6]) - int(item[5]),
                    "gene_type": item[3],
                    "transcripts": gene_detail_dict[item[1]]["transcripts"],
                    "gen_id": gen_id
                }
                result = self.col_find_one("sg_species_version_detail", {"genome_version_id": genome_id_, "gene_id": item[1]})
                if not result:
                    gene_detail.append(gene_detail_)
        if gene_detail:
            self.col_insert_data("sg_species_version_detail", gene_detail)
        else:
            print "gene在以前已经存在，gene_detail为空！"
        if len(data_list) == 0:
            self.bind_object.logger.info("gene_total注释列表为空！{}".format(pop_summary_path))
        else:
            self.col_insert_data("sg_region_anno_detail", data_list)

    def sg_region_anno_go_stat(self, anno_id, go_stat_detail, go_summary_path):
        """
        关联区域go分类统计表
        go.stat.detail
        go_summary_path: go summary结果pop.2.enrich
        """
        anno_id = self.check_objectid(anno_id)
        self.check_exists(go_stat_detail)
        self.check_exists(go_summary_path)
        data_list, value_list, go_id_dic = [], [], {}
        with open(go_summary_path, "r") as fr:
            lines = fr.readlines()[1:]
            for line in lines:
                item = line.strip().split("\t")
                go_id_dic[item[0]] = {"cate": item[2], "level2": item[-1]}
        with open(go_stat_detail, "r") as r:
            lines = r.readlines()
            if len(lines) == 0:
                return
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[0] in go_id_dic.keys():
                    level2_id = go_id_dic[item[0]]["level2"]
                    category = go_id_dic[item[0]]["cate"]
                else:
                    category, level2_id = "--", "--"
                insert_data = {
                    "anno_id": anno_id,
                    "go_id": item[0],
                    "des": item[1],
                    "eff_variant": int(item[2]),
                    "total_variant": int(item[3]),
                    "tran_list": item[-1].split(";"),
                    "category": category,
                    "level2_id": level2_id
                }
                data_list.append(insert_data)
                value_list.append((item[0] + ":" + item[1], int(item[2]), item[0]))
            if len(value_list) == 0 or len(data_list) == 0:
                self.bind_object.logger.info("anno_go注释列表为空！{}".format(go_stat_detail))
            else:
                self.col_insert_data("sg_region_anno_go_stat", data_list)
                sort_value = sorted(value_list, key=lambda value: value[1], reverse=True)  # 将数据按照value值倒序排列
                bar_id = self.sg_bar(self.task_id, anno_id, "GO分类统计柱形图", "", 2, location="region_anno_go_class")
                for tup in sort_value:  # 按倒序顺序存入数据库，方便前端取数据，categories对应
                    self.sg_bar_detail(bar_id, tup[2], tup[1], "false", "2", tooltip=tup[0])

    def sg_region_anno_kegg_stat(self, anno_id, kegg_stat_detail, pathway_dir, origin_path):
        """
        kegg分类统计表
        pathway_dir: 对象存储的pathway_dir的路径
        origin_path: 集群的pathway_dir的路径
        """
        anno_id = self.check_objectid(anno_id)
        self.check_exists(kegg_stat_detail)
        data_list, value_list = [], []
        with open(kegg_stat_detail, "r") as r:
            lines = r.readlines()
            if len(lines) == 0:
                return
            for line in lines[1:]:
                item = line.strip().split("\t")
                pdf_path = os.path.join(pathway_dir, item[0] + ".pdf")
                png_path = os.path.join(pathway_dir, item[0] + ".png")
                if not os.path.exists(os.path.join(origin_path, item[0] + ".pdf")) or not os.path.exists(os.path.join(origin_path, item[0] + ".png")):
                    pdf_path = "kegg数据库中没有找到对应的图片"
                    png_path = "kegg数据库中没有找到对应的图片"
                insert_data = {
                    "anno_id": anno_id,
                    "ko_id": item[0],
                    "des": item[1],
                    "eff_variant": int(item[2]),
                    "total_variant": int(item[3]),
                    "pdf_path": pdf_path,
                    "png_path": png_path,
                    "tran_list": item[-1],
                    "pathway_id": item[0].replace("ko", "pathway")
                }
                data_list.append(insert_data)
                value_list.append((item[0] + ":" + item[1], int(item[2]), item[0]))
            if len(value_list) == 0 or len(data_list) == 0:
                self.bind_object.logger.info("anno_kegg注释列表为空！{}".format(kegg_stat_detail))
            else:
                self.col_insert_data("sg_region_anno_kegg_stat", data_list)
                sort_value = sorted(value_list, key=lambda value: value[1], reverse=True)  # 将数据按照value值倒序排列
                bar_id = self.sg_bar(self.task_id, anno_id, "KEGG分类统计柱形图", "", 2, location="region_anno_kegg_class")
                for tup in sort_value:  # 按倒序顺序存入数据库，方便前端取数据，categories对应
                    self.sg_bar_detail(bar_id, tup[2], tup[1], "false", "2", tooltip=tup[0])

    def sg_region_anno_eggnog_stat(self, anno_id, eggnog_stat_detail):
        """
        eggnog分类统计表
        """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(eggnog_stat_detail)  # 检查文件是否存在
        data_list, categories = [], []  # x轴
        value_dict = defaultdict(list)
        with open(eggnog_stat_detail, "r") as r:
            lines = r.readlines()
            if len(lines) == 0:
                return
            for line in lines[1:]:
                item = line.strip().split("\t")
                categories.append(item[0])
                value_dict[item[0]] = [item[0] + ":" + item[1], int(item[2])]
                insert_data = {
                    "anno_id": anno_id,
                    "eggnog_id": item[0],
                    "des": item[1],
                    "eff_variant": int(item[2]),
                    "total_variant": int(item[3]),
                    "tran_list": item[-1]
                }
                data_list.append(insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("anno_eggnog注释列表为空！{}".format(eggnog_stat_detail))
            else:
                self.col_insert_data("sg_region_anno_eggnog_stat", data_list)
                categories.sort()  # 将X轴按字母顺序排序
                bar_id = self.sg_bar(self.task_id, anno_id, "EggNog分类统计柱形图", "", 3,
                                     location="region_anno_eggnog_class")
                for x in categories:   # 将value按照X轴排序后的顺序依次存放
                    self.sg_bar_detail(bar_id, value_dict[x][0], value_dict[x][1], "false", "3", categories=x)

    def sg_region_anno_pfam_stat(self, anno_id, pfam_stat_detail):
        """
        pfam分类统计表
        """
        anno_id = self.check_objectid(anno_id)
        self.check_exists(pfam_stat_detail)
        data_list, categories = [], []  # x轴
        value_dict = defaultdict(list)
        with open(pfam_stat_detail, "r") as r:
            lines = r.readlines()
            if len(lines) == 0:
                return
            for line in lines[1:]:
                item = line.strip().split("\t")
                categories.append(item[0])
                value_dict[item[0]] = [item[0] + ":" + item[1], int(item[2])]
                insert_data = {
                    "anno_id": anno_id,
                    "pfam_acc": item[0],
                    "pfam_anno": item[1],
                    "eff_variant": int(item[2]),
                    "total_variant": int(item[3]),
                    "tran_list": item[-1]
                }
                data_list.append(insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("anno_pfam注释列表为空！{}".format(pfam_stat_detail))
            else:
                self.col_insert_data("sg_region_anno_pfam_stat", data_list)
                categories.sort()  # 将X轴按字母顺序排序
                bar_id = self.sg_bar(self.task_id, anno_id, "pfam分类统计柱形图", "", 3,
                                     location="region_anno_pfam_class")
                for x in categories:   # 将value按照X轴排序后的顺序依次存放
                    self.sg_bar_detail(bar_id, value_dict[x][0], value_dict[x][1], "false", "3", categories=x)

    def check_split_info(self, strings):
        """
        该函数用于对字符串进行split，然后判断split后的列表的长度，长度如果为3就返回3个值，小于4的--串代替
        gene_name|GeneID|Genbank|transcript_id|protein
        :param strings:
        :return:
        """
        temp = strings.strip().split("|")
        try:
            temp1 = temp[0]
        except:
            temp1 = "--"
        try:
            temp2 = temp[1]
            if temp2 == "":
                temp2 = "--"
        except:
            temp2 = "--"
        try:
            temp3 = temp[2]
            if temp3 == "":
                temp3 = "--"
        except:
            temp3 = "--"
        try:
            temp4 = temp[3]
            if temp4 == "":
                temp4 = "--"
        except:
            temp4 = "--"
        try:
            temp5 = temp[4]
            if temp5 == "":
                temp5 = "--"
        except:
            temp5 = "--"
        return temp1, temp2, temp3, temp4, temp5

    def add_sg_anno_params(self, task_id, project_sn, types, params, desc, analysis_name, region_path_list, region_path):
        """
        添加注释表的主表
        :param task_id:
        :param project_sn:
        :param types:
        :param params:
        :param desc:
        :param analysis_name:
        :param region_path_list:
        :param region_path:
        :return:
        """
        data_list = []
        insert_data = {
            "task_id": task_id,
            "project_sn": project_sn,
            "type": types,
            "params": params,
            "desc": desc,
            "analysis_name": analysis_name,
            "region_path_list": region_path_list,
            "region_path": region_path
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_anno_params", data_list)
        return main_id

    def update_info(self, coll, main_id, update_dict):
        """
        更新sg_region_anno表的pop_summary
        """
        main_id = self.check_objectid(main_id)
        self.db[coll].update({"main_id": main_id}, {"$set": update_dict})


if __name__ == "__main__":
    a = RegionAnno(None)
    project_sn = "wgs_v2"
    task_id = "sanger_85433"
    species_version_id = "5ca50d9b17b2bf5bb702dc2d"
    # anno_id = a.add_sg_region_anno(project_sn, task_id)
    pop_summary_path = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/05.annovar/anno_count/pop.summary"
    # a.add_sg_region_anno_detail(anno_id, pop_summary_path, species_version_id)
    go_stat_detail = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/05.annovar/anno_count/pop.go.stat"
    go_summary_path = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/05.annovar/go_summary/pop.2.enrich"
    # pathway_dir = "s3://rerewrweset/files/wgs_test/wgs_test/wgs_test/hongdong/dadou/05.annovar/kegg_anno/pathway_dir/"
    kegg_stat_detail = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/05.annovar/anno_count/pop.kegg.stat"
    eggnog_stat_detail = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/05.annovar/anno_count/pop.eggnog.stat"
    pfam_stat_detail = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/05.annovar/anno_count/pop.pfam.stat"
    # origin_path = "/mnt/ilustre/users/sanger-dev/workspace/20190409/WgsV2_sanger_85433/output/05.annovar/kegg_anno/pathway_dir"
    anno_id = "5ccf99c717b2bf3f7608b212"
    pathway_dir = 's3://common/files/m_188/188_5cb3e5db67c6a/tsg_34072/workflow_results/05.annovar/kegg_anno/pathway_dir/'
    origin_path = '/mnt/ilustre/users/sanger-dev/workspace/20190505/WgsV2_tsg_34072/output/05.annovar/kegg_anno/pathway_dir'
    kegg_stat_detail = "/mnt/ilustre/users/sanger-dev/workspace/20190505/WgsV2_tsg_34072/output/05.annovar/anno_count/pop.kegg.stat"
    # a.sg_region_anno_go_stat(anno_id, go_stat_detail, go_summary_path)
    a.sg_region_anno_kegg_stat(anno_id, kegg_stat_detail, pathway_dir, origin_path)
    # a.sg_region_anno_eggnog_stat(anno_id, eggnog_stat_detail)
    # a.sg_region_anno_pfam_stat(anno_id, pfam_stat_detail)
