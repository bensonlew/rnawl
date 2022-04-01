# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'

from api_base import ApiBase
import datetime
import os
import re
from collections import defaultdict
import math


class Anno(ApiBase):
    def __init__(self, bind_object):
        """
        wgs中所有的注释相关的导表
        __author__ = HONGDONG
        __last_modify__ = 20180419
        :param bind_object:
        """
        super(Anno, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_anno(self, task_id, project_sn, member_id=None, params=None, name=None):
        """
        基因注释主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param params:
        :param name:
        :return:
        """
        name = name if name else "origin_anno"
        params = params if params else ""
        main_id = self.add_main_table("sg_anno", task_id, project_sn, params, name, "基因注释主表", member_id)
        return main_id

    def add_sg_anno_detail(self, file_path, anno_id, task_id, species_version_id):
        """
        基因注释细节表---这里也进行了基因详情页的导表
        第三列： gene_name|GeneID|Genbank|transcript_id|protein
        :param file_path:
        :param species_version_id:
        :param anno_id:
        :param task_id:
        :return:
        """
        gene_detail = []
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        genome_id_ = self.check_objectid(species_version_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        # genome_id_ = self.add_sg_genome(species, genome_version, desc)
        with open(file_path, 'r') as r:
            data = r.readlines()
            for m in data:
                line = m.strip().split('\t')
                if re.match(r'^##', line[0]):
                    # info = line[0].strip().split(' ')  # 用于获取基因组信息
                    # species = info[0][2:]
                    # genome_version = info[1]
                    continue
                elif re.match(r'^#', line[0]):
                    continue
                if line[5] in ['chr1', 'chr4', 'chr5']:    # 由于pop_summary中第五列有字符串，需要黄总监那边处理
                    continue
                try:
                    insert_data = {
                        "anno_id": anno_id,
                        # "gene_id": line[1].strip().split(";")[0],
                        "gene_name": line[0],
                        "gene_id_name": str(line[1] + ":" + line[0]),
                        "gene_id": line[1],
                        "chr": line[4],
                        "start": int(line[5]),
                        "end": int(line[6]),
                        "high": int(line[7]),
                        "moderate": int(line[8]),
                        "low": int(line[9]),
                        "modifier": int(line[10]),
                        "nr_id": line[11],
                        "nr_desc": line[12],
                        "uniprot_id": line[13],
                        "uniprot_desc": line[14],
                        "go_term": line[17],
                        "go_desc": line[18],
                        "ko_id": line[15],
                        "ko_desc": line[16],
                        "egg": line[19],
                        "egg_desc": line[20]
                    }
                except:
                    continue
                data_list.append(insert_data)
                gene_name, gen_id, gbank, tran_id, prote = self.check_split_info(line[2])
                gene_detail_ = {
                    "genome_version_id": genome_id_,
                    "gene_id": line[1],
                    "gene_name": line[0],
                    "genebank": gbank,
                    "chr": line[4],
                    "start": int(line[5]),
                    "end": int(line[6]),
                    "protein": prote,
                    "gene_length": int(line[6])-int(line[5]),
                    "gene_type": line[3],
                    "transcripts": tran_id,
                    "gen_id": gen_id
                }
                if not self.col_find_one("sg_species_version_detail", {"genome_version_id": genome_id_, "gene_id": line[1]}):
                    gene_detail.append(gene_detail_)
        self.col_insert_data("sg_anno_detail", data_list)
        # self.add_sg_version(task_id, genome_id_)
        if gene_detail:
            self.col_insert_data("sg_species_version_detail", gene_detail)
        else:
            print "gene在以前已经存在，gene_detail为空！"

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
            temp1 = '--'
        try:
            temp2 = temp[1]
            if temp2 == "":
                temp2 = '--'
        except:
            temp2 = '--'
        try:
            temp3 = temp[2]
            if temp3 == "":
                temp3 = '--'
        except:
            temp3 = '--'
        try:
            temp4 = temp[3]
            if temp4 == "":
                temp4 = '--'
        except:
            temp4 = '--'
        try:
            temp5 = temp[4]
            if temp5 == "":
                temp5 = '--'
        except:
            temp5 = '--'
        return temp1, temp2, temp3, temp4, temp5

    def add_sg_version(self, task_id, genome_id):
        """
        添加基因组版本的主表
        :param task_id:
        :param genome_id:
        :return:
        """
        genome_id = self.check_objectid(genome_id)
        self.col_insert_data("sg_version", [{"created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                                             "task_id": task_id, "genome_id": genome_id}])

    def add_sg_genome(self, species, genome_version, desc):
        """
        添加基因信息主表
        :param species:
        :param genome_version:
        :param desc:
        :return:
        """
        result = self.col_find_one("sg_genome", {"species": species, "genome_version": genome_version})
        if result:
            __id = result["_id"]
        else:
            __id = self.col_insert_data("sg_genome", [{
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "species": species,
                "genome_version": genome_version, "desc": desc}])
        return __id

    def update_specimen_type(self, task_id, wp, mp, wb, mb):
        """
        更新参与分析的样本，增加type为1，便于关联区域详情筛选功能样本展示正确
        """
        results = self.db["sg_specimen"].find({"task_id": task_id})
        update_dict = {"type": 1}
        if results:
            if wp:
                query_dict = {"task_id": task_id, "old_name": wp}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)
            if mp:
                query_dict = {"task_id": task_id, "old_name": mp}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)
            if wb:
                query_dict = {"task_id": task_id, "old_name": wb}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)
            if mb:
                query_dict = {"task_id": task_id, "old_name": mb}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)

    def sg_anno_go_stat(self, anno_id, go_stat_detail):
        """
        go分类统计表
        pop.go.final.stat.detail
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(go_stat_detail)  # 检查文件是否存在
        result = self.col_find_one("sg_anno", {"main_id": anno_id})
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
                    self.col_insert_data("sg_anno_go_stat", data_list)
                    sort_value = sorted(value_list, key=lambda value: value[1], reverse=True)  # 将数据按照value值倒序排列
                    # for tup in sort_value:
                    #     categories.append(tup[0])
                    bar_id = self.sg_bar(task_id, anno_id, 'GO分类统计柱形图', "", 2, location='anno_go_class')
                    for tup in sort_value:  # 按倒序顺序存入数据库，方便前端取数据，categories对应
                        self.sg_bar_detail(bar_id, tup[2], tup[1], "false", "2", tooltip=tup[0])

    def sg_anno_kegg_stat(self, anno_id, kegg_stat_detail, origin_path, pathway_dir):
        """
        关联区域kegg分类统计表
        pop.kegg.final.stat.detail
        origin_path png图片所在output路径
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(kegg_stat_detail)  # 检查文件是否存在
        result = self.col_find_one("sg_anno", {"main_id": anno_id})
        task_id = result["task_id"]
        # pathway_dir = "/rerewrweset" + pathway_dir.split("rerewrweset")[1]
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
                    if not os.path.exists(os.path.join(origin_path, ko_id + '.pdf')) \
                            or not os.path.exists(os.path.join(origin_path, ko_id + '.png')):
                        pdf_path = "kegg数据库中没有找到对应的图片"
                        png_path = "kegg数据库中没有找到对应的图片"
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
                    self.col_insert_data("sg_anno_kegg_stat", data_list)
                    sort_value = sorted(value_list, key=lambda value: value[1], reverse=True)  # 将数据按照value值倒序排列
                    # for tup in sort_value:
                    #     categories.append(tup[0])
                    bar_id = self.sg_bar(task_id, anno_id, 'KEGG分类统计柱形图', "", 2, location='anno_kegg_class')
                    for tup in sort_value:  # 按倒序顺序存入数据库，方便前端取数据，categories对应
                        self.sg_bar_detail(bar_id, tup[2], tup[1], "false", '2', tooltip=tup[0])

    def sg_anno_eggnog_stat(self, anno_id, eggnog_stat_detail):
        """
        关联区域eggnog分类统计表
        pop.eggnog.final.stat.detail
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(eggnog_stat_detail)  # 检查文件是否存在
        result = self.col_find_one("sg_anno", {"main_id": anno_id})
        task_id = result["task_id"]
        data_list = []
        categories = []  # x轴
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
                    value_dict[eggnog_id] = [eggnog_id + ':' + des, int(tmp[2])]  # 每个图例对应的值
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
                    self.col_insert_data("sg_anno_eggnog_stat", data_list)
                    categories.sort()  # 将X轴按字母顺序排序
                    bar_id = self.sg_bar(task_id, anno_id, 'EggNog分类统计柱形图', "", 3,
                                         location='anno_eggnog_class')
                    for x in categories:  # 将value按照X轴排序后的顺序依次存放
                        self.sg_bar_detail(bar_id, value_dict[x][0], value_dict[x][1], "false", "3", categories=x)


if __name__ == "__main__":
    project_sn = '188_5b48676c687ce'
    task_id = 'tsg_31253'
    member_id = None
    t = Anno(None)
    t.project_type = "dna_gmap"
    anno_id = t.add_sg_anno(task_id, project_sn, member_id, "origin_params")
    t.add_sg_anno_detail("/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5b48676c687ce/tsg_30967/workflow_results/05.annovar/anno_count/pop.summary", anno_id, task_id, "5ad05fd6a4e1af4315161715")
    # t.sg_anno_go_stat(anno_id, "/mnt/ilustre/users/sanger-dev/workspace/20180417/Single_annovar_module_20180416/Annovar/output/go_anno/pop.go.final.stat.detail")
    # t.sg_anno_kegg_stat(anno_id, "/mnt/ilustre/users/sanger-dev/workspace/20180417/Single_annovar_module_20180416/Annovar/output/kegg_anno/pop.kegg.final.stat.detail", "/mnt/ilustre/users/sanger-dev/workspace/20180417/Single_annovar_module_20180416/Annovar/output/kegg_anno/pathway_dir")
    # t.sg_anno_eggnog_stat(anno_id, "/mnt/ilustre/users/sanger-dev/workspace/20180417/Single_annovar_module_20180416/Annovar/output/eggnog_anno/pop.eggnog.final.stat.detail")
