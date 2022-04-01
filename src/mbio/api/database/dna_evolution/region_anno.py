# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0816

from api_base import ApiBase
from collections import defaultdict
import os
import datetime
import json


class RegionAnno(ApiBase):
    """
    群体进化，基因注释接口导表
    """
    def __init__(self, bind_object):
        super(RegionAnno, self).__init__(bind_object)
        self._project_type = "dna_evolution"

    def add_sg_region_anno(self, project_sn, task_id, params=None, name=None):
        """
        基因注释主表
        """
        if params:
            new_params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            new_params = "null"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": new_params,
            "name": name if name else "origin_sg_region_anno",
            "desc": "基因注释接口",
        }
        main_id = self.db['sg_region_anno'].insert_one(insert_data).inserted_id
        self.update_db_record("sg_region_anno", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_region_anno_detail(self, anno_id, pop_summary_path, species_version_id):
        """
        基因注释细节表
        pop.summary
        """
        gene_detail = []
        anno_id = self.check_objectid(anno_id)   # 检查id是否是OBjectID
        genome_id_ = self.check_objectid(species_version_id)  # 检查id是否是OBjectID
        self.check_exists(pop_summary_path)   # 检查文件是否存在
        data_list = []
        gene_detail_dict = {}
        with open(pop_summary_path, 'r') as fr:
            lines = fr.readlines()[1:]
            for line in lines:
                if line.startswith("#") or line.startswith("@"):
                    pass
                else:
                    tmp = line.strip().split('\t')
                    gene_data_list = tmp[2].strip().split('|')
                    if tmp[1] in gene_detail_dict.keys():
                        if tmp[0] == '' or tmp[0] == '--' or gene_data_list[3] == '' or gene_data_list[3] == '--':
                            transcrit_id_list = [tmp[1]]
                            new_transcrit_id_list = gene_detail_dict[tmp[1]]["transcripts"] + transcrit_id_list
                        else:
                            transcrit_id_list = gene_data_list[3].strip().split(';')
                            new_transcrit_id_list = gene_detail_dict[tmp[1]]["transcripts"] + transcrit_id_list
                        if int(tmp[5]) < gene_detail_dict[tmp[1]]["start"]:
                            new_start = int(tmp[5])
                        else:
                            new_start = gene_detail_dict[tmp[1]]["start"]
                        if int(tmp[6]) > gene_detail_dict[tmp[1]]["end"]:
                            new_end = int(tmp[6])
                        else:
                            new_end = gene_detail_dict[tmp[1]]["end"]
                        gene_detail_dict[tmp[1]] = {"start": new_start, "end": new_end, "transcripts": new_transcrit_id_list}
                    else:
                        if len(gene_data_list) != 5:  # 解决奇葩的gff文件导致的gene_data_list不为5列的报错
                            print line
                            continue
                        if tmp[0] == '' or tmp[0] == '--' or gene_data_list[3] == '' or gene_data_list[3] == '--':
                            transcrit_id_list = [tmp[1]]
                        else:
                            transcrit_id_list = gene_data_list[3].strip().split(';')
                        gene_detail_dict[tmp[1]] = {"start": int(tmp[5]),"end": int(tmp[6]),"transcripts": transcrit_id_list}
        with open(pop_summary_path, 'r') as r:
            lines = r.readlines()[1:]
            for line in lines:
                if line.startswith("#") or line.startswith("@"):
                    pass
                else:
                    tmp = line.strip().split('\t')
                    gene_data_list = tmp[2].strip().split('|')
                    if len(gene_data_list) != 5:
                        continue
                    if gene_data_list[3] == '' or gene_data_list[3] == '--':
                        transcrit_id_list = [tmp[1]]
                    else:
                        transcrit_id_list = gene_data_list[3].strip().split(';')
                        # transcrit_id = transcrit_id_list[0]
                    for transcrit_id in transcrit_id_list:
                        insert_data = {
                            "anno_id": anno_id,
                            "gene_name": tmp[0],
                            "gene_id": tmp[1],
                            "gene_name_id": tmp[1] + ":" + tmp[0],
                            "transcrit_id": transcrit_id,
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
                    gene_name, gen_id, gbank, tran_id, prote = self.check_split_info(tmp[2])
                    gene_detail_ = {
                        "genome_version_id": genome_id_,
                        "gene_id": tmp[1],
                        "gene_name": tmp[0],
                        "genebank": gbank,
                        "chr": tmp[4],
                        "start": gene_detail_dict[tmp[1]]["start"],
                        "end": gene_detail_dict[tmp[1]]["end"],
                        "protein": prote,
                        "gene_length": int(tmp[6]) - int(tmp[5]),
                        "gene_type": tmp[3],
                        "transcripts": gene_detail_dict[tmp[1]]["transcripts"],
                        "gen_id": gen_id
                    }
                    if not self.col_find_one("sg_species_version_detail",
                                             {"genome_version_id": genome_id_, "gene_id": tmp[1]}):
                            gene_detail.append(gene_detail_)
        if gene_detail:
            self.col_insert_data("sg_species_version_detail", gene_detail)
        else:
            print "gene在以前已经存在，gene_detail为空！"
        if len(data_list) == 0:  # add by hongdong 20180409
            self.bind_object.logger.info("gene_total注释列表为空！{}".format(pop_summary_path))
        else:
            self.col_insert_data("sg_region_anno_detail", data_list)

    def sg_region_anno_go_stat(self, task_id, anno_id, go_stat_detail, go_summary_path, pop_summary=None):
        """
        关联区域go分类统计表
        /mnt/ilustre/users/sanger-dev/workspace/20180817/Single_region_anno_0817_220180817/RegionAnno/output/select_region.go.stat.detail
        go_summary_path: go summary结果pop.2.enrich
        """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(go_stat_detail)  # 检查文件是否存在
        data_list = []
        value_list = []
        go_id_dic = {}
        gene_dict = {}
        if pop_summary:
            gene_dict = self.get_genes(pop_summary, 'go')
        with open(go_summary_path, 'r') as fr:
            lines = fr.readlines()[1:]
            for line in lines:
                tmp = line.strip().split('\t')
                level2_id = [tmp[0]]
                desc = [tmp[2]]
                go_id = tmp[6]
                go_id_list = go_id.strip().split(',')
                for i in go_id_list:
                    if i in go_id_dic.keys():
                        new_desc = go_id_dic[i]["desc"] + desc
                        new_level2_id = go_id_dic[i]["level2_id"] + level2_id
                        go_id_dic[i] = {"desc": new_desc, "level2_id": new_level2_id}
                    else:
                        go_id_dic[i] = {"desc": desc, "level2_id": level2_id}
        with open(go_stat_detail, 'r') as r:
            lines = r.readlines()
            if len(lines) == 0:  # add by wzy 20180330  如果为空，则只导主表
                pass
            else:
                for line in lines[1:]:
                    tmp = line.strip().split('\t')
                    go_id = tmp[0].strip("\"")
                    des = tmp[1].strip("\"")
                    if go_id in go_id_dic.keys():
                        level2_desc = go_id_dic[go_id]["desc"]
                        level2_id = go_id_dic[go_id]["level2_id"]
                    else:
                        level2_desc = "--"
                        level2_id = "--"
                    insert_data = {
                        "anno_id": anno_id,
                        "go_id": go_id,
                        "des": des,
                        "eff_variant": int(tmp[2]),
                        "total_variant": int(tmp[3]),
                        # "transcript_id": tmp[6],
                        "level2_desc": level2_desc,
                        "level2_id": level2_id
                    }
                    if pop_summary:
                        if go_id in gene_dict.keys():
                            trans = gene_dict[go_id]
                        else:
                            trans = []
                        insert_data.update({"gene_id": trans})
                    else:
                        insert_data.update({"gene_id": tmp[6]})
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

    def sg_region_anno_kegg_stat(self, task_id, anno_id, kegg_stat_detail, origin_path, pathway_dir):
        """
        kegg分类统计表
        /mnt/ilustre/users/sanger-dev/workspace/20180817/Single_region_anno_0817_220180817/RegionAnno/output/select_region.kegg.stat.detail
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(kegg_stat_detail)  # 检查文件是否存在
        # pathway_dir = "/rerewrweset" + pathway_dir.split("rerewrweset")[1]
        data_list = []
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
                    self.col_insert_data("sg_region_anno_kegg_stat", data_list)
                    sort_value = sorted(value_list, key=lambda value: value[1], reverse=True)  # 将数据按照value值倒序排列
                    bar_id = self.sg_bar(task_id, anno_id, 'KEGG分类统计柱形图', "", 2, location='region_anno_kegg_class')
                    for tup in sort_value:  # 按倒序顺序存入数据库，方便前端取数据，categories对应
                        self.sg_bar_detail(bar_id, tup[2], tup[1], "false", '2', tooltip=tup[0])

    def sg_region_anno_eggnog_stat(self, task_id, anno_id, eggnog_stat_detail, pop_summary=None):
        """
        eggnog分类统计表
        /mnt/ilustre/users/sanger-dev/workspace/20180817/Single_region_anno_0817_220180817/RegionAnno/output/select_region.eggnog.stat.detail
         """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(eggnog_stat_detail)  # 检查文件是否存在
        data_list = []
        categories = []   # x轴
        value_dict = defaultdict(list)
        gene_dict = {}
        if pop_summary:
            gene_dict = self.get_genes(pop_summary, 'egg')
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
                    value_dict[eggnog_id] = [eggnog_id + ':' + des, int(tmp[2])]
                    insert_data = {
                        "anno_id": anno_id,
                        "eggnog_id": eggnog_id,
                        "des": des,
                        "eff_variant": int(tmp[2]),
                        "total_variant": int(tmp[3])
                    }
                    if pop_summary:
                        if ':'.join([eggnog_id, des]) in gene_dict.keys():
                            trans = gene_dict[':'.join([eggnog_id, des])]
                        else:
                            trans = []
                        insert_data.update({"gene_id": trans})
                    else:
                        trans = tmp[6].split(",")
                        insert_data.update({"gene_id": trans})
                    data_list.append(insert_data)
                if len(data_list) == 0:
                    self.bind_object.logger.info("anno_eggnog注释列表为空！{}".format(eggnog_stat_detail))
                else:
                    self.col_insert_data("sg_region_anno_eggnog_stat", data_list)
                    categories.sort()  # 将X轴按字母顺序排序
                    bar_id = self.sg_bar(task_id, anno_id, 'EggNog分类统计柱形图', "", 3,
                                         location='region_anno_eggnog_class')
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

    def get_genes(self, pop_summary, types='go'):
        """
        获取基因与go以及基因与egg之间对应关系
        GO:0005801,GO:0006904,GO:0009636,GO:0051223,GO:1901998,GO:1902902
        V,W:Defense mechanisms;Extracellular structures
        要处理成V：Defense mechanisms   W：Extracellular structures
        :param pop_summary:
        :param types:
        :return:
        """
        tmp_ = defaultdict(list)
        with open(pop_summary, 'r') as r:
            lines = r.readlines()[1:]
            for line in lines:
                temp = line.strip().split("\t")
                if types == 'go':
                    if temp[17] != '--':
                        for m in temp[17].split(','):
                            tmp_[m].append(temp[1])
                else:
                    if temp[19] != '--':
                        id_ = temp[19].split(':')[0].split(',')
                        des_id_ = temp[19].split(':')[1].split(';')
                        if len(id_) != len(des_id_):
                            continue
                        for i in range(0, len(id_)):
                            tmp_[':'.join([id_[i], des_id_[i]])].append(temp[1])
        return tmp_


if __name__ == "__main__":
    a = RegionAnno(None)
    # project_sn = 'dna_evolution_cui'
    # params = "{\"sg_anno_parmas_id\":\"5adecf95a4e1af4fa5c55351\"}"
    # anno_id = a.add_sg_region_anno(project_sn, task_id, params)
    # summary_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/4.evolution/1.gwas_region/1.data/pop.summary"
    # # summary_path接口task_id查表(tool不支持查表)，传给tool
    # a.add_sg_region_anno_detail(anno_id, summary_path)
    # go_path = '/mnt/ilustre/users/sanger-dev/workspace/20180817/Single_region_anno_0817_220180817/RegionAnno/output/select_region.go.stat.detail'
    # kegg_path = '/mnt/ilustre/users/sanger-dev/workspace/20180817/Single_region_anno_0817_220180817/RegionAnno/output/select_region.kegg.stat.detail'
    # eggnog_path = '/mnt/ilustre/users/sanger-dev/workspace/20180817/Single_region_anno_0817_220180817/RegionAnno/output/select_region.eggnog.stat.detail'
    # pathway_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180817/Single_region_anno_0817_220180817/RegionAnno/output/pathway_dir'
    # a.sg_region_anno_go_stat(task_id, anno_id, go_path)
    # a.sg_region_anno_kegg_stat(task_id, anno_id, kegg_path, pathway_dir)
    # a.sg_region_anno_eggnog_stat(task_id, anno_id, eggnog_path)
    pop_summary_path = "/mnt/lustre/users/sanger/workspace/20181228/Evolution_sanger_149673/Annovar/AnnoCount/output/pop.summary"
    a.add_sg_region_anno_detail("5bab1e5777b3f3b113740aaa", pop_summary_path, "5bab1e5777b3f3b113740aaa")
    # task_id = 'test_liu'
    # anno_id = "5bab1e5777b3f3b113740124"
    # pop_summary = "/mnt/ilustre/users/sanger-dev/workspace/20181112/Evolution_tsg_1112/Annovar/AnnoCount/output/pop.summary"
    # go_stat_detail = "/mnt/ilustre/users/sanger-dev/workspace/20181112/Evolution_tsg_1112/Annovar/AnnoCount/output/pop.go.stat"
    # egg_stat_detail = "/mnt/ilustre/users/sanger-dev/workspace/20181112/Evolution_tsg_1112/Annovar/AnnoCount/output/pop.eggnog.stat"
    # go_summary_path =  "/mnt/ilustre/users/sanger-dev/workspace/20181112/Evolution_tsg_1112/GoSummary/output/pop.2.enrich"
    # a.sg_region_anno_go_stat(task_id, anno_id, go_stat_detail, go_summary_path, pop_summary)
    # a.sg_region_anno_eggnog_stat(task_id, anno_id, egg_stat_detail, pop_summary)
