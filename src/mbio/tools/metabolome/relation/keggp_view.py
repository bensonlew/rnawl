# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from collections import defaultdict
from bson import ObjectId
import time
from mbio.packages.metabolome.kegg_mark import mark_kegg
from mbio.packages.metabolome.dump_mongo_data import dump_mongo_data, get_db
from mbio.packages.metabolome.dump_mongo_data import expdiff_group_samples, get_samples_map
from mbio.packages.metabolome.get_data import dump_trans_data, dump_group_samples
from biocluster.config import Config
import pandas as pd
import json


class KeggpViewAgent(Agent):
    def __init__(self, parent):
        super(KeggpViewAgent, self).__init__(parent)
        options = [
            {"name": "metab_set", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "gene_set", "type": "string"},
            {"name": "metab_annot", "type": "string"},  # 代谢集文件
            {"name": "gene_annot", "type": "string"},  # 代谢集总览表
            {"name": "metab_diff", "type": "string"},
            {"name": "gene_diff", "type": "string"},
            {"name": "relation_task_id", "type": "string"},
            {"name": "relation_proj_type", "type": "string"},
            {"name": "gene_group", "type": "string"},
            {"name": "metab_group", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps("ipath")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.ipath.start()
        self.step.update()

    def stepfinish(self):
        self.step.ipath.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(KeggpViewAgent, self).end()


class KeggpViewTool(Tool):
    def __init__(self, config):
        super(KeggpViewTool, self).__init__(config)
        self.task_id = "_".join(self.config.current_workflow_id.split('_')[:2])
        self.db = get_db()
        self.ref_db = get_db(True)
        relation_info = self.db["sg_relation_analysis"].find_one({"task_id": self.task_id,
                                                                 "relate_task_id": self.option("relation_task_id")})
        if relation_info and "relate_db_version" in relation_info:
            self.relation_db_v = relation_info["relate_db_version"]
        else:
            self.relation_db_v = self.config.DBVersion
        task_info = self.db["sg_task"].find_one({"task_id": self.task_id})
        db_v = ""
        if "database" in task_info:
            db_v = json.loads(task_info["database"])["kegg"]
        if db_v in ["v2021.09.18"]:
            self.kegg_db = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/v2021.09.18/"
            self.html_db = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/v2021.09.18/html"
            self.mongo_org = "kegg_v2021.09.18_organisms"
        elif db_v in ["v94.2"]:
            self.html_db = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/html/"
            self.kegg_db = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/"
            self.mongo_org = "kegg_v94.2_organisms"
        else:
            self.html_db = self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/map_html/"
            self.kegg_db = self.html_path = self.config.SOFTWARE_DIR + "/database/KEGG/"
            self.mongo_org = "kegg_organisms"

        self.path_img = os.path.join(self.output_dir, "pathway_img")
        self.level_out = os.path.join(self.output_dir, "level.txt")
        if not os.path.exists(self.path_img):
            os.mkdir(self.path_img)
        self.cmp2map = self.get_cmp_map()
        self.k2map = self.get_k_map()
        self.map_des = self.get_map_des()
        self.metab_pathways = self.get_metab_pathways()
        self.ori_ipath_input = self.work_dir + '/ori_ipath_input.xls'
        self.sets = []

    def run(self):
        """
        运行
        :return:
        """
        super(KeggpViewTool, self).run()
        self.get_kegg_marks()
        self.logger.info("ipath运行完毕")
        # self.set_output()
        self.end()

    def get_cmp_map(self):
        cmp_db = os.path.join(self.kegg_db, "compound_br08001.txt")
        cmp_map = defaultdict(set)
        with open(cmp_db, 'r') as r:
            for l in r:
                l = l.strip().split('\t')
                l[4] = l[4].replace('map', '')
                cmp_map[l[0]].update(l[4].split(';'))
        return cmp_map

    def get_k_map(self):
        k_db = os.path.join(self.kegg_db, "ko_info.txt")
        k_map = defaultdict(set)
        with open(k_db, 'r') as r:
            for l in r:
                l = l.strip().split('\t')
                l[7] = l[7].replace('ko', '')
                k_map[l[0]].update(l[7].split(';'))
        return k_map

    def get_map_des(self):
        map_db = os.path.join(self.kegg_db, "pathway_des.txt")
        map_des = defaultdict(list)
        with open(map_db, 'r') as r:
            for l in r:
                l = l.strip().split('\t')
                map_des[l[0]] = l[1:]
        return map_des

    def get_metab_pathways(self):
        org_info = self.db["anno_keggp"].find_one({"name" : "AnnoKeggp_Origin", "task_id": self.task_id})
        if not org_info:
            return set()
        orgs = json.loads(org_info["params"])["organism"].split(';')
        maps = set()
        if len(orgs) == 2:
            rets = self.ref_db[self.mongo_org].find_one({"second_category": orgs[1]})
            if rets:
                maps = set(rets["map_list"].split(';'))
        else:
            rets = self.ref_db[self.mongo_org].find_one({"first_category": orgs[0]})
            if rets:
                for one in self.ref_db[self.mongo_org].find({"first_category": orgs[0]}):
                    maps.update(one["map_list"].split(';'))

        return maps

    def fc_change(self):
        '''
        检查比较组对应关系，确定上下调是否需要反转
        '''
        trans_sp2metab_sp = get_samples_map(self.task_id, self.option("relation_task_id"))

        gene_groups = dump_group_samples(self.option("relation_proj_type"),
                                          self.option("relation_task_id"),
                                          self.option("gene_group"),
                                          db_version=self.relation_db_v)
        self.logger.info("gene groups ori: {}".format(gene_groups))
        for g in gene_groups:
            gene_groups[g] = [trans_sp2metab_sp[s] if s in trans_sp2metab_sp else s for s in gene_groups[g]]
        metab_groups = expdiff_group_samples(self.option("metab_diff"))

        metab_gnames = self.option("metab_group").split("_vs_")
        gene_gnames = self.option("gene_group").split("|")

        # 产品确认 代谢中 A_vs_B 和转录中的 B|A 表示的上下调关系一致
        # 两个项目中的分组需要是包含或被包含的关系
        metab_sps = set(metab_groups[metab_gnames[0]])
        gene_sps = set(gene_groups[gene_gnames[1]])
        gene_sps2 = set(gene_groups[gene_gnames[0]])
        self.logger.info("metab samples: {}".format(metab_sps))
        self.logger.info("gene samples: {}".format(gene_sps))

        if metab_sps <= gene_sps or metab_sps > gene_sps:
            fc_change = False
        elif metab_sps <= gene_sps2 or metab_sps > gene_sps2:
            fc_change = True
        else:
            print("metab_groups: {}".format(metab_groups))
            print("gene_groups: {}".format(gene_groups))
            print("trans_sp2metab_sp: {}".format(trans_sp2metab_sp))
            self.set_error("代谢和转录中使用的分组样本需要是包含或被关系, 请重新选项差异分组")
        return fc_change

    def get_kegg_marks(self):
        '''
        根据基因集、代谢集，和注释信息获取需要在kegg html中高亮的信息
        '''
        fc_change = self.fc_change()
        metab_set, gene_set, gene_diff = self.get_set()
        metab_anno, metab_name, gene_anno = self.get_anno()
        gene_name, _ = dump_trans_data("", self.option("relation_proj_type"), self.option("relation_task_id"),
                                       "gene_name", db_version=self.relation_db_v)
        if gene_name is not None:  # 判断是否有gene name
            gene_name = {x[1]["gene_id"]: x[1]["gene_name"] for x in gene_name.iterrows() if x[1]["gene_name"] != '-'}
        metab_diff = dump_mongo_data("", self.option("metab_diff"), "exp_diff")
        metab_diff = metab_diff[metab_diff["metab_id"].isin(metab_set)]
        if 'diff_group' in metab_diff.columns:
            metab_diff = metab_diff[metab_diff["diff_group"] == self.option("metab_group")]

        color_mark = {}  # metab/gene id to color
        map_anno = defaultdict(list)  # {map_id: [gene_id, k, metab_id, compound]}
        to_marks = defaultdict(set)
        def get_color(g_id, fc, cl_dic, gset, anno_dic, reverse=False):
            if g_id not in gset or g_id not in anno_dic:
                return
            up_down = ["#ff0000", "#00ff00"]
            mix = "#0000ff"
            if reverse:
                up_down.reverse()
            for a_id in re.split(r";|,", anno_dic[g_id]):
                try:
                    fc = float(fc)
                    if fc > 1:
                        cl = up_down[0]
                    else:
                        cl = up_down[1]
                except:
                    if fc == "up":
                        cl = up_down[0]
                    else:
                        cl = up_down[1]
                if a_id in cl_dic and cl_dic[a_id] != cl:
                    cl = mix
                cl_dic[a_id] = cl

        metab_diff.agg(lambda x: get_color(x["metab_id"], x["fc"], color_mark, metab_set, metab_anno), axis=1)
        gene_diff.agg(lambda x: get_color(x["gene_id"], x["regulate"], color_mark, gene_set, gene_anno, fc_change), axis=1)

        for one_metab in metab_set:
            cmps = re.split(r";|,", metab_anno[one_metab])
            for cmp in cmps:
                if cmp not in self.cmp2map:
                    continue
                map_ids = self.cmp2map[cmp]
                for map_id in map_ids:
                    if map_id not in map_anno:
                        map_anno[map_id] = [set(), set(), set(), set()]
                    map_anno[map_id][2].add(one_metab)
                    map_anno[map_id][3].add(cmp)

        for one_gene in gene_set:
            if one_gene not in gene_anno:
                continue
            ks = re.split(r";|,", gene_anno[one_gene])
            for k in ks:
                if k not in self.k2map:
                    continue
                map_ids = self.k2map[k]
                for map_id in map_ids:
                    if self.metab_pathways and map_id not in self.metab_pathways:
                        continue
                    if map_id not in map_anno:
                        map_anno[map_id] = [set(), set(), set(), set()]
                    map_anno[map_id][0].add(one_gene)
                    map_anno[map_id][1].add(k)
        url_cl_dic = {"#ff0000": "%09red", "#00ff00": "%09green", "#0000ff": "%09blue"}
        with open(self.level_out, 'w') as w:
            w.write("pathway_id\tdecription\tsecond_category\tfirst_category\tmetab_list\tmetab_name_list\tmetab_count\tgene_list\tgene_name_list\tgene_count\thyperlink\tpdf\n")
            for map_id in map_anno:
                hyperlink = "http://www.genome.jp/dbget-bin/show_pathway?map{}+{}"
                if map_id not in self.map_des:
                    continue
                gene_list, k_list, metab_list, cmp_list = map_anno[map_id]
                if len(gene_list) == 0 or len(metab_list) == 0:  # 产品要求，任意为0及不需相关结果
                    continue
                if gene_name is not None:  #判断是否有基因名称
                    g_name_list = [gene_name[g] if g in gene_name else g for g in gene_list]
                else:
                    g_name_list = gene_list
                metab_name_list = [metab_name[m] if m in metab_name else m for m in metab_list]
                mark_to = {}
                colors = []
                for a in k_list | cmp_list:
                    if a not in color_mark:
                        continue
                    url_cl = url_cl_dic[color_mark[a]]
                    colors.append(a + url_cl)
                    mark_to[a] = color_mark[a]
                h_link = hyperlink.format(map_id, "+".join(colors))
                pdf = mark_kegg(mark_to, map_id, self.path_img, self.html_db)
                pathway_id = "map" + map_id
                out_fmt = "\t".join(["{}"] * 10) + '\n'
                w.write(out_fmt.format(pathway_id, "\t".join(self.map_des[map_id]), ";".join(metab_list), ";".join(metab_name_list),
                                       len(metab_list), ';'.join(gene_list), ';'.join(g_name_list),
                                       len(gene_list), h_link, pdf))

    def get_set(self):
        metab_mongo = self.db
        metab_set_main = metab_mongo["metab_set"].find_one({"$or": [{"_id": ObjectId(self.option("metab_set"))},
                                                                    {"mai_id": ObjectId(self.option("metab_set"))}]})
        metab_set = metab_mongo["metab_set_detail"].find_one({"set_id": ObjectId(self.option("metab_set"))})["set_list"]
        df, name = dump_trans_data("", self.option("relation_proj_type"), self.option("relation_task_id"),
                                         "geneset", self.option("gene_set"),
                                         db_version=self.relation_db_v)
        gene_diff = df
        gene_set = df["gene_id"].tolist()
        self.sets = [metab_set_main["name"], name]
        return metab_set, gene_set, gene_diff

    def get_anno(self):
        metab_anno = dump_mongo_data("", self.option("metab_annot"), 'anno_overview')
        metab_info = {}
        metab_name = {}
        metab_anno.agg(lambda x: metab_info.update({x["metab_id"]: x["compound_id"]}), axis=1)
        metab_anno.agg(lambda x: metab_name.update({x["metab_id"]: x["metab"]}), axis=1)

        gene_anno, _ = dump_trans_data("", self.option("relation_proj_type"), self.option("relation_task_id"),
                                       "kegg_anno", self.option("gene_annot"),
                                       db_version=self.relation_db_v)

        gene_info = {}
        gene_anno.agg(lambda x: gene_info.update({x["gene_id"]: x["K_lists"]}), axis=1)

        return metab_info, metab_name, gene_info
