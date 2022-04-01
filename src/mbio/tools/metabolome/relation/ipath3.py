# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from collections import defaultdict
from bson import ObjectId
import time
import json
from mbio.packages.metabolome.ipath3 import Ipath
from mbio.packages.metabolome.ipath_mark import mark_ipath, id_to_anno_type
from mbio.packages.metabolome.dump_mongo_data import dump_mongo_data, get_db
from mbio.packages.metabolome.get_data import dump_trans_data
from biocluster.config import Config
import pandas as pd


class Ipath3Agent(Agent):
    def __init__(self, parent):
        super(Ipath3Agent, self).__init__(parent)
        options = [
            {"name": "metab_set", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "gene_set", "type": "string"},
            {"name": "metab_annot", "type": "string"},  # 代谢集文件
            {"name": "gene_annot", "type": "string"},  # 代谢集总览表
            {"name": "relation_task_id", "type": "string"},
            {"name": "relation_proj_type", "type": "string"},
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
        super(Ipath3Agent, self).end()


class Ipath3Tool(Tool):
    def __init__(self, config):
        super(Ipath3Tool, self).__init__(config)
        self.ipath_db = self.config.SOFTWARE_DIR + "/database/IPATH3_metab"
        self.ori_ipath_input = self.output_dir + '/ori_ipath_input.xls'
        self.sets = []

    def run(self):
        """
        运行
        :return:
        """
        super(Ipath3Tool, self).run()
        self.get_ipath_marks()
        mark_ipath(self.marks, self.output_dir, self.ipath_db)
        self.svg2png()
        self.logger.info("ipath运行完毕")
        self.end()

    def get_ipath_marks(self):
        '''
        根据基因集、代谢集，和注释信息获取需要在ipath中高亮的信息
        '''
        all_ipath = pd.read_csv(self.ipath_db +'/all_map',sep='\t', header=None)
        has_map = all_ipath[0].tolist()
        to_marks = defaultdict(list)
        metab_set, gene_set = self.get_set()
        metab_anno, gene_anno = self.get_anno()
        anno2metab = defaultdict(set)
        anno2gene = defaultdict(set)
        db_type = ["Antibiotics", "Metabolism", "Microbial_metabolism", "Secondary_metabolites"]
        with open(self.ori_ipath_input, 'wb') as w2, open("set_to_ipath.txt", 'w') as sw:
            color = "#ff0000"
            to_marks[color] = [self.sets[0], set()]
            for one_metab in metab_set:
                name, compounds = metab_anno[one_metab]
                for each in re.split(r";|,", compounds):
                    if each in has_map and (each != '-'):
                        to_marks[color][1].add(each)
                        anno2metab[each].add(one_metab)
                        w2.write("{}\t{}\t#ff0000\tW20\n".format(name, each))

            color = "#00ff00"
            to_marks[color] = [self.sets[1], set()]
            for one_gene in gene_set:
                if one_gene not in gene_anno:
                    continue
                ks = gene_anno[one_gene]
                for each in re.split(r";|,", ks):
                    if each in has_map and (each != '-'):
                        to_marks[color][1].add(each)
                        anno2gene[each].add(one_gene)
                        w2.write("{}\t{}\t#00ff00\tW20\n".format(one_gene, each))
            metab_anno_type = id_to_anno_type(anno2metab, self.ipath_db)
            gene_anno_type = id_to_anno_type(anno2gene, self.ipath_db)
            for t in db_type:
                sw.write("metab_set\t{}\t{}\t{}\n".format(t, self.sets[0], ";".join(metab_anno_type[t])))
                sw.write("gene_set\t{}\t{}\t{}\n".format(t, self.sets[1], ";".join(gene_anno_type[t])))
        self.marks = to_marks

    def svg2png(self):
        rsvg = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/figsave/bin/rsvg-convert")
        for f in os.listdir(self.output_dir):
            if not f.endswith('.svg'):
                continue
            svg = os.path.join(self.output_dir, f)
            png = svg.replace("svg", 'png')
            rsvg_cmd = "{} -f png -o {} {}".format(rsvg, png, svg)
            os.system(rsvg_cmd)

    def get_set(self):
        metab_mongo = get_db()
        metab_set_main = metab_mongo["metab_set"].find_one({"$or": [{"_id": ObjectId(self.option("metab_set"))},
                                                                    {"mai_id": ObjectId(self.option("metab_set"))}]})
        metab_set = metab_mongo["metab_set_detail"].find_one({"set_id": ObjectId(self.option("metab_set"))})["set_list"]
        gene_set, name = dump_trans_data("", self.option("relation_proj_type"), self.option("relation_task_id"),
                                         "geneset", self.option("gene_set"))
        gene_set = gene_set["gene_id"].tolist()
        self.sets = [metab_set_main["name"], name]
        return metab_set, gene_set

    def get_anno(self):
        metab_anno = dump_mongo_data("", self.option("metab_annot"), 'anno_overview')
        metab_info = {}
        metab_anno.agg(lambda x: metab_info.update({x["metab_id"]: [x["metab"], x["compound_id"]]}), axis=1)

        gene_anno, _ = dump_trans_data("", self.option("relation_proj_type"), self.option("relation_task_id"),
                                       "kegg_anno", self.option("gene_annot"))

        gene_info = {}
        gene_anno.agg(lambda x: gene_info.update({x["gene_id"]: x["K_lists"]}), axis=1)

        return metab_info, gene_info

    def set_output(self):
        svg = ['Antibiotics.svg', 'Metabolism.svg','Microbial_metabolism.svg', 'Secondary_metabolites.svg']
        png = ['Antibiotics.png', 'Metabolism.png','Microbial_metabolism.png', 'Secondary_metabolites.png']

        all_files = ['ori_ipath_input.xls']
        all_files.extend(svg)
        all_files.extend(png)

        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

