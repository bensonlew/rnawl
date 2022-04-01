# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from collections import defaultdict
import gridfs
from bson import ObjectId
import time
from mbio.packages.denovo_rna.express.kegg_regulate import KeggRegulate
from biocluster.config import Config
import subprocess
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer


class KeggClassAgent(Agent):
    """
    Kegg分类统计分析，主要用于基因集的重运行步骤
    version v1.0.1
    author: shijin
    last_modify: 2017.8.16
    """
    def __init__(self, parent):
        super(KeggClassAgent, self).__init__(parent)
        options = [
            {"name": "geneset_kegg", "type": "string"},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "kegg_table", "type": "infile", "format": "annotation.kegg.kegg_table"},
            {"name": "geneset_id", "type": "string"},
            {"name": "background_links", "type": "string", "default": ""}  # 底图的地址信息，add_info
        ]
        self.add_option(options)
        self.step.add_steps("kegg_regulate")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kegg_regulate.start()
        self.step.update()

    def stepfinish(self):
        self.step.kegg_regulate.finish()
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
        self._cpu = 11
        self._memory = '10G'

    def end(self):
        super(KeggClassAgent, self).end()


class KeggClassTool(Tool):
    def __init__(self, config):
        super(KeggClassTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python = '/program/Python/bin/'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.map_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/map5.r"
        self.db_path = self.config.SOFTWARE_DIR + "/database/KEGG/xml/"
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.map_dict = {}
        self.geneset_list = list()


    def run(self):
        """
        运行
        :return:
        """
        super(KeggClassTool, self).run()
        # self.get_kegg_pics()
        self.get_kegg_png()
        self.get_dicts()
        self.get_kegg_stat_tmp()
        if self.option("background_links"):
            self.get_background_info()
        self.generate_kegg_stat()
        pathways = self.output_dir + '/pathways'
        if not os.path.exists(pathways):
            os.mkdir(pathways)
        self.generate_ko_txt_dir()
        self.generate_new_pics()
        self.end()

    def get_kegg_pics(self):
        geneset_id = self.option("geneset_id")
        self.mongo_db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        #self.mongo_db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
        fs = gridfs.GridFS(self.mongo_db)
        annotation_collection = self.mongo_db["sg_annotation_kegg"]
        geneset_collection = self.mongo_db["sg_geneset"]
        if "," in geneset_id:
            geneset_id = geneset_id.split(",")[0]
        result = geneset_collection.find_one({"_id": ObjectId(geneset_id)})
        task_id = result["task_id"]
        anno_type = result["type"]
        main_id = annotation_collection.find_one({"task_id":task_id})["_id"]
        kegg_level_collection = self.mongo_db["sg_annotation_kegg_level"]
        results = kegg_level_collection.find({"kegg_id":main_id, "seq_type":"all", "anno_type":anno_type})
        anno_path = self.work_dir + "/png"
        if not os.path.exists(anno_path):
            os.mkdir(anno_path)
        self.logger.info("开始导出kegg png格式图片")
        start = time.time()
        for result in results:
            graph_id = result["graph_png_id"]
            pathway_id = result["pathway_id"]
            with open(anno_path + "/" + pathway_id + ".png", "w") as fw:
                content = fs.get(graph_id).read()
                fw.write(content)
        end = time.time()
        self.logger.info("png格式图片生成完成，耗时{}".format(end - start))

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            os.system('ln -s {} {}'.format(path.rstrip("/"), to_path.rstrip("/")))
        else:
            try:
                transfer = MultiFileTransfer()
                transfer.add_download(path, to_path)
                transfer.perform()
            except:
                self.set_error('file can not find {}'.format(path))
        return to_path


    def get_kegg_png(self):
        geneset_id = self.option("geneset_id")
        mongo_db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        kegg_collection = mongo_db["sg_annotation_kegg"]
        geneset_collection = mongo_db["sg_geneset"]
        if "," in geneset_id:
            geneset_id = geneset_id.split(",")[0]
        result = geneset_collection.find_one({"_id": ObjectId(geneset_id)})
        if not result:
            self.set_error('geneset with main_id: {} was not found'.format(geneset_id))
        task_id = result["task_id"]
        geneset_type = result["type"]
        kegg_main_table = kegg_collection.find_one({"task_id": task_id, 'status': "end"})
        if 'graph_dir' not in kegg_main_table:
            self.get_kegg_pics()
            return
        kegg_annot_output = kegg_main_table['graph_dir']
        if geneset_type == 'transcript':
            # tsanger_27356 > workflow_results > Annotation > TransAnnotation > KEGG >
            kegg_pathways = os.path.join(kegg_annot_output, 'TransAnnotation', 'KEGG', 'alltrans_pathway')
        else:
            kegg_pathways = os.path.join(kegg_annot_output, 'GeneAnnotation', 'KEGG', 'allgene_pathway')
        if not kegg_pathways.endswith("/"):
            kegg_pathways += "/"

        anno_path = self.work_dir + "/png/"

        self.download_s3_file(kegg_pathways, anno_path)

    def get_dicts(self):
        ko_genes, path_ko = self.option('kegg_table').get_pathway_koid()
        geneset_ko = defaultdict(set)
        regulate_gene = {}
        with open(self.option("geneset_kegg"), "r") as f:
            for line in f:
                line = line.strip().split("\t")
                regulate_gene[line[0]] = line[1].split(",")
                geneset_ko[line[0]] = []
                self.geneset_list.append(line[0])
                for key in ko_genes.keys():
                    for gene in regulate_gene[line[0]]:
                        if gene in ko_genes[key]:
                            geneset_ko[line[0]].append(key)
        self.category = geneset_ko  # 基因集与ko的对应关系
        self.ko_genes= ko_genes  # ko与基因的对应关系
        self.path_ko = path_ko  # path小ko与基因大Ko的对应关系
        self.geneset_gene = regulate_gene  # 基因集与基因的对应关系

    def get_kegg_stat_tmp(self):
        self.logger.info("生成kegg_stat_tmp.xls文件")
        KeggRegulate().get_regulate_table(ko_gene=self.ko_genes, path_ko=self.path_ko,
                                          regulate_gene=self.geneset_gene,
                                          output= self.work_dir + '/kegg_stat_tmp.xls',
                                          gene_list = self.geneset_list)
        self.logger.info("生成kegg_stat_tmp.xls文件完毕")

    def generate_ko_txt_dir(self):
        # catgory=self.geneset_gene
        kegg_regulate=self.output_dir + '/kegg_stat.xls'
        out_dir=self.output_dir
        if not os.path.exists(out_dir + "/ko"):
            os.mkdir(out_dir + "/ko")
        f = open(kegg_regulate, "r")
        f.readline()
        for line in f:
            tmp = line.strip().split("\t")
            path = tmp[0]
            gene_num1 = tmp[2]
            gene1_list = []
            gene2_list = []
            if len(self.geneset_gene) == 1:
                if gene_num1 and not gene_num1 == "0" and not gene_num1.startswith("http"):
                    gene1_list_tmp = [x.split("(")[1] if not ")" in x else x.strip(")").split("(")[1] for x in tmp[3].split(");")]
                    for gene in gene1_list_tmp:
                        if gene.find(";") != 1:
                            gene1_list.append(gene)
                        else:
                            gene1_list.extend(gene.split(";"))
            else:
                if gene_num1 and not gene_num1 == "0" and not gene_num1.startswith("http"):
                    gene1_list_tmp = [x.split("(")[1] if not ")" in x else x.strip(")").split("(")[1] for x in tmp[3].split(");")]
                    for gene in gene1_list_tmp:
                        if gene.find(";") != 1:
                            gene1_list.append(gene)
                        else:
                            gene1_list.extend(gene.split(";"))
                gene_num2 = tmp[4]
                if gene_num2 and not gene_num2.startswith("http") and not gene_num2 == "0":
                    gene2_list_tmp = [x.split("(")[1] if not ")" in x else x.strip(")").split("(")[1] for x in tmp[5].split(");")]
                    for gene in gene2_list_tmp:
                        if gene.find(";") != 1:
                            gene2_list.append(gene)
                        else:
                            gene2_list.extend(gene.split(";"))
            gene_list = []
            gene_list.extend(gene1_list)
            gene_list.extend(gene2_list)
            gene_list_unrepeat = list(set(gene_list))
            color_dict = {}
            for gene in gene_list_unrepeat:
                color_dict[gene] = []
                if len(self.geneset_gene) == 1:
                    if gene in gene1_list:
                        color_dict[gene].append("#0000cd")  # 蓝色
                elif len(self.geneset_gene) ==2:
                    if gene in gene1_list and gene in gene2_list:
                        color_dict[gene].append("#0000cd,#ff0000")  # 蓝色,大红
                    elif gene in gene1_list:
                        color_dict[gene].append("#0000cd")  # 蓝色
                    elif gene in gene2_list:
                        color_dict[gene].append("#ff0000")  # 大红
                else:
                    pass
            with open(out_dir + "/ko/" + path, "w") as fw:
                fw.write("#KO\tbg\tfg\n")
                for key in color_dict.keys():
                    if len(color_dict[key]) != 0:
                        str_ = key + "\tNA\t" + ",".join(color_dict[key]) + "\n"
                    else:
                        str_ = key + "\tNA\t" + "NA" + "\n"
                    fw.write(str_)
        self.logger.info("ko文件生成完毕")
        return


    def generate_new_pics(self):
        path_ko=self.path_ko
        kos_path=self.output_dir + "/ko"
        out_dir=self.output_dir
        path_list = path_ko.keys()
        cmd1_list = []
        cmd2_list = []
        for path in path_list:
            ko_path = kos_path + "/" + path
            if os.path.exists(ko_path):
                ko = path.replace("map", "ko")
                cmd = "{} {} {} {} {} {} {}".format(self.r_path, self.map_path, path,
                                                 ko_path, out_dir + "/pathways/" + path + ".png",
                                                 self.db_path + ko + ".xml",
                                                 self.work_dir + "/png/" + path + ".png")
                cmd1_list.append(cmd)
                pdf_path = out_dir + "/pathways/" + path + ".pdf"
                cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                      + out_dir + "/pathways/" + path + ".png" + ' ' + pdf_path
                cmd2_list.append(cmd)
            else:
                db_png_path =  self.work_dir + "/png/" + path + ".png"
                os.link(db_png_path, out_dir + "/pathways/" + path + ".png")
                cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' \
                      + db_png_path + ' ' + out_dir + "/pathways/" + path + ".pdf"
                cmd1_list.append(cmd)
        self.logger.info("开始生成新kegg图片")
        with open(self.work_dir + "/cmd1.list", "w") as fw:
            for i in range(len(cmd1_list)):
                fw.write(cmd1_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd1.list")
        cmd1_obj = self.add_command("cmd1", cmd, ignore_error=True).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("cmd1 list执行成功")
        with open(self.work_dir + "/cmd2.list", "w") as fw:
            for i in range(len(cmd2_list)):
                fw.write(cmd2_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd2.list")
        cmd2_obj = self.add_command("cmd2", cmd, ignore_error=True).run()
        self.wait(cmd2_obj)
        if cmd2_obj.return_code == 0:
            self.logger.info("cmd2 list执行成功")
        self.logger.info("kegg图片生成完毕")


    def get_background_info(self):
        background_links = self.option("background_links")
        if not os.path.exists(background_links):
            self.logger.info("不存在背景信息文件，程序退出")
            return
        with open(background_links) as r:
            r.readline()
            for line in r:
                tmp = line.strip().split("\t")
                pathway = tmp[0]
                links = tmp[1].split("?")[1].split("/")
                links.pop(0)  # 去除掉map
                dict = {}
                for link in links:
                    ko = link.split("%09")[0]  # K06101
                    color = link.split("%09")[1]  # tomato
                    dict[ko] = color
                self.map_dict[pathway] = dict  # self.map_dict["map05340"] = {"K10887": "yellow"}

    def generate_kegg_stat(self):
        with open(self.work_dir + '/kegg_stat_tmp.xls', "r") as file, open(self.output_dir + '/kegg_stat.xls', "w") as fw:
            line = file.readline()
            fw.write(line)
            for line in file:
                tmp = line.strip().split("\t")
                links = tmp[-1]
                pathway = tmp[0]
                links_ko = links.split("?")[1].split("/")
                links_ko.pop(0)  # 去除掉map
                lnk = links.split("?")[0] + "?map=" + pathway + "&multi_query="
                ko_tmp = [x.split("%09")[0] for x in links_ko]
                if pathway in self.map_dict:  # 含有背景色
                    for ko in self.map_dict[pathway].keys():  # 对背景色中的所有项进行循环
                        # self.logger.info(ko)
                        if ko == "":
                            continue
                        if ko in ko_tmp:  # 基因ko显著富集
                            if self.get_color(ko):
                                lnk += ko + "+{},{}%0d%0a".format(self.map_dict[pathway][ko], self.get_color(ko))  # 有背景色,前景色
                            else:
                                lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                        else:
                            lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                else:  # 只标边框
                    for ko in ko_tmp:
                        if ko == "":
                            continue
                        if self.get_color(ko):
                            lnk += ko + "+{},{}%0d%0a".format("white", self.get_color(ko))  # 只有背景色
                        else:
                            pass
                tmp[-1] = lnk
                str_ = "\t".join(tmp) + "\n"
                fw.write(str_)

    def get_color(self, ko):
        """
        根据基因的ko号，获取颜色
        :param ko: 基因的ko号
        :return:
        """
        if len(self.category) == 1:
            if ko in self.category[self.category.keys()[0]]:
                return "blue"
            else:
                return False
        elif len(self.category) == 2:
            lst = self.geneset_list
            # lst = list(self.category.keys())  # 基因集列表
            # lst.sort()
            if ko in self.category[lst[0]] and ko in self.category[lst[1]]:
                return "pink"
            elif ko in self.category[lst[0]]:
                return "blue"
            elif ko in self.category[lst[1]]:
                return "red"
            else:
                return False
