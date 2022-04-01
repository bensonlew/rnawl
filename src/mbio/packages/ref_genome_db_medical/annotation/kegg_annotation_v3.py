# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
from biocluster.config import Config
import re
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis_t import KGMLCanvas
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from reportlab.lib import colors
import collections
import json
from itertools import islice
import subprocess
import gridfs
import os
import sys
from mbio.packages.rna.annot_config import AnnotConfig

class KeggAnnotation(object):
    def __init__(self, kegg_version="2017"):
        """
        设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections2
        """
        # self.client = Config().biodb_mongo_client
        # self.mongodb = self.client.sanger_biodb
        print "version is {}".format(kegg_version)
        self.client = Config().get_mongo_client(mtype="ref_rna", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("ref_rna", ref=True)]  # 20171101 by zengjing 数据库连接方式修改
        self.gene_coll = self.mongodb.kegg_gene_v1
        self.ko_coll = self.mongodb.kegg_ko_v1
        self.png_coll = self.mongodb.kegg_pathway_png_v1
        self.path = collections.defaultdict(str)
        self.ko2gene = dict()
        self.kegg_version = kegg_version
        print "version is {}".format(kegg_version)

        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=kegg_version)
        self.kegg_json = self.kegg_files_dict["br08901.json"]
        self.kegg_pathway = self.kegg_files_dict["pathway"]
        self.kegg_ko = self.kegg_files_dict["ko_des"]
        self.gene_ko = self.kegg_files_dict["ko_genes.list"]
        '''

        if self.kegg_version == "202003":
            self.kegg_json = Config().SOFTWARE_DIR + "/database/Annotation/other{}/kegg{}/br08901.json".format(kegg_version[0:4], kegg_version)
            self.kegg_pathway = Config().SOFTWARE_DIR + "/database/Annotation/other{}/kegg{}/pathway".format(kegg_version[0:4], kegg_version)
            self.kegg_ko = Config().SOFTWARE_DIR + "/database/Annotation/other{}/kegg{}/ko_des".format(kegg_version[0:4], kegg_version)
            self.gene_ko = Config().SOFTWARE_DIR + "/database/Annotation/other{}/kegg{}/ko_genes.list".format(kegg_version[0:4], kegg_version)
        elif self.kegg_version == "201909":
            self.kegg_json = Config().SOFTWARE_DIR + "/database/Annotation/other2019/kegg201909/br08901.json"
            self.kegg_pathway = Config().SOFTWARE_DIR + "/database/Annotation/other2019/kegg201909/pathway"
            self.kegg_ko = Config().SOFTWARE_DIR + "/database/Annotation/other2019/kegg201909/ko_des"
            self.gene_ko = Config().SOFTWARE_DIR + "/database/Annotation/other2019/kegg201909/ko_genes.list"

        else:
            self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
            self.kegg_pathway = Config().SOFTWARE_DIR + "/database/KEGG/pathway"
            self.kegg_ko = Config().SOFTWARE_DIR + "/database/KEGG/ko_des"
            self.gene_ko = Config().SOFTWARE_DIR + "/database/KEGG/ko_genes.list"

        '''
        self.ko2name = dict()
        self.filter_path = list()

        self.gloabl = ["map01100", "map01110", "map01120", "map01130", "map01200", "map01210", "map01212", "map01230", "map01220"]
        print  "\n".join([self.kegg_json, self.kegg_pathway, self.kegg_ko, self.gene_ko])
        # self.map_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/script/map4.r"
        # self.r_path = "/mnt/ilustre/users/sanger-dev/app/program/R-3.3.3/bin/Rscript"

    def get_filter_path(self, filter_ko):
        filter_list = list()
        with open(filter_ko, 'rb') as ko_f:
            for line in ko_f:
                if line.startswith("ko"):
                    filter_list.append(line.strip().replace("ko", "map"))
                else:
                    filter_list.append(line.strip())
        self.filter_path = list(set(filter_list))

    def get_ko2path(self):
        '''
        获取ko pathway 对应关系
        '''
        ko2path = dict()
        with open(self.kegg_pathway, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                path = cols[0].split(":")[1]
                ko = cols[1].split(":")[1]
                if path in self.gloabl:
                    continue
                if ko in ko2path:
                    if path in ko2path[ko]:
                        pass
                    else:
                        ko2path[ko].add(path)
                else:
                    ko2path[ko] = set([path])
        self.ko2path = ko2path
        return ko2path

    def get_ko2name(self):
        '''
        获取ko 描述对应关系
        '''
        ko2des = dict()
        with open(self.kegg_ko, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                ko = cols[0].split(":")[1]
                if ";" in cols[1]:
                    name = cols[1].split(";")[0]
                else:
                    name = ""
                if ko in ko2des:
                    pass
                else:
                    ko2des[ko] = name
        self.ko2name = ko2des
        return ko2des

    def get_gene2ko(self):
        '''
        获取gene ko对应关系
        '''
        gene2ko = dict()
        with open(self.gene_ko, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                ko = cols[0].split(":")[1]
                gene = cols[1]
                if gene in gene2ko:
                    gene2ko[gene].add(ko)
                else:
                    gene2ko[gene] = set([ko])
        self.gene2ko = gene2ko
        return gene2ko

    def get_kegg_class(self):
        '''
        获取kegg分类注释
        '''
        map2class = dict()
        with open(self.kegg_json, "rb") as f:
            root = json.load(f)
        classI = root['children']
        for classI_child in classI:
            classI_name = classI_child['name']
            for classII_child in classI_child['children']:
                classII_name = classII_child['name']
                for classIII_child in classII_child['children']:
                    classIII_name = classIII_child['name']
                    names = classIII_name.split("  ")
                    path_id = "map" + names[0]
                    name = names[1]
                    map2class[path_id] = [classI_name, classII_name, name]
        self.map2class = map2class
        return map2class



    def pathSearch(self, blast_xml, kegg_table, taxonomy=None):
        # 输入blast比对的xml文件
        """
        输入blast比对的xml文件(Trinity_vs_kegg.xml)，输出kegg_table.xls
        """
        ko2path = self.get_ko2path()
        self.ko2path = ko2path
        gene2ko = self.get_gene2ko()
        ko2name = self.get_ko2name()
        path_db = self.get_keggdb_paths()
        map2class = self.get_kegg_class()
        self.map2class = map2class
        # print map2class
        tablefile = open(kegg_table, "wb")
        ko_list = []
        if taxonomy:
            kofile = open(taxonomy, "rb").readlines()
            for line in kofile:
                line = line.strip()
                ko_list.append(line)
        self.ko_list = ko_list
        tablefile.write('#Query\tKO_ID (Gene id)\tKO_name (Gene name)\tHyperlink\tPaths\tKegg_genes\n')

        if blast_xml == None:
            pass
        else:
            raise Exception("not support ")
        print "pathSearch finished!"

    def merge_known(self, known_ko_annot, kegg_table):
        # 合并已知注释
        print "合并 known_ko_annot {}".format(known_ko_annot)
        def map_kos2paths(kos):
            paths = list()
            for ko in kos.split(";"):
                if ko in self.ko2path:
                    for path in self.ko2path[ko]:
                        if path in self.map2class:
                            if path.replace("map", "ko") in self.ko_list:
                                paths.append(path)
                    # paths.extend(list(self.ko2path[ko]))
            return ";".join([path[3:] for path in set(paths)])
        known_ko = dict()
        annot_ko = dict()
        self.ko2gene = dict()
        with open(known_ko_annot, 'r') as f_known:
            for line in f_known:
                cols = line.strip("\n").split("\t")
                kegg_genes = cols[1]
                if kegg_genes in self.gene2ko:
                    kos = ";".join(self.gene2ko[kegg_genes])
                    ko_name = ";".join([self.ko2name.get(x, "") for x in kos.split(";")])
                    link = "http://www.genome.jp/dbget-bin/www_bget?{}".format(kegg_genes)
                    paths = list()
                    '''
                    for x in kos.split(";"):
                        if x in self.ko2path:
                            paths.extend(list(self.ko2path[x]))
                    # (paths.extend(list(self.ko2path[x])) for x in kos.split(";") if x in self.ko2path)
                    '''
                    # path = ";".join(list(set(paths)))
                    path = map_kos2paths(kos)
                    known_ko[cols[0]] = "\t".join([cols[0], kos, ko_name, link, path, kegg_genes])

        header = ""
        with open(kegg_table, 'r') as f_annot:
            header = f_annot.readline()
            for line in f_annot:
                cols = line.strip().split("\t")
                annot_ko[cols[0]] = line.strip("\n")

        with open(kegg_table, 'w') as w_annot:
            w_annot.write(header)
            for tran_id in list(set(known_ko.keys() + annot_ko.keys())):
                if tran_id in known_ko.keys():
                    ko_annot = known_ko[tran_id].split("\t")
                    ko = ko_annot[1] if len(ko_annot) >= 4 else str()
                    for i in ko.split(";"):
                        if self.ko2gene.has_key(i):
                            self.ko2gene[i] += "|" + tran_id
                        else:
                            self.ko2gene[i] = "accession: " + tran_id

                    ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
                    #result = self.ko_coll.find_one({"ko_id": ko})
                    if ko in self.ko2name:
                        name = self.ko2name[ko]
                    elif ";" in ko:
                        name = ";".join([self.ko2name.get(a_ko, '') for a_ko in ko.split(";")])
                    else:
                        name = ""
                        print "{} do not have name".format(ko)
                    if len(ko_annot) >= 5:
                        if len(self.filter_path) !=0:
                            maps = ";".join(["map" + x for x in ko_annot[4].split(";") if "map" + x not in self.gloabl and "map" + x in self.filter_path])
                        else:
                            maps = ";".join(["map" + x for x in ko_annot[4].split(";") if "map" + x not in self.gloabl])
                    else:
                        maps = ""

                    if len(ko_annot) >= 6:
                        kegg_genes = ko_annot[5]
                    else:
                        kegg_genes = ""
                    w_annot.write("\t".join([
                        tran_id,
                        ko,
                        name,
                        ko_hlink,
                        maps,
                        kegg_genes
                    ]) + "\n")
                else:
                    pass
                '''
                    w_annot.write(annot_ko[tran_id] + "\n")
                    kos = annot_ko[tran_id].split("\t")[1]
                    for i in kos.split(";"):
                        if self.ko2gene.has_key(i):
                            self.ko2gene[i] += "|" + tran_id
                        else:
                            self.ko2gene[i] = "accession: " + tran_id
                '''


    def pathSearch_upload(self, kegg_ids, kegg_table, taxonomy=None):
        # 输入blast比对的xml文件
        """
        输入基因/转录本id对应的K编号文件(kegg.list)，输出kegg_table.xls
        """
        tablefile = open(kegg_table, "wb")
        ko_list = []
        if taxonomy:
            kofile = open(taxonomy, "rb").readlines()
            for line in kofile:
                line = line.strip()
                ko_list.append(line)
        tablefile.write('#Query\tKO_ID (Gene id)\tKO_name (Gene name)\tHyperlink\tPaths\n')
        kegg = open(kegg_ids, "rb").readlines()
        for line in kegg:
            ko, ko_name, ko_hlink, path = [], [], [], []
            line = line.strip().split("\t")
            query = line[0]
            kos = line[1].split(";")
            for ko_id in kos:
                ko.append(ko_id)
                # result = self.ko_coll.find_one({"ko_id": ko_id})
                if 1:
                    ko_name.append(self.ko2name[ko_id])
                    pids = self.ko2path[ko_id]
                    for index, i in enumerate(pids):
                        if ko_list:
                            if i in ko_list:
                                map_id = re.sub("ko", "map", i)
                                class1 = self.map2class[map_id]
                                if map_id not in self.gloabl:
                                    self.path[map_id] = class1[2]  # 对应pathway的definition
                                    if len(self.filter_path) !=0:
                                        if map_id in self.filter_path:
                                            path.append(map_id)
                                    else:
                                        path.append(map_id)
                        else:
                            map_id = re.sub("ko", "map", i)
                            if map_id not in self.gloabl:
                                self.path[map_id] = self.path[map_id] = class1[2]
                                if len(self.filter_path) !=0:
                                    if map_id in self.filter_path:
                                        path.append(map_id)
                                else:
                                    path.append(map_id)
                else:
                    print "没有在kegg_ko数据库找到%s" % ko_id
            ko = ';'.join(ko)
            ko_name = ';'.join(ko_name)
            ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
            path = ';'.join(path)
            if not path:
                path = '\t'
            if ko:
                if ko_name:
                    tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\n')
                else:
                    print "没有在kegg_ko数据库找到%s" % ko_id
        print "pathSearch finished!"

    def get_keggdb_paths(self):
        with open(self.kegg_json, "rb") as f:
            root = json.load(f)
        classI = root['children']
        classII = []
        for i in classI:
            classII.extend(i['children'])
        classIII = []
        for i in classII:
            classIII.extend(i['children'])
        db_paths = ["map" + str(i['name']).split(" ")[0] for i in classIII]
        return db_paths

    def pathTable(self, r_path, map_path, kegg_table, pathway_path, pidpath, link_bgcolor, png_bgcolor, pathwaydir, image_magick):
        """
        根据pathSearch生成的kegg_table.xls统计pathway的信息，输入文件为kegg_table.xls,输出文件为pathway_table.xls,pid.txt
        """
        if not os.path.exists(pathwaydir):
            os.makedirs(pathwaydir)
        path_table_xls = open(pathway_path, "wb")  # 输出文件path_table.xls
        pid_txt = open(pidpath, "wb")  # 输出文件pid.txt
        header_line = "Pathway" + "\t" + "First Category" + "\t" + "Second Category" + "\t" + "Pathway_definition" + "\t" + "num_of_seqs" + "\t" + "seqs_kos/gene_list" + "\t" + "pathway_imagename" + "\t" + "Hyperlink" + "\n"
        path_table_xls.write(header_line)
        path_table = collections.defaultdict(list)
        kegg_table = islice(open(kegg_table), 1, None)  # 打开kegg_table.xls
        kegg = [i.strip('\n').split('\t') for i in kegg_table]
        table = [(i[0] + '(' + i[1] + ')', i[4]) for i in kegg]
        for i in table:
            if i[1] == "":
                continue
            for path in i[1].split(';'):
                if len(self.filter_path) !=0 and path in self.filter_path:
                    path_table[path].append(i[0])
                else:
                    path_table[path].append(i[0])
                print "path is {}".format(path)


        if "/T/" in pathway_path:
            level = "T"
        else:
            level = "G"


        path_db = self.get_keggdb_paths()

        sorted_paths = sorted(path_table.keys(), key=lambda x:path_db.index(x))
        for key in sorted_paths:
            if key:
                pid = key
                # pid = re.sub("map", "ko", key)
                # print "key is" + key
                if key in self.path:
                    definition = self.path[key][2]
                else:
                    definition = ""

                # koids = [i.split('(')[1][0:-1] for i in path_table[key]]
                koids = []
                for i in path_table[key]:
                    for j in i.split('(')[1][0:-1].split(';'):
                        koids.append(j)
                koids = set(koids)
                koid_str = ';'.join(koids)
                ko_color = []
                fgcolor = "NA"
                kos_path = os.path.join(os.getcwd(), key + level + '_' +  "KOs.txt")
                with open(kos_path, "w") as w:
                    w.write("#KO\tbg\tfg\n")
                    for k in koids:
                        ko_color.append(k + "%09" + link_bgcolor)
                        if png_bgcolor.startswith("#"):
                            ko_bgcolor = png_bgcolor
                        else:
                            ko_bgcolor = "#" + png_bgcolor

                        w.write(k + "\t" + ko_bgcolor + "\t" + fgcolor + "\n")
                png_path = pathwaydir + '/' + key + ".png"
                pdf_path = pathwaydir + '/' + key + ".pdf"
                html_path = pathwaydir + '/' + key + ".html"
                self.get_pic(r_path, map_path, key, kos_path, png_path, html_path, png_bgcolor)
                # if image_magick:
                #     cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
                #     try:
                #         subprocess.check_output(cmd, shell=True)
                #     except subprocess.CalledProcessError:
                #         print '图片格式pdf转png出错'
                link = 'http://www.genome.jp/dbget-bin/show_pathway?' + key + '/' + '/'.join(ko_color)
                pid_txt.write(key + '\t' + koid_str + '\n')
                # result = self.ko_coll.find_one({"pathway_id": {"$in": [pid]}})  # 找到对应的集合
                if key in self.map2class:
                    pids = self.map2class[key]  # 找到对应的pid列表
                    layer_1st = pids[0]
                    layer_2nd = pids[1]
                    definition = pids[2]
                    '''
                    for index, i in enumerate(pids):
                        if i == pid:
                            category = result["pathway_category"][index]  # [pathway_index]#找到pid对应的层级信息
                            layer_1st = category[0]  # 找到第一层
                            layer_2nd = category[1]  # 找到第二层
                    '''
                    num_of_seqs = len(path_table[key])
                    geneids = [j.split('(')[0] for j in path_table[key]]
                    genes = ';'.join(geneids)
                    path_image = key + '.png'
                    line = key + '\t' + layer_1st + '\t' + layer_2nd + '\t' + definition + "\t" + str(num_of_seqs) + "\t" + genes + "\t" + path_image + "\t" + link
                    path_table_xls.write(line + '\n')
            else:
                print "key==None，该基因没有对应的pathway！"
        print "pathTable finished!!!"

    def get_pic(self, r_path, map_path, path, kos_path, png_path, html_path, png_bgcolor):
        """
        画通路图
        """
        # fs = gridfs.GridFS(self.mongodb)
        map_id = re.sub("ko", "map", path)
        pid = re.sub("map", "ko", path)
        map_html = KeggHtml(version=self.kegg_version)
        map_html.color_bg[0] = png_bgcolor
        map_html.run(self.html_path + '/' + map_id + '.html', html_path, path + '.png', self.ko2gene)
        ko_list = [set(self.ko2gene.keys())]
        html_mark = html_path + '.mark'
        map_html.run_html_mark(self.html_path + '/' + map_id + '.html', html_mark, path + '.png', self.ko2gene, ko_list)

        '''
        with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
            result = self.png_coll.find_one({"pathway_id": pid})
            if result:
                kgml_id = result['pathway_ko_kgml']
                png_id = result['pathway_map_png']
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
        '''
        cmd = "{} {} {} {} {} {} {}".format(r_path, map_path, path, kos_path, png_path, self.html_path + '/' + map_id + '.kgml', self.html_path + '/' + map_id + '.png')
        try:
            subprocess.check_output(cmd, shell=True)
            print "{} ".format(cmd)
        except subprocess.CalledProcessError:
            print "{} 画图出错".format(cmd)
            # print "{}画图出错".format(path)
            os.system("cp {} {}".format("pathway.png", png_path))

    def keggLayer(self, pathway_table, layerfile):
        """
        输入pathway_table.xls，获取分类信息文件
        """
        f = open(pathway_table)
        d = {}
        ko = {}
        for record in islice(f, 1, None):
            iterm = record.strip('\n').split('\t')
            ko_list = iterm[5].split(";")
            pid = iterm[0] #  = re.sub("map", "ko", iterm[0])
            # result = self.ko_coll.find_one({"pathway_id": {"$in": [pid]}})  # 找到对应的集合
            if pid in self.map2class:
                # pids = result["pathway_id"]  # 找到对应的pid列表
                layer = True
                layer_1st = self.map2class[pid][0]  # 找到第一层
                layer_2nd = self.map2class[pid][1]  # 找到第二层
                '''
                for index, i in enumerate(pids):
                    if i == pid:
                        category = result["pathway_category"][index]  # [pathway_index]#找到pid对应的层级信息
                        layer_1st = category[0]  # 找到第一层
                        layer_2nd = category[1]  # 找到第二层
                        layer = True
                '''
                if layer:
                    if ko.has_key(layer_1st):
                        if ko[layer_1st].has_key(layer_2nd):
                            for k in ko_list:
                                ko[layer_1st][layer_2nd].append(k)
                        else:
                            ko[layer_1st][layer_2nd] = ko_list
                    else:
                        ko[layer_1st] = {}
                        ko[layer_1st][layer_2nd] = ko_list
        # 按官网顺序排序

        with open(layerfile, "w+") as k:
            for i in ["Metabolism","Genetic Information Processing","Environmental Information Processing","Cellular Processes","Organismal Systems", "Human Diseases","Drug Development"]:
                if ko.has_key(i):
                    for j in ko[i]:
                        ko[i][j] = list(set(ko[i][j]))
                        line = i + "\t" + j + "\t" + str(len(ko[i][j])) + "\t" + ';'.join(ko[i][j]) + "\n"
                        k.write(line)
                else:
                    pass


    def getPic(self, pidpath, pathwaydir, image_magick=None):
        """
        输入文件pid.txt，输出文件夹pathways，作图
        image_magick:将pdf转为png的软件目录(/mnt/ilustre/users/sanger-dev/app/program/ImageMagick/bin/convert)
        """
        fs = gridfs.GridFS(self.mongodb)
        f = open(pidpath)
        if not os.path.exists(pathwaydir):
            os.makedirs(pathwaydir)
        for i in f:
            if i:
                i = i.strip('\n').split('\t')
                pid = i[0]
                koid = i[1].split(';')
                l = []
                kgml_path = os.path.join(os.getcwd(), "pathway.kgml")
                png_path = os.path.join(os.getcwd(), "pathway.png")
                if os.path.exists(kgml_path) and os.path.exists(png_path):
                    os.remove(kgml_path)
                    os.remove(png_path)
                with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
                    result = self.png_coll.find_one({"pathway_id": pid})
                    if result:
                        kgml_id = result['pathway_ko_kgml']
                        png_id = result['pathway_ko_png']
                        k.write(fs.get(kgml_id).read())
                        p.write(fs.get(png_id).read())
                p_kgml = KGML_parser.read(open("pathway.kgml"))
                p_kgml.image = png_path
                for ortholog in p_kgml.orthologs:
                    for g in ortholog.graphics:
                        g.bgcolor = colors.Color(alpha=0)
                for ko in koid:
                    for degree in p_kgml.entries.values():
                        if re.search(ko, degree.name):
                            l.append(degree.id)
                    for n in l:
                        for graphic in p_kgml.entries[n].graphics:
                            graphic.fgcolor = '#CC0000'
                canvas = KGMLCanvas(p_kgml, import_imagemap=True, label_compounds=True,
                                    label_orthologs=False, label_reaction_entries=False,
                                    label_maps=False, show_maps=False, draw_relations=False, show_orthologs=True,
                                    show_compounds=False, show_genes=False,
                                    show_reaction_entries=False)
                pdf = pathwaydir + '/' + pid + '.pdf'
                png = pathwaydir + '/' + pid + '.png'
                canvas.draw(pdf)
                if image_magick:
                    cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + pdf + ' ' + png
                    try:
                        subprocess.check_output(cmd, shell=True)
                    except subprocess.CalledProcessError:
                        print '图片格式pdf转png出错'

        print "getPic finished!!!"

    def run(self, r_path, map_path, blast_xml, kegg_ids, kegg_table, pidpath, pathwaydir, pathway_table, layerfile, taxonomy=None, link_bgcolor="green", png_bgcolor= "#00CD00", image_magick=None, html_path="/mnt/ilustre/users/sanger-dev/app/database/KEGG/kegg_2017-05-01/kegg/pathway/map", known_ko=None):
        """blast_xml存在对比对到kegg库的xml文件进行kegg注释统计，kegg_ids存在，对客户上传的kegg注释文件进行kegg注释统计"""

        if not png_bgcolor.startswith("#"):
            png_bgcolor = "#" + png_bgcolor
        if taxonomy != None:
            self.get_filter_path(taxonomy)
        # if blast_xml:
        self.pathSearch(blast_xml, kegg_table, taxonomy)
        if kegg_ids:
            self.pathSearch_upload(kegg_ids, kegg_table, taxonomy)
        # self.getPic(pidpath, pathwaydir, image_magick)
        if known_ko and os.path.isfile(known_ko):
            self.merge_known(known_ko, kegg_table)
        else:
            pass





        self.html_path = html_path
        self.pathTable(r_path, map_path, kegg_table, pathway_table, pidpath, link_bgcolor, png_bgcolor, pathwaydir, image_magick)
        self.keggLayer(pathway_table, layerfile)

if __name__ == '__main__':

    if sys.argv[3] == "None":
        sys.argv[3] = None
    if sys.argv[4] == "None":
        sys.argv[4] = None
    if sys.argv[10] == "None":
        sys.argv[10] = None
    if sys.argv[11] == "None":
        sys.argv[11] = "green"
    if sys.argv[12] == "None":
        sys.argv[12] = "#00CD00"
    if sys.argv[13] == "None":
        sys.argv[13] = None
    if len(sys.argv) > 15 and sys.argv[15] != "None":
        known_ko = sys.argv[15]
    else:
        known_ko = None

    if len(sys.argv) > 16 and sys.argv[16] != "None":
        kegg_version = sys.argv[16]
    else:
        kegg_version = ''

    kegg_anno = KeggAnnotation(kegg_version)

    print sys.argv[0], sys.argv[12]
    kegg_anno.run(r_path=sys.argv[1], map_path=sys.argv[2], blast_xml=sys.argv[3], kegg_ids=sys.argv[4], kegg_table=sys.argv[5],
                  pidpath=sys.argv[6], pathwaydir=sys.argv[7], pathway_table=sys.argv[8],
                  layerfile=sys.argv[9], taxonomy=sys.argv[10], link_bgcolor=sys.argv[11], png_bgcolor=sys.argv[12], image_magick=sys.argv[13], html_path=sys.argv[14], known_ko=known_ko)

    # python kegg_annotation_v2.py /mnt/ilustre/users/sanger-dev/app/program/R-3.3.3/bin/Rscript
    #        /mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/script/map4.r
    #        /mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/taxonomy/mouse/new/anno_stat/blast/gene_kegg.xml
    #        None kegg_table.xls pid.txt pathways pathway_table.xls kegg_layer.xls kegg_taxonomy.xls None None None
