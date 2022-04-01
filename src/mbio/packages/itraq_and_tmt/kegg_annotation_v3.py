# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
from biocluster.config import Config
import re
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis_t import KGMLCanvas
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from reportlab.lib import colors
import collections
from itertools import islice
import subprocess
import gridfs
import os
import sys
import tarfile
import shutil
import json
import glob
from mbio.packages.rna.annot_config import AnnotConfig


class KeggAnnotation(object):
    def __init__(self, kegg_version="2017", abr=None):
        """
        设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections
        """
        # self.client = Config().biodb_mongo_client
        # self.mongodb = self.client.sanger_biodb
        '''
        self.client = Config().get_mongo_client(mtype="ref_rna", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("ref_rna", ref=True)]  # 20171101 by zengjing 数据库连接方式修改
        self.gene_coll = self.mongodb.kegg_gene_v1
        self.ko_coll = self.mongodb.kegg_ko_v1
        self.png_coll = self.mongodb.kegg_pathway_png_v1
        self.path = collections.defaultdict(str)
        self.gloabl = ["01100", "01110", "01120", "01130", "01200", "01210", "01212", "01230", "01220"]
        '''
        self.path = collections.defaultdict(str)

        self.kegg_version = kegg_version
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=kegg_version)
        self.kegg_json = self.kegg_files_dict["br08901.json"]
        self.kegg_pathway = self.kegg_files_dict["pathway"]
        self.kegg_ko = self.kegg_files_dict["ko_des"]
        self.gene_ko = self.kegg_files_dict["ko_genes.list"]
        # self.gene_path = self.kegg_files_dict["species_path"] + "/{}".format(abr)

        self.gloabl = ["map01100", "map01110", "map01120", "map01130", "map01200", "map01210", "map01212", "map01230", "map01220"]
        self.species_path = []
        self.ko2gene = dict()
        self.species_abr = abr
        self.filter_path = []
        self.spegene2gene = dict()
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

    def pathSearch(self, blast_xml, kegg_table, taxonomy=None, species_abr=None):
        ko2path = AnnotConfig().get_kegg_ko2path(self.kegg_version)
        gene2ko = AnnotConfig().get_kegg_gene2ko(self.kegg_version)
        gene2path = AnnotConfig().get_kegg_gene2path(species_abr, version=self.kegg_version)
        ko2name = AnnotConfig().get_kegg_ko2name(self.kegg_version)
        path_db = AnnotConfig().get_keggdb_paths(self.kegg_version)
        map2class = AnnotConfig().get_kegg_map2class(version=self.kegg_version)
        self.spegene2ko = gene2ko
        self.map2class = map2class

        print "ko2path is {}".format({k:v for k,v in  ko2path.items()[:5]})
        print "ko2name is {}".format({k:v for k,v in  ko2path.items()[:5]})
        print "gene2ko is {}".format({k:v for k,v in  gene2ko.items()[:5]})
        print "gene2path is {}".format({k:v for k,v in  gene2path.items()[:5]})
        print "path_db is {}".format([p for p in path_db[:5]])
        print "map2class is {}".format({k:v for k,v in map2class.items()[:5]})


        # 输入blast比对的xml文件
        """
        输入blast比对的xml文件(Trinity_vs_kegg.xml)，输出kegg_table.xls
        """
        tablefile = open(kegg_table, "wb")
        ko_list = []
        if taxonomy:
            kofile = open(taxonomy, "rb").readlines()
            for line in kofile:
                line = line.strip()
                ko_list.append(line)
        tablefile.write('#Query\tKO_ID (Gene id)\tKO_name (Gene name)\tHyperlink\tPaths\tGene_ID\n')

        docment = ET.parse(blast_xml)
        root = docment.getroot()
        iterns = root.find('BlastOutput_iterations')
        for itern in iterns:
            query = itern.find('Iteration_query-def').text.split()[0]
            iter_hits = itern.find('Iteration_hits')
            hits = iter_hits.findall('Hit')
            if len(hits) > 0:
                mark = 0
                ko, ko_name, ko_hlink, path, gene = [], [], [], [], []
                for hit in hits:
                    mark += 1
                    if mark == 6:
                        break
                    gid = hit.find('Hit_id').text
                    if self.kegg_version >= "2020":
                        gid = species_abr + ":" + gid.split(":")[0]
                    gene.append(gid)
                    if gid in self.spegene2gene:
                        self.spegene2gene[gid] += "," + query
                    else:
                        self.spegene2gene[gid] = query
                    # 在数据库中寻找该gene id对应的ko信息，可能改基因并不能在数据库中查到对应的信息
                    koids = gene2ko.get(gid, set())
                    paths = gene2path.get(gid, set())

                    # koids = self.gene_coll.find({"gene_id": gid})

                    for item in koids:
                        ko.append(item)
                        # result = ko2path[item]
                        if item in ko2path:
                            # print item
                            # ko_name.append(result['ko_name'])
                            pids = paths
                            # print "pids is {}".format(pids)
                            for index, i in enumerate(pids):
                                if 0:
                                    if i in ko_list:
                                        map_id = re.sub(species_abr, "map", i)
                                        if map_id not in self.gloabl:
                                            self.path[map_id] = map2class[map_id]
                                            path.append(map_id)  # 对应pathway的definition
                                else:
                                    print "map_id is {}".format(i)
                                    map_id = re.sub("ko", "map", i)
                                    if map_id not in self.gloabl:
                                        if map_id.replace(self.species_abr, "map") in map2class:
                                            #print "path is {}".format(path_list)
                                            self.path[map_id] = map2class[map_id.replace(self.species_abr, "map")]
                                            if len(self.filter_path) !=0:
                                                if map_id in self.filter_path:
                                                    path.append(i)
                                            else:
                                                path.append(i)

                        else:
                            pass

                # ko_name = ';'.join([ko2name[x] if x in ko2name else "" for x in ko])
                ko_name = ';'.join([ko2name[x] if x in ko2name else "" for x in set(ko)])

                ko = ';'.join(list(set(ko)))
                gene = ';'.join(list(set(gene)))

                for i in ko.split(";"):
                    if ko == "":
                        continue
                    if self.ko2gene.has_key(i):
                        self.ko2gene[i] += "|" + query
                    else:
                        self.ko2gene[i] = "accession: " + query


                ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
                gene_hlink = 'https://www.kegg.jp/dbget-bin/www_bget?{}'.format(gene)
                path = ';'.join(list(set(path)))
                if not path:
                    path = ''
                if gene:
                    if ko_name:
                        tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\t' + gene + '\n')
                    else:
                        tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\t' + gene + '\n')
                        print "没有在kegg_ko数据库找到%s" % ko_name
                else:
                    print "%s没有在数据库kegg_gene找到相应的koid" % gid
            else:
                print "没有找到在该query下对应的基因信息！"  #kgml文件中该query没有找到对应的基因
        print "pathSearch finished!"

    def get_ko2spegene(self):
        '''
        获取ko 和spe 基因对应关系
        '''
        ko2spegene = dict()
        for k, vs in self.spegene2ko.items():
            for v in vs:
                if v in ko2spegene:
                    ko2spegene[v].add(k)
                else:
                    ko2spegene[v] = set([k])
        return ko2spegene
    def get_ko_spegene2gene(self):
        '''
        获取单物种ko 基因对应关系
        '''
        ko2spegene = self.get_ko2spegene()
        ko_spegene2gene = dict()

        for k,g in self.ko2gene.items():
            spe_genes = set(self.spegene2gene.keys()).intersection(ko2spegene[k])
            for spe_gene in spe_genes:
                if k in ko_spegene2gene:
                    ko_spegene2gene[k] += "\\n" + spe_gene + ":" + self.spegene2gene[spe_gene]
                else:
                    ko_spegene2gene[k] = "\\n" + spe_gene + ":" + self.spegene2gene[spe_gene]
        return ko_spegene2gene



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
                result = self.ko_coll.find_one({"ko_id": ko_id})
                if result:
                    ko_name.append(result['ko_name'])
                    pids = result['pathway_id']
                    for index, i in enumerate(pids):
                        if ko_list:
                            if i in ko_list:
                                map_id = re.sub("ko", "map", i)
                                if map_id not in self.gloabl:
                                    self.path[map_id] = result['pathway_category'][index]  # 对应pathway的definition
                                    path.append(map_id)
                        else:
                            map_id = re.sub("ko", "map", i)
                            if map_id not in self.gloabl:
                                self.path[map_id] = result['pathway_category'][index]
                                path.append(map_id)
                else:
                    print "没有在kegg_ko数据库找到%s" % ko_id
            ko = ';'.join(ko)
            ko_name = ';'.join(ko_name)

            ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
            path = ';'.join(path)
            if not path:
                path = ''
            if ko:
                if ko_name:
                    tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\n')
                else:
                    print "没有在kegg_ko数据库找到%s" % ko_id
        print "pathSearch finished!"


    def pathTable(self, r_path, map_path, kegg_table, pathway_path, pidpath, link_bgcolor, png_bgcolor, pathwaydir, image_magick, species_abr, species_path_dir ):
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
        path_table2 = collections.defaultdict(list)
        kegg_table = islice(open(kegg_table), 1, None)  # 打开kegg_table.xls
        kegg = [i.strip('\n').split('\t') for i in kegg_table]
        table = [(i[0] + '(' + i[1] + ')', i[4], i[0] + '(' + i[5] + ')') for i in kegg]
        for i in table:
            for path in i[1].split(';'):
                path_table[path].append(i[0])
                path_table2[path].append(i[2])
        self.ko_spegene2gene =  self.get_ko_spegene2gene()

        for key in path_table:
            if key:
                pid = re.sub(species_abr, "ko", key)
                definition = self.path[key][2]
                # koids = [i.split('(')[1][0:-1] for i in path_table[key]]
                koids = []
                geneids = []
                for i in path_table[key]:
                    for j in i.split('(')[1][0:-1].split(';'):
                        koids.append(j)
                for i in path_table2[key]:
                    for j in i.split('(')[1][0:-1].split(';'):
                        geneids.append(j)
                koids = set(koids)
                koid_str = ';'.join(koids)
                geneids = set(geneids)
                geneid_str = ';'.join(geneids)
                ko_color = []
                gene_color = []
                fgcolor = "NA"
                kos_path = os.path.join(os.getcwd(), "KOs.txt")
                with open(kos_path, "w") as w:
                    w.write("#KO\tbg\tfg\n")
                    for k in koids:
                        ko_color.append(k + "%09" + link_bgcolor)
                        # gene_color.append(k + "%09" + link_bgcolor)
                        w.write(k + "\t" + png_bgcolor + "\t" + fgcolor + "\n")
                    for gene in geneids:
                        gene_color.append(gene + "%09" + link_bgcolor)
                        png_path = pathwaydir + '/' + key + ".png"

                pdf_path = pathwaydir + '/' + key + ".pdf"
                html_path = pathwaydir + '/' + key + ".html"
                self.get_pic(r_path, map_path, key, kos_path, png_path, html_path, species_abr, species_path_dir, png_bgcolor)
                if image_magick:
                    cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
                    try:
                        subprocess.check_output(cmd, shell=True)
                    except subprocess.CalledProcessError:
                        print '图片格式pdf转png出错'
                link = 'http://www.genome.jp/dbget-bin/show_pathway?' + key + '/' + '/'.join(gene_color)
                pid_txt.write(key + '\t' + koid_str + '\t' + geneid_str + '\n')
                # result = self.ko_coll.find_one({"pathway_id": {"$in": [pid]}})  # 找到对应的集合
                if key.replace(species_abr, "map") in self.map2class:
                    pids = self.map2class[key.replace(species_abr, "map")]  # 找到对应的pid列表
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

    def get_pic(self, r_path, map_path, path, kos_path, png_path, html_path, species_abr, species_path_dir, png_bgcolor):
        """
        画通路图
        """
        map_id = re.sub(self.species_abr, "map", path)
        pid = re.sub(self.species_abr, "ko", path)
        species_path_id = path

        map_html = KeggHtml(self.kegg_version)
        map_html.color_bg[0] = png_bgcolor
        # map_html.run(self.html_path + '/' + map_id + '.html', html_path, species_path_id + '.png', self.ko2gene)
        ko_list = [set(self.ko2gene.keys())]
        html_mark = html_path + '.mark'
        ko_spegene2gene =  self.ko_spegene2gene
        map_html.run_html_mark(self.html_path + '/' + map_id + '.html', html_mark, path + '.png', ko_spegene2gene, ko_list)
        print "+++ {} {} {} +++\n".format(self.html_path + '/' + map_id + '.html', html_path, species_path_id + '.png')

        png_dir = os.path.join(species_path_dir, path + '.png')
        os.system("cp {} {}".format(png_dir, "pathway.png"))
        kgml = self.html_path + '/map{}.kgml'.format(path[-5:])
        cmd = "{} {} {} {} {} {} {}".format(r_path, map_path, path, kos_path, png_path, kgml, "pathway.png")
        print "*** " + cmd
        try:
            subprocess.check_output(cmd, shell=True)
            if os.path.exists(png_path):
                pass
            else:
                print "{}画图出错".format(cmd)
                os.system("cp {} {}".format(self.html_path + '/' + map_id + '.png', png_path))
        except subprocess.CalledProcessError:
            print "{}画图出错".format(path)
            os.system("cp {} {}".format("pathway.png", png_path))

    def keggLayer(self, pathway_table, layerfile, species_abr):
        """
        输入pathway_table.xls，获取分类信息文件
        """
        f = open(pathway_table)
        d = {}
        ko = {}
        for record in islice(f, 1, None):
            iterm = record.strip('\n').split('\t')
            ko_list = iterm[5].split(";")
            pid = re.sub(species_abr, "map", iterm[0])
            # result = self.ko_coll.find_one({"pathway_id": {"$in": [pid]}})  # 找到对应的集合
            if pid in self.map2class:
                # pids = result["pathway_id"]  # 找到对应的pid列表
                layer = True
                layer_1st = self.map2class[pid][0]  # 找到第一层
                layer_2nd = self.map2class[pid][1]  # 找到第二层

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

    def get_species_path(self, species_abr, species_path_dir):
        """
        根据物种数据库路径获得该物种相关的通路
        """
        spe_paths = []
        path_file = os.path.join(species_path_dir, species_abr + '.list')
        if os.path.exists(path_file):
            with open (path_file, "r") as f:
                for line in f.readlines():
                    path_all = line.strip().split('\t')[0]
                    spe_paths.append(re.findall(r"\d*$", path_all)[0])
            self.species_path = list(set(spe_paths))
        else:
            files = glob.glob(species_path_dir + "/" +species_abr + "*.png")
            for sfile in files:
                name = os.path.basename(sfile).split(".")[0]
                spe_paths.append(re.findall(r"\d*$", name)[0])
            self.species_path = list(set(spe_paths))




    def run(self, r_path, map_path, blast_xml, kegg_ids, kegg_table, pidpath, pathwaydir, pathway_table, layerfile, taxonomy=None, link_bgcolor="green", png_bgcolor= "#00CD00", image_magick=None, species_abr="map", species_path_dir="/mnt/ilustre/users/sanger-dev/app/database/KEGG/kegg_2017-05-01/kegg/pathway/map", html_path="/mnt/ilustre/users/sanger-dev/app/database/KEGG/kegg_2017-05-01/kegg/pathway/map"):
        """blast_xml存在对比对到kegg库的xml文件进行kegg注释统计，kegg_ids存在，对客户上传的kegg注释文件进行kegg注释统计"""
        print "kegg run params {}".format(locals())
        self.get_species_path(species_abr, species_path_dir)
        self.html_path = html_path
        if blast_xml:
            self.pathSearch(blast_xml, kegg_table, taxonomy, species_abr)
        if kegg_ids:
            self.pathSearch_upload(kegg_ids, kegg_table, taxonomy)
        # self.getPic(pidpath, pathwaydir, image_magick)
        self.pathTable(r_path, map_path, kegg_table,
                       pathway_table, pidpath, link_bgcolor,
                       png_bgcolor, pathwaydir, image_magick,
                       species_abr, species_path_dir)
        self.keggLayer(pathway_table, layerfile, species_abr)


if __name__ == '__main__':
    kegg_anno = KeggAnnotation(sys.argv[-1], abr=sys.argv[14])
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
    print sys.argv[0], sys.argv[12]
    kegg_anno.run(r_path=sys.argv[1], map_path=sys.argv[2], blast_xml=sys.argv[3], kegg_ids=sys.argv[4], kegg_table=sys.argv[5],
                  pidpath=sys.argv[6], pathwaydir=sys.argv[7], pathway_table=sys.argv[8],
                  layerfile=sys.argv[9], taxonomy=sys.argv[10], link_bgcolor=sys.argv[11], png_bgcolor=sys.argv[12], image_magick=sys.argv[13], species_abr=sys.argv[14], species_path_dir=sys.argv[15],html_path=sys.argv[16])

    # python kegg_annotation_v2.py /mnt/ilustre/users/sanger-dev/app/program/R-3.3.3/bin/Rscript
    #        /mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/script/map4.r
    #        /mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/taxonomy/mouse/new/anno_stat/blast/gene_kegg.xml
    #        None kegg_table.xls pid.txt pathways pathway_table.xls kegg_layer.xls kegg_taxonomy.xls None None None
