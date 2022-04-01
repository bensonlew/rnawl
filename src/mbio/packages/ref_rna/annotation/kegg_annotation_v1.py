# -*- coding: utf-8 -*-
import xml.etree.ElementTree as ET
import pymongo
from biocluster.config import Config
import re
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import collections
from itertools import islice
import gridfs
from IPython.display import Image, HTML
from pymongo import MongoClient
import os
import sys


class KeggAnnotation(object):
    def __init__(self):
        """
        设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections
        """
        #self.client = MongoClient(Config().MONGO_BIO_URI)
        self.mongodb = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)]
        self.gene_coll = self.mongodb.kegg_gene
        self.ko_coll = self.mongodb.kegg_ko
        self.png_coll = self.mongodb.kegg_pathway_png
        self.path = collections.defaultdict(str)

    def pathSearch(self, blast_xml, kegg_table, taxonomy):
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
        tablefile.write('#Query\tKO_ID (Gene id)\tKO_name (Gene name)\tHyperlink\tPaths\n')
        docment = ET.parse(blast_xml)
        root = docment.getroot()
        iterns = root.find('BlastOutput_iterations')
        for itern in iterns:
            query = itern.find('Iteration_query-def').text.split()[0]
            iter_hits = itern.find('Iteration_hits')
            hits = iter_hits.findall('Hit')
            if len(hits) > 0:
                hit = hits[0]
                gid = hit.find('Hit_id').text
                # 在数据库中寻找该gene id对应的ko信息，可能改基因并不能在数据库中查到对应的信息
                koids = self.gene_coll.find({"gene_id": gid})
                ko, ko_name, ko_hlink, path = [], [], [], []
                for item in koids:
                    ko.append(item["koid"])
                    result = self.ko_coll.find_one({"ko_id": item["koid"]})
                    if result:
                        ko_name.append(result['ko_name'])
                        pids = result['pathway_id']
                        for index, i in enumerate(pids):
                            if ko_list:
                                if i in ko_list:
                                    self.path[i] = result['pathway_category'][index]  # 对应pathway的definition
                                    path.append('path:' + i)
                            else:
                                self.path[i] = result['pathway_category'][index]
                                path.append('path:' + i)
                    else:
                        print "没有在kegg_ko数据库找到%s" % item["koid"]
                ko = ';'.join(ko)
                ko_name = ';'.join(ko_name)
                ko_hlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + ko
                path = ';'.join(path)
                if not path:
                    path = '\t'
                if ko:
                    # tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\n')
                    if ko_name:
                        tablefile.write(query + '\t' + ko + '\t' + ko_name + '\t' + ko_hlink + '\t' + path + '\n')
                    else:
                        print "没有在kegg_ko数据库找到%s" % item["koid"]
                else:
                    print "%s没有在数据库kegg_gene找到相应的koid" % gid
            else:
                print "没有找到在该query下对应的基因信息！"  #kgml文件中该query没有找到对应的基因
        print "pathSearch finished!"

    def pathTable(self, kegg_table, pathway_path, pidpath):
        """
        根据pathSearch生成的kegg_table.xls统计pathway的信息，输入文件为kegg_table.xls,输出文件为pathway_table.xls,pid.txt
        """
        path_table_xls = open(pathway_path, "wb")  # 输出文件path_table.xls
        pid_txt = open(pidpath, "wb")  # 输出文件pid.txt
        header_line = "Pathway" + "\t" + "Pathway_definition" + "\t" + "num_of_seqs" + "\t" + "seqs_kos/gene_list" + "\t" + "pathway_imagename" + "\n"
        path_table_xls.write(header_line)
        path_table = collections.defaultdict(list)
        kegg_table = islice(open(kegg_table), 1, None)  # 打开kegg_table.xls
        kegg = [i.strip('\n').split('\t') for i in kegg_table]
        table = [(i[0] + '(' + i[1] + ')', i[4]) for i in kegg]
        for i in table:
            for path in i[1].split(';'):
                path_table[path].append(i[0])
                # definition = self.path[path.split(':')[1]]#pid对应的名字
        for key in path_table:
            if key:
                pid = key.split(':')[1]
                definition = self.path[pid][2]
                koids = [i.split('(')[1][0:-1] for i in path_table[key]]
                koid_str = ';'.join(koids)
                pid_txt.write(pid + '\t' + koid_str + '\n')
                try:
                    num_of_seqs = len(path_table[key])
                    if num_of_seqs > 1:
                        genes = ';'.join(path_table[key])
                    else:
                        genes = path_table[key][0]
                    path_image = key.split(":")[1] + '.png'
                    line = key + '\t' + definition + "\t" + str(num_of_seqs) + "\t" + genes + "\t" + path_image
                    path_table_xls.write(line + '\n')
                except Exception as e:
                    print e
            else:
                print "key==None，该基因没有对应的pathway！"
        print "pathTable finished!!!"

    def getPic(self, pidpath, pathwaydir):
        """
        输入文件pid.txt，输出文件夹pathways，作图
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
                    # kgml_id = colln.find_one({"pathway_id": pid})['pathway_ko_kgml']
                    # png_id = colln.find_one({"pathway_id": pid})['pathway_ko_png']
                    result = self.png_coll.find_one({"pathway_id": pid})
                    if result:
                        kgml_id = result['pathway_ko_kgml']
                        png_id = result['pathway_ko_png']
                        k.write(fs.get(kgml_id).read())
                        p.write(fs.get(png_id).read())
                p_kgml = KGML_parser.read(open("pathway.kgml"))
                p_kgml.image = png_path
                for ko in koid:
                    for degree in p_kgml.entries.values():
                        if re.search(ko, degree.name):
                            l.append(degree.id)
                    for n in l:
                        for graphic in p_kgml.entries[n].graphics:
                            graphic.fgcolor = '#CC0000'
                    canvas = KGMLCanvas(p_kgml, import_imagemap=True)
                    canvas.draw(pathwaydir + '/' + pid + '.pdf')
        print "getPic finished!!!"

    def keggLayer(self, pathway_table, layerfile, taxonomyfile):
        """
        输入pathway_table.xls，获取分类信息文件
        """
        f = open(pathway_table)
        d = {}
        ko = {}
        for record in islice(f, 1, None):
            iterm = record.strip('\n').split('\t')
            seqnum = int(iterm[2])
            ko_list = iterm[3]
            pid = iterm[0].split(":")[1]
            result = self.ko_coll.find_one({"pathway_id": {"$in": [pid]}})  # 找到对应的集合
            if result:
                pids = result["pathway_id"]  # 找到对应的pid列表
                for index, i in enumerate(pids):
                    category = result["pathway_category"][index]
                    # [pathway_index]#找到pid对应的层级信息
                    layer_1st = category[0]  # 找到第一层
                    layer_2nd = category[1]  # 找到第二层
                if d.has_key(layer_1st):
                    if d[layer_1st].has_key(layer_2nd):
                        d[layer_1st][layer_2nd] += seqnum
                        ko[layer_1st][layer_2nd] = ko[layer_1st][layer_2nd] + ";" + ko_list
                    else:
                        d[layer_1st][layer_2nd] = seqnum
                        ko[layer_1st][layer_2nd] = ko_list
                else:
                    d[layer_1st] = {}
                    ko[layer_1st] = {}
                    d[layer_1st][layer_2nd] = seqnum
                    ko[layer_1st][layer_2nd] = ko_list
        with open(layerfile, "w+") as k, open(taxonomyfile, "w+") as t:
            for i in d:
                n = 0
                doc = ''
                if i in ko:
                    for j in d[i]:
                        if j in ko[i]:
                            line = i + "\t" + j + "\t" + str(d[i][j]) + "\t" + ko[i][j] + "\n"
                            k.write(line)
                    n += d[i][j]
                    doc += '--' + j + '\t' + str(d[i][j]) + '\n'
                t.write(i + '\t' + str(n) + '\n')
                t.write(doc)

if __name__ == '__main__':
    kegg_anno = KeggAnnotation()
    kegg_anno.pathSearch(sys.argv[1])
    print sys.argv[1]
    kegg_anno.pathTable()
    kegg_anno.getPic()
    kegg_anno.keggLayer()
