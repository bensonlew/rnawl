# -*- coding: utf-8 -*-
import sys
from pymongo import MongoClient
from Bio.Blast import NCBIXML
from bson.objectid import ObjectId
import os
import argparse
from biocluster.config import Config
import gridfs
from reportlab.lib import colors
import re
import subprocess
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from IPython.display import Image, HTML
from mbio.packages.metagenomic.kegg_html import KeggHtml
import pandas as pd
import lxml.html


class KeggPathwayImg(object):
    def __init__(self):
        self.client = Config().get_mongo_client(mtype="bacgenome", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("bacgenome", ref=True)]
        self.png_coll = self.mongodb.kegg_pathway_png_v1
        self.fs = gridfs.GridFS(self.mongodb)

    def run_img(self, pathway_file, outdir, png_file, html_db):
        if not os.path.exists(pathway_file):
            raise Exception("{}文件不存在！".format(pathway_file))
        pathwaydir = outdir + "/pathway_img"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(pathwaydir):
            os.mkdir(pathwaydir)
        #enzymes = pd.read_table(enzyme_file, index_col=0).index.tolist() if enzyme_file else []
        enzymes = []
        png_have_color = True  # 旧版宏基因组，细菌基因组等
        if html_db:
            png_have_color = False
        with open(pathway_file, "r") as f1:
            head = f1.next().strip().split("\t")
            for line in f1:
                line = line.strip().split("\t")
                ko = line[5].split(";")[0]
                KOs = line[6]
                path_KO = KOs.split(";")
                l = []
                result = self.png_coll.find_one({"pathway_id": ko})
                if result:
                    kgml_id = result['pathway_ko_kgml']
                    if png_have_color:
                        png_id = result['pathway_ko_png']
                    else:
                        png_id = result["pathway_ko_png"]
                    kgml_path = os.path.join(outdir, "pathway.kgml")
                    png_path = os.path.join(outdir, "pathway.png")
                    filter_path = os.path.join(outdir, "filter.xls")
                    html_path = os.path.join(pathwaydir, ko + '.html')
                    png_base_name = ko
                    with open(kgml_path, "w+") as k, open(png_path, "w+") as p:
                        k.write(self.fs.get(kgml_id).read())
                        p.write(self.fs.get(png_id).read())
                    p_kgml = KGML_parser.read(open(kgml_path))
                    p_kgml.image = png_path
                    for ortholog in p_kgml.orthologs:
                        l = []
                        for g in ortholog.graphics:
                             g.bgcolor = colors.Color(alpha=0)
                    for each in path_KO:
                        for degree in p_kgml.entries.values():
                            if re.search(each, degree.name):
                                l.append(degree.id)
                        for n in l:
                            for graphic in p_kgml.entries[n].graphics:
                                graphic.fgcolor = '#CC0000'
                    canvas = KGMLCanvas(p_kgml, import_imagemap=True, label_compounds=True, label_orthologs=False,
                                        label_reaction_entries=False, label_maps=True, show_maps=False,
                                        draw_relations=True, show_orthologs=True, show_compounds=False, show_genes=True,
                                        show_reaction_entries=False)
                    pdf = pathwaydir + '/' + ko + '.pdf'
                    canvas.draw(pdf)
                    if png_file == 1 or png_file == "1":  # 微生物基因组
                        png = pathwaydir + '/' + ko + '.png'
                        # canvas.draw(pdf)
                        software_dir= Config().SOFTWARE_DIR
                        image_magick = software_dir + "/program/ImageMagick/bin/convert"
                        cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + pdf + ' ' + png
                        try:
                            subprocess.check_output(cmd, shell=True)
                        except subprocess.CalledProcessError:
                            raise Exception('图片格式pdf转png出错')
                    elif png_file:
                        png = pathwaydir + '/' + ko + '.png'
                        db_file = os.path.join(html_db, ko.replace("ko", "map") + '.html')
                        map_html = KeggHtml()
                        map_html.color_bg[0] = "#FFFF00"
                        self.get_html(db_file, html_path, png_base_name)
                        ko_list = [set(path_KO + enzymes)]
                        ko_filted, enzyme_filted = map_html.run_mg_mark(html_path, html_path + ".mark", png_base_name, ko_list)
                        with open("tmp_filted", "w") as f4:
                            f4.write("\n".join(list(ko_list[0])))
                        self.get_pic(ko, "tmp_filted", png, html_db)
                    if enzymes:
                        with open (filter_path, "w") as f3:
                            f3.write("\t".join(ko_filted) + "\n" + "\t".join(enzyme_filted) + "\n")
                else:
                    raise Exception("{} id 找不到pathway图片".format(ko))

    def get_html(self, db_file, html_path, pic_name):
        if not os.path.isfile(db_file):
            print "不存在html图片{}".format(db_file)
            return
        html = lxml.html.parse(db_file)
        root = html.getroot()
        body = root.find('body')
        img = body.find('img')
        if img == None:
            img = body.find_class('image')[0].find('img')
        try:
            img.set('src', pic_name)
        except:
            raise Exception("%s文件没有img" % db_file)
        html.write(html_path)

    # self.get_pic(ko, pathway_file, png, html_db)
    def get_pic(self, ko, ko_file, png, html_db):
        map_id = re.sub("ko", "map", ko)
        r_path = Config().SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        map_path = Config().PACKAGE_DIR + "/metagenomic/pathway_pic.r"
        cmd = "{} {} {} {} {} {} {}".format(r_path, map_path, ko, ko_file, png, html_db + '/' + map_id + ".kgml",
                                html_db + '/' + map_id + '.png')
        print cmd
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "{}画图出错".format(ko)
            os.system("cp {} {}".format("pathway.png", png))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', metavar='[pathway_file]',
                        help='Input pathway_file contain ko_id and KO_list like K0001;K0002;K0003')
    parser.add_argument('-o', metavar='[output dir]', required=True, help='output dir')
    parser.add_argument('-png_file', metavar='[png_file]', required=False, default=None,  help='PNG_file')
    parser.add_argument('-html', metavar='[html_db]', required=False, default=None, help='database html file')
    parser.add_argument('-enzyme_file', metavar='[enzyme_file]', required=False, default=None, help='enzyme profile file')
    args = parser.parse_args()
    pathway_file = args.p
    outdir = args.o
    png_file = args.png_file
    run_img = KeggPathwayImg()
    html_db = args.html
    run_img.run_img(pathway_file, outdir, png_file, html_db)
