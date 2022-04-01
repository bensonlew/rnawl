# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
import os
import re
import gridfs
import subprocess
import sys
from biocluster.config import Config
from mbio.packages.ref_rna_v2.kegg_html import KeggHtml
from biocluster.config import Config
import json

class MergeKeggPathway(object):
    def __init__(self):
        self.mongodb = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)]
        # self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.png_coll = self.mongodb.kegg_pathway_png_v1
        self.parafly = Config().SOFTWARE_DIR +  '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'
        self.png_modify_cmd = list()
        self.png2pdf = list()
        self.kegg_json = Config().SOFTWARE_DIR + "/database/KEGG/br08901.json"
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

    def get_pic(self, map_path, r_path, path, kos_path, png_path, raw_html_path):

        """
        画通路图
        map_math : R脚本
        path： 通路名如 map00710
        """
        ''' 改为由文件获得图片和kgml
        fs = gridfs.GridFS(self.mongodb)
        pid = re.sub("map", "ko", path)
        with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
            result = self.png_coll.find_one({"pathway_id": pid})
            if result:
                kgml_id = result['pathway_ko_kgml']
                png_id = result['pathway_map_png']
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
        '''
        with open(raw_html_path + "/" + path + ".kgml", "rb") as kgml:
            if len(kgml.readlines()) == 0:
                cmd = "cp {}/{}.png {}".format(raw_html_path, path, png_path)
            else:
                cmd = "{} {} {} {} {} {} {}".format(r_path, map_path, path, kos_path, png_path,
                                                    raw_html_path + "/" + path + ".kgml",
                                                    raw_html_path + "/" + path + ".png")
        self.png_modify_cmd.append(cmd)
        '''
        print cmd

        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "{}画图出错".format(path)
            os.system("cp {} {}".format("pathway.png", png_path))
        '''
    def get_ko2gene(self, all_gene_ko):
        '''
        获取kogene基因转换列表
        '''
        ko2gene = dict()
        with open(all_gene_ko, 'rb') as f:
            f.readline()
            for line in f.readlines():
                lines = line.strip().split("\t")
                for ko in lines[1].split(";"):
                    if ko2gene.has_key(ko):
                        if lines[0] in ko2gene[ko].split(","):
                            pass
                        else:
                            ko2gene[ko] = ko2gene[ko] + "," + lines[0]
                    else:
                        ko2gene[ko] = lines[0]
        return ko2gene

    def merge_kegg_pathway(self, map_path, r_path, image_magick, r_level_path, n_level_path, all_level_path, all_pathways, all_gene_ko, all_html_path):
        """
        r_level_path: pathway_table.xls(ref)
        n_level_path: pathway_table.xls(new)
        """
        self.ko2gene = self.get_ko2gene(all_gene_ko)
        if not os.path.exists(all_pathways):
            os.makedirs(all_pathways)
        path_def = {}
        r_path_list = {}
        n_path_list = {}
        pathway_table = open(all_level_path, "wb")
        pathway_table.write("Pathway\tFirst Category\tSecond Category\tPathway_definition\tnum_of_seqs\tseqs_kos/gene_list\tpathway_imagename\tHyperlink\n")
        with open(r_level_path, "rb") as r, open(n_level_path, "rb") as n:
            lines = r.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                if line[1] == "Metabolism" and line[2] == "Global and overview maps":
                    pass
                else:
                    path = line[0] + "|||" + line[1] + "|||" + line[2] + "|||" + line[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + line[0]
                    k_cols = line[7].split(sp)
                    k_ids = k_cols[1].split("/")
                    k_list = []
                    for k in k_ids:
                        if k != "":
                            k_id = k.split('%09')[0]
                            k_list.append(k_id)
                    seqlist = line[5].split(";")
                    path_def[line[0]] = path
                    r_path_list[line[0]] = []
                    r_path_list[line[0]].append(seqlist)
                    r_path_list[line[0]].append(k_list)
            for item in items[1:]:
                item = item.strip().split("\t")
                if item[1] == "Metabolism" and item[2] == "Global and overview maps":
                    pass
                else:
                    path = item[0] + "|||" + item[1] + "|||" + item[2] + "|||" + item[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + item[0]
                    k_cols = item[7].split(sp)
                    k_ids = k_cols[1].split("/")
                    k_list = []
                    for k in k_ids:
                        if k != "":
                            k_id = k.split('%09')[0]
                            k_list.append(k_id)
                    seqlist = item[5].split(";")
                    n_path_list[item[0]] = []
                    n_path_list[item[0]].append(seqlist)
                    n_path_list[item[0]].append(k_list)
                    if item[0] not in path_def:
                        path_def[item[0]] = path
        path_db = self.get_keggdb_paths()
        sorted_paths = sorted(path_def.keys(), key=lambda x:path_db.index(x))
        for map_id in sorted_paths:
            link = []
            r_kos, n_kos, b_kos = [], [], []
            ref, new, both = [], [], []
            paths = path_def[map_id].split("|||")
            first_category = paths[1]
            second_category = paths[2]
            pathway_definition = paths[3]
            try:
                seq_list = r_path_list[map_id][0]
                r_ko = r_path_list[map_id][1]
            except:
                seq_list = []
                r_ko = []
            try:
                for q in n_path_list[map_id][0]:
                    seq_list.append(q)
                n_ko = n_path_list[map_id][1]
            except:
                n_ko = []
            seq_list = list(set(seq_list))
            all_r_kos = [r_c.split("%09")[0] for r_c in r_ko]
            all_n_kos = [n_c.split("%09")[0] for n_c in n_ko]
            '''
            for r_c in r_ko:
                # r = r_c.split("%09yellow")[0]
                r = r_c.split("%09")[0]
                if r_c not in n_ko:
                    r_id = r + '%09' + 'yellow'
                    link.append(r_id)
                    r_kos.append(r)
                else:
                    b_id = r + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(r)
            for n_c in n_ko:
                # n = n_c.split("%09green")[0]
                n = n_c.split("%09")[0]
                if n not in r_kos:
                    n_id = n + '%09' + 'green'
                    link.append(n_id)
                    n_kos.append(n)
                else:
                    b_id = n + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(n)
            '''
            all_r_kos = set(all_r_kos)
            all_n_kos = set(all_n_kos)
            b_kos = list(all_r_kos & all_n_kos)
            n_kos = list(all_n_kos - all_r_kos)
            r_kos = list(all_r_kos - all_n_kos)
            link = [ko + '%09' + 'tomoto' for ko in b_kos]
            link += [ko + '%09' + 'green' for ko in n_kos]
            link += [ko + '%09' + 'yellow' for ko in r_kos]

            # b_kos = list(set(b_kos))
            # n_kos = list(set(n_kos))
            # r_kos = list(set(r_kos))
            link = 'http://www.genome.jp/kegg-bin/show_pathway?' + map_id + '/' + '/'.join(link)
            pathway_table.write(map_id + "\t" + first_category + "\t" + second_category + "\t" + pathway_definition + "\t" + str(len(seq_list)) + "\t" + ";".join(seq_list) + "\t" + map_id + ".png" + "\t" + link + "\n")
            png_path = all_pathways + '/' + map_id + ".png"
            pdf_path = all_pathways + '/' + map_id + ".pdf"
            html_path = all_pathways + '/' + map_id + ".html"
            raw_html_path = all_html_path + "/" + map_id + ".html"
            map_html = KeggHtml()
            map_html.run(raw_html_path, html_path, map_id + '.png', self.ko2gene)
            kos_path = all_pathways + '/' + map_id + ".KOs.txt"
            with open(kos_path, "w") as w:
                w.write("#KO\tbg\tfg\n")
                for k in n_kos:
                    w.write(k + "\t" + "#00CD00" + "\t" + "NA" + "\n")
                for k in r_kos:
                    w.write(k + "\t" + "#FFFF00" + "\t" + "NA" + "\n")
                for k in b_kos:
                    w.write(k + "\t" + "#FFFF00,#00CD00" + "\t" + "NA" + "\n")
            self.get_pic(map_path, r_path, map_id, kos_path, png_path, all_html_path)
            cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
            self.png2pdf.append(cmd)

        # 将生成图片串行改为并行
        with open("cmd_getpic.sh", "wb") as f:
            f.write("\n".join(self.png_modify_cmd))

        cmd = "{} -c {} -CPU {} -v -shuffle".format(self.parafly, "cmd_getpic.sh", 10)
        print cmd
        os.system(cmd)

        with open("cmd_png2pdf.sh", "wb") as f:
            f.write("\n".join(self.png2pdf))
        cmd = "{} -c {} -CPU {} -v -shuffle".format(self.parafly, "cmd_png2pdf.sh", 10)
        print cmd
        os.system(cmd)
        # print cmd
        # try:
        #     subprocess.check_output(cmd, shell=True)
        # except subprocess.CalledProcessError:
        #     print '图片格式pdf转png出错'

if __name__ == "__main__":
    MergeKeggPathway().merge_kegg_pathway(map_path=sys.argv[1], r_path=sys.argv[2], image_magick=sys.argv[3], r_level_path=sys.argv[4], n_level_path=sys.argv[5], all_level_path=sys.argv[6], all_pathways=sys.argv[7], all_gene_ko=sys.argv[8], all_html_path=sys.argv[9])
    # test = MergeKeggPathway()
    # map_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/script/map4.r"
    # r_path = "/mnt/ilustre/users/sanger-dev/app/program/R-3.3.3/bin/Rscript"
    # image_magick = "/mnt/ilustre/users/sanger-dev/app/program/ImageMagick/bin/convert"
    # r_level_path = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/Ensemble_release_36/Annotation/anno_stat/kegg_stat/gene_pathway_table.xls"
    # n_level_path = "/mnt/ilustre/users/sanger-dev/workspace/20170822/NewAnno_arab_test_anno_8/RefAnnotation1/output/anno_stat/kegg_stat/gene_pathway_table.xls"
    # all_level_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/annotation/script/all_pathway_table.xls"
    # all_pathways = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/annotation/script/all_pathways"
    # test.merge_kegg_pathway(map_path, r_path, image_magick, r_level_path, n_level_path, all_level_path, all_pathways)
