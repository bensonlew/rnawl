# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
import os
import re
import gridfs
import subprocess
import sys
from biocluster.config import Config


class MergeKeggPathway(object):
    def __init__(self):
        self.mongodb = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)]
        # self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.png_coll = self.mongodb.kegg_pathway_png_v1

    def get_pic(self, map_path, r_path, path, kos_path, png_path):
        """
        画通路图
        """
        fs = gridfs.GridFS(self.mongodb)
        pid = re.sub("map", "ko", path)
        with open("pathway.kgml", "w+") as k, open("pathway.png", "w+") as p:
            result = self.png_coll.find_one({"pathway_id": pid})
            if result:
                kgml_id = result['pathway_ko_kgml']
                png_id = result['pathway_map_png']
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
        cmd = "{} {} {} {} {} {} {}".format(r_path, map_path, path, kos_path, png_path, "pathway.kgml", "pathway.png")
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print "{}画图出错".format(path)
            os.system("cp {} {}".format("pathway.png", png_path))

    def merge_kegg_pathway(self, map_path, r_path, image_magick, r_level_path, n_level_path, all_level_path, all_pathways):
        """
        r_level_path: pathway_table.xls(ref)
        n_level_path: pathway_table.xls(new)
        """
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
                    k_ids = k_cols[1].split("%09yellow")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
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
                    k_ids = k_cols[1].split("%09green")
                    k_list = []
                    for k in k_ids:
                        if k.startswith('/'):
                            k_id = k.split('/')[1]
                            k_list.append(k_id)
                    seqlist = item[5].split(";")
                    n_path_list[item[0]] = []
                    n_path_list[item[0]].append(seqlist)
                    n_path_list[item[0]].append(k_list)
                    if item[0] not in path_def:
                        path_def[item[0]] = path
        for map_id in path_def:
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
            for r_c in r_ko:
                r = r_c.split("%09yellow")[0]
                if r_c not in n_ko:
                    r_id = r + '%09' + 'yellow'
                    link.append(r_id)
                    r_kos.append(r)
                else:
                    b_id = r + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(r)
            for n_c in n_ko:
                n = n_c.split("%09green")[0]
                if n_c not in r_ko:
                    n_id = n + '%09' + 'green'
                    link.append(n_id)
                    n_kos.append(n)
                else:
                    b_id = n + '%09' + 'tomato'
                    link.append(b_id)
                    b_kos.append(n)
            link = list(set(link))
            b_kos = list(set(b_kos))
            n_kos = list(set(n_kos))
            r_kos = list(set(r_kos))
            link = 'http://www.genome.jp/kegg-bin/show_pathway?' + map_id + '/' + '/'.join(link)
            pathway_table.write(map_id + "\t" + first_category + "\t" + second_category + "\t" + pathway_definition + "\t" + str(len(seq_list)) + "\t" + ";".join(seq_list) + "\t" + map_id + ".png" + "\t" + link + "\n")
            png_path = all_pathways + '/' + map_id + ".png"
            pdf_path = all_pathways + '/' + map_id + ".pdf"
            kos_path = os.path.join(os.getcwd(), "KOs.txt")
            with open(kos_path, "w") as w:
                w.write("#KO\tbg\tfg\n")
                for k in n_kos:
                    w.write(k + "\t" + "#00CD00" + "\t" + "NA" + "\n")
                for k in r_kos:
                    w.write(k + "\t" + "#FFFF00" + "\t" + "NA" + "\n")
                for k in b_kos:
                    w.write(k + "\t" + "#FFFF00,#00CD00" + "\t" + "NA" + "\n")
            self.get_pic(map_path, r_path, map_id, kos_path, png_path)
            cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print '图片格式pdf转png出错'

if __name__ == "__main__":
    MergeKeggPathway().merge_kegg_pathway(map_path=sys.argv[1], r_path=sys.argv[2], image_magick=sys.argv[3], r_level_path=sys.argv[4], n_level_path=sys.argv[5], all_level_path=sys.argv[6], all_pathways=sys.argv[7])
    # test = MergeKeggPathway()
    # map_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/script/map4.r"
    # r_path = "/mnt/ilustre/users/sanger-dev/app/program/R-3.3.3/bin/Rscript"
    # image_magick = "/mnt/ilustre/users/sanger-dev/app/program/ImageMagick/bin/convert"
    # r_level_path = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/Ensemble_release_36/Annotation/anno_stat/kegg_stat/gene_pathway_table.xls"
    # n_level_path = "/mnt/ilustre/users/sanger-dev/workspace/20170822/NewAnno_arab_test_anno_8/RefAnnotation1/output/anno_stat/kegg_stat/gene_pathway_table.xls"
    # all_level_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/annotation/script/all_pathway_table.xls"
    # all_pathways = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/annotation/script/all_pathways"
    # test.merge_kegg_pathway(map_path, r_path, image_magick, r_level_path, n_level_path, all_level_path, all_pathways)
