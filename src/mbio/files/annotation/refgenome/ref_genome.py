# -*- coding: utf-8 -*-
# __author__ = 'shijin'

import os
from biocluster.iofile import File
import sys
from biocluster.config import Config

class RefGenomeFile(File):
    def __init__(self,path):
        self.workflow_path = path
        self.ref_anno_path = os.path.join(path, "RefAnnotation")
        self.transcript_kegg_path = os.path.join(self.ref_anno_path + "/KeggAnnotation/output/kegg_table.xls")
        self.transcript_go_path = os.path.join(self.ref_anno_path + "/GoAnnotation/output/blast2go.annot")
        self.transcript_cog_path = os.path.join(self.ref_anno_path + "/String2cogv9/cog_list.xls")
        self.gene_kegg_path = os.path.join(self.ref_anno_path + "/RefAnnoStat/kegg_stat/gene_kegg_table.xls")
        self.gene_go_path = os.path.join(self.ref_anno_path + "/RefAnnoStat/go_stat/gene_blast2go.annot")
        self.gene_cog_path = os.path.join(self.ref_anno_path + "/RefAnnoStat/go_stat/gene_cog_list.xls")
        self.path_list = []
        self.path_list.append(self.transcript_cog_path)
        self.path_list.append(self.transcript_go_path)
        self.path_list.append(self.transcript_kegg_path)
        self.path_list.append(self.gene_kegg_path)
        self.path_list.append(self.gene_go_path)
        self.path_list.append(self.gene_cog_path)

    def check(self):
        for path in self.path_list:
            if not os.path.exists(path):
                print "{} doesnot exist".format(path)
                sys.exit(1)

    def get_taxon(self):
        base_name = os.path.basename(self.workflow_path)
        tmp = base_name.split("_")
        base_name = tmp[1] + "_" + tmp[2]
        self.taxon = base_name
        print self.taxon

    def find_ref_genome_path(self):
        self.get_taxon()
        with open(Config().SOFTWARE_DIR + "/database/refGenome/scripts/list_1_1") as file:
            for line in file:
                line = line.strip()
                tmp = line.split("\t")
                taxon = tmp[0]
                if taxon.lower() == self.taxon:
                    print taxon + " has been found"
                    path = tmp[1]
        self.ref_genome_path = path

    def mkdir(self):
        self.ref_anno_path = os.path.split(self.ref_genome_path)[0] + "/anno"
        if not os.path.exists(self.ref_anno_path):
            os.mkdir(self.ref_anno_path)
        else:
            for file in os.listdir(self.ref_anno_path):
                path = os.path.join(self.ref_anno_path, file)
                os.remove(path)
        print self.ref_anno_path


    def move_cog(self):
        pass

    def move_kegg(self):
        t_id_list = []
        ko_id_list = []
        with open(self.transcript_kegg_path, "r") as r:
            r.readline()
            for line in r:
                tmp = line.split("\t")
                t_id = tmp[0]
                ko_id = tmp[1]
                t_id_list.append(t_id)
                ko_id_list.append(ko_id)
        with open(self.ref_anno_path + "/kegg_info", "w") as w:
            for i in range(len(t_id_list)):
                str =  t_id_list[i] + "\t" + "transcript" + "\t" + ko_id_list[i] + "\n"
                w.write(str)
        g_id_list = []
        g_ko_list = []
        with open(self.gene_kegg_path, "r") as r:
            r.readline()
            for line in r:
                tmp = line.split("\t")
                g_id = tmp[0]
                g_ko_id = tmp[1]
                g_id_list.append(g_id)
                g_ko_list.append(g_ko_id)
        with open(self.ref_anno_path + "/kegg_info", "a") as a:
            for i in range(len(g_id_list)):
                str = g_id_list[i] + "\t" + "gene" + "\t"+ g_ko_list[i] + "\n"
                a.write(str)
        print "ko表构建完成"

    def move_go(self):
        t_id_list = []
        go_id_list = []
        with open(self.transcript_go_path, "r") as r:
            # r.readline()
            for line in r:
                tmp = line.split("\t")
                t_id = tmp[0]
                go_id = tmp[1]
                t_id_list.append(t_id)
                go_id_list.append(go_id)
        with open(self.ref_anno_path + "/go_info", "w") as w:
            for i in range(len(t_id_list)):
                str =  t_id_list[i] + "\t" + "transcript" + "\t" + go_id_list[i] + "\n"
                w.write(str)
        g_id_list = []
        g_go_list = []
        with open(self.gene_go_path, "r") as r:
            # r.readline()
            for line in r:
                tmp = line.split("\t")
                g_id = tmp[0]
                g_ko_id = tmp[1]
                g_id_list.append(g_id)
                g_go_list.append(g_ko_id)
        with open(self.ref_anno_path + "/go_info", "a") as a:
            for i in range(len(g_id_list)):
                str =  g_id_list[i] + "\t" + "gene" + "\t" + g_go_list[i] + "\n"
                a.write(str)

    def run(self):  # 默认使用run方法
        self.find_ref_genome_path()
        self.mkdir()
        self.move_kegg()
        self.move_go()

