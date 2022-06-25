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
from mbio.packages.rna.annot_config import AnnotConfig
import json

class MergeKeggPathway(object):
    def __init__(self, kegg_version=None, kegg_species=None):
        self.kegg_version = kegg_version
        # self.mongodb = Config().get_mongo_client(mtype='ref_rna', ref=True)[Config().get_mongo_dbname('ref_rna', ref=True)]
        # self.png_coll = self.mongodb.kegg_pathway_png_v1
        self.parafly = os.path.join(Config().SOFTWARE_DIR, 'program/parafly-r2013-01-21/bin/bin/ParaFly')
        self.png_modify_cmd = list()
        self.png2pdf = list()
        self.kegg_json = os.path.join(Config().SOFTWARE_DIR, 'database/KEGG/br08901.json')
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=kegg_version)
        self.kegg_json = self.kegg_files_dict["br08901.json"]
        self.gene_ko = self.kegg_files_dict["ko_genes.list"]
        try:
            genome_abr, self.species_path_dir = AnnotConfig().get_kegg_spe_path_dir(version=kegg_version, species=kegg_species)
        except:
            self.set_error("该版本不支持单物种注释")

    def get_keggdb_paths(self):
        with open(self.kegg_json, 'rb') as f:
            root = json.load(f)
        classI = root['children']
        classII = []
        for i in classI:
            classII.extend(i['children'])
        classIII = []
        for i in classII:
            classIII.extend(i['children'])
        db_paths = [str(i['name']).split(' ')[0] for i in classIII]

        return db_paths

    def get_pic(self, map_path, r_path, path, kos_path, png_path, raw_html_path):
        path = "map" + path[-5:]
        spe_html_path = self.species_path_dir
        with open(raw_html_path + '/' + path + '.kgml', 'rb') as kgml:
            if len(kgml.readlines()) == 0:
                cmd = 'cp {}/{}.png {}'.format(raw_html_path, path, png_path)
            else:
                cmd = '{} {} {} {} {} {} {}'.format(
                    r_path,
                    map_path,
                    path,
                    kos_path,
                    png_path,
                    raw_html_path + '/' + path + '.kgml',
                    spe_html_path + '/' + path.replace("map", self.species_abr) + '.png'
                )
        self.png_modify_cmd.append(cmd)

    def get_ko2gene(self, all_gene_ko):
        ko2gene = dict()
        self.spegene2gene = dict()
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
                if len(lines) >= 6 and lines[5] != "":
                    for spe_gene in lines[5].split(";"):
                        if spe_gene in self.spegene2gene:
                            self.spegene2gene[spe_gene].add(lines[0])
                        else:
                            self.spegene2gene[spe_gene] = set([lines[0]])
        return ko2gene

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
                if gene.startswith(self.species_abr):
                    if gene in gene2ko:
                        gene2ko[gene].add(ko)
                    else:
                        gene2ko[gene] = set([ko])
        self.gene2ko = gene2ko
        return gene2ko

    def get_ko2spegene(self):
        '''
        获取ko 和spe 基因对应关系
        '''
        self.spegene2ko = self.get_gene2ko()
        ko2spegene = dict()
        for k, vs in self.spegene2ko.items():
            for v in vs:
                if v in ko2spegene:
                    ko2spegene[v].add(k)
                else:
                    ko2spegene[v] = set([k])
        return ko2spegene
    def get_ko_spegene2gene(self, all_gene_ko):
        '''
        获取单物种ko 基因对应关系
        '''
        # self.ko2gene = self.get_ko2gene(all_gene_ko)
        ko2spegene = self.get_ko2spegene()
        print "ko2spegene is {}".format({k:v for k,v in ko2spegene.items()[:5]})
        print "ko2gene is {}".format({k:v for k,v in self.ko2gene.items()[:5]})
        print "spegene2gene is {}".format({k:v for k,v in self.spegene2gene.items()[:5]})
        ko_spegene2gene = dict()

        for k,g in self.ko2gene.items():
            if k in ko2spegene:
                spe_genes = set(self.spegene2gene.keys()).intersection(ko2spegene[k])
                for spe_gene in spe_genes:
                    if k in ko_spegene2gene:
                        ko_spegene2gene[k] += "\\n" + spe_gene + ":" + ",".join(list(self.spegene2gene[spe_gene]))
                    else:
                        ko_spegene2gene[k] = "\\n" + spe_gene + ":" + ",".join(list(self.spegene2gene[spe_gene]))
        print "ko_spegene2gene is {}".format({k:v for k,v in ko_spegene2gene.items()[:5]})
        return ko_spegene2gene

    def merge_kegg_pathway(self, map_path, r_path, image_magick, r_level_path, n_level_path, all_level_path,
                           all_pathways, all_gene_ko, all_html_path,
                           kegg_species
    ):
        print "merge_kegg_pathway {}".format(locals())
        self.species_abr = kegg_species
        self.ko2gene = self.get_ko2gene(all_gene_ko)
        self.ko_spegene2gene =  self.get_ko_spegene2gene(all_gene_ko)

        prefix = all_gene_ko.split("/")[-1]
        if not os.path.exists(all_pathways):
            os.makedirs(all_pathways)
        path_def = {}
        r_path_list = {}
        n_path_list = {}
        pathway_table = open(all_level_path, 'wb')
        pathway_table.write('Pathway\tFirst Category\tSecond Category\tPathway_definition\tnum_of_seqs\tseqs_kos/gene_list\tpathway_imagename\tHyperlink\n')
        with open(r_level_path, 'rb') as r, open(n_level_path, 'rb') as n:
            lines = r.readlines()
            items = n.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if line[1] == 'Metabolism' and line[2] == 'Global and overview maps':
                    pass
                else:
                    path = line[0] + '|||' + line[1] + '|||' + line[2] + '|||' + line[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + line[0]
                    k_cols = line[7].split(sp)
                    k_ids = k_cols[1].split('/')
                    k_list = []  # 单物种此处为基因
                    ko_list = []
                    for k in k_ids:
                        if k != '':
                            k_id = k.split('%09')[0]
                            k_list.append(k_id)
                            if k_id in self.spegene2ko:
                                ko_list += list(self.spegene2ko[k_id])
                    seqlist = line[5].split(';')
                    path_def[line[0]] = path
                    r_path_list[line[0]] = []
                    r_path_list[line[0]].append(seqlist)
                    r_path_list[line[0]].append(k_list)
                    r_path_list[line[0]].append(ko_list)
            for item in items[1:]:
                item = item.strip().split('\t')
                if item[1] == 'Metabolism' and item[2] == 'Global and overview maps':
                    pass
                else:
                    path = item[0] + '|||' + item[1] + '|||' + item[2] + '|||' + item[3]
                    sp = 'http://www.genome.jp/dbget-bin/show_pathway?' + item[0]
                    k_cols = item[7].split(sp)
                    k_ids = k_cols[1].split('/')
                    k_list = []
                    ko_list = []
                    for k in k_ids:
                        if k != '':
                            k_id = k.split('%09')[0]
                            k_list.append(k_id)
                            if k_id in self.spegene2ko:
                                ko_list += list(self.spegene2ko[k_id])
                    seqlist = item[5].split(';')
                    n_path_list[item[0]] = []
                    n_path_list[item[0]].append(seqlist)
                    n_path_list[item[0]].append(k_list)
                    n_path_list[item[0]].append(ko_list)
                    if item[0] not in path_def:
                        path_def[item[0]] = path
        path_db = self.get_keggdb_paths()
        paths = path_def.keys()
        paths_select = [path for path in paths if path[-5:] in path_db]
        sorted_paths = sorted(paths_select, key=lambda x:path_db.index(x[-5:]))
        for map_id in sorted_paths:
            link = []
            r_kos, n_kos, b_kos = [], [], []
            ref, new, both = [], [], []
            paths = path_def[map_id].split('|||')
            first_category = paths[1]
            second_category = paths[2]
            pathway_definition = paths[3]
            try:
                seq_list = r_path_list[map_id][0]
                r_k = r_path_list[map_id][1]
                r_ko = r_path_list[map_id][2]
            except:
                seq_list = []
                r_k = []
                r_ko = []
            try:
                for q in n_path_list[map_id][0]:
                    seq_list.append(q)
                n_k = n_path_list[map_id][1]
                n_ko = n_path_list[map_id][2]
            except:
                n_k = []
                n_ko = []
            # link
            seq_list = list(set(seq_list))
            all_r_ks = [r_c.split('%09')[0] for r_c in r_k]
            all_n_ks = [n_c.split('%09')[0] for n_c in n_k]
            all_r_ks = set(all_r_ks)
            all_n_ks = set(all_n_ks)
            b_ks = list(all_r_ks & all_n_ks)
            n_ks = list(all_n_ks - all_r_ks)
            r_ks = list(all_r_ks - all_n_ks)
            link = [k + '%09' + 'tomato' for k in b_ks]
            link += [k + '%09' + 'green' for k in n_ks]
            link += [k + '%09' + 'yellow' for k in r_ks]
            link = 'http://www.genome.jp/kegg-bin/show_pathway?' + map_id + '/' + '/'.join(link)

            # png图片
            all_r_kos = [r_c.split('%09')[0] for r_c in r_ko]
            all_n_kos = [n_c.split('%09')[0] for n_c in n_ko]
            all_r_kos = set(all_r_kos)
            all_n_kos = set(all_n_kos)
            b_kos = list(all_r_kos & all_n_kos)
            n_kos = list(all_n_kos - all_r_kos)
            r_kos = list(all_r_kos - all_n_kos)

            pathway_table.write(map_id + '\t' + first_category + '\t' + second_category + '\t' + pathway_definition + '\t' + str(len(seq_list)) + '\t' + ';'.join(seq_list) + '\t' + map_id + '.png' + '\t' + link + '\n')
            png_path = all_pathways + '/' + map_id + '.png'
            pdf_path = all_pathways + '/' + map_id + '.pdf'

            html_path = all_pathways + '/' + map_id + '.html'
            raw_html_path = all_html_path + '/' + map_id.replace(self.species_abr, "map") + '.html'
            map_html = KeggHtml(self.kegg_version)
            map_html.run(raw_html_path, html_path, map_id + '.png', self.ko2gene)
            map_html.run_html_mark(raw_html_path, html_path+ '.mark', map_id + '.png', self.ko_spegene2gene, [all_r_kos, all_n_kos])
            kos_path = all_pathways + '/' + map_id + '.KOs.txt'
            with open(kos_path, 'w') as w:
                w.write('#KO\tbg\tfg\n')
                for k in n_kos:
                    w.write(k + '\t' + '#00CD00' + '\t' + 'NA' + '\n')
                for k in r_kos:
                    w.write(k + '\t' + '#FFFF00' + '\t' + 'NA' + '\n')
                for k in b_kos:
                    w.write(k + '\t' + '#FFFF00,#00CD00' + '\t' + 'NA' + '\n')
            self.get_pic(map_path, r_path, map_id, kos_path, png_path, all_html_path)
            cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + png_path + ' ' + pdf_path
            self.png2pdf.append(cmd)
        with open(prefix + 'cmd_getpic.sh', 'wb') as f:
            f.write('\n'.join(self.png_modify_cmd))
        cmd = '{} -c {} -CPU {} -v -shuffle'.format(self.parafly, prefix +  'cmd_getpic.sh', 10)
        print cmd
        os.system(cmd)
        with open(prefix + 'cmd_png2pdf.sh', 'wb') as f:
            f.write('\n'.join(self.png2pdf))
        cmd = '{} -c {} -CPU {} -v -shuffle'.format(self.parafly, prefix +  'cmd_png2pdf.sh', 10)
        print cmd
        # os.system(cmd)

if __name__ == '__main__':
    if len(sys.argv) <11 or sys.argv[10] == "None":
        MergeKeggPathway().merge_kegg_pathway(
            map_path=sys.argv[1],
            r_path=sys.argv[2],
            image_magick=sys.argv[3],
            r_level_path=sys.argv[4],
            n_level_path=sys.argv[5],
            all_level_path=sys.argv[6],
            all_pathways=sys.argv[7],
            all_gene_ko=sys.argv[8],
            all_html_path=sys.argv[9]
        )
    else:
        MergeKeggPathway(kegg_version=sys.argv[10], kegg_species=sys.argv[11]).merge_kegg_pathway(
                map_path=sys.argv[1],
                r_path=sys.argv[2],
                image_magick=sys.argv[3],
                r_level_path=sys.argv[4],
                n_level_path=sys.argv[5],
                all_level_path=sys.argv[6],
                all_pathways=sys.argv[7],
                all_gene_ko=sys.argv[8],
                all_html_path=sys.argv[9],
                kegg_species=sys.argv[11]
            )
