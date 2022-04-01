# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis_t import KGMLCanvas
from biocluster.config import Config
import gridfs
import re
import os
from reportlab.lib import colors
import subprocess


class KeggRegulate(object):
    def __init__(self):
        project_type = "ref_rna"
        self.mong_db = Config().get_mongo_client(mtype=project_type, ref=True)[Config().get_mongo_dbname(project_type, ref=True)]

    def get_kgml_and_png(self, pathway_id, kgml_path, png_path):
        collection = self.mong_db['kegg_pathway_png']
        result = collection.find_one({'pathway_id': pathway_id})
        if result:
            kgml_id = result['pathway_ko_kgml']
            png_id = result['pathway_ko_png']
            with open(kgml_path, 'wb') as k, open(png_path, 'wb') as p:
                fs = gridfs.GridFS(self.mong_db)
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
            return True
        else:
            return False

    def get_regulate_table(self, ko_gene, path_ko, regulate_gene, output, gene_list=None):
        """
        生成kegg调控统计表
        ko_gene:ko对应的gene信息:{'ko1': [gene1,gene2], ...,'ko2': [gene1,gene2]}
        path_ko:path对应的ko信息:{'pathway1': [ko1,ko2], ...,'pathway2': [ko1,ko2]}
        output:输出结果的路径
        regulate_dict：gene调控信息:{'up': [gene1,gene2], 'down': [gene1,gene2]}
        """
        colors_w = ['red', 'yellow', 'blue', "green", 'purple', 'pink']
        same_color = "orange"
        new_path_ko = {}
        with open(output, 'wb') as w:
            # w.write('Pathway_id\tKo_ids\tup_numbers\tdown_numbers\tup_genes\tdown_genes\n')
            # modified by qindanhua add 7 line 支持两个以上的基因集统计
            if gene_list:
                genelist_names = gene_list
            else:
                genelist_names = regulate_gene.keys()
            w.write('Pathway_id\tKo_ids\t')
            for gn in genelist_names:
                w.write("{}_numbers\t{}_genes\t".format(gn, gn))
            w.write("\n")

            for path in path_ko:
                ko_ids = set(path_ko[path])
                # up_genes = []
                # down_genes = []
                write_dict = {}
                for gn in genelist_names:
                    write_dict[gn] = []
                for ko in ko_ids:
                    # print ko
                    genes = set(ko_gene[ko])
                    for gn in genelist_names:
                        geneset = set(regulate_gene[gn])
                        same_gene = genes & geneset
                        # print same_genes
                        if len(same_gene) > 0:
                            # print same_gene
                            for sg in same_gene:
                                write_dict[gn].append('{}({})'.format(sg, ko))
                count = 0
                link = 'http://www.genome.jp/kegg-bin/show_pathway?' + path
                # print link
                for gn in genelist_names:
                    count += len(write_dict[gn])
                if count > 0:
                    w.write('{}\t{}\t'.format(path, ";".join(ko_ids)))
                    set_color = {}
                    for n, gn in enumerate(genelist_names):
                        geneko = write_dict[gn]
                        gene_num = len(set([x.split("(")[0] for x in  write_dict[gn]]))
                        w.write("{}\t{}\t".format(gene_num, ";".join(write_dict[gn])))
                        # w.write("{}\t{}\t".format(len(write_dict[gn]), ";".join(write_dict[gn])))
                        for gk in geneko:
                            gko = gk.split('(')[-1][:-1]
                            set_color[gko] = same_color if gko in set_color else colors_w[n]
                    for k in set_color:
                        color_gk = '/' + k + '%09' + set_color[k]
                        link += color_gk
                        # print link
                    w.write('{}'.format(link))
                    w.write("\n")
                    new_path_ko[path] = path_ko[path]
        return new_path_ko

    def get_regulate_table2(self, ko_gene, path_ko, regulate_gene, gene_kegg_gene, output):
        """
        生成kegg调控统计表, 特定物种
        ko_gene:ko对应的gene信息:{'ko1': [gene1,gene2], ...,'ko2': [gene1,gene2]}
        path_ko:path对应的ko信息:{'pathway1': [ko1,ko2], ...,'pathway2': [ko1,ko2]}
        output:输出结果的路径
        regulate_dict：gene调控信息:{'up': [gene1,gene2], 'down': [gene1,gene2]}
        """
        colors_w = ['red', 'yellow', 'blue', "green", 'purple', 'pink']
        same_color = "orange"
        new_path_ko = {}
        with open(output, 'wb') as w:
            # w.write('Pathway_id\tKo_ids\tup_numbers\tdown_numbers\tup_genes\tdown_genes\n')
            # modified by qindanhua add 7 line 支持两个以上的基因集统计
            genelist_names = regulate_gene.keys()
            w.write('Pathway_id\tKo_ids\t')
            for gn in genelist_names:
                w.write("{}_numbers\t{}_genes\t".format(gn, gn))
            w.write("\n")

            for path in path_ko:
                ko_ids = set(path_ko[path])
                # up_genes = []
                # down_genes = []
                write_dict = {}
                for gn in genelist_names:
                    write_dict[gn] = []
                for ko in ko_ids:
                    # print ko
                    genes = set(ko_gene[ko])
                    for gn in genelist_names:
                        geneset = set(regulate_gene[gn])
                        same_gene = genes & geneset
                        # print same_genes
                        if len(same_gene) > 0:
                            # print same_gene
                            for sg in same_gene:
                                kegg_genes = gene_kegg_gene[sg]
                                kegg_genes_rename = kegg_genes[0].replace(',', ';')
                                write_dict[gn].append('{}({})'.format(sg, kegg_genes_rename))
                count = 0
                link = 'http://www.genome.jp/kegg-bin/show_pathway?' + path
                # print link
                for gn in write_dict:
                    count += len(write_dict[gn])
                if count > 0:
                    w.write('{}\t{}\t'.format(path, ";".join(ko_ids)))
                    set_color = {}
                    for n, gn in enumerate(write_dict):
                        geneko = write_dict[gn]
                        gene_num = len(set([x.split("(")[0] for x in  write_dict[gn]]))
                        w.write("{}\t{}\t".format(gene_num, ";".join(write_dict[gn])))
                        # w.write("{}\t{}\t".format(len(write_dict[gn]), ";".join(write_dict[gn])))
                        for gk in geneko:
                            gko = gk.split('(')[-1][:-1]
                            set_color[gko] = same_color if gko in set_color else colors_w[n]
                    for k in set_color:
                        color_gk = '/' + k + '%09' + set_color[k]
                        link += color_gk
                        # print link
                    w.write('{}'.format(link))
                    w.write("\n")
                    new_path_ko[path] = path_ko[path]
        return new_path_ko

    def get_pictrue(self, path_ko, out_dir, regulate_dict=None, image_magick=None):
            """
            传入path_ko统计信息，生成pathway绘图文件夹
            path_ko：path对应的ko信息:{'pathway': [ko1,ko2], ...,'pathway': [ko1,ko2]}
            out_dir:输出结果的目录
            regulate_dict：ko调控信息:{'up': [ko1,ko2], 'up_down': [ko1,ko2],'down': [ko1,ko2]}
            """
            colors_l = ['#CC0000', '#EEEE00', '#388E3C', "#EE6AA7", '#9B30FF', '#7B68EE']
            same_colors = "#FF9800"
            for path in path_ko:
                # path = path.replace("map","ko")
                koid = path_ko[path]
                path = path.replace("map","ko")
                l = {}
                kgml_path = out_dir + '/kgml_tmp.kgml'
                png_path = out_dir + '/{}.png'.format(path)
                result = self.get_kgml_and_png(pathway_id=path, kgml_path=kgml_path, png_path=png_path)
                if result:
                    pathway = KGML_parser.read(open(kgml_path))
                    pathway.image = png_path
                    for ko in koid:
                        for degree in pathway.entries.values():
                            if re.search(ko, degree.name):
                                l[degree.id] = ko
                    for ortholog in pathway.orthologs:
                        for g in ortholog.graphics:
                            g.bgcolor = colors.Color(alpha=0)
                    if not regulate_dict == None:
                        for theid in l:
                            for graphic in pathway.entries[theid].graphics:
                                # modified by qindanhua 20170602 适应基因集的修改，输入的字典名称根据基因集名称变化，不限制于上下调基因
                                same_count = 0
                                for n, gs in enumerate(regulate_dict):
                                    if l[theid] in regulate_dict[gs]:
                                        same_count += 1
                                        graphic.fgcolor = colors_l[n]
                                if same_count > 1:
                                    graphic.fgcolor = same_colors
                    canvas = KGMLCanvas(pathway, import_imagemap=True, label_compounds=True,
                                    label_orthologs=False, label_reaction_entries=False,
                                    label_maps=False, show_maps=False, draw_relations=False, show_orthologs=True,
                                    show_compounds=False, show_genes=False,
                                    show_reaction_entries=False)
                    pdf = out_dir + '/' + path + '.pdf'
                    png = out_dir + '/' + path + '.png'
                    canvas.draw(pdf)
                    os.remove(kgml_path)
                    os.remove(png_path)
                    if image_magick:
                        cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + pdf + ' ' + png
                        try:
                            subprocess.check_output(cmd, shell=True)
                        except subprocess.CalledProcessError:
                            print '图片格式pdf转png出错'


if __name__ == "__main__":
    kegg = KeggRegulate()
    from mbio.files.annotation.kegg.kegg_table import KeggTableFile
    a = KeggTableFile()
    # a.set_path("/mnt/ilustre/users/sanger-dev/workspace/20170802/GenesetKegg_demo_0711_5922_9801/gene_kegg_table.xls")
    pathways = "/mnt/ilustre/users/sanger-dev/workspace/20170727/GenesetKegg_demo_0711_8933_190/output/pathways"
    # a.get_regulate_table
    # ko_genes, path_ko = a.get_pathway_koid()
    # for key in path_ko.keys():
    #     print key + "\t" + ";".join(path_ko[key])
    #     print "***************"
    # kegg.get_pictrue(path_ko=path_ko, out_dir=pathways)
