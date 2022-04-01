# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.api.database.base import Base
import re
import gridfs
from reportlab.lib import colors
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import shutil
import subprocess
import lxml.html
from biocluster.config import Config
from mbio.packages.metabolome.anno_keggp import KeggHtml
# 参数需要传入anno_overview的id，通过to_file生成anno_overview_ko表

class Metabsetp(Base):
    def __init__(self, image_magick, html_path):
        super(Metabsetp, self).__init__()
        self._project_type = 'metabolome'
        self.ko_data = ""
        self.seta_ko_dic = {}
        self.setb_ko_dic = {}
        self.image_magick = image_magick
        self.html_path = html_path

    def one_set_run(self, set_list, stat_path):
        # to_file table format ["metab_id", "pathway_id", "description", "first_category", "second_category", "compound_id", "hyperlink"]
        category_dic = {}
        ko_dic = {}
        for indexs in self.ko_data.index:
            metab_id = self.ko_data['metab_id'][indexs]
            fc = self.ko_data['first_category'][indexs]
            sc = self.ko_data['second_category'][indexs]
            # path_id = self.ko_data['pathway_id'][indexs]
            path_id = indexs
            metab_list = metab_id.split(";")
            fil_metab = set(metab_list) & set(set_list)
            if len(fil_metab) != 0:
                ko_dic[path_id] = list(fil_metab)
                if category_dic.has_key(fc):
                    if category_dic[fc].has_key(sc):
                        category_dic[fc][sc].append(path_id)
                    else:
                        category_dic[fc][sc] = [path_id]
                else:
                    category_dic[fc] = {sc: [path_id]}
        self.set_stat_output(category_dic, ko_dic, stat_path)
        return ko_dic

    def set_stat_output(self, category_dic, ko_dic, stat_path):
        stat_file = open(stat_path, "w")
        stat_file.write("first_category\tsecond_category\tmetab_id\tcount\n")
        for fc in sorted(category_dic.keys()):
            for sc in sorted(category_dic[fc].keys()):
                if sc == "Global and overview maps":
                    continue  # 不对此二级分类作图 modified by GHD @ 20181205
                ko_list = category_dic[fc][sc]
                metab_set = set()
                for ko in ko_list:
                    metab_set = metab_set | set(ko_dic[ko])
                metab_list = list(metab_set)
                metab_id = ';'.join(metab_list)
                count = len(metab_list)
                stat_file.write("%s\t%s\t%s\t%s\n" % (fc,sc,metab_id,count))
        stat_file.close()

    def set_level_output(self, level_path, outdir, version=2):
        # hyperlink_list = []
        if self.setb_ko_dic == {}:
            level_data = self.ko_data[self.ko_data.index.isin(self.seta_ko_dic.keys())]
            # level_data = self.ko_data[self.seta_ko_dic.keys()]
            level_data["count"] = 0
            for indexs in level_data.index:
                hyperlink = "http://www.genome.jp/dbget-bin/show_pathway?" + indexs
                seta_ko_list = self.seta_ko_dic[indexs]
                metab_list, compound_list = self.cal_com(seta_ko_list, level_data['metab_id'][indexs], level_data['compound_id'][indexs])
                # level_data['metab_id'][indexs] = ";".join(metab_list)
                # level_data['compound_id'][indexs] = ";".join(compound_list)
                level_data.at[indexs, 'metab_id'] = ";".join(metab_list)
                level_data.at[indexs, 'compound_id'] = ';'.join(compound_list)
                level_data.at[indexs, "count"] = len(metab_list)
                for c in compound_list:
                    hyperlink += "+" + c + "%09red"   #green
                # hyperlink_list.append(hyperlink)
                # level_data['hyperlink'][indexs] = hyperlink
                level_data.at[indexs, 'hyperlink'] = hyperlink
                self.run_img(indexs, compound_list, outdir)

        # version 2.0
        elif self.setb_ko_dic != {}  and  version==2:
            bingji = set(self.seta_ko_dic.keys()) | set(self.setb_ko_dic.keys())
            level_data = self.ko_data[self.ko_data.index.isin(list(bingji))]
            level_data["count"] = 0
            level_data["compound_id_b"] = '-'
            level_data["metab_id_b"] = '-'
            level_data["count_b"] = 0  ##上述程序顺序不要调整
            for indexs in level_data.index:
                hyperlink = "http://www.genome.jp/dbget-bin/show_pathway?" + indexs
                if indexs not in self.seta_ko_dic.keys():
                    seta_ko_set = set([])
                else:
                    seta_ko_set = set(self.seta_ko_dic[indexs])

                if indexs not in self.setb_ko_dic.keys():
                    setb_ko_set = set([])
                else:
                    setb_ko_set = set(self.setb_ko_dic[indexs])


                # jiaoji_set = seta_ko_set & setb_ko_set
                # in_a_set = seta_ko_set - setb_ko_set
                # in_b_set = setb_ko_set - seta_ko_set
                #metab_list,compound_list = self.cal_com(jiaoji_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])

                #metab_list_a,compound_list_a = self.cal_com(seta_ko_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])
                #metab_list_b,compound_list_b = self.cal_com(setb_ko_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])

                #metab_in_a,compound_in_a = self.cal_com(in_a_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])
                #metab_in_b,compound_in_b = self.cal_com(in_b_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])

                #20200717 有的代谢物不一样，但compund id是一样的，所以要不能用上面的方法，根据代谢物id来判断交集
                metab_list_a,compound_list_a = self.cal_com(seta_ko_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])
                metab_list_b,compound_list_b = self.cal_com(setb_ko_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])

                compound_list = list(set(compound_list_a) & set(compound_list_b))
                compound_in_a = []
                #metab_in_a = []
                compound_in_b = []
                #metab_in_b = []

                for cid in range(len(compound_list_a)):
                    compound = compound_list_a[cid]
                    #meta =  metab_list_a[cid]
                    if compound not in compound_list:
                        compound_in_a.append(compound)
                        #metab_in_a.append(meta)

                for cid in range(len(compound_list_b)):
                    compound = compound_list_b[cid]
                    #meta =  metab_list_b[cid]
                    if compound not in compound_list:
                        compound_in_b.append(compound)
                        #metab_in_b.append(meta)




                for c in compound_list:
                    hyperlink += "+" + c + "%09blue"  #
                for c in compound_in_a:
                    hyperlink += "+" + c + "%09red"
                for c in compound_in_b:
                    hyperlink += "+" + c + "%09green"
                # level_data['metab_id'][indexs] = ";".join(metab_list)
                # level_data['compound_id'][indexs] = ";".join(compound_list)
                # level_data['hyperlink'][indexs] = hyperlink
                level_data.at[indexs, 'metab_id'] = ";".join(metab_list_a)
                level_data.at[indexs, 'compound_id'] = ";".join(compound_list_a)
                level_data.at[indexs, 'hyperlink'] = hyperlink
                level_data.at[indexs, "count"] = len(metab_list_a)

                level_data.at[indexs, 'metab_id_b'] = ";".join(metab_list_b)
                level_data.at[indexs, 'compound_id_b'] = ";".join(compound_list_b)
                level_data.at[indexs, "count_b"] = len(metab_list_b)
                print 'jioaji:'
                print compound_list
                print 'a:'
                print compound_in_a
                print 'b：'
                print compound_in_b

                self.run_img(indexs, compound_in_a, outdir,compound_list , compound_in_b)



        else:   #交集  v1版
            jiaoji = set(self.seta_ko_dic.keys()) & set(self.setb_ko_dic.keys())
            level_data = self.ko_data[self.ko_data.index.isin(list(jiaoji))]
            for indexs in level_data.index:
                hyperlink = "http://www.genome.jp/dbget-bin/show_pathway?" + indexs
                seta_ko_set = set(self.seta_ko_dic[indexs])
                setb_ko_set = set(self.setb_ko_dic[indexs])
                jiaoji_set = seta_ko_set & setb_ko_set
                in_a_set = seta_ko_set - setb_ko_set
                in_b_set = setb_ko_set - seta_ko_set
                metab_list,compound_list = self.cal_com(jiaoji_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])
                metab_in_a,compound_in_a = self.cal_com(in_a_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])
                metab_in_b,compound_in_b = self.cal_com(in_b_set, level_data['metab_id'][indexs], level_data['compound_id'][indexs])
                for c in compound_list:
                    hyperlink += "+" + c + "%09green"
                for c in compound_in_a:
                    hyperlink += "+" + c + "%09orange"
                for c in compound_in_b:
                    hyperlink += "+" + c + "%09blue"
                # level_data['metab_id'][indexs] = ";".join(metab_list)
                # level_data['compound_id'][indexs] = ";".join(compound_list)
                # level_data['hyperlink'][indexs] = hyperlink
                level_data.at[indexs, 'metab_id'] = ";".join(metab_list)
                level_data.at[indexs, 'compound_id'] = ";".join(compound_list)
                level_data.at[indexs, 'hyperlink'] = hyperlink
                level_data.at[indexs, "count"] = len(metab_list)

                self.run_img(indexs, compound_in_a, outdir,compound_list , compound_in_b)
                # self.run_img(indexs, metab_in_a, metab_in_b, jiaoji, img_path)
                # hyperlink_list.append()

        if os.path.exists(outdir + "_tmp"):
            shutil.rmtree(outdir + "_tmp")
        level_data.to_csv(level_path, sep="\t")

    def run_img_v1(self, ko, path_C, outdir, path_C_a=None, path_C_b=None):
        # print "pathway: %s" % ko
        if ko.startswith("map"):
            ko = ko.replace("map","ko")
        map = ko.replace("ko", "map")
        png_coll = self.ref_db['kegg_pathway_png_v1']
        result = png_coll.find_one({"pathway_id": ko})
        if result:
            kgml_id = result['pathway_ko_kgml']
            png_id = result['pathway_ko_png']
            # kgml_path = os.path.join(outdir + "_tmp", "pathway.kgml")
            # png_path = os.path.join(outdir + "_tmp", "pathway.png")
            kgml_path = os.path.join(outdir + "_tmp", ko + ".kgml")
            png_path = os.path.join(outdir + "_tmp", ko + ".png")
            # pathway_dir = os.path.join(outdir, "pathway_img")
            pathway_dir = outdir
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            if not os.path.exists(outdir + "_tmp"):
                os.mkdir(outdir + "_tmp")
            fs = gridfs.GridFS(self.ref_db)  # 这样写是否可行？
            with open(kgml_path, "w+") as k, open(png_path, "w+") as p:
                k.write(fs.get(kgml_id).read())
                p.write(fs.get(png_id).read())
            p_kgml = KGML_parser.read(open(kgml_path))
            p_kgml.image = png_path
            for ortholog in p_kgml.orthologs:
                for g in ortholog.graphics:
                     g.bgcolor = colors.Color(alpha=0)
            for each in path_C:
                self.change_color(each, "#00FF00", p_kgml)
            if path_C_a:
                for each in list(path_C_a):
                    self.change_color(each, "#FFA500", p_kgml)
            if path_C_b:
                for each in list(path_C_b):
                    self.change_color(each, "#0000FF", p_kgml)
            for each in path_C:  # 共有的颜色覆盖掉独有分组颜色
                self.change_color(each, "#00FF00", p_kgml)
            canvas = KGMLCanvas(p_kgml, import_imagemap=True, label_compounds=False, label_orthologs=False,
                                label_reaction_entries=False, label_maps=True, show_maps=False,
                                draw_relations=True, show_orthologs=True, show_compounds=True, show_genes=True,
                                show_reaction_entries=False)
            pdf = pathway_dir + '/' + ko + '.pdf'
            png = pathway_dir + '/' + ko + '.png'
            html_path = pathway_dir + '/' + ko + '.html'
            db_file = os.path.join(self.html_path, map + '.html')
            png_base_name = ko + '.png'
            self.get_html(db_file, html_path, png_base_name)
            # if os.path.isfile(db_file):
            #     os.link(db_file, html_path)
            # else:
            #     print "db_file: %s not exists" % db_file
            # map_html = KeggHtml()
            # map_html.run("/mnt/ilustre/users/sanger-dev/app/database/KEGG/kegg_2017-05-01/kegg/pathway/map/" + ko + '.html', html_path, ko + '.png', self.)
            # self.get_pic(r_path, map_path, ko, kos_path, png_path, html_path)
            canvas.draw(pdf)
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + pdf + ' ' + png
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print "圖片格式pdf轉png出錯"
            # try:
            #     canvas.draw(pdf)
            # except Exception,e:
            #     import pickle
            #     with open("pk_debug", "wb") as pk_file:
            #         pickle.dump(target, pk_file)
            #     raise Exception(e)
            # shutil.rmtree(outdir + "_tmp")
        else:
            # raise Exception("{} id 找不到pathway图片".format(ko))
            print "{} id 找不到pathway图片".format(ko)

    def get_html(self, db_file, html_path, pic_name):
        if not os.path.isfile(db_file):
            print "不存在html图片{}".format(db_file)
            return
        html = lxml.html.parse(db_file)
        root = html.getroot()
        body = root.find('body')
        img = body.find('img')
        if img == None:
            try:
                img = body.find_class('image-block')[0].find('img')
            except:
                print "未找到对应的img值"
                return
        try:
            img.set('src', pic_name)
        except:
            print("%s文件没有img" % db_file)
        html.write(html_path)

    def change_color(self, member, color, p_kgml):
        l = []
        for degree in p_kgml.entries.values():
            if re.search(member, degree.name):
                l.append(degree.id)
        for n in l:
            for graphic in p_kgml.entries[n].graphics:
                graphic.bgcolor = color

    def cal_com(self, ref_ko, metab_str, compound_str, rm_repeat=True):
        if len(ref_ko) == 0:
            return [],[]
        metab_list = metab_str.split(";")
        compound_list = compound_str.split(";")
        return_metab = []
        return_compound = []
        conv = {}
        for index,i in enumerate(metab_list):
            if i in ref_ko:
                return_metab.append(i)
                return_compound.append(compound_list[index])
        # return_metab = ';'.join(return_metab)
        # return_compound = ';'.join(return_compound)
        if rm_repeat:
            return list(set(return_metab)), list(set(return_compound))
        else:
            return return_metab, return_compound

    def run(self, set_table, overview_ko, output_dir):
        self.ko_data = pd.read_table(overview_ko, index_col="pathway_id")
        set_table_list = []
        with open(set_table) as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")[1]
                line = line.split(',')
                set_table_list.append(line)
        # set_table = set_table.split(",")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        if len(set_table_list) == 1:
            stat_path = os.path.join(output_dir, "stat.xls")
            self.seta_ko_dic = self.one_set_run(set_table_list[0], stat_path)
        elif len(set_table_list) == 2:
            stat_path1 = os.path.join(output_dir, "stat1.xls")
            stat_path2 = os.path.join(output_dir, "stat2.xls")
            print "record set1: %s" % set_table[0]
            print "record set2: %s" % set_table[1]
            self.seta_ko_dic = self.one_set_run(set_table_list[0], stat_path1)
            self.setb_ko_dic = self.one_set_run(set_table_list[1], stat_path2)
        else:
            raise Exception("不能用超过两个基因集，现在个数%s" % len(set_table_list))
        level_path = os.path.join(output_dir, "level.xls")
        img_path = os.path.join(output_dir, "pathway_img")
        self.set_level_output(level_path, img_path)

    #20190722
    def run_img(self, ko, path_C, outdir, path_C_jiao=None, path_C_b=None):
        html_db = Config().SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007/html/"
        map_pre = ko
        if ko.startswith("map"):
            ko = ko.replace("map","ko")
        else:
            map_pre = ko.replace("ko",'map')

        kgml_ori_path = os.path.join(html_db, map_pre +'.kgml')
        png_ori_path = os.path.join(html_db, map_pre +'.png')
        html_ori_path =  os.path.join(html_db, map_pre +'.html')

        if  os.path.exists(kgml_ori_path) and os.path.exists(png_ori_path) and os.path.exists(html_ori_path) and os.path.getsize(kgml_ori_path) != 0 :
            pathway_dir = outdir
            #pathway_dir = os.path.join(outdir, "pathway_img")
            kgml_path = os.path.join(outdir, "pathway.kgml")
            png_path = os.path.join(pathway_dir, "%s.png" % ko)

            # if not os.path.exists(outdir):
            #     os.mkdir(outdir)
            if not os.path.exists(pathway_dir):
                os.mkdir(pathway_dir)
            if os.path.exists(kgml_path):
                os.remove(kgml_path)
            os.link(kgml_ori_path, kgml_path)
            if os.path.exists(png_path):
                os.remove(png_path)
            os.link(png_ori_path, png_path)

            p_kgml = KGML_parser.read(open(kgml_path))
            p_kgml.image = png_path
            for ortholog in p_kgml.orthologs:
                for g in ortholog.graphics:
                     g.bgcolor = colors.Color(alpha=0)
            compound_color_dic = {}


            if path_C:
                for each in list(path_C):
                    self.change_color(each, "#FF0000", p_kgml)  # 红色
                    compound_color_dic[each] = '#FF0000'

            if path_C_jiao:
                for each in path_C_jiao:
                    self.change_color(each, "#0000FF", p_kgml)  # 蓝色
                    compound_color_dic[each] = '#0000FF'

            if path_C_b:
                for each in list(path_C_b):
                    self.change_color(each, "#00FF00", p_kgml)  # 绿色
                    compound_color_dic[each] = '#00FF00'
            canvas = KGMLCanvas(p_kgml, import_imagemap=True, label_compounds=False, label_orthologs=False,
                                label_reaction_entries=False, label_maps=True, show_maps=False,
                                draw_relations=True, show_orthologs=True, show_compounds=True, show_genes=True,
                                show_reaction_entries=False)
            pdf = pathway_dir + '/' + ko + '.pdf'

            canvas.draw(pdf)
            os.remove(kgml_path)

            """ produce *.mark   """
            html_path = os.path.join(pathway_dir, ko + '.html')
            self.get_html(html_ori_path, html_path, ko)
            all_c = list(path_C)
            if path_C_jiao:
                all_c.extend(list(path_C_jiao))
            if path_C_b:
                all_c.extend(list(path_C_b))
            if os.path.exists(html_path):
                map_html = KeggHtml(html_path, html_path + ".mark", map_pre, all_c)
                map_html.get_html_compound_mark(compound_color_dic=compound_color_dic)

        else:
            # raise Exception("{} id 找不到pathway图片".format(ko))
            print "{} id 找不到pathway图片".format(ko)