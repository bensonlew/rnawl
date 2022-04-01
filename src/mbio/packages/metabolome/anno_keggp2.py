# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base
import pandas as pd
import os
import subprocess
from biocluster.config import Config
import lxml.html
import re
from reportlab.lib import colors
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

class AnnoKeggp(Base):
    def __init__(self):
        super(AnnoKeggp, self).__init__()
        self._project_type = 'metabolome'
        self.compound_set = set()  # compound_list 改为compound_set by ghd @ 20191014
        self.c2m = {}
        self.stat_dic = {}

    def run(self, table, organisms, level_path, stat_path, color='green',html_db=None):
        """
        table: 所用的数据详情表
        organisms:查询的相关物种XXX;xxx或XXX
        level_path:输出的层级结果
        stat_path:输出的统计结果
        color:绘图用的颜色
        """
        data = pd.read_table(table)
        data.apply(self.process, axis=1, args=())
        table_ko_list = self.find_compound_ko()
        if organisms not in  ["False", "All"]:
            ko_list = self.get_ko_list(organisms)
            ko_set = set(ko_list) & set(table_ko_list)
            print "check ko_list: %s" % set(ko_list)
            print "check table_ko_list: %s" % set(table_ko_list)
        else:
            ko_set = set(table_ko_list)
        ko_result = self.find_kegg(list(ko_set))
        # level_data = pd.DataFrame(columns=["pathway_id", "first_category", "second_category", "description", "compound_id", "metab_id", "count", "hyperlink"])
        level_dict = {
            "pathway_id": [],
            "first_category": [],
            "second_category": [],
            "description": [],
            "compound_id": [],
            "metab_id": [],
            "count": [],
            "hyperlink": []
        }
        for one in ko_result:
            compound_list = one['compound_id'].split(';')
            result_compound = set(compound_list) & set(self.c2m.keys())
            compound_list = list(result_compound)
            count = len(result_compound)
            compound_str = ';'.join(compound_list)
            # metab_list = []
            metab_set = set()  # metab_list改为metab_set by ghd @20191015
            for one_list in compound_list:
                # metab_list.append(self.c2m[one_list])
                metab_set.update(self.c2m[one_list])
            metab_list = list(metab_set)
            metab_str = ';'.join(metab_list)
            level_dict['pathway_id'].append(one['pathway_id'])
            level_dict['first_category'].append(one['first_category'])
            level_dict['second_category'].append(one['second_category'])
            level_dict['description'].append(one['discription'])
            level_dict['compound_id'].append(compound_str)
            level_dict['metab_id'].append(metab_str)
            level_dict['count'].append(count)
            level_dict['hyperlink'].append(self.get_link(one['pathway_id'], compound_list, color=color))  # 需确认方法
            self.run_img(one['pathway_id'], compound_list, os.path.dirname(level_path),html_db=html_db)
            if self.stat_dic.has_key(one['first_category']):
                if self.stat_dic[one['first_category']].has_key(one['second_category']):
                    self.stat_dic[one['first_category']][one['second_category']] += metab_list
                else:
                    self.stat_dic[one['first_category']][one['second_category']] = metab_list
            else:
                self.stat_dic[one['first_category']] = {
                    one['second_category']: metab_list
                }
        level_data = pd.DataFrame(data=level_dict,
                                  columns=["pathway_id", "first_category", "second_category", "description",
                                           "compound_id", "metab_id", "count", "hyperlink"])
        level_data.sort_values(by='count', ascending=False, inplace=True)
        level_data.to_csv(level_path, index=False, sep='\t')
        stat_file = open(stat_path, "w")
        stat_file.write("first_category\tsecond_category\tmetab_id\tcount\n")
        for fc in sorted(self.stat_dic.keys()):
            for sc in sorted(self.stat_dic[fc].keys()):
                if sc == "Global and overview maps":
                    continue  # 不对此二级分类作图 modified by GHD @ 20181205
                metab_list = list(set(self.stat_dic[fc][sc]))
                count = len(metab_list)
                metab_id = ';'.join(metab_list)
                stat_file.write("%s\t%s\t%s\t%s\n" % (fc, sc, metab_id, count))
        stat_file.close()

    def process(self, one_result):
        """
        将一行中的数据存储到类变量self.compound_list/self.c2m中
        此方法兼容了metab含多个compound id 情况 by ghd @20191015
        :param one_result:
        :return:
        """
        # self.compound_list.append(one_result['KEGG Compound ID'])
        self.compound_set.update(one_result["KEGG Compound ID"].split(";"))
        for compound_id in one_result["KEGG Compound ID"].split(";"):
            if self.c2m.has_key(compound_id):
                self.c2m[compound_id].append(one_result["metab_id"])
            else:
                self.c2m[compound_id] = [one_result['metab_id']]

    def get_ko_list(self, organisms):
        """
        根据kegg_organisms查询相关物种的通路信息
        :param organisms: 物种名称 XXX;xxx
        :return:
        """
        tmp_list = organisms.split(";")
        ref_org_db = self.ref_db["kegg_organisms"]
        ko_list = []
        if len(tmp_list) == 1:
            ko_set = set()
            results = ref_org_db.find({"first_category": tmp_list[0]})
            for result in results:
                tmp_ko_list = result['map_list'].split(";")
                ko_set = ko_set | set(tmp_ko_list)
            ko_list = list(ko_set)
        elif len(tmp_list) == 2:
            if tmp_list[0] == 'All':
                result = ref_org_db.find_one({"second_category": tmp_list[1]})
            else:
                result = ref_org_db.find_one({"first_category": tmp_list[0], "second_category": tmp_list[1]})
            ko_list = result['map_list'].split(";")
        ko_list = map(lambda x: 'map' + x, ko_list)
        return ko_list

    def find_compound_ko(self):
        """
        查询kegg_compound参考库，获取在self.compound_list self.compound_set中的化合物相关的所有通路id
        :return:
        """
        ref_kegg_compound = self.ref_db["kegg_compound"]
        ko_list = []
        result = ref_kegg_compound.find({"entry": {"$in": list(self.compound_set)}}) # 将compound_list 改为compound_set by ghd @ 20191014
        for one in result:
            if one["pathway"] != "-":
                ko_list += one["pathway"].split(";")
        return ko_list

    def find_kegg(self, ko_list=None):
        """
        根据ko列表，查询kegg_pathway_level数据库
        :param ko_list:
        :return:
        """
        ref_kegg_db = self.ref_db["kegg_pathway_level"]
        result = ref_kegg_db.find({"pathway_id": {"$in": ko_list}}, no_cursor_timeout=True)
        return result

    def get_link(self, ko, compound_list, color="green"):
        """
        获取网页版的通路图片，对关注的化合物进行着色
        :param ko:  kegg pathway id
        :param compound_list:一组kegg compound id
        :param color:标识的颜色
        :return: link = 'http://www.genome.jp/dbget-bin/show_pathway?' + ko + '/' + '/'.join(ko_color)
        """
        link = 'http://www.genome.jp/dbget-bin/show_pathway?' + ko
        for i in compound_list:
            link += '+' + i + '%09' + color
        return link

    def run_img_v1(self, ko, path_C, outdir):
        """
        绘制本地pdf图片
        :param ko: kegg pathway id
        :param path_C: compound in list format
        :param outdir: output directory
        :return:
        """
        import re
        import gridfs
        from reportlab.lib import colors
        from Bio.KEGG.KGML import KGML_parser
        from Bio.Graphics.KGML_vis import KGMLCanvas
        # print "pathway: %s" % ko
        if ko.startswith("map"):
            ko = ko.replace("map","ko")
        png_coll = self.ref_db['kegg_pathway_png_v1']
        result = png_coll.find_one({"pathway_id": ko})
        if result:
            kgml_id = result['pathway_ko_kgml']
            png_id = result['pathway_ko_png']
            kgml_path = os.path.join(outdir, "pathway.kgml")
            png_path = os.path.join(outdir, "pathway.png")
            pathway_dir = os.path.join(outdir, "pathway_img")

            if not os.path.exists(outdir):
                os.mkdir(outdir)
            if not os.path.exists(pathway_dir):
                os.mkdir(pathway_dir)
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
                l = []
                for degree in p_kgml.entries.values():
                    if re.search(each, degree.name):
                        l.append(degree.id)
                for n in l:
                    for graphic in p_kgml.entries[n].graphics:
                        # print "id: %s\tbefore color: %s" % (n,graphic.fgcolor)
                        graphic.bgcolor = "#00FF00"
            canvas = KGMLCanvas(p_kgml, import_imagemap=True, label_compounds=False, label_orthologs=False,
                                label_reaction_entries=False, label_maps=True, show_maps=False,
                                draw_relations=True, show_orthologs=True, show_compounds=True, show_genes=True,
                                show_reaction_entries=False)
            pdf = pathway_dir + '/' + ko + '.pdf'
            canvas.draw(pdf)
            os.remove(kgml_path)
            os.remove(png_path)
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
            img = body.find_class('image')[0].find('img')
        try:
            img.set('src', pic_name)
        except:
            raise Exception("%s文件没有img" % db_file)
        html.write(html_path)


    def run_img(self, ko, path_C, outdir,html_db=None):  #20190612
        """
        绘制本地pdf图片
        :param ko: kegg pathway id
        :param path_C: compound in list format
        :param outdir: output directory
        :return:
        """
        # print "pathway: %s" % ko
        if not html_db:
            raise Exception('html_db params is None')
        map_pre = ko
        if ko.startswith("map"):
            ko = ko.replace("map","ko")
        else:
            map_pre = ko.replace("ko",'map')

        kgml_ori_path = os.path.join(html_db, map_pre +'.kgml')
        png_ori_path = os.path.join(html_db, map_pre +'.png')
        html_ori_path =  os.path.join(html_db, map_pre +'.html')

        if  os.path.exists(kgml_ori_path) and os.path.exists(png_ori_path) and os.path.exists(html_ori_path) and os.path.getsize(kgml_ori_path) !=0:
            pathway_dir = os.path.join(outdir, "pathway_img")
            kgml_path = os.path.join(outdir, "pathway.kgml")
            png_path = os.path.join(pathway_dir, "%s.png" % ko)

            if not os.path.exists(outdir):
                os.mkdir(outdir)
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
            for each in path_C:
                l = []
                for degree in p_kgml.entries.values():
                    if re.search(each, degree.name):
                        l.append(degree.id)
                for n in l:
                    for graphic in p_kgml.entries[n].graphics:
                        print ("id: %s\tbefore color: %s" % (n, graphic.fgcolor))
                        graphic.bgcolor = "#FF0000" #red   #"#00FF00" green
            canvas = KGMLCanvas(p_kgml, import_imagemap=True, label_compounds=False, label_orthologs=False,
                                label_reaction_entries=False, label_maps=True, show_maps=False,
                                draw_relations=True, show_orthologs=True, show_compounds=True, show_genes=True,
                                show_reaction_entries=False)
            pdf = pathway_dir + '/' + ko + '.pdf'
            canvas.draw(pdf)  #20200420
            os.remove(kgml_path)

            """ produce *.mark   """
            html_path = os.path.join(pathway_dir, ko + '.html')
            self.get_html(html_ori_path, html_path, ko)
            map_html = KeggHtml(html_path, html_path + ".mark", map_pre, path_C)
            map_html.get_html_compound_mark()

        else:
            # raise Exception("{} id 找不到pathway图片".format(ko))
            print ("{} id 找不到pathway参考数据".format(ko))

#zouguanqing
class KeggHtml:
    def __init__(self, map_file, out_file,id, hit_compound_list):
        self.map_file = map_file
        self.map_id = id
        self.out_file = out_file
        self.compounds = hit_compound_list

    def read_html(self):
        """
        读取html文件
        :return: areas
        """
        html = lxml.html.parse(self.map_file)
        root = html.getroot()
        body = root.find('body')

        # 修改图片链接
        img = body.find('img')
        try:
            img.set('src', self.map_id.replace('map','ko'))
        except:
            print "{} has no image src".format(self.map_id)

        # 修改鼠标点击矩形区域链接
        maps = body.find('map')
        areas = maps.findall('area')
        return html,areas


    def get_html_compound_mark(self,compound_color_dic=None):
        """
        获取标签信息，导入用于前端作图
        """
        html,areas = self.read_html()
        with open(self.out_file, 'w') as mark:
            mark.write('ko\tquery\tshape\tcoords\ttitle\thref\tcolor\n')
            for area in areas:
                area_dict=dict(area.items())
                if area_dict['shape'] in ['circle']:
                    kos = area_dict['href'].split('?')[-1].split('+')
                    hit_cs = list(set(kos) & set(self.compounds))
                    color =''
                    if compound_color_dic:
                        for c in hit_cs:
                            if c in compound_color_dic.keys():
                                color = compound_color_dic[c]
                                break
                    # print "change***".join(change_kos)
                    if len(hit_cs) > 0:
                        mark.write("\t".join([self.map_id,','.join(hit_cs), area_dict['shape'], area_dict['coords'], area_dict['title'], area_dict['href'],color]) + "\n")
                    else:
                        mark.write("\t".join([self.map_id, '', area_dict['shape'], area_dict['coords'],  area_dict['title'], area_dict['href'],color]) + "\n")
                elif area_dict['shape'] in ['rect']:
                    mark.write("\t".join([self.map_id, "",area_dict['shape'], area_dict['coords'],  area_dict['title'] , area_dict['href'],'']) + "\n")

