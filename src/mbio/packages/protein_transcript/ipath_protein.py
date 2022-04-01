#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/8 17:19
@file    : ipath_protein.py
@author  : yitong.feng 照搬刘彬旭老师的做法
@contact: yitong.feng@majorbio.com
"""


import os
import sys
import xml.etree.ElementTree as ET
import sys

class Ipath(object):
    def __init__(self, kegg_table, protein_list, prefix=''):
        self.db_dir = "/mnt/ilustre/centos7users/yitong.feng/script/ipath/IPATH"
        self.maps = ['Metabolic_pathways', 'Regulatory_pathways', 'Biosynthesis_of_secondary_metabolities']
        self.K = dict()
        self.legend = []
        self.kegg_table = kegg_table
        self.protein_list = protein_list.split(',')
        if len(self.protein_list) > 2:
            print('传入的list文件数量大于2，不好意思做不了')
        self.prefix = prefix
        self.p2K = dict()
        self.out = ''

    def set_legend(self):
        '''
        添加图片标签
        '''
        if len(self.protein_list) == 2:
            n = 1
            for pro in self.protein_list:
                pro = os.path.basename(pro)
                if 'up' in pro:
                    self.legend.append(self.prefix + '_up')
                elif 'down' in pro:
                    self.legend.append(self.prefix + '_down')
                else:
                    self.legend.append(self.prefix + '_' + str(n))
                    n += 1
        else:
            self.legend.append(self.prefix)

        if len(self.legend) == 2:
            self.out = self.prefix + '_' + self.legend[0].split('_')[-1] + '_' + self.legend[1].split('_')[
                -1] + '_'
        else:
            self.out = self.prefix + '_'

    def set_p2K_dict(self):
        with open(self.kegg_table, 'r') as kr:
            _ = kr.readline()
            for line in kr:
                if line.strip():
                    line = line.strip().split('\t')
                    K2p = line[6]
                    if u'(' in K2p:
                        for i in K2p.split(');'):
                            K, ps = i.split('(')
                            if K not in self.p2K:
                                self.p2K[K] = list()
                            self.p2K[K] += ps.split(';')

    def get_input_files(self):
        proteinset2 = list()
        if len(self.protein_list) == 2:
            with open(self.protein_list[0]) as p1, open(self.protein_list[1]) as p2:
                proteinset1 = p1.read().split('\n')
                proteinset2 = p2.read().split('\n')
        else:
            with open(self.protein_list[0]) as p1:
                proteinset1 = p1.read().split('\n')

        out1 = self.out + 'ipath_input.xls'
        out2 = self.out + 'protein_ipath_input.xls'
        with open(out1, 'w') as w1, open(out2, 'w') as w2:
            for K, ps in self.p2K.items():
                p1 = list(set(ps).intersection(set(proteinset1)))
                p2 = list(set(ps).intersection(set(proteinset2)))
                if p1 and p2:
                    w1.write("{}\t#0000ff\tW15\n".format(K))
                    w2.write("{}\t{}\t#0000ff\tW15\n".format(";".join(set(p1 + p2)), K))
                elif p1:
                    w1.write("{}\t#ff0000\tW15\n".format(K))
                    w2.write("{}\t{}\t#ff0000\tW15\n".format(";".join(set(p1)), K))
                elif p2:
                    w1.write("{}\t#00ff00\tW15\n".format(K))
                    w2.write("{}\t{}\t#00ff00\tW15\n".format(";".join(set(p2)), K))
                else:
                    pass
        return out1

    def get_K_color_width(self, color_width_file):
        '''
        获得KO和线条长度颜色信息
        '''
        with open(color_width_file, 'r') as f:
            for line in f.readlines():
                k, color, width = line.strip().split('\t')
                width=width.replace('W','')
                if self.K.has_key(k):
                    pass
                else:
                    self.K.update({k:{"color":color, "width":width}})

    def map_ko_to_svg(self, K, id_file):
        '''
        将KO的颜色长度信息映射至图形svg
        '''
        svg_label = dict()
        with open(id_file, 'r') as f:
            for line in f.readlines():
                label, first_K, all_K = line.strip().split('\t')
                all_ele = [x.strip().split(" ")[0] for x in all_K.split(",")]
                for ele in all_ele:
                    if K.has_key(ele):
                        if svg_label.has_key(label):
                            if K[ele]["color"] == "#0000ff":
                                svg_label.update({label:{"color":K[ele]["color"],"width":K[ele]["width"]}})
                            # pass
                        else:
                            svg_label.update({label:{"color":K[ele]["color"],"width":K[ele]["width"]}})
                    else:
                        pass
        return svg_label

    def set_svg_legend(self, element):
        '''
        修改图例
        '''
        leg_col = element.find('{http://www.w3.org/2000/svg}path')
        text=element.find('{http://www.w3.org/2000/svg}text')
        text.text = self.legend[0]
        text_x = int(text.attrib['x'])
        text_font_size = int(text.attrib['style'].split(";")[0].strip()[-4:-2])
        if len(self.legend) == 2:
            keg_col2= leg_col.copy()
            text2 = text.copy()
            keg_col3= leg_col.copy()
            text3 = text.copy()

            text2_x = text_x  + int(text_font_size * 0.6 * len(self.legend[0])) + 40
            keg_col2.attrib['d'] = 'M{},30L{},30'.format(text2_x , text2_x+100)
            text2.attrib['x'] = str(text2_x + 150)
            text2.attrib['fill'] = "#00ff00"
            keg_col2.attrib['stroke'] = "#00ff00"

            text3_x = text2_x + 150 + int(text_font_size * 0.6 * len(self.legend[1])) + 40
            keg_col3.attrib['d'] = 'M{},30L{},30'.format(text3_x , text3_x+100)
            text3.attrib['x'] = str(text3_x + 150)

            keg_col3.attrib['stroke'] = "#0000ff"
            text3.attrib['fill'] = "#0000ff"

            text2.text = self.legend[1]
            # print "legend 2 is {}".format(self.legend[1])
            text3.text = 'both'

            element.insert(2, keg_col2)
            element.insert(3, text2)
            element.insert(4, keg_col3)
            element.insert(5, text3)

    def map_file(self):
        '''
        修改svg图形改变其中的线条颜色和宽度
        '''
        for path_map in self.maps:
            map_file = os.path.join(self.db_dir, path_map + ".raw.svg")
            id_file = os.path.join(self.db_dir, path_map + ".map")
            svg_lable = self.map_ko_to_svg(self.K, id_file)
            out_file = self.out + path_map + ".svg"
            xml = ET.parse(map_file)
            root = xml.getroot()
            layer_a = root.findall("{http://www.w3.org/2000/svg}g")
            for one_of_a in layer_a:
                if one_of_a.attrib['id'] == 'title':
                    continue
                if one_of_a.attrib['id'] == 'legend':
                    self.set_svg_legend(one_of_a)
                    continue
                layer_b = one_of_a.findall("{http://www.w3.org/2000/svg}g")
                for one_of_b in layer_b:
                    # 修改标签

                    # print self.legend
                    if one_of_b.attrib['id'] == 'title':
                        continue
                    if one_of_b.attrib['id'] == 'legend':
                        self.set_svg_legend(one_of_b)
                        continue

                    layer_c = one_of_b.findall("{http://www.w3.org/2000/svg}g")
                    for one_of_c in layer_c:
                        m_id = "no_id"
                        if one_of_c.attrib.has_key('id'):
                            m_id = one_of_c.attrib['id']
                        else:
                            pass
                        layer_d = one_of_c.findall("{http://www.w3.org/2000/svg}path")
                        layer_d2 = one_of_c.findall("{http://www.w3.org/2000/svg}ellipse")
                        for one_of_d in layer_d + layer_d2:
                            # print m_id,
                            attr_list = ['style', 'stroke-width' , 'stroke','fill']
                            if svg_lable.has_key(m_id) and one_of_d.attrib.has_key("stroke-width"):
                                one_of_d.attrib["stroke-width"] = svg_lable[m_id]["width"]
                                one_of_d.attrib["style"] =  ' opacity: 0.8; '
                            else:
                                pass

                            if svg_lable.has_key(m_id) and one_of_d.attrib.has_key("stroke"):
                                one_of_d.attrib["stroke"] = svg_lable[m_id]["color"]
                                if svg_lable[m_id]["color"] == "#0000ff":
                                    one_of_d.attrib["style"] =  ' opacity: 0.8; '
                                elif len(self.legend) == 2:
                                    one_of_d.attrib["style"] =  ' opacity: 0.4; '
                            else:
                                pass

            xml.write(out_file)

    def run(self):
        self.set_legend()
        self.set_p2K_dict()
        file = self.get_input_files()
        self.get_K_color_width(file)
        self.map_file()
        self.convert()

    def convert(self):
        for file in os.listdir('.'):
            if file.endswith('.svg') and file.startswith(self.out):
                pdf = file.split('.svg')[0] + '.pdf'
                cmd = 'convert %s %s'%(file, pdf)
                os.system(cmd)
                png = file.split('.svg')[0] + '.png'
                cmd = 'convert -flatten -quality 100 -density 130 -background white %s %s' % (pdf, png)
                os.system(cmd)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        exit('USAGE: %s protein_list pathway_table prefix'%sys.argv[0])
    protein_list = sys.argv[1]
    pathway_table = sys.argv[2]
    prefix = sys.argv[3]
    Ipath1 = Ipath(pathway_table, protein_list, prefix)
    Ipath1.run()