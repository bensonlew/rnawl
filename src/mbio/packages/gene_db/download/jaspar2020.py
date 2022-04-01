# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET
import lxml.html
import chardet
import json

class Jaspar(object):
    def __init__(self):
        pass

    def get_html(self, ma_id):
        '''
        设置原始html路径
        '''
        page = 'http://jaspar.genereg.net/matrix/{}/'.format(ma_id)
        return page
        # htmlCode = page.read()
        # return htmlCode


    def get_matrix(self, ma_id):
        '''
        重置html文件中的标签属性
        '''
        html_url = self.get_html(ma_id)
        print html_url
        html = lxml.html.parse(html_url)
        root = html.getroot()
        body = root.find('body')

        detail_dict = dict()

        matrix_detail = body.get_element_by_id("matrix-detail")
        for tr in matrix_detail.getchildren():
            title = tr.getchildren()[0].getchildren()[0].text
            value = tr.getchildren()[1].text
            if title in ["Species:", "Validation:",
                         "Uniprot ID:",	"Pazar TF:",
                         "TFBSshape ID:", "TFencyclopedia IDs:",	"Source:", "Comment:"]:
                eles = list()
                for sub_ele in tr.getchildren()[1].getchildren():
                    eles.append(str(sub_ele.text).strip())
                detail_dict[title] = ";".join(eles)

            else:
                detail_dict[title] = str(value).strip()
        return detail_dict

    def para(self, profile_list, out):
        eles = ["Name:", "Matrix ID: ", "Class:","Family:",
                "Collection:", "Taxon:", "Species:", "Data Type:",
                "Validation:",	"Uniprot ID:",	"Pazar TF:",
                "TFBSshape ID:", "TFencyclopedia IDs:",	"Source:", "Comment:"]
        with open(profile_list, 'r') as f:
            profile_json = json.load(f)
        with open(out, 'w') as out_fo:
            out_fo.write("\t".join(eles) + "\n")
            for mid in profile_json:
                ma_id = mid
                if not ma_id.startswith("MA"):
                    continue
                ma_dict = self.get_matrix(ma_id)
                print ma_dict
                value_list = [ma_dict.get(x, "-") for x in eles]
                out_fo.write("\t".join(value_list) + "\n")

if __name__ == "__main__":
    Jaspar = Jaspar()
    profile_list = sys.argv[1]
    out_file = sys.argv[2]
    Jaspar.para(profile_list, out_file)
