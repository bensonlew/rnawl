#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2020/5/6 12:33
@file    : get_protein_information_from_picture.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

from biocluster.config import Config
from PIL import Image
import pandas as pd
import numpy as np
import math
import glob
import os
# import pytesseract
import re


crop_black = True
train_pictures_dir = Config().PACKAGE_DIR + '/tool_lab/train_pictures'
# test_picture = 'Y:\\蛋白代谢事业部\\蛋白与代谢实验部\\客户项目信息\\2019年项目\\2019-10\\MJ20191031086-常二梅-侧柏-愈伤组织-iTRAQ-10个\\搜库结果\\新数据库搜库结果\\Ck25CEM_bjhb1_second_Result.PNG'
# test_picture = 'Y:\\蛋白代谢事业部\\蛋白与代谢实验部\\客户项目信息\\2019年项目\\2019-10\\MJ20191031086-常二梅-侧柏-愈伤组织-iTRAQ-10个\\搜库结果\\老师提供旧转录组的数据库\\Ck25CEM_bjhb_Result.JPG'
# test_picture = 'Y:\\蛋白代谢事业部\\蛋白与代谢实验部\\客户项目信息\\2019年项目\\2019-12\\MJ20191231020-冯胜勇-小鼠肝脏-TMT-6个\\搜库结果\\Db10FSY_bjhb1_Result.PNG'
# test_picture = 'Y:\\蛋白代谢事业部\\蛋白与代谢实验部\\客户项目信息\\2019年项目\\2019-12\\MJ20191206067-陈润-耻垢分枝杆菌-TMT-9个\\搜库结果\\Cl09CR_bjhb1_Result.JPG'
# test_picture = ''
def trim(im, cut=2):
    x1 = 0
    y1 = 0
    x2, y2 = im.size
    x_t, y_t = im.size
    b = False
    for x in range(0, x_t):
        if b:
            x1 = x
            break
        for y in range(0, y_t):
            pix = im.getpixel((x, y))
            if pix != 255:
                b = True
                break
    b = False
    for y in range(0, y_t):
        if b:
            y1 = y
            break
        for x in range(0, x_t):
            pix = im.getpixel((x, y))
            if pix != 255:
                b = True
                break
    b = False
    for x in range(x_t - 1, x1, -1):
        if b:
            x2 = x
            break
        for y in range(y_t - 1, y1, -1):
            pix = im.getpixel((x, y))
            if pix != 255:
                b = True
                break
    b = False
    for y in range(y_t - 1, y1, -1):
        if b:
            y2 = y
            break
        for x in range(x_t - 1, x1, -1):
            pix = im.getpixel((x, y))
            if pix != 255:
                b = True
                break
    return im.crop((x1-cut, y1-cut, x2+cut, y2+cut))



def convert_to_blackwhite(im):
    im.convert("P")
    x, y = im.size
    if crop_black:
        im = im.crop((0+2, 0+2, x-2, y-2))
    im2 = Image.new("P", im.size, 255)
    for x in range(im.size[1]):
        for y in range(im.size[0]):
            pix = im.getpixel((y, x))
            if pix[0] < 60:
                im2.putpixel((y, x), 0)
    return im2

def get_draw_pos(im):
    for x in range(im.size[1]):
        for y in range(im.size[0]):
            pix = im.getpixel((y, x))
            if pix != 255:
                return (y, x)


def split_image_to_single(im):
    s_ims = list()
    inletter = False
    foundletter = False
    start = 0
    end = 0
    letters = []

    for y in range(im.size[0]):
        for x in range(im.size[1]):
            pix = im.getpixel((y, x))
            if pix != 255:
                inletter = True
        if foundletter == False and inletter == True:
            foundletter = True
            start = y

        if foundletter == True and inletter == False:
            foundletter = False
            end = y
            letters.append((start, end))

        inletter = False

    for letter in letters:
        im_ = im.crop((letter[0], 0, letter[1], im.size[1]))
        s_ims.append(im_)

    return s_ims


def get_all_possibility():
    all_poss_d = list()
    for png in glob.glob(os.path.join(train_pictures_dir, '*.png')):
        possiable = os.path.basename(png).split('.png')[0].strip().split('+')[0]
        image = Image.open(png)
        image.convert('P')
        image = trim(image, 2)
        all_poss_d.append({
            'possiable': possiable,
            'image': image
        })
    return pd.DataFrame(all_poss_d)


class VectorCompare:
    def magnitude(self, concordance):
        total = 0
        for word, count in concordance.items():
            total = total + count ** 2
        return math.sqrt(total)

    def relation(self, concordance1, concordance2):
        topvalue = 0
        for word, count in concordance1.items():
            if word in concordance2:
                topvalue = topvalue + count * concordance2[word]
        return topvalue / ((self.magnitude(concordance1) * self.magnitude(concordance2)) + 0.00001)


#返回的是一個字典 key是count，（count是自加的，代辦像素的個數，value 是 像素數字 圖片的8bit一個像素，正好是0-255之間的元素）
def buildvector(im):
    d1 = dict()
    count = 0
    for i in im.getdata():
        d1[count] = i
        count += 1
    return d1

def image_similarity_vectors_via_numpy(image1, image2):
    image1 = image1.convert('L')
    image2 = image2.convert('L')
    images = [image1, image2]
    vectors = []
    norms = []
    for image in images:
        vector = []
        for pixel_tuple in image.getdata():
            vector.append(np.average(pixel_tuple))
        vectors.append(vector)
        # linalg=linear（线性）+algebra（代数），norm则表示范数
        # 求图片的范数？？
        norms.append(np.linalg.norm(vector, 2))
    a, b = vectors
    a_norm, b_norm = norms
    # dot返回的是点积，对二维数组（矩阵）进行计算
    res = np.dot(a / a_norm, b / b_norm)
    return res


def get_most_simility(df, im):
    def cal_rela_(row):
        t_im = im.resize(row['image'].size, Image.ANTIALIAS)
        im_d = buildvector(t_im)
        v = VectorCompare()
        return v.relation(buildvector(row['image']), im_d)
    def cal_rela(row):
        t_im = im.resize(row['image'].size, Image.ANTIALIAS)
        return image_similarity_vectors_via_numpy(t_im, row['image'])
    df['rela'] = df.apply(cal_rela, axis=1)
    df = df.sort_values(by='rela', ascending=False)
    df = df.reset_index(drop=True)
    return df.iloc[0, :]['possiable']


def get_chars_from_picture(pic_files):
    char = ''
    im = Image.open(pic_files)
    # print(pytesseract.image_to_string(im))
    im = convert_to_blackwhite(im)
    # im.show()
    # print(pytesseract.image_to_string(im))
    try:
        im = trim(im, cut=1)
    except:
        im = trim(im, cut=1)
    # print(pytesseract.image_to_string(im))
    # im.show()
    split_ims = split_image_to_single(im)
    # pos = get_draw_pos(split_ims[0])
    all_p_df = get_all_possibility()
    n = 30
    for s_m in split_ims:
        s_m = trim(s_m, 2)
        # s_m.save(os.path.join('train_pictures', str(n) + '.png'))
        n += 1
        c_ = get_most_simility(all_p_df.copy(), s_m)
        char += c_
    return char


def export_p_info(chars, dir_):
    chars = chars.lower()
#     '4530___5884Proteins----4530ProteinGroups----45733PePtideGroups----109192PSMs----302401MS___MSSPectrumInfo'
    res = re.search('(?P<protein>\d+_*\d*).*?proteins----(?P<proteingroup>\d+).*?proteingrou\ws----(?P<pepgroup>\d+).*?pe\wtidegrou\ws----(?P<psm>\d+).*?psms----(?P<ms>\d+).*?ms', chars)
    res_d = res.groupdict()
    try:
        tmp = res_d['protein'].split('_')
        res_d['protein'] = tmp[-1]
        res_d['proteingroup'] = tmp[0]
    except Exception:
        pass
    p_f = os.path.join(dir_, 'Protein_information.xls')
    with open(p_f, 'w') as pw:
        pw.write('Total Spectrum\tIdentified Spectrum\tPeptide number\tProtein number\tProtein group number\n')
        # pw.write(res_d['ms'] + '\t' + res_d['psm'] + '\t' + res_d['pepgroup'] + '\t' + res_d['protein'] + '\t' + res_d['proteingroup'] + '\n')
        pw.write('{ms}\t{psm}\t{pepgroup}\t{protein}\t{proteingroup}\n'.format(**res_d))
    with open(p_f) as pr:
        print(pd.read_csv(pr, sep='\t'))

def run():
    chars = get_chars_from_picture(picture)
    print(chars)
    export_p_info(chars, output)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Extract the data from the protein picture')
    parser.add_argument('-picture_path', type=str,
                        help='picture_path')
    parser.add_argument('-if_black', type=str,
                        help='if have black border')
    parser.add_argument('-output', type=str, help='output path')
    args = parser.parse_args()
    global picture, crop_black, output
    picture = args.picture_path
    output = args.output
    if args.if_black == 'no':
        crop_black = False
    try:
        run()
    except Exception:
        crop_black = not crop_black
        run()
