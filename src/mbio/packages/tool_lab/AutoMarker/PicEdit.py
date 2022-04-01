#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 2021/2/26 13:51
# @Author  : U make me wanna surrender my soul
import faulthandler
faulthandler.enable()
import sys
import cv2
import numpy as np
import tensorflow as tf
import xlrd
from PIL import ImageFont, ImageDraw, Image
from collections import defaultdict
import argparse
import os


# --------------Model preparation----------------
# Path to frozen detection graph. This is the actual model that is used for the object detection.
# PATH_TO_CKPT = r'inference_graph/frozen_inference_graph.pb'

def main(PATH_TO_CKPT, input_file, input_table, out_file, font_path):
    np.set_printoptions(threshold=np.inf)
    # Load a (frozen) Tensorflow model into memory
    detection_graph = tf.Graph()
    with detection_graph.as_default():
        od_graph_def = tf.compat.v1.GraphDef()
        with tf.compat.v1.gfile.GFile(PATH_TO_CKPT, 'rb') as fid:
            serialized_graph = fid.read()
            od_graph_def.ParseFromString(serialized_graph)
            tf.import_graph_def(od_graph_def, name='')

    image_tensor = detection_graph.get_tensor_by_name('image_tensor:0')
    # Each box represents a part of the image where a particular
    # object was detected.
    gboxes = detection_graph.get_tensor_by_name('detection_boxes:0')
    # Each score represent how level of confidence for each of the objects.
    # Score is shown on the result image, together with the class label.
    gscores = detection_graph.get_tensor_by_name('detection_scores:0')
    gclasses = detection_graph.get_tensor_by_name('detection_classes:0')
    gnum_detections = detection_graph.get_tensor_by_name('num_detections:0')

    if not os.path.exists(input_file):
        print('请先创建该文件夹，并将未编辑的胶图移动至该文件夹')
        sys.exit()
    if not os.path.exists(input_table):
        print('请检查当前目录下是否存在该文件：{input_table}'.format(input_table=args.t))
        sys.exit()
    if not os.path.exists(out_file):
        os.mkdir(out_file)

    with detection_graph.as_default():
        with tf.compat.v1.Session(graph=detection_graph) as sess:
            font = ImageFont.truetype(font_path, 13)
            Dict_marker = defaultdict(dict)
            data = xlrd.open_workbook(input_table)  # 在此处修改输入表格
            table = data.sheets()[0]
            nrows = table.nrows
            if (nrows + 1) % 6 != 0:
                print("您的表格记录行数有误，请核对后再运行")
                sys.exit(0)
            ncols = table.ncols
            if ncols % 4 != 0:
                print("您的表格记录行数有误，请核对后再运行")
                sys.exit(0)
            for n in range(int((nrows + 2) / 6)):
                imgFile = str(table.row_values(6 * n)[0])
                for i in range(1, 5):
                    list1 = table.row_values(i + 6 * n)
                    list1 = [k for k in list1 if k != '']
                    if len(list1) % 4 != 0:
                        print("您的表格记录有误，请核对后再运行--" + "__".join([str(a) for a in list1]))
                        sys.exit(0)
                    Dict_marker[imgFile][i] = list1  # 循环添加元素到创建的空字典
            for key in Dict_marker:
                tif = str(input_file + '/' + key + '.tif')
                jpg = tif + ".jpg"
                print(tif)
                if not os.path.exists(tif):
                    print("该文件夹没有此图片！！ " + tif)
                    continue
                img = cv2.imread(tif)  # 读入胶图文件
                cv2.imwrite(out_file + '/' + jpg, img)  # 创建胶图副本
                # img_marker = cv2.imread("marker.png", 0)  # 读入marker图
                line_num = int((img.shape[0]) / 4)  # 计算胶图四等分 cv2.imread如果是三通道读取图片，则图片的高度在前，长度在后
                # if k == 1:
                last_split = img[0:(line_num + 50), 0:img.shape[1]]  # 取胶图图片的第一张
                # print(last_split)
                # last_split = img[((k - 1) * line_num):(k * line_num), 0:img.shape[1]]
                last_img, last_coor = detect_image_objects(last_split, sess, gboxes, gscores, gclasses, gnum_detections, image_tensor)
                try:
                    mx = int(np.sort(last_coor, axis=0)[-1][1])  # 取每一行的最后一个marker的横坐标
                    # print(mx)
                    my = int(np.sort(last_coor, axis=0)[-1][0])  # 同上纵坐标
                    if my < 20:
                        my = 20
                    # print(my)
                except IndexError:
                    print('未找到{pic}的marker，请检查图片'.format(pic=input_file + '/' + key + '.tif'))

                axisX = 0
                space_bak = 0
                # 构造列表，存放每个客户的点样数量信息，对应marker表内每行的数据
                for k in range(1, 5):
                    # print(k)
                    strq = Dict_marker[key][k]
                    marker_line_list = []
                    for num in range(int(len(strq) / 4)):
                        numMin = strq[1 + num * 4]
                        numMax = strq[2 + num * 4]
                        if isinstance(numMin, str) or isinstance(numMax, str):
                            print("您的数据记录此行有误，请核对" + "_".join([str(a) for a in strq]))
                            break
                        else:
                            numMin = int(numMin)
                            numMax = int(numMax)
                        for j in range(numMin, numMax + 1):
                            marker_line_list.append(str(j))
                        if strq[num * 4 + 3] == "CK" or strq[num * 4 + 3] == "ck":
                            marker_line_list.append("CK")
                        marker_line_list.append("M")

                    # 对原始文件进行分块操作
                    if k < 4:
                        img_split = img[((k - 1) * line_num):(k * line_num + 50), 0:img.shape[1]]
                    else:
                        img_split = img[((k - 1) * line_num):(k * line_num), 0:img.shape[1]]
                    # no_less_img 对应的图像对象
                    no_less_img, Match_coor = detect_image_objects(img_split, sess, gboxes, gscores, gclasses, gnum_detections, image_tensor)  # 返回匹配目标的像素矩阵
                    # print(Match_coor)
                    if Match_coor.size == 0:
                        if axisX == 0:
                            print("Sorry,你的图片第" + str(k) + "行匹配有误！！请重新截图marker！！ " + key + ".tif")
                            continue
                        else:
                            y = axisY + line_num
                            x = axisX
                    if Match_coor.size > 0:
                        lst = (Match_coor[:, 1]).tolist()  # 取x位置，转为list
                        lst = sorted(set(lst))  # 对list排序
                        y = int(np.mean(Match_coor, axis=0)[0] + (k - 1) * line_num)
                        x_marker = judge(lst)
                        if len(x_marker) < 2:
                            # print(axisX)
                            if axisX == 0:
                                print("Sorry,你的图片此行匹配有误！！请重新截图marker！！ " + key + ".tif")
                                continue
                            else:
                                x = axisX
                                space = space_bak
                        else:
                            x = int(x_marker[0])  # 删除 +15
                            len1 = x_marker[-1] - x_marker[0]
                        if len1 > 200:
                            try:
                                space = len1 / int(len(marker_line_list))
                            except ZeroDivisionError:
                                print('您输入的数据表记录有误！请检查该文件！')
                            # print(space)
                        else:
                            if space_bak != 0:
                                space = space_bak

                        if abs(x - axisX) > 50 and axisX != 0:
                            x = axisX
                            space = space_bak
                    axisY = y
                    axisX = x
                    space_bak = space

                    img_pi = Image.fromarray(img)
                    draw = ImageDraw.Draw(img_pi)
                    draw.text((x + 20, y - 20), "M", font=font, fill=(0, 255, 255))  # -15 -10 -20  修改第一个Marker的位置
                    img = np.array(img_pi)
                    m = 0
                    ck_marker_count = [x]
                    # print('ck:', ck_marker_count)
                    for j in marker_line_list:
                        m = m + 1
                        # 调整数字位置
                        if str.isdigit(j) and int(j) < 10:  # 个位数
                            pos = m * space + x + 15  # 10
                            # print('个位数:', pos)
                        elif str.isdigit(j) and 10 <= int(j) < 100:  # 两位数
                            pos = m * space + x + 10  # 2
                            # print('两位数:', pos)
                        elif str.isdigit(j) and int(j) >= 100:
                            pos = m * space + x + 10
                        else:  # 三位数
                            pos = m * space + x  # no
                            # print('三位数:', pos)
                        if j == 'CK' or j == 'M' or m == len(marker_line_list):
                            ck_marker_count.append(pos)
                            # print('DLorCK:', ck_marker_count)
                        if j == 'M':
                            img = drawText(img, pos + 10, y - 20, str(j), font_path)  # y - 15 修改尾数Marker的位置
                        else:
                            img = drawText(img, pos, y - 20, str(j), font_path)  # 修改数字位置
                        cv2.imwrite(out_file + '/' + key + '.jpg', img)
                cv2.imwrite(out_file + '/' + jpg, img)
                print("即将进行下一张胶图编辑...")

    if os.path.exists('add_text.jpg'):
        os.remove('add_text.jpg')
    print('程序运行完毕')

def detect_image_objects(image, sess, gboxes, gscores, gclasses, gnum_detections, image_tensor):
    # Expand dimensions since the model expects images to have
    # shape: [1, None, None, 3]
    image_np_expanded = np.expand_dims(image, axis=0)

    # Actual detection.
    (boxes, scores, classes, num_detections) = sess.run(
        [gboxes, gscores, gclasses, gnum_detections],
        feed_dict={image_tensor: image_np_expanded})

    # Visualization of the results of a detection.
    boxes = np.squeeze(boxes)
    scores = np.squeeze(scores)
    height, width = image.shape[:2]
    match = np.empty([0, 2])
    for u in range(boxes.shape[0]):
        if scores is None or scores[u] > 0.95:  # 设置显示阈值
            x_list = []
            y_list = []
            ymin, xmin, ymax, xmax = boxes[u]
            ymin = int(ymin * height)
            xmin = int(xmin * width)
            text_x = np.max((0, xmin))
            text_y = np.max((0, ymin))
            x_list.append(text_x)
            y_list.append(text_y)
            loc = zip(y_list, x_list)
            for k in loc:
                match = np.append(match, [k], axis=0)

    return image, match  # 返回点的坐标


def drawText(img_1, x_1, y_1, text, font_path):
    img_pil = Image.fromarray(img_1)
    draw_ = ImageDraw.Draw(img_pil)
    font_text = None
    if isinstance(text, str):
        font_text = ImageFont.truetype(font_path, 13)
    else:
        if int(text) < 100:
            font_text = ImageFont.truetype(font_path, 13)
        elif int(text) >= 100:
            font_text = ImageFont.truetype(font_path, 13)
    draw_.text((x_1, y_1), text, font=font_text, fill=(0, 255, 255))
    bk_img = np.array(img_pil)
    return bk_img


# 删除多余的marker
def judge(inlist):
    b = []
    for po in range(len(inlist)):
        if po == len(inlist) - 1:
            pass
        else:
            if inlist[po + 1] - inlist[po] >= 55:
                b.append(inlist[po])
                b.append(inlist[po + 1])
    if b[-1] - b[-2] <= 55:
        del b[-1]
    new_list = []
    for pb in b:
        if pb not in new_list:
            new_list.append(pb)
    return new_list


# 在图片的最后一个marker的后面添加ck
def add_CK_and_M(img_ck, mx, my, font_path):
    font_ck = ImageFont.truetype(font_path, 13)
    ck_pp = Image.fromarray(img_ck)
    draw_ck = ImageDraw.Draw(ck_pp)
    ckx = mx + 20
    cky = my + 2
    draw_ck.text((mx + 70, my + 2), "CK", font=font_ck, fill=(0, 255, 255))
    draw_ck.text((ckx, cky), "M", font=font_ck, fill=(0, 255, 255))
    return ck_pp


def add_only_CK(img_ck, Mx, My, font_path):
    font_ck = ImageFont.truetype(font_path, 13)
    ck_pp = Image.fromarray(img_ck)
    draw_ck = ImageDraw.Draw(ck_pp)
    ckx = Mx + 25
    cky = My + 50
    draw_ck.text((ckx, cky), "CK", font=font_ck, fill=(0, 255, 255))
    return ck_pp

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="输入图片存放文件夹", type=str, required=True)
    parser.add_argument("-t", help="输入对应编号表", type=str, required=True)
    parser.add_argument("-p", help="训练模型的路径", type=str, required=True)
    parser.add_argument("-f", help="字体存放路径", type=str, required=True)
    parser.add_argument("-o", help="输入结果存放文件夹", type=str, required=True)
    args = parser.parse_args()
    PATH_TO_CKPT = args.p
    input_file = args.i
    input_table = args.t
    out_file = args.o
    font_path = args.f
    main(PATH_TO_CKPT, input_file, input_table, out_file, font_path)