#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2021/4/30 11:24
@file    : find_sample_rawdatas.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""



import os
import pandas as pd

# 本脚本采用强行匹配的方法，规则加的越多越容易报错，运行时间也不一定就能省了
# 有三个原始数据存放的路径，会先通过合同号确定是什么项目，存在那个位置
# 找到在什么位置之后便通过传来的样本列表，找到每个样本对应的原始数据的绝对路径

class not_find_MJ(Exception): pass  #没有发现项目会报这个错误


def get_project_dir(MJ, PM=''):
    """
    :param MJ: 美吉的合同编号
    :param PM: 任务号，如果没有便为空
    :return: 项目所在的绝对路径
    """
    raw_dirs = [
        '/mnt/clustre/centos7users/protein/Metabonomics', #代谢项目
        '/mnt/clustre/centos7users/protein/Proteomics/', #蛋白项目
        '/mnt/clustre/centos7users/protein/Metabonomics/GCMS/', #GC项目
    ]
    year = MJ[2: 6]
    month = MJ[6: 8]
    MJ_DIR = ''
    for dir_ in raw_dirs:
        dir_ = os.path.join(dir_, year, year + month)
        if os.path.exists(dir_):
            for md in os.listdir(dir_):
                if MJ in md:
                    MJ_DIR = os.path.join(dir_, md)
                    if not PM:
                        return MJ_DIR
                    for pm in os.listdir(MJ_DIR):
                        if PM in pm:
                            return os.path.join(MJ_DIR, pm)
    if not MJ_DIR:
        raise not_find_MJ("找不到合同号对应的路径")


def find_sample_rawdatas(MJ, PM='', samples=[]):
    """
    :param MJ: 美吉的合同编号
    :param PM: 如果没有便为空
    :param samples: 样本列表，改名前的
    :return: 返回值有两个，一个是项目的绝对路径，另一个是一个pandas表格，第一列是样本名，第二列是样本对应的绝对路径
    """
    mj_dir = get_project_dir(MJ, PM)
    all_files = list()
    for root, dirs, files in os.walk(mj_dir):
        for file in files:
            all_files.append(os.path.join(root, file))
    all_files = pd.Series(all_files)
    match_df = pd.DataFrame()
    match_df['sample'] = samples
    def match(sample):
        paths = all_files.loc[all_files.str.contains('_%s\.' % sample)]
        s_num2path = dict()
        for p in paths:
            if not 'GCMS' in p:
                s_num = len(os.path.basename(p).split('_'))
                if s_num not in s_num2path:
                    s_num2path[s_num] = set()
                s_num2path[s_num].add(p)
            else:
                p = p.split('/')
                for n, d in enumerate(p):
                    if d.endswith('.D'):
                        s_num = len(d.split('_'))
                        if s_num not in s_num2path:
                            s_num2path[s_num] = set()
                        s_num2path[s_num].add('/'.join(p[0: n + 1]))
                        break
                else:
                    if 1 not in s_num2path:
                        s_num2path[1] = set()
                    s_num2path[1].add('/'.join(p))
        if not s_num2path:
            return ''
        min_ = min(s_num2path.keys())
        return ';'.join(list(s_num2path[min_]))
    match_df['paths'] = match_df['sample'].apply(match)
    return mj_dir, match_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="通过合同号和样本列表找到样本对应的原始文件路径")
    parser.add_argument("-MJ", type=str, required=True, help="项目的合同号")
    parser.add_argument("-PM", type=str, default='', help="项目的任务号，不区分的话填空，默认为空")
    parser.add_argument("-samples", type=str, default='', help="样本列表，以分号分隔")
    parser.add_argument("-out", type=str, default='', help="输出文件内容")
    parser.add_argument("-out_dir", type=str, default='', help="项目对应的文件夹")
    args = parser.parse_args()

    samples = args.samples.split(';')
    dir_, df = find_sample_rawdatas(args.MJ, args.PM, samples)
    print(dir_)
    # df.to_csv('%s.xls' % args.MJ, sep='\t', index=False)
    df.to_csv(args.out, sep='\t', index=False)
    with open(args.out_dir, "wb") as w:
        w.write(dir_ + "\n")
