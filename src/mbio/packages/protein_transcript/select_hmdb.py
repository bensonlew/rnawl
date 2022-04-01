#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/8 9:43
@file    : select_hmdb.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


def select_hmdb(hmdb, filters, out):
    with open(hmdb, 'rb') as hr,\
    open(out, 'wb') as ow:
        filters = filters.split(';')
        while '' in filters:
            filters.remove('')
        ow.write(hr.readline())
        for line in hr:
            for filter in filters:
                if filter.encode('utf-8') not in line:
                    break
            else:
                ow.write(line)
                # ow.write(line.encode('utf-8'))


if __name__ == '__main__':
    filters = input("请输入你想要筛选的字段，如果有多个，字段之间要用';'分隔\nfor example:Saliva;Blood;Urine;Endogenous;Glycine max\n")
    out = input("请输入筛选后生成文件的名字：")
    if not out:
        out = 'selected.xls'
    hmdb = 'Y:\\蛋白代谢事业部\\蛋白与代谢生信分析部\\hmdb\\hmdb_metabolites.xls'
    select_hmdb(hmdb, filters, out)