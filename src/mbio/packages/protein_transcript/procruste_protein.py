#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/1 14:35
@file    : procruste_protein.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

# !/usr/bin/env python

__contact__ = 'Edward-Wang@foxmail.com'

import sys
sys.path.insert(1, '/mnt/ilustre/users/hongming.wang/scripts/Procrustes')
from treatData import *
import pandas as pd
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import os
# import rpy2.robjects as robj
# import rpy2.robjects.pandas2ri
# from plot import procrustes_points


def twoDimPic(map_file, group_header):
    reference, trans, importance, merged = get_rpy_input('procrustes_results/pcoa1_transformed_reference.txt',
                                                         'procrustes_results/pcoa2_transformed_q1.txt', map_file)
    p_value = get_p_value('procrustes_results/procrustes_results.txt')
    p_value = float(p_value)
    importance = map(lambda x: str(round(x * 100, 2)), importance)
    importance = ','.join(importance)
    merged.to_csv('merge', index=True, sep='\t')
    r_cmd = Template(filename='/mnt/ilustre/users/ting.kuang/ALL-SCRIPT/procruste.r')
    r_cmd = r_cmd.render(importance=importance,
                         pvalue=p_value,
                         group_header=group_header,
                         )
    r_path = 'procruste.r'
    with open(r_path, 'w') as rw:
        rw.write(r_cmd)
    os.system('Rscript procruste.r')



def main(distance1, distance2, map_file=None, group_header=None, ThreeDim=False):
    perform_pcoa(distance1, distance2)
    perform_proc('pcoa1.txt', 'pcoa2.txt')
    if ThreeDim:
        perform_emperer('procrustes_results', map_file)
    else:
        twoDimPic(map_file, group_header)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Perform Procrustes Analysis!')
    parser.add_argument('-d1', '--distance1', required=True, help='The first distance matrix')
    parser.add_argument('-d2', '--distance2', required=True, help='The second distance matrix')
    parser.add_argument('-m', '--map', default=None, help="map file; header should be '#SampleID\tGroup'")
    parser.add_argument('-g', '--group', default=None, help="Header in map file such as Group")
    parser.add_argument('-TD', '--ThreeDimentions', action='store_true',
                        help='add this if you want to draw 3 dimensional figure. map file should be provided.')
    args = parser.parse_args()

    main(args.distance1, args.distance2, map_file=args.map, group_header=args.group, ThreeDim=args.ThreeDimentions)
