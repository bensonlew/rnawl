#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = "fengyitong 2018-12-20"

from __future__ import print_function
import sys
import os
import time
import random

abelt = dict(
    a='''

     /\\
    /  \\
   / /\ \\
  / ____ \\
 /_/    \_\\
           ''',
    b='''
  ____
 |  _ \\
 | |_) |
 |  _ <
 | |_) |
 |____/
        ''',
    c='''
    _____
  / ____|
 | |
 | |
 | |____
  \_____|
         ''',
    d='''
  _____
 |  __ \\
 | |  | |
 | |  | |
 | |__| |
 |_____/
         ''',
    e='''
  ______ 
 |  ____|
 | |__   
 |  __|  
 | |____ 
 |______|
        ''',
    f='''
  ______
 |  ____|
 | |__
 |  __|
 | |
 |_|
         ''',
    g='''
    _____
  / ____|
 | |  __
 | | |_ |
 | |__| |
  \_____|
         ''',
    h='''
  _    _
 | |  | |
 | |__| |
 |  __  |
 | |  | |
 |_|  |_|
         ''',
    i='''
  _____
 |_   _|
   | |
   | |
  _| |_
 |_____|
         ''',
    j='''
       _
      | |
      | |
  _   | |
 | |__| |
  \____/
         ''',
    k='''
  _  __
 | |/ /
 | ' /
 |  <
 | . \\
 |_|\_\\
       ''',
    l='''
  _
 | |
 | |
 | |
 | |
 |_|
       ''',
    m='''
  __  __
 |  \/  |
 | \  / |
 | |\/| |
 | |  | |
 |_|  |_|
         ''',
    n='''
  _   _
 | \ | |
 |  \| |
 | . ` |
 | |\  |
 |_| \_|
        ''',
    o='''
   ____
  / __ \\
 | |  | |
 | |  | |
 | |__| |
  \____/
         ''',
    p='''
  _____
 |  __ \\
 | |__) |
 |  ___/
 | |
 |_|
         ''',
    q='''
   ____
  / __ \\
 | |  | |
 | |  | |
 | |__| |
  \___\_\\
         ''',
    r='''
  _____
 |  __ \\
 | |__) |
 |  _  /
 | | \ \\
 |_|  \_\\
         ''',
    s='''
   _____
  / ____|
 | (___
  \___ \\
  ____) |
 |_____/
         ''',
    t='''
  _______
 |__   __|
    | |
    | |
    | |
    |_|
          ''',
    u='''
  _    _
 | |  | |
 | |  | |
 | |  | |
 | |__| |
  \____/
         ''',
    v='''
 __      __
 \ \    / /
  \ \  / / 
   \ \/ /  
    \  /   
     \/    
           ''',
    w='''
 __          __
 \ \        / /
  \ \  /\  / /
   \ \/  \/ /
    \  /\  /
     \/  \/
               ''',
    x='''
 __   __
 \ \ / /
  \ V /
   > <
  / . \\
 /_/ \_\\
        ''',
    y='''
 __     __
 \ \   / /
  \ \_/ /
   \   /
    | |
    |_|
          ''',
    z='''
  ______
 |___  /
    / /
   / /
  / /__
 /_____|
        '''
)

info_type = dict(
    info = 'INFO :  {}运行成功',
    warning1= 'WARNING : {}的参数没有默认值，请确认确实不需要默认值？',
    warning2= 'WARNING : *************************************************************************************还有子对象没有运行完成,请确认此情况是否程序本意？如非故意如此, 此问题可能会导致后续程序抛出EventStopError异常或严重的逻辑错误！请注意on_rely对象为数组时，只对on_rely语句执行时刻的已有数组成员有效，后期再添加的数组成员不能自动成为Rely对象!************************************************************************************',
    debug1='DEBUG : {}开始运行',
    debug2='DEBUG : {}正在运行中，流程保持监控',
    debug3='DEBUG : {}已经运行完成，正在检查并开启下一步骤的运行',
    error='ERROR : {}出现未知错误，运行失败！！！'
)

def color_print(info):
    if u'INFO' in info:
        print('\033[1;32;40m', end='')
        print(str(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) + '\t' + info)
        print('\033[0m', end='')
    if u'DEBUG' in info:
        print('\033[0;34;40m', end='')
        print(str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))) + '\t' + info)
        print('\033[0m', end='')
    if u'WARNING' in info:
        print('\033[4;33;40m', end='')
        print(str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))) + '\t' + info)
        print('\033[0m', end='')
    if u'ERROR' in info:
        print('\033[7;35;40m', end='')
        print(str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))) + '\t' + info)
        print('\033[0m', end='')

def print_steps(step, error =0):
    color_print(info_type['debug1'].format(step))
    # time.sleep(random.randint(0, 2))
    time.sleep(random.uniform(0,0.1))
    if random.random() < 0.8:
        color_print(info_type['warning1'].format(step))
        # time.sleep(random.randint(0, 1))
        time.sleep(random.uniform(0,0.1))
    color_print(info_type['debug2'].format(step))
    if random.random() < 0.1:
        color_print(info_type['warning2'].format(step))
        # time.sleep(random.randint(0, 3))
        time.sleep(random.uniform(0,0.1))
    color_print(info_type['debug3'].format(step))
    # time.sleep(random.randint(0, 1))
    time.sleep(random.uniform(0,0.1))
    if error:
        color_print(info_type['error'].format(step))
    else:
        color_print(info_type['info'].format(step))
    time.sleep(0.3)


def print_flower_letter(letter,success = 1):
    if success:
        print('\033[5;32;40m')
    else:
        print('\033[5;31;40m')
    print_info = [abelt[x] for x in letter.lower()]
    length_line = [len(x.strip().split('\n')) for x in print_info]
    for n in range(max(length_line) + 1):
        info = ' '
        for i in print_info:
            i = i.strip('\n').split('\n')
            try:
                info += i[n].replace('\n', ' ') + ' ' * (
                            len(i[-1].replace('\n', ' ')) - len(i[n].replace('\n', ' '))) + ' '
            except:
                info += ' ' * len(i[-1])
        print(info)
    print('\033[0m')

first_p = ['Prokrna工作流    INFO : 开始更新步骤信息...',
'Prokrna工作流    DEBUG : 初始化RPC服务器..!',
'FileTransfer    DEBUG : 开始启动新进程20178 ....',
'FileTransfer    DEBUG : 开始传输文件',
'FileTransfer    DEBUG : 所有传输线程结束请队列为空，进程20178等待30秒 ....',
'FileTransfer    INFO : 所有文件传输完成']

steps_list = ["rock_index", "filecheck", "rna_qc", "mapping", "express", "diffexpress", "snp_rna",
                     "map_qc", "annot_mapdb", "annot_orfpfam", "annot_filter",  "annot_class", "exp_pca",
                     "exp_corr", "exp_venn", "qc_stat_before", "qc_stat_after", "rockhopper", "srna", "promote"
                     ]

t = 0.1
for p in first_p:
    color_print(p)
    t += 0.1
    time.sleep(t)

for step in steps_list:
    print_steps(step)

print_steps('output', error=1)

print_flower_letter('error', success=0)
# print_flower_letter('failed', success=0)

