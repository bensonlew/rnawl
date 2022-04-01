# -*- coding: utf-8 -*-

import os
import sys
from concurrent.futures import ThreadPoolExecutor

if len(sys.argv) != 4:
    exit('USAGE: %s stat_dir count_dir output' %sys.argv[0])

stat_dir = sys.argv[1]
count_dir = sys.argv[2]
out = sys.argv[3]
all_srcipts = list()

if not os.path.exists(count_dir):
    exit('the %s not exists'% count_dir)

def steal_user_scripts(user):
    scripts_list = list()
    dir_ = os.path.join(out, user)
    print(dir_)
    if not os.path.exists(dir_):
        os.mkdir(dir_)
    with open(os.path.join(stat_dir, user+'_file_counts.txt'),'r') as u_r:
        _ = u_r.readline()
        n=0
        for line in u_r.readlines():
            line = line.strip().split('\t')
            if len(line) > 1:
                file = line[0]
                if file.endswith('.py') or file.endswith('.pl') or file.endswith('.sh'):
                    if not u'__' in file and not 'qsub' in file and not 'MJ' in file and not os.path.basename(file) in all_srcipts:
                        n += 1
                        if not os.path.basename(file) in all_srcipts:
                            scripts_list.append(os.path.basename(file))
                            all_srcipts.append(os.path.basename(file))
                            target = os.path.join(dir_, os.path.basename(file))
                            try:
                                os.link(file, target)
                            except:
                                pass
                        else:
                            target = os.path.join(dir_, str(n) + os.path.basename(file))
                            try:
                                os.link(file, target)
                            except:
                                pass


users = [x for x in os.listdir(count_dir) if os.path.isdir(os.path.join(count_dir, x))]

with ThreadPoolExecutor(35) as pool:
    pool.map(steal_user_scripts, users)

