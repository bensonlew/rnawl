## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
# 20180719

import argparse
import os
from collections import defaultdict
from collections import OrderedDict
import sys
import pandas as pd
from biocluster.config import Config
from mako.template import Template


if len(sys.argv) != 6:
    exit('USAGE: python %s group_list trimPairFq.list special memory' % sys.argv[0])
group_list = sys.argv[1]
trim_list = sys.argv[2]
special = sys.argv[3]
memory = sys.argv[4]
work_dir = sys.argv[5]

s2g = defaultdict(list)
java_path = Config().SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin/java'
rock_path = Config().SOFTWARE_DIR + '/bioinfo/rna_pro/Rockhopper.jar'

group_f = pd.read_table(group_list, sep='\t', header = 0)
samples = group_f['#sample'].tolist()
samples = [str(x) for x in samples]
groups = group_f['group'].tolist()
groups = [str(x) for x in groups]
g2str = OrderedDict()
for z in zip(samples, groups):
    s2g[z[1]].append(z[0])

with open(trim_list, 'r') as trim_r:
    for line in trim_r.readlines():
        if line.strip():
            line = line.strip().split('\t')
            sample = line[0]
            for group in s2g:
                if sample in s2g[group]:
                    if not group in g2str:
                        g2str[group] = ''
                        g2str[group] += line[1] + '%' + line[2] + ','
                    else:
                        g2str[group] += line[1] + '%' + line[2] + ','

group_str = ','.join(g2str.keys())
trims = [x.strip(',') for x in g2str.values()]
trim_str = ' '.join(trims)

index_str = ''
for index in os.listdir('rock_index'):
    ptt_path = os.path.join(work_dir, 'rock_index/{}/{}.ptt'.format(index,index))
    if len(open(ptt_path,'r').readlines()) > 3:
        index_str += 'rock_index/' + index + ','

if special.lower() == 'true':
    bash_info = u"""${java} -Xms${memory}000m -Xmx${memory}000m -cp ${rock_path} Rockhopper -s true -rf -p 20 -z 0.2 -g ${index_str} ${trim_str} -L ${group_str} -v true
"""
else:
    bash_info = u"""${java} -Xms${memory}000m -Xmx${memory}000m -cp ${rock_path} Rockhopper -s true -p 20 -z 0.2 -g ${index_str} ${trim_str} -L ${group_str} -v true
"""

f = Template(bash_info)
bash_info = f.render(rock_path=rock_path,
                     index_str=index_str.strip(','),
                     trim_str=trim_str,
                     group_str=group_str,
                     java = java_path,
                     memory = memory
                     )

with open('run_rockhopper.bash', 'w') as rock_sh:
    rock_sh.write(bash_info)

os.system('bash run_rockhopper.bash')