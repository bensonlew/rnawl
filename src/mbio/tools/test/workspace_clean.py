# -*- coding: utf-8 -*-
# __author__ = 'shijin'

import os

for dir in os.listdir("/mnt/ilustre/users/sanger/workspace"):
    if dir == "tmp":
        continue
    dir_num = int(dir)
    if dir_num <= 20161130:
        os.system("mv /mnt/ilustre/users/sanger/workspace/{}  /mnt/ilustre/users/sanger/bak".format(dir))