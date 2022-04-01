# -*- coding: utf-8 -*-
import sys
import os

"""
tmp = sys.argv[1]
print tmp
print "***********************"
tmp = sys.argv[2]
print tmp
print "***********************"
merge_file = sys.argv[3]
print merge_file
print "***********************"
dirs = []
dirs.append(sys.argv[1])
dirs.append(sys.argv[2])
filt = []
for path in dirs:
    print path
    if os.path.exists(path):
        with open(path, "rb") as f:
            # lines = f.readlines()
            for line in f:
                if line not in filt:
                    filt.append(line)
    else:
        raise Exception("{}文件不存在".format(path))
with open(merge_file, "wb") as w:
    for line in filt:
        w.write(line)
        """
        
os.system("cat {} {} > {}".format(sys.argv[1], sys.argv[2], sys.argv[3]))