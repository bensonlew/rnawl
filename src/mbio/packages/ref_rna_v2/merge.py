# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import sys
import os

class Merge(object):
    def __init__(self):
        pass
    def merge2file(self, file1, file2, fileto, withhead = True):
        '''
        合并两个文件 cat
        file1 : 文件1路径
        file2： 文件2路径
        fileto: 结果文件路径
        '''
        print "**{} **{} **{}".format(file1, file2, fileto)
        if withhead == True:
            os.system("cat {} > {}".format(file1, fileto))
            os.system("sed '1d' {} >> {}".format(file2, fileto))
        else:
            os.system("cat {} {} > {}".format(file1, file2, fileto))
