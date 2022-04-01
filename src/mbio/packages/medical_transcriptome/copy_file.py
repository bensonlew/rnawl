# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import sys
import os
import shutil

class CopyFile(object):
    def __init__(self):
        pass

    def linkdir(self, olddir, newdir, mode='link'):
        """
        移动一个目录到另一个目录
        ; olddir 需要移动的目录参数
        ；newdir 需要移动的目的位置
        """
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在{}'.format(olddir))
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            else:
                self.linkdir(oldfiles[i], newfiles[i])

    def linkfile(self, oldname, newname, mode='link'):
        """
        移动一个目文件到另一个文件
        ; oldname 需要移动的目录参数
        ；newname 需要移动的目的位置
        """
        if not os.path.exists(oldname):
            raise Exception('需要移动到output文件不存在{}'.format(oldname))
        if os.path.exists(newname):
            os.remove(newname)
        os.link(oldname, newname)
