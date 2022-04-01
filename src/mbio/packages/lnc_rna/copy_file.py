# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import sys
import os
import shutil
import glob

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
            if os.path.exists(olddir):
                self.linkfile(olddir, newdir, mode=mode)
            else:
                print '需要移动到output目录的文件夹不存在{}'.format(olddir)
        else:
            if os.path.exists(newdir):
                pass
                # shutil.rmtree(newdir):
            else:
                os.makedirs(newdir)
            allfiles = os.listdir(olddir)
            oldfiles = [os.path.join(olddir, i) for i in allfiles]
            newfiles = [os.path.join(newdir, i) for i in allfiles]
            for i in range(len(allfiles)):
                if os.path.isfile(oldfiles[i]):
                    self.linkfile(oldfiles[i], newfiles[i], mode=mode)
                else:
                    self.linkdir(oldfiles[i], newfiles[i], mode=mode)

    def linkdir_test(self, olddir, newdir, mode='link'):
        """
        移动一个目录到另一个目录
        ; olddir 需要移动的目录参数
        ；newdir 需要移动的目的位置
        """
        if not os.path.isdir(olddir):
            if os.path.exists(olddir):
                self.linkfile(olddir, newdir, mode=mode)
            else:
                print '需要移动到output目录的文件夹不存在{}'.format(olddir)
        else:
            if os.path.exists(newdir):
                pass
                # shutil.rmtree(newdir):
            else:
                os.makedirs(newdir)
            allfiles = os.listdir(olddir)
            oldfiles = [os.path.join(olddir, i) for i in allfiles]
            if os.path.basename(olddir) in ['diffexpress_t', 'exp_pca_t']:
                newfiles = [os.path.join(newdir, 'T_'+i) for i in allfiles]
            elif os.path.basename(olddir) in ['diffexpress', 'exp_pca']:
                newfiles = [os.path.join(newdir, 'G_' + i) for i in allfiles]
            else:
                newfiles = [os.path.join(newdir, i) for i in allfiles]
            for i in range(len(allfiles)):
                if os.path.isfile(oldfiles[i]):
                    self.linkfile(oldfiles[i], newfiles[i], mode=mode)
                else:
                    self.linkdir(oldfiles[i], newfiles[i], mode=mode)

    def linkfile(self, oldname, newname, mode='link'):
        """
        移动一个目文件到另一个文件
        ; oldname 需要移动的目录参数
        ；newname 需要移动的目的位置
        """

        if not os.path.exists(oldname):
            print '需要移动到output目录的文件夹不存在{}'.format(oldname)
        else:
            if os.path.exists(newname):
                os.remove(newname)
            if os.path.exists(os.path.dirname(newname)):
                pass
            else:
                os.makedirs(os.path.dirname(newname))
            if mode == "link":
                os.link(oldname, newname)
            else:
                os.system("cp {} {}".format(oldname, newname))

    def renameafile(self, old, new):
        print "rename"
        print old
        print new
        old_dir = os.path.dirname(old)
        new_dir = os.path.dirname(new)
        if os.path.exists(new_dir):
            os.rename(old, new)
        else:
            if new_dir == "":
                print "new dir empty"
            else:
                os.makedirs(new_dir)
                os.rename(old, new)

    def renamefile(self, files, rename):
        """
        移动文件，重命名文件
        ; files 文件通配符
        ；rename 两个元素的列表，表示替换的字符串
        """
        if type(rename) == list:
            file_list = glob.glob(files)
            for afile in file_list:
                print rename
                print afile

                new_file = os.path.basename(afile).replace(rename[0], rename[1])
                self.renameafile(afile, os.path.join(os.path.dirname(afile), new_file))

    def remove(self, files):
        """
        删除文件
        ; files 需要删除的文件通配符
        """
        file_list = glob.glob(files)
        for afile in file_list:
            os.remove(afile)
