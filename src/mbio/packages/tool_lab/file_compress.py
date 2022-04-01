# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import gzip
import os
import tarfile
import zipfile
import rarfile
import subprocess
from biocluster.config import Config

class FileCompress(object):

    def __init__(self,uncompress,out_dir):
        self.ucomp = uncompress
        self.out_dir = out_dir


    def uncompress_file(self,uncompress,out_dir):
        if uncompress.endswith("gz"):
            self.un_gz(uncompress,out_dir)
        elif uncompress.endswith("tgz") or uncompress.endswith("tar.gz"):
            self.un_tgz(uncompress,out_dir)
        elif uncompress.endswith("rar"):
            self.un_rar(uncompress,out_dir)
        elif uncompress.endswith("tar"):
            self.un_tar(uncompress,out_dir)
        elif uncompress.endswith("un_zip"):
            self.un_zip(uncompress, out_dir)
        else:
            raise Exception("不可以用这样子的压缩文件")



    # gz
    # 因为gz一般仅仅压缩一个文件，全部常与其它打包工具一起工作。比方能够先用tar打包为XXX.tar,然后在压缩为XXX.tar.gz
    # 解压gz，事实上就是读出当中的单一文件


    def un_gz(self,file_name,out_dir):
        """ungz zip file"""
        f_name = os.path.join(out_dir,file_name.replace(".gz", ""))
        #获取文件的名称，去掉
        g_file = gzip.GzipFile(file_name)
        #创建gzip对象
        open(f_name, "w+").write(g_file.read())
        #gzip对象用read()打开后，写入open()建立的文件里。
        g_file.close()
        #关闭gzip对象

    def un_tgz(self, file_name, out_dir):
        jobs = "tar -xzvf {} -C {}".format(file_name,out_dir)
        subprocess.call(jobs, shell=True)

    def un_rar(self,file_name,out_dir):
        rar_path = Config().SOFTWARE_DIR + "/program/rar/rar/"
        jobs = "{}rar x {} -C {}".format(rar_path, file_name, out_dir)
        subprocess.call(jobs, shell=True)


    # tar
    # XXX.tar.gz解压后得到XXX.tar，还要进一步解压出来。
    # 注：tgz与tar.gz是同样的格式，老版本号DOS扩展名最多三个字符，故用tgz表示。
    # 因为这里有多个文件，我们先读取全部文件名称。然后解压。例如以下：
    # 注：tgz文件与tar文件同样的解压方法。
    def un_tar(self,file_name,out_dir):
           # untar zip file"""
        tar = tarfile.open(file_name)
        names = tar.getnames()
        if os.path.isdir(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        #因为解压后是很多文件，预先建立同名目录
        for name in names:
            tar.extract(name, file_name + "_files/")
        tar.close()

    # zip
    # 与tar类似，先读取多个文件名称，然后解压。例如以下：

    def un_zip(self,file_name,out_dir):
        """unzip zip file"""
        zip_file = zipfile.ZipFile(file_name)
        if os.path.isdir(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        for names in zip_file.namelist():
            zip_file.extract(names,out_dir)
        zip_file.close()


    # rar
    # 由于rar通常为window下使用，须要额外的Python包rarfile。
    #
    # 可用地址： http://sourceforge.net/projects/rarfile.berlios/files/rarfile-2.4.tar.gz/download
    #
    # 解压到Python安装文件夹的/Scripts/文件夹下，在当前窗体打开命令行,
    #
    # 输入Python setup.py install
    #
    # 安装完毕。

    # def un_rar(file_name):
    #     """unrar zip file"""
    #     rar = rarfile.RarFile(file_name)#待解压文件
    #     if os.path.isdir(file_name + "_files"):
    #         pass
    #     else:
    #         os.mkdir(file_name + "_files")
    #     os.chdir(file_name + "_files")
    #     rar.extractall()#解压指定目录
    #     rar.close()

