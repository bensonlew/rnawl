# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import Directory
import os
from biocluster.core.exceptions import FileError
from mbio.files.align.blast.blast_table import BlastTableFile


class BlastTableDirFile(Directory):
    """
    bam文件夹格式
    """
    def __init__(self):
        super(BlastTableDirFile, self).__init__()

    def check(self):
        if super(BlastTableDirFile, self).check():
            self.get_info()
            return True


    def get_info(self):
        files = os.listdir(self.path)
        full_files = []
        files_obj = []
        basenames = []
        if not len(files):
            raise FileError('文件夹为空，请检查！')
        for f in files:
            base_name = os.path.splitext(f)[0]
            blasttable = BlastTableFile()
            path = os.path.join(self.path, f)
            blasttable.set_path(path)
            blasttable.check()
            full_files.append(path)
            files_obj.append(blasttable)
            basenames.append(base_name)
        self.set_property('files_num', len(files))  # 文件数量
        self.set_property('files', full_files)  # 文件路径
        self.set_property("file_objs", files_obj)  # 文件对象
        self.set_property("basenames", basenames)  # 文件名去后缀
