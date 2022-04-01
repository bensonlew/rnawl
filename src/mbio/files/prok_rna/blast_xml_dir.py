# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.core.exceptions import FileError
from biocluster.iofile import Directory
from mbio.files.align.blast.blast_xml import BlastXmlFile
import os


class BlastXmlDirFile(Directory):
    """
    bam文件夹格式
    """
    def __init__(self):
        super(BlastXmlDirFile, self).__init__()

    def check(self):
        if super(BlastXmlDirFile, self).check():
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
            blastxml = BlastXmlFile()
            path = os.path.join(self.path, f)
            blastxml.set_path(path)
            blastxml.check()
            full_files.append(path)
            files_obj.append(blastxml)
            basenames.append(base_name)
        self.set_property('files_num', len(files))  # 文件数量
        self.set_property('files', full_files)  # 文件路径
        self.set_property("file_objs", files_obj)  # 文件对象
        self.set_property("basenames", basenames)  # 文件名去后缀

    def convert2table(self, out_dir):
        outfiles = []
        for f, name in zip(self.prop['file_objs'], self.prop['basenames']):
            out_file = os.path.join(out_dir, name + '.xls')
            f.convert2table(out_file)
            outfiles.append(out_file)
        return outfiles

    def convert2blast6default(self, out_dir):
        outfiles = []
        for f, name in zip(self.prop['file_objs'], self.prop['basenames']):
            out_file = os.path.join(out_dir, name + '.xls')
            f.convert2blast6default(out_file)
            outfiles.append(out_file)
        return outfiles

if __name__ == '__main__':
    f_dir = "C:\\Users\\sheng.he.MAJORBIO\\Desktop\\kegg_diamond"
    a = BlastXmlDirFile()
    a.set_path(f_dir)
    a.check()
    a.convert2blast6default(f_dir)
