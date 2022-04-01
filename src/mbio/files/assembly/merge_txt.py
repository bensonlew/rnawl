# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""txt格式文件类"""

from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class MergeTxtFile(File):
    """
    txt类,所有需要merge的所有gtf文件的绝对路径
    """
    def __init__(self):
        super(MergeTxtFile, self).__init__()
        
    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(MergeTxtFile, self).check():
            super(MergeTxtFile, self).get_info()
            with open(self.prop["path"], "r") as f:
                for line in f:
                    if line.startswith("/"):
                        if line.endswith(".gtf\n"):
                            return True
                        else:
                            raise FileError("%s行格式错误，必须是gtf文件", variables = (line), code = "41300201")
                    else:
                        raise FileError("文件格式错误，文本路径为gtf文件的绝对路径", code = "41300202")
        else:
            raise FileError("文件格式错误", code = "41300203")
# if __name__ == "__main__":
#     a = MergeTxtFile()
#     # a.set_path('/mnt/ilustre/users/sanger-dev/workspace/20170316/Single_assembly_module_tophat_stringtie_gene1/Assembly/assembly_gtf.txt')
#     # a.set_path('tielog_1.txt')
#     a.check()