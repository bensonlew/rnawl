# -*- coding: utf-8 -*-
# __author__ = 'xuan.zhou'

"""宏基因非冗余基因集生成的reads_profile.txt"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError
from biocluster.core.exceptions import FileError


class ProfileTableFile(File):
    """
    txt类
    """

    def __init__(self):
        super(ProfileTableFile, self).__init__()

    #
    # def check(self):
    # 	"""
    # 	检测文件是否为空，数据是否存在
    # 	:return:
    # 	"""
    # 	if super(ProfileTableFile, self).check():
    # 		return True
    # 	else:
    # 		raise FileError("文件格式错误")

    def check(self):
        with open(self.prop['path'], 'r') as r:
            n = 0
            number = 0
            for line in r:
                n += 1
                if n == 1:
                    number = len(line.rstrip('\r\n').split("\t"))
                else:
                    if len(line.rstrip('\r\n').split("\t")) != number:
                        raise FileError('文件中信息不全，请检查', code="44001401")
        return True

    #def get_top(self, n, top_header):
        """
        get top n profile
        :param n: top profile number
        :param top_header: header named used to calculated
        :return: change file with its top n
        """
        #file_path = self.prop['path']
        #tmp_file = self.prop['path'] + '.tmp'
        #data = pd.read_table(file_path, )

