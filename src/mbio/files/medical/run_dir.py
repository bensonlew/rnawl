# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# creat at 20171111

from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import os


class RunDirFile(Directory):
    '''
    base call file：下机数据文件夹<run folder>/Data/Intensities/BaseCalls/
    '''
    def __init__(self):
        super(RunDirFile, self).__init__()

    def get_info(self):
        '''
        获取文件属性
        '''
        super(RunDirFile, self).get_info()

    def check(self):
        '''
        检查文件格式
        '''
        if super(RunDirFile, self).check():
            self.get_info()
            # flist = os.listdir(self.prop['path'])
            # if len(flist) == 0:
            #     raise FileError("文件夹为空")
            # if not os.path.isdir(self.prop['path'] + "/Data/Intensities/BaseCalls/"):
            #     raise FileError("BaseCalls 文件夹不存在")

            # file_numbers = ['330', '324', '322', '320', '294', '256', '177', '172', '171', '168', '167', '166', '146',
            #                 '122', '84']       # 从以往下机BaseCalls文件夹L001 - L004中统计出来的文件个数
            # basecall = self.prop['path'] + "/Data/Intensities/BaseCalls/"
            # for i in ["L001", "L002", "L003", "L004"]:
            #     lane_dir = os.path.join(basecall, i)
            #     if not os.path.exists(lane_dir):
            #         raise FileError("{}:文件夹不存在！".format(lane_dir))
            #     count = len(os.listdir(lane_dir))
            #     if str(count) not in file_numbers:
            #         raise FileError("{}文件夹下文件数为{}，与预先设定的文件数({})"
            #                         "不符".format(i, str(count), '|'.join(file_numbers)))
        else:
            raise FileError("文件格式错误")
        return True
