# -*- coding: utf-8 -*-
# __author__ = 'hongdong.xuan'

"""用于拆分表split_check的检查主要是检查barcode"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError
import re


class SplitCheckFile(File):
    """
    SplitCheck类
    """

    def __init__(self):
        super(SplitCheckFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(SplitCheckFile, self).get_info()
        samples = self.get_sample_list()
        self.set_property("wq_sample_list", samples[0])
        self.set_property("ws_sample_list", samples[1])

    def get_sample_list(self):
        wq_sample_list = []
        ws_sample_list = []
        with open(self.prop['path'], 'r') as r:
            for line in r:
                line = line.strip().split('\t')
                results = re.match('WQ([0-9]*)-.*', line[3])
                results_ws = re.match('WS(.*)', line[3])
                if results:
                    wq_sample_list.append('WQ' + results.group(1))
                elif results_ws:
                    ws_sample_list.append(results_ws.group(0))
                else:
                    pass
        return list(set(wq_sample_list)), list(set(ws_sample_list))

    def sample_type_check(self):
        with open(self.prop['path'], 'r') as r:
            for line in r:
                line = line.strip().split('\t')
                results = re.match('WQ([0-9]*)-.*', line[3])
                results_ws = re.match('WS(.*)', line[3])
                if results:
                    if 'F' in line[3] or 'M' in line[3]:
                        if line[2] != 'dcpt':
                            raise FileError('父本与母本样本{}类型不正确，应该是dcpt!'.format(line[3]))
                    elif 'S' in line[3]:
                        if line[2] != 'pt':
                            raise FileError('胎儿样本{}类型不正确，应该是dcpt!'.format(line[3]))
                if results_ws:
                    if line[2] != 'nipt':
                        raise FileError('产筛样本{}类型不正确，应该是nipt!'.format(line[3]))
        return True

    def get_barcode_info(self):
        """
        获取样本的barcode，并进行检查
        :return:
        """
        barcode_list = []
        with open(self.prop['path'], 'r') as r:
            for line in r:
                line = line.strip().split('\t')


    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(SplitCheckFile, self).check():
            # self.sample_type_check()
            return True
        else:
            raise FileError("文件格式错误")

if __name__ == '__main__':
    path = "/mnt/ilustre/tsanger-data/rerewrweset/files/m_5950/20171016/170928_TPNB500180_0112_AHMWJKAFXX_1508133981.xls"
    a = SplitCheckFile()
    a.set_path(path)
    a.check()
