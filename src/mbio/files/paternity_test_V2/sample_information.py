# -*- coding: utf-8 -*-
# __author__ = 'chenhongyu'

'''PT_NIPT 生产表、反馈表内容检查'''

from biocluster.iofile import File
import subprocess
from biocluster.config import Config
import os
import re
from biocluster.core.exceptions import FileError
from collections import defaultdict
import codecs
from chardet.universaldetector import UniversalDetector

class SampleInformationFile(File):
    '''

    '''
    def __init__(self):
        '''

        '''
        super(SampleInformationFile, self).__init__()
        self.header_split = list()
        self.content_split_array = list()
        self.table_reformed = list()



        #########################################定义一些常量############################################################
        # 建库类型缩写
        self.libraryTypePre = ["FA", "NIPT", "WQ",  "YCZL", "HLY", "QA", "CT", "XA"]
        # libraryType列表为上机表中的所有的建库类型
        self.libraryType = [u"NIPT文库构建实验流程", u"WQ gDNA建库实验流程_多重", u"WQ cfDNA建库实验流程_杂捕", u"FA gDNA建库实验流程_杂捕", u"FA ctDNA建库实验流程_杂捕",
                       u"YCZL 建库实验流程", u"QA 建库实验流程", u"CT 建库实验流程", u"HLY 建库实验流程", u"FA gDNA建库实验流程", u"FA cfDNA建库实验流程_杂捕",
                       u"RXA 建库实验流程", u"snp 建库实验流程", u"XA RNA建库实验流程", u"small RNA建库实验流程", u"XA DNA建库实验流程", u"甲基化 建库实验流程"]
        # analyLibraryType列表为需要分析的测序项目的检录类型，为libraryType的子列表
        self.analyLibraryType = [u"NIPT文库构建实验流程", u"WQ gDNA建库实验流程_多重", u"WQ cfDNA建库实验流程_杂捕", u"FA gDNA建库实验流程_杂捕",
                            u"FA ctDNA建库实验流程_杂捕", u"YCZL 建库实验流程", u"QA 建库实验流程", u"CT 建库实验流程", u"HLY 建库实验流程", u"FA gDNA建库实验流程",
                            u"FA cfDNA建库实验流程_杂捕", u"XA RNA建库实验流程", u"XA DNA建库实验流程"]
        # 样本类型及其缩写
        self.snDNAType = {u"全血": "QX", u"血浆": "XJ", u"蜡块": "SL", u"石蜡": "SL", u"石蜡切片": "SL", u"石蜡切片(白片)": "SL", u"胸腹水": "XS",
                     u"手术标本": "XZ", u"穿刺标本": "XZ", u"穿刺样本": "XZ", u"组织标本": "XZ", u"新鲜组织": "XZ", u"蜡卷": "SL"}
        self.analysisType = [u"生产", u"测试", u"生产返工"]
        # 合并后的上机表的表头
        self.all_title = [u'序号', u'收样日期', u'样本姓名', u'建库类型', u'产品线', u'样本编号', u'内部订单编号', u'亲本', u'样本类型', u'文库名称',
                          u'上机测序量', u'index序列', u'片段大小', u'下机数据量', u'OT率', u'去重比率', u'去重平均深度', u'深度覆盖率', u'下机质控结果',
                          u'系统实验结果', u'分析类型', u'备注']

        ########################################定义存储错误条目的数据结构##############################################
        self.library_sn_not_match = list()
        self.library_illegal = list()
        self.ctDNA_sample_type_error = list()
        self.test_sample_without_test_prefix = list()
        self.FXLX_error = list()


    def get_info(self):
        '''
        获取文件路径、basename
        :return:
        '''
        super(SampleInformationFile, self).get_info()

    def read_table(self):
        '''
        根据传入的文件的路径，获取文件内容：表头存入self.header_split中，内容存入self.content_split_array，

        :return:
        '''

        encoding = self.detect_encoding(self.prop['path'])

        with codecs.open(self.prop['path'], 'r', encoding) as FIN:

            lines = FIN.readlines()
            if len(lines) < 1:
                raise FileError("文件为空")
            first_line_split = self.rm_extra_space(lines[0].split(','))
            if first_line_split[0] == '':                                        #如果是线下上机表（这个判断方式需要后续维护）
                if len(lines) < 4:
                    raise FileError("线下上机表内容为空")
                self.header_split = self.rm_extra_space(lines[2].rstrip('\n').split(','))
                for index in range(3, len(lines)):
                    each_line_split = lines[index].rstrip('\n').split(',')
                    self.content_split_array.append(self.rm_extra_space(each_line_split))
            else:                                                                                        #如果是线上上机表
                if len(lines) < 2:
                    raise FileError("线上上机表内容为空")
                self.header_split = self.rm_extra_space(lines[0].rstrip('\n').split(','))
                for index in range(1, len(lines)):
                    self.content_split_array.append(self.rm_extra_space(lines[index].rstrip('\n').split(',')))

            temp = list()
            for i in self.content_split_array:
                if ''.join(i) != '':
                    temp.append(i)
            self.content_split_array = temp

    def rm_extra_space(self, line_split):
        '''
        分割后数据中多余的空白
        :param line_split:
        :return:
        '''
        for i in range(0, len(line_split)):
            line_split[i] = line_split[i].strip()
        return line_split

    def get_index(self, title):
        '''
        根据列名返回对应的索引
        :param title:
        :return:
        '''
        for i in range(0, len(self.header_split)):
            if self.header_split[i] == title:
                return i

    def check(self):
        '''
        检测文件是否满足要求
        :return:
        '''
        if super(SampleInformationFile, self).check():                               #如果文件路径已设置，获取文件内容，按列名获取index
            self.get_info()                                                          #获取文件路径与basename
            if os.path.getsize(self.prop['path']) < 100*1024*1024:                   #上机表文件应该在100kb左右，这里设置了100M的限制
                self.read_table()
                self.index_library_type = self.get_index(u'建库类型')
                self.index_sn = self.get_index(u'内部订单编号')
                self.index_parent = self.get_index(u'亲本')
                self.index_sample_type = self.get_index(u'样本类型')
            else:
                raise FileError("上机表文件超过100M，请检查文件或修改程序")
        self.checker_one()
        self.checher_two()
        self.checher_three()
        self.table_reform()
        #self.set_property('header_split', self.header_split)
        #self.set_property('content_split_array', self.content_split_array)
        return True



    def checker_one(self):
        '''
        检查线上上机表和线下上机表的建库类型和样本名称是否能核对上，并检查建库类型是否合法
        :return:
        '''

        for i in range(0, len(self.content_split_array)):
            library_type = self.content_split_array[i][self.index_library_type]
            sn = self.content_split_array[i][self.index_sn]
            parent = self.content_split_array[i][self.index_parent]
            if library_type in self.analyLibraryType:
                if library_type == u'NIPT文库构建实验流程':
                    if 'WS' in sn:
                        pass
                    else:
                        self.library_sn_not_match.append('\t'.join(self.content_split_array[i]))

                if library_type == u'WQ cfDNA建库实验流程_杂捕':
                    if 'WQ' in sn and 'S' in parent:
                        pass
                    else:
                        self.library_sn_not_match.append('\t'.join(self.content_split_array[i]))

                if library_type == u'WQ gDNA建库实验流程_多重':
                    if 'WQ' in sn and ('F' in parent or 'M' in parent):
                        pass
                    else:
                        self.library_sn_not_match.append('\t'.join(self.content_split_array[i]))

            else:
                self.library_illegal.append('\t'.join(self.content_split_array[i]))
        self.set_property('library_sn_not_match', self.library_sn_not_match)         # 将捕捉到的错误添加到属性里，传递给tool
        self.set_property('library_illegal', self.library_illegal)

    def checher_two(self):
        '''
        检查线上上机表和线下上机表中FA样品的样本类型是否和样本类型库一致
        :return:
        '''
        for i in range(0, len(self.content_split_array)):
            library_type = self.content_split_array[i][self.index_library_type]
            sample_type = self.content_split_array[i][self.index_sample_type]
            sn = self.content_split_array[i][self.index_sn]

            if 'FA' in library_type:
                if sample_type in self.snDNAType.keys():
                    if len(sn.split('-')) > 1:
                        abbreviation = self.snDNAType[sample_type]
                        if abbreviation != self.snDNAType[sample_type]:
                            self.ctDNA_sample_type_error.append('\t'.join(self.content_split_array[i]))
                    else:

                        self.content_split_array[i][self.index_sn] = '-'.join([sn, self.snDNAType[sample_type]])

                else:
                    self.ctDNA_sample_type_error.append('\t'.join(self.content_split_array[i]))
        self.set_property('ctDNA_sample_type_error', self.ctDNA_sample_type_error)

    def checher_three(self):
        '''
        线下上机表的"生产/测试"列，必须要有生产或者测试信息
        线下上机表的测试样本，"生产/测试"列必须填写"测试"，同时样本名称前必须加"TEST"或者"test"
        :return:
        '''
        if u'分析类型' in self.header_split:                                    #线下上机表才有“分析类型”这个字段
            index_FXLX = self.get_index(u'分析类型')
            for i in range(0, len(self.content_split_array)):
                FXLX = self.content_split_array[i][index_FXLX]
                sn = self.content_split_array[i][self.index_sn]
                if FXLX in self.analysisType:
                    if FXLX == u'测试':
                        if sn.startswith('TEST') or sn.startswith('test'):
                            pass
                        else:
                            self.test_sample_without_test_prefix.append('\t'.join(self.content_split_array[i]))
                else:
                    self.FXLX_error.append('\t'.join(self.content_split_array[i]))
        self.set_property('test_sample_without_test_prefix', self.test_sample_without_test_prefix)
        self.set_property('FXLX_error', self.FXLX_error)

    def table_reform(self):
        '''
        统一内容模板，将转化后的上机表传到tool里
        时间复杂度有点儿高，先这样吧。
        :return:
        '''
        self.table_reformed.append(self.all_title)                                   # 先存入标题行
        for i in range(0, len(self.content_split_array)):
            line = list()
            for title in self.all_title:                                             # self.all_title是两个上机表标题的并集
                if title in self.header_split:                                       # self.header_split是读取的表中含有的标题
                    index = self.get_index(title)
                    line.append(self.content_split_array[i][index])
                else:
                    line.append('')
            self.table_reformed.append(line)

        self.set_property('table_reformed', self.table_reformed)

    def detect_encoding(self, path):
        detector = UniversalDetector()
        with open(path, 'r') as FIN:
            for line in FIN:
                detector.feed(line)
                if detector.done:
                    break
            detector.close()
        return detector.result['encoding']




