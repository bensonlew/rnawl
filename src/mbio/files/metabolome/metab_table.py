# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import defaultdict
import pandas as pd


class MetabTableFile(File):
    """
    定义上传的metab_table格式文件
    检查表头名称，个数，顺序
    """
    def __init__(self):
        super(MetabTableFile, self).__init__()
        self.data = None

    def get_info(self):
        """
        获取文件属性
        """
        super(MetabTableFile, self).get_info()
        info = self.get_file_info()
        # self.set_property("sample_number", len(info[0]))
        # self.set_property("sample", info[0])
        # self.set_property("group_scheme", info[1])
        # self.set_property("is_empty", info[2])

    def get_file_info(self):
        """
        获取metab_table文件的信息
        """
        self.data = pd.read_csv(self.prop['path'], sep='\t')
        try:
            self.data['']
        except:
            pass

    def format_check(self):
        ## 改成检查固定列
        head_list = ["Metabolite", "Mode", "KEGG Compound ID"]
        for li in head_list:
            if li not in self.data.columns.tolist():
                raise Exception('缺失固定的列 %s' %li)

        # head_list = ["Metabolite", "Mode", "Formula", "m/z", "RT (min)", "KEGG Compound ID", "HMDB_ID", "CAS number"]
        # if self.data.columns[:8].tolist() == head_list:
        #     pass
        # else:
        #     raise Exception("上传数据表表头名称不对应")
        #     raise FileError("上传数据表表头名称不对应", code="44700601")
        if True in self.data.isnull().any().tolist():
            raise Exception("表格中有缺失值，需检查")
            raise FileError("上传数据表中有缺失值，需检查，并用'-'填补", code="44700602")

    def check(self):
        if super(MetabTableFile, self).check():
            self.get_info()
            self.format_check()

