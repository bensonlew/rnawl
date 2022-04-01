# -*- coding: utf-8 -*-
# __author__ = 'xuting'

import numpy as np
import pandas as pd
import time
from memory_profiler import profile

def print_subtype_range(type_list):
    if type(type_list) == str:
        if type_list in set(["int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"]):
            print (np.iinfo(type_list))
        elif type_list in set(["float16", "float32", "float64"]):
            print (np.finfo(type_list))
        else:
            print ("%s不在检测范围" % type_list)
    elif type(type_list) == list:
        for i in type_list:
            print_subtype_range(i)
    else:
        raise Exception("type_list必须为字符或列表")

def mem_usage(pandas_obj):
    """
    统计一个pandas表一共占用的内存大小
    """
    if isinstance(pandas_obj, pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
    else:
        usage_b = pandas_obj.memory_usage(deep=True)
    usage_mb = usage_b / 1024.0 ** 2
    return "{:03.2f} MB".format(usage_mb)

class LargeData(object):
    def __init__(self, data):
        super(LargeData, self).__init__()
        if isinstance(data, pd.DataFrame):
            self.data = data
        else:
            raise Exception("数据类型不符合：%s" % type(data))
        self.str_col = list()
        self.int_col = list()
        self.float_col = list()

    def print_memory_usage(self, type_list):
        # print self.data.info()
        print self.data.get_dtype_counts()
        # print self.data.memory_usage(deep=True)
        memory_usage_all = self.data.memory_usage(deep=True).sum() / 1024.0 ** 2
        print "All memory usage for all columns: {:03.2f} MB".format(memory_usage_all)
        for dtype in type_list:
            if dtype not in set(["float", "int", "object"]):
                print "ERROR DEBUG: %s is not in type_list" % dtype
                continue
            selected_dtype = self.data.select_dtypes(include=[dtype])
            mean_usage_b = selected_dtype.memory_usage(deep=True).mean()
            mean_usage_mb= mean_usage_b / 1024.0 ** 2
            sum_usage_b = selected_dtype.memory_usage(deep=True).sum()
            sum_usage_b = sum_usage_b / 1024.0 ** 2
            print "Average memory usage for {} columns: {:03.2f} MB".format(dtype, mean_usage_mb)
            print "Sum memory usage for {} columns: {:03.2f} MB".format(dtype, sum_usage_b)

    def select_data(self, dtype, col_list):
        if len(col_list) == 0:
            data = self.data.select_dtypes(include=[dtype])
            head = data.columns
        else:
            data = self.data[col_list]
            head = col_list
        return data, head

    def convert_int(self, astype="unsigned", inplace=False):
        """
        astype: integer, signed, unsigned   coverted to possible min
                int8, int16, int32, int64 forced converted to this type, may lose value or raise Exception error
                unit8, unit16, unit32, unit64  also forced
        inplace: if True, change obj data
                 if False, return a new data
        """
        # int_data = self.data.select_dtypes(include=['int'])
        # int_head = int_data.columns
        int_data,int_head = self.select_data('int', self.int_col)
        if astype in set(["integer", "signed", "unsigned"]):
            converted_int = int_data.apply(pd.to_numeric, downcast=astype)
        else:
            converted_int = int_data.astype(astype)
        new_data = self.data.apply(lambda x: converted_int.loc[:, x.name] if x.name in int_head else x)
        if inplace:
            self.data = new_data
            return self
        else:
            return new_data

    # @profile(stream = open("alogwrite","a+"))
    def convert_float(self, astype="float", inplace=False):
        """
        astype: float converted to possible min
                float16, float32, float64 fored converted to this type, may lose value or raise Exception error
        inplace: if True, change obj data
                 if False, return a new data
        """
        # float_data = self.data.select_dtypes(include=['float'])
        # float_head = float_data.columns
        float_data, float_head = self.select_data('float', self.float_col)
        if astype == "float":
            converted_float = float_data.apply(pd.to_numeric, downcast="float")
        else:
            converted_float = float_data.astype(astype)
        new_data = self.data.apply(lambda x: converted_float.loc[:, x.name] if x.name in float_head else x)
        if inplace:
            self.data = new_data
            return self
        else:
            return new_data

    def convert_str(self, inplace=False):
        str_data,str_head = self.select_data('object', self.str_col)
        # if len(str_head) == 0:
        #     if inplace:
        #         return
        #     else:
        #         return self.data
        new_str_head = list()  # need to check 0.5
        # for col in str_head:
        #     num_uniq = float(len(self.data[col].unique()))
        #     num_total_uniq = len(self.data[col])
        #     if num_uniq / num_total_uniq < 0.5:
        #         # self.data.loc[:, col] = self.data[col].astype('category')
        #         new_str_head.append(col)
        num_total = float(len(str_data))
        # str_data.apply(lambda x: x.name)
        str_data.apply(lambda x: new_str_head.append(x.name) if len(x.unique()) / num_total < 0.5 else 0)
        if len(new_str_head) == 0 or new_str_head == [None]:
            if inplace:
                return self
            else:
                return self.data
        str_data.loc[:, new_str_head] = str_data[new_str_head].astype('category')
        new_data = self.data.apply(lambda: str_data[new_str_head] if x.name in new_str_head else x)
        if inplace:
            self.data = new_data
            return self
        else:
            return new_data

def table_parser(table_reader, chunk_size=100000, ignore_index=False):
    """
    对数据分片导入
    table_reader: TextFileReader
    chunk_size: 分片大小
    ignore_index: 是否取消索引列
    """
    loop = True
    chunks = list()
    while loop:
        try:
            chunk = table_reader.get_chunk(chunk_size)
            new_chunk = LargeData(chunk).convert_float(inplace=True).convert_int(inplace=True).data
            chunks.append(new_chunk)
        except StopIteration:
            loop = False
    df = pd.concat(chunks, ignore_index=ignore_index)
    return df

def test(table):
    start1 = time.time()
    data = pd.read_table(table, index_col=0)
    end1 = time.time()
    spend_time = end1 - start1
    print mem_usage(data)
    print "cost time(s): %s" % spend_time
    # data = pd.read_table(table)
    for i in range(50000, 500001, 50000):
        # chunk_size=100000
        print "Chunck Size = %s" % i
        start2 = time.time()
        reader = pd.read_table(table, index_col=0, iterator=True)
        data = table_parser(reader, chunk_size=i)
        end2 = time.time()
        spend_time = end2 - start2
        print mem_usage(data)
        print "cost time(s): %s" % spend_time
    """
    type_list = ["float", "int", "object"]
    obj = LargeData(data)
    obj.print_memory_usage(type_list)
    obj.convert_float(inplace=True)
    obj.print_memory_usage(type_list)
    obj.convert_int(inplace=True)
    obj.print_memory_usage(type_list)
    obj.convert_str(inplace=True)
    obj.print_memory_usage(type_list)
    """

if __name__ == '__main__':
    test("reads_number_relative.xls")