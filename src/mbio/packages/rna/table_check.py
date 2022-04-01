# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# lastmodied: 20190311

import re
import os
import sys
import codecs
import xlrd

class Table(object):
    def __init__(self, path=None):
        self.file_path = path
        self.file_type = "txt"

    def set_path(self, path=None):
        '''
        设置ipath数据库路径，包括下载的svg文件和map图形KO对应关系
        '''
        self.file_path = path

    def check_xls(self):
        '''
        检查excel文件
        '''
        print "check xls"
        try:
            print self.file_path
            xlrd.open_workbook(self.file_path)
        except Exception as e:
            print e
            return False
        else:
            return True

    def check_unicode(self):
        try:
            f = codecs.open(self.file_path, encoding='utf-8', errors='strict')
            for line in f:
                pass
            print "Valid utf-8"
        except UnicodeDecodeError:
            print "invalid utf-8"
            return False
        else:
            return True

    def get_xls_sheet1(self, sheet_path):
        if self.file_type == "xls":
            xls_obj =  xlrd.open_workbook(self.file_path)
            table1 = xls_obj.sheets()[0]
            with open(sheet_path, 'w') as sheet_f:
                for row in range(table1.nrows):
                    sheet_f.write("\t".join(map(str, table1.row_values(row))) + "\n")
        else:
            raise Exception("文件不是excel格式")


    def check(self):
        if self.check_unicode():
            self.file_type = "txt"
            return True
        elif self.check_xls():
            self.file_type = "xls"
            return True
        else:
            return False
            # raise Exception("无法识别文件格式")


if __name__ == "__main__":
    t = Table()
    t.set_path(sys.argv[1])
    t.check()
    t.get_xls_sheet1(sys.argv[2])
