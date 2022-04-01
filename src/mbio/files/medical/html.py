# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# creat at 20171114

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from HTMLParser import HTMLParser


class HtmlFile(File):
    '''
    html文件，获取表格内容
    '''
    def __init__(self):
        super(HtmlFile, self).__init__()
        self.parser = MyHTMLParser()
        self.tab_list = None

    def get_info(self):
        '''
        获取文件属性
        '''
        super(HtmlFile, self).get_info()
        with open(self.prop['path'], "r") as f:
            mycontent = f.read()
            self.parser.feed(mycontent)
        self.tab_list = self.parser.tab_list
        # self.set_property("tab_list", self.parser.tab_list)

    def check(self):
        '''
        检查文件格式
        '''
        if super(HtmlFile, self).check():
            self.get_info()
        else:
            raise FileError("文件格式错误")
        return True


class MyHTMLParser(HTMLParser):

    def __init__(self):
        HTMLParser.__init__(self)
        self.tab_list = list()

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.tab = []
        if tag == 'tr':
            self.row = []

    def handle_endtag(self, tag):
        if tag == 'table':
            self.tab_list.append(self.tab)
        if tag == 'tr':
            self.tab.append(self.row)

    def handle_data(self, data):
        if data.isspace() is False:
            self.row.append(data)


if __name__ == '__main__':
    path = ""
    a = HtmlFile()
    a.set_path(path)
    a.check()
    print a.tab_list
