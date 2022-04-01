"""
format and check table used as input file
with random forest tool
version 0.0.0 
author: yingnn
date: 2017.10.16

"""
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from biocluster.config import Config


class TableRfFile(File):
    """
    format and check table used as input file
    with random forest tool
    
    """
    def __init__(self):
        super(TableRfFile, self).__init__()

    def check(self):
        if super(TableRfFile, self).check():
            self.get_info()
            return True
            
    def get_info(self):
        if 'path' in self.prop.keys():
            self.new_path = self.prop['path']
            self.set_property("new_table", self.new_path)
        else:
            raise FileError("something wrong with file path")


if __name__ == "__main__":
    import argparse
    ap= argparse.ArgumentParser(description='rf test')
    ap.add_argument('file', help='input file')
    ap.add_argument('-p', '--prefix', help='out prefix', default='1')
    args= ap.parse_args()
    a = TableRfFile()
    a.set_path(args.file)
    a.check()

