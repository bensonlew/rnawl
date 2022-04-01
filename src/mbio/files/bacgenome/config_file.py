# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171214
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class ConfigFileFile(File):
    """
    定义高通量数据拆分的输入config_file文件的格式
    只有一列
    """
    def __init__(self):
        super(ConfigFileFile, self).__init__()

    def check(self):
        super(ConfigFileFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("%s文件不存在", variables=(self.prop["path"]), code="41400401")
        self.get_info()
        #self.set_property('samples', self.samples)
        #self.set_property('fastqs', self.fastqs)

    def get_info(self):
        with open(self.prop["path"], "r") as f:
            for line in f:
                item = line.strip().split()
                if len(item) != 1:
                    raise FileError("config 文件不正确 ", code="41400402")


if __name__ == "__main__":
    a = ConfigFileFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ncRNA/nc_qc_after_list.txt")
    a.check()
    print a.prop["path"]
