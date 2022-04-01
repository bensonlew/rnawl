# -*_ coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modified: 20180329
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class ListFileFile(File):
    """
    定义高通量数据拆分的输入list_file文件的格式
    第一列为文件名称，可为绝对路径或者相对于list_file下的相对路径
    第二列为样本名称
    第三列为双端的l/r/s
    列间用空格或\t分隔
    """
    def __init__(self):
        super(ListFileFile, self).__init__()
        self.samples = {}

    def check(self):
        super(ListFileFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("%s文件不存在", variables=(self.prop["path"]), code="41400601")

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(ListFileFile, self).get_info()
        self.get_list_info()
        self.set_property('samples', self.samples)

    def get_list_info(self):
        base_dir = os.path.dirname(self.prop["path"])
        with open(self.prop["path"], "r") as f:
            for line in f:
                item = line.strip().split()
                fp = item[0]  # 通过接口投递时，检查不到文件，测试时先不检查
                if len(item) == 3:
                    if item[1] not in self.samples.keys():
                        self.samples[item[1]] = ['','','']
                    if item[2] == "l":
                        #self.samples[item[1]].insert(0, fp)
                        self.samples[item[1]][0] = fp
                    elif item[2] == 'r':
                        self.samples[item[1]][1] = fp  #zouguanqing 20181102
                    else:
                        self.samples[item[1]][2] = fp



if __name__ == "__main__":
    a = ListFileFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ncRNA/nc_qc_after_list.txt")
    a.check()
    print a.prop["path"]
