# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171214
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class ListFileFile(File):
    """
    第一列分析名称，第二列批次样本名，第三列左端序列，第四列右端序列
    第二列批次样本名不能重复
    """
    def __init__(self):
        super(ListFileFile, self).__init__()
        self.samples = {}
        self.fastqs = []

    def check(self):
        super(ListFileFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("%s文件不存在", variables=(self.prop["path"]), code="45500105")
        # self.get_info()
        # self.set_property('samples', self.samples)
        # self.set_property('fastqs', self.fastqs)

    def get_info(self):
        base_dir = os.path.dirname(self.prop["path"])
        with open(self.prop["path"], "r") as f:
            for line in f:
                item = line.strip().split()
                fp = item[0]  # 通过接口投递时，检查不到文件，测试时先不检查
                self.fastqs.append(fp)
                if len(item) == 2:
                    if item[1] in self.samples.keys():
                        raise FileError("%s样本%s重复，请检查", variables=(self.prop["path"], item[1]), code="45500106")
                    else:
                        self.samples[item[1]] = []
                        self.samples[item[1]].append(fp)
                if len(item) == 3:
                    if item[1] in self.samples.keys():
                        if item[2] == "l":
                            self.samples[item[1]].insert(0, fp)
                        else:
                            self.samples[item[1]].append(fp)
                    else:
                        self.samples[item[1]] = []
                        self.samples[item[1]].append(fp)
        for s in self.samples.keys():
            if len(self.samples[s]) > 2:
                self.set_error('%s里样本%s重名，请改样本名或分开运行！', variables=(self.prop["path"], s), code="45500103")


if __name__ == "__main__":
    a = ListFileFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ncRNA/nc_qc_after_list.txt")
    a.check()
    print a.prop["path"]
