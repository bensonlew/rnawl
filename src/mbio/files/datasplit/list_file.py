# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171214
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class ListFileFile(File):
    """
    定义高通量数据拆分的输入list_file文件的格式
    第一列为文件名称，可为绝对路径或者相对于list_file下的相对路径
    第二列为样本名称
    第三列为双端的l/r，若为单端则只有两列
    列间用空格或\t分隔
    """
    def __init__(self):
        super(ListFileFile, self).__init__()
        self.samples = {}
        self.fastqs = []

    def check(self):
        super(ListFileFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("%s文件不存在", variables=(self.prop["path"]), code="41800201")
        self.get_info()
        self.set_property('samples', self.samples)
        self.set_property('fastqs', self.fastqs)

    def get_info(self):
        base_dir = os.path.dirname(self.prop["path"])
        with open(self.prop["path"], "r") as f:
            for line in f:
                item = line.strip().split()
                fp = item[0]  # 通过接口投递时，检查不到文件，测试时先不检查
                self.fastqs.append(fp)
                if len(item) == 2:
                    if item[1] in self.samples.keys():
                        raise FileError("{}样本{}重复，请检查".format(self.prop["path"], item[1]))
                    else:
                        self.samples[item[1]] = []
                        self.samples[item[1]].append(fp)
                if len(item) >= 3:
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
                raise Exception('{}里样本{}重名，请改样本名或分开运行！'.format(self.prop["path"], s))


if __name__ == "__main__":
    a = ListFileFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ncRNA/nc_qc_after_list.txt")
    a.check()
    print a.prop["path"]
