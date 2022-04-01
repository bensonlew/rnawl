# -*- coding: utf-8 -*-

import os
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class PopSummaryFile(File):
    """
    wentian.liu
    """

    def __init__(self):
        super(PopSummaryFile, self).__init__()

    def is_file(self, file_path):
        """
        检查是否是文件是否存在
        :param file_path:
        :return:
        """
        if not os.path.isfile(file_path) or not os.path.exists(file_path):
            raise FileError("原始文件中不存在{}文件！".format(file_path))

    def is_dir(self, dir_path):
        """
        检查文件夹是否存在
        :param dir_path:
        :return:
        """
        if not os.path.isdir(dir_path) or not os.path.exists(dir_path):
            raise FileError("原始文件中不存在{}路径！".format(dir_path))

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(PopSummaryFile, self).check():
            with open(self.path, "r") as r:
                data = r.readlines()
                for m in data:
                    temp = m.strip().split("\t")
                    # print len(temp)
                    if re.match(r'^#', temp[0]):
                        if len(temp) not in [21, 25]:
                            raise FileError("pop_summary的列数不对--应该是21列或者25列！", code="41500108")
                        else:
                            continue
                    else:
                        # print temp[-2]
                        if temp[19] == "--":
                            continue
                        else:
                            tem1 = temp[19].split(':')
                            if len(tem1) == 2:
                                if len(tem1[0].split(",")) == len(tem1[1].split(';')):
                                    pass
                                elif len(tem1[0].split(";")) == len(tem1[1].split(';')):
                                    pass
                                else:
                                    raise FileError("pop_summary文件EggNOG-ID列：%s命名不正确！".format(temp[19]), variables=(temp[19]), code="41500109")
                            # else:
                            #     raise FileError("pop_summary文件EggNOG-ID列：%s命名不正确！".format(temp[19]), variables=(temp[19]), code="41500110")

if __name__ == "__main__":
    a = PopSummaryFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/workspace/20190329/Single_annovar2/Annovar/output/anno_count/pop.summary")
    print a.check()
    print "检查通过"