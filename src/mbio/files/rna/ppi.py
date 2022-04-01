# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
# time: 2017.04.20


from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class PpiFile(File):
    """
    ppi文件夹格式
    """
    def __init__(self):
        super(PpiFile, self).__init__()

    def check(self):
        super(PpiFile, self).check()
        with open(self.prop["path"], "r") as r:
            # line = r.readlines()[0]
            line1 = r.readline()
            line2 = r.readline()
            if not line1 and not line2:
                return True
            else:
                if line1:
                    if str(line1.strip()) != "gene_id":
                        raise FileError("基因列表中第一行第一个字段必须为gene_id")
            return True

if __name__ == "__main__":
    a = PpiFile()
    a.set_path("/mnt/ilustre/users/sanger-test/workspace/20170704/Refrna_small_test_5/Express/output/diff/genes_diff/network_A_vs_B")
    a.check()

