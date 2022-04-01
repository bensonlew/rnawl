# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.21

import os
import re
import xlrd
import xlwt
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class TraitFile(File):
    """
    检查性状文件，生成xls文件，用于性状分析（原因：性状分析的脚本读取的必须是xls文件）
    """
    def __init__(self):
        super(TraitFile, self).__init__()
        self.specimen_ids, self.trait_list = [], []
        self.excel = ""

    def check(self):
        if super(TraitFile, self).check():
            if not os.path.exists(self.path):
                raise FileError("文件:%s不存在，请检查", variables=(self.path), code="44801201")
            self.read_excel()
            if len(self.specimen_ids) != len(list(set(self.specimen_ids))):
                raise FileError("性状文件里样本不能重复，请检查", code="44801202")
            if len(self.trait_list) != len(list(set(self.trait_list))):
                raise FileError("性状文件里性状不能重复，请检查", code="44801203")
            self.set_property('is_excel', self.excel)
            self.set_property('specimen_ids', self.specimen_ids)
            print self.trait_type
            self.set_property('trait_type', self.trait_type)

    def read_excel(self):
        """
        读取excel文件
        """
        trait_value = []
        try:
            data = xlrd.open_workbook(self.path)
            table = data.sheets()[0]
            self.excel = True
            nrows = table.nrows  # 行数
            head = table.row_values(0)
            self.trait_list = head[1:]
            for i in xrange(1, nrows):
                row_values = table.row_values(i)
                trait_value = list(set(trait_value).union(set(row_values)))
                self.specimen_ids.append(row_values[0])
        except:
            self.excel = False
            with open(self.path, "r") as f:
                lines = f.readlines()
                if re.search(",", lines[0]):
                    head = lines[0].strip().split(",")
                    self.trait_list = head[1:]
                    for line in lines[1:]:
                        item = line.strip().split(",")
                        trait_value = list(set(trait_value).union(set(item[1:])))
                        self.specimen_ids.append(item[0])
                else:
                    head = lines[0].strip().split("\t")
                    self.trait_list = head[1:]
                    for line in lines[1:]:
                        item = line.strip().split("\t")
                        trait_value = list(set(trait_value).union(set(item[1:])))
                        self.specimen_ids.append(item[0])
        if len(trait_value) == 2 and set(trait_value).issubset(set(["0", "1"])):
            self.trait_type = "btl"
        else:
            self.trait_type = "qtl"

    def get_xls_trait(self, new_trait):
        """
        生成xls格式的trait文件
        """
        wbk = xlwt.Workbook()
        sheet = wbk.add_sheet('Sheet1', cell_overwrite_ok=True)
        with open(self.path, "r") as f:
            lines = f.readlines()
            for i in range(len(lines)):
                item = lines[i].strip().split("\t")
                for j in range(len(item)):
                    if item[j]:
                        value = item[j]
                    else:
                        value = "NaN"
                    sheet.write(i, j, value)
        wbk.save(new_trait)


if __name__ == "__main__":
    a = TraitFile()
    # a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/trit/trit3.xls")
    a.set_path("/mnt/ilustre/users/sanger-dev/workspace/20180801/Single_tsg_31236_QtlAnalysis_0801090753005038/trait_dir/GT.txt")
    a.check()
    # print a.prop["is_excel"]
    # a.get_xls_trait("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/trit/trit3_new.xls")
