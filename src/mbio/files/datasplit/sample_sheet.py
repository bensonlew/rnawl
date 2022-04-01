# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171218
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os


class SampleSheetFile(File):
    """
    定义数据拆分bcl2fastq的--sample-sheet文件
    第一行为[Data]
    # 第二行表头：单端"Lane\tSample_ID\tample_Plate\tSample_Well\tI7_Index_ID\tindex\tSample_Project\tDesription\n"
    # 第二行表头：双端"Lane\tSample_ID\tSample_Plate\tSample_Well\tI7_Index_ID\tindex\tI5_Index_ID\tindex\tSample_Project\tDesription\n"
    表格中表格中Lane、Sample_ID、Sample_Name、index、Sample_Project五列内容为必填
    Sample_ID内容需在整个表格中唯一(用于区分fastq.gz文件)
    """
    def __init__(self):
        super(SampleSheetFile, self).__init__()

    def check(self):
        super(SampleSheetFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("{}文件不存在".format())
        sample_ids, lane_index = [], {}
        with open(self.prop["path"], "r") as f:
            lines = f.readlines()
            # if lines[0].startswith("[Data]"):
            #    raise FileError("{}第一行要为[Data]".format(self.prop["path"]))
            # if lines[1] == "Lane\tSample_ID\tample_Plate\tSample_Well\tI7_Index_ID\tindex\tSample_Project\tDesription\n":
            #     pass
            # elif lines[1] == "Lane\tSample_ID\tSample_Plate\tSample_Well\tI7_Index_ID\tindex\tI5_Index_ID\tindex\tSample_Project\tDesription\n":
            #     pass
            # else:
            #     raise FileError("{}第二行表头错误，请检查".format(self.prop["path"]))
            if len(lines) < 3:
                raise FileError("没有拆分的文库，请检查")
            for line in lines[2:]:
                item = line.strip().split(",")
                if item[0] not in lane_index.keys():
                    lane_index[item[0]] = []
                # if item[6] in lane_index[item[0]]:
                #     raise FileError("lane:%s 里的index: %s重复，请检查" % (item[0], item[6]))
                # lane_index[item[0]].append(item[6])
                if len(item) > 8:
                    if item[6]+":"+item[8] in lane_index[item[0]]:
                        raise FileError("lane:%s 里的index: %s重复，请检查" % (item[0], item[6]+":"+item[8]))
                    lane_index[item[0]].append(item[6]+":"+item[8])
                else:
                    if item[6] in lane_index[item[0]]:
                        raise FileError("lane:%s 里的index: %s重复，请检查" % (item[0], item[6]))
                    lane_index[item[0]].append(item[6])
                if not item[0]:
                    raise FileError("Lane不能为空")
                if not item[1]:
                    raise FileError("Sample_ID不能空")
                if not item[6]:
                    raise FileError("index不能为空")
                if item[1] in sample_ids:
                    raise FileError("Sample_ID：{}重复，请检查".format(item[1]))
                else:
                    sample_ids.append(item[1])

if __name__ == "__main__":
    a = SampleSheetFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/tsanger/workspace/20200722/LibrarySplit_CF3-20200717mNova_20200722_112109/1.sample_sheet.csv")
    a.check()
