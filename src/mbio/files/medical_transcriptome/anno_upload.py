# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2017.04.13
from biocluster.iofile import File
import re
from biocluster.core.exceptions import FileError
from collections import defaultdict


class AnnoUploadFile(File):
    """
    定义客户上传的kegg/go注释文件
    文件每行需为三列,以\t分隔
    第一列为基因或转录本id,第二列标注第一列id是基因或转录本，kegg注释文件第三列为该基因或转录本注释到KEGG数据库中的gene ID,以";"分隔
    go注释文件第三列为该基因或转录本注释到的GO term，以";"分隔
    """
    def __init__(self):
        super(AnnoUploadFile, self).__init__()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(AnnoUploadFile, self).check():
            super(AnnoUploadFile, self).get_info()
            with open(self.prop["path"], "r") as f:
                lines = f.readlines()
                for line in lines:
                    try:
                        line = line.strip().split("\t")
                        if len(line) == 3:
                            pass
                        else:
                            raise FileError("文件每行需为三列", code = "43700101")
                        if line[1] not in ["transcript", "gene"]:
                            raise FileError("文每行第二列要为gene或transcript，便于区分是基因注释还是转录本注释", code = "43700102")
                        items = line[2].split(";")
                        for item in items:
                            if re.match(r"K|GO.+$", item):
                                pass
                            else:
                                raise FileError("文件第三列应为注释到KEGG数据库中的gene ID", code = "43700103")
                    except:
                        raise FileError("文件每行不是以\t分隔", code = "43700104")

    def get_transcript_anno(self, outdir):
        with open(self.prop["path"], "rb") as f, open(outdir, "wb") as w:
            with open(self.prop["path"], "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    if line[1] == "transcript":
                        w.write(line[0] + "\t" + line[2] + "\n")

    def get_gene_anno(self, outdir):
        with open(self.prop["path"], "rb") as f, open(outdir, "wb") as w:
            with open(self.prop["path"], "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    if line[1] == "gene":
                        w.write(line[0] + "\t" + line[2] + "\n")
