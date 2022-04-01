# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171214
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os
import json


class LibraryParamsFile(File):
    """
    定义高通量数据拆分
    """
    def __init__(self):
        super(LibraryParamsFile, self).__init__()

    def check(self):
        super(LibraryParamsFile, self).check()
        if not os.path.exists(self.prop["path"]):
            raise FileError("{}文件不存在".format(self.prop["path"]))
        self.dump_json()

    def dump_json(self):
        """
        解析json文件
        """
        f = open(self.prop["path"], "rb")
        try:
            json_dict = json.loads(f.read())
        except:
            raise FileError("json: %s格式不正确" % self.prop["path"])
        # if "dna" in json_dict.keys():
        #     for dna_info in json_dict["dna"]:
        #         for lib in dna_info["samples"]:
        #             sample_list = dna_info["samples"][lib]
        #             new_sample_list = list(set(sample_list))
        #             if len(sample_list) != len(new_sample_list):
        #                 raise FileError("dna里文库：%s对应的样本有重复的，请检查" % lib)
        # if "meta" in json_dict.keys():
        #     meta_info = json_dict["meta"]
        #     for i in range(len(meta_info)):
        #         lib_info = meta_info[i]["lib_info"]
        #         samples = []
        #         lib_len = []
        #         lib_name = {}
        #         for l in lib_info:
        #             lib = l[l.keys()[0]]
        #             if lib["lib_insert_size"] in lib_len:
        #                 for s in lib["samples"]:
        #                     if s not in samples:
        #                         lib_name[lib] = lib["lib_path"]
        # return json_dict

    # def get_info(self):
    #     """
    #     获取文件属性
    #     """
    #     super(LibraryParamsFile, self).get_info()
    #     self.json_dict = self.dump_json()
    #     types = self.json_dict.keys()
    #     for ty in types:
    #         if ty not in ["microbial_genome", "rna", "ncrna", "mirna", "meta", "metagenomic", "dna"]:
    #             raise FileError("项目类型{}不在范围内".format(ty))
    #     self.params_info = {}


if __name__ == "__main__":
    a = LibraryParamsFile()
    b = {'library_split': [{'lane1': {'sample_sheet': '/mnt/ilustre/users/isanger/sanger-dev/workspace/20181213/Bcl2fastq_5c1202eaa209543324304426_20181213_151116/sample_sheet.csv', 'data_path': '/mnt/clustre/upload/hiseq/hiseq4000/2018/BCL/20180917sXten', 'bases_mask': 'y151,i6nn,y151', 'barcode_mismatch': 0}}]}
    with open("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/isanger/datasplit_v2/module_test/library_split.json", "w") as w:
        w.write(json.dumps(b) + "\t")
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/isanger/datasplit_v2/module_test/library_split.json")
    a.check()
    # a.dump_json()
