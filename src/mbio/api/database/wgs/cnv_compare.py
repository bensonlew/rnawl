# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 20180409

from api_base import ApiBase


class CnvCompare(ApiBase):
    def __init__(self, bind_object):
        """
        WGS项目导表
        """
        super(CnvCompare, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_cnv_diff_result(self, compare_id, file_path):
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        title = open(file_path, 'r').readlines()[0].strip().split("\t")
        with open(file_path, 'r') as r:
            lines = r.readlines()[1:]
            for line in lines:
                tmp = line.strip().split('\t')
                insert_data = {
                    "compare_id": compare_id,
                    "chr": tmp[0],
                    "pos1": tmp[1],
                    "pos2": tmp[2],
                    "len": tmp[3],
                    "type": tmp[4],
                    title[5]: tmp[5],
                    title[6]: tmp[6],
                    title[7]: tmp[7],
                    title[8]: tmp[8],
                    "gene_num": tmp[9],
                    "gene": tmp[10]
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_cnv_compare_detail", data_list)

if __name__ == "__main__":
    a = CnvCompare(None)
    a.add_cnv_diff_result("5acadac2a4e1af431d80b3d8", "/mnt/ilustre/users/sanger-dev/workspace/20180409/CnvDiff_wgs_wgs_test_0409144110_5548/output/cnv_diff/A8-10_Lands.xls")