# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.07

from api_base import ApiBase
import datetime
import os


class PrimerDesign(ApiBase):
    """
    引物设计
    """
    def __init__(self, bind_object):
        super(PrimerDesign, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def add_sg_primer_detail(self, primer_id, primer_result, primer_num, download_file=None):
        """
        sg_primer_detail
        primer_result: variation.result
        """
        primer_id = self.check_objectid(primer_id)
        self.check_exists(primer_result)
        if download_file:
            self.update_db_record("sg_primer", {"_id": primer_id}, {"download_path": download_file})
        detail_title = ["primer_id", "marker_id", "chr", "pos", "type"]
        for num in range(primer_num):
            detail_title.append(("forward_primer" + str(num + 1)))
            detail_title.append(("forward_tm" + str(num + 1)))
            detail_title.append(("forward_gc" + str(num + 1)))
            detail_title.append(("forward_len" + str(num + 1)))
            detail_title.append(("reverse_primer" + str(num + 1)))
            detail_title.append(("reverse_tm" + str(num + 1)))
            detail_title.append(("reverse_gc" + str(num + 1)))
            detail_title.append(("reverse_len" + str(num + 1)))
            detail_title.append(("product_size" + str(num + 1)))
            detail_title.append(("start" + str(num + 1)))
            detail_title.append(("end" + str(num + 1)))
        self.update_db_record("sg_primer", {"_id": primer_id}, {"detail_title": detail_title})
        data_list = []
        with open(primer_result, "r") as f:
            # lines = f.readlines()
            # for line in lines[1:]:
            f.next()
            for line in f:
                item = line.strip().split("\t")
                if len(item) < 9 + 11 * (primer_num):
                    continue
                line_length = len(item)
                insert_data = {
                    "primer_id": primer_id,
                    "marker_id": str(item[0]) + "_" + str(item[1]),
                    "chr": item[0],
                    "pos": int(item[1]),
                    "type": item[3],
                    # "marker_start": item[7],
                    # "marker_end": item[8],
                }
                if line_length > 9:
                    for num in range(primer_num):
                        insert_data[("forward_primer" + str(num + 1))] = item[(9 + 11 * num)]
                        insert_data[("forward_tm" + str(num + 1))] = item[(10 + 11 * num)]
                        insert_data[("forward_gc" + str(num + 1))] = item[(11 + 11 * num)]
                        insert_data[("forward_len" + str(num + 1))] = item[(12 + 11 * num)]
                        insert_data[("reverse_primer" + str(num + 1))] = item[(13 + 11 * num)]
                        insert_data[("reverse_tm" + str(num + 1))] = item[(14 + 11 * num)]
                        insert_data[("reverse_gc" + str(num + 1))] = item[(15 + 11 * num)]
                        insert_data[("reverse_len" + str(num + 1))] = item[(16 + 11 * num)]
                        insert_data[("product_size" + str(num + 1))] = item[(17 + 11 * num)]
                        insert_data[("start" + str(num + 1))] = item[(18 + 11 * num)]
                        insert_data[("end" + str(num + 1))] = item[(19 + 11 * num)]
                    data_list.append(insert_data)
        self.col_insert_data("sg_primer_detail", data_list)

if __name__ == "__main__":
    a = PrimerDesign(None)
    primer_id = "5c7ce9ac17b2bf0f24111999"
    primer_result = "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20200219/PrimerDesign_majorbio_238342_0219143052649791_1362/output/primer_design/variation.result"
    download_file = "/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20191021/PrimerDesign_i-sanger_210305_1021113322352019_1857/output/primer_design/new_variant_result.xls"
    primer_num = 3
    a.add_sg_primer_detail(primer_id, primer_result, primer_num)
