# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200506

from api_base import ApiBase

# 注意：
# 这个脚本并没有被使用！！！！！！


class PrimerDesign(ApiBase):
    """
    引物设计
    """
    def __init__(self, bind_object):
        super(PrimerDesign, self).__init__(bind_object)

    def add_sg_primer_detail(self, primer_id, primer_result):
        """
        sg_primer_detail
        primer_result: variation.result
        """
        primer_id = self.check_objectid(primer_id)
        self.check_exists(primer_result)
        detail_title = ",".join(["sequence", "forward_primer", "forward_tm", "forward_gc", "forward_size", "forward_self_end_th",
                        "forward_self_any_th",
                        #"forward_hairpin_th", "forward_end_stability", "forward_penalty",
                        "reverse_primer", "reverse_tm", "reverse_gc", "reverse_size", "reverse_self_end_th",
                        "reverse_self_any_th",
                        #"reverse_hairpin_th", "reverse_end_stability", "reverse_penalty",
                        "produce_size", "start", "end", "pair_penalty", "pair_compl_any_th", "pair_compl_end_th"])
        self.update_db_record("sg_primer", {"_id": primer_id}, {"detail_title": detail_title})
        data_list = []
        with open(primer_result, "r") as f:
            lines = f.readlines()
            for line in lines[1:-1]:
                item = line.strip().split("\t")
                insert_data = {
                    "primer_id": primer_id,
                    "sequence": item[24],
                    "forward_primer": item[0],
                    "forward_tm": item[1],
                    "forward_gc": item[2],
                    "forward_size": item[3],
                    "forward_self_end_th": item[4],
                    "forward_self_any_th": item[5],
                    # "forward_hairpin_th": item[6],
                    # "forward_end_stability": item[7],
                    # "forward_penalty": item[8],
                    "reverse_primer": item[9],
                    "reverse_tm": item[10],
                    "reverse_gc": item[11],
                    "reverse_size": item[12],
                    "reverse_self_end_th": item[13],
                    "reverse_self_any_th": item[14],
                    # "reverse_hairpin_th": item[15],
                    # "reverse_end_stability": item[16],
                    # "reverse_penalty": item[17],
                    "produce_size": item[18],
                    "start": int(item[19]),
                    "end": int(item[20]) + 1 - int(item[12]),
                    "pair_penalty": item[21],
                    "pair_compl_any_th": item[22],
                    "pair_compl_end_th": item[23]}
                data_list.append(insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("{}文件为空！".format(primer_result))
            else:
                self.col_insert_data("sg_primer_detail", data_list)

if __name__ == "__main__":
    a = PrimerDesign(None)
    primer_id = "5c7ce9ac17b2bf0f24111999"
    primer_result = "/mnt/ilustre/users/sanger-dev/workspace/20200506/PrimerDesign_tsg_3421_0506162802184639_4732/PrimerDesign/output/variation.result"
    a.add_sg_primer_detail(primer_id, primer_result)