# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.05

from api_base import ApiBase
import datetime
import os


class SvCall(ApiBase):
    """
    sv统计导表
    """
    def __init__(self, bind_object):
        super(SvCall, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def add_sg_sv_call(self, project_sn, task_id, params=None, name=None):
        """
        sg_indel_call
        添加主表，add by hongdong
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_sv_call",
            "params": params if params else "null",
            "desc": "SV检测主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_sv_call", data_list)
        self.update_db_record("sg_sv_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_sv_call_stat(self, main_id, file):
        """
        sv数据统计表
        file:stat.txt
        """
        self.check_exists(file)
        call_id = self.check_objectid(main_id)
        insert_data_list = []
        with open(file, "r")as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                insert_data = {
                    "call_id": call_id,
                    "specimen_id": tmp[0],
                    "del": tmp[1],
                    "ins": tmp[2],
                    "dup": tmp[3],
                    "inv": tmp[4],
                    "bnd": tmp[5]
                }
                insert_data_list.append(insert_data)
        self.col_insert_data("sg_sv_call_stat", insert_data_list)

    def add_sg_sv_len(self, main_id, task_id, file_path):
        """
        sv长度分布图
        file_path:sv_stat_v2.output_dir路径
        导入output中的XXX.sv.length.txt文件
        """
        origin_id = self.check_objectid(main_id)
        categories = ["0-1k", "1-2k", "2-3k", "3-4k", "4-5k", "5-6k", "6-7k", "7-8k", "8-9k", "9-10k", ">10k"]
        file_list = os.listdir(file_path)
        for file in file_list:
            if file.endswith(".sv.length.txt"):
                temp = file.strip().split(".")
                sample = temp[0]
                with open(os.path.join(file_path, file), "r")as fr:
                    lines = fr.readlines()
                    for line in lines[1:]:
                        tmp = line.strip().split("\t")
                        location = ""
                        if tmp[0] == "DEL":
                            location = "sv_del_length"
                        elif tmp[0] == "INS":
                            location = "sv_ins_length"
                        elif tmp[0] == "DUP":
                            location = "sv_dup_length"
                        elif tmp[0] == "INV":
                            location = "sv_inv_length"
                        elif tmp[0] == "BND":
                            location = "sv_bnd_length"
                        # list_value = tmp[1:]
                        list_value = map(int, tmp[1:])
                        bar_id = self.sg_bar(task_id, origin_id, sample, categories, 1, location, "",
                                             {"type": "sample", "title": sample})
                        self.sg_bar_detail(bar_id, sample, list_value)

if __name__ == "__main__":
    a = SvCall(None)
    main_id = "5cb4481d17b2bf62a0a4cb99"
    file = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/output/sv_stat_v2/stat.txt"
    # a.add_sg_sv_call_stat(main_id, file)
    # file_len_1 = "/mnt/ilustre/users/sanger-dev/workspace/20190307/Single_test_sv_stat_v2_2019030709091416/SvStatV2/output/SRR5739119.sv.length.txt"
    task_id = "sanger_85433"
    # a.add_sg_sv_len(main_id, file_len_1, task_id, "SRR5739119")
    file_len_2 = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/output/sv_stat_v2"
    a.add_sg_sv_call_stat(main_id, file)
    # a.add_sg_sv_len(main_id, task_id, file_len_2)