# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.07

from api_base import ApiBase
import datetime
import os
import re


class Transgene(ApiBase):
    """
    转基因
    """
    def __init__(self, bind_object):
        super(Transgene, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def add_transgene(self, main_id, file, insert_file, ref_file, merge_sam):
        """
        插入位点统计表、插入位点详情表
        file:279_1.insert.xls
        insert_file:输入的外源片段序列。
        ref_file:参考基因组。
        merge_sam:merge.sam
        如果没结果时会给出报错。
        """
        self.check_exists(file)
        self.check_exists(insert_file)
        self.check_exists(ref_file)
        transgenosis_id = self.check_objectid(main_id)
        insert_data_list = []
        insert1_data_list = []
        merge_lines_list = []
        with open(insert_file, "r")as f:
            lines = f.readlines()
            header = re.findall(r">(.+?)\n", lines[0])[0]
            seq = lines[1]
        with open(ref_file, "r")as fref:
            ref_lines = fref.readlines()
        with open(merge_sam, "r")as frm:
            merge_lines = frm.readlines()
            for line in merge_lines:
                if line.startswith("@"):
                    pass
                else:
                    tmp = line.strip().split("\t")
                    if tmp[6] == "insert":
                        merge_lines_list.append(line)
        with open(file, "r")as fr:
            lines = fr.readlines()
            chr_dict = {}
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                if tmp[4] in chr_dict.keys():
                    chr_dict[tmp[4]] += 1
                else:
                    chr_dict[tmp[4]] = 1
                insert = seq[(int(tmp[1]) - 1):int(tmp[2])]
                ref_sign = ">" + tmp[4] + "\n"
                start_num = 0
                ref_seq = ""
                for ref_line in ref_lines:
                    if ref_line.startswith(">"):
                        if ref_line == ref_sign:
                            start_num = 1
                        else:
                            start_num = 0
                    else:
                        if start_num == 1:
                            ref_seq += ref_line.strip("\n")
                ref_start = tmp[5].strip().split("-")[0].strip()
                ref_end = tmp[5].strip().split("-")[1].strip()
                if int(ref_start) < 2001:
                    sequence_f = ref_seq[0:int(ref_start)]
                else:
                    sequence_f = ref_seq[(int(ref_start) - 2001):int(ref_start)]
                if int(ref_end) > (len(ref_seq) - 2001):
                    sequence_b = ref_seq[(int(ref_end) - 1):len(ref_seq)]
                else:
                    sequence_b = ref_seq[(int(ref_end) - 1):(int(ref_end) + 1999)]
                insert_data = {
                    "transgenosis_id": transgenosis_id,
                    "insertion_id": tmp[0],
                    "chr": tmp[4],
                    "pos": tmp[5],
                    "insert": insert,
                    "start": int(tmp[1]),
                    "end": int(tmp[2]),
                    "mapped_reads": tmp[6],
                    "sequence_f": sequence_f,
                    "sequence_b": sequence_b
                }
                detail_id = self.db["sg_transgenosis_detail"].insert_one(insert_data).inserted_id
                left_data_list = []
                right_data_list = []
                pos_list = tmp[5].strip().split(" ")
                left_pos = pos_list[0]
                right_pos = pos_list[2]
                for line in merge_lines_list:
                    temp = line.strip().split("\t")
                    if (int(left_pos) - 150) <= temp[3] <= int(left_pos):
                        left_data_list.append(len(temp[9]))
                    elif int(right_pos) <= temp[3] <= (int(right_pos) + 150):
                        right_data_list.append(len(temp[9]))
                insert_pic_data = {
                    "detail_id": detail_id,
                    "left_data": left_data_list,
                    "right_data": right_data_list,
                    "upstream_seq_num": left_pos,
                    "downstream_seq_num": right_pos,
                    "chr_id": tmp[4]
                }
                self.db["sg_transgenosis_pic"].insert_one(insert_pic_data)
                # insert_data_list.append(insert_data)
        # if insert_data_list:
        #     self.col_insert_data("sg_transgenosis_detail", insert_data_list)
        # else:
        #     self.bind_object.logger.info("转基因插入位点详情表结果为空！")
        for i in chr_dict.keys():
            insert_data1 = {
                "transgenosis_id": transgenosis_id,
                "chr": chr_dict[i],
                "chr_id": i,
                "header": header
            }
            insert1_data_list.append(insert_data1)
        if insert1_data_list:
            self.col_insert_data("sg_transgenosis_stat", insert1_data_list)
        else:
            self.bind_object.logger.info("转基因插入位点统计表结果为空！")
            # print "转基因插入位点统计表结果为空！"

if __name__ == "__main__":
    a = Transgene(None)
    main_id = "5c7ce9ac17b2bf0f241118b1"
    # file = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/ref_wgs_v2/transgenosis/output/279_1.insert.xls"
    file = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/ref_wgs_v2/transgenosis/output/279_1.insert.xls"
    insert_file = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/ref_wgs_v2/transgenosis/data/insert.fa"
    ref_file = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/ref_wgs_v2/transgenosis/data/ref.fa"
    merge_sam = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/ref_wgs_v2/transgenosis/sam/merge.sam"
    a.add_transgene(main_id, file, insert_file, ref_file, merge_sam)
