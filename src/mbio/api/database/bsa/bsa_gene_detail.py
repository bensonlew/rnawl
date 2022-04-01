# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.26

from api_base import ApiBase
import datetime
import os
import re


class BsaGeneDetail(ApiBase):
    """
    BSA基因详情页导表
    """
    def __init__(self, bind_object):
        super(BsaGeneDetail, self).__init__(bind_object)
        self._project_type = "bsa"

    def add_sg_gene_seq(self, main_id, seq_path):
        """
        基因详情页基因序列表
        seq_path: seq.fa
        """
        main_id = self.check_objectid(main_id)   # 检查id是否是OBjectID
        self.check_exists(seq_path)   # 检查文件是否存在
        seq = ''
        with open(seq_path, "r") as f:
            lines = f.readlines()
            try:
                seq = lines[0]
            except:
                pass
            length = len(seq)
        if seq:
            update_dict = {
                "seq": seq,
                "length": length,
                "status": "end"
            }
            self.update_db_record("sg_gene_seq", {"_id": main_id}, update_dict, upsert=True, multi=True)
        else:
            self.update_db_record("sg_gene_seq", {"_id": main_id}, {"status": "end"}, upsert=True, multi=True)

    def add_sg_gene_index_detail(self, index_id, index_path):
        """
        基因详情页基因型频率详情表
        index_path: filter.result.index
        """
        index_id = self.check_objectid(index_id)   # 检查id是否是OBjectID
        self.check_exists(index_path)   # 检查文件是否存在
        data_list = []
        title_list = ["chrom", "pos", "type", "ref"]  # 表头
        with open(index_path, "r") as f:
            lines = f.readlines()
            chrom = ""
            header = lines[0].strip().split("\t")
            keys = {}
            for i in range(len(header[4:-5])):
                s = header[4:-5][i]
                if re.search(r"(.+)-GT", s):
                    k = s.split("-GT")[0] + "_gt"
                elif re.search(r"(.+)-AD", s):
                    k = s.split("-AD")[0] + "_ad"
                else:
                    k = s
                keys[k] = i + 4
                title_list.append(k)
            for line in lines[1:]:
                if not line.startswith("@"):
                    item = line.strip().split("\t")
                    insert_data = {
                        "index_id": index_id,
                        "chrom": item[0],
                        "pos": int(item[1]),
                        "type": item[2],
                        "ref": item[3],
                        "annotation": item[-5],
                        "high": int(item[-4]),
                        "moderate": int(item[-3]),
                        "low": int(item[-2]),
                        "modifier": int(item[-1])
                    }
                    for k in keys.keys():
                        m = re.match(r"(.+)_(ad)", k)
                        n = re.match(r"(.+)_(gt)", k)
                        if m:
                            insert_data[k] = item[keys[k]]
                            s = m.group(1)
                            s_ad = 0
                            for i in item[keys[k]].split(","):
                                s_ad += int(i)
                            insert_data[s + "_dp"] = s_ad
                        elif n:
                            insert_data[k] = item[keys[k]]
                        else:
                            insert_data[k] = round(float(item[keys[k]]), 4)
                    data_list.append(insert_data)
        if data_list:
            title_list.append("annotation")
            self.col_insert_data("sg_gene_index_detail", data_list)
            self.update_db_record("sg_gene_index", {"_id": index_id}, {"index_title": title_list, "status": "end"})
        else:
            self.update_db_record("sg_gene_index", {"_id": index_id}, {"status": "end"})
        # self.col_insert_data("sg_gene_index_detail", data_list)

    def update_sg_gene_status(self, main_id, collection):
        """
        更新sg_gene_seq/sg_gene_index表的状态为failed
        """
        main_id = self.check_objectid(main_id)   # 检查id是否是OBjectID
        self.update_db_record(collection, {"_id": main_id}, {"status": "failed"})


if __name__ == "__main__":
    a = BsaGeneDetail(None)
    index_id = "5a9635e8a4e1af45055c7abb"
    index_path = "/mnt/ilustre/users/sanger-dev/workspace/20180302/Single_bsa_test_new_0302112710_377/GeneDetail/output/filter.result.index"
    # a.add_sg_gene_index_detail(index_id, index_path)
    a.update_sg_gene_status(main_id="5abb615fa4e1af62b9199001", collection="sg_gene_seq")
