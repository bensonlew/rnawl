# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.26

from biocluster.api.database.base import Base, report_check
from api_base import ApiBase
import datetime
import json
import re


class BsaRegion(ApiBase):
    """
    BSA关联区域定位及详情导表
    """
    def __init__(self, bind_object):
        super(BsaRegion, self).__init__(bind_object)
        self._project_type = "bsa"

    def add_sg_region(self, task_id, project_sn, slidingwin_id, wp, mp, wb, mb, name=None, params=None):
        """
        关联区域过滤主表
        """
        slidingwin_id = self.check_objectid(slidingwin_id)  # 检查id是否是OBjectID
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "slidingwin_id": slidingwin_id,
            "name": name if name else "origin_filter_analysis",
            "params": params,
            "desc": "关联区域过滤主表",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "wp": wp,
            "mp": mp,
            "wb": wb,
            "mb": mb
        }
        main_id = self.db['sg_region'].insert_one(insert_data).inserted_id
        self.db["sg_region"].update({"_id": main_id}, {"$set": {"main_id": main_id}}, upsert=True, multi=True)
        return main_id

    def add_sg_region_variant(self, region_id, gene_total_path, variant_total_path):
        """
        关联区域变异位点统计表
        gene_total_path: region.threshold.gene.total
        variant_total_path: region.threshold.variant.total
        """
        region_id = self.check_objectid(region_id)   # 检查id是否是OBjectID
        self.check_exists(gene_total_path)   # 检查文件是否存在
        self.check_exists(variant_total_path)   # 检查文件是否存在
        result = self.col_find_one("sg_region", {"main_id": region_id})
        task_id = result["task_id"]
        data_list = []
        start, end, snp_number, indel_number, eff_snp, eff_indel, gene_num, eff_gene = 0, 0, 0, 0, 0, 0, 0, 0
        with open(variant_total_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("@"):
                    item = line.strip().split("\t")
                    if start > int(float(item[1])):
                        start = int(float(item[1]))
                    if end < int(float(item[2])):
                        end = int(float(item[2]))
                    snp_number += int(item[3])
                    indel_number += int(item[5])
                    eff_snp += int(item[4])
                    eff_indel += int(item[6])
                    insert_data = {
                        "region_id": region_id,
                        "chr": item[0].split("@")[1],
                        "start": int(float(item[1])),
                        "end": int(float(item[2])),
                        "region_name": item[0].split("@")[1] + ":" + str(int(float(item[1]))) + "-" + str(int(float(item[2]))),
                        "snp_number": int(item[3]),
                        "indel_number": int(item[5]),
                        "eff_snp": int(item[4]),
                        "eff_indel": int(item[6])
                    }
                    data_list.append(insert_data)
        self.col_insert_data("sg_region_variant", data_list)
        with open(gene_total_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("@"):
                    item = line.strip().split("\t")
                    gene_num += int(item[3])
                    eff_gene += int(item[4])
                    query_dict = {
                        "region_id": region_id,
                        "chr": item[0].split("@")[1],
                        "region_name": item[0].split("@")[1] + ":" + str(int(float(item[1]))) + "-" + str(int(float(item[2])))
                    }
                    update_dict = {
                        "gene_num": int(item[3]),
                        "eff_gene": int(item[4])
                    }
                    self.update_db_record("sg_region_variant", query_dict, update_dict, "false",
                                          upsert=True, multi=True)
        # 这一条是总的total的值
        insert_data = {
            "region_id": region_id,
            "chr": "total",
            "start": "--",
            "end": "--",
            "region_name": "total:" + str(start) + "-" + str(end),
            "snp_number": snp_number,
            "indel_number": indel_number,
            "eff_snp": eff_snp,
            "eff_indel": eff_indel,
            "gene_num": gene_num,
            "eff_gene": eff_gene
        }
        self.db['sg_region_variant'].insert_one(insert_data)

    # @report_check
    def add_sg_region_vcf(self, region_id, vcf_total_path):
        """
        关联区域变异位点详情表--lasted modified by hongdong@20180820 修改导表bug
        vcf_total_path: region.threshold.vcf.total
        """
        region_id = self.check_objectid(region_id)   # 检查id是否是OBjectID
        self.check_exists(vcf_total_path)   # 检查文件是否存在
        data_list = []
        title_list = ["chrom", "pos", "type", "ref"]  # 表头
        variant_list = self.find_region_variant(region_id)
        with open(vcf_total_path, "r") as f:
            lines = f.readlines()
            # chrom = ""
            header = lines[0].strip().split("\t")
            keys = {}
            for i in range(len(header[4:-5])):
                s = header[4:-5][i]
                # k = s.split("-")[0] + "_" + s.split("-")[1].lower()
                if re.search(r"(.+)-GT", s):
                    k = s.split("-GT")[0] + "_gt"
                elif re.search(r"(.+)-AD", s):
                    k = s.split("-AD")[0] + "_ad"
                else:
                    k = s
                keys[k] = i + 4
                title_list.append(k)
            for line in lines[1:]:
                item = line.strip().split("\t")
                # if chrom != item[0]:
                #     chrom = item[0]
                #     query_dic = {
                #         "region_id": region_id,
                #         "chr": chrom
                #     }
                #     result = self.col_find_one("sg_region_variant", query_dic)
                #     variant_id = result["_id"]

                # if item[0] == chrom:
                variant_id = self.check_variant(variant_list, item[0], int(item[1]))
                insert_data = {
                    "variant_id": variant_id,
                    "chrom": item[0],
                    "pos": int(float(item[1])),
                    "type": item[2],
                    "ref": item[3],
                    "annotation": item[-5],
                    "high": int(item[-4]),
                    "moderate": int(item[-3]),
                    "low": int(item[-2]),
                    "modifier": int(item[-1])
                }
                for k in keys.keys():
                    m = re.match(r"(.+)_(gt)", k)
                    n = re.match(r"(.+)_(ad)", k)
                    if m:
                        insert_data[k] = item[keys[k]]
                    if n:
                        insert_data[k] = item[keys[k]]
                        s = n.group(1)
                        s_ad = 0
                        for i in item[keys[k]].split(","):
                            s_ad += int(i)
                        insert_data[s + "_dp"] = s_ad
                data_list.append(insert_data)
        # self.col_insert_data("sg_region_vcf", data_list)
        print "开始导入vcf文件"
        title_list.append("annotation")
        print title_list
        self.gevet_insert_data("sg_region_vcf", data_list, 1000000)
        self.update_db_record("sg_region", {"_id": region_id}, {"vcf_title": title_list})

    def find_region_variant(self, region_id):
        """
        根据region_id查找所有的变异区域
        :param region_id:
        :return: list
        """
        variant = []
        result = self.col_find("sg_region_variant", {"region_id": region_id})
        if result.count() != 0:
            for m in result:
                if m['chr'] != "total":
                    variant.append({"chr": m['chr'], "start": m['start'], "end": m['end'], "variant_id": m['_id']})
        return variant

    def check_variant(self, variant_list, chrs, pos):
        """
        输入pos位点检查 在那个区域
        :param variant_list:
        :param pos:
        :param chrs:
        :return:
        """
        variant_id = '123456789012345678901234'
        for m in variant_list:
            if chrs == m['chr']:
                if m['end'] > pos >= m['start']:
                    variant_id = m['variant_id']
                else:
                    continue
            else:
                continue
        return variant_id

    def add_sg_region_index(self, region_id, variant_total_path):
        """
        关联区域基因型频率详情表-- lasted modified by hongdong@20180820 修改导表bug
        variant_total_path: region.threshold.variant.total
        """
        region_id = self.check_objectid(region_id)   # 检查id是否是OBjectID
        self.check_exists(variant_total_path)   # 检查文件是否存在
        data_list = []
        title_list = ["chrom", "pos", "type", "ref"]  # 表头
        variant_list = self.find_region_variant(region_id)
        with open(variant_total_path, "r") as f:
            lines = f.readlines()
            # chrom = ""
            header = lines[1].strip().split("\t")
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
            for line in lines[2:]:
                if not line.startswith("@"):
                    # print line
                    item = line.strip().split("\t")
                    # if chrom != item[0]:
                    #     chrom = item[0]
                    #     query_dic = {
                    #         "region_id": region_id,
                    #         "chr": chrom
                    #     }
                    #     result = self.col_find_one("sg_region_variant", query_dic)
                    #     variant_id = result["_id"]
                    variant_id = self.check_variant(variant_list, item[0], int(item[1]))
                    insert_data = {
                        "variant_id": variant_id,
                        "chrom": item[0],
                        "pos": int(float(item[1])),
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
        self.col_insert_data("sg_region_index", data_list)
        title_list.append("annotation")
        print title_list
        self.update_db_record("sg_region", {"_id": region_id}, {"index_title": title_list})

    def update_sg_region_quantile(self, region_id, index_quantile_file):   # add by wzy,20180305
        """
        更新主表中的阈值及index值
        """
        region_id = self.check_objectid(region_id)  # 检查id是否是OBjectID
        self.check_exists(index_quantile_file)
        with open(index_quantile_file, 'r') as r:
            data = r.readlines()[1]
            temp = data.strip().split("\t")
            quantile = temp[0]
            index = temp[1][:6]
        self.update_db_record("sg_region", {"_id": region_id}, {"index": float(index), "quantile": float(quantile)})


if __name__ == "__main__":
    a = BsaRegion(None)
    # project_sn = 'bsa_test_new'
    # task_id = 'bsa_test_new'
    # slidingwin_id = "5a8fd9d5a4e1af0e64fca5a2"
    # region_id = a.add_sg_region(task_id, project_sn, slidingwin_id, name=None, params=None)
    # region_id = "5a97bc36a4e1af07498d49e1"
    # gene_total_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/region.threshold.gene.total"
    # variant_total_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/region.threshold.variant.total"
    # vcf_total_path = "/mnt/ilustre/users/sanger-dev/workspace/20180301/SlidingwinFilterAnalysis_bsa_test_new_0301163918_129/output/region_vcf/region.threshold.vcf.total"
    # a.add_sg_region_variant(region_id, gene_total_path, variant_total_path)
    a.add_sg_region_vcf("5b7a20a3ec02cc65f3b721a6",
                        "/mnt/ilustre/users/sanger/workspace/20180819/Bsa_i-sanger_107929/RegionAnalysis/output/reg"
                        "ion_vcf/region.threshold.vcf.total")
    a.add_sg_region_index("5b7a20a3ec02cc65f3b721a6",
                          "/mnt/ilustre/users/sanger/workspace/20180819/Bsa_i-sanger_107929"
                          "/RegionAnalysis/output/region_variant/region.threshold.variant.total")
    print "111"
    # a.add_sg_region_vcf(region_id, vcf_total_path)
    # a.add_sg_region_index(region_id, variant_total_path)
