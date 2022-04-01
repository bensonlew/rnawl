# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from api_base import ApiBase
import datetime
import re
from collections import defaultdict


class Slidingwin(ApiBase):
    def __init__(self, bind_object):
        """
        BSA项目将标记筛选的结果（滑窗结果）导入数据库中
        __author__ = wangzhaoyue
        lasted modified by hd 20180612
        :param bind_object:
        """
        super(Slidingwin, self).__init__(bind_object)
        self._project_type = "bsa"    # 连接到测试的mongo集群 192.168.10.16

    def get_task_info(self, task_id):
        """
        获取任务id相关信息
        :param task_id:
        :return:
        """
        return self.col_find_one("sg_task", {"task_id": task_id})

    def add_sg_slidingwin(self, task_id, project_sn, member_id,  index_path, variant_path, results_path, slid_path, mb,
                          wp=None, mp=None, wb=None, params=None, name=None):
        """
        添加 标记筛选和分析主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param wp: 野生型亲本名称
        :param mp: 突变型亲本名称
        :param wb: 野生型混池名称
        :param mb: 突变型混池名称，必填
        :param params:
        :param name:
        :return:
        """
        if not mb:
            self.set_error("标记筛选和分析主表中缺少突变型混池名称信息", code="51500701")
            raise Exception("标记筛选和分析主表中缺少突变型混池名称信息!")
        name = name if name else "origin_slidingwin"
        params = params if params else ""
        main_id = self.add_main_table("sg_slidingwin", task_id, project_sn, params, name, "标记筛选和分析主表", member_id)
        self.update_db_record("sg_slidingwin", {"_id": main_id},
                              {"wp": wp, "mp": mp, "wb": wb, "mb": mb, "calc_index_path": index_path,
                               "calc_variant_path": variant_path, "slidingwin_result_path": results_path,
                               "slid_result_path": slid_path})
        return main_id

    def update_sg_slidingwin_path(self, slidingwin_id, index_path, variant_path, results_path, slid_path):
        """
        更新 标记筛选和分析主表结果路径
        """
        slidingwin_id = self.check_objectid(slidingwin_id)   # 检查id是否是OBjectID
        # self.check_exists(index_path)   # 检查文件是否存在
        # self.check_exists(variant_path)   # 检查文件是否存在
        # self.check_exists(results_path)   # 检查文件是否存在
        # self.check_exists(slid_path)  # 检查文件是否存在
        self.update_db_record("sg_slidingwin", {"_id": slidingwin_id},
                              {"calc_index_path": index_path, "calc_variant_path": variant_path,
                               "slidingwin_result_path": results_path, "slid_result_path": slid_path})
        self.bind_object.logger.info("更新sg_slidingwin路径成功")

    def add_sg_slidingwin_stat(self, stat_file_path, slidingwin_id):
        """
        样本信息比对细节表
        :param slidingwin_id: 主表id
        :param stat_file_path: 统计信息文件
        :return:
        """
        slidingwin_id = self.check_objectid(slidingwin_id)   # 检查id是否是OBjectID
        self.check_exists(stat_file_path)   # 检查文件是否存在
        data_list = []
        with open(stat_file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                insert_data = {
                    "slidingwin_id": slidingwin_id,
                    "chr": line[0],
                    # "gene_number": int(line[1]),
                    "snp_number": int(line[1]),
                    "indel_number": int(line[2]),
                    "eff_snp": int(line[3]),
                    "eff_indel": int(line[4])
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_slidingwin_stat", data_list)

    def add_sg_manhattan(self, slidingwin_id, results_path, ref_chrlist=None):
        """
        曼哈顿图主表
        :param slidingwin_id: 主表id
        :param results_path: sliding-win.result
        :param ref_chrlist: ref_chrlist   # 不需要
        :return:
        """
        slidingwin_id = self.check_objectid(slidingwin_id)  # 检查id是否是OBjectID
        self.check_exists(results_path)  # 检查文件是否存在
        main_data = [{
            "origin_id": slidingwin_id,
            "name": "曼哈顿图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": "region_mahattan"
        }]
        main_id = self.col_insert_data("sg_manhattan", main_data)

        chr_list = []
        sca_list = []
        other_list = []
        final_list = []
        insert_data = []
        with open(results_path)as fr:
            lines = fr.readlines()
            col = lines[1].strip().split("\t")
            chr_data = defaultdict(list)  # {染色体：[pos],[index]]}或{染色体：[pos],[index],[index2],[delta_index]}
            if len(col) == 9:   # 只有一个index
                for line in lines[1:]:
                    tmp = line.strip().split("\t")
                    chr_name = tmp[0].strip("\"")
                    if chr_name not in chr_data.keys():
                        chr_data[chr_name] = [[float(tmp[2])], [float(tmp[3][:6])]]  # 有科学计数法，需要取浮点型在变成字符串取值
                    else:
                        chr_data[chr_name][0].append(float(tmp[2]))
                        chr_data[chr_name][1].append(float(tmp[3][:6]))
            else:  # index.index2,delta_index
                for line in lines[1:]:
                    tmp = line.strip().split("\t")
                    chr_name = tmp[0].strip("\"")
                    if chr_name not in chr_data.keys():
                        chr_data[chr_name] = [[float(tmp[2])], [float(tmp[3][:6])], [float(tmp[4][:6])], [float(tmp[5][:6])]]
                    else:
                        chr_data[chr_name][0].append(float(tmp[2]))
                        chr_data[chr_name][1].append(float(tmp[3][:6]))
                        chr_data[chr_name][2].append(float(tmp[4][:6]))
                        chr_data[chr_name][3].append(float(tmp[5][:6]))
        for i in chr_data:   # change by wzy 20180403 对染色体排序
            if i.startswith("chr"):
                num = i.strip().split("chr")[-1]
                try:
                    chr_list.append(int(num))
                except:
                    chr_list.append(num)
            elif i.startswith("sca"):
                num = i.strip().split("sca")[-1]
                try:
                    sca_list.append(int(num))
                except:
                    sca_list.append(num)
            else:
                other_list.append(i)
        chr_list.sort()  # add by wzy 20180403 对染色体排序
        sca_list.sort()
        for i in chr_list:
            final_list.append("chr" + str(i))
        for i in sca_list:
            final_list.append("sca" + str(i))
        for i in other_list:
            final_list.append(i)
        for key in final_list:
            data = {
                "manhattan_id":  main_id,
                "name": str(key),
                "pos_len": len(chr_data[key][0]),
                "x_categories": chr_data[key][0],
                "index_value": chr_data[key][1],
            }
            if len(chr_data[key]) == 4:  # 当有index2，delta_index时，需要加进去
                data.update({"index2_value": chr_data[key][2], "delta_value": chr_data[key][3]})
            insert_data.append(data)
        self.col_insert_data("sg_manhattan_detail", insert_data)
        self.update_db_record("sg_manhattan", {"_id": main_id}, {"chr_list": final_list})

    def add_sg_circos(self, slidingwin_id, gene_num_path, chr_num_path_, chr_num_path, snp_path, indel_path, delta_path):
        """
        导入circos的所有的数据入库
        :param slidingwin_id: 滑窗主表id
        :param gene_num_path: gene.num.csv
        :param chr_num_path_: circos.chrlist 远程磁盘路径
        :param snp_path: snp.win.csv
        :param indel_path: indel.win.csv
        :param delta_path: sliding.win.csv
        :param chr_num_path  真实路径
        :return:
        """
        slidingwin_id = self.check_objectid(slidingwin_id)  # 检查id是否是OBjectID
        # self.check_exists(gene_num_path)  # 检查文件是否存在
        self.check_exists(chr_num_path)  # 检查文件是否存在
        # self.check_exists(snp_path)  # 检查文件是否存在
        # self.check_exists(indel_path)  # 检查文件是否存在
        # self.check_exists(delta_path)  # 检查文件是否存在
        chr_list = []
        with open(chr_num_path, "r") as f:
            for line in f:
                m = re.match(r'.*\{\"id\"\:\"(.+)\",\"label\".+}', line)
                if m:
                    chr_list.append(m.group(1))
        # print chr_list
        data = [{
            "origin_id": slidingwin_id,
            "name": "circos图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "gene_num_path": gene_num_path,
            "chr_num_path": chr_num_path_,
            "snp_path": snp_path,
            "indel_path": indel_path,
            "delta_path": delta_path,
            "chr_list": chr_list,
            "location": "slidingwin_circos"
        }]
        self.col_insert_data("sg_circos", data)

if __name__ == "__main__":
    # project_sn = 'bsa_test'
    task_id = 'bsa_test'
    # member_id = 'm_test'
    # index_path = '/mnt/ilustre/users/sanger-dev/workspace/20180223/Single_test_module_slidingwin_analysis_wzy3/SlidingwinAnalysis/output/index-calc.result.index'
    # variant_path = '/mnt/ilustre/users/sanger-dev/workspace/20180223/Single_test_module_slidingwin_analysis_wzy3/SlidingwinAnalysis/output/index-calc.result.variant'
    # results_path = '/mnt/ilustre/users/sanger-dev/workspace/20180223/Single_test_module_slidingwin_analysis_wzy3/SlidingwinAnalysis/output/sliding-win.result'
    t = Slidingwin(None)
    # main_id = t.add_sg_slidingwin(task_id, project_sn, member_id, index_path, variant_path, results_path, "ZJU_co", "B23XC")
    # stat_file = '/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/bsa/final.stat'
    # t.add_sg_slidingwin_stat(stat_file, main_id)
    slidingwin_id = "5a8fd9d5a4e1af0e64fca5a2"
    # results_path = '/mnt/ilustre/users/sanger-dev/workspace/20180301/Single_test_module_slidingwin_analysis_ninanjie_wzy/SlidingwinAnalysis/output/sliding-win.result'
    # manhattan_id = "5a97c4b6a4e1af043ea1189f"
    t.add_sg_manhattan(slidingwin_id, "/mnt/ilustre/users/sanger-dev/workspace/20180316/Bsa_tsg_28536/SlidingwinAnalysis/output/sliding-win.result")
    # region_id = "5a93a887a4e1af1d48f392c5"
    # t.add_sg_manhattan_value(manhattan_id, region_id, 0.999, 0.454545454)
    # gene_num_path = '/mnt/ilustre/tsanger-data/rerewrweset/BSA_files/circos/tsg_28714/gene.num.csv'
    # chr_num_path = '/mnt/ilustre/tsanger-data/rerewrweset/BSA_files/circos/tsg_28714/circos.chrlist'
    # snp_path = '/mnt/ilustre/tsanger-data/rerewrweset/BSA_files/circos/tsg_28714/snp.win.csv'
    # indel_path = '/mnt/ilustre/tsanger-data/rerewrweset/BSA_files/circos/tsg_28714/indel.win.csv'
    # delta_path = '/mnt/ilustre/tsanger-data/rerewrweset/BSA_files/circos/tsg_28714/sliding.win.csv'
    # t.add_sg_circos(slidingwin_id, gene_num_path, chr_num_path, snp_path, indel_path, delta_path)
