# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 20180412

from api_base import ApiBase
import datetime
from collections import defaultdict
import gc


class SnpIndelCompare(ApiBase):
    def __init__(self, bind_object):
        """
        WGS snp indel 比较分析接口的导表
        """
        super(SnpIndelCompare, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_snp_indel_compare_stat(self, task_id, compare_id, file_path, types="snp"):
        """
        页面snp/indel差异统计表  indel/SNP差异功能统计直方图
        :param compare_id:
        :param task_id:
        :param file_path:  Ann.stat
        :param types:
        :return:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        x_categories = []
        value = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines:
                tmp = line.strip().split('\t')
                x_categories.append(tmp[0])
                value.append(int(tmp[1]))
                insert_data = {
                    "compare_id": compare_id,
                    "type": tmp[0],
                    "num": int(tmp[1])
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            if types == "snp":
                self.col_insert_data("sg_snp_compare_stat", data_list)
            else:
                self.col_insert_data("sg_indel_compare_stat", data_list)
        main_id = self.sg_bar(task_id, compare_id, "", x_categories, 1, location="{}_function_stat".format(types))
        self.sg_bar_detail(main_id, "{}_function_stat".format(types), value, is_show_log="true", types='1')

    def add_sg_snp_indel_compare_eff_stat(self, task_id, compare_id, file_path, types="snp"):
        """
        snp/indel差异功效统计表 和 snp/Indel差异功效统计直方图
        :param compare_id:
        :param file_path:  Eff.stat
        :param types:
        :param task_id:
        :return:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        x_categories = []
        value = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines:
                tmp = line.strip().split('\t')
                x_categories.append(tmp[0])
                value.append(int(tmp[1]))
                insert_data = {
                    "compare_id": compare_id,
                    "type": tmp[0],
                    "num": int(tmp[1])
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            if types == "indel":
                self.col_insert_data("sg_indel_compare_eff_stat", data_list)
            else:
                self.col_insert_data("sg_snp_compare_eff_stat", data_list)
        main_id = self.sg_bar(task_id, compare_id, "", x_categories, 1, location="{}_eff_stat".format(types))
        self.sg_bar_detail(main_id, "{}_eff_stat".format(types), value, is_show_log="true", types='1')

    def add_sg_snp_indel_compare_detail(self, compare_id, file_path, types="snp", analysis_type="sample", download_path=None):
        """
        差异snp/indel详情表--样本样本比较,组与组比较，组件比较
        :param compare_id:
        :param file_path:  pop.variant
        :param types:
        :param analysis_type:sample/one_group/two_group
        :return:
        """
        if analysis_type not in ['sample', 'one_group', "two_group"]:
            raise Exception("{}不合法！, 必须为sample/one_group/two_group其中一种！".format(analysis_type))
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        title = open(file_path, 'r').readlines()[0].strip().split("\t")
        print title
        # self.bind_object.logger.info(file_path)
        with open(file_path, 'r') as r:
            # lines = r.readlines()[1:]
            for line in r:
                if line.startswith("#"):
                    continue
                tmp = line.strip().split('\t')
                if analysis_type == 'one_group':
                    insert_data = {
                        "compare_id": compare_id,
                        "chr": tmp[0],
                        "pos": int(tmp[1]),
                        "ref": tmp[2],
                        "alt": tmp[3],
                        "average": float(tmp[4]),
                        "miss": float(tmp[5]),
                        "maf": float(tmp[6]),
                        "frequence": tmp[7]
                    }
                    try:
                        insert_data["annotion"] = tmp[8]
                    except:
                        insert_data["annotion"] = ""
                elif analysis_type == 'two_group':
                    insert_data = {
                        "compare_id": compare_id,
                        "chr": tmp[0],
                        "pos": int(tmp[1]),
                        "ref": tmp[2],
                        "alt": tmp[3],
                        title[4]: float(tmp[4]),
                        title[5]: float(tmp[5]),
                        title[6]: float(tmp[6]),
                        title[7]: tmp[7],
                        title[8]: float(tmp[8]),
                        title[9]: float(tmp[9]),
                        title[10]: float(tmp[10]),
                        title[11]: tmp[11]
                    }
                    try:
                        insert_data["annotion"] = tmp[12]
                    except:
                        insert_data["annotion"] = ""
                else:
                    insert_data = {
                        "compare_id": compare_id,
                        "chr": tmp[0],
                        "pos": int(tmp[1]),
                        "ref": tmp[2],
                        "alt": tmp[3],
                        title[4]: tmp[4],
                        title[5]: tmp[5]
                    }
                    try:
                        insert_data["annotion"] = tmp[-1]
                    except:
                        insert_data["annotion"] = ""
                    if len(title) == 11:
                        insert_data.update({
                            title[6]: tmp[6],
                            title[7]: tmp[7],
                            title[8]: int(tmp[8]),
                            title[9]: int(tmp[9])
                        })
                    elif len(title) == 8:
                        insert_data.update({
                            title[6]: int(tmp[6]),
                        })
                    else:
                        raise Exception("文件{}表头列数不正确！".format(file_path))
                data_list.append(insert_data)
                if len(data_list) == 100000:
                    if types == "indel":
                        self.col_insert_data("sg_indel_compare_detail", data_list)
                    else:
                        self.col_insert_data("sg_snp_compare_detail", data_list)
                    data_list = []
        if types == "indel":
            self.update_db_record("sg_indel_compare", {"main_id": compare_id}, {"diff_variant": file_path})
            if download_path:
                self.update_db_record("sg_indel_compare", {"main_id": compare_id}, {"download_path": download_path})
        else:
            self.update_db_record("sg_snp_compare", {"main_id": compare_id}, {"diff_variant": file_path})
            if download_path:
                self.update_db_record("sg_snp_compare", {"main_id": compare_id}, {"download_path": download_path})
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            if types == "indel":
                self.col_insert_data("sg_indel_compare_detail", data_list)
            else:
                self.col_insert_data("sg_snp_compare_detail", data_list)

    def add_snp_indel_compare(self, task_id, project_sn, types="snp", analysis_type="sample"):
        """
        snp/indel主表
        :param task_id:
        :param project_sn:
        :param types:
        :param analysis_type:sample/one_group/two_group
        :return:
        """
        collection = "sg_indel_compare" if types == "indel" else "sg_snp_compare"
        main_id = self.add_main_table(collection, task_id, project_sn, " ", "origin_{}_compare".format(types),
                                      "{}比较分析主表".format(types), "", is_update="true")
        self.update_db_record(collection, {"_id": main_id}, {"type": analysis_type})
        return main_id

    def add_sg_manhattan(self, compare_id, results_path, types="snp"):
        """
        曼哈顿图导表
        results_path: diff.variant
        """
        compare_id = self.check_objectid(compare_id)  # 检查id是否是OBjectID
        self.check_exists(results_path)  # 检查文件是否存在
        location = "snp_manhattan" if types == "snp" else "indel_manhattan"
        main_data = [{
            "origin_id": compare_id,
            "name": "曼哈顿图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": location,
            "type": 1
        }]
        main_id = self.col_insert_data("sg_manhattan", main_data)
        chr_lists = []
        insert_data = []
        with open(results_path)as fr:
            lines = fr.readlines()
            chr_data = defaultdict(list)  # {染色体：[pos],[index]]}或{染色体：[pos],[index],[index2],[delta_index]}
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                if tmp[0] not in chr_data.keys():
                    chr_data[tmp[0]] = [[], [], [], []]
                chr_data[tmp[0]][0].append(float(tmp[1]))
                chr_data[tmp[0]][1].append(float(tmp[6]))
                chr_data[tmp[0]][2].append(float(tmp[10]))
                delta = abs(float(tmp[6]) - float(tmp[10]))
                chr_data[tmp[0]][3].append(float(delta))
        for key in chr_data:
            chr_lists.append(str(key).strip("\""))
            chrom = str(key).strip("\"")
            data = {
                "manhattan_id":  main_id,
                "chrom": chrom,
                "pos_len": len(chr_data[key][0]),
                "name": "x_categories",
                "value": chr_data[key][0],
            }
            insert_data.append(data)
            data = {
                "manhattan_id":  main_id,
                "chrom": chrom,
                "name": "index_value",
                "value": chr_data[key][1],
            }
            insert_data.append(data)
            if len(chr_data[key]) == 4:  # 当有index2，delta_index时，需要加进去
                data = {
                    "manhattan_id":  main_id,
                    "chrom": chrom,
                    "name": "index2_value",
                    "value": chr_data[key][2],
                }
                insert_data.append(data)
                data = {
                    "manhattan_id":  main_id,
                    "chrom": chrom,
                    "name": "delta_value",
                    "value": chr_data[key][3],
                }
                insert_data.append(data)
        if len(insert_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(results_path))
        else:
            self.col_insert_data("sg_manhattan_detail", insert_data)
            chr_list, sca_list, other_list, final_list = [], [], [], []
            for i in chr_lists:   # 对染色体排序
                if i.startswith("chr"):
                    num = i.strip().split("chr")[-1]
                    chr_list.append(int(num))
                elif i.startswith("sca"):
                    num = i.strip().split("sca")[-1]
                    sca_list.append(int(num))
                else:
                    other_list.append(i)
            chr_list.sort()
            sca_list.sort()
            for i in chr_list:
                final_list.append("chr" + str(i))
            for i in sca_list:
                final_list.append("sca" + str(i))
            for i in other_list:
                final_list.append(i)
            self.update_db_record("sg_manhattan", {"_id": main_id}, {"chr_list": final_list})
            if types == "snp":
                self.update_db_record("sg_snp_compare", {"_id": compare_id}, {"chr_list": final_list})
            else:
                self.update_db_record("sg_indel_compare", {"_id": compare_id}, {"chr_list": final_list})

    def add_sg_distribution_new(self, compare_id, results_path, types="snp"):
        """
        染色体分布图导表
        results_path: win.stat.xls
        """
        compare_id = self.check_objectid(compare_id)  # 检查id是否是OBjectID
        self.check_exists(results_path)  # 检查文件是否存在
        location = "snp_distribution" if types == "snp" else "indel_distribution"
        main_data = [{
            "origin_id": compare_id,
            "name": "染色体分布图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": location,
            "type": 1
        }]
        main_id = self.col_insert_data("sg_distribution", main_data)
        chr_lists = []
        insert_data = []
        win_data = {}
        start, end, s_max = 0, 0, 0
        with open(results_path, "r")as f:
            lines = f.readlines()
            for l in range(len(lines)):
                item = lines[l].strip().split("\t")
                if item[0] not in win_data.keys():
                    chr_lists.append(item[0])
                    win_data[item[0]] = []
                    if s_max > end:
                        end = s_max
                win_data[item[0]].append(int(item[2]))
                s_max = int(item[1])
        chr_list, sca_list, other_list, final_list = [], [], [], []
        for i in chr_lists:   # 对染色体排序
            if i.startswith("chr"):
                num = i.strip().split("chr")[-1]
                chr_list.append(int(num))
            elif i.startswith("sca"):
                num = i.strip().split("sca")[-1]
                sca_list.append(int(num))
            else:
                other_list.append(i)
        chr_list.sort()
        sca_list.sort()
        for i in chr_list:
            final_list.append("chr" + str(i))
        for i in sca_list:
            final_list.append("sca" + str(i))
        for i in other_list:
            final_list.append(i)
        for chr in final_list:
            data = {
                "distribution_id": main_id,
                "name": chr,
                "data_list": win_data[chr]
            }
            insert_data.append(data)
        if len(insert_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(results_path))
        else:
            self.col_insert_data("sg_distribution_detail", insert_data)
            self.update_db_record("sg_distribution", {"_id": main_id}, {"chr_list": final_list, "start": start, "end": end})
        if types == "snp":
            self.update_db_record("sg_snp_compare", {"main_id": compare_id}, {"chr_list": final_list})
        else:
            self.update_db_record("sg_indel_compare", {"main_id": compare_id}, {"chr_list": final_list})

    def add_sg_distribution(self, compare_id, results_path, s_step=10000, types="snp"):
        """
        染色体分布图导表
        results_path: win.stat
        s_step: 滑窗步长
        """
        compare_id = self.check_objectid(compare_id)  # 检查id是否是OBjectID
        self.check_exists(results_path)  # 检查文件是否存在
        # self.bind_object.logger.info("开始导染色体分布图，进行滑窗计算")
        # self.bind_object.logger.info(results_path)
        location = "snp_distribution" if types == "snp" else "indel_distribution"
        main_data = [{
            "origin_id": compare_id,
            "name": "染色体分布图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": location,
            "type": 1
        }]
        main_id = self.col_insert_data("sg_distribution", main_data)
        chr_lists = []
        insert_data = []
        win_data = {}
        with open(results_path, "r")as f:
            lines = f.readlines()
            for l in range(len(lines)):
                item = lines[l].strip().split("\t")
                if item[0] not in win_data.keys():
                    chr_lists.append(item[0])
                    win_data[item[0]] = []
                    i = 0
                    x = 0
                    step = s_step
                if int(item[1]) <= step:
                    x += 1
                else:
                    n = int(item[1]) / s_step
                    win_data[item[0]].append(x)
                    i += 1
                    if n > i:
                        m = i
                        for j in range(m, n):
                            step += s_step
                            x = 0
                            win_data[item[0]].append(x)
                            i += 1
                    x = 1
                    step += s_step
                if l < len(lines) - 1:
                    tmp = lines[l+1].strip().split("\t")
                    if item[0] != tmp[0]:
                        win_data[item[0]].append(x)
            try:
                win_data[item[0]].append(x)
            except:
                self.bind_object.logger.info("{}文件为空！".format(results_path))
        start, end = 0, 0
        chr_list, sca_list, other_list, final_list = [], [], [], []
        for i in chr_lists:   # 对染色体排序
            if i.startswith("chr"):
                num = i.strip().split("chr")[-1]
                chr_list.append(int(num))
            elif i.startswith("sca"):
                num = i.strip().split("sca")[-1]
                sca_list.append(int(num))
            else:
                other_list.append(i)
        chr_list.sort()
        sca_list.sort()
        for i in chr_list:
            final_list.append("chr" + str(i))
        for i in sca_list:
            final_list.append("sca" + str(i))
        for i in other_list:
            final_list.append(i)
        for chr in final_list:
            s_max = len(win_data[chr]) * s_step
            if s_max > end:
                end = s_max
            data = {
                "distribution_id": main_id,
                "name": chr,
                "data_list": win_data[chr]
            }
            insert_data.append(data)
        if len(insert_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(results_path))
        else:
            self.col_insert_data("sg_distribution_detail", insert_data)
            self.update_db_record("sg_distribution", {"_id": main_id}, {"chr_list": final_list, "start": start, "end": end})
        if types == "snp":
            self.update_db_record("sg_snp_compare", {"main_id": compare_id}, {"chr_list": final_list})
        else:
            self.update_db_record("sg_indel_compare", {"main_id": compare_id}, {"chr_list": final_list})

    def add_sg_heatmap(self, compare_id, results_path, types="snp"):
        """
        SNP/indel热图导表
        results_path: diff.matrix
        """
        compare_id = self.check_objectid(compare_id)  # 检查id是否是OBjectID
        self.check_exists(results_path)  # 检查文件是否存在
        location = "snp_distribution" if types == "snp" else "indel_distribution"
        main_data = [{
            "origin_id": compare_id,
            "name": "{}热图".format(types),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": location,
            "type": 1
        }]
        main_id = self.col_insert_data("sg_heatmap", main_data)
        insert_data = []
        columns = []
        with open(results_path)as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                columns.append(item[0])
                value = []
                for s in item[1:]:
                    value.append(int(s))
                data = {
                    "heatmap_id": main_id,
                    "name": item[0],
                    "data_list": value
                }
                insert_data.append(data)
        if len(insert_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(results_path))
        else:
            self.col_insert_data("sg_heatmap_detail", insert_data)
            self.update_db_record("sg_heatmap", {"_id": main_id}, {"columns": columns})

    def update_download_path(self, compare_id, types, download_path):
        """
        更新主表sg_indel_compare/sg_snp_compare的download_path，用于下载表格的链接和下游venn的分析
        """
        compare_id = self.check_objectid(compare_id)  # 检查id是否是OBjectID
        if types == "indel":
            self.update_db_record("sg_indel_compare", {"main_id": compare_id}, {"download_path": download_path})
        else:
            self.update_db_record("sg_snp_compare", {"main_id": compare_id}, {"download_path": download_path})

if __name__ == "__main__":
    a = SnpIndelCompare(None)
    task_id = "wgs_test"
    project_sn = "wgs_test"
    # compare_id = "5b03bf9fa4e1af1482e33207"
    # b1 = a.add_snp_indel_compare(task_id, project_sn)
    # b2 = a.add_snp_indel_compare(task_id, project_sn, "indel")
    # a.add_sg_snp_indel_compare_stat(task_id, b2, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/Ann.stat", "indel")
    # a.add_sg_snp_indel_compare_stat(task_id, b1, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/Ann.stat")
    # a.add_sg_snp_indel_compare_eff_stat(task_id, b1, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/Eff.stat")
    # a.add_sg_snp_indel_compare_eff_stat(task_id, b2, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/Eff.stat", "indel")
    # a.add_sg_snp_indel_compare_detail(b1, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/pop.variant")
    # a.add_sg_snp_indel_compare_detail(compare_id, "/mnt/ilustre/users/sanger-dev/workspace/20180523/Single_tsg_29900_0523141120277624_4024/DoubleGroupCompare/output/diff.variant", "snp", "two_group")
    compare_id = "5b0648e9a4e1af1f77f28b06"
    results_path = "/mnt/ilustre/users/sanger-dev/workspace/20180524/Single_tsg_29900_0524130857538605_2512/DoubleGroupCompare/output/diff.variant.chr_1"
    a.add_sg_manhattan(compare_id, results_path, types="indel")
    # results_path = "/mnt/ilustre/users/sanger-dev/workspace/ 20180414/Single_double_group_compare/DoubleGroupCompare/diff/win.stat"
    # compare_id = "5b03bf9fa4e1af1482e33207"
    # results_path = "/mnt/ilustre/users/sanger-dev/workspace/20180523/Single_tsg_29900_0523175002121620_8874/SingleGroupCompare/output/win.stat"
    # a.add_sg_distribution(compare_id, results_path, types="snp")
