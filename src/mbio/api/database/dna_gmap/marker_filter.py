# -*- coding: utf-8 -*-
# __author__ = 'qing_mei'
# modified 20180619
# modified 20180707
# api

from api_base import ApiBase
import datetime
import re
import json


class MarkerFilter(ApiBase):
    def __init__(self, bind_object):
        """
        dna_gmap项目导表, 遗传标记筛选模块导表
        """
        super(MarkerFilter, self).__init__(bind_object)
        self._project_type = "dna_gmap"

    def add_sg_marker(self, project_sn, task_id, type, params=None, detail_info_path=None, name=None,
                      filter_marker_path=None, chr_list=None):
        """
        sg_marker(遗传标记筛选主表)
        20180621补进去chr_list字段；下面读如文件时更新
        path = filtered.marker_path # binmarker用
        这里缺少chr_list字段；已经在controller内导入主表，这里无须修改
        """
        if params:
            new_params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            new_params = "null"
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "type": type,
            "params": new_params,
            "name": name if name else "origin_sg_marker",
            "desc": "遗传标记筛选主表",
            "detail_info_path": detail_info_path,
            "filtered_marker_path": filter_marker_path,
            "chr_list": chr_list
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_marker", data_list)
        self.update_db_record("sg_marker", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_marker_detail(self, task_id, marker_id, file_path, types):
        """
        1.1的图和表
        :sg_marker_detail  标记分类类型分布表
        根据主表的sg_marker的type的字段来确定sg_marker_detail的类型，如SNP/InDel。
        :sg_bar(柱形图主表)
        :sg_bar_detail(柱形图细节表)
        file_path:/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/1.20180619.test.data/pop.filtered.marker.stat.stat.xls
        """
        call_id_ = self.check_objectid(marker_id)  # 检查id是否是ObjectID
        self.check_exists(file_path)
        with open(file_path, 'r') as r:
            data = r.readlines()    # 文件全部为一个list。list[0]是第一行
            type_list = []
            snpindel_list = {}
            name = "标记分类类型分布表"
            data_list = []
            for m in data:
                line = m.strip().split('\t')
                type_list.append(line[0])
                insert_data = {
                    "task_id": task_id,
                    "marker_id": call_id_,
                    "marker_type": line[0],     # 是type
                }
                if types.lower() == 'snp':
                    if types.lower() not in snpindel_list.keys():
                        snpindel_list['snp'] = []
                    insert_data['snp_num'] = int(line[1])
                    snpindel_list['snp'].append(int(line[1]))
                elif types.lower() == 'indel':
                    if types.lower() not in snpindel_list.keys():
                        snpindel_list['indel'] = []
                    insert_data['indel_num'] = int(line[2])
                    snpindel_list['indel'].append(int(line[2]))
                else:
                    if len(snpindel_list.keys()) == 0:
                        snpindel_list['snp'] = []
                        snpindel_list['indel'] = []
                    insert_data['snp_num'] = int(line[1])
                    insert_data['indel_num'] = int(line[2])
                    snpindel_list['snp'].append(int(line[1]))
                    snpindel_list['indel'].append(int(line[2]))
                data_list.append(insert_data)
            if len(data_list) != 0 and len(type_list) != 0:
                self.col_insert_data("sg_marker_detail", data_list)
                curve_id = self.sg_bar(task_id, marker_id, name, type_list, 1, "markertype_bar", "")
                keys_sort = sorted(snpindel_list.keys(), reverse=True)
                for type in keys_sort:
                    if re.search('snp', type.lower()):  # re.I
                        i = 'SNP Number'
                        # if len(snpindel_list[type]) != 0:
                        self.sg_bar_detail(curve_id, i, snpindel_list[type])
                    if re.search('indel', type.lower()):
                        i = 'InDel Number'
                        # if len(snpindel_list[type]) != 0:
                        self.sg_bar_detail(curve_id, i, snpindel_list[type])
            else:
                self.bind_object.logger.info("sg_marker_detail:{}文件为空".format(file_path))  # 这步没有用 by hd

    def add_sg_distribution(self, task_id, origin_id, results_path, step=100000):
        """
        sg_distribution(染色体分布图主表)
        sg_distribution_detail(染色体分布图细节表)
        origin_id : 主表id
        start end是最长的chr的起始和终止
        # chr1_48700      lmxll   1371    20      23
        results_path = /mnt/ilustre/users/caixia.tian/Develop/10.Gmap/other/matrix/12.marker 过滤后的marker
        """
        origin_id = self.check_objectid(origin_id)  # 检查id是否是OBjectID
        self.check_exists(results_path)  # 检查文件是否存在
        main_data = [{
            "task_id": task_id,
            "origin_id": origin_id,
            "name": "标记筛选染色体分布图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "type": 1,
            "location": "marker_distribution",
        }]
        main_id = self.col_insert_data("sg_distribution", main_data)
        insert_data = []
        chr_list = []   # 存全部chr from "chr_3578"
        chr_allpos = {}
        start = 0
        end = 0
        with open(results_path)as f:
            lines = f.readlines()[1:]
            for line in lines:
                items = line.strip().split("\t")
                item = items[0].split("_")
                if item[0] not in chr_allpos.keys():
                    chr_list.append(item[0])
                    chr_allpos[item[0]] = []  # 定义成list
                try:
                    chr_allpos[item[0]].append(int(item[1]))  # 存chr下的all位置
                except:
                    if len(item) == 3 and item[1] == 'random':
                        chr_allpos[item[0]].append(int(item[2]))
                    else:
                        pass
                try:
                    location = int(item[1])
                except:
                    if len(item) == 3 and item[1] == 'random':
                        location = int(item[2])
                    else:
                        continue
                if start == 0:
                    start = location
                if location > end:
                    end = location  # 存最大值
                if location < start:
                    start = location    # 存最小值
        self.update_db_record("sg_marker", {"_id": origin_id}, {"chr_list": chr_list})
        win_data = {}
        # for chr_one in chr_allpos.keys():
        for chr_one in chr_list:
            s_step = step
            start_ = start
            win_data[chr_one] = []
            i = 0
            chr_allpos[chr_one].sort()  # 数字排序
            for pos_true in chr_allpos[chr_one]:
                # str 和 int 易错！
                if pos_true >= start_ + s_step:
                    # win_data[chr_one].append(i)
                    while pos_true >= start_ + s_step:
                        win_data[chr_one].append(i)
                        i = 0
                        # s_step += step
                        start_ += step
                    # i += 1
                    # win_data[chr_one].append(i)
                else:
                    i += 1
            data = {
                "distribution_id": main_id,
                "name": chr_one,
                "value": win_data[chr_one]
            }
            insert_data.append(data)
        if len(insert_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(results_path))
        else:
            self.col_insert_data("sg_distribution_detail", insert_data)
            self.update_db_record("sg_distribution", {"_id": main_id}, {"chr_list": chr_list,
                                                                        "start": start, "end": end})

    def add_sg_heatmap(self, task_id, heatnap_id, results_path, path=None,
                       type_subtype=None, params=None, number1=800000, number2=400000):
        """
        每20万来拆sg_heatmap的主表的rows----number2
        每100万来拆sg_heatmap_detail的主表的value---number1
        ！！sg_heatmap(热图主表)
        ！！sg_heatmap_detail(热图细节表)
        :location用于区分一个项目内的sg_distribution表的字段
        ：heatnap_id
        path:/mnt/ilustre/users/caixia.tian/Develop/10.Gmap/test2.xls
        ！！ sg_subtype_matrix 分型矩阵主表
        :type_subtype 分型矩阵主表是原始生成的还是是存在记录的参数的flag ["all", "part"]
        数据模式
                sample1 sample2 sample3 (横坐标ID-对应表结构：rows)
        chr_id1 10  12  13
        chr_id2 20  21  23
        chr_id3 14  15  16
        (纵坐标ID-columns) (每一列的数字-data_list)
        """
        type_subtype = 'all'    # 分型矩阵接口删掉，这个参数无用，默认all
        heatnap_id = self.check_objectid(heatnap_id)  # 检查id是否是OBjectID
        self.check_exists(results_path)  # 检查文件是否存在
        location = "marker_subtype_matrix"
        main_data = {
            "task_id": task_id,
            "origin_id": heatnap_id,
            "name": "marker_filter热图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": location,
            "type": 1
        }
        main_id1 = self.db["sg_heatmap"].insert_one(main_data).inserted_id
        main_data = [{
            "task_id": task_id,
            "origin_id": heatnap_id,
            "type": type_subtype,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": "Matrix_{}".format(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")),
            "path": path if path else results_path,     # rere的path
        }]
        if len(main_data) == 0:
            self.bind_object.logger.info("sg_subtype_matrix导主表失败")
        else:
            main_id = self.col_insert_data("sg_subtype_matrix", main_data)  # 导表外加取main_id = _id
            self.update_db_record("sg_subtype_matrix", {"_id": main_id}, {"main_id": main_id})
        insert_data = []
        rows = []
        sample = []
        with open(results_path)as f:
            lines = f.readlines()
            headmap_data = {}
            for num in range(len(lines)):
                row = lines[num].strip().split("\t")
                if num == 0:
                    sample = row[1:]    # 储存横坐标样品ID
                    continue
                rows.append(row[0])     # 存chrids   # 要排序
                for ii in range(len(sample)):
                    if sample[ii] not in headmap_data.keys():
                        headmap_data[sample[ii]] = []
                    headmap_data[sample[ii]].append(int(row[ii + 1]))
            for column in headmap_data.keys():
                data = {
                    "value": headmap_data[column],     # 为样本上该列的数据
                    "heatmap_id": main_id1,
                    "name": column,     # 存主表的横坐标
                }
                insert_data.append(data)
        if len(insert_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(results_path))
        else:
            # 细节表value<100万的时候，不分开导细节表---[主表标记的个数==细节表value的个数]
            if len(rows) <= number1:
                self.col_insert_data("sg_heatmap_detail", insert_data)
            else:
                num = 0
                for k in insert_data:
                    self.split_sg_heatmap_detail(k['value'], k['heatmap_id'], k['name'], int(number1))
                    num += 1
                print("标记筛选插入{}条sg_heatmap成功".format(num))
            # 标记个数<20万的时候，不分开导主表
            if len(rows) <= number2:
                self.update_db_record("sg_heatmap", {"_id": main_id1}, {"rows": rows, "columns": sample},
                                      is_show_log="false")
            else:
                # 此处主表做分开导表，只有第一个导进去的主表的_id关联细节表
                self.split_sg_heatmap(rows, task_id, heatnap_id, int(number2), sample, main_id1)

    def split_sg_heatmap_detail(self, list, heatmap_id, name, number):
        """
        用list按照每number份拆后，导表：sg_heatmap
        """
        num = 0
        split_list = []
        while list[num:num + number]:
            split_list.append(list[num:num + number])
            num += number
        # print(len(split_list))
        insert_data = []
        for i in split_list:
            data = {
                "value": i,     # 为样本上该列的数据
                "heatmap_id": heatmap_id,
                "name": name,     # 存主表的横坐标
            }
            insert_data.append(data)
        self.col_insert_data("sg_heatmap_detail", insert_data)
        return True

    def split_sg_heatmap(self, list, task_id, heatnap_id, number, sample, main_id):
        """
        sg_heatmap_detail细节表拆
        """
        num = 0
        split_list = []
        while list[num:num + number]:
            split_list.append(list[num:num + number])
            num += number
        insert_data = []
        self.update_db_record("sg_heatmap", {"_id": main_id}, {"rows": split_list[0], "columns": sample},
                              is_show_log="false")
        for i in split_list[1:]:
            main_data = {
                "task_id": task_id,
                "origin_id": heatnap_id,
                "name": "marker_filter热图",
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "location": "marker_subtype_matrix",
                "type": 1,
                "rows": i,
                "columns": sample
            }
            insert_data.append(main_data)
        self.col_insert_data("sg_heatmap", insert_data, is_show_log="False")

    def add_sg_marker_stat(self, task_id, marker_id, file_path):
        """
        1.4
        sg_marker_stat  标记数据统计表
        file_path:/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/1.20180619.test.data/marker_info.xls
        /mnt/ilustre/users/sanger-dev/workspace/20180730/Gmap_tsg_31260/output/02.marker_filter/marker_info.xls
        该数据有的行数少于9列的就跳过了，为什么少的话，崔青美回来看下
        """
        call_id_ = self.check_objectid(marker_id)  # 检查id是否是ObjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                try:
                    pos = int(line[2])
                except:
                    pos = int(line[0].split('_')[-1])
                insert_data = {
                    "task_id": task_id,
                    "marker_id": call_id_,
                    "locus": line[0],
                    "chr": line[1],
                    "pos": pos,
                    "variant_type": str(line[3]),
                    "average_depth": float(line[4]),
                    "miss_tatio": float(line[5]),
                    "x2": float(line[6]),
                    "df": int(line[7]),
                    "signif": float(line[8]),
                    "class": line[9]
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_marker_stat", data_list)

    def add_sg_marker_child_stat(self, task_id, marker_id, file_path):
        """
        1.5
        sg_marker_child_stat(子代数据统计表)
        file_path:/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/1.20180619.test.data/sub.generation.xls
        """
        call_id_ = self.check_objectid(marker_id)  # 检查id是否是ObjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                insert_data = {
                    "marker_id": call_id_,
                    "task_id": task_id,
                    "specimen_id": line[0],
                    "average_depth": float(line[1]),
                    "miss_tatio": float(line[2])
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_marker_child_stat", data_list)
            print("sg_marker_child_stat导表成功")

    def get_chrlist_repath(self, task_id):
        """
        获取sg_task内ref_chrlist文件，读取rere路径，在workflow里存chr_list[]
        return：更新sg_marker[chr_list]
        """
        result = self.col_find_one("sg_task", {"task_id": task_id})
        return result['ref_chrlist']

    # def get_update_chrlist(self, main_id, chrlist, filtered_marker_path, detail_info_path):
    #     self.update_db_record("sg_marker", {"_id": main_id}, {"chr_list": chrlist,
    #                                                           "filtered_marker_path": filtered_marker_path,
    #                                                           "detail_info_path": detail_info_path})

    def get_update_chrlist(self, main_id, filtered_marker_path, detail_info_path):
        self.update_db_record("sg_marker", {"_id": main_id}, {"filtered_marker_path": filtered_marker_path,
                                                              "detail_info_path": detail_info_path})

    # def get_child_list(self, generationid):
    #     """
    #     generationid = 5b33874fdde3eea81100002a
    #     """
    #     generationid = self.check_objectid(generationid)
    #     result = self.col_find_one("sg_child_list", {"main_id": generationid})
    #     return result

    def add_sg_subtype_matrix(self, task_id, results_path):
        main_data = [{
            "task_id": task_id,
            "origin_id": "",
            "type": "all",
            "params": "",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": "origin_matrix",
            "path": results_path,
        }]
        if len(main_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(results_path))
        else:
            main_id = self.col_insert_data("sg_subtype_matrix", main_data)  # 导表外加取main_id = _id
            self.update_db_record("sg_subtype_matrix", {"_id": main_id}, {"main_id": main_id})

    # def get_sub_matrix(self, matrixid):
    #     """
    #     matrixid = 5b338ce7dde3eec01100002b
    #     """
    #     matrixid = self.check_objectid(matrixid)
    #     result = self.col_find_one("sg_subtype_matrix", {"main_id": matrixid})
    #     return result

    # def get_sg_marker(self, master_id):
    #     """
    #     main_id = 5b331bc8dde3eec01100002a
    #     """
    #     master_id = self.check_objectid(master_id)
    #     result = self.col_find_one("sg_marker", {"main_id": master_id})
    #     return result


if __name__ == "__main__":
    a = MarkerFilter(None)  # 继承父类，需要赋值。否则self.bind_object.logger.info报错
    member_id = ""
    member_type = 1
    cmd_id = 1
    project_sn = 'test_test_gmap_test'
    task_id = 'test_test_tsanger_30729'
    type = 'snp,indel'
    # detail_info_path = "detail_info_path"
    # marker_upload_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/data/filter.marker_upload.xls"
    # params = "{\"pid\":\"14w_27\",\"mid\":\"16S15\",\"type\":\"snp,indel\",\"pdep\":10,\"odep\":1,\"popt\":\"F2\",\"miss_tatio\":0.3,\"signif\":0.05,\"matrixid\":1,\"generationid\":1,\"marker_upload\":\"/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/data/filter.marker_upload.xls\"}",
    # marker_id = a.add_sg_marker(project_sn, task_id, type, params=params, detail_info_path=detail_info_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180711/MarkerFilter_Marker_tsg_30925_0711163508640764/MarkerFilter/output/pop.filtered.gtype.stat"
    # # a.add_sg_marker_detail(task_id, marker_id, file_path, type)
    # file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180724/Gmap_tsg_31167/MarkerFilter1/output/pop.filtered.marker"
    a.add_sg_distribution(task_id, "5c075ea8326c8c2d6f8b4567", "/mnt/lustre/users/sanger/workspace/20181204/Gmap_sanger_143999/output/02.marker_filter/pop.filtered.marker")
    print "ok"
    # file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180711/MarkerFilter_Marker_tsg_30925_0711163508640764/MarkerFilter/output/pop.filtered.matrix"
    # file_path = '/mnt/ilustre/users/sanger-dev/workspace/20180724/Gmap_tsg_31167/MarkerFilter1/output/pop.filtered.matrix'
    # type_subtype = "all"    # part
    # a.add_sg_heatmap(task_id, marker_id, file_path, path="rere")
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/1.20180619.test.data/marker_info.xls"
    # file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180708/Single_markerfilter_tool8_20180708/MarkerFilter/output/marker_info.xls_3"
    # a.add_sg_marker_stat(task_id, marker_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/1.20180619.test.data/sub.generation.xls"
    # a.add_sg_marker_child_stat(task_id, marker_id, file_path)

# 导表：
## sg_marker
## sg_marker_detail
# sg_bar
# sg_bar_detail
# sg_distribution
# sg_distribution_detail
# sg_heatmap
# sg_heatmap_detail
## sg_subtype_matrix 分型矩阵主表
## sg_marker_stat
## sg_marker_child_stat
