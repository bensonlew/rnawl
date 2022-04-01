# -*- coding: utf-8 -*-
# __author__ = 'qing_mei'
# modified 20180704
# api

from api_base import ApiBase
import datetime


class QtlExport(ApiBase):
    def __init__(self, bind_object):
        """
        dna_gmap项目导表, 遗传标记筛选模块导表
        """
        super(QtlExport, self).__init__(bind_object)
        self._project_type = "dna_gmap"

    def add_sg_qtl_export(self, project_sn, task_id, path, params, name=None, desc=None):
        """
        sg_qtl_export(格式化软件导出主表)
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": params,
            "name": name if name else "origin_sg_qtl_export",
            "desc": "格式化软件导出主表",
            "export_path": path,
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_qtl_export", data_list)
        self.update_db_record("sg_qtl_export", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_qtl_export_detail(self, qtl_export_id, taskid, sg_lg_id, sg_feature_file_id):
        """
        :sg_qtl_export_detail (格式化软件导出细节表)
        file_path:/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/data/QTl_result/Qtl_export.xls
        无用file_path
        """
        data_list = []
        qtl_export_id = self.check_objectid(qtl_export_id)  # 检查id是否是ObjectID
        tuple_num = self.get_lg_num(taskid, sg_lg_id)
        insert_data = {
            "qtl_export_id": qtl_export_id,
            "marker_number": tuple_num[2],
            "individual_number": tuple_num[0],
            "population_type": self.get_popt_type(taskid),
            "chromosome_number": tuple_num[1],
            "phenotype_number": self.get_trit_num(taskid, sg_feature_file_id)
        }
        data_list.append(insert_data)
        if len(data_list) != 0:
            self.col_insert_data("sg_qtl_export_detail", data_list)
        else:
            raise Exception("sg_qtl_export_detail:{}文件为空！".format(file_path))
            self.bind_object.logger.info("sg_qtl_export_detail:{}文件为空".format(file_path))

    def get_lg_num(self, taskid, sg_lg_id):
        """
        计算样本个数，lg的染色体个数；
        """
        marker_num = 0
        sg_lg_id = self.check_objectid(sg_lg_id)
        # result = self.col_find_one("sg_lg", {"task_id": taskid})
        result = self.col_find_one("sg_lg", {"_id": sg_lg_id})
        indi_num = len(result['specimen_ids'])
        # objectid_id = result['main_id']
        # objectid_id = self.check_objectid(objectid_id)
        result = self.col_find("sg_lg_stat", {"lg_id": sg_lg_id})    # 返回多个子表，为list//
        lg_num = 0
        for i in result:
            lg_num += 1
            marker_num += i['marker_num']
        return indi_num, lg_num, marker_num

    def get_popt_type(self, taskid):
        """
        根据task_id从sg_task里取出群体类型
        """
        result = self.col_find_one("sg_task", {"task_id": taskid})
        poptype = result['poptype']
        return poptype

    def get_trit_num(self, taskid, sg_feature_file_id):
        """
        从sg_feature_file内取性状个数
        """
        sg_feature_file_id = self.check_objectid(sg_feature_file_id)
        # result = self.col_find_one("sg_feature_file", {"task_id": taskid})
        result = self.col_find_one("sg_feature_file", {"_id": sg_feature_file_id})
        trit_dict = result['feature_dict']
        trit_keys = trit_dict.keys()
        trit_num = len(trit_keys)
        return trit_num

    def updata_sg_qtl_export(self, main_id, export_path):
        main_id = self.check_objectid(main_id)
        # result = self.col_find_one("sg_qtl_export", {"main_id": main_id})
        self.update_db_record("sg_qtl_export", {"_id": main_id}, {"export_path": export_path})


    # def get_mark_lg_path(self, sg_lg_id):
    #     """
    #     以下均给tool返回结果用
    #     返回sg_lg的结果：
    #         CP群体返回sexAver_loc_path、sexAver_map_path
    #         NPCP群体返回total_csv_path、total_map_path、total_loc_path
    #     taskid =
    #     """
    #     sexAver_map = ''
    #     sexAver_loc = ''
    #     sg_lg_id = self.check_objectid(sg_lg_id)
    #     print(sg_lg_id)
    #     result = self.col_find_one("sg_lg", {"main_id": sg_lg_id})
    #     print(result)
    #     if 'sexAver_map' in result:   # F1
    #         sexAver_loc = result['sexAver_loc_path']     # total.sexAver.loc
    #         sexAver_map = result['sexAver_map_path']     # total.sexAver.map
    #         return sexAver_loc, sexAver_map
    #     else:   # NO F1
    #         total_csv = result['total_csv_path']        # total.csv 7.6号修改为total.
    #         total_loc = result['total_loc_path']
    #         total_map = result['total_map_path']
    #         return total_csv, total_loc, total_map

    # def get_nocp_trit(self, sg_feature_file_id):
    #     """
    #     传sg_feature的id，得到sg_feature的样本id和sg_feature_file_detail的性状
    #     """
    #     sg_feature_file_id = self.check_objectid(sg_feature_file_id)
    #     master_result = self.col_find_one("sg_feature", {"_id": sg_feature_file_id})
    #     print("^^^^^^^^^^^^^^")
    #     print(master_result)
    #     print(sg_feature_file_id)
    #     file_result = self.col_find_one('sg_feature_file', {"task_id": master_result['task_id']})
    #     detail_result = self.col_find("sg_feature_file_detail", {"feature_id": master_result['main_id']})
    #     feature_dict = file_result['feature_dict']  # feature_dict['trait1']还原本来的名字；对应关系
    #     # sample_list = master_result['specimen_list']
    #     data_dict = {}
    #     sample_list = []
    #     print("$$$$$$$$$$$$")
    #     print(detail_result)    # <pymongo.cursor.Cursor object at 0x7fd0c8184710> 多数据
    #     print(feature_dict)
    #     for record in detail_result:
    #         print(record)
    #         sample = record['sample']
    #         if sample:
    #             if sample not in sample_list:
    #                 sample_list.append(sample)
    #         for traits in record:     # 读每一条表记录
    #             if traits in feature_dict.keys():   # keys里是trait1
    #                 trait = feature_dict[traits]
    #                 if sample not in data_dict.keys():
    #                     data_dict[sample] = {}
    #                 if trait not in data_dict[sample].keys():
    #                     data_dict[sample][trait] = ''
    #                 data_dict[sample][trait] = str(float(record[traits]))
    #     return data_dict

    # def get_cp_trit(self, sg_feature_file_id):
    #     """
    #     demo_file_path = /mnt/ilustre/users/qingmei.cui/newmdt/sanger/6.gmap/QTl_convert/CP/trit.xls
    #     备注：测试先用，等曾静那边搞定CP的运行程序后，确定数据样子。
    #     返回的数据不能是路径，是路径的话，要用workflow包装后然后传给tool
    #     """
    #     trit_path = "/mnt/ilustre/users/qingmei.cui/newmdt/sanger/6.gmap/QTl_convert/CP/trit.xls"
    #     return trit_path


if __name__ == "__main__":
    a = QtlExport(None)  # 继承父类，需要赋值。否则self.bind_object.logger.info报错
    member_id = ""
    member_type = 1
    cmd_id = 1
    project_sn = 'gmap_test'
    task_id = 'tsanger_30729'
    file_path = "data/rerewrweset/files/m_188/info_1520321510.list"     # 压缩文件路径
    params = "{\"sg_lg_id\": \"5b3b0230a4e1af0d52ce16ba\",\"sg_feature_file_id\":\"5b2c68aaa4e1af2ac8a1db82\", \"type\": [\"rqtl\", \"mapqtl\", \"qtlcart\"]}"  # controller导主表
    marker_id = a.add_sg_qtl_export(project_sn, task_id, file_path, params)
    a.add_sg_qtl_export_detail(marker_id, file_path)

# 导表：
# sg_qtl_export (格式化软件导出主表)
# sg_qtl_export_detail (格式化软件导出细节表)

# marker_number: 多个表，sg_lg_stat['marker_num']的和；
# individual_number: 单个表，sg_lg['specimen_ids']，类型：list；
# population_type: 单个表，sg_task['poptype']；
# chromosome_number: 多个表，sg_lg_stat['lg'] keys的和；
# phenotype_number: 单个表，sg_feature_file['feature_dict']，类型：dict, 需要将dict的value去掉-1，-2来去重。
