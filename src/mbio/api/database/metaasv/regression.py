# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from mainapp.libs.param_pack import group_detail_sort, param_pack
from mbio.packages.metaasv.common_function import find_group_name
from biocluster.api.database.base import Base, report_check
import json
import datetime
from bson.son import SON
from bson.objectid import ObjectId
from types import StringTypes

class Regression(Base):
    def __init__(self, bind_object):
        super(Regression, self).__init__(bind_object)
        self._project_type = "metaasv"

    def add_regression(self,asv_id,name=None,params=None,group_id=None,spname_spid=None):
        task_id = self.bind_object.sheet.id
        if spname_spid and params:
            if group_id not in [None, "All", "all", "ALL"]:
                ## 调用common模块，功能将导入的分组方案返回group_detail
                group_detail = find_group_name(task_id)
            else:
                group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
            params['level_id'] = int(9)
            params = param_pack(params)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "asv_id": asv_id,
            "name":  self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "status": "end",
            "desc": "回归分析",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["regression"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_resgression_site_detail(self,file_path,table_id,MDS_PC):
        """
        回归分析的data表
        :param file_path:
        :param table_id:
        :param MDS_PC:
        :return:
        """
        data_list = []
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！")
        with open(file_path, 'rb') as r:
            n = 0
            type = ''
            for line in r:
                line = line.strip('\n')
                line_data = line.split('\t')
                if line_data[0] == 'name' :##标题行的数据
                    if MDS_PC in ["Sobs","Shannon","Simpson","Ace","Chao"]:
                        type = MDS_PC
                    else:
                        n += 1
                        type = MDS_PC+str(n)
                else:##非标题行的数据
                    data = [("regression_id", table_id),
                            ("specimen", line_data[0].strip("\"")),
                            ("x", line_data[2]),
                            ("y", line_data[1]),
                            ("estimators", type)]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["regression_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        try:
            main_collection = self.db["regression"]
            scatter_data = {
                "scatter_data": {"name":"specimen",
                        "category": "",
                        "condition": {"type":"estimators"}
                }
            }
            scatter_data_json = json.dumps(scatter_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": ObjectId(table_id)},{"$set": {"status": 'end',
                                                                         "scatter_data":scatter_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("更新主表regression信息出错")
        else:
            self.bind_object.logger.info("更新主表regression：%s信息成功!" % file_path)
        return data_list, table_id

    def add_resgression_message_detail(self,file_path,table_id,MDS_PC):
        """
        回归分析的message
        :param file_path:
        :param table_id:
        :param MDS_PC:
        :return:
        """
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            n = 0
            lab_v = ''
            for line in r:
                if i == 0:
                    i = 1
                else:
                    n += 1
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    line_d = ''
                    # 增加fvalue、pvalue、adj_r >>>
                    if float('%.3f' % float(line_data[-1])) > 0:
                        line_d ='y='+str('%.3f' % float(line_data[-2]))+'x+'+ str('%.3f' % float(line_data[-1]))
                    elif float('%.3f' % float(line_data[-1]))< 0:
                        line_d = 'y=' + str('%.3f' % float(line_data[-2])) + 'x' + str('%.3f' % float(line_data[-1]))
                    elif float('%.3f' % float(line_data[-1])) == 0.000:
                        line_d = 'y=' + str('%.3f' % float(line_data[-2])) + 'x'
                    des = []
                    for i in range(4):
                        if line_data[i] == 'NA':
                            des.append('NA')
                        else:
                            des.append('%.6f' % float(line_data[i]))
                    data = [("regression_id", table_id),
                            ("fvalue", des[0]),
                            ("pvalue", des[1]),
                            ("R_2", des[2]),
                            ("adj_r", des[3]),
                            ("xmin", line_data[4]),
                            ("ymin", line_data[5]),
                            ("xmax", line_data[6]),
                            ("ymax", line_data[7]),
                            ("k", line_data[8]),
                            ("b", line_data[9]),
                            ("lineFormula",line_d)
                         ]
                    if MDS_PC in ["Sobs","Shannon","Simpson","Ace","Chao"]:
                        data.append(("estimators", MDS_PC))
                    else:
                        data.append(("estimators", MDS_PC+str(n)))
                    data_son = SON(data)
                    lab_v += MDS_PC+str(n)+','
                    data_list.append(data_son)
                    i = 0
        try:
            collection = self.db["regression_line"]
            collection.insert_many(data_list)
            self.bind_object.logger.info("导入regression_line成功！")
            main_collection = self.db["regression"]
            settled_params = {"software" : "R-3.3.1 (stat)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            table_data = {"table_data": ["estimators", "fvalue", "pvalue", "R_2", "adj_r", "xmin", "ymin", "xmax", "ymax", "k", "b"]}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            line_data = {"line_data": {"name":"estimators","condition": {"type":"estimators"}}}
            line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
            text_data = {"line_data": {"name":"pvalue","condition": {"type":"estimators"}}}
            text_data_json = json.dumps(text_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": ObjectId(table_id)},
                                   {"$set": {"status": 'end',
                                            "main_id": table_id,
                                             "settled_params": settled_params_json,
                                             "table_data": table_data_json,
                                             "line_data": line_data_json,
                                             "text_data": text_data_json}})
            if MDS_PC in ["Sobs","Shannon","Simpson","Ace","Chao"]:
                lab_v = MDS_PC+"1"
            main_collection.update({"_id": ObjectId(table_id)},
                                   {"$set": {"lab": lab_v[:-1]}})

        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list