# -*- coding: utf-8 -*-
# __author__ = 'gaohao guanqing'


from biocluster.api.database.base import Base, report_check
import json
import datetime
from bson.son import SON
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.param_pack import group_detail_sort

class Regression(Base):


    def __init__(self, bind_object):
        super(Regression, self).__init__(bind_object)
        self._project_type = "meta"

    def add_regression(self,tax_anno_id,func_anno_id,geneset_id,tax_level,func_level,name=None,params=None):
        task_id = self.bind_object.sheet.id
        group_detail = {'All': [str(i) for i in spname_spid.values()]}
        params['group_detail'] = group_detail_sort(group_detail)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "geneset_id": geneset_id,
            "tax_anno_id": tax_anno_id,
            "func_anno_id":func_anno_id,
            "name":  self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "tax_level": tax_level,
            "func_level": func_level,
            "status": "end",
            "desc": "回归分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["regression"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_resgression_site_detail(self,file_path,table_id,MDS_PC):    
        data_list = []
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="51006801")
        with open(file_path, 'rb') as r:
            ##guanqing.zou 20180518 实现从一文件读多个PC 水平的坐标值，以type字段进行区分
            n = 0
            type = ''
            for line in r:
                line = line.strip('\n')
                line_data = line.split('\t')
                if line_data[0] == 'name' :
                    # by houshuang 20191009 >>>
                    if MDS_PC in ["Sobs","Shannon","Simpson","Ace","Chao"]:
                        type = MDS_PC
                    else:
                        n += 1
                        type = MDS_PC+str(n)
                else:
                    data = [("environmental_regression_id", table_id), ("sample_name", line_data[
                        0].strip("\"")), ("X_PCA", line_data[1]), ("Y_factor", line_data[2]),("type", type)]    
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_environmental_regression_curve"]  
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51006802")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    def add_resgression_message_detail(self,file_path,table_id,MDS_PC):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="51006801")
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
                    # by houshuang 20191009 增加fvalue、pvalue、adj_r >>>
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
                    data = [("environmental_regression_id", table_id), ("fvalue", des[0]), ("pvalue", des[1]),
                            ("R_2", des[2]), ("adj_r", des[3]), ("xmin", line_data[4]), ("ymin", line_data[5]),
                            ("xmax", line_data[6]), ("ymax", line_data[7]), ("k", line_data[8]), ("b", line_data[9]),
                            ("lineFormula",line_d)
                         ]
                    if MDS_PC in ["Sobs","Shannon","Simpson","Ace","Chao"]:
                        data.append(("type", MDS_PC))
                    else:
                        data.append(("type", MDS_PC+str(n)))
                    # <<<
                    data_son = SON(data)
                    lab_v += MDS_PC+str(n)+','
                    data_list.append(data_son)
                    i = 0   #guanqing.zou 20180518  实现隔一行取一行信息
        try:
            collection = self.db["sg_environmental_regression_line"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_environmental_regression"]
            main_collection.update({"_id": ObjectId(table_id)},
                                   {"$set": {"status": 'end'}})
            # by houshuang 20191009 >>>
            if MDS_PC in ["Sobs","Shannon","Simpson","Ace","Chao"]:
                lab_v = MDS_PC+"1"
            main_collection.update({"_id": ObjectId(table_id)},
                                   {"$set": {"lab": lab_v[:-1]}})
            main_collection.update({"_id": ObjectId(table_id)},
                                   {"$set": {"version": "4.0"}})
            # <<<
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51006802")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list








