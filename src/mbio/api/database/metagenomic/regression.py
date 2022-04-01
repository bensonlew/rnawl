# -*- coding: utf-8 -*-
# __author__ = 'gaohao'


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
        self._project_type = "metagenomic"

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

    def add_resgression_site_detail(self,file_path,table_id, level=None):
        data_list = []
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="52802801")
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("regression_id", table_id), ("sample_name", line_data[
                        0].strip("\"")), ("X_func", line_data[1]), ("Y_taxon", line_data[2])]
                    if level:
                        data.append(('name',level))
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["regression_curve"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入file_path信息出错", code="52802802")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    def add_resgression_message_detail(self,file_path,table_id,level=None):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="52802803")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip("\n")
                    line_data = line.split('\t')
                    line_d = ''
                    #if float('%.3f' % float(line_data[6])) > 0:
                   #     line_d ='y='+str('%.3f' % float(line_data[5]))+'x+'+ str('%.3f' % float(line_data[6]))
                   # elif float('%.3f' % float(line_data[6]))< 0:
                   #     line_d = 'y=' + str('%.3f' % float(line_data[5])) + 'x' + str('%.3f' % float(line_data[6]))
                    line_data[6] = line_data[6].strip()
                    if line_data[6][0] == '-' :
                        line_d = 'y=' + self.format(line_data[5]) + 'x' +  self.format(line_data[6])
                    else:
                        line_d = 'y=' + self.format(line_data[5]) + 'x+' +  self.format(line_data[6])
                    if line_data[0] == 'NA':
                        des= 'NA'
                    else:
                        des ='%.2f' % float(line_data[0])
                    data = [("regression_id", table_id), ("R_2", self.format(des)), ("xmin", self.format(line_data[1])),
                            ("ymin", self.format(line_data[2])), ("xmax", self.format(line_data[3])), ("ymax", self.format(line_data[4])), ("regression_equation",line_d),
                            ('k', self.format(line_data[5])),('b', self.format(line_data[6])),
                            ('raw_r2', self.format(line_data[7])), ('f_value',self.format(line_data[8])),('p_value',self.format(line_data[9]))
                         ]
                    if level:
                        data.append(('name',level))
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["regression_line"]
            collection.insert_many(data_list)
            main_collection = self.db["regression"]
            main_collection.update({"_id": ObjectId(table_id)},
                                   {"$set": {"status": 'end','version':'1'}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入file_path信息出错", code="52802804")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    def format(self, va):
        v = va
        e_ = ''
        if not isinstance(va, str):
            v = str(va)
        if 'e' in v:
            tmp = v.split('e')
            v = tmp[0]
            e_ = tmp[1]
        spv = v.split('.')
        if len(spv) == 2:
            num = 1
            tmp = spv[1][0]
            for j in spv[1][1:]:
                tmp += j
                num +=1
                if num > 3 and int(tmp) !=0:
                    break
            v = spv[0]+'.'+ tmp
        if e_ != '':
            v = v+'e'+e_
        return v










