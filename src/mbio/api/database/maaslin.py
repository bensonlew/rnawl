# -*- coding: utf-8 -*-

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import datetime
import re
import pandas as pd


class Maaslin(Base):

    def __init__(self, bind_object):
        super(Maaslin, self).__init__(bind_object)
        self._project_type = 'meta'


    @report_check
    def add_site_detail(self, file_path, table_id=None):
        self.bind_object.logger.info('start insert mongo')

        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="51008701")

        data_list = []
        #sp_mat = re.compile('__[dkpcofgs]*__')
        with open(file_path, 'rb') as r:
            head = r.readline()
            head = head.strip()
            sps = head.split('\t')[2:]
            #sps = [ sp_mat.split(i)[-1] for i in sps]
            for line in r:
                line = line.strip()
                spline = line.split('\t')
                for index,e in enumerate(spline[2:],0):
                    data = {
                        'maaslin_id' : table_id,
                        'spe' : spline[0],
                        'x' : spline[1],
                        'y' : e,
                        'feature' : sps[index],
                        'type' : 'specimen'
                    }
                    data_list.append(data)
        try:
            collection = self.db["sg_maaslin_line"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)


    @report_check
    def add_data_detail(self, file_path, table_id=None):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51008702")

        data_list = []
        line_data = pd.read_table(file_path,sep='\t',header=0)

        for index in range(len(line_data)):
            data = [
                ("maaslin_id", table_id),
                ("variable",line_data['Variable'][index]),
                ("feature",line_data['Feature'][index]),
                ("coefficient",line_data['Coefficient'][index]),
                ("n",line_data['N'][index]),
                ("n_not",line_data['N not 0'][index]),
                ("pvalue",line_data['P-value'][index]),
                ("qvalue",line_data['Q-value'][index])
            ]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["sg_maaslin_detail"]
            collection.insert_many(data_list)
            main_collection = self.db['sg_maaslin']
            #main_collection.update({"_id":table_id},{"$set": {"main_id":table_id}})

        except Exception as e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    @report_check
    def add_line_detail(self, file_path, table_id):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, StringTypes):
                table_id = ObjectId(table_id)
        else:
            self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51008703")

        data_list = []
        line_data = pd.read_table(file_path,sep='\t',header=0)
        #sp_mat = re.compile(';\s*')
        for index in range(len(line_data)):
            formula = "y={0}x+{1}".format(round(line_data['K'][index],4), round(line_data['b'][index],4))
            #spe = re.sub(";\s*","__",line_data['spe'][index])
            #$feature = sp_mat.split(line_data['spe'][index])[-1]
            feature = line_data['spe'][index]
            data = [
                ("maaslin_id", table_id),
                ("feature",  feature),
                ("xmin", line_data['xmin'][index]),
                ("ymin", line_data['ymin'][index]),
                ("xmax", line_data['xmax'][index]),
                ("ymax", line_data['ymax'][index]),
                ("formula",formula),
                ("type", "line")
            ]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["sg_maaslin_line"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list
