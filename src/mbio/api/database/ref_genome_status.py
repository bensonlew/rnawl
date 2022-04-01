# -*- coding: utf-8 -*-
# __author__ = 'shijin'
from biocluster.api.database.base import Base, report_check
import re
from biocluster.config import Config



class GenomeInfo(Base):
    def __init__(self, bind_object):
        super(GenomeInfo, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + "_ref_rna"

    @report_check
    def add_genome_info(self, file_path, major=True):
        self.bind_object.logger.info("add_genome_info start")
        data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "species_name": "",
                    "ref_anno_version": "",
                    "hyperlink": ""
                }
        try:
            collection = self.db["sg_specimen"]
            result = collection.insert_many(data)
            self.sample_table_ids = result.inserted_ids[:]
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
        else:
            if major == True:
                self.add_genome_info_detail(file_path)
            self.bind_object.logger.info("导入样品信息数据成功:%s" % result.inserted_ids)


    def add_genome_info_detail(self, file_path):

        data_list = []
        with open(file_path, 'r') as f:
            l = f.readline()
            if not l.strip().endswith("Gene_Number"):
                raise Exception("文件%s格式不正确，请选择正确的样品信息文件" % file_path)
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                data = {
                    "project_sn": self.bind_object.sheet.project_sn,
                    "task_id": self.bind_object.sheet.id,
                    "specimen_name": line_data[0],
                    "read_number": int(line_data[1]),
                    "base_number": int(line_data[2]),
                    "average_length": float(line_data[3]),
                    "min_length": int(line_data[4]),
                    "max_length": int(line_data[5]),
                    "is_initial": 1
                }
                data_list.append(data)
        try:
            collection = self.db["sg_specimen"]
            result = collection.insert_many(data_list)
            self.sample_table_ids = result.inserted_ids[:]
        except Exception, e:
            self.bind_object.logger.error("导入样品信息数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样品信息数据成功:%s" % result.inserted_ids)
        return self.sample_table_ids
