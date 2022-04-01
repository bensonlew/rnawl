# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
from biocluster.api.database.base import Base, report_check
import datetime
import os
from bson import ObjectId
from bson.son import SON
from biocluster.config import Config

class HmdbMap(Base):
    def __init__(self, bind_object):
        super(HmdbMap, self).__init__(bind_object)
        sanger_type, sanger_path = self.bind_object._sheet.output.split(':')
        sanger_prefix = Config().get_netdata_config(sanger_type)
        self.work_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.work_dir = self.bind_object.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'hmdb'
        else:
            self._db_name = 'hmdb'
        self._project_type = 'hmdb'
        self.check()

    @report_check
    def run(self):
        """
        运行函数
        """
        main_id = self.insert_table(self.work_dir + '/annotation_result.xls', 'Map分析结果表')
        self.insert_detail(main_id, self.work_dir + '/annotation_result.xls')

    def insert_table(self, path, name):
        self.bind_object.logger.info('开始导入table表')
        return_id = self.db['cgc_map'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name=name,
            demo=0,
            status='end',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            position=path
        )).inserted_id
        self.bind_object.logger.info('table表导入结束')
        return return_id

    def insert_detail(self, main_id, table_path):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        data_list = list()
        with open(table_path, "r") as file:
            lines = file.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = [
                    ("map_id", main_id),
                    ("gene_id",line[0]),
                    ("value", line[1]),
                    ("taxon_domain", line[2]),
                    ("taxon_kingdom", line[3]),
                    ("taxon_phylum", line[4]),
                    ("taxon_class", line[5]),
                    ("taxon_order", line[6]),
                    ("taxon_family", line[7]),
                    ("taxon_genus", line[8]),
                    ("taxon_species", line[9]),
                    ("egg_anno", line[10]),
                    ("egg_desc", line[11]),
                    ("egg_fun2", line[12]),
                    ("egg_fun1", line[13]),
                    ("kegg_anno", line[14]),
                    ("kegg_desc", line[15]),
                    ("kegg_pathway", line[16]),
                    ("kegg_fun2", line[18]),
                    ("kegg_fun1", line[17])
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["cgc_map_detail"]
        try:
            detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.logger.info("map detail 表导入失败")
        else:
            self.bind_object.logger.info("map detail 表导入成功")


    def check(self):
        """
        检查文件格式是否正确
        """
        pass