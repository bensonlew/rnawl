# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20190902

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re, os
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class Hcluster(Base):
    def __init__(self, bind_object):
        super(Hcluster, self).__init__(bind_object)
        self._project_type = "bac_comparative"

    @report_check
    def add_hcluster(self, params=None, name=None):
        """
        pan_genomes的主表；params中记录了group_detail等字段
        :param params:
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Hcluster分析主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Hcluster_Origin",
        }
        collection = self.db["hcluster"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    def add_hcluster_data(self, inserted_id, hcluster_data=None):
        """
        导表的cluster表
        :param inserted_id: 主表id
        :param pca_data: hcluster的dir文件夹
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        if hcluster_data == None:
            self.bind_object.set_error("没有生成正确的文件夹hcluster_data: %s" % (hcluster_data))

        tree_path = os.path.join(hcluster_data, 'hcluster.tre')
        site_path = os.path.join(hcluster_data, 'hcluster_matrix.xls')
        try:
            self.insert_table_detail(new_inserted_id, tree_path, 'tree')
            self.insert_table_detail(new_inserted_id, site_path, 'specimen')
        except Exception, e:
            self.bind_object.set_error("更新hcluster_detail出错:%s" % (e))


    def insert_table_detail(self, inserted_id, table_path, type):
        """
        插入详情表
        :param inserted_id: 主表id
        :param table_path: 插入主表信息
        :return:
        """
        if type in ["tree"]:
            with open(table_path, 'r') as f:
                lines = f.readlines()
                line = lines[0].strip()
                tree = line
                raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):[0-9])', tree)
                sample_list = [i[1].split("; ")[-1].strip() for i in raw_samp]
            main_collection = self.db["hcluster"]
            try:
                main_collection.update_one({'_id': inserted_id}, {'$set': {'tree': tree,
                                                                           "spe_sort": ",".join(sample_list)}})
            except Exception, e:
                self.bind_object.set_error("更新hcluster主表出错：%s" % (e))
        else:
            with open(table_path, "r") as f:
                lines = f.readlines()
                header = lines[0].strip().split("\t")
                data_list = []
                specimen_list = header[0:]

                num = 1
                sp_dict = {}
                new_specimen_list = []
                for sp_name in specimen_list:
                    spname = "{}:vs{}".format(sp_name, str(num))
                    new_specimen_list.append(spname)
                    sp_dict[sp_name] = "vs{}".format(str(num))
                    num += 1

                for line in lines[1:]:
                    line = line.strip().split("\t")
                    sample_name = line[0]
                    data = {
                        "hcluster_id": inserted_id,
                        "specimen_id": sample_name,
                        "type": type,
                        }
                    for sp in specimen_list:
                        sp_index = specimen_list.index(sp)
                        sp_replace = sp_dict[sp]
                        data[sp_replace] = float(line[sp_index+1])
                    data_son = SON(data)
                    data_list.append(data_son)
            try:
                collection = self.db["hcluster_detail"]
                collection.insert_many(data_list)
                self.bind_object.logger.info("导入hcluster_detail详情表成功！")
            except Exception, e:
                self.bind_object.set_error("导入hcluster_detail结果表出错:%s" % (e))
            if type in ["specimen"]:
                try:
                    main_collection = self.db["hcluster"]
                    main_collection.update_one({'_id': inserted_id}, {'$set': {'spe_list': "|".join(new_specimen_list), 'main_id':inserted_id}})
                except Exception, e:
                    self.bind_object.set_error("导入hcluster_detail结果表出错:%s" % (e))