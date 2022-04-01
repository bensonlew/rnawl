# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.14

from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
import datetime
import os


class TraitAnnovar(Base):
    def __init__(self, bind_object):
        """
        性状分析导表
        """
        super(TraitAnnovar, self).__init__(bind_object)
        self._project_type = "dna_gmap"

    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        :param id_:
        :return:
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        :param file_path:
        :return:
        """
        if not os.path.exists(file_path):
            raise Exception("文件或文件夹{}不存在！".format(file_path))

    def add_sg_feature(self, project_sn, task_id, params=None, name=None):
        """
        sg_feature
        """
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_feature_annovar",
            "params": params if params else "null",
            "desc": "性状分析主表"
        }
        main_id = self.db["sg_feature"].insert_one(insert_data).inserted_id
        self.db["sg_feature"].update_one({"_id": main_id}, {"$set": {"main_id": main_id}})
        return main_id

    def add_sg_feature_annovar(self, task_id, feature_id, anno_trait):
        """
        sg_feature_annovar
        anno_trait: annovar.trait.xls
        """
        feature_id = self.check_objectid(feature_id)
        self.check_exists(anno_trait)
        data_list = []
        categories = []
        # trait_value = []
        with open(anno_trait, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                categories.append(item[0])
                # trait_value.append(int(item[1])-1)
                insert_data = {
                    "feature_id": feature_id,
                    "trait": item[0],
                    "df": int(item[1])
                }
                try:
                    insert_data["f_value"] = float(item[4])
                    insert_data["p_value"] = float(item[5])
                    insert_data["ms"] = float(item[2])
                    insert_data["variance"] = float(item[3])
                except:
                    insert_data["f_value"] = item[4]
                    insert_data["p_value"] = item[5]
                    insert_data["ms"] = item[2]
                    insert_data["variance"] = item[3]
                data_list.append(insert_data)
        if data_list:
            self.db["sg_feature_annovar"].insert_many(data_list)
            # bar_id = self.add_sg_bar(feature_id, task_id, categories)
            # self.add_sg_bar_detail(bar_id, "trait", trait_value)
        else:
            self.bind_object.logger.info("annovar.trait.xls结果为空，请检查")

    def add_sg_feature_bar(self, feature_id, task_id, bar_dir):
        feature_id = self.check_objectid(feature_id)
        self.check_exists(bar_dir)
        for f in os.listdir(bar_dir):
            f_ = os.path.join(bar_dir, f)
            trit = f.split(".bar.txt")[0]
            categories, value = [], []
            with open(f_, "r") as f1:
                lines = f1.readlines()
                for line in lines:
                    item = line.strip().split("\t")
                    categories.append(item[0])
                    value.append(int(item[1]))
            bar_id = self.add_sg_bar(feature_id, task_id, categories, trit)
            self.add_sg_bar_detail(bar_id, trit, value)

    def add_sg_bar(self, feature_id, task_id, categories, name):
        """
        sg_bar
        """
        insert_data = {
             "origin_id": feature_id,
             "task_id": task_id,
             "location": "feature_bar",
             "type": 1,
             "categories": categories,
             "name": name,
             "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        bar_id = self.db["sg_bar"].insert_one(insert_data).inserted_id
        return bar_id

    def add_sg_bar_detail(self, bar_id, name, values):
        """
        sg_bar_detail
        """
        bar_id = self.check_objectid(bar_id)
        insert_data = {
             "bar_id": bar_id,
             "name": name,
             "value": values
        }
        self.db["sg_bar_detail"].insert_one(insert_data)

    def add_sg_feature_samples(self, feature_id, anno_sample):
        """
        sg_feature_samples
        anno_sample: annovar.sample.xls
        """
        feature_id = self.check_objectid(feature_id)
        self.check_exists(anno_sample)
        data_list, specimen_list = [], []
        trait_dict = {}
        with open(anno_sample, "r") as f:
            lines = f.readlines()
            head = lines[0].strip().split("\t")
            for i in range(1, len(head)):
                trait_dict["trait" + str(i)] = head[i]
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "feature_id": feature_id,
                    "sample": item[0],
                }
                specimen_list.append(item[0])
                for i in range(1, len(head)):
                    try:
                        insert_data["trait" + str(i)] = float(item[i])
                    except:
                        insert_data["trait" + str(i)] = item[i]
                data_list.append(insert_data)
        specimen_list = list(set(specimen_list))
        if data_list:
            self.db["sg_feature_samples"].insert_many(data_list)
        else:
            self.bind_object.logger.info("annovar.sample.xls文件为空，请检查")
        self.db["sg_feature"].update({"main_id": feature_id}, {"$set": {"feature_dict": trait_dict, "specimen_list": specimen_list}})


if __name__ == "__main__":
    a = TraitAnnovar(None)
    project_sn = "gmap_test"
    task_id = "tsanger_30729"
    anno_trait = "/mnt/ilustre/users/sanger-dev/workspace/20180614/Single_trait_annovar/TraitAnnovar/output/annovar.trit.xls"
    anno_sample = "/mnt/ilustre/users/sanger-dev/workspace/20180614/Single_trait_annovar/TraitAnnovar/output/annovar.sample.xls"
    feature_id = a.add_sg_feature(project_sn, task_id)
    a.add_sg_feature_annovar(task_id, feature_id, anno_trait)
    a.add_sg_feature_samples(feature_id, anno_sample)
